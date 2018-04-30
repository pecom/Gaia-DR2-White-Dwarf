import numpy as np
import scipy.stats as sp
from scipy.optimize import least_squares
from WhiteDwarf import WhiteDwarf
import matplotlib
matplotlib.use('TkAgg')
from lmfit import minimize, Parameters, Minimizer
import matplotlib.pyplot as plt
import corner
import time
import os

#USING CGS UNITS FOR STELLAR MASS, RADIUS, AND DISTANCES

bergZeroes = [3.684e-9, 6.548e-9, 3.804e-9, 2.274e-9, 1.119e-9, 3.106e-10, 1.143e-10, 4.206e-11, 1.1436e-8, 4.9804e-9, 2.8638e-9, 1.9216e-9, 1.3343e-9, 3.699e-9]
bergeronOrder = ["U", "B", "V", "R", "I", "J", "H", "K", "u", "g", "r", "i", "z", "y"]
bergeronIndices = [5+a for a in range(len(bergeronOrder))]
wavelengthDictionary = {"U": 3971, "B": 4481, "V": 5423, "R": 6441, "I": 8071, "J": 12350, "H": 16620, "K": 21590, "u": 3146, "g": 4670, "r": 6156, "i": 7471, "z": 8918, "y": 5456}
pcInCm = 3.086e18
universalG = 6.672e-8
solarMassinGram = 1.989e33
solarRadiusinCm = 6.96e10
whiteDwarfs = []

def readFluxFrom(fileName):
    """ Reads flux from a given file name.
    Each line of the file should refer to one star and the data should be put in the order of:
    name, parallax (in arcseconds), parallax error (arcseconds), f_λ (erg cm^-2 s^-1 Å^-1), f_λ error, f_λ, f_λ error, ..., [Filter band, Filter band, ....]
    Check whiteDwarfList for an example
    """
    f = open(fileName, 'r')
    for a in f:
        star = WhiteDwarf()
        roughData = a.split(', ')
        noNewLine = []
        for a in roughData:
            a = a.replace('\n', '')
            noNewLine.append(a)
        name = noNewLine[0]
        parallax = float(noNewLine[1])
        paraErr = float(noNewLine[2])
        star.set_name(name)
        star.set_para(parallax, paraErr)
        noName = noNewLine[3:]
        numbers = []
        endLabel = []
        dictionary = False
        for a in noName:
            if a[0] == '[':
                numbers.extend(noName[:noName.index(a)])
                endLabel.extend(noName[noName.index(a):])
                break
        flux = []
        error = []
        labels = []
        for a in range(int(len(numbers)/2)):
            flux.append(float(numbers[2*a]))
            error.append(float(numbers[2*a+1]))
        for a in endLabel:
            labels.append(a.translate(str.maketrans('[]', '  ')).replace(' ', ''))
        star.set_fluxes(np.array(flux), np.array(error))
        star.set_order(labels)
        whiteDwarfs.append(star)

def interpolateScale(start, fin, scale):
    """Simple linear interpolation given a start point, end point, and the amount to scale"""
    return (np.array(fin) - np.array(start))*scale + np.array(start)

def interpolatePoint(findT, findG, tRange, gRange, berg, capture):
    """Calculates the model fluxes for a given effective temp and log(g)
    Also requires a list of possible temperature values, log(g) values, the Bergeron model, and a list of filters to capture (all of these are lists)
    """
    teff = list(tRange)
    logg = list(gRange)
    teff.append(findT)
    logg.append(findG)
    sortedTemp = sorted(teff, key = float)
    sortedLogg = sorted(logg, key = float)
    x = [findT, findG]
    TeffInterval = [sortedTemp[sortedTemp.index(x[0]) - 1], sortedTemp[sortedTemp.index(x[0]) + 1]]
    loggInterval = [sortedLogg[sortedLogg.index(x[1]) - 1], sortedLogg[sortedLogg.index(x[1]) + 1]]

    for a in berg:              #x-axis is teff and y-axis is logg, the names should then make sense
        if a[0] == TeffInterval[0] and a[1] == loggInterval[0]:
            bottomLeft = a
        if a[0] == TeffInterval[1] and a[1] == loggInterval[0]:
            bottomRight = a
        if a[0] == TeffInterval[0] and a[1] == loggInterval[1]:
            topLeft = a
        if a[0] == TeffInterval[1] and a[1] == loggInterval[1]:
            topRight = a

    teffScale = (findT - TeffInterval[0])/(TeffInterval[1] - TeffInterval[0])
    loggScale = (findG - loggInterval[0])/(loggInterval[1] - loggInterval[0])
    bottomPoints = [careAbout(bottomLeft, capture), careAbout(bottomRight, capture)]
    topPoints = [careAbout(topLeft, capture), careAbout(topRight, capture)]
    bottomG = interpolateScale(bottomPoints[0], bottomPoints[1], teffScale)
    topG = interpolateScale(topPoints[0], topPoints[1], teffScale)
    foundIt = interpolateScale(bottomG, topG, loggScale)

    return foundIt

def careAbout(bergRow, capture):
    """ Returns the fluxes of bergeron's table we are interested in """
    return [bergRow[bergeronIndices[i]] for i in capture]

def residual(params, temps, loggs, logg, berg, capture, data, eps):
    """ Residual function for the Levenberg-Marquardt fit
    The fit parameters are effective temperature and radius/distance
    Calls interpolatePoint() so it has a lot of those arguments in this function's arguments
    """
    teff = params['teff']
    rd = params['rd']

    # teff = params[0]
    # rd = params[1]

    model = np.array(interpolatePoint(teff, logg, temps, loggs, berg, capture))*(rd**2)
    return (model - data)/eps

def findTheFit(wdStar, logg, bergTable, teffR, loggR, against):
    """ Finds the model fit for a given white Dwarf
    Uses LMFit package to do Levenberg-Marquardt fit.
    Calls interpolatePoint() so it has a lot of those arguments in this function's arguments
    """
    print("Fitting for logg: %f" % logg)
    captureList = []
    for a in wdStar.labels:
        captureList.append(bergeronOrder.index(a))

    params = Parameters()
    params.add('teff', value = 1.569e4, min = 1.5e3, max = 1.2e5, vary = True)
    params.add('rd', value = 3e-12, min = 0, max = 1, vary = True)

    tatooine = Minimizer(residual, params, fcn_args = (teffR, loggR, logg, bergTable, captureList, against, wdStar.fluxErr))
    locMinima = tatooine.leastsq(xtol = 3e-16, ftol = 3e-16)

    #locMinima = minimize(residual, params, method = "leastsq", args = (teffR, loggR, logg, bergTable, captureList, against, wdStar.fluxErr))

    print(locMinima.params.pretty_print())
    print("Chi sqr is %f: " % locMinima.chisqr)
    return [locMinima.params['teff'].value, locMinima.params['rd'].value, locMinima.chisqr, locMinima.params['teff'].stderr, locMinima.params['rd'].stderr]

def startFits(wdStar):
    """ Calls findTheFit() for each white dwarf 5 times
    For each iteration a log(g) is calculated based on the mass and radius derived from the previous iteration
    The first iteration assumes log(g) = 8
    The temperature, radius, mass, log(g), and chi-squared are stored into the WhiteDwarf object
    """

    print("****************************************************")
    print("Finding the fit for %s" % wdStar.name)
    print("****************************************************")
    logg = 8.0
    f = open('bergeronFlux.csv', 'r')
    bergTable = []
    for a in f:
        bergTable.append(list(map(float, a.split(', '))))
    teffR = []
    loggR = [7.0, 7.5, 8.0, 8.5, 9.0, 9.5]
    for a in bergTable:
        if not (a[0] in teffR):
            teffR.append(a[0])

    mrBerg = massRadiusTable()[0]
    massOrder = massRadiusTable()[1]
    for a in range(5):
        temp, rd, chiVal, teffErr, rdErr = findTheFit(wdStar, logg, bergTable, teffR, loggR, wdStar.observedFluxes)
        print("Temp is %f" % temp)
        wdDistance = (1/wdStar.parallax)*pcInCm
        radius = rd*wdDistance
        print("Radius is %r" % (radius*1.e-5))
        smolMass = getMassFromRT(radius, temp, mrBerg, massOrder)
        print("Mass is %f" % (smolMass/solarMassinGram))
        print("\n")
        logg = np.log10(universalG*smolMass/(radius**2))

        wdStar.set_chi(chiVal)
        wdStar.set_Teff(temp)
        wdStar.set_logg(logg)
        wdStar.set_mass(smolMass)
        wdStar.set_radius(radius)

        calculateErrors(wdStar, rd, rdErr, teffErr)

def fluxGraph(wdStar):
    """ Creates a graph of observed fluxes with error bar in blue vs wavelength (Å) with the model points overlayed in green
    The effective temperature, log(g), radius, mass and chi-squared are put on the graph at the point (2000, 0)
    The images are saved under 'images/' + WhiteDwarfName + '_flux.png'
    """

    rd = wdStar.radius*wdStar.parallax/pcInCm
    teff = wdStar.T_eff
    logg = wdStar.logg

    print(teff)
    print(logg)

    captureList = []
    for a in wdStar.labels:
        captureList.append(bergeronOrder.index(a))

    print(captureList)

    f = open('bergeronFlux.csv', 'r')
    bergTable = []
    for a in f:
        bergTable.append(list(map(float, a.split(', '))))
    teffR = []
    loggR = [7.0, 7.5, 8.0, 8.5, 9.0, 9.5]
    for a in bergTable:
        if not (a[0] in teffR):
            teffR.append(a[0])



    model = np.array(interpolatePoint(wdStar.T_eff, wdStar.logg, teffR, loggR, bergTable, captureList))*(rd**2)
    actual = wdStar.observedFluxes
    errorBar = wdStar.fluxErr

    xticks = [wavelengthDictionary[a] for a in wdStar.labels]
    print(xticks)

    obi = "Teff: %f Log(g): %f \n" % (wdStar.T_eff, wdStar.logg)
    wan = "Radius: %f Mass: %f  \n" % (wdStar.radius/solarRadiusinCm, wdStar.mass/solarMassinGram)
    kenobi = "Chi-squared: %f" % wdStar.chisqr

    plt.figure()
    plt.errorbar(xticks, actual, yerr = errorBar, fmt = 'o')
    plt.scatter(xticks, model, c = 'g')
    plt.xlim(0, 22000)
    plt.annotate(xy = (0, 0), fontsize = 8, s = obi + wan + kenobi)
    plt.savefig("images/" + wdStar.name + "_flux.png")
    plt.close()

def massRadiusTable():
    """ Grabs the necessary data for the getMassFromRT() function
    Requires the bergeron mass radius files to be under the BergeronFiles directory in either a Thin or Thick folder
    """
    fileNames = [f for f in os.listdir("BergeronFiles/Thick")]
    mrBerg = []
    massOrder = []
    for name in fileNames:
        f = open("BergeronFiles/Thick/" + name, 'r')
        tempMeanwhile = []
        for l in f:
            tempMeanwhile.append(list(map(float, l.split())))
        useTheForce = [tempMeanwhile[3*a] + tempMeanwhile[3*a + 1] + tempMeanwhile[3*a + 2] for a in range(int(len(tempMeanwhile)/3))]
        mrBerg.append(useTheForce)
        massOrder.append(int(name[3:6]))
        f.close()
    return mrBerg, massOrder

def getMassFromRT(radius, teff, mrBerg, massOrder):
    """ Returns a mass given a radius and effective temperature """
    endPoints = []
    for b in mrBerg:
        for a in range(len(b)):
            tempPoints = []
            if b[a][1] <= teff:
                tempPoints.append(b[a])
                tempPoints.append(b[a - 1])
                endPoints.append(tempPoints)
                break

    interpolatedPoints = []
    for c in range(len(endPoints)):
        a = endPoints[c]
        scaleFactor = (teff - a[0][1])/(a[1][1] - a[0][1])
        interpolatedPoints.append([massOrder[c], interpolateScale(a[0][2], a[1][2], scaleFactor)])

    radiiList = []

    for a in interpolatedPoints:
        radiiList.append([a[0], np.sqrt(universalG*(a[0]/100)*solarMassinGram/(10**(a[1])))])
    massInterpolation = []

    for r in range(len(radiiList)):
        radii = radiiList[r][1]
        if radius > radii:
            massInterpolation.append(radiiList[r])
            massInterpolation.append(radiiList[r - 1])
            break

    radiusScale = (radius - massInterpolation[0][1])/(massInterpolation[1][1] - massInterpolation[0][1])
    finalMass = interpolateScale(massInterpolation[0][0], massInterpolation[1][0], radiusScale)
    return finalMass*1e-2*solarMassinGram

def monteFit(logg, bergTable, teffR, loggR, against, err, labels):
    """ Fitting function for the Monte-Carlo calls
    Effectively the same, just less printing
    """

    captureList = []
    for a in labels:
        captureList.append(bergeronOrder.index(a))

    params = Parameters()
    params.add('teff', value = 1.569e4, min = 1.5e3, max = 1.2e5, vary = True)
    params.add('rd', value = 3e-12, min = 0, max = 1, vary = True)

    dagobah = Minimizer(residual, params, fcn_args = (teffR, loggR, logg, bergTable, captureList, against, err))
    monteMin = dagobah.leastsq(xtol = 3e-16, ftol = 3e-16)

    return [monteMin.params['teff'].value, monteMin.params['rd'].value]

def monteGame(bergTable, mrBerg, massOrder, teffR, loggR, against, err, labels, parallax):
    """ Same as startFits() but with no printing """

    logg = 8.0
    r2d2 = [0, 0, 0, 0]
    for a in range(5):
        temp, rd = monteFit(logg, bergTable, teffR, loggR, against, err, labels)
        wdDistance = (1/parallax)*pcInCm
        radius = rd*wdDistance
        smolMass = getMassFromRT(radius, temp, mrBerg, massOrder)
        logg = np.log10(universalG*smolMass/(radius**2))


        r2d2[0] = temp
        r2d2[1] = logg
        r2d2[2] = smolMass/solarMassinGram
        r2d2[3] = radius/solarRadiusinCm

    return r2d2

def monteCarloTime(wdStar):
    """ Does a Monte-Carlo run of n=1000 to create a corner plot
    The images are saved as 'images/' + WhiteDwarfName + '_corner.png'
    """

    observed = wdStar.observedFluxes
    errors = wdStar.fluxErr
    checkLog = wdStar.logg
    arcsec = wdStar.parallax
    arcErr = wdStar.paraerr


    f = open('bergeronFlux.csv', 'r')
    bergTable = []
    for a in f:
        bergTable.append(list(map(float, a.split(', '))))
    teffR = []
    loggR = [7.0, 7.5, 8.0, 8.5, 9.0, 9.5]
    for a in bergTable:
        if not (a[0] in teffR):
            teffR.append(a[0])

    mrBerg = massRadiusTable()[0]
    massOrder = massRadiusTable()[1]

    nowThisisPodRacing = 500

    darthPlagueis = []
    start = time.time()
    for n in range(nowThisisPodRacing):
        print("On Monte-Carlo run: %i" % n)
        randCoeff = np.random.standard_normal((1, len(observed)))
        print("got the coeffs")
        randPara = np.random.standard_normal()*arcErr + arcsec
        randFluxMeasures = np.array(observed) + np.array(errors)*randCoeff
        darthPlagueis.append(monteGame(bergTable, mrBerg, massOrder, teffR, loggR, randFluxMeasures, errors, wdStar.labels, arcsec) +  [1/randPara] + randFluxMeasures.tolist()[0])
    end = time.time()
    print("%i runs took %f seconds" % (nowThisisPodRacing, end - start))
    print("Building corner plot . . .")
    cornerData = darthPlagueis
    cornerLabel = ["Teff", "Log(g)", "Mass", "Radius", "Distance"] + wdStar.labels
    midichlorians = corner.corner(cornerData, labels = cornerLabel, show_titles = True, use_math_text = True, range = [1.0 for a in range(len(cornerLabel))])
    maceWindu = corner.quantile([a[3] for a in cornerData], [.25, .5, .75])
    print("Quantiles?: ")
    print(maceWindu)
    midichlorians.savefig('images/' + wdStar.name + '_corner.png')

def magnitudeToFlux(obsMag, obsErr, label):
    """ Converts observed magnitudes and errors into fluxes f_λ (erg cm^-2 s^-1 Å^-1)
    The label argument is a list of the filter bands used and is the same label in the WhiteDwarf object and file.
    Example label: ['J', 'B', 'v', 'z']
    """
    flux_zeroes = [bergZeroes[bergeronOrder.index(i)] for i in label]
    flux_mags = np.array([10**(obsMag[a]*-.4)*flux_zeroes[a] for a in range(len(obsMag))])
    flux_errs = [(10**(obsErr[a]*.4) - 1)*flux_mags[a] for a in range(len(obsMag))]
    distFlux = flux_mags
    stringo = ""
    for a in range(len(distFlux)):
        stringo += str(distFlux[a]) + ", " + str(flux_errs[a]) + ", "

    print(distFlux)
    print(flux_errs)
    print(stringo)

def calculateErrors(wdStar, rdVal, rdErr, tErr):
    distanceErr = (wdStar.paraerr)/(wdStar.parallax)*wdStar.distance
    radiusErr = np.sqrt((rdErr/rdVal)**2 + (distanceErr/wdStar.distance)**2)*wdStar.radius
    radiusEndPoints = [wdStar.radius - radiusErr, wdStar.radius + radiusErr]
    massPoints = [getMassFromRT(radiusEndPoints[0], wdStar.T_eff, massRadiusTable()[0], massRadiusTable()[1]), getMassFromRT(radiusEndPoints[1], wdStar.T_eff, massRadiusTable()[0], massRadiusTable()[1])]
    massErr = np.abs(massPoints[0] - massPoints[1])
    logErr = wdStar.logg*np.sqrt((massErr/wdStar.mass)**2 + 4*(radiusErr/wdStar.radius)**2)

    wdStar.set_Derr(distanceErr)
    wdStar.set_Rerr(radiusErr)
    wdStar.set_Merr(massErr)
    wdStar.set_Lgerr(logErr)
    wdStar.set_Terr(tErr)

    print("Here are the errors")
    print("Distance: %f | Radius: %f | Mass: %f | Log(g): %f | Temperature: %f" % (wdStar.disterr/wdStar.distance, wdStar.Rerr/wdStar.radius, wdStar.Merr/wdStar.mass, wdStar.Lgerr/wdStar.logg, wdStar.T_err/wdStar.T_eff))

def mainSequence():
    """ Main start function.
    Grabs the white dwarfs from the file 'whiteDwarfList', calculates the fit and then creates the flux plot for each one
    """
    print("how dhi")
    readFluxFrom("whiteDwarfList")
    for a in whiteDwarfs:
        startFits(a)
        fluxGraph(a)

def hayashiTrack():
    """ Monte-Carlo start function
    Grabs the white dwarfs from the file 'whiteDwarfList' and then proceeds with Monte-Carlo on each one
    """
    readFluxFrom("whiteDwarfList")
    for a in whiteDwarfs:
        monteCarloTime(a)




#magnitudeToFlux([17.5], [.02], ["g"])

#print(getMassFromRT(758963351.3, 4860, massRadiusTable()[0], massRadiusTable()[1]))

# WD 0011-721, .05277, .00110, 3.25267377581e-15, 6.04717203796e-17, 2.56324701534e-15, 4.76543199408e-17, 1.693675037e-15, 3.14877698504e-17, 6.42985902934e-16, 1.80140847477e-17, 2.9515233979e-16, 1.10765997373e-17, 1.14780805545e-16, 6.5215594766e-18, [V, R, I, J, H, K]
# WD 0326-273, .0413, .00131, 1.36848647093e-14, 6.44947804023e-16, 9.05295705839e-15, 4.26652721729e-16, 4.79545794325e-15, 2.26002970104e-16, 1.60030017007e-15, 1.54394073696e-16, 6.51697763217e-16, 5.63233855959e-17, 2.42030037644e-16, 2.82844162595e-17, [V, R, I, J, H, K]
# WD 0034-602, .04141, .00157, 1.38072971258e-14, 6.50717866045e-16, 8.87647446892e-15, 4.18335353531e-16, 4.79502041295e-15, 2.25982349937e-16, 2.33792237043e-15, 1.10182886774e-16, 5.54883040749e-16, 2.08238814836e-17, 1.73000050696e-16, 9.8294319744e-18, 6.54439144653e-17, 5.65603111961e-18, [B, V, R, I, J, H, K]


mainSequence()

#hayashiTrack()
