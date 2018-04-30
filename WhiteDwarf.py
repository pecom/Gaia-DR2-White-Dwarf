class WhiteDwarf:
    """An object that holds all the needed information on a White Dwarf for aNewHope.py"""
    #TEFF/RADIUS/MASS/LOGG are all in CGS

    def __init__(self):
        self.name = ""
        self.observedMagnitudes = []
        self.magnitudeErr = []
        self.observedFluxes = []
        self.fluxErr = []
        self.labels = []
        self.T_eff = 0
        self.T_err = 0
        self.radius = 0
        self.Rerr = 0
        self.mass = 0
        self.Merr = 0
        self.logg = 0
        self.Lgerr = 0
        self.parallax = 0
        self.paraerr = 0
        self.distance = 0
        self.disterr = 0
        self.chisqr = 0

    def set_name(self, name):
        self.name = name

    def set_fluxes(self, fluxes, err):
        self.observedFluxes = fluxes
        self.fluxErr = err

    def set_magnitudes(self, magni, err):
        self.observedMagnitudes = magni
        self.magnitudeErr = err

    def set_Teff(self, teff):
        self.T_eff = teff

    def set_Terr(self, err):
        self.T_err = err

    def set_radius(self, radius):
        self.radius = radius

    def set_Rerr(self, err):
        self.Rerr = err

    def set_mass(self, mass):
        self.mass = mass

    def set_Merr(self, err):
        self.Merr = err

    def set_logg(self, logg):
        self.logg = logg

    def set_Lgerr(self, err):
        self.Lgerr = err

    def set_order(self, label):
        self.labels = label

    def set_para(self, arcsec, error):
        self.parallax = arcsec
        self.distance = 1/arcsec*3.086e18
        self.paraerr = error

    def set_Derr(self, err):
        self.disterr = err

    def set_chi(self, chi):
        self.chisqr = chi
