import subprocess
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import ICRS
import math


sed_freq = []
sed_flux = []
sed_eflux = []
sed_filter = []

#   14:03:57.185 -15:01:10.52
def grabFromVizier():
    coords = input("Please put the desired RA DEC: ")
    f = open("VizierData/" + coords.replace(' ', '').replace('.', '').replace(':', '_') + ".vot", 'w+')
    subprocess.run(["vizquery", "-mime=vot", "-c=" + coords, "-phot", "-c.rs=5", "-out=_sed4"], stdout = f)
    f.close()

def parseTable(filename):
    t = Table.read(filename, format='votable')
    sed_freq = t['sed_freq']
    sed_flux = t['sed_flux']
    sed_eflux = t['sed_eflux']
    sed_filter = t['sed_filter']
    return t

def gaiaToJ2000(ra, dec, pmra, pmdec, ref):
    years = ref - 2000
    J2000 = [ra - pmra/3.6e6*year, dec-pmdec/3.6e6*year]
    return J2000



#grabFromVizier()
#print(parseTable('EC14012.vot'))
