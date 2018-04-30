# Gaia-DR2-White-Dwarf

## BergeronSolver.py
This script is supposed to solve for some of the basic properties of a white dwarf using Bergeron's models. The output is either the solved white dwarf with: distance, radius, mass, Log(g), and temperature (mainSequence()) or a monte-carlo run of n=1000 with corner plots generated (hayashiTrack())

### Requirements
1. bergeronFlux.csv (included)
2. Bergeron's Thick and Thin cooling models to be located under a BergeronFiles directory.
3. An images directory to hold the graphs (included)
4. numpy, scipy, WhiteDwarf (included), matplotlib, lmfit, and corner
5. whiteDwarfList (example included)


## grabPhotometry.py
Simply grabs photometry from Vizier given an RA and DEC. The output is a .vot file with the photometry

### Requirements
1. The vizquery tool needs to be installed on the machine
2. astropy
3. A VizierData directory to hold the .vot files
