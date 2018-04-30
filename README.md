# Gaia-DR2-White-Dwarf

#BergeronSolver.py
This script is supposed to solve for some of the basic properties of a white dwarf using Bergeron's models
The files required are: bergeronFlux.csv (included) and Bergeron's Thick and Thin cooling models to be located under a BergeronFiles directory
The modules required are: numpy, scipy, WhiteDwarf (included), matplotlib, lmfit, and corner
The output is either the solved white dwarf with: distance, radius, mass, Log(g), and temperature

#grabPhotometry.py
This file requires the vizquery tool to be installed on the machine
The modules requried are: astropy
The output is a .vot file with the photometry
