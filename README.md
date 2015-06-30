# UCSFsaxs
GUI for SAXS analysis, with smearing correction for line-collimated beam.

This Python program uses PyQt4, Numpy, Scipy, Matplotlib and Pyqtgraph (which is included here in a folder).
It also uses a routine written in Fortran90 to speed up the smearing matrix calculation. You need to compile the trans_smear.f90 using f2py.
Install f2py (which should be included in Numpy)

>>> f2py -c trans_smear.f90 -m trans_smear

And you should have 'trans_smear.so' to import into Python.

Then, run the program via python:
>>> python main_saxsgui.py

The algorithm implemented here is from:
Rambo, Robert P., and John A. Tainer. "Accurate assessment of mass, models and resolution by small-angle scattering." Nature 496.7446 (2013): 477-481.
Hansen, Steen. "Bayesian estimation of hyperparameters for indirect Fourier transformation in small-angle scattering." Journal of applied crystallography 33.6 (2000): 1415-1421.


-Daniel Elnatan (Agard Lab @ UCSF)
