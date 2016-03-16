# UCSFsaxs
GUI for SAXS analysis, with smearing correction for line-collimated beam.

This Python program uses PyQt4, Numpy, Scipy, Matplotlib and Pyqtgraph (which is included here in a folder. I added support to display error bars in log-scale because as of this time, Pyqtgraph doesn't support it).
It also uses a couple routines written in Fortran90 to speed up the smearing matrix calculation and automatic determination of most linear region for Guinier analysis. You need to compile 'trans_smear.f90' & 'autorg_de.f90' using f2py. If you used pip to install numpy you should have it.
<pre>
f2py -c trans_smear.f90 -m trans_smear
f2py -c autorg_de.f90 -m autorg_de
</pre>
And you should have 'trans_smear.so' and 'autorg_de.so' to import into Python.

Then, run the program via python:
<pre>
python main_saxsgui.py
</pre>

The algorithms implemented here are directly taken from:<br>
<p>For Molecular Weight estimation:<br>
Rambo, Robert P., and John A. Tainer. "Accurate assessment of mass, models and resolution by small-angle scattering." Nature 496.7446 (2013): 477-481. 

<p>For Bayesian IFT:<br>
Hansen, Steen. "Bayesian estimation of hyperparameters for indirect Fourier transformation in small-angle scattering." Journal of applied crystallography 33.6 (2000): 1415-1421.


-Daniel Elnatan (Agard Lab @ UCSF)
