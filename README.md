# Delta-Sigma
A stand alone version of the Delta Sigma code currently used in CosmoSIS that doesn't require the use of the CosmoSIS data block.

# Compilation
To run with a python interface, from Delta-Sigma/ run:

$make SHARED=yes

# Running
For an exmaple of how to run, you can use:

$python pymain.py

If you want to change cosmologies then you have to supply a different power spectrum and change the paths around in pymain.py