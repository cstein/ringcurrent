# ringcurrent
A python implementation of [*Definitive Benchmark Study of Ring Current Effects on Amide Proton Chemical Shifts*](http://dx.doi.org/10.1021/ct2002607) by [Anders S. Christensen](https://github.com/andersx), Stephan Sauer and [Jan H Jensen](https://github.com/jhjensen2)

The implementation uses the point-dipole approximation to evaluate how aromatic rings changes the chemical shift of the backbone amide proton.

It is possible to search for other atoms (based on the atom name) and as such use the same model for alpha carbon hydrogens for instance.

## Requirements
  * [Open Babel](http://openbabel.org) with [python bindings](http://openbabel.org/docs/current/UseTheLibrary/PythonInstall.html)
  * Numpy

## Examples
You can test it out by running

    python rc.py examples/gb3_00000.pdb
