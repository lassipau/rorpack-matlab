# The RORPack (for Matlab) project

The Robust Output Regulation Package is an open-source Matlab package for controller design and simulation for robut output tracking and disturbance rejection for linear partial differential equations. The package includes a set of examples (folder 'examples/') on simulation of different types of controlled PDE systems on 1D and 2D spatial domains.

## Requirements

The package works on Matlab (R2020a and later, though earlier versions may be fine too). The controller design routines ObserverBasedROMRC and DualObserverBasedROMRC require the Control System Toolbox, and a couple of the examples make use of the (free) 3rd party package Chebfun (www.chebfun.org).

## Installation 

 1. Copy the package contents to a folder on your hard drive
 2. Add the folder "RORPack" and its subfolders ("controllers", "simulation", "visualization") to the Matlab PATH.
 5. Test to see that everything works by trying an example, e.g. running the file 'examples/Wave1DMain.m'


## Documentation

The mathematical documentation and general introduction to the package is contained in the folder ./"Introduction to RORPack". Additional documentation and lists of parameters of the routines are included in the beginnings of the code files and can be shown with the commands of the form 'help ObserverBasedRC'. The use of the routines for controller design, closed-loop simulations and visualization of the results are also illustrated in the included example files in the folder "examples/".


## Disclaimer

The purpose of RORPack is to serve as a tool to illustrate the theory of robust output regulation for distributed parameter systems and it should not (yet) be considered as a serious controller design software. The developers of the software do not take any responsibility and are not liable for any damage caused through use of this software. In addition, at its present stage, the library is not optimized in terms of numerical accuracy. Any helpful comments on potential improvements of the numerical aspects will be greatly appreciated!


## License

See LICENSE.txt for licensing information of the RORPack library.


