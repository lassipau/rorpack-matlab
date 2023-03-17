# Changelog for RORPack for Matlab

The main changes to "RORPack for Matlab" are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)

## [Unreleased]

### Added 

-

### Changed

- A new stabilization option 'full_K' in ObserverBasedRC
- Both DualObserverBasedRC and ObserverBasedRC can now be used also for systems with more inputs than outputs, i.e., dimU > dimY.

### Removed

- 

### Fixed

- ObserverBasedROMRC: A mistake in the stabilization parameter dimensions.
- Fixed typos in the example files using LowGainRC (added missing complex unit to computation of Pvals)

### Bugs found

- The approximate transfer function values used by LowGainRC in some of the examples are computed at P(w_k) instead of at P(iw_k)
- Computation of the initial state is not incorrectly in ConstrHeat2DCase1, ConstrHeat2DCase2, and ConstrHeat2DCase3.

## [v1.0.0] - 2021-11-08

### Added 

- Animation of the solution of the controlled PDE for the Timoshenko beam example case.

### Fixed

- Heat2DCase2 animation works for M != N now.
- Small fix to ObserverBasedROMRC.

## [v0.9.0] - 2021-07-07

### Added 

- RORPack for Matlab. RORPack is a software package for controller design and simulation of robust output tracking and disturbance rejection for linear partial diﬀerential equation (PDE) models.
- Examples that make use of RORPack and illustrate its functionality. These examples include 1-dimensional heat, wave and beam equations as well as 2-dimensional heat and wave equations.
- Documentation that gives some background on the theory that RORPack is used to approximate and simulate. Also includes descriptions of the typical structure of an example case, of the implemented controllers and of all the example cases.


[unreleased]: https://github.com/lassipau/rorpack-matlab/tree/dev
[v1.0.0]: https://github.com/lassipau/rorpack-matlab/releases/tag/v1.0.0
[v0.9.0]: https://github.com/lassipau/rorpack-matlab/releases/tag/v0.9.0

