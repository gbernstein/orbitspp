# orbitspp
Update of the Bernstein & Khushalani orbit-fitting software for outer solar system bodies.  Improvements include:
* All object-oriented code in C++
* Improve algorithms to avoid wayward solutions for very distant objects and high-e orbits.
* Use of multiprocessing for large-scale use [not yet]
* Switch to currently supported libraries:
** NASA/JPL `cspice` code for planetary ephemerides
** `Eigen` for linear algebra

## Prerequisites

* A C++-11-compliant C++ compiler, e.g. gcc version 5 or higher.
* NASA/JPL `cspice` code for planetary ephemerides:
https://naif.jpl.nasa.gov/naif/toolkit.html
* `cfitsio` FITS interface, also from NASA
https://heasarc.gsfc.nasa.gov/fitsio/
* `Eigen` template-based C++ linear algebra library:
http://eigen.tuxfamily.org/index.php?title=Main_Page or mirrored on github
* Three other repos, `astrometry`, `gbutil` and
  `gbfits`. https://github.com/gbernstein/astrometry,
  https://github.com/gbernstein/gbfits,
  https://github.com/gbutil/

## Installation

* First, install the three external prereqs as per their instructions,
  being sure to use the compiler you will use for this code.  I found
  the installation for CSPICE to be nonstandard and klunky.
  Before running the build script I had to set up the CXX and CC
  environment variables and `setenv TKCOMPILER gcc-6` to use my gcc-6
  compiler too.  For the Mac, I had to alter several scripts replace ar, ranlib calls with
  `libtool -static -o $LIBRARY.a *.o`.  Then I needed to change the
  output library name to `libspice.a` to meet the usual conventions so
  the linker can find it.
* The makefiles for my repos take their cues from environment
variables.  The csh version for my setup looks like:
```
setenv G6_DIR ~/G6
# These variables should point to the places you've installed the
# other packages.  The .h files for each (except EIGEN, which is
# header-only) should be in /include subdirectories of each of these.
setenv EIGEN_DIR $G6_DIR/include
setenv SPICE_DIR $G6_DIR/cspice
setenv CFITSIO_DIR $G6_DIR/cfitsio
setenv GBUTIL_DIR ~/CODE/gbutil
setenv GBFITS_DIR ~/CODE/gbfits
setenv ASTROMETRY_DIR ~/CODE/astrometry

# Give your compilation commands here
setenv CXX "g++-6 -fopenmp -march=haswell"
setenv CC "gcc-6 -fopenmp -march=haswell"
setenv CXXFLAGS "-O"

# Use whatever command you need to point to these packages libraries
# for runtime linking, e.g. maybe LD_LIBRARY_PATH here.
setenv DYLD_LIBRARY_PATH $CFITSIO_DIR/lib:$SPICE_DIR/lib

# These variables are needed during execution of the program.
# The orbit-fitting code uses them as default paths to important
# files.
# Note that this repo's /data directory contains the needed files
# for our TNO tasks.
# This one gives the locations of CTIO and other observatories:
setenv ORBIT_OBSERVATORIES $HOME/CODE/orbitspp/data/observatories.dat
# This one points to the master data file for SPICE:
setenv SPICE_KERNEL $HOME/CODE/orbitspp/data/de430.tm
# This is data on all the DES exposures
setenv DES_EXPOSURE_TABLE $HOME/CODE/orbitspp/data/alldes.exposure.positions.fits
```

* With these set up, you should just be able to type `make` in the
  root directory of this repository.  It will automatically execute
  the makefiles in the other three repos. `make clean` will likewise
  clean up the other ones too.
* Executables will end up in the _orbitspp/bin_ directory.  Put that
in your path as needed.
* `make tests` will build some test programs and install them in
  _orbitspp/testbin_.  I haven't yet documented the proper usage of
  these tests though.
* The SPICE kernel file in data/de430.tm tells SPICE where to find
its data files.  The PATH_VALUES line should be changed to match the
directory where you keep such things, e.g. the /data subdir of this
repo. Either edit this file, or (preferentially, to keep it from
getting rewritten every time you pull the repo) make another version
with a new name and make sure that the SPICE_KERNEL environment
variable points to it.


## Usage

The CFITSIO software is a bit picky about data types, so beware of
these details:  in the input FITS files to the _Bulk*cpp_ programs,
the `ORBITID` and `EXPNUM` columns should 4-byte integers (FITS type
"J") and the `OBJID` column should be 8-byte (FITS type "K"). The
`MJD` and other coordinate columns should be double-precision ("D").




