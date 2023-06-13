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
  https://github.com/gbernstein/gbutil/

## Installation

* First, install the three external prereqs as per their instructions,
  being sure to use the compiler you will use for this code.  I found
  the installation for CSPICE to be nonstandard and clunky.
  Before running the build script I had to set up the CXX and CC
  environment variables and `setenv TKCOMPILER gcc-6` to use my gcc-6
  compiler too.  For the Mac, I had to alter several scripts replace ar, ranlib calls with
  `libtool -static -o $LIBRARY.a *.o`.  Then I needed to change the
  output library name to `libcspice.a` to meet the usual conventions so
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

Apologies that documentation is currently minimal.  Most of the 
executables will give you brief help if called with no arguments or
with `-h` as an argument. 

Once you've built the
executables, some simple tests/demos can be done as follows:

`% bin/MpcFit junk < testdata/pluto14yr.ast > junk.fit`

will run the `MpcFit` program on the fictitious observations of Pluto
stored in the specified `.ast` file. The stderr output will be various
diagnostics of the fit, including the residual errors.  Two new files
named `junk.aei` and `junk.abg` will be created with ASCII renditions'
of the results of the fit, in element space and in the BK 
sort-of-cartesian system.  The elements in `junk.aei` should make
sense for Pluto and be close (but not exactly match) the elements in the 
saved file `testdata/pluto14yr.aei`.  The `testdata/pluto14yr.ast` is
an example of the MPC-like format that can be used as input to this 
program to fit a single object's observations.  There are also 
some example files for 1992 QB1 = Albion and for some fictitious
Planet 9 observations.

With the `junk.fit` program that results from std output of 
`MpcFit`, one can predict future positions (or past) with 
a command such as

`% bin/PrintEphemeris -startdate 2022.0 -enddate 2023.0 -stepdate 30 < junk.fit`

The various `bin/Bulk*` executables are written to process larger
amounts of data for multiple objects, using input and output in
FITS tables.  Their help messages give some details about formats.

### Some gotchas

The code only knows about observatories whose MPC codes appear
in the file `data/observatories.dat`.  You'll need to add the
relevant info about new sites to that file.

The CFITSIO software is a bit picky about data types, so beware of
these details:  in the input FITS files to the _Bulk*cpp_ programs,
the `EXPNUM` column should 4-byte integers (FITS type
"J") and the `ORBITID` and `OBJID` columns should be 8-byte (FITS type "K"). The
`MJD` and other coordinate columns should be double-precision ("D").



