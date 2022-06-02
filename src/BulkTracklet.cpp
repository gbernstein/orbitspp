// Run orbit fits given sets of observations stored in FITS tables
#include "Pset.h"
#include "Fitter.h"
#include "Elements.h"
#include "StringStuff.h"
#include "FTable.h"
#include "FitsTable.h"
#include "Astrometry.h"
#include "Exposures.h"

#include <iostream>

const string usage =
  "BulkFit: produce orbit fits for sets of DES observations.  The input data are\n"
  "linked sets of detections, identified by common ID number.  Outputs are \n"
  "fit results for each orbit ID.\n"
  "Usage:\n"
  "  BulkFit [parameter file...] [-<key> <value>...]\n"
  "  where any parameter file(s) given will be scanned first, then parameter key/value\n"
  "  pairs on cmd line will be read and override file values.\n"
  "  Program options are listed below.\n"
  "\n"
  "Input detections can be provided either in a binary FITS table or as\n"
  "ASCII file at stdin. Detections for given orbit should be contiguous\n"
  "in the input.  FITS file should have columns for TRACK_ID, EXPNUM (=DECam exposure\n"
  "number, in which case exposure file must be provided) or MJD, then for RA, DEC\n"
  "(in degrees), and SIGMA (uncertainty per coordinate in arcsec). Header must contain\n"
  "entries for RA0, DEC0, MJD0, X0, Y0, Z0 specifying reference frame orientation,\n"
  "reference time, and origin.  All detections should be in vicinity of RA0,DEC0\n"
  "\n"
  "If no FITS input file is given, data will be read from stdin.  First non-comment\n"
  "line should give the reference frame info.  Successive lines have format\n"
  "  <orbit id> <expnum or mjd>  <ra>  <dec>  <sigma>\n"
  "\n"
  "When expnum is used, atmospheric turbulence contribution is added to measurement\n"
  "errors if it is known.\n"
  "\n"
  "Output FITS binary table will the reference frame info in header, and contain\n"
  "columns for \n"
  " orbit ID,\n"
  " a FLAG that is set for orbit-fitting failures,\n"
  " the 6-element ABG array specifying best-fit orbit,\n"
  " the ABGCOV 6x6 covariance matrix, flattened.\n"
  "If no FITS output file is given, orbit fit results will be sent to stdout.\n"
  "\n"
  "If the 'residuals' parameter is set to true/T/1, then the residuals of the data\n"
  "to the best-fit orbit will be printed to stdout (including for aborted fits).";

using namespace std;
using namespace orbits;

// Some possible fitting failure modes
const double MAX_LINEAR_CHISQ_PER_PT = 5e5; //**100.;
const int LINEAR_CHISQ_TOO_HIGH = 1;
const double MAX_LINEAR_ORBIT_KE = 10.;  // Max ratio of |KE/PE| before quitting
const int LINEAR_KE_TOO_HIGH = 2;
const double MAX_LINEAR_GAMMA = 0.5;
const int LINEAR_GAMMA_TOO_HIGH = 4;
const int NONCONVERGENCE = 8;
const int INSUFFICIENT_OBSERVATIONS = 16;

int main(int argc,
	 char *argv[])
{
  try {
    string ephemerisPath;
    string exposurePath;
    string observationPath;
    string orbitPath;
    double bindingFactor;
    double gamma0;
    double dGamma;
    int obscode;
    bool showResiduals;
    Pset parameters;
   
    {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMember("observationFile",&observationPath, def,
			   "FITS file holding observation list (null=>ASCII/stdin)", "");
      parameters.addMember("orbitFile",&orbitPath, def,
			   "FITS file for output orbit fits (null=>ASCII/stdout)", "");
      parameters.addMember("ephemerisFile",&ephemerisPath, def,
			   "SPICE file (null=>environment", "");
      parameters.addMember("exposureFile",&exposurePath, def,
			   "FITS file holding DECam exposure info", "");
      parameters.addMember("bindingFactor",&bindingFactor, def | low,
			   "Chisq penalty at unbinding", 4., 0.);
      parameters.addMember("gamma0",&gamma0, def,
			   "Center of gamma prior (<=0 for none)", 0.);
      parameters.addMember("dGamma",&dGamma, def | lowopen,
			   "Width (sigma) of gamma prior", 0.01, 0.);
      parameters.addMember("obscode",&obscode, def | low,
			   "Observatory code (default: CTIO)", 807, 0);
      parameters.addMember("residuals",&showResiduals, def,
			   "Print fit residuals to stdout? (false)", false);
    }
    parameters.setDefault();

    if (argc<2 || string(argv[1])=="-h" || string(argv[1])=="--help") {
      cout << usage << endl;
      parameters.dump(cerr);
      exit(1);
    }
    
    {
      // Read any parameter files
      int nPositional=0;
      for (int iarg=1; iarg < argc && argv[iarg][0]!='-'; iarg++) {
	ifstream ifs(argv[iarg]);
	if (!ifs) {
	  cerr << "Can't open parameter file " << argv[iarg] << endl;
	  cerr << usage << endl;
	  exit(1);
	}
	try {
	  parameters.setStream(ifs);
	} catch (std::runtime_error &m) {
	  cerr << "In file " << argv[iarg] << ":" << endl;
	  quit(m,1);
	}
	nPositional++;
      }
      // And read arguments from the remaining command line entries
      parameters.setFromArguments(argc, argv);
    }

    // Read the ephemeris
    Ephemeris eph(ephemerisPath);

    // Open the input frame if we have it, and set reference Frame.
    bool observationsFromFile = !observationPath.empty();
    bool useExpnum; // Set true if input will have expnum instead of MJD.
    Frame frame;
    vector<LONGLONG> idIn;
    vector<int> expnum1;
    vector<int> expnum2;
    vector<double> mjd1;
    vector<double> mjd2;
    vector<double> ra1;
    vector<double> ra2;
    vector<double> dec1;
    vector<double> dec2;
    vector<Matrix44> sigma;
    
    if (observationsFromFile) {
      FITS::FitsTable ft(observationPath,FITS::ReadOnly, 1);
      img::FTable obsTable = ft.extract();

      // Get reference frame info from header
      double ra0, dec0, mjd0, x0, y0, z0;
      obsTable.header()->getValue("RA0",ra0);
      obsTable.header()->getValue("DEC0",dec0);
      obsTable.header()->getValue("MJD0",mjd0);
      astrometry::SphericalICRS pole(ra0*DEGREE, dec0*DEGREE);
      astrometry::Orientation orient(pole);
      orient.alignToEcliptic();  // We'll do this by default.
      if ((obsTable.header()->getValue("X0",x0)
	   && obsTable.header()->getValue("Y0",y0)
	   && obsTable.header()->getValue("Z0",z0))) {
	frame = Frame(astrometry::CartesianICRS(x0,y0,z0),
		      orient, eph.mjd2tdb(mjd0));
      } else {
	// Put origin at position of observatory at MJD0 if not given
	frame.tdb0 = eph.mjd2tdb(mjd0);
	frame = Frame(eph.observatory(obscode, frame.tdb0),
		      orient, frame.tdb0);
      }

      // Extract data from the table
      obsTable.readCells(idIn, "ORBITID");
      // Are we using expnum or MJD?
      auto colnames = obsTable.listColumns();
      useExpnum = (std::find(colnames.begin(),colnames.end(),"EXPNUM1")!=colnames.end());
      if (useExpnum) {
	obsTable.readCells(expnum1, "EXPNUM1");
	obsTable.readCells(expnum2, "EXPNUM2");

	// Read the exposure table later.
      } else {
	// We'll use MJD's
	obsTable.readCells(mjd1,"MJD1");
	obsTable.readCells(mjd2,"MJD2");
	
      }

      obsTable.readCells(ra1, "RA1");
      obsTable.readCells(ra2, "RA2");
      obsTable.readCells(dec1,"DEC1");
      obsTable.readCells(dec2,"DEC2");

      obsTable.readCells(sigma,"SIGMA"); 

    } else {

      // Info will come from stdin.  Read the reference frame from first
      // non-comment line.
      string buffer;
      stringstuff::getlineNoComment(cin, buffer);
      std::istringstream iss(buffer);
      double ra0, dec0, mjd0;
      iss >> ra0 >> dec0 >> mjd0;
      frame.orient.setPole(astrometry::SphericalICRS(ra0*DEGREE,dec0*DEGREE));
      // Align the frame to ecliptic
      frame.orient.alignToEcliptic();
      frame.tdb0 = eph.mjd2tdb(mjd0);
      // Put origin at position of observatory at MJD0
      frame.origin = eph.observatory(obscode, frame.tdb0);
      // Expnum vs MJD will be determined line by line.
      useExpnum = false;
    }

    // Create vectors to hold output if we are writing FITS
    vector<LONGLONG> idOut;
    vector<int> fitFlags;  // Record info on fitting
    vector<double> chisq;
    vector<int> dof;
    vector<vector<double>> abg;
    vector<vector<double>> abgCov;
    vector<vector<double>> abgInvCov;

    // Are we writing to a FITS table or stdout?
    bool orbitsToFile = !orbitPath.empty();

    // If writing to stdout, print out the reference frame
    if (!orbitsToFile) {
      Frame::writeHeader(cout);
      frame.write(cout);
    }
    
    // If we read the exposure table it'll go here:
    ExposureTable* exposureTable = nullptr;
    
    // And we will build this index from expnum to table row number:
    std::map<int,int> expnumIndex;
      
    // Keep next line of input data, either from stdin or next row of table
    string inputLine;
    long inputRow = 0;
    bool done = false;  // set flag when no more input observations
    if (observationsFromFile) {
      done = inputRow==idIn.size();
    } else {
      if (!stringstuff::getlineNoComment(cin, inputLine))
	done = true;
    }

    while (!done) {
      // Set up a new fitter
      Fitter fit(eph, Gravity::GIANTS);
      if (bindingFactor > 0.)
	fit.setBindingConstraint(bindingFactor);
      if (gamma0 > 0.)
	fit.setGammaConstraint(gamma0, dGamma);
      
      // Read input observations until no more match orbit_id of first
      bool isFirst = true;
      LONGLONG idThis;
      while (true) {
	double raThis1, raThis2, decThis1, decThis2, sigmaThis;
	int expnumThis1, expnumThis2;
	double mjdThis1, mjdThis2;
	if (observationsFromFile) {
	  // Read new data from table
	  if (isFirst) {
	    idThis = idIn[inputRow];
	    isFirst = false;
	  } else {
	    if (idIn[inputRow]!=idThis)
	      break; // Done getting data for this orbit
	  }
	  raThis1 = ra1[inputRow];
	  raThis2 = ra2[inputRow];
	  decThis1 = dec1[inputRow];
	  decThis2 = dec2[inputRow];
	  sigmaThis = sigma[inputRow];
	  if (useExpnum) {
	    expnumThis1 = expnum1[inputRow];
	    expnumThis2 = expnum2[inputRow];
	  } else {
	    mjdThis1 = mjd1[inputRow];
	    mjdThis2 = mjd2[inputRow];
	  }
	  ++inputRow;
	  if (inputRow >= idIn.size())
	    done = true; // This is last input observation.
	} else {
	  // Read data from input line.
	  std::istringstream iss(inputLine);
	  LONGLONG idTest;
	  iss >> idTest >> mjdThis >> raThis >> decThis >> sigmaThis;
	  if (isFirst) {
	    idThis = idTest;
	    isFirst = false;
	  } else {
	    if (idTest!=idThis)
	      break; // Done getting data for this orbit
	  }
	  // Is input expnum or MJD?
	  useExpnum = mjdThis1 > 100000.; // Expnums are >100,000, MJD's are ~50,000
	  if (useExpnum) {
	    expnumThis1 = static_cast<int> (round(mjdThis1));
	    expnumThis2 = static_cast<int> (round(mjdThis2));
	  }
	  // Read next input line (if any)
	  done = !getlineNoComment(cin, inputLine);
	}

	// Build Observation
	Tracklet track;
	track.radec1 = astrometry::SphericalICRS(raThis1*DEGREE, decThis1*DEGREE);
	track.radec2 = astrometry::SphericalICRS(raThis2*DEGREE, decThis2*DEGREE);
	
	track.cov = sigmaThis * ARCSEC * ARCSEC;
      
	if (useExpnum) {
	  // Read exposure table if needed and not here
	  if (!exposureTable) {
	    exposureTable = new ExposureTable(exposurePath);
	  }

	  double mjdThis1, mjdThis2;
	  astrometry::CartesianICRS xyzThis1, xyzThis2;
	  astrometry::SphericalICRS radecThis1, radecThis2;
	  // Get xE and mjd from the DECam table
	  if (!exposureTable->observingInfo(expnumThis1, mjdThis1, xyzThis1, radecThis1)) {
	    /**/cerr << "Missing exposure info for " << expnumThis1 << endl;
	    continue;
	  }
	  if (!exposureTable->observingInfo(expnumThis2, mjdThis2, xyzThis2, radecThis2)) {
	    /**/cerr << "Missing exposure info for " << expnumThis2 << endl;
	    continue;
	  }

	    
	  for (int i=0; i<3; i++){
	  	track.observer1[i] = xyzThis1[i];
	  	track.observer2[i] = xyzThis2[i];
	  }

	  track.tdb1 = eph.mjd2tdb(mjdThis1);
		track.tdb2 = eph.mjd2tdb(mjdThis2);

	  if (exposureTable->isAstrometric(expnumThis1)) {
	  	Matrix44 atm;
	  	for (int i = 0; i < 4, i++){
	  		for (int j = 0; i < 4, i++){
	  			atm(i,j) = 0;
	  		}
	  	}
	  	atm(0,0) = exposureTable->atmosphereCov(expnumThis1)(0,0);
	  	atm(1,1) = exposureTable->atmosphereCov(expnumThis1)(1,1);
	  	atm(0,1) = atm(1,0) = exposureTable->atmosphereCov(expnumThis1)(1,0);
	  	atm(2,2) = exposureTable->atmosphereCov(expnumThis2)(0,0);
	  	atm(3,3) = exposureTable->atmosphereCov(expnumThis2)(1,1);
	  	atm(2,3) = atm(3,2) = exposureTable->atmosphereCov(expnumThis2)(1,0);

	    track.cov += atm;
	  }
	} else {
	  // Input data is an MJD.  Use ephemeris to get observatory posn.
	  track.tdb1 = eph.mjd2tdb(mjdThis1);
	  track.tdb2 = eph.mjd2tdb(mjdThis2);
	  track.observer1 = eph.observatory(obscode, track.tdb1);
	  track.observer2 = eph.observatory(obscode, track.tdb2);
	  
	  // Note no atmospheric term added to covariance here.
	}

	fit.addObservation(track);
      } // end collecting observations for this orbit


      // Set up Fitter and do fit
      int errorCode = 0;
      if (fit.nObservations() < 3) {
	// Do not attempt to fit this orbit
	errorCode = INSUFFICIENT_OBSERVATIONS;
	if (orbitsToFile) {
	  idOut.push_back(idThis);
	  fitFlags.push_back(errorCode);
	  chisq.push_back(0.);
	  dof.push_back(0);
	  vector<double> zero(6,0.);
	  abg.push_back(zero);
	  vector<double> zzero(36,0.);
	  abgCov.push_back(zzero);
	  abgInvCov.push_back(zzero);
	} else {
	  cout << "# ------------------------" << endl;
	  cout << "OrbitID " << idThis
	       << " ERROR " << errorCode
	       << endl;
	}
	continue;
      }
      fit.setFrame(frame);
      try {
	fit.setLinearOrbit();
	fit.setLinearOrbit(); // Another iteration
	auto abg = fit.getABG(true);
	if (fit.getChisq() / (2*fit.nObservations()) > MAX_LINEAR_CHISQ_PER_PT) {
	  errorCode = LINEAR_CHISQ_TOO_HIGH;
	} else if (abg[ABG::G] > MAX_LINEAR_GAMMA) {
	  errorCode = LINEAR_GAMMA_TOO_HIGH;
	} else if ( (abg[ABG::ADOT]*abg[ABG::ADOT] + abg[ABG::BDOT]*abg[ABG::BDOT])
		    / (2*GM*pow(abs(abg[ABG::G]),3.)) > MAX_LINEAR_ORBIT_KE) {
	  errorCode = LINEAR_KE_TOO_HIGH;
	} else {
	  fit.newtonFit();
	}
      } catch (std::runtime_error& e) {
	errorCode = NONCONVERGENCE;
      }

      // Print fitting residuals if requested
      if (showResiduals) {
	cout << "# Fitting results for " << idThis << endl;
	fit.printResiduals(cout);
      }
      
      // Output results
      if (orbitsToFile) {
	// save to output vectors
	idOut.push_back(idThis);
	fitFlags.push_back(errorCode);
	// No useful chisq if no convergence:
	chisq.push_back( errorCode == NONCONVERGENCE ? 0. : fit.getChisq()); 
	dof.push_back(fit.getDOF());
	if (errorCode>0) {
	  // Failed fit just gets empty results
	  vector<double> zero(6,0.);
	  abg.push_back(zero);
	  vector<double> zzero(36,0.);
	  abgCov.push_back(zzero);
	  abgInvCov.push_back(zzero);
	} else {
	  auto a = fit.getABG();
	  vector<double> v(6);
	  for (int i=0; i<6; i++) v[i] = a[i];
	  abg.push_back(v);
	  vector<double> vv(36);
	  ABGCovariance c = fit.getInvCovarABG().inverse();
	  for (int i=0; i<6; i++)
	    for (int j=0; j<6; j++) vv[6*i+j] = c(i,j);
	  abgCov.push_back(vv);
	  
	  vector<double> vinv(36);

		ABGCovariance cinv = fit.getInvCovarABG();
		for (int i=0; i<6; i++)
			for (int j=0; j<6; j++) vinv[6*i+j] = cinv(i,j);
		abgInvCov.push_back(vinv);
	}

      } else {
	// Print results to stdout
	cout << "# ------------------------" << endl;
	if (errorCode > 0) {
	  cout << "OrbitID " << idThis
	       << " ERROR " << errorCode;
	  if (errorCode!=NONCONVERGENCE)
	    cout << " chisq " << fit.getChisq()
		 << " DOF " << fit.getDOF();
	  cout << endl;
	} else {
	  auto a = fit.getABG();
	  cout << "OrbitID " << idThis
	       << " chisq " << fit.getChisq()
	       << " DOF " << fit.getDOF()
	       << endl; /*
	       << " KE/|PE| " << (a[ABG::ADOT]*a[ABG::ADOT] +
				  a[ABG::BDOT]*a[ABG::BDOT] +
				  a[ABG::GDOT]*a[ABG::GDOT]) / (GM * pow(a[ABG::G],3))
				  << endl;*/
	  cout << a << endl;
	  ABGCovariance abgCov = fit.getInvCovarABG().inverse();
	  abgCov.write(cout);
	}
      }

    } // End orbit loop
    
    if (orbitsToFile) {
      // Write the output table
      FITS::FitsTable ft(orbitPath, FITS::Create + FITS::OverwriteFile);
      img::FTable out = ft.use();
      // Write reference frame information into header
      double ra0, dec0;
      frame.orient.getPole().getLonLat(ra0,dec0);
      out.header()->replace("RA0",ra0/DEGREE);
      out.header()->replace("DEC0",dec0/DEGREE);
      out.header()->replace("TDB0",frame.tdb0);
      out.header()->replace("PA0",frame.orient.getPA()/DEGREE);
      out.header()->replace("X0",frame.origin[0]);
      out.header()->replace("Y0",frame.origin[1]);
      out.header()->replace("Z0",frame.origin[2]);
      
      out.addColumn(idOut,"ORBITID");
      out.addColumn(fitFlags,"FLAGS");
      out.addColumn(chisq,"CHISQ");
      out.addColumn(dof,"DOF");
      out.addColumn(abg,"ABG",6);
      out.addColumn(abgCov,"ABGCOV",36);
      out.addColumn(abgInvCov,"ABGINVCOV",36);
    }
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
