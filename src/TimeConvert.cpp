#include <vector>
#include <iostream>

#include "Std.h"
#include "FitsTable.h"
#include "Ephemeris.h"
#include "Pset.h"

using namespace std;

string usage =
  "Convert time entries between formats / time scales.\n"
  "usage: TimeConvert [filename:column] [filename:column] [-i format] [-o format]\n"
  "       <filename:column> entries give source / destination as either '-' for\n"
  "                         stdin/stdout, or a combination of FITS table filename\n"
  "                         and column\n"
  "       [-i/o format]     give the format of input/output, among these:\n"
  "                         jd,mjd:   (Modified) Julian Date (UTC)\n"
  "                         tdb:      Years since J2000. (TDB)\n"
  "                         ymd:      YYYY MM DD.DDDDD (UTC)\n"
  "                         ymdhms:   YYYY MM DD HH:MM:SS.SSS (UTC)\n"
  "                         Defaults to tdb in, mjd out.";



enum Format {JD, MJD, TDB, YMD, YMDHMS};
Format getFormat(const string s) {
  // Make uppercase version of string
  string ss = "";
  for (auto& c : s) ss+= std::toupper((unsigned char) c);
  if (ss=="JD")
    return JD;
  else if (ss=="MJD")
    return MJD;
  else if (ss=="TDB")
    return TDB;
  else if (ss=="YMD")
    return YMD;
  else if (ss=="YMDHMS")
    return YMDHMS;
  else
    throw std::runtime_error("Unknown time format " + s);
}

bool isStringFormat(const Format f) {
  if (f==YMD || f==YMDHMS)
    return true;
  else
    return false;
}
  
int
main(int argc, char *argv[]) {

  try {
    Pset parameters;
    string inFormat;
    string outFormat;
    string ephemerisPath;

    {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMember("i",&inFormat, def,
			   "input time format", "tdb");
      parameters.addMember("o",&outFormat, def,
			   "output time format", "mjd");
      parameters.addMember("ephemerisFile",&ephemerisPath, def,
			   "SPICE file (null=>environment", "");
    }
    parameters.setDefault();

    if (argc>1 && (string(argv[1])=="-h" || string(argv[1])=="--help")) {
      cout << usage << endl;
      parameters.dump(cerr);
      exit(1);
    }
    
    int positionalArgs = parameters.setFromArguments(argc, argv);
    
    string inFilename;
    string outFilename;
    string inColumnName;
    string outColumnName;
    
    // Parse the input and output destinations and formats
    auto inF = getFormat(inFormat);
    auto outF = getFormat(outFormat);
    
    string arg = positionalArgs > 1 ? argv[1] : "-";
    if (arg=="-") {
      // Do nothing, inFilename stays null
    } else {
      auto fc = stringstuff::split(arg,':');
      if (fc.size()!=2) {
	throw std::runtime_error("Bad filename:column spec: " + arg);
      }
      inFilename = fc.front();
      fc.pop_front();
      inColumnName = fc.front();
      fc.pop_front();
    }
    arg = positionalArgs > 2 ? argv[2] : "-";
    if (arg=="-") {
      // Do nothing, outFilename stays null
    } else {
      auto fc = stringstuff::split(arg,':');
      if (fc.size()!=2) {
	throw std::runtime_error("Bad filename:column spec: " + arg);
      }
      outFilename = fc.front();
      fc.pop_front();
      outColumnName = fc.front();
      fc.pop_front();
    }
    
    orbits::Ephemeris ephem;
    // vectors of I/O values
    vector<double> vd;
    vector<string> vs;

    // Acquire inputs
    if (!inFilename.empty()) {
      // Read from a table
      FITS::FitsTable ff(inFilename,FITS::ReadOnly,1);
      auto table = ff.use();
      if (isStringFormat(inF)) {
	table.readCells(vs,inColumnName);
      } else {
	table.readCells(vd,inColumnName);
      }
    } else {
      // Read from stdin
      string buffer;
      while (stringstuff::getlineNoComment(cin,buffer)) {
	// Get just the first number of needed fields of the input line
	auto fields = stringstuff::split(buffer);
	string ss = fields.front();
	fields.pop_front();
	if (inF==YMD) {
	  if (fields.size()<2) {
	    throw std::runtime_error("Insufficient YMD info at " + buffer);
	  }
	  for (int i=0; i<2; i++) {
	    ss += " " + fields.front();
	    fields.pop_front();
	  }
	  vs.push_back(ss);
	} else if (inF==YMDHMS) {
	  if (fields.size()<3) {
	    throw std::runtime_error("Insufficient YMDHMS info at " + buffer);
	  }
	  for (int i=0; i<3; i++) {
	    ss += " " + fields.front();
	    fields.pop_front();
	  }
	  vs.push_back(ss);
	} else {
	  // Convert the first field into double
	  vd.push_back(strtod(ss.c_str(),nullptr));
	}
      } // End of input
    }

    // Prepare arrays we'll need
    int n = isStringFormat(inF) ? vs.size() : vd.size();
    vector<astrometry::UT> vut(n);
    if (isStringFormat(outF)) {
      if (vs.empty())
	vs.resize(n);
    } else {
      if (vd.empty())
	vd.resize(n);
    }

    // Convert all inputs to UTC
    for (int i=0; i<n; i++) {
      switch (inF) {
      case JD:
	vut[i].set(vd[i]);
	break;
      case MJD:
	vut[i].setMJD(vd[i]);
	break;
      case TDB:
	vut[i].set(ephem.tdb2jd(vd[i]));
	break;
      default:
	// String-based ones
	std::istringstream iss(vs[i]);
	iss >> vut[i];
      }
    }

    // Now convert to outputs
    for (int i=0; i<n; i++) {
      switch (outF) {
      case JD:
	vd[i] = vut[i].getJD();
	break;
      case MJD:
	vd[i] = vut[i].getMJD();
	break;
      case TDB:
	vd[i] = ephem.jd2tdb(vut[i].getJD());
	break;
      case YMD:
	{
	  std::ostringstream oss;
	  oss << std::setprecision(7);
	  vut[i].writeYMD(oss);
	  vs[i] = oss.str();
	}
	break;
      case YMDHMS:
	{
	  std::ostringstream oss;
	  oss << std::setprecision(2);
	  vut[i].writeYMDHMS(oss);
	  vs[i] = oss.str();
	}
	break;
      default:
	throw runtime_error("Logic error reaching default in switch");
      }
    }

    // Output
    if (!outFilename.empty()) {
      // Write to a table
      FITS::FitsTable ff(outFilename,FITS::Create,1);
      auto table = ff.use();
      // Table length must match array length
      if (table.nrows() > 0 && n != table.nrows()) {
	FormatAndThrow<std::runtime_error>() << "Input length does not match output table size, "
					     << n
					     << " vs "
					     << table.nrows();
      }
      // Blow away old column if it exists
      if (table.hasColumn(outColumnName))
	table.eraseColumn(outColumnName);
      if (isStringFormat(inF)) {
	table.addColumn(vs,outColumnName);
      } else {
	table.addColumn(vd,outColumnName);
      }
    } else {
      // Write to stdout
      if (isStringFormat(outF)) {
	for (auto s : vs)
	  cout << s << endl;
      } else {
	// Set output string to about 0.01s precision
	if (outF==TDB)
	  cout << std::fixed << std::setprecision(9) << std::setw(12);
	else if (outF==JD)
	  cout << std::fixed << std::setprecision(7) << std::setw(15);
	else if (outF==MJD)
	  cout << std::fixed << std::setprecision(7) << std::setw(13);

	for (auto d : vd)
	  cout << d << endl;
      }
    }

  } catch (std::runtime_error& m) {
    quit(m);
  }

  exit(0);
}
