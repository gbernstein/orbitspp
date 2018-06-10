// Ephemeris interface through CSPICE package //
// Note that spice uses TDB as seconds since J2000 TDB.
// Our code uses year as unit of time.  So our code
// will express TDB as Julian years past J2000 TDB, and
// convert on input/output of these codes.

#include "Ephemeris.h"
#include "AstronomicalConstants.h"
#include "StringStuff.h"
#include "Std.h"
#include <cstring>

/* / cspice files:
#include "SpiceZdf.h"
#include "SpiceCK.h"
#include "SpiceZpr.h"
*/
#include "SpiceUsr.h"

using namespace orbits;
using namespace std;

extern "C" SpiceDouble spd_c(void);

void checkSpice() {
  // Routine to check spice error status and throw exception if bad
  if (!failed_c()) return;
  int msglen = 2048;
  char msg[msglen];
  getmsg_c("SHORT",msglen, msg);
  string err = msg;
  getmsg_c("EXPLAIN",msglen,msg);
  if (strlen(msg)>0)
    err = err + "\n" + msg;
  getmsg_c("LONG",msglen,msg);
  if (strlen(msg)>0)
    err = err + "\n" + msg;

  // Reset error status in case someone decides to catch this
  // exception and proceed
  reset_c();
  
  throw SpiceError(err);
}

bool
Ephemeris::init = false;

// Name of the environment variable that can hold kernel filename:
const string SPICE_ENVIRON = "SPICE_KERNEL";
// Name of environment variable giving path to observatories file
const string OBS_ENVIRON = "ORBIT_OBSERVATORIES";

double Ephemeris::earthRadius;
double Ephemeris::earthFlatten;

Ephemeris::Ephemeris(const string& kernelFile) {
  if (init)
    throw std::runtime_error("Cannot create two instances of Ephemeris"
			     " under cspice");
  // Tell CSPICE not to abort on errors, instead go into "fall-through" mode
  const int STRING_LENGTH = 256;
  SpiceChar msg[STRING_LENGTH];
  strncpy(msg, "RETURN",STRING_LENGTH);
  erract_c("SET",1,msg);
  // Shut off error message dumping
  strncpy(msg, "NULL",STRING_LENGTH);
  errdev_c("SET",1,msg);

  
  if (kernelFile.empty())  {
    // Find kernel file name from environment variable
    char *kpath = getenv(SPICE_ENVIRON.c_str());
    if (kpath == NULL) 
      throw std::runtime_error("No path given for SPICE kernel file");
    furnsh_c(kpath);
  } else {
    furnsh_c(kernelFile.c_str());
  }

  /* Get radii and flattening of Earth */
  int n;
  double abc[3];
  bodvcd_c(EARTH, "RADII", 3, &n, abc );
  checkSpice();
  earthRadius = abc[0];
  earthFlatten =  ( abc[0]-abc[2]) / abc[0];

  init = true;
}
  
State
Ephemeris::state(int body,
		 double tdb) const {
  double s[6];
  double lighttime;
  spkgeo_c(body, tdb/SECOND, "J2000", SSB, s, &lighttime);
  checkSpice();
  State out;
  out.tdb = tdb;
  // Convert spice's km, s to AU, yr
  out.x[0] = s[0] * (1000. / AU);
  out.x[1] = s[1] * (1000. / AU);
  out.x[2] = s[2] * (1000. / AU);
  out.v[0] = s[3] * (1000. / AU / SECOND);
  out.v[1] = s[4] * (1000. / AU / SECOND);
  out.v[2] = s[5] * (1000. / AU / SECOND);
  return out;
}
  
astrometry::CartesianICRS
Ephemeris::position(int body,
		    double tdb) const {
  return state(body,tdb).x;
}

double
Ephemeris::utc2tdb(string utc) const {
  // Use spice string parsing method, which does all the leapseconds
  double tdb;
  str2et_c(utc.c_str(), &tdb);
  checkSpice();
  return tdb*SECOND;
}

double
Ephemeris::jd2tdb(double jd) const {
  // Assuming here that input JD is in UTC
  // turn into seconds past J2000 = 2000 JAN 01 12:00:00
  double utc = (jd - j2000_c()) * spd_c(); // (latter call is seconds per day)
  double delta;
  // Let spice calculate ET - UTC
  deltet_c(utc, "UTC", &delta);
  checkSpice();
  // Convert to years
  return (utc + delta) * SECOND;
}

double
Ephemeris::mjd2tdb(double mjd) const {
  return jd2tdb(mjd + MJD0);
}

double
Ephemeris::utc2tdb(const astrometry::UT& utc) const {
  return jd2tdb(utc.getJD());
}

astrometry::CartesianICRS
Ephemeris::observatory(double lon,
		       double lat,
		       double elev,
		       double tdb) const {
  // Convert lat and lon to cartesian:
  double abc[3], geo[3];
  double rotate[3][3];
  double lighttime;
  // Convert geodetic coordinates to Earth-fixed cartesian:
  georec_c( lon, lat, elev, earthRadius, earthFlatten, geo);
  // Get rotation matrix for Earth into J2000 
  pxform_c( "ITRF93", "J2000",  tdb/SECOND, rotate);
  // Apply rotate to get geo->obs vector in ICRS
  mxv_c(rotate, geo, abc);
  // Get geocenter position
  double s[6];
  spkgeo_c(EARTH, tdb/SECOND, "J2000", SSB, s, &lighttime);
  checkSpice();
  astrometry::CartesianICRS out;
  for (int i=0; i<3; i++)
    out[i] = (abc[i] + s[i]) * (1000. / AU);
  return out;
}

astrometry::CartesianICRS
Ephemeris::observatory(int obsid,
		       double tdb) const {
  
  // Read observatory file if we don't have it
  if (obsTable.empty()) {
    // Find kernel file name from environment variable
    char *kpath = getenv(OBS_ENVIRON.c_str());
    if (kpath == NULL) 
      throw std::runtime_error("No path given for observatory info file");
    // Now open the file named by the environment variable
    std::ifstream ifs(kpath);
    if (!ifs) {
      string fname(kpath);
      FormatAndThrow<std::runtime_error>() << "Could not open observatory info file " << fname;
    }
    string buffer;
    while (stringstuff::getlineNoComment(ifs, buffer)) {
      int id;
      string lonstring;
      string latstring;
      double elev;
      istringstream iss(buffer);
      if (!(iss >> id >> lonstring >> latstring >> elev))
	throw std::runtime_error("Bad observatory spec: <" + buffer + ">");

      ObsInfo oi;
      oi.lon = -astrometry::dmsdeg(lonstring) * DEGREE; // Table gives degrees W lon, we want rad E
      oi.lat = astrometry::dmsdeg(latstring) * DEGREE;
      oi.elev = elev / 1000.;   // cspice wants km; MPC uses meters
      obsTable[id] = oi;
    }
    ifs.close();
  }

  // Look up our id
  if (obsTable.count(obsid)==0) 
    FormatAndThrow<std::runtime_error>() << "Unknown observatory id " << obsid;
  auto obs = obsTable[obsid];
  // Call with lon/lat/elev
  return observatory( obs.lon, obs.lat, obs.elev, tdb);
}

Observation
orbits::mpc2Observation(const MPCObservation& mpc, const Ephemeris& ephem) {
  Observation out;
  /**/cerr << "In mpc2";
  out.radec = mpc.radec;
  out.tdb = ephem.mjd2tdb(mpc.mjd);
  /**/cerr << " got tdb " << out.tdb;
  out.cov(1,1) = out.cov(0,0) = pow(mpc.sigma*ARCSEC,2);
  out.cov(0,1) = out.cov(1,0) = 0.;
  out.observer = ephem.observatory(mpc.obscode, out.tdb);
  /**/cerr << " and earth at " << out.observer << endl;
  return out;
};
