// Functions to produce ASCII files in the format produced by the legacy C orbfit code.
#include "Legacy.h"
#include "Astrometry.h"
#include <cstdio>
#include <vector>

using namespace astrometry;
using namespace orbits;
using namespace std;

void
orbits::writeOldAEI(string filename, const Elements& elements,
	    const ElementCovariance& cov,
	    double tdb0,
	    Ephemeris& eph)
{
  // The order of elements in cov matrices:
  vector<int> indices={Elements::A, Elements::E, Elements::I,
		       Elements::LAN, Elements::AOP, Elements::TOP};
  auto fptr=std::fopen(filename.c_str(), "w");

  fprintf(fptr,"# Osculating elements at epoch %.1f:\n",
	  eph.tdb2jd(tdb0));
  fprintf(fptr,"#    a            e       i      Node   Arg of Peri   Time of Peri\n");
  fprintf(fptr,"%12.6f  %9.6f  %8.3f %8.3f  %8.3f %11.3f\n",
	  elements[Elements::A],
	  elements[Elements::E],
	  elements[Elements::I]/DEGREE,
	  elements[Elements::LAN]/DEGREE,
	  elements[Elements::AOP]/DEGREE,
	  eph.tdb2jd(elements[Elements::TOP]));
  fprintf(fptr,"+-%10.6f  %9.6f  %8.3f %8.3f  %8.3f %11.3f\n",
	  sqrt(cov(Elements::A,Elements::A)),
	  sqrt(cov(Elements::E,Elements::E)),
	  sqrt(cov(Elements::I,Elements::I))/DEGREE,
	  sqrt(cov(Elements::LAN,Elements::LAN))/DEGREE,
	  sqrt(cov(Elements::AOP,Elements::AOP))/DEGREE,
	  sqrt(cov(Elements::TOP,Elements::TOP))/DAY);
  fprintf(fptr,"# covariance matrix:\n");
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      fprintf(fptr,"%11.4e ",cov(indices[i],indices[j]));
    }
    fprintf(fptr,"\n");
  }
  fclose(fptr);
}

void
orbits::writeOldABG(string filename, const ABG& abg,
		    const ABGCovariance& abgcov,
		    const Frame& frame,
		    Ephemeris& eph)
{
  // ABG ordering for C code:
  vector<int> indices={ABG::A, ABG::ADOT,
		       ABG::B, ABG::BDOT,
		       ABG::G, ABG::GDOT};
  auto fptr=std::fopen(filename.c_str(), "w");
  std::fprintf(fptr,"# Exact a, adot, b, bdot, g, gdot:\n");
  std::fprintf(fptr,"%11.8f %11.8f %11.8f %11.8f %11.8f %11.8f\n",
	       abg[indices[0]],
	       abg[indices[1]],
	       abg[indices[2]],
	       abg[indices[3]],
	       abg[indices[4]],
	       abg[indices[5]]);

  fprintf(fptr, "# Covariance matrix: \n");
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      fprintf(fptr,"%11.4e ",abgcov(indices[i],indices[j]));
    }
    fprintf(fptr,"\n");
  }
 
  /* Print out information on the coordinate system */
  fprintf(fptr,"# Orbital reference frame:\n");
  fprintf(fptr,"#     lat0       lon0       xBary     yBary      zBary   JD0\n");
  double lat0, lon0, xBary, yBary, zBary, tdb0;
  // Old code used reference frames in ecliptic coords
  SphericalEcliptic pole(frame.orient.getPole());
  pole.getLonLat(lat0, lon0);
  double pa = frame.orient.getPA(); // Not written, assumed ecliptic orient!
  xBary = frame.origin.getVector()[0];
  yBary = frame.origin.getVector()[1];
  zBary = frame.origin.getVector()[2];

  fprintf(fptr,"%12.7f %12.7f %10.7f %10.7f %10.7f  %.6f\n",
	  lat0/DEGREE,lon0/DEGREE,xBary,yBary,zBary,eph.tdb2jd(frame.tdb0));
	
  std::fclose(fptr);
}
