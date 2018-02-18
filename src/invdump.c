/* Give listing of eq & ec coords corresponding to points on invariable plane*/
/* Usage:  invdump lonmin lonmax step 
 * all in degrees. Default step is 1 degree.
 * 1/20/02 gmb
 */
#include "orbfit.h"
#include <stdio.h>

extern double 
dmsdeg(char *string);
extern double 
hmsdeg(char *string);
extern void
degdms(double degr,
       char *outbuff);
extern void
deghms(double degr,
       char *outbuff);

int
main(int argc,
     char *argv[])
{
  double lonmin, lonmax, lonstep;
  double lon, ra, dec, lon_ecl, lat_ecl;
  char rastring[64], decstring[64];
  if (argc<3) {
    fprintf(stderr,"Missing command-line arguments\n");
    exit(1);
  }
  lonmin = atof(argv[1]);
  lonmax = atof(argv[2]);
  if (argc>=4)
    lonstep = atof(argv[3]);
  else
    lonstep = 1.;

  printf("# Invariable              J2000                Ecliptic  \n",
	 "# lat   lon       RA              Dec        lat     lon \n");
  for (lon=lonmin; lon<=lonmax; lon+=lonstep) {
    inv_to_eq(0., lon*DTOR, &ra, &dec);
    while (ra<0) ra+=TPI;
    deghms(ra/DTOR, rastring);
    degdms(dec/DTOR, decstring);
    eq_to_ec(ra, dec, &lat_ecl, &lon_ecl, NULL);
    while (lon_ecl<0) lon_ecl+=TPI;
    printf("%6.2f %6.2f %12s %12s %6.2f %7.2f\n",0., lon, rastring,
	   decstring, lat_ecl/DTOR, lon_ecl/DTOR);
  }
}
