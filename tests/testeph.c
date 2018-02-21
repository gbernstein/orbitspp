/* testeph.c: see if I'm getting RA & Dec of planets correct
 * 1/20/02 gmb
 */
#include "orbfit.h"
#include "ephem_types.h"

char   *help[] = {
  " usage: testeph obscode",
    0
};

void
print_help(void)
{
  int	i;
  for (i = 0; help[i] != 0; i++)
    fprintf (stderr, "%s\n", help[i]);
  exit (1);
}

extern void 
deghms(double degr,
       char *outbuff);
extern void 
degdms(double degr,
       char *outbuff);

int
main(int argc, char *argv[])
{
  PBASIS p;
  OBSERVATION	futobs;
  struct date_time dt;
  char	inbuff[256],rastring[20],decstring[20];
  double lat_ecl,lon_ecl,lat_inv, lon_inv;
  double ra,dec;
  double yr,mo,day,hr,mn,ss, jd;
  double obs[3], xyz[3];
  int i, nfields;

  if (argc<2 || *argv[1]=='^') print_help();
  futobs.obscode = atof(argv[1]);

  /* echo the command line to output */
  printf("#");
  for (i=0; i<argc; i++) printf(" %s",argv[i]);
  {
#include <time.h>
    time_t timettt;
    time(&timettt);
    /* note that ctime returns string with newline at end */
    printf("\n#---%s",ctime(&timettt));
  }

  fprintf (stderr,"Enter JD's or Y M D ... of observations, -1 to quit:\n");
  printf("#Body   Inv. lat/lon    RA   (J2000)    Dec         Ecliptic lat/lon\n");
  while ( fgets_nocomment(inbuff,255,stdin,NULL)!=NULL) {
    nfields=sscanf(inbuff,"%lf %lf %lf %lf %lf %lf",
		   &yr,&mo,&day,&hr,&mn,&ss);
    if (nfields==0 ) {
      fprintf(stderr,"Error on time spec:\n->%s\n",inbuff);
      exit(1);
    } else if (yr<0.) {
      /*done*/
      exit(0);
    } else if (nfields==1 || nfields==2) {
      /* Got a JD. (probably...)*/
      jd = yr;
      futobs.obstime = (yr-jd0)*DAY;
    } else {
      dt.y = yr;
      dt.mo = mo;
      dt.d = day;
      if (nfields>=4) dt.h = hr;  else dt.h=0.;
      if (nfields>=5) dt.mn = mn; else dt.mn=0.;
      if (nfields>=6) dt.s = ss;  else dt.s=0.;
      jd = date_to_jd(dt);
      futobs.obstime = (jd-jd0)*DAY;
    }

    /* Get observatory coordinates in ICRS:*/
    earth_ssbary(jd, futobs.obscode, &obs[0], &obs[1], &obs[2]);

    /* Apparent Location of the Sun: (ignore time delay)*/
    bodycenter_ssbary(jd, xyz, SUN, NULL);
    for (i=0; i<3; i++) xyz[i]-=obs[i];
    /* Find line of site in ICRS */
    {
      double r;
      r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
      dec = asin(xyz[2] / r);
      ra = atan2(xyz[1],xyz[0]);
      if (ra<0) ra+=TPI;
    }
    deghms(ra/DTOR, rastring);
    degdms(dec/DTOR, decstring);
    eq_to_ec(ra, dec, &lat_ecl, &lon_ecl, NULL);
    eq_to_inv(ra, dec, &lat_inv, &lon_inv);
    if (lon_inv<0.) lon_inv+=TPI;
    if (lon_ecl<0.) lon_ecl+=TPI;
    printf("Sun:  %6.2f %6.2f %12s %12s %6.2f %7.2f\n",
	   lat_inv/DTOR, lon_inv/DTOR, rastring,
	   decstring, lat_ecl/DTOR, lon_ecl/DTOR);

    /* Apparent Location of the Moon: (ignore time delay)*/
    bodycenter_ssbary(jd, xyz, MOON, NULL);
    for (i=0; i<3; i++) xyz[i]-=obs[i];
    /* Find line of site in ICRS */
    {
      double r;
      r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
      dec = asin(xyz[2] / r);
      ra = atan2(xyz[1],xyz[0]);
      if (ra<0) ra+=TPI;
    }
    deghms(ra/DTOR, rastring);
    degdms(dec/DTOR, decstring);
    eq_to_ec(ra, dec, &lat_ecl, &lon_ecl, NULL);
    eq_to_inv(ra, dec, &lat_inv, &lon_inv);
    if (lon_inv<0.) lon_inv+=TPI;
    if (lon_ecl<0.) lon_ecl+=TPI;
    printf("Moon: %6.2f %6.2f %12s %12s %6.2f %7.2f\n",
	   lat_inv/DTOR, lon_inv/DTOR, rastring,
	   decstring, lat_ecl/DTOR, lon_ecl/DTOR);

    /* Apparent Location of Neptune: (ignore time delay)*/
    bodycenter_ssbary(jd, xyz, NEPTUNE, NULL);
    for (i=0; i<3; i++) xyz[i]-=obs[i];
    /* Find line of site in ICRS */
    {
      double r;
      r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
      dec = asin(xyz[2] / r);
      ra = atan2(xyz[1],xyz[0]);
      if (ra<0) ra+=TPI;
    }
    deghms(ra/DTOR, rastring);
    degdms(dec/DTOR, decstring);
    eq_to_ec(ra, dec, &lat_ecl, &lon_ecl, NULL);
    eq_to_inv(ra, dec, &lat_inv, &lon_inv);
    if (lon_inv<0.) lon_inv+=TPI;
    if (lon_ecl<0.) lon_ecl+=TPI;
    printf("Nept: %6.2f %6.2f %12s %12s %6.2f %7.2f\n",
	   lat_inv/DTOR, lon_inv/DTOR, rastring,
	   decstring, lat_ecl/DTOR, lon_ecl/DTOR);

  }
  exit(0);
}
