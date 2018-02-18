/* Use Horizons to get positions for SS body barycenters (Nept & Pluto),
 * then run these positions through our routines to get apparent positions,
 * then compare these positions with what Horizons gives.
 *  This should check my machinery for generating positions from 3d data to see
 * if this is the source of my residuals.
 * 6/9/00 gmb
 */
#include "orbfit.h"

void
body3d(double t,	/* time is in years here */
	int body,
	double *x, double *y, double *z);
void
body2d(int body,
      OBSERVATION *obs,
      double *x,
      double *y);

char   *help[] = {
    " check_posn:  Fit KBO orbit to observed astrometric RA/DEC",
    " usage:  check_posn bodycode <obsfile>",
    "  bodycode    is 8 for Neptune, 9 for Pluto barycenters.",
    "  obsfile     contains one observation per line, with format",
    "              JD  RA Dec error obscode",
    "  output is a file containing position residuals.",
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

OBSERVATION obsarray[MAXOBS];
int	nobs;

int
main(int argc, char *argv[])
{
  PBASIS p;
  ORBIT  orbit;
  XVBASIS xv;

  double x,y,t,*a,**covar;
  double chisq;
  int	i,j,bodycode;

  if (argc!=3) print_help();

  bodycode = atoi(argv[1]);
  if (bodycode>9 || bodycode<=0) {
    fprintf(stderr,"Bad body spec for JPL files: %s\n",argv[1]);
    print_help();
  }
  if (read_radec(obsarray,argv[2],&nobs)) {
    fprintf(stderr, "Error reading data file %s\n",argv[2]);
    exit(1);
  }
  
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

  /* Print out information on the coordinate system */
  printf("#     lat0       lon0       xBary     yBary      zBary   JD0\n");
  printf("%11.6f  %11.6f  %9.6f %9.6f %9.6f  %.5f\n",
	 lat0/DTOR,lon0/DTOR,xBary,yBary,zBary,jd0);

  /* Dump residuals to stdout */
  fprintf(stdout,"#Best fit orbit gives:\n");
  fprintf(stdout,"#obs  time        x      x_resid       y   y_resid  elong:\n");
  for (i=0; i<nobs; i++) {
    double x1, y1, z1, x2, y2, z2;
    body2d(bodycode, &obsarray[i], &x, &y);
    body3d(obsarray[i].obstime, bodycode,  &x1,&y1,&z1);
    x1 -= obsarray[i].xe; y1 -= obsarray[i].ye; z1 -= obsarray[i].ze; 
    body3d(obsarray[i].obstime+DAY, bodycode,  &x2,&y2,&z2);
    x2 -= x1; y2 -= y1; z2 -= z1;
    earth3d(obsarray[i].obstime+DAY, obsarray[i].obscode,  &x1,&y1,&z1);
    x2 -= x1; y2 -= y1; z2 -= z1;
    x2 *= 1.5e8/86400; y2 *= 1.5e8/86400; z2 *= 1.5e8/86400;
    fprintf(stdout,"%3d %9.5f %9.2f %6.3f %9.2f %6.3f %6.2f %6.2f %6.2f\n",
	    i, obsarray[i].obstime,
	    obsarray[i].thetax/ARCSEC, (obsarray[i].thetax-x)/ARCSEC,
	    obsarray[i].thetay/ARCSEC, (obsarray[i].thetay-y)/ARCSEC,
	    x2, y2, z2);
  }

  exit(0);
} 

/* Project KBO position onto tangent plane at a given time */
void
body2d(int body,
       OBSERVATION *obs,
       double *x,
       double *y)
{
  double        xk,yk,zk;
  double        xe,ye,ze;
  double        invz;
  int           i;
  double distance; 

  /* Get the Earth position if not already calculated*/
  if (obs->xe < -900.) {
    earth3d(obs->obstime, obs->obscode, &xe,&ye,&ze);
    obs->xe = xe;
    obs->ye = ye;
    obs->ze = ze;
  } else {
    xe = obs->xe;
    ye = obs->ye;
    ze = obs->ze;
  }

  /* Get preliminary KBO position */
  body3d(obs->obstime,body,
	&xk,&yk,&zk);
  /*fprintf(stderr," body3d gives %.4f %4.f %.4f\n",xk,yk,zk);*/

  /* At this point one should account for the light-travel time
     delay between observed time and kbo orbit position time.
     Calculate distance & time delay using naive position.
  */ 
  distance=sqrt( (xk-xe)*(xk-xe) + (yk-ye)*(yk-ye)
		 + (zk-ze)*(zk-ze) );

  body3d(obs->obstime-distance/SPEED_OF_LIGHT, body,
	&xk,&yk,&zk);
 
  /*{ 
    double dx,dy,dz;
    dx = xk-0.5*(xe+xBary);
    dy = yk-0.5*(ye+yBary);
    dz = zk-0.5*(ze+zBary);
    distance=sqrt( dx*dx + dy*dy+dz*dz);
    }*/

  distance=sqrt( (xk-xe)*(xk-xe) + (yk-ye)*(yk-ye)
    + (zk-ze)*(zk-ze) );
  body3d(obs->obstime-distance/SPEED_OF_LIGHT, body,
  &xk,&yk,&zk);
 
  invz = 1./(zk-ze);
  *x = (xk - xe)*invz;
  *y = (yk - ye)*invz;
 
  /* Try reversing the calculation to get back the equatorial barycentric coords*/
  zk = sqrt( ( (xk-xe)*(xk-xe) + (yk-ye)*(yk-ye)
	       + (zk-ze)*(zk-ze)) /
	     (*x * *x + *y * *y + 1) );
  xk = *x * zk;
  yk = *y * zk; /* now have x/y/z from telescope in projected system*/
  xk += xe; /* this is vector to object from initial point*/
  yk += ye; 
  zk += ze; 
  xk -= xBary; /* get to barycentric center*/
  yk -= yBary;
  zk -= zBary;
  {/* rotate to equatorial:*/
    double xec, yec, zec, xxx[3];
    xyz_proj_to_ec(xk, yk, zk, &xec, &yec, &zec, lat0, lon0, NULL);
    xyz_ec_to_eq(xec, yec, zec, &xk, &yk, &zk, NULL);
    
    /* compare to original: */
    bodycenter_ssbary((obs->obstime-distance/SPEED_OF_LIGHT)/DAY+jd0, xxx, body);
    xxx[0]-= xk;
    xxx[1]-= yk;
    xxx[2]-= zk;
    /*fprintf(stderr," Orig/reversed x,y,z: %f %f %f %f %f %f\n",
      xk,xxx[0],yk,xxx[1],zk,xxx[2]);*/
  }

  return;
}
/* Function to return xyz coords of a JPL ephemeris body in
 * standard coordinate system.
 */
void
body3d(double t,	/* time is in years here */
	int body,
	double *x, double *y, double *z)
{
  double xxx[3];
  double xec, yec, zec; 

  /* get observatory posn wrt barycenter */
  bodycenter_ssbary(t/DAY+jd0, xxx, body);
  /*fprintf(stderr,"bodycenter for %d: %f %f %f\n",body,xxx[0],xxx[1],xxx[2]);*/
  /* convert to tangent-point coord system */
  /* via ecliptic */
  xyz_eq_to_ec(xxx[0], xxx[1], xxx[2], &xec, &yec, &zec,NULL);
  xyz_ec_to_proj(xec, yec, zec, x, y, z, lat0, lon0, NULL);

  /* Translate to our origin */
  *x += xBary;
  *y += yBary;
  *z += zBary;

  return;
}
