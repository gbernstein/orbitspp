/* Compare Horizons orbit to my integrated orbit.
 * 6/9/00 gmb
 */
#include "orbfit.h"

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

  double x,y,t,*a,**covar;
  double chisq;
  int	i,j,bodycode;
  double xk[3], vk[3], xjpl, yjpl, zjpl, vxyz[3];

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

  /* Set up the orbit parameters using Horizons initial pos/velocity*/
  body3d(obsarray[0].obstime, bodycode,  &xjpl, &yjpl, &zjpl, vxyz);
  p.a = xjpl / zjpl;
  p.adot = vxyz[0] / zjpl;
  p.b = yjpl / zjpl;
  p.bdot = vxyz[1] / zjpl;
  p.g = 1./zjpl;
  p.gdot = vxyz[2] / zjpl;

  /* Calculate orbit, show residuals */
  fprintf(stdout,"time     x    x_resid    y   y_resid    z    z_resid:\n");
  for (i=0; i<nobs; i++) {
    double xg, yg, zg;
    kbo3d(&p, obsarray[i].obstime, xk, vk, NULL, NULL, NULL);
    body3d(obsarray[i].obstime, bodycode,  &xjpl,&yjpl,&zjpl, NULL);
    /* get a total acceleration term: */
    xg = xjpl - p.a/p.g - p.adot*obsarray[i].obstime/p.g;
    yg = yjpl - p.b/p.g - p.bdot*obsarray[i].obstime/p.g;
    zg = zjpl - 1./p.g - p.gdot*obsarray[i].obstime/p.g;
    fprintf(stdout,
	    "%9.5f %12.8f %.3g %.3g %12.8f  %.3g %.3g %12.8f  %.3g %.3g\n",
	    obsarray[i].obstime,
	    xk[0], xk[0]-xjpl, xg, xk[1], xk[1]-yjpl, yg, xk[2], 
	    xk[2]-zjpl, zg);

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
	&xk,&yk,&zk, NULL);
  /*fprintf(stderr," body3d gives %.4f %4.f %.4f\n",xk,yk,zk);*/

  /* At this point one should account for the light-travel time
     delay between observed time and kbo orbit position time.
     Calculate distance & time delay using naive position.
  */ 
  distance=sqrt( (xk-xe)*(xk-xe) + (yk-ye)*(yk-ye)
		 + (zk-ze)*(zk-ze) );

  body3d(obs->obstime-distance/SPEED_OF_LIGHT, body,
	&xk,&yk,&zk, NULL);
 
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
  &xk,&yk,&zk, NULL);
 
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
    bodycenter_ssbary((obs->obstime-distance/SPEED_OF_LIGHT)/DAY+jd0, xxx, body, NULL);
    xxx[0]-= xk;
    xxx[1]-= yk;
    xxx[2]-= zk;
    /*fprintf(stderr," Orig/reversed x,y,z: %f %f %f %f %f %f\n",
      xk,xxx[0],yk,xxx[1],zk,xxx[2]);*/
  }

  return;
}
