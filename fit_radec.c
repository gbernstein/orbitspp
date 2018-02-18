/* 	$Id: fit_radec.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: fit_radec.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* Fit orbit to RA/DEC file.  First
 * preliminary fit, THEN do a full Marquandt-Levenberg fit.
 * Use sigma matrix of full fit to predict positions and
 * errors using the full nonlinear formulae.  
 * usage:  orbtest3 <datafile> <# params>
 *  datafile contains list of observations (e.g. output of orbtest1)
 *  # params is either 5 (don't fit to gdot) or 6 (fit all params)
 * 8/9/99 gmb
 */
#include "orbfit.h"

char   *help[] = {
    " fit_radec:  Fit KBO orbit to observed astrometric RA/DEC",
    " usage:  fit_radec [-m mpc_error] [-j JPL_file] [-o observatory_file] [-v]",
    "  mpc_error  is arcsec uncertainty to apply to MPC-format",
    "             observations.  Default is 0.2",
    "  JPL_file   is binary ephemeris file.  Default is binEphem.405, or",
    "             a file specified by environment variable ORBIT_EPHEMERIS",
    "  observatory_file   is file with observatory site data.  Default is",
    "             observatories.dat, or a file specified by environment",
    "             variable ORBIT_OBSERVATORIES",
    "  stdin      contains one observation per line, with format",
    "             JD  RA Dec error obscode",
    "  stdout     is a file containing best-fit alpha, beta, etc.",
    "             and its covariance matrix.",
    "  Residuals are dumped to stderr.",
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


int
fit_observations(OBSERVATION obsarray[],
		 int nobs,
		 PBASIS *p,
		 double **covar,
		 double *chisq,
		 int *dof,
		 FILE *logfile);

OBSERVATION obsarray[MAXOBS];
int	nobs;

int
main(int argc, char *argv[])
{
  PBASIS p;
  ORBIT  orbit;
  XVBASIS xv;

  double **covar;
  double chisq;
  int i, dof;
  double rmserr, maxerr;

  covar = dmatrix(1,6,1,6);

  {
    int iarg=1;
    if (argc>1 && *argv[1]=='^') print_help();
    if (read_options(&iarg, argc, argv)) print_help();
  }

  if (read_radec(obsarray,NULL,&nobs)) {
    fprintf(stderr, "Error reading input observations\n");
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

  printf("# Fitting %d observations\n",nobs);

  /* Call subroutine to do the actual fitting: */
  fit_observations(obsarray, nobs, &p, covar, &chisq, &dof,stdout);

  printf("# Chi-squared of fit: %.2f DOF: %d\n",chisq,dof);
  printf("# Exact a, adot, b, bdot, g, gdot:\n");
  printf("%11.8f %11.8f %11.8f %11.8f %11.8f %11.8f\n",p.a,p.adot,p.b,
	p.bdot, p.g, p.gdot);
  pbasis_to_bary(&p, &xv, NULL);

  orbitElements(&xv, &orbit);
  printf("# a=%lf AU,e=%lf,i=%lf deg\n",orbit.a, orbit.e, orbit.i);
  {
    double d, dd;
    d = sqrt(xBary*xBary + yBary*yBary + pow(zBary-1/p.g,2.));
    dd = d*d*sqrt(covar[5][5]);
    printf("# Barycentric distance %.3f+-%.3f\n",d,dd);
  }

  /* Print the covariance matrix to stdout */
  printf("# Covariance matrix: \n");
  print_matrix(stdout,covar,6,6);

  /* Print out information on the coordinate system */
  printf("#     lat0       lon0       xBary     yBary      zBary   JD0\n");
  printf("%12.7f %12.7f %10.7f %10.7f %10.7f  %.6f\n",
	 lat0/DTOR,lon0/DTOR,xBary,yBary,zBary,jd0);

  /* Dump residuals to stderr */
  fprintf(stderr,"Best fit orbit gives:\n");
  fprintf(stderr,"obs  time        x      x_resid       y   y_resid\n");
  rmserr = 0.;
  maxerr = 0.;
  for (i=0; i<nobs; i++) {
    double x,y,dx,dy,errsq;
    kbo2d(&p, &obsarray[i], &x, NULL, &y, NULL);
    dx = (obsarray[i].thetax-x)/ARCSEC;
    dy = (obsarray[i].thetay-y)/ARCSEC;
    fprintf(stderr,"%3d %9.4f %10.3f %7.3f %10.3f %7.3f\n",
	    i, obsarray[i].obstime,
	    obsarray[i].thetax/ARCSEC, dx,
	    obsarray[i].thetay/ARCSEC, dy);
    errsq = dx*dx + dy*dy;
    rmserr += errsq;
    if (rmserr > maxerr) maxerr = rmserr;
  }
  rmserr = sqrt(rmserr/(nobs-1));
  maxerr = sqrt(maxerr);
  printf("# RMS Error: %.3f  Maximum err: %.3f (arcsec)\n",rmserr, maxerr);

  free_dmatrix(covar,1,6,1,6);

  exit(0);
} 

