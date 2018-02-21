/**** use amoeba fitter to see where we end up ****/
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


#define  MAXIT	100	/*maximum mrqmin iterations*/
#define	 CHITOL 0.0001	/*tolerance for chisq minimum*/

OBSERVATION obsarray[MAXOBS];
int	nobs;

char   *help[] = {
    " fit_radec:  Fit KBO orbit to observed astrometric RA/DEC",
    " usage:  fit_radec <obsfile>",
    "  obsfile     contains one observation per line, with format",
    "              JD  RA Dec error obscode",
    "  output is a file containing best-fit alpha, beta, etc.",
    "         and its covariance matrix.",
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

double
amoebafunc(double *a) {
  int i;
  PBASIS	params;
  double x,y, chisq=0., dx, dy;

  params.a = a[1];
  params.adot = a[2];
  params.b = a[3];
  params.bdot = a[4];
  params.g = a[5];
  params.gdot = a[6];

  for (i=0; i<nobs; i++) {
    kbo2d(&params, &obsarray[i], &x, NULL, &y, NULL);
    dx = (obsarray[i].thetax - x)/(obsarray[i].dthetax);
    dy = (obsarray[i].thetay - y)/(obsarray[i].dthetay);
    chisq += dx*dx+dy*dy;
  }
  return chisq;
}

void
amoeba(double **p, double y[], int ndim, double ftol,
       double (*funk)(double []), int *nfunk);
void 
mrqmin_orbit(OBSERVATION obsarray[], int ndata, double a[], int ia[],
		  int ma, double **covar, double **alpha, double *chisq,
		  double *alamda);


int
main(int argc, char *argv[])
{
  PBASIS p;
  ORBIT  orbit;
  XVBASIS xv;

  double x,y,t,*a,**covar, **alpha, alambda;
  double chisq;
  int	ma=6,*ia,i,j;

  int nfunc;
  double *yy, **pp;
  a     = dvector(1,6);
  ia    = ivector(1,6);
  alpha = dmatrix(1,6,1,6);
  covar = dmatrix(1,6,1,6);

  if (argc!=2) print_help();

  if (read_radec(obsarray,argv[1],&nobs)) {
    fprintf(stderr, "Error reading data file %s\n",argv[1]);
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

  prelim_fit(obsarray,nobs,&p,covar);

  fprintf(stderr,"Preliminary a, adot, b, bdot, g, gdot:\n");
  fprintf(stderr,"%lf %lf %lf %lf %lf %lf\n",p.a,p.adot,p.b,
	p.bdot, p.g, p.gdot);  

  /* Prepare for the amoeba */
  pp = dmatrix(1,7,1,6);
  yy = dvector(1,7);
  for (j=1; j<=7; j++) {
    pp[j][1] = p.a;
    pp[j][2] = p.adot;
    pp[j][3] = p.b;
    pp[j][4] = p.bdot;
    pp[j][5] = p.g;
    pp[j][6] = p.gdot;
  }
  pp[1][1] += 1e-5;
  pp[2][2] += 1e-4;
  pp[3][3] += 1e-5;
  pp[4][4] += 1e-4;
  pp[5][5] += 1e-4;
  pp[6][6] += 1e-3;
  for (i=1; i<=7; i++) yy[i] = amoebafunc(pp[i]);

  amoeba(pp, yy, 6, 0.0001, amoebafunc, &nfunc);

  /* restart: */
  pp[1][1] += 1e-5;
  pp[2][2] += 1e-4;
  pp[3][3] += 1e-5;
  pp[4][4] += 1e-4;
  pp[5][5] += 1e-4;
  pp[6][6] += 1e-3;

  for (i=1; i<=7; i++) yy[i] = amoebafunc(pp[i]);
  amoeba(pp, yy, 6, 0.0001, amoebafunc, &nfunc);
 
  p.a = pp[1][1];
  p.adot = pp[1][2];
  p.b = pp[1][3];
  p.bdot = pp[1][4];
  p.g = pp[1][5];
  p.gdot = pp[1][6];
 
  printf("# Exact a, adot, b, bdot, g, gdot:\n");
  printf("%10.7f %lf %10.7f %lf %lf %lf\n",p.a,p.adot,p.b,
	p.bdot, p.g, p.gdot);
  pbasis_to_bary(&p, &xv, NULL);

  orbitElements(&xv, &orbit);
  printf("# a=%lf AU,e=%lf,i=%lf deg\n",orbit.a, orbit.e, orbit.i);

  /* use mrqmin to get covariance matrix */
  for (i=1; i<=6; i++) ia[i]=1;

  a[1]=p.a;
  a[2]=p.adot;
  a[3]=p.b;
  a[4]=p.bdot;
  a[5]=p.g;
  a[6]=p.gdot;
  alambda=0.;
  /*  mrqmin_orbit(obsarray, nobs, a, ia,
      ma, covar, alpha, &chisq, &alambda); */

  /* Print the covariance matrix to stdout */
  /*  printf("# Covariance matrix: \n");
      print_matrix(stdout,covar,6,6); */

  /* Print out information on the coordinate system */
  printf("#     lat0       lon0       xBary     yBary      zBary   JD0\n");
  printf("%11.6f  %11.6f  %9.6f %9.6f %9.6f  %.5f\n",
	 lat0/DTOR,lon0/DTOR,xBary,yBary,zBary,jd0);

  /* Dump residuals to stderr */
  fprintf(stderr,"Best fit orbit gives:\n");
  fprintf(stderr,"obs  time        x      x_resid       y   y_resid\n");
  for (i=0; i<nobs; i++) {
    kbo2d(&p, &obsarray[i], &x, NULL, &y, NULL);
    fprintf(stderr,"%3d %9.4f %9.2f %6.2f %9.2f %6.2f\n",
	    i, obsarray[i].obstime,
	    obsarray[i].thetax/ARCSEC, (obsarray[i].thetax-x)/ARCSEC,
	    obsarray[i].thetay/ARCSEC, (obsarray[i].thetay-y)/ARCSEC);
  }


  free_ivector(ia,1,6);
  free_dmatrix(covar,1,6,1,6);

  exit(0);
} 

