/* plutinos.c: plot distance of a Plutino and longitude relative 
 * to Neptune at various orbit phases.
 * 1/21/02 gmb
 */

#include "orbfit.h"

/*suffix for orbit-fit files*/
#define SUFFIX ".abg"
#define OBSCODE_GRND 807	/*Where the followup observation will be*/
#define PHASE_STEPS 100	   /*number of phases for each object's orbit*/

char   *help[] = {
  "plutinos.c: plot distance of a Plutino and longitude relative ",
  "  to Neptune at time steps around orbit.",
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
main(int argc, char *argv[])
{
  PBASIS p;
  ORBIT o, onept;
  XVBASIS xv, xvnept;

  FILE *abgfile;
  char	inbuff[256],objname[256],fname[256];

  double lat,lon, latnept, lonnept, distance;
  int i,nfields;
  double ra, dec, period, jd, peri;
  int iphase;
  double **covar;
  covar = dmatrix(1,6,1,6);

  if (argc>1) print_help();

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

  printf("#   ID    phase distance  /peri  Neptune longitude\n");

  /* Get Neptune: */
  if (read_abg("nept.abg",&p,covar)) {
    fprintf(stderr, "Error reading alpha/beta/gamma file %s\n",fname);
    exit(1);
  }

  /* look at several points around the orbital phase */
  /* get orbital elements for this one */
  pbasis_to_bary(&p, &xv, NULL);
  orbitElements(&xv, &onept);

  /* Loop through input containing names of objects */
  while ( fgets_nocomment(inbuff,255,stdin,NULL)!=NULL) {
    sscanf(inbuff," %s ",objname);
    strcpy(fname,objname);
    strcat(fname,SUFFIX);

    if (read_abg(fname,&p,covar)) {
      fprintf(stderr, "Error reading alpha/beta/gamma file %s\n",fname);
      exit(1);
    }

    /* look at several points around the orbital phase */
    /* get orbital elements for this one */
    pbasis_to_bary(&p, &xv, NULL);
    orbitElements(&xv, &o);

    period = pow(o.a, 1.5)/DAY;

    for (iphase=0; iphase<2*PHASE_STEPS; iphase++) {
      /* Advance jd through 2 periods of the orbit. */
      jd = o.T + iphase*period/PHASE_STEPS;

      elements_to_xv(&onept, jd, &xv);
      distance = sqrt(xv.x*xv.x + xv.y*xv.y + xv.z*xv.z);
      dec = asin(xv.z / distance);
      ra = atan2(xv.y,xv.x);
      eq_to_ec(ra, dec, &latnept, &lonnept, NULL);

      elements_to_xv(&o, jd, &xv);
      distance = sqrt(xv.x*xv.x + xv.y*xv.y + xv.z*xv.z);
      if (iphase==0) peri=distance;
      dec = asin(xv.z / distance);
      ra = atan2(xv.y,xv.x);
      eq_to_ec(ra, dec, &lat, &lon, NULL);

      lon -= lonnept;
      while (lon > PI) lon -= TPI;
      while (lon < -PI) lon += TPI;
   
      printf("%10s %4.2f %5.2f %5.3f %6.1f\n",
	     objname, ((double) iphase)/PHASE_STEPS, 
	     distance, distance/peri, lon/DTOR);
    } /*phase loop */
  } /*object loop*/
  exit(0);
}
