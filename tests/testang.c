/* test coordinate transformation routines */
#include "orbfit.h"

int
main(int argc, char **argv) {
  double	ra0, dec0, lat0, lon0;
  double	lat,lon,ra,dec,xp, yp;
  double	**deriv1, **deriv2, **prod;

  deriv1 = dmatrix(1,2,1,2);
  deriv2 = dmatrix(1,2,1,2);
  prod = dmatrix(1,2,1,2);

  printf("Enter ra0, dec0 for projection: (in degrees) ");
  while (scanf("%lf %lf", &ra0, &dec0)==2) {

    ra0 *= PI / 180.;
    dec0 *= PI / 180.;
    
    /* convert eq to ecliptic */
    eq_to_ec(ra0, dec0, &lat0, &lon0, deriv1);
    
    printf("Ecliptic coords are: lon %.3f lat %.3f\n",
	   lon0*180./PI, lat0*180./PI);
    printf("derivative matrix:\n");
    print_matrix(stdout,deriv1,2,2);
    printf(" det=%f, cos dec / cos lat = %f\n",
	   deriv1[1][1]*deriv1[2][2]-deriv1[1][2]*deriv1[2][1],
	   cos(dec0)/cos(lat0));

    /* go back to eq */
    ec_to_eq(lat0,lon0,&ra0, &dec0, deriv2);
    
    printf("Equatorial are ra %.3f dec %.3f\n",
	   ra0*180./PI, dec0*180./PI);
    printf("derivative matrix:\n");
    print_matrix(stdout,deriv2,2,2);
    printf(" det=%f, cos dec / cos lat = %f\n",
	   deriv2[1][1]*deriv2[2][2]-deriv2[1][2]*deriv2[2][1],
	   cos(lat0)/cos(dec0));

    matrix_multiply(deriv1,deriv2,prod,2,2,2,2);
    printf("Product matrix:\n");
    print_matrix(stdout,prod,2,2);

  }

  exit(0);
}
