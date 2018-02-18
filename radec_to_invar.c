/* Transform J2000 equatorial coordinates to invariable lat & lon*/
#include "orbfit.h"
extern double 
dmsdeg(char *string);
extern double 
hmsdeg(char *string);

int
main(int argc,
     char *argv[])
{
  char name[64], rastring[64], decstring[64];
  double ra,dec,lat,lon;
  while (scanf("%s %s %s",name,rastring,decstring)==3) {
    ra = hmsdeg(rastring)*DTOR;
    dec = dmsdeg(decstring)*DTOR;
    eq_to_inv(ra,dec,&lat,&lon);
    printf("%10s %6.2f %6.2f\n",name,lat/DTOR,lon/DTOR);
  }
}
