/* Subroutine to determine partial derivative matrix of orbital params
 * w.r.t. barycentric phase space params. Incredible algebra by Bharat.
 * Adapted for new Orbits classes in 2018.
 */

#include "LinearAlgebra.h"
#include "OrbitTypes.h"
#include "AstronomicalConstants.h"
#include "Elements.h"

namespace orbits {
  Matrix66
  aei_derivs( const orbits::State& xv, bool heliocentric)
  {
    double x,y,z,xdot,ydot,zdot;
    double mu= heliocentric ? GM : SolarSystemGM; /* AU^3/yr^2 */

    // State vectors are ICRS, in AU, yr units.
    // Convert state vector from ICRS to ecliptic
    astrometry::Matrix33 partials;
    astrometry::CartesianEcliptic x_ecl(xv.x);
    astrometry::CartesianEcliptic v_ecl(xv.v, partials);

    x = x_ecl[0];
    y = x_ecl[1];
    z = x_ecl[2];
    xdot = v_ecl[0];
    ydot = v_ecl[1];
    zdot = v_ecl[2];

    // Common quantities:
    double x2=x*x;
    double y2=y*y;
    double z2=z*z;

    double rr=x2+y2+z2;
    double r = sqrt(rr);
    double r3 = pow(rr,1.5);

    double vx2=xdot*xdot;
    double vy2=ydot*ydot;
    double vz2=zdot*zdot;
    double vv = vx2+vy2+vz2;
    double xdotv = x*xdot + y*ydot + z*zdot;

    double E2=2/r - vv/mu;	//-2*energy/mu
    double ee=-(1/r) + vv/mu; 
    double E2sq=pow(E2,2);
    double E2cube=pow(E2,3);

    // Create output
    Matrix66 derivs;
    /* ax is del a by del x */
    derivs(0,0) = (2*x)/(r3*E2sq);
    derivs(0,1) = (2*y)/(r3*E2sq);
    derivs(0,2) = (2*z)/(r3*E2sq);
    derivs(0,3) = (2*xdot)/(mu*E2sq);
    derivs(0,4) = (2*ydot)/(mu*E2sq);
    derivs(0,5) = (2*zdot)/(mu*E2sq);
 
    /* Partials for e now */

    derivs(1,0) = 
      (2*(-(vx2/mu) + 
	  x2/r3 - 
	  1/r + vv/mu)*
       (-((xdot*xdotv)/mu) + 
	x*ee)
       + 2*(-((xdot*ydot)/mu) + 
	    (x*y)/r3)*
       (-((ydot*xdotv)/mu) + 
	y*ee)
       + 2*((x*z)/
	    r3 - 
	    (xdot*zdot)/mu)*(-((zdot*xdotv)/
			       mu) + z*ee)
       )/(2.*sqrt(pow(-((xdot*xdotv)/
			mu) + x*ee
		      ,2) + pow(-((ydot*xdotv)/mu) + 
				y*ee
				,2) + pow(-((zdot*xdotv)/mu) + 
					  z*ee
					  ,2)));
 
    derivs(1,1) = 
      (2*(-((xdot*ydot)/mu) + (x*y)/
	  r3)*
       (-((xdot*xdotv)/mu) + 
	x*ee)
       + 2*(-(vy2/mu) + 
	    y2/r3 - 
	    1/r +  vv/mu)*
       (-((ydot*xdotv)/mu) + 
	y*ee)
       + 2*((y*z)/
	    r3 - 
	    (ydot*zdot)/mu)*(-((zdot*xdotv)/
			       mu) + z*ee)
       )/(2.*sqrt(pow(-((xdot*xdotv)/
			mu) + x*ee
		      ,2) + pow(-((ydot*xdotv)/mu) + 
				y*ee
				,2) + pow(-((zdot*xdotv)/mu) + 
					  z*ee
					  ,2)));
 
    derivs(1,2) = 
      (2*((x*z)/r3 - 
	  (xdot*zdot)/mu)*(-((xdot*xdotv)/
			     mu) + x*ee)
       + 2*((y*z)/
	    r3 - 
	    (ydot*zdot)/mu)*(-((ydot*xdotv)/
			       mu) + y*ee)
       + 2*(z2/r3 - 
	    1/r - 
	    vz2/mu + 
	    vv/mu)*
       (-((zdot*xdotv)/mu) + 
	z*ee)
       )/(2.*sqrt(pow(-((xdot*xdotv)/
			mu) + x*ee
		      ,2) + pow(-((ydot*xdotv)/mu) + 
				y*ee
				,2) + pow(-((zdot*xdotv)/mu) + 
					  z*ee
					  ,2)));
 
    derivs(1,3) = 
      (2*((x*xdot)/mu - (xdotv)/mu)*
       (-((xdot*xdotv)/mu) + 
	x*ee)
       + 2*((2*xdot*y)/mu - (x*ydot)/mu)*
       (-((ydot*xdotv)/mu) + 
	y*ee)
       + 2*((2*xdot*z)/mu - (x*zdot)/mu)*
       (-((zdot*xdotv)/mu) + 
	z*ee)
       )/(2.*sqrt(pow(-((xdot*xdotv)/
			mu) + x*ee
		      ,2) + pow(-((ydot*xdotv)/mu) + 
				y*ee
				,2) + pow(-((zdot*xdotv)/mu) + 
					  z*ee
					  ,2)));
 
    derivs(1,4) = 
      (2*(-((xdot*y)/mu) + (2*x*ydot)/mu)*
       (-((xdot*xdotv)/mu) + 
	x*ee)
       + 2*((y*ydot)/mu - (xdotv)/mu)*
       (-((ydot*xdotv)/mu) + 
	y*ee)
       + 2*((2*ydot*z)/mu - (y*zdot)/mu)*
       (-((zdot*xdotv)/mu) + 
	z*ee)
       )/(2.*sqrt(pow(-((xdot*xdotv)/
			mu) + x*ee
		      ,2) + pow(-((ydot*xdotv)/mu) + 
				y*ee
				,2) + pow(-((zdot*xdotv)/mu) + 
					  z*ee
					  ,2)));
 
    derivs(1,5) = 
      (2*(-((xdot*z)/mu) + (2*x*zdot)/mu)*
       (-((xdot*xdotv)/mu) + 
	x*ee)
       + 2*(-((ydot*z)/mu) + (2*y*zdot)/mu)*
       (-((ydot*xdotv)/mu) + 
	y*ee)
       + 2*((z*zdot)/mu - (xdotv)/mu)*
       (-((zdot*xdotv)/mu) + 
	z*ee)
       )/(2.*sqrt(pow(-((xdot*xdotv)/
			mu) + x*ee
		      ,2) + pow(-((ydot*xdotv)/mu) + 
				y*ee
				,2) + pow(-((zdot*xdotv)/mu) + 
					  z*ee
					  ,2)));
 
    /* Partials of i now */

    derivs(2,0) = 
      -((-((-(xdot*y) + x*ydot)*
	   (2*ydot*(-(xdot*y) + x*ydot) - 
	    2*zdot*(xdot*z - x*zdot)))/
	 (2.*pow(pow(-(xdot*y) + x*ydot,2) + 
		 pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2),1.5)) + 
	 ydot/sqrt(pow(-(xdot*y) + x*ydot,2) + 
		   pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
		   ))/
	sqrt(1 - pow(-(xdot*y) + x*ydot,2)/
	     (pow(-(xdot*y) + x*ydot,2) + 
	      pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
	      )));
 
    derivs(2,1) = 
      -((-((-(xdot*y) + x*ydot)*
	   (-2*xdot*(-(xdot*y) + x*ydot) + 
	    2*zdot*(-(ydot*z) + y*zdot)))/
	 (2.*pow(pow(-(xdot*y) + x*ydot,2) + 
		 pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2),1.5)) - 
	 xdot/sqrt(pow(-(xdot*y) + x*ydot,2) + 
		   pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
		   ))/
	sqrt(1 - pow(-(xdot*y) + x*ydot,2)/
	     (pow(-(xdot*y) + x*ydot,2) + 
	      pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
	      )));
 
    derivs(2,2) = 
      ((-(xdot*y) + x*ydot)*(2*xdot*(xdot*z - x*zdot) - 
			     2*ydot*(-(ydot*z) + y*zdot)))/
      (2.*pow(pow(-(xdot*y) + x*ydot,2) + 
	      pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2),
	      1.5)*sqrt(1 - pow(-(xdot*y) + x*ydot,2)/
			(pow(-(xdot*y) + x*ydot,2) + 
			 pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
			 )));
 
    derivs(2,3) = 
      -((-((-(xdot*y) + x*ydot)*
	   (-2*y*(-(xdot*y) + x*ydot) + 2*z*(xdot*z - x*zdot)))/
	 (2.*pow(pow(-(xdot*y) + x*ydot,2) + 
		 pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2),1.5)) - 
	 y/sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
		))/
	sqrt(1 - pow(-(xdot*y) + x*ydot,2)/
	     (pow(-(xdot*y) + x*ydot,2) + 
	      pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
	      )));
 
    derivs(2,4) = 
      -((-((-(xdot*y) + x*ydot)*
	   (2*x*(-(xdot*y) + x*ydot) - 2*z*(-(ydot*z) + y*zdot))
	   )/
	 (2.*pow(pow(-(xdot*y) + x*ydot,2) + 
		 pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2),1.5)) + 
	 x/sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
		))/
	sqrt(1 - pow(-(xdot*y) + x*ydot,2)/
	     (pow(-(xdot*y) + x*ydot,2) + 
	      pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
	      )));
 
    derivs(2,5) = 
      ((-(xdot*y) + x*ydot)*(-2*x*(xdot*z - x*zdot) + 
			     2*y*(-(ydot*z) + y*zdot)))/
      (2.*pow(pow(-(xdot*y) + x*ydot,2) + 
	      pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2),
	      1.5)*sqrt(1 - pow(-(xdot*y) + x*ydot,2)/
			(pow(-(xdot*y) + x*ydot,2) + 
			 pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2)
			 )));
 

    /* Partials of capital Omega (long of asc node) now */

    derivs(3,0) = 
      -(((zdot*(xdot*z - x*zdot)*(-(xdot*z) + x*zdot))/
	 pow(pow(xdot*z - x*zdot,2) + 
	     pow(-(ydot*z) + y*zdot,2),1.5) + 
	 zdot/sqrt(pow(xdot*z - x*zdot,2) + 
		   pow(-(ydot*z) + y*zdot,2)))/
	sqrt(1 - pow(-(xdot*z) + x*zdot,2)/
	     (pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2))
	     ));
 
    derivs(3,1) = 
      (zdot*(-(xdot*z) + x*zdot)*(-(ydot*z) + y*zdot))/
      (pow(pow(xdot*z - x*zdot,2) + 
	   pow(-(ydot*z) + y*zdot,2),1.5)*
       sqrt(1 - pow(-(xdot*z) + x*zdot,2)/
	    (pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2))
	    ));
 
    derivs(3,2) = 
      -((-((-(xdot*z) + x*zdot)*
	   (2*xdot*(xdot*z - x*zdot) - 
	    2*ydot*(-(ydot*z) + y*zdot)))/
	 (2.*pow(pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2),1.5)) - 
	 xdot/sqrt(pow(xdot*z - x*zdot,2) + 
		   pow(-(ydot*z) + y*zdot,2)))/
	sqrt(1 - pow(-(xdot*z) + x*zdot,2)/
	     (pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2))
	     ));
 
    derivs(3,3) = 
      -((-((z*(xdot*z - x*zdot)*(-(xdot*z) + x*zdot))/
	   pow(pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2),1.5)) - 
	 z/sqrt(pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2)))/
	sqrt(1 - pow(-(xdot*z) + x*zdot,2)/
	     (pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2))
	     ));
 
    derivs(3,4) = 
      -((z*(-(xdot*z) + x*zdot)*(-(ydot*z) + y*zdot))/
	(pow(pow(xdot*z - x*zdot,2) + 
	     pow(-(ydot*z) + y*zdot,2),1.5)*
	 sqrt(1 - pow(-(xdot*z) + x*zdot,2)/
	      (pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2)))));
 
    derivs(3,5) = 
      -((-((-(xdot*z) + x*zdot)*
	   (-2*x*(xdot*z - x*zdot) + 2*y*(-(ydot*z) + y*zdot)))/
	 (2.*pow(pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2),1.5)) + 
	 x/sqrt(pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2)))/
	sqrt(1 - pow(-(xdot*z) + x*zdot,2)/
	     (pow(xdot*z - x*zdot,2) + pow(-(ydot*z) + y*zdot,2))
	     ));
 

    /* partials of small omega (arg of per) now */

    derivs(4,0) = 
      -((-(((-(xdot*z) + x*zdot)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) + 
	    (-(ydot*z) + y*zdot)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee))*
	   (2*(-(vx2/mu) + 
	       x2/
	       pow(rr,
		   1.5) - 1/r + vv/mu)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) + 
	    2*(-((xdot*ydot)/mu) + 
	       (x*y)/
	       pow(rr,
		   1.5))*(-((ydot*xdotv)/
			    mu) + y*
			  (ee)) + 
	    2*((x*z)/
	       pow(rr,
		   1.5) - (xdot*zdot)/mu)*
	    (-((zdot*xdotv)/mu) + 
	     z*ee)))/
	 (2.*sqrt(pow(xdot*z - x*zdot,2) + 
		  pow(-(ydot*z) + y*zdot,2))*
	  pow(pow(-((xdot*xdotv)/
		    mu) + x*ee,2) + 
	      pow(-((ydot*xdotv)/mu) + 
		  y*ee,2) + 
	      pow(-((zdot*xdotv)/mu) + 
		  z*ee,2),1.5)) + 
	 ((-((xdot*ydot)/mu) + 
	   (x*y)/
	   r3)*
	  (-(ydot*z) + y*zdot) + 
	  (-(xdot*z) + x*zdot)*
	  (-(vx2/mu) + 
	   x2/r3\
	   - 1/r +  vv/mu
	   ) + zdot*(-((xdot*xdotv)/
		       mu) + x*ee))/
	 (sqrt(pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	  sqrt(pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) + 
	 (zdot*(xdot*z - x*zdot)*
	  ((-(xdot*z) + x*zdot)*
	   (-((xdot*xdotv)/mu) + 
	    x*ee) + 
	   (-(ydot*z) + y*zdot)*
	   (-((ydot*xdotv)/mu) + 
	    y*ee)))/
	 (pow(pow(xdot*z - x*zdot,2) + 
	      pow(-(ydot*z) + y*zdot,2),1.5)*
	  sqrt(pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))))/
	sqrt(1 - pow((-(xdot*z) + x*zdot)*
		     (-((xdot*xdotv)/mu) + 
		      x*ee) + 
		     (-(ydot*z) + y*zdot)*
		     (-((ydot*xdotv)/mu) + 
		      y*ee),2)/
	     ((pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	      (pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2)))));
 
    derivs(4,1) = 
      -((-(((-(xdot*z) + x*zdot)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) + 
	    (-(ydot*z) + y*zdot)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee))*
	   (2*(-((xdot*ydot)/mu) + 
	       (x*y)/
	       pow(rr,
		   1.5))*(-((xdot*xdotv)/
			    mu) + x*
			  (ee)) + 
	    2*(-(vy2/mu) + 
	       y2/
	       pow(rr,
		   1.5) - 1/r + vv/mu)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee) + 
	    2*((y*z)/
	       pow(rr,
		   1.5) - (ydot*zdot)/mu)*
	    (-((zdot*xdotv)/mu) + 
	     z*ee)))/
	 (2.*sqrt(pow(xdot*z - x*zdot,2) + 
		  pow(-(ydot*z) + y*zdot,2))*
	  pow(pow(-((xdot*xdotv)/
		    mu) + x*ee,2) + 
	      pow(-((ydot*xdotv)/mu) + 
		  y*ee,2) + 
	      pow(-((zdot*xdotv)/mu) + 
		  z*ee,2),1.5)) + 
	 ((-((xdot*ydot)/mu) + 
	   (x*y)/
	   r3)*
	  (-(xdot*z) + x*zdot) + 
	  (-(ydot*z) + y*zdot)*
	  (-(vy2/mu) + 
	   y2/r3\
	   - 1/r + 
	   vv/mu
	   ) + zdot*(-((ydot*xdotv)/
		       mu) + y*ee))/
	 (sqrt(pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	  sqrt(pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) - 
	 (zdot*(-(ydot*z) + y*zdot)*
	  ((-(xdot*z) + x*zdot)*
	   (-((xdot*xdotv)/mu) + 
	    x*ee) + 
	   (-(ydot*z) + y*zdot)*
	   (-((ydot*xdotv)/mu) + 
	    y*ee)))/
	 (pow(pow(xdot*z - x*zdot,2) + 
	      pow(-(ydot*z) + y*zdot,2),1.5)*
	  sqrt(pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))))/
	sqrt(1 - pow((-(xdot*z) + x*zdot)*
		     (-((xdot*xdotv)/mu) + 
		      x*ee) + 
		     (-(ydot*z) + y*zdot)*
		     (-((ydot*xdotv)/mu) + 
		      y*ee),2)/
	     ((pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	      (pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2)))));
 
    derivs(4,2) = 
      -((-(((-(xdot*z) + x*zdot)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) + 
	    (-(ydot*z) + y*zdot)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee))*
	   (2*((x*z)/
	       pow(rr,
		   1.5) - (xdot*zdot)/mu)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) + 
	    2*((y*z)/
	       pow(rr,
		   1.5) - (ydot*zdot)/mu)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee) + 
	    2*(z2/
	       pow(rr,
		   1.5) - 1/r - 
	       vz2/mu + 
	       vv/mu)*
	    (-((zdot*xdotv)/mu) + 
	     z*ee)))/
	 (2.*sqrt(pow(xdot*z - x*zdot,2) + 
		  pow(-(ydot*z) + y*zdot,2))*
	  pow(pow(-((xdot*xdotv)/
		    mu) + x*ee,2) + 
	      pow(-((ydot*xdotv)/mu) + 
		  y*ee,2) + 
	      pow(-((zdot*xdotv)/mu) + 
		  z*ee,2),1.5)) + 
	 ((-(xdot*z) + x*zdot)*
	  ((x*z)/
	   r3\
	   - (xdot*zdot)/mu) + 
	  (-(ydot*z) + y*zdot)*
	  ((y*z)/
	   r3\
	   - (ydot*zdot)/mu) - 
	  xdot*(-((xdot*xdotv)/mu) + 
		x*ee) - 
	  ydot*(-((ydot*xdotv)/mu) + 
		y*ee))/
	 (sqrt(pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	  sqrt(pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) - 
	 ((2*xdot*(xdot*z - x*zdot) - 
	   2*ydot*(-(ydot*z) + y*zdot))*
	  ((-(xdot*z) + x*zdot)*
	   (-((xdot*xdotv)/mu) + 
	    x*ee) + 
	   (-(ydot*z) + y*zdot)*
	   (-((ydot*xdotv)/mu) + 
	    y*ee)))/
	 (2.*pow(pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2),1.5)*
	  sqrt(pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))))/
	sqrt(1 - pow((-(xdot*z) + x*zdot)*
		     (-((xdot*xdotv)/mu) + 
		      x*ee) + 
		     (-(ydot*z) + y*zdot)*
		     (-((ydot*xdotv)/mu) + 
		      y*ee),2)/
	     ((pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	      (pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2)))));
 
    derivs(4,3) = 
      -((-(((-(xdot*z) + x*zdot)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) + 
	    (-(ydot*z) + y*zdot)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee))*
	   (2*((x*xdot)/mu - (xdotv)/mu)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) + 
	    2*((2*xdot*y)/mu - (x*ydot)/mu)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee) + 
	    2*((2*xdot*z)/mu - (x*zdot)/mu)*
	    (-((zdot*xdotv)/mu) + 
	     z*ee)))/
	 (2.*sqrt(pow(xdot*z - x*zdot,2) + 
		  pow(-(ydot*z) + y*zdot,2))*
	  pow(pow(-((xdot*xdotv)/
		    mu) + x*ee,2) + 
	      pow(-((ydot*xdotv)/mu) + 
		  y*ee,2) + 
	      pow(-((zdot*xdotv)/mu) + 
		  z*ee,2),1.5)) + 
	 (((2*xdot*y)/mu - (x*ydot)/mu)*(-(ydot*z) + y*zdot) + 
	  (-(xdot*z) + x*zdot)*
	  ((x*xdot)/mu - (xdotv)/mu) - 
	  z*(-((xdot*xdotv)/mu) + 
	     x*ee))/
	 (sqrt(pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	  sqrt(pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) - 
	 (z*(xdot*z - x*zdot)*
	  ((-(xdot*z) + x*zdot)*
	   (-((xdot*xdotv)/mu) + 
	    x*ee) + 
	   (-(ydot*z) + y*zdot)*
	   (-((ydot*xdotv)/mu) + 
	    y*ee)))/
	 (pow(pow(xdot*z - x*zdot,2) + 
	      pow(-(ydot*z) + y*zdot,2),1.5)*
	  sqrt(pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))))/
	sqrt(1 - pow((-(xdot*z) + x*zdot)*
		     (-((xdot*xdotv)/mu) + 
		      x*ee) + 
		     (-(ydot*z) + y*zdot)*
		     (-((ydot*xdotv)/mu) + 
		      y*ee),2)/
	     ((pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	      (pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2)))));
 
    derivs(4,4) = 
      -((-(((-(xdot*z) + x*zdot)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) + 
	    (-(ydot*z) + y*zdot)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee))*
	   (2*(-((xdot*y)/mu) + (2*x*ydot)/mu)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) + 
	    2*((y*ydot)/mu - (xdotv)/mu)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee) + 
	    2*((2*ydot*z)/mu - (y*zdot)/mu)*
	    (-((zdot*xdotv)/mu) + 
	     z*ee)))/
	 (2.*sqrt(pow(xdot*z - x*zdot,2) + 
		  pow(-(ydot*z) + y*zdot,2))*
	  pow(pow(-((xdot*xdotv)/
		    mu) + x*ee,2) + 
	      pow(-((ydot*xdotv)/mu) + 
		  y*ee,2) + 
	      pow(-((zdot*xdotv)/mu) + 
		  z*ee,2),1.5)) + 
	 ((-((xdot*y)/mu) + (2*x*ydot)/mu)*
	  (-(xdot*z) + x*zdot) + 
	  (-(ydot*z) + y*zdot)*
	  ((y*ydot)/mu - (xdotv)/mu) - 
	  z*(-((ydot*xdotv)/mu) + 
	     y*ee))/
	 (sqrt(pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	  sqrt(pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) + 
	 (z*(-(ydot*z) + y*zdot)*
	  ((-(xdot*z) + x*zdot)*
	   (-((xdot*xdotv)/mu) + 
	    x*ee) + 
	   (-(ydot*z) + y*zdot)*
	   (-((ydot*xdotv)/mu) + 
	    y*ee)))/
	 (pow(pow(xdot*z - x*zdot,2) + 
	      pow(-(ydot*z) + y*zdot,2),1.5)*
	  sqrt(pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))))/
	sqrt(1 - pow((-(xdot*z) + x*zdot)*
		     (-((xdot*xdotv)/mu) + 
		      x*ee) + 
		     (-(ydot*z) + y*zdot)*
		     (-((ydot*xdotv)/mu) + 
		      y*ee),2)/
	     ((pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	      (pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2)))));
 
    derivs(4,5) = 
      -((-(((-(xdot*z) + x*zdot)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) + 
	    (-(ydot*z) + y*zdot)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee))*
	   (2*(-((xdot*z)/mu) + (2*x*zdot)/mu)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) + 
	    2*(-((ydot*z)/mu) + (2*y*zdot)/mu)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee) + 
	    2*((z*zdot)/mu - (xdotv)/mu)*
	    (-((zdot*xdotv)/mu) + 
	     z*ee)))/
	 (2.*sqrt(pow(xdot*z - x*zdot,2) + 
		  pow(-(ydot*z) + y*zdot,2))*
	  pow(pow(-((xdot*xdotv)/
		    mu) + x*ee,2) + 
	      pow(-((ydot*xdotv)/mu) + 
		  y*ee,2) + 
	      pow(-((zdot*xdotv)/mu) + 
		  z*ee,2),1.5)) + 
	 ((-(xdot*z) + x*zdot)*
	  (-((xdot*z)/mu) + (2*x*zdot)/mu) + 
	  (-(ydot*z) + y*zdot)*
	  (-((ydot*z)/mu) + (2*y*zdot)/mu) + 
	  x*(-((xdot*xdotv)/mu) + 
	     x*ee) + 
	  y*(-((ydot*xdotv)/mu) + 
	     y*ee))/
	 (sqrt(pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	  sqrt(pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) - 
	 ((-2*x*(xdot*z - x*zdot) + 2*y*(-(ydot*z) + y*zdot))*
	  ((-(xdot*z) + x*zdot)*
	   (-((xdot*xdotv)/mu) + 
	    x*ee) + 
	   (-(ydot*z) + y*zdot)*
	   (-((ydot*xdotv)/mu) + 
	    y*ee)))/
	 (2.*pow(pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2),1.5)*
	  sqrt(pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))))/
	sqrt(1 - pow((-(xdot*z) + x*zdot)*
		     (-((xdot*xdotv)/mu) + 
		      x*ee) + 
		     (-(ydot*z) + y*zdot)*
		     (-((ydot*xdotv)/mu) + 
		      y*ee),2)/
	     ((pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	      (pow(-((xdot*xdotv)/mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2)))));
 
 
    /* partials of T (time from periapse) now */

    derivs(5,0) = 
      -((((xdotv)*
	  sqrt(pow(-(xdot*y) + x*ydot,2) + 
	       pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	  (E2)*(-2*(-(vx2/mu) + 
		    x2/
		    pow(rr,
			1.5) - 1/r + 
		    vv/mu)*
		(-((xdot*xdotv)/mu) + 
		 x*ee) - 
		2*(-((xdot*ydot)/mu) + 
		   (x*y)/
		   pow(rr,
		       1.5))*(-((ydot*xdotv)/
				mu) + y*
			      (ee)) - 
		2*((x*z)/
		   pow(rr,
		       1.5) - (xdot*zdot)/mu)*
		(-((zdot*xdotv)/mu) + 
		 z*ee)))/
	 (2.*mu*pow(1 - 
		    pow(-((xdot*xdotv)/mu) + 
			x*ee,2) - 
		    pow(-((ydot*xdotv)/mu) + 
			y*ee,2) - 
		    pow(-((zdot*xdotv)/mu) + 
			z*ee,2),1.5)) + 
	 (2*x*xdotv*
	  sqrt(pow(-(xdot*y) + x*ydot,2) + 
	       pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2)))/
	 (mu*r3*
	  sqrt(1 - pow(-((xdot*xdotv)/
			 mu) + x*
		       (ee),2) - 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) - 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) - 
	 ((xdotv)*
	  (2*ydot*(-(xdot*y) + x*ydot) - 
	   2*zdot*(xdot*z - x*zdot))*
	  (E2))/
	 (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		     pow(xdot*z - x*zdot,2) + 
		     pow(-(ydot*z) + y*zdot,2))*
	  sqrt(1 - pow(-((xdot*xdotv)/
			 mu) + x*
		       (ee),2) - 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) - 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) - 
	 (xdot*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		    pow(xdot*z - x*zdot,2) + 
		    pow(-(ydot*z) + y*zdot,2))*
	  (E2))/
	 (mu*sqrt(1 - pow(-((xdot*
			     (xdotv))/mu) + 
			  x*ee,2) - 
		  pow(-((ydot*xdotv)/mu) + 
		      y*ee,2) - 
		  pow(-((zdot*xdotv)/mu) + 
		      z*ee,2))) + 
	 (-((xdotv)*
	    sqrt(pow(-(xdot*y) + x*ydot,2) + 
		 pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2))*
	    (E2)*
	    (2*(-(vx2/mu) + 
		x2/r3 - 
		1/r
		+ vv/mu)*
	     (-((xdot*xdotv)/mu) + 
	      x*ee) + 
	     2*(-((xdot*ydot)/mu) + 
		(x*y)/
		r3)*
	     (-((ydot*xdotv)/mu) + 
	      y*ee) + 
	     2*((x*z)/
		r3 - (xdot*zdot)/mu)*
	     (-((zdot*xdotv)/mu) + 
	      z*ee)))/
	  (2.*mu*sqrt(1 - 
		      pow(-((xdot*xdotv)/
			    mu) + 
			  x*ee,2) - 
		      pow(-((ydot*xdotv)/
			    mu) + 
			  y*ee,2) - 
		      pow(-((zdot*xdotv)/
			    mu) + 
			  z*ee,2))*
	   pow(pow(-((xdot*xdotv)/
		     mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/
		     mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/
		     mu) + 
		   z*ee,2),1.5)) - 
	  ((xdotv)*
	   sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2))*
	   (E2)*
	   (-2*(-(vx2/mu) + 
		x2/
		pow(rr,
		    1.5) - 
		1/r + vv/mu)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) - 
	    2*(-((xdot*ydot)/mu) + 
	       (x*y)/
	       pow(rr,
		   1.5))*
	    (-((ydot*xdotv)/mu) + 
	     y*ee) - 
	    2*((x*z)/
	       pow(rr,
		   1.5) - (xdot*zdot)/mu)*
	    (-((zdot*xdotv)/mu) + 
	     z*ee)))/
	  (2.*mu*pow(1 - 
		     pow(-((xdot*xdotv)/
			   mu) + 
			 x*ee,2) - 
		     pow(-((ydot*xdotv)/
			   mu) + 
			 y*ee,2) - 
		     pow(-((zdot*xdotv)/
			   mu) + 
			 z*ee,2),1.5)*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) - 
	  (2*x*xdotv*
	   sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2)))/
	  (mu*pow(rr,
		  1.5)*sqrt(1 - 
			    pow(-((xdot*xdotv)/
				  mu) + 
				x*ee,2) - 
			    pow(-((ydot*xdotv)/
				  mu) + 
				y*ee,2) - 
			    pow(-((zdot*xdotv)/
				  mu) + 
				z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) + 
	  ((xdotv)*
	   (2*ydot*(-(xdot*y) + x*ydot) - 
	    2*zdot*(xdot*z - x*zdot))*
	   (E2))/
	  (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		      pow(xdot*z - x*zdot,2) + 
		      pow(-(ydot*z) + y*zdot,2))*
	   sqrt(1 - pow(-((xdot*
			   (xdotv))/mu) + 
			x*ee,2) - 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) - 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) + 
	  (xdot*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		     pow(xdot*z - x*zdot,2) + 
		     pow(-(ydot*z) + y*zdot,2))*
	   (E2))/
	  (mu*sqrt(1 - pow(-((xdot*
			      (xdotv))/mu) + 
			   x*ee,2) - 
		   pow(-((ydot*xdotv)/
			 mu) + 
		       y*ee,2) - 
		   pow(-((zdot*xdotv)/
			 mu) + 
		       z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))))/
	 sqrt(1 - (pow(xdotv,2)*
		   (pow(-(xdot*y) + x*ydot,2) + 
		    pow(xdot*z - x*zdot,2) + 
		    pow(-(ydot*z) + y*zdot,2))*
		   E2sq)/
	      (pow(mu,2)*(1 - 
			  pow(-((xdot*xdotv)/
				mu) + 
			      x*ee,2) - 
			  pow(-((ydot*xdotv)/
				mu) + 
			      y*ee,2) - 
			  pow(-((zdot*xdotv)/
				mu) + 
			      z*ee,2))*
	       (pow(-((xdot*xdotv)/mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2)))))/
	sqrt(mu*E2cube
	     )) - (3*mu*x*E2sq*
		   (-(((xdotv)*
		       sqrt(pow(-(xdot*y) + x*ydot,2) + 
			    pow(xdot*z - x*zdot,2) + 
			    pow(-(ydot*z) + y*zdot,2))*
		       (E2))/
		      (mu*sqrt(1 - pow(-((xdot*
					  (xdotv))/mu) + 
				       x*ee,2) - 
			       pow(-((ydot*xdotv)/mu) + 
				   y*ee,2) - 
			       pow(-((zdot*xdotv)/mu) + 
				   z*ee,2)))) + 
		    asin(((xdotv)*
			  sqrt(pow(-(xdot*y) + x*ydot,2) + 
			       pow(xdot*z - x*zdot,2) + 
			       pow(-(ydot*z) + y*zdot,2))*
			  (E2))/
			 (mu*sqrt(1 - pow(-((xdot*
					     (xdotv))/mu) + 
					  x*ee,2) - 
				  pow(-((ydot*xdotv)/mu) + 
				      y*ee,2) - 
				  pow(-((zdot*xdotv)/mu) + 
				      z*ee,2))*
			  sqrt(pow(-((xdot*xdotv)/
				     mu) + x*
				   (ee),2) + 
			       pow(-((ydot*xdotv)/mu) + 
				   y*ee,2) + 
			       pow(-((zdot*xdotv)/mu) + 
				   z*ee,2))))))/
      (r3*
       pow(mu*E2cube
	   ,1.5));
 
    derivs(5,1) = 
      -((((xdotv)*
	  sqrt(pow(-(xdot*y) + x*ydot,2) + 
	       pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	  (E2)*(-2*(-((xdot*ydot)/mu) + 
		    (x*y)/
		    pow(rr,
			1.5))*(-((xdot*xdotv)/
				 mu) + x*
			       (ee)) - 
		2*(-(vy2/mu) + 
		   y2/
		   pow(rr,
		       1.5) - 1/r + 
		   vv/mu)*
		(-((ydot*xdotv)/mu) + 
		 y*ee) - 
		2*((y*z)/
		   pow(rr,
		       1.5) - (ydot*zdot)/mu)*
		(-((zdot*xdotv)/mu) + 
		 z*ee)))/
	 (2.*mu*pow(1 - 
		    pow(-((xdot*xdotv)/mu) + 
			x*ee,2) - 
		    pow(-((ydot*xdotv)/mu) + 
			y*ee,2) - 
		    pow(-((zdot*xdotv)/mu) + 
			z*ee,2),1.5)) + 
	 (2*y*xdotv*
	  sqrt(pow(-(xdot*y) + x*ydot,2) + 
	       pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2)))/
	 (mu*r3*
	  sqrt(1 - pow(-((xdot*xdotv)/
			 mu) + x*
		       (ee),2) - 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) - 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) - 
	 ((xdotv)*
	  (-2*xdot*(-(xdot*y) + x*ydot) + 
	   2*zdot*(-(ydot*z) + y*zdot))*
	  (E2))/
	 (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		     pow(xdot*z - x*zdot,2) + 
		     pow(-(ydot*z) + y*zdot,2))*
	  sqrt(1 - pow(-((xdot*xdotv)/
			 mu) + x*
		       (ee),2) - 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) - 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) - 
	 (ydot*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		    pow(xdot*z - x*zdot,2) + 
		    pow(-(ydot*z) + y*zdot,2))*
	  (E2))/
	 (mu*sqrt(1 - pow(-((xdot*
			     (xdotv))/mu) + 
			  x*ee,2) - 
		  pow(-((ydot*xdotv)/mu) + 
		      y*ee,2) - 
		  pow(-((zdot*xdotv)/mu) + 
		      z*ee,2))) + 
	 (-((xdotv)*
	    sqrt(pow(-(xdot*y) + x*ydot,2) + 
		 pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2))*
	    (E2)*
	    (2*(-((xdot*ydot)/mu) + 
		(x*y)/
		r3)*
	     (-((xdot*xdotv)/mu) + 
	      x*ee) + 
	     2*(-(vy2/mu) + 
		y2/r3 - 
		1/r
		+ vv/mu)*
	     (-((ydot*xdotv)/mu) + 
	      y*ee) + 
	     2*((y*z)/
		r3 - (ydot*zdot)/mu)*
	     (-((zdot*xdotv)/mu) + 
	      z*ee)))/
	  (2.*mu*sqrt(1 - 
		      pow(-((xdot*xdotv)/
			    mu) + 
			  x*ee,2) - 
		      pow(-((ydot*xdotv)/
			    mu) + 
			  y*ee,2) - 
		      pow(-((zdot*xdotv)/
			    mu) + 
			  z*ee,2))*
	   pow(pow(-((xdot*xdotv)/
		     mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/
		     mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/
		     mu) + 
		   z*ee,2),1.5)) - 
	  ((xdotv)*
	   sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2))*
	   (E2)*
	   (-2*(-((xdot*ydot)/mu) + 
		(x*y)/
		pow(rr,
		    1.5))*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) - 
	    2*(-(vy2/mu) + 
	       y2/
	       pow(rr,
		   1.5) + ee)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee) - 
	    2*((y*z)/
	       pow(rr,
		   1.5) - (ydot*zdot)/mu)*
	    (-((zdot*xdotv)/mu) + 
	     z*ee)))/
	  (2.*mu*pow(1 - 
		     pow(-((xdot*xdotv)/
			   mu) + 
			 x*ee,2) - 
		     pow(-((ydot*xdotv)/
			   mu) + 
			 y*ee,2) - 
		     pow(-((zdot*xdotv)/
			   mu) + 
			 z*ee,2),1.5)*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) - 
	  (2*y*xdotv*
	   sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2)))/
	  (mu*pow(rr,
		  1.5)*sqrt(1 - 
			    pow(-((xdot*xdotv)/
				  mu) + 
				x*ee,2) - 
			    pow(-((ydot*xdotv)/
				  mu) + 
				y*ee,2) - 
			    pow(-((zdot*xdotv)/
				  mu) + 
				z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) + 
	  ((xdotv)*
	   (-2*xdot*(-(xdot*y) + x*ydot) + 
	    2*zdot*(-(ydot*z) + y*zdot))*
	   (E2))/
	  (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		      pow(xdot*z - x*zdot,2) + 
		      pow(-(ydot*z) + y*zdot,2))*
	   sqrt(1 - pow(-((xdot*
			   (xdotv))/mu) + 
			x*ee,2) - 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) - 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) + 
	  (ydot*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		     pow(xdot*z - x*zdot,2) + 
		     pow(-(ydot*z) + y*zdot,2))*
	   (E2))/
	  (mu*sqrt(1 - pow(-((xdot*
			      (xdotv))/mu) + 
			   x*ee,2) - 
		   pow(-((ydot*xdotv)/
			 mu) + 
		       y*ee,2) - 
		   pow(-((zdot*xdotv)/
			 mu) + 
		       z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))))/
	 sqrt(1 - (pow(xdotv,2)*
		   (pow(-(xdot*y) + x*ydot,2) + 
		    pow(xdot*z - x*zdot,2) + 
		    pow(-(ydot*z) + y*zdot,2))*
		   E2sq)/
	      (pow(mu,2)*(1 - 
			  pow(-((xdot*xdotv)/
				mu) + 
			      x*ee,2) - 
			  pow(-((ydot*xdotv)/
				mu) + 
			      y*ee,2) - 
			  pow(-((zdot*xdotv)/
				mu) + 
			      z*ee,2))*
	       (pow(-((xdot*xdotv)/mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2)))))/
	sqrt(mu*E2cube
	     )) - (3*mu*y*E2sq*
		   (-(((xdotv)*
		       sqrt(pow(-(xdot*y) + x*ydot,2) + 
			    pow(xdot*z - x*zdot,2) + 
			    pow(-(ydot*z) + y*zdot,2))*
		       (E2))/
		      (mu*sqrt(1 - pow(-((xdot*
					  (xdotv))/mu) + 
				       x*ee,2) - 
			       pow(-((ydot*xdotv)/mu) + 
				   y*ee,2) - 
			       pow(-((zdot*xdotv)/mu) + 
				   z*ee,2)))) + 
		    asin(((xdotv)*
			  sqrt(pow(-(xdot*y) + x*ydot,2) + 
			       pow(xdot*z - x*zdot,2) + 
			       pow(-(ydot*z) + y*zdot,2))*
			  (E2))/
			 (mu*sqrt(1 - pow(-((xdot*
					     (xdotv))/mu) + 
					  x*ee,2) - 
				  pow(-((ydot*xdotv)/mu) + 
				      y*ee,2) - 
				  pow(-((zdot*xdotv)/mu) + 
				      z*ee,2))*
			  sqrt(pow(-((xdot*xdotv)/
				     mu) + x*
				   (ee),2) + 
			       pow(-((ydot*xdotv)/mu) + 
				   y*ee,2) + 
			       pow(-((zdot*xdotv)/mu) + 
				   z*ee,2))))))/
      (r3*
       pow(mu*E2cube
	   ,1.5));

    derivs(5,2) = 
      -((((xdotv)*
	  sqrt(pow(-(xdot*y) + x*ydot,2) + 
	       pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	  (E2)*(-2*((x*z)/
		    pow(rr,
			1.5) - (xdot*zdot)/mu)*
		(-((xdot*xdotv)/mu) + 
		 x*ee) - 
		2*((y*z)/
		   pow(rr,
		       1.5) - (ydot*zdot)/mu)*
		(-((ydot*xdotv)/mu) + 
		 y*ee) - 
		2*(z2/
		   pow(rr,
		       1.5) - 1/r - 
		   vz2/mu + 
		   vv/mu)*
		(-((zdot*xdotv)/mu) + 
		 z*ee)))/
	 (2.*mu*pow(1 - 
		    pow(-((xdot*xdotv)/mu) + 
			x*ee,2) - 
		    pow(-((ydot*xdotv)/mu) + 
			y*ee,2) - 
		    pow(-((zdot*xdotv)/mu) + 
			z*ee,2),1.5)) + 
	 (2*z*xdotv*
	  sqrt(pow(-(xdot*y) + x*ydot,2) + 
	       pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2)))/
	 (mu*r3*
	  sqrt(1 - pow(-((xdot*xdotv)/
			 mu) + x*
		       (ee),2) - 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) - 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) - 
	 ((xdotv)*
	  (2*xdot*(xdot*z - x*zdot) - 
	   2*ydot*(-(ydot*z) + y*zdot))*
	  (E2))/
	 (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		     pow(xdot*z - x*zdot,2) + 
		     pow(-(ydot*z) + y*zdot,2))*
	  sqrt(1 - pow(-((xdot*xdotv)/
			 mu) + x*
		       (ee),2) - 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) - 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) - 
	 (zdot*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		    pow(xdot*z - x*zdot,2) + 
		    pow(-(ydot*z) + y*zdot,2))*
	  (E2))/
	 (mu*sqrt(1 - pow(-((xdot*
			     (xdotv))/mu) + 
			  x*ee,2) - 
		  pow(-((ydot*xdotv)/mu) + 
		      y*ee,2) - 
		  pow(-((zdot*xdotv)/mu) + 
		      z*ee,2))) + 
	 (-((xdotv)*
	    sqrt(pow(-(xdot*y) + x*ydot,2) + 
		 pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2))*
	    (E2)*
	    (2*((x*z)/
		r3 - (xdot*zdot)/mu)*
	     (-((xdot*xdotv)/mu) + 
	      x*ee) + 
	     2*((y*z)/
		r3 - (ydot*zdot)/mu)*
	     (-((ydot*xdotv)/mu) + 
	      y*ee) + 
	     2*(z2/r3 - 
		1/r
		- vz2/mu + 
		vv/mu)*
	     (-((zdot*xdotv)/mu) + 
	      z*ee)))/
	  (2.*mu*sqrt(1 - 
		      pow(-((xdot*xdotv)/
			    mu) + 
			  x*ee,2) - 
		      pow(-((ydot*xdotv)/
			    mu) + 
			  y*ee,2) - 
		      pow(-((zdot*xdotv)/
			    mu) + 
			  z*ee,2))*
	   pow(pow(-((xdot*xdotv)/
		     mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/
		     mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/
		     mu) + 
		   z*ee,2),1.5)) - 
	  ((xdotv)*
	   sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2))*
	   (E2)*
	   (-2*((x*z)/
		pow(rr,
		    1.5) - (xdot*zdot)/mu)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) - 
	    2*((y*z)/
	       pow(rr,
		   1.5) - (ydot*zdot)/mu)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee) - 
	    2*(z2/
	       pow(rr,
		   1.5) - 
	       1/r\
	       - vz2/mu + 
	       vv/mu)*
	    (-((zdot*xdotv)/mu) + 
	     z*ee)))/
	  (2.*mu*pow(1 - 
		     pow(-((xdot*xdotv)/
			   mu) + 
			 x*ee,2) - 
		     pow(-((ydot*xdotv)/
			   mu) + 
			 y*ee,2) - 
		     pow(-((zdot*xdotv)/
			   mu) + 
			 z*ee,2),1.5)*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) - 
	  (2*z*xdotv*
	   sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2)))/
	  (mu*pow(rr,
		  1.5)*sqrt(1 - 
			    pow(-((xdot*xdotv)/
				  mu) + 
				x*ee,2) - 
			    pow(-((ydot*xdotv)/
				  mu) + 
				y*ee,2) - 
			    pow(-((zdot*xdotv)/
				  mu) + 
				z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) + 
	  ((xdotv)*
	   (2*xdot*(xdot*z - x*zdot) - 
	    2*ydot*(-(ydot*z) + y*zdot))*
	   (E2))/
	  (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		      pow(xdot*z - x*zdot,2) + 
		      pow(-(ydot*z) + y*zdot,2))*
	   sqrt(1 - pow(-((xdot*
			   (xdotv))/mu) + 
			x*ee,2) - 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) - 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) + 
	  (zdot*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		     pow(xdot*z - x*zdot,2) + 
		     pow(-(ydot*z) + y*zdot,2))*
	   (E2))/
	  (mu*sqrt(1 - pow(-((xdot*
			      (xdotv))/mu) + 
			   x*ee,2) - 
		   pow(-((ydot*xdotv)/
			 mu) + 
		       y*ee,2) - 
		   pow(-((zdot*xdotv)/
			 mu) + 
		       z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))))/
	 sqrt(1 - (pow(xdotv,2)*
		   (pow(-(xdot*y) + x*ydot,2) + 
		    pow(xdot*z - x*zdot,2) + 
		    pow(-(ydot*z) + y*zdot,2))*
		   E2sq)/
	      (pow(mu,2)*(1 - 
			  pow(-((xdot*xdotv)/
				mu) + 
			      x*ee,2) - 
			  pow(-((ydot*xdotv)/
				mu) + 
			      y*ee,2) - 
			  pow(-((zdot*xdotv)/
				mu) + 
			      z*ee,2))*
	       (pow(-((xdot*xdotv)/mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2)))))/
	sqrt(mu*E2cube
	     )) - (3*mu*z*E2sq*
		   (-(((xdotv)*
		       sqrt(pow(-(xdot*y) + x*ydot,2) + 
			    pow(xdot*z - x*zdot,2) + 
			    pow(-(ydot*z) + y*zdot,2))*
		       (E2))/
		      (mu*sqrt(1 - pow(-((xdot*
					  (xdotv))/mu) + 
				       x*ee,2) - 
			       pow(-((ydot*xdotv)/mu) + 
				   y*ee,2) - 
			       pow(-((zdot*xdotv)/mu) + 
				   z*ee,2)))) + 
		    asin(((xdotv)*
			  sqrt(pow(-(xdot*y) + x*ydot,2) + 
			       pow(xdot*z - x*zdot,2) + 
			       pow(-(ydot*z) + y*zdot,2))*
			  (E2))/
			 (mu*sqrt(1 - pow(-((xdot*
					     (xdotv))/mu) + 
					  x*ee,2) - 
				  pow(-((ydot*xdotv)/mu) + 
				      y*ee,2) - 
				  pow(-((zdot*xdotv)/mu) + 
				      z*ee,2))*
			  sqrt(pow(-((xdot*xdotv)/
				     mu) + x*
				   (ee),2) + 
			       pow(-((ydot*xdotv)/mu) + 
				   y*ee,2) + 
			       pow(-((zdot*xdotv)/mu) + 
				   z*ee,2))))))/
      (r3*
       pow(mu*E2cube
	   ,1.5));
 
    derivs(5,3) = 
      -((((xdotv)*
	  sqrt(pow(-(xdot*y) + x*ydot,2) + 
	       pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	  (E2)*(-2*((x*xdot)/mu - 
		    (xdotv)/mu)*
		(-((xdot*xdotv)/mu) + 
		 x*ee) - 
		2*((2*xdot*y)/mu - (x*ydot)/mu)*
		(-((ydot*xdotv)/mu) + 
		 y*ee) - 
		2*((2*xdot*z)/mu - (x*zdot)/mu)*
		(-((zdot*xdotv)/mu) + 
		 z*ee)))/
	 (2.*mu*pow(1 - 
		    pow(-((xdot*xdotv)/mu) + 
			x*ee,2) - 
		    pow(-((ydot*xdotv)/mu) + 
			y*ee,2) - 
		    pow(-((zdot*xdotv)/mu) + 
			z*ee,2),1.5)) + 
	 (2*xdot*xdotv*
	  sqrt(pow(-(xdot*y) + x*ydot,2) + 
	       pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2)))/
	 (pow(mu,2)*sqrt(1 - 
			 pow(-((xdot*xdotv)/mu) + 
			     x*ee,2) - 
			 pow(-((ydot*xdotv)/mu) + 
			     y*ee,2) - 
			 pow(-((zdot*xdotv)/mu) + 
			     z*ee,2))) - 
	 ((xdotv)*
	  (-2*y*(-(xdot*y) + x*ydot) + 2*z*(xdot*z - x*zdot))*
	  (E2))/
	 (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		     pow(xdot*z - x*zdot,2) + 
		     pow(-(ydot*z) + y*zdot,2))*
	  sqrt(1 - pow(-((xdot*xdotv)/
			 mu) + x*
		       (ee),2) - 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) - 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) - 
	 (x*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		 pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2))*
	  (E2))/
	 (mu*sqrt(1 - pow(-((xdot*
			     (xdotv))/mu) + 
			  x*ee,2) - 
		  pow(-((ydot*xdotv)/mu) + 
		      y*ee,2) - 
		  pow(-((zdot*xdotv)/mu) + 
		      z*ee,2))) + 
	 (-((xdotv)*
	    sqrt(pow(-(xdot*y) + x*ydot,2) + 
		 pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2))*
	    (E2)*
	    (2*((x*xdot)/mu - 
		(xdotv)/mu)*
	     (-((xdot*xdotv)/mu) + 
	      x*ee) + 
	     2*((2*xdot*y)/mu - (x*ydot)/mu)*
	     (-((ydot*xdotv)/mu) + 
	      y*ee) + 
	     2*((2*xdot*z)/mu - (x*zdot)/mu)*
	     (-((zdot*xdotv)/mu) + 
	      z*ee)))/
	  (2.*mu*sqrt(1 - 
		      pow(-((xdot*xdotv)/
			    mu) + 
			  x*ee,2) - 
		      pow(-((ydot*xdotv)/
			    mu) + 
			  y*ee,2) - 
		      pow(-((zdot*xdotv)/
			    mu) + 
			  z*ee,2))*
	   pow(pow(-((xdot*xdotv)/
		     mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/
		     mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/
		     mu) + 
		   z*ee,2),1.5)) - 
	  ((xdotv)*
	   sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2))*
	   (E2)*
	   (-2*((x*xdot)/mu - 
		(xdotv)/mu)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) - 
	    2*((2*xdot*y)/mu - (x*ydot)/mu)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee) - 
	    2*((2*xdot*z)/mu - (x*zdot)/mu)*
	    (-((zdot*xdotv)/mu) + 
	     z*ee)))/
	  (2.*mu*pow(1 - 
		     pow(-((xdot*xdotv)/
			   mu) + 
			 x*ee,2) - 
		     pow(-((ydot*xdotv)/
			   mu) + 
			 y*ee,2) - 
		     pow(-((zdot*xdotv)/
			   mu) + 
			 z*ee,2),1.5)*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) - 
	  (2*xdot*xdotv*
	   sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2)))/
	  (pow(mu,2)*sqrt(1 - 
			  pow(-((xdot*xdotv)/
				mu) + 
			      x*ee,2) - 
			  pow(-((ydot*xdotv)/
				mu) + 
			      y*ee,2) - 
			  pow(-((zdot*xdotv)/
				mu) + 
			      z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) + 
	  ((xdotv)*
	   (-2*y*(-(xdot*y) + x*ydot) + 
	    2*z*(xdot*z - x*zdot))*
	   (E2))/
	  (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		      pow(xdot*z - x*zdot,2) + 
		      pow(-(ydot*z) + y*zdot,2))*
	   sqrt(1 - pow(-((xdot*
			   (xdotv))/mu) + 
			x*ee,2) - 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) - 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) + 
	  (x*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		  pow(xdot*z - x*zdot,2) + 
		  pow(-(ydot*z) + y*zdot,2))*
	   (E2))/
	  (mu*sqrt(1 - pow(-((xdot*
			      (xdotv))/mu) + 
			   x*ee,2) - 
		   pow(-((ydot*xdotv)/
			 mu) + 
		       y*ee,2) - 
		   pow(-((zdot*xdotv)/
			 mu) + 
		       z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))))/
	 sqrt(1 - (pow(xdotv,2)*
		   (pow(-(xdot*y) + x*ydot,2) + 
		    pow(xdot*z - x*zdot,2) + 
		    pow(-(ydot*z) + y*zdot,2))*
		   E2sq)/
	      (pow(mu,2)*(1 - 
			  pow(-((xdot*xdotv)/
				mu) + 
			      x*ee,2) - 
			  pow(-((ydot*xdotv)/
				mu) + 
			      y*ee,2) - 
			  pow(-((zdot*xdotv)/
				mu) + 
			      z*ee,2))*
	       (pow(-((xdot*xdotv)/mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2)))))/
	sqrt(mu*E2cube
	     )) - (3*xdot*E2sq*
		   (-(((xdotv)*
		       sqrt(pow(-(xdot*y) + x*ydot,2) + 
			    pow(xdot*z - x*zdot,2) + 
			    pow(-(ydot*z) + y*zdot,2))*
		       (E2))/
		      (mu*sqrt(1 - pow(-((xdot*
					  (xdotv))/mu) + 
				       x*ee,2) - 
			       pow(-((ydot*xdotv)/mu) + 
				   y*ee,2) - 
			       pow(-((zdot*xdotv)/mu) + 
				   z*ee,2)))) + 
		    asin(((xdotv)*
			  sqrt(pow(-(xdot*y) + x*ydot,2) + 
			       pow(xdot*z - x*zdot,2) + 
			       pow(-(ydot*z) + y*zdot,2))*
			  (E2))/
			 (mu*sqrt(1 - pow(-((xdot*
					     (xdotv))/mu) + 
					  x*ee,2) - 
				  pow(-((ydot*xdotv)/mu) + 
				      y*ee,2) - 
				  pow(-((zdot*xdotv)/mu) + 
				      z*ee,2))*
			  sqrt(pow(-((xdot*xdotv)/
				     mu) + x*
				   (ee),2) + 
			       pow(-((ydot*xdotv)/mu) + 
				   y*ee,2) + 
			       pow(-((zdot*xdotv)/mu) + 
				   z*ee,2))))))/
      pow(mu*E2cube,
	  1.5);
 
 
    derivs(5,4) = 
      -((((xdotv)*
	  sqrt(pow(-(xdot*y) + x*ydot,2) + 
	       pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	  (E2)*(-2*(-((xdot*y)/mu) + (2*x*ydot)/mu)*
		(-((xdot*xdotv)/mu) + 
		 x*ee) - 
		2*((y*ydot)/mu - (xdotv)/mu)*
		(-((ydot*xdotv)/mu) + 
		 y*ee) - 
		2*((2*ydot*z)/mu - (y*zdot)/mu)*
		(-((zdot*xdotv)/mu) + 
		 z*ee)))/
	 (2.*mu*pow(1 - 
		    pow(-((xdot*xdotv)/mu) + 
			x*ee,2) - 
		    pow(-((ydot*xdotv)/mu) + 
			y*ee,2) - 
		    pow(-((zdot*xdotv)/mu) + 
			z*ee,2),1.5)) + 
	 (2*ydot*xdotv*
	  sqrt(pow(-(xdot*y) + x*ydot,2) + 
	       pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2)))/
	 (pow(mu,2)*sqrt(1 - 
			 pow(-((xdot*xdotv)/mu) + 
			     x*ee,2) - 
			 pow(-((ydot*xdotv)/mu) + 
			     y*ee,2) - 
			 pow(-((zdot*xdotv)/mu) + 
			     z*ee,2))) - 
	 ((xdotv)*
	  (2*x*(-(xdot*y) + x*ydot) - 
	   2*z*(-(ydot*z) + y*zdot))*
	  (E2))/
	 (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		     pow(xdot*z - x*zdot,2) + 
		     pow(-(ydot*z) + y*zdot,2))*
	  sqrt(1 - pow(-((xdot*xdotv)/
			 mu) + x*
		       (ee),2) - 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) - 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) - 
	 (y*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		 pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2))*
	  (E2))/
	 (mu*sqrt(1 - pow(-((xdot*
			     (xdotv))/mu) + 
			  x*ee,2) - 
		  pow(-((ydot*xdotv)/mu) + 
		      y*ee,2) - 
		  pow(-((zdot*xdotv)/mu) + 
		      z*ee,2))) + 
	 (-((xdotv)*
	    sqrt(pow(-(xdot*y) + x*ydot,2) + 
		 pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2))*
	    (E2)*
	    (2*(-((xdot*y)/mu) + (2*x*ydot)/mu)*
	     (-((xdot*xdotv)/mu) + 
	      x*ee) + 
	     2*((y*ydot)/mu - 
		(xdotv)/mu)*
	     (-((ydot*xdotv)/mu) + 
	      y*ee) + 
	     2*((2*ydot*z)/mu - (y*zdot)/mu)*
	     (-((zdot*xdotv)/mu) + 
	      z*ee)))/
	  (2.*mu*sqrt(1 - 
		      pow(-((xdot*xdotv)/
			    mu) + 
			  x*ee,2) - 
		      pow(-((ydot*xdotv)/
			    mu) + 
			  y*ee,2) - 
		      pow(-((zdot*xdotv)/
			    mu) + 
			  z*ee,2))*
	   pow(pow(-((xdot*xdotv)/
		     mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/
		     mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/
		     mu) + 
		   z*ee,2),1.5)) - 
	  ((xdotv)*
	   sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2))*
	   (E2)*
	   (-2*(-((xdot*y)/mu) + (2*x*ydot)/mu)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) - 
	    2*((y*ydot)/mu - 
	       (xdotv)/mu)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee) - 
	    2*((2*ydot*z)/mu - (y*zdot)/mu)*
	    (-((zdot*xdotv)/mu) + 
	     z*ee)))/
	  (2.*mu*pow(1 - 
		     pow(-((xdot*xdotv)/
			   mu) + 
			 x*ee,2) - 
		     pow(-((ydot*xdotv)/
			   mu) + 
			 y*ee,2) - 
		     pow(-((zdot*xdotv)/
			   mu) + 
			 z*ee,2),1.5)*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) - 
	  (2*ydot*xdotv*
	   sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2)))/
	  (pow(mu,2)*sqrt(1 - 
			  pow(-((xdot*xdotv)/
				mu) + 
			      x*ee,2) - 
			  pow(-((ydot*xdotv)/
				mu) + 
			      y*ee,2) - 
			  pow(-((zdot*xdotv)/
				mu) + 
			      z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) + 
	  ((xdotv)*
	   (2*x*(-(xdot*y) + x*ydot) - 
	    2*z*(-(ydot*z) + y*zdot))*
	   (E2))/
	  (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		      pow(xdot*z - x*zdot,2) + 
		      pow(-(ydot*z) + y*zdot,2))*
	   sqrt(1 - pow(-((xdot*
			   (xdotv))/mu) + 
			x*ee,2) - 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) - 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) + 
	  (y*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		  pow(xdot*z - x*zdot,2) + 
		  pow(-(ydot*z) + y*zdot,2))*
	   (E2))/
	  (mu*sqrt(1 - pow(-((xdot*
			      (xdotv))/mu) + 
			   x*ee,2) - 
		   pow(-((ydot*xdotv)/
			 mu) + 
		       y*ee,2) - 
		   pow(-((zdot*xdotv)/
			 mu) + 
		       z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))))/
	 sqrt(1 - (pow(xdotv,2)*
		   (pow(-(xdot*y) + x*ydot,2) + 
		    pow(xdot*z - x*zdot,2) + 
		    pow(-(ydot*z) + y*zdot,2))*
		   E2sq)/
	      (pow(mu,2)*(1 - 
			  pow(-((xdot*xdotv)/
				mu) + 
			      x*ee,2) - 
			  pow(-((ydot*xdotv)/
				mu) + 
			      y*ee,2) - 
			  pow(-((zdot*xdotv)/
				mu) + 
			      z*ee,2))*
	       (pow(-((xdot*xdotv)/mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2)))))/
	sqrt(mu*E2cube
	     )) - (3*ydot*E2sq*
		   (-(((xdotv)*
		       sqrt(pow(-(xdot*y) + x*ydot,2) + 
			    pow(xdot*z - x*zdot,2) + 
			    pow(-(ydot*z) + y*zdot,2))*
		       (E2))/
		      (mu*sqrt(1 - pow(-((xdot*
					  (xdotv))/mu) + 
				       x*ee,2) - 
			       pow(-((ydot*xdotv)/mu) + 
				   y*ee,2) - 
			       pow(-((zdot*xdotv)/mu) + 
				   z*ee,2)))) + 
		    asin(((xdotv)*
			  sqrt(pow(-(xdot*y) + x*ydot,2) + 
			       pow(xdot*z - x*zdot,2) + 
			       pow(-(ydot*z) + y*zdot,2))*
			  (E2))/
			 (mu*sqrt(1 - pow(-((xdot*
					     (xdotv))/mu) + 
					  x*ee,2) - 
				  pow(-((ydot*xdotv)/mu) + 
				      y*ee,2) - 
				  pow(-((zdot*xdotv)/mu) + 
				      z*ee,2))*
			  sqrt(pow(-((xdot*xdotv)/
				     mu) + x*
				   (ee),2) + 
			       pow(-((ydot*xdotv)/mu) + 
				   y*ee,2) + 
			       pow(-((zdot*xdotv)/mu) + 
				   z*ee,2))))))/
      pow(mu*E2cube,
	  1.5);
 
    derivs(5,5) = 
      -((((xdotv)*
	  sqrt(pow(-(xdot*y) + x*ydot,2) + 
	       pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2))*
	  (E2)*(-2*(-((xdot*z)/mu) + (2*x*zdot)/mu)*
		(-((xdot*xdotv)/mu) + 
		 x*ee) - 
		2*(-((ydot*z)/mu) + (2*y*zdot)/mu)*
		(-((ydot*xdotv)/mu) + 
		 y*ee) - 
		2*((z*zdot)/mu - (xdotv)/mu)*
		(-((zdot*xdotv)/mu) + 
		 z*ee)))/
	 (2.*mu*pow(1 - 
		    pow(-((xdot*xdotv)/mu) + 
			x*ee,2) - 
		    pow(-((ydot*xdotv)/mu) + 
			y*ee,2) - 
		    pow(-((zdot*xdotv)/mu) + 
			z*ee,2),1.5)) + 
	 (2*zdot*xdotv*
	  sqrt(pow(-(xdot*y) + x*ydot,2) + 
	       pow(xdot*z - x*zdot,2) + 
	       pow(-(ydot*z) + y*zdot,2)))/
	 (pow(mu,2)*sqrt(1 - 
			 pow(-((xdot*xdotv)/mu) + 
			     x*ee,2) - 
			 pow(-((ydot*xdotv)/mu) + 
			     y*ee,2) - 
			 pow(-((zdot*xdotv)/mu) + 
			     z*ee,2))) - 
	 ((xdotv)*
	  (-2*x*(xdot*z - x*zdot) + 2*y*(-(ydot*z) + y*zdot))*
	  (E2))/
	 (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		     pow(xdot*z - x*zdot,2) + 
		     pow(-(ydot*z) + y*zdot,2))*
	  sqrt(1 - pow(-((xdot*xdotv)/
			 mu) + x*
		       (ee),2) - 
	       pow(-((ydot*xdotv)/mu) + 
		   y*ee,2) - 
	       pow(-((zdot*xdotv)/mu) + 
		   z*ee,2))) - 
	 (z*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		 pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2))*
	  (E2))/
	 (mu*sqrt(1 - pow(-((xdot*
			     (xdotv))/mu) + 
			  x*ee,2) - 
		  pow(-((ydot*xdotv)/mu) + 
		      y*ee,2) - 
		  pow(-((zdot*xdotv)/mu) + 
		      z*ee,2))) + 
	 (-((xdotv)*
	    sqrt(pow(-(xdot*y) + x*ydot,2) + 
		 pow(xdot*z - x*zdot,2) + 
		 pow(-(ydot*z) + y*zdot,2))*
	    (E2)*
	    (2*(-((xdot*z)/mu) + (2*x*zdot)/mu)*
	     (-((xdot*xdotv)/mu) + 
	      x*ee) + 
	     2*(-((ydot*z)/mu) + (2*y*zdot)/mu)*
	     (-((ydot*xdotv)/mu) + 
	      y*ee) + 
	     2*((z*zdot)/mu - 
		(xdotv)/mu)*
	     (-((zdot*xdotv)/mu) + 
	      z*ee)))/
	  (2.*mu*sqrt(1 - 
		      pow(-((xdot*xdotv)/
			    mu) + 
			  x*ee,2) - 
		      pow(-((ydot*xdotv)/
			    mu) + 
			  y*ee,2) - 
		      pow(-((zdot*xdotv)/
			    mu) + 
			  z*ee,2))*
	   pow(pow(-((xdot*xdotv)/
		     mu) + 
		   x*ee,2) + 
	       pow(-((ydot*xdotv)/
		     mu) + 
		   y*ee,2) + 
	       pow(-((zdot*xdotv)/
		     mu) + 
		   z*ee,2),1.5)) - 
	  ((xdotv)*
	   sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2))*
	   (E2)*
	   (-2*(-((xdot*z)/mu) + (2*x*zdot)/mu)*
	    (-((xdot*xdotv)/mu) + 
	     x*ee) - 
	    2*(-((ydot*z)/mu) + (2*y*zdot)/mu)*
	    (-((ydot*xdotv)/mu) + 
	     y*ee) - 
	    2*((z*zdot)/mu - 
	       (xdotv)/mu)*
	    (-((zdot*xdotv)/mu) + 
	     z*ee)))/
	  (2.*mu*pow(1 - 
		     pow(-((xdot*xdotv)/
			   mu) + 
			 x*ee,2) - 
		     pow(-((ydot*xdotv)/
			   mu) + 
			 y*ee,2) - 
		     pow(-((zdot*xdotv)/
			   mu) + 
			 z*ee,2),1.5)*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) - 
	  (2*zdot*xdotv*
	   sqrt(pow(-(xdot*y) + x*ydot,2) + 
		pow(xdot*z - x*zdot,2) + 
		pow(-(ydot*z) + y*zdot,2)))/
	  (pow(mu,2)*sqrt(1 - 
			  pow(-((xdot*xdotv)/
				mu) + 
			      x*ee,2) - 
			  pow(-((ydot*xdotv)/
				mu) + 
			      y*ee,2) - 
			  pow(-((zdot*xdotv)/
				mu) + 
			      z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) + 
	  ((xdotv)*
	   (-2*x*(xdot*z - x*zdot) + 
	    2*y*(-(ydot*z) + y*zdot))*
	   (E2))/
	  (2.*mu*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		      pow(xdot*z - x*zdot,2) + 
		      pow(-(ydot*z) + y*zdot,2))*
	   sqrt(1 - pow(-((xdot*
			   (xdotv))/mu) + 
			x*ee,2) - 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) - 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))) + 
	  (z*sqrt(pow(-(xdot*y) + x*ydot,2) + 
		  pow(xdot*z - x*zdot,2) + 
		  pow(-(ydot*z) + y*zdot,2))*
	   (E2))/
	  (mu*sqrt(1 - pow(-((xdot*
			      (xdotv))/mu) + 
			   x*ee,2) - 
		   pow(-((ydot*xdotv)/
			 mu) + 
		       y*ee,2) - 
		   pow(-((zdot*xdotv)/
			 mu) + 
		       z*ee,2))*
	   sqrt(pow(-((xdot*xdotv)/
		      mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2))))/
	 sqrt(1 - (pow(xdotv,2)*
		   (pow(-(xdot*y) + x*ydot,2) + 
		    pow(xdot*z - x*zdot,2) + 
		    pow(-(ydot*z) + y*zdot,2))*
		   E2sq)/
	      (pow(mu,2)*(1 - 
			  pow(-((xdot*xdotv)/
				mu) + 
			      x*ee,2) - 
			  pow(-((ydot*xdotv)/
				mu) + 
			      y*ee,2) - 
			  pow(-((zdot*xdotv)/
				mu) + 
			      z*ee,2))*
	       (pow(-((xdot*xdotv)/mu) + 
		    x*ee,2) + 
		pow(-((ydot*xdotv)/
		      mu) + 
		    y*ee,2) + 
		pow(-((zdot*xdotv)/
		      mu) + 
		    z*ee,2)))))/
	sqrt(mu*E2cube
	     )) - (3*zdot*E2sq*
		   (-(((xdotv)*
		       sqrt(pow(-(xdot*y) + x*ydot,2) + 
			    pow(xdot*z - x*zdot,2) + 
			    pow(-(ydot*z) + y*zdot,2))*
		       E2)/
		      (mu*sqrt(1 - pow(-((xdot*
					  (xdotv))/mu) + 
				       x*ee,2) - 
			       pow(-((ydot*xdotv)/mu) + 
				   y*ee,2) - 
			       pow(-((zdot*xdotv)/mu) + 
				   z*ee,2)))) + 
		    asin(((xdotv)*
			  sqrt(pow(-(xdot*y) + x*ydot,2) + 
			       pow(xdot*z - x*zdot,2) + 
			       pow(-(ydot*z) + y*zdot,2))*
			  (E2))/
			 (mu*sqrt(1 - pow(-((xdot*
					     (xdotv))/mu) + 
					  x*ee,2) - 
				  pow(-((ydot*xdotv)/mu) + 
				      y*ee,2) - 
				  pow(-((zdot*xdotv)/mu) + 
				      z*ee,2))*
			  sqrt(pow(-((xdot*xdotv)/
				     mu) + x*
				   (ee),2) + 
			       pow(-((ydot*xdotv)/mu) + 
				   y*ee,2) + 
			       pow(-((zdot*xdotv)/mu) + 
				   z*ee,2))))))/
      pow(mu*E2cube,
	  1.5);
 
    // Need to include derivatives between ecliptic and ICRS states!
    Matrix66 eclDerivs(0.);
    eclDerivs.subMatrix(0,3,0,3) = partials;
    eclDerivs.subMatrix(3,6,3,6) = partials;

    // Appears that the derivatives have a sign error when the angles are >180
    orbits::Elements el = orbits::getElements(xv);
    if (el[Elements::LAN] >=PI) {
      // Flip signs of these derivatives
      for (int j=0; j<6; j++) 
	derivs(Elements::LAN,j) *= -1.;
    }
    if (el[Elements::AOP] >=PI) {
      // Flip signs of these derivatives
      for (int j=0; j<6; j++) 
	derivs(Elements::AOP,j) *= -1.;
    }
    
    return derivs * eclDerivs;
  }

} //namespace ephemeris
