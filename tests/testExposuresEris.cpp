// Test interaction of Fitter and Exposure list using Eris data
#include "Fitter.h"
#include "Elements.h"
#include "StringStuff.h"
#include "Astrometry.h"
#include "Exposures.h"

#include <iostream>

using namespace std;
using namespace orbits;

const string horizons =
  "2455197.500000000  24.3185230  -4.6068506 2010-Jan-01 00:00 \n"
  "2455287.500000000  24.7241933  -4.2193035 2010-Apr-01 00:00 \n"
  "2455378.500000000  25.4398624  -3.9956041 2010-Jul-01 00:00 \n"
  "2455470.500000000  25.1221594  -4.2500940 2010-Oct-01 00:00 \n"
  "2455562.500000000  24.4455242  -4.3499116 2011-Jan-01 00:00 \n"
  "2455652.500000000  24.8469113  -3.9656684 2011-Apr-01 00:00 \n"
  "2455743.500000000  25.5653737  -3.7384094 2011-Jul-01 00:00 \n"
  "2455835.500000000  25.2523014  -3.9891975 2011-Oct-01 00:00 \n"
  "2455927.500000000  24.5728771  -4.0926595 2012-Jan-01 00:00 \n"
  "2456018.500000000  24.9792191  -3.7070478 2012-Apr-01 00:00 \n"
  "2456109.500000000  25.6942614  -3.4815958 2012-Jul-01 00:00 \n"
  "2456201.500000000  25.3736579  -3.7316665 2012-Oct-01 00:00 \n"
  "2456293.500000000  24.6976939  -3.8331897 2013-Jan-01 00:00 \n"
  "2456383.500000000  25.1018474  -3.4531670 2013-Apr-01 00:00 \n"
  "2456474.500000000  25.8195634  -3.2241968 2013-Jul-01 00:00 \n"
  "2456566.500000000  25.5034309  -3.4706100 2013-Oct-01 00:00 \n"
  "2456658.500000000  24.8243616  -3.5758303 2014-Jan-01 00:00 \n"  //** reference epoch
  "2456748.500000000  25.2240337  -3.1992256 2014-Apr-01 00:00 \n"
  "2456839.500000000  25.9443482  -2.9668135 2014-Jul-01 00:00 \n"
  "2456931.500000000  25.6325330  -3.2095212 2014-Oct-01 00:00 \n"
  "2457023.500000000  24.9504765  -3.3183358 2015-Jan-01 00:00 \n"
  "2457113.500000000  25.3456162  -2.9451975 2015-Apr-01 00:00 \n"
  "2457204.500000000  26.0684795  -2.7093391 2015-Jul-01 00:00 \n"
  "2457296.500000000  25.7611459  -2.9482341 2015-Oct-01 00:00 \n"
  "2457388.500000000  25.0760324  -3.0606177 2016-Jan-01 00:00 \n"
  "2457479.500000000  25.4760508  -2.6861656 2016-Apr-01 00:00 \n"
  "2457570.500000000  26.1956671  -2.4519383 2016-Jul-01 00:00 \n"
  "2457662.500000000  25.8805935  -2.6901367 2016-Oct-01 00:00 \n"
  "2457754.500000000  25.1987255  -2.8006069 2017-Jan-01 00:00 \n"
  "2457844.500000000  25.5968593  -2.4315435 2017-Apr-01 00:00 \n"
  "2457935.500000000  26.3192351  -2.1938270 2017-Jul-01 00:00 \n"
  "2458027.500000000  26.0086409  -2.4281391 2017-Oct-01 00:00 \n";
/*
  "2458119.500000000  25.3239165  -2.5420374 2018-Jan-01 00:00 \n"
  "2458209.500000000  25.7177163  -2.1764158 2018-Apr-01 00:00 \n"
  "2458300.500000000  26.4429386  -1.9351945 2018-Jul-01 00:00 \n"
  "2458392.500000000  26.1370448  -2.1655122 2018-Oct-01 00:00 \n"
  "2458484.500000000  25.4494369  -2.2828180 2019-Jan-01 00:00 \n";
*/

int main(int argc,
	 char *argv[])
{
  try {
    Ephemeris ephem;

    Fitter fit(ephem, Gravity::GIANTS);
    fit.setBindingConstraint(1.);

    {
      std::istringstream iss(horizons);
      string buffer;
      while (stringstuff::getlineNoComment(iss, buffer)) {
	  std::istringstream line(buffer);
	  double jd, ra, dec;
	  line >> jd >> ra >> dec;
	  Observation obs;
	  obs.radec = astrometry::SphericalICRS(ra*DEGREE, dec*DEGREE);
	  obs.tdb = ephem.jd2tdb(jd);
	  obs.observer = ephem.observatory(807, obs.tdb);
	  obs.cov(0,0) = obs.cov(1,1) = pow(0.1*ARCSEC,2.);
	  obs.cov(1,0) = obs.cov(0,1) = 0.;
	  fit.addObservation(obs);
	}
    }
	  
    fit.chooseFrame(16);
    fit.setLinearOrbit();
    fit.newtonFit();
    cerr << "Chisq of Eris fit: " << fit.getChisq() << endl;

    Frame::writeHeader(cerr);
    Frame frame = fit.getFrame();
    frame.write(cerr);

    auto abg = fit.getABG();
    ABG::writeHeader(cerr);
    abg.write(cerr);
      
    // Calculate orbital elements:
    cerr << "Elements: " << endl;
    Elements::writeHeader(cerr);
    cerr << fit.getElements() << endl;
    {
      auto cov =  fit.getElementCovariance();
      cov.write(cerr);
    }

    // Find exposures that potentially contain Eris
    double gamma0 = 1./90.;
    double dGamma = gamma0*0.1;
    auto possibleExposures = orbits::selectExposures(frame, ephem,
						     gamma0, dGamma, 0.5);


    cerr << "Possible exposures: " << possibleExposures.size() << endl;
    // Now make a tree
    DESTree tree(possibleExposures, ephem, frame, gamma0);

    cerr << "Tree has " << tree.countNodes() << " nodes. Start search" << endl;
    // Search the tree for Eris hits
    auto hits = tree.find(fit);
			      
    {
      // Check all input exposures to see if orbit is in them, note if any are missing from possible list.
      int n = possibleExposures.size();
      DVector tobs(n);
      DVector x(n);
      DVector y(n);
      DMatrix earth(n,3);
      for (int i=0; i<n; i++) {
	tobs[i] = possibleExposures[i]->tdb - frame.tdb0;
	earth.row(i) = possibleExposures[i]->earth.transpose();
      }
      fit.predict(tobs,earth,&x,&y);
	
      int nCross=0;
      int nMiss=0;
      for (int i=0; i<n; i++) {
	const Exposure& expo = *possibleExposures[i];
	double ra,dec;
	astrometry::Gnomonic gn(x[i],y[i],frame.orient);
	astrometry::SphericalICRS(gn).getLonLat(ra,dec);
	bool hit=false;
	for (auto e : hits)
	  if (e->expnum == expo.expnum) {
	    hit = true;
	    break;
	  }
	double rad = hypot(expo.axis[0]-x[i],
			   expo.axis[1]-y[i]);
	if (rad >1.1*DEGREE) continue;
	nCross++;
	if (!hit) {
	  nMiss++;
	  cout << "MISSED " 
	       << expo.expnum << " " << tobs[i]
	       << " " << expo.axis[0]/DEGREE << " " << expo.axis[1]/DEGREE
	       << " Eris " << x[i]/DEGREE << " " << y[i]/DEGREE
	       << endl;
	}
      }
      cout << "True hits: " << nCross << " of " << hits.size()
	   << ", " << nMiss << " misses" << endl;
    }

    // Examine transient lists from each possible exposure
    {
      // Check all input exposures to see if orbit is in them, note if any are missing from possible list.
      int n = possibleExposures.size();
      DVector tobs(n);
      DVector x(n);
      DVector y(n);
      DVector covxx(n);
      DVector covxy(n);
      DVector covyy(n);
      DMatrix earth(n,3);
      for (int i=0; i<n; i++) {
	tobs[i] = possibleExposures[i]->tdb - frame.tdb0;
	earth.row(i) = possibleExposures[i]->earth.transpose();
      }
      fit.predict(tobs,earth,&x,&y,&covxx,&covxy,&covyy);
      int totalCounts = 0;
      //**double CHISQ_THRESHOLD=9.21; // 99% point of chisq
      double CHISQ_THRESHOLD=50; //** ?? try more
      for (int i=0; i<n; i++) {
	DVector chisq = possibleExposures[i]->chisq(x[i], y[i], covxx[i], covyy[i], covxy[i]);
	int counts = (chisq.array() < CHISQ_THRESHOLD).count();
	totalCounts += counts;
	if (counts>0) {
	  cout << possibleExposures[i]->expnum << " detects: ";
	  for (int j=0; j<chisq.size(); j++) {
	    if (chisq[j] < CHISQ_THRESHOLD) {
	      cout << possibleExposures[i]->id[j]
		   << " " << chisq[j]
		   << " rate " << sqrt(covxx[i]*covyy[i]-covxy[i]*covxy[i])
		                       *possibleExposures[i]->detectionDensity*CHISQ_THRESHOLD; 
	    }
	  }
	  cout << endl;
	}
      }
      cout << "Total detections: " << totalCounts << endl;
    }
    
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}

