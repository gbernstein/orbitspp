// Various conversions for OrbitTypes

#include "OrbitTypes.h"
#include "StringStuff.h"

using namespace orbits;

astrometry::DMatrix
ABG::getXYZ(const astrometry::DVector& t) const {
  // Calculate inertial XYZ positions at an array of times, return n x 3 matrix
  astrometry::DMatrix xyz(t.size(), 3, 0.);
  xyz.col(0).array() = (*this)[ADOT]*t.array() + (*this)[A];
  xyz.col(1).array() = (*this)[BDOT]*t.array() + (*this)[B];
  xyz.col(2).array() = (*this)[GDOT]*t.array() + 1.;
  xyz *= 1./(*this)[G];
  return xyz;
}

Matrix66
ABG::getStateDerivatives() const {
  Matrix66 out(0.);
  double invG = 1. / (*this)[G];
  out(0, A) = invG;
  out(0, G) = -(*this)[A]*invG*invG;
  out(1, B) = invG;
  out(1, G) = -(*this)[B]*invG*invG;
  out(2, G) = -invG*invG;
  out(3, ADOT) = invG;
  out(3, G) = -(*this)[ADOT]*invG*invG;
  out(4, BDOT) = invG;
  out(4, G) = -(*this)[BDOT]*invG*invG;
  out(5, GDOT) = invG;
  out(5, G) = -(*this)[GDOT]*invG*invG;
  return out;
}
  
void
ABG::writeTo(std::ostream& os, int precision) const {
  // Write ABG on one line
  stringstuff::StreamSaver ss(os);  // Cache stream state
  os << std::fixed << std::showpos << std::setprecision(precision);
  for (int i=0; i<6; i++)
    os << (*this)[i] << " ";
  return;  // Stream state restored on destruction of ss
}

astrometry::Vector3
Frame::toICRS(const astrometry::Vector3& x,
	      bool isVelocity) const {
  astrometry::Vector3 out = orient.m().transpose() * x;
  if (!isVelocity) {
    out += origin.getVector();
  }
  return out;
}

astrometry::Vector3
Frame::fromICRS(const astrometry::Vector3& x,
		bool isVelocity) const {
  if (isVelocity)
    return orient.m() * x;
  else
    return orient.m() * (x - origin.getVector());
}

astrometry::DMatrix
Frame::toICRS(const astrometry::DMatrix& x,
	      bool isVelocity) const {
  astrometry::DMatrix out = orient.m().transpose() * x;
  if (!isVelocity)
    out.colwise() += origin.getVector();
  return out;
}

astrometry::DMatrix
Frame::fromICRS(const astrometry::DMatrix& x,
		bool isVelocity) const {
  if (isVelocity) {
    return orient.m() * x;
  } else {
    astrometry::DMatrix tmp = x.colwise() - origin.getVector();
    return orient.m() * tmp;
  }
}

MPCObservation::MPCObservation(const string& line) {
  auto fields = stringstuff::split(line);
  auto wordPtr = fields.begin();
  if (fields.size() < 5) {
    // No way this is enough info
    throw std::runtime_error("Insufficient info on MPCObservation: " + line);
  }
  {
    // First field is indicator of time
    double mjd_in;
    bool fail = false;
    string word = *wordPtr;
    istringstream iss(word);
    iss >> mjd_in;
    if (!iss.eof() || iss.fail()) {
      throw std::runtime_error("Bad date/time in MPCObservation: " + line);
    }

    if (mjd_in < 10000.) {
      // This is probably a year.
      // Try YYYY MM DD.DDDD
      int yr = static_cast<int> (floor(mjd_in));
      int mo;
      float day;
      if (abs(mjd_in-yr) > 1e-7) {
	// Was not integral year
	throw std::runtime_error("Bad date/time in MPCObservation: " + line);
      }
      ++wordPtr;
      string moWord = *wordPtr;
      istringstream iss1(moWord);
      iss1 >> mo;
      if (!iss1.eof() || iss1.fail() || mo<1 || mo>12) {
	throw std::runtime_error("Bad date/time in MPCObservation: " + line);
      }

      ++wordPtr;
      string dayWord = *wordPtr;
      istringstream iss2(dayWord);
      iss2 >> day;
      if (!iss2.eof() || iss2.fail() || day<1.0 || day >=32.) {
	throw std::runtime_error("Bad date/time in MPCObservation: " + line);
      }

      // Steal skycalc algorithm for ymd -> jd.  Only valid 1900-2100
      if  ( yr<=1900 || yr>=2100) {
	throw runtime_error("Year out of 1900-2100 range for MJDObservation: " + line);
      }
      if(mo <= 2) {
	yr--;
	mo += 13;
      } else {
	mo += 1;
      }

      double jdint = floor(365.25*yr);  /* truncates */
      double inter = floor(30.6001*mo);
      mjd = jdint + inter + day + 1720982 - 0.5 - MJD0;
      
    } else if (mjd_in < 300000) {
      // It's an MJD
      mjd = mjd_in;
    } else {
      // It's a JD
      mjd = mjd_in - MJD0;
    }
    ++wordPtr;
  }
  // Now get RA, Dec, sigma, obscode
  if (wordPtr==fields.end())
    throw std::runtime_error("Insufficient info on MPCObservation: " + line);
  string raWord = *wordPtr;

  ++wordPtr;
  if (wordPtr==fields.end())
    throw std::runtime_error("Insufficient info on MPCObservation: " + line);
  string decWord = *wordPtr;

  ++wordPtr;
  if (wordPtr==fields.end())
    throw std::runtime_error("Insufficient info on MPCObservation: " + line);
  string sigWord = *wordPtr;

  ++wordPtr;
  if (wordPtr==fields.end())
    throw std::runtime_error("Insufficient info on MPCObservation: " + line);
  string obscodeWord = *wordPtr;

  {
    istringstream iss( raWord + " " + decWord);
    iss >> radec;
    if (!iss.eof() || iss.fail()) {
      throw std::runtime_error("Bad RA/Dec in MPCObservation: " + line);
    }
  }
  {
    istringstream iss(sigWord);
    iss >> sigma;
    if (!iss.eof() || iss.fail()) {
      throw std::runtime_error("Bad sigma in MPCObservation: " + line);
    }
  }
  {
    istringstream iss(obscodeWord);
    iss >> obscode;
    if (!iss.eof() || iss.fail()) {
      throw std::runtime_error("Bad obscode in MPCObservation: " + line);
    }

  }
}

