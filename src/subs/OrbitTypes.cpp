// Various conversions for OrbitTypes

#include "OrbitTypes.h"
#include "StringStuff.h"

#include <string>
#include <iostream>
#include <iomanip>

using namespace orbits;

//********************************
// Here is some nice code to center strings in output fields
// From https://stackoverflow.com/questions/14861018/center-text-in-fixed-width-field-with-stream-manipulators-in-c
template<typename charT, typename traits = std::char_traits<charT> >
class center_helper {
    std::basic_string<charT, traits> str_;
public:
    center_helper(std::basic_string<charT, traits> str) : str_(str) {}
    template<typename a, typename b>
    friend std::basic_ostream<a, b>& operator<<(std::basic_ostream<a, b>& s, const center_helper<a, b>& c);
};

template<typename charT, typename traits = std::char_traits<charT> >
center_helper<charT, traits> centered(std::basic_string<charT, traits> str) {
    return center_helper<charT, traits>(str);
}

// redeclare for std::string directly so we can support anything that implicitly converts to std::string
center_helper<std::string::value_type, std::string::traits_type> centered(const std::string& str) {
    return center_helper<std::string::value_type, std::string::traits_type>(str);
}

template<typename charT, typename traits>
std::basic_ostream<charT, traits>& operator<<(std::basic_ostream<charT, traits>& s, const center_helper<charT, traits>& c) {
    std::streamsize w = s.width();
    if (w > c.str_.length()) {
        std::streamsize left = (w + c.str_.length()) / 2;
        s.width(left);
        s << c.str_;
        s.width(w - left);
        s << "";
    } else {
        s << c.str_;
    }
    return s;
}

// Order to read/write components:
const std::vector<int> abgOrder = {ABG::A, ABG::B, ABG::G,
				   ABG::ADOT, ABG::BDOT, ABG::GDOT};
const std::vector<string> abgNames={"alpha","beta","gamma","alphadot","betadot","gammadot"};
const std::vector<string> elementNames={"a","e","i","LAN","AoP","ToP"}; // Matches order of enum


//********************************

DMatrix
ABG::getXYZ(const DVector& t) const {
  // Calculate inertial XYZ positions at an array of times, return n x 3 matrix
  DMatrix xyz(t.size(), 3, 0.);
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
  
std::ostream& 
ABG::write(std::ostream& os, int precision) const {
  // Write ABG on one line.  All ABG's are +-0.xxx in practice.
  stringstuff::StreamSaver ss(os);  // Cache stream state
  os << std::fixed << std::showpos << std::setprecision(precision);
  for (int i=0; i<6; i++)
    os << std::setw(precision+3) << (*this)[abgOrder[i]] << " ";
  return os;  // Stream state restored on destruction of ss
}

std::ostream& 
ABG::writeHeader(std::ostream& os, int precision) {
  os << "#";
  for (auto& s : abgNames)
    os << std::setw(precision+3) << centered(s) << " ";
  os << endl;
  return os;
}

std::ostream&
orbits::operator<<(std::ostream& os, const ABG& abg) {return abg.write(os);}
std::istream&
orbits::operator>>(std::istream& is, ABG& abg) {return abg.read(is);}
std::ostream&
orbits::operator<<(std::ostream& os, const Elements& el) {return el.write(os);}
std::istream&
orbits::operator>>(std::istream& is, Elements& el) {return el.read(is);}

std::istream& 
ABG::read(std::istream& is) {
  for (int i=0; i<6; i++)
    is >> (*this)[i];
  return is;  // Stream state restored on destruction of ss
}

std::ostream& 
Elements::write(std::ostream& os, int precision) const {
  // Write Elements on one line, in degrees.  Aim to have total
  // number of digits = precision, but without going into scientific notation
  // like std::defaultfloat.
  stringstuff::StreamSaver ss(os);  // Cache stream state
  os << std::fixed << std::noshowpos;
  os << std::setprecision(precision- ((*this)[A]<10. ? 3 : 4))
     << std::setw(precision+1) << (*this)[A] << " ";
  os << std::setprecision(precision) << std::setw(precision+2) 
     << (*this)[E] << " ";
  // Inclination is 0-180 degrees
  os << std::setprecision(precision) << std::setw(precision+4)
     << (*this)[I]/DEGREE << " ";
  // Angles 0-360 degrees
  for (int i=Elements::LAN; i<=Elements::AOP; i++) {
    double degrees = (*this)[i]/DEGREE;
    if (degrees<0) degrees+=360.;
    os << std::setprecision(precision) << std::setw(precision+4)
       << degrees << " ";
  }
  // Time of perihelion - in TDB years post J2000.  Since TNO at opposition
  // moves about 1e-6 degrees in 1e-7 years, add 1 digit.
  os << std::setprecision(precision+1) << std::setw(precision+4)
     << (*this)[Elements::TOP];
  
  return os;  // Stream state restored on destruction of ss
}

std::ostream& 
Elements::writeHeader(std::ostream& os, int precision) {
  os << "#"
     << std::setw(precision+1) << centered(elementNames[Elements::A]) << " "
     << std::setw(precision+2) << centered(elementNames[Elements::E]) << " "
     << std::setw(precision+4) << centered(elementNames[Elements::I]) << " "
     << std::setw(precision+4) << centered(elementNames[Elements::LAN]) << " "
     << std::setw(precision+4) << centered(elementNames[Elements::AOP]) << " "
     << std::setw(precision+4) << centered(elementNames[Elements::TOP])
     << endl;
  return os;
}

std::istream& 
Elements::read(std::istream& is) {
  is >> (*this)[Elements::A]
     >> (*this)[Elements::E]
     >> (*this)[Elements::I]
     >> (*this)[Elements::LAN]
     >> (*this)[Elements::AOP]
     >> (*this)[Elements::TOP];
  (*this)[Elements::I] *= DEGREE;
  (*this)[Elements::LAN] *= DEGREE;
  (*this)[Elements::AOP] *= DEGREE;
  return is;  // Stream state restored on destruction of ss
}


Vector3
Frame::toICRS(const Vector3& x,
	      bool isVelocity) const {
  Vector3 out = orient.m().transpose() * x;
  if (!isVelocity) {
    out += origin.getVector();
  }
  return out;
}

Vector3
Frame::fromICRS(const Vector3& x,
		bool isVelocity) const {
  if (isVelocity)
    return orient.m() * x;
  else
    return orient.m() * (x - origin.getVector());
}

DMatrix
Frame::toICRS(const DMatrix& x,
	      bool isVelocity) const {
  DMatrix out = orient.m().transpose() * x;
  if (!isVelocity)
    out.colwise() += origin.getVector();
  return out;
}

DMatrix
Frame::fromICRS(const DMatrix& x,
		bool isVelocity) const {
  if (isVelocity) {
    return orient.m() * x;
  } else {
    DMatrix tmp = x.colwise() - origin.getVector();
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

std::ostream&
orbits::writeCovariance6(std::ostream& os, const Matrix66& m, int precision) {
  stringstuff::StreamSaver ss(os);
  DVector sd = m.diagonal().cwiseSqrt();
  os << "# Standard deviations: " << endl;
  os << std::scientific << std::setprecision(precision) << std::noshowpos;
  for (int i=0; i<6; i++) 
    os << sd[i] << " ";
  os << endl;
  // Calculate correlation matrix
  os << "# Correlation matrix:" << endl;
  sd = sd.cwiseInverse();
  Matrix66 out = sd.asDiagonal() * m * sd.asDiagonal();
  os << std::fixed << std::showpos;
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++)
      os << std::setw(precision+3) << out(i,j) << " ";
    os << endl;
  }
  return os;
}

Matrix66
orbits::readCovariance6(std::istream& is) {
  string buffer;
  stringstuff::getlineNoComment(is, buffer);
  Vector6 sd;
  {
    std::istringstream iss(buffer);
    for (int i=0; i<6; i++)
      iss >> sd[i];
  }
  Matrix66 cov;
  for (int i=0; i<6; i++) {
    stringstuff::getlineNoComment(is, buffer);
    std::istringstream iss(buffer);
    for (int j=0; j<6; j++)
      iss >> cov(i,j);
  }
  return sd.asDiagonal() * cov * sd.asDiagonal();
}
  
