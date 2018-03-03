// Various conversions for OrbitTypes

#include "OrbitTypes.h"
#include "StringStuff.h"

using namespace orbits;

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

