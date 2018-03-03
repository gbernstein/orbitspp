// Various conversions for OrbitTypes

#include "OrbitTypes.h"
#include "StringStuff.h"

using namespace orbits;

MPCObservation::MPCObservation(const string& line) {
  auto fields = stringstuff::split(line,' ');
  // Clean & strip the list
  for (auto i = fields.begin(); i!=fields.end(); ) {
    stringstuff::stripWhite(*i);
    if (i->size()==0) {
      i = fields.erase(i);
    } else {
      ++i;
    }
  }
  auto wordPtr = fields.begin();
  if (fields.size() < 4) {
    // No way this is enough info
    throw std::runtime_error("Insufficient info on MPCObservation: " + line);
  }
  {
    // First field is indicator of time
    double jd;
    bool fail = false;
    string word = *wordPtr;
    istringstream iss(word);
    iss >> jd;
    if (!iss.eof() || iss.fail()) {
      throw std::runtime_error("Bad date/time in MPCObservation: " + line);
    }

    if (jd < 10000.) {
      // This is probably a year.
      // Try YYYY MM DD.DDDD
      int yr = static_cast<int> (floor(jd));
      int mo;
      float day;
      if (abs(jd-yr) > 1e-7) {
	// Was no integral year
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
      mjd = jdint + inter + day + 1720982 - 0.5;
      
    } else if (jd < 300000) {
      // It's an MJD
      mjd = jd;
    } else {
      // It's a JD
      mjd = jd - MJD0;
    }
    wordPtr++;
  }
  // Now get RA, Dec, obscode
  if (wordPtr==fields.end())
    throw std::runtime_error("Insufficient info on MPCObservation: " + line);
  ++wordPtr;
  string raWord = *wordPtr;

  if (wordPtr==fields.end())
    throw std::runtime_error("Insufficient info on MPCObservation: " + line);
  ++wordPtr;
  string decWord = *wordPtr;

  if (wordPtr==fields.end())
    throw std::runtime_error("Insufficient info on MPCObservation: " + line);
  ++wordPtr;
  string obscodeWord = *wordPtr;

  {
    istringstream iss( raWord + " " + decWord);
    iss >> radec;
    if (!iss.eof() || iss.fail()) {
      throw std::runtime_error("Bad RA/Dec in MPCObservation: " + line);
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

