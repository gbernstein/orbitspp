# orbitspp
Update of the Bernstein & Khushalani orbit-fitting software for outer solar system bodies.  Improvements to include:
* Use NASA/JPL `cspice` code for planetary ephemerides
* All object-oriented code in C++
* Use `Eigen` for linear algebra
* Improve algorithm performance to avoid wayward solutions for very distant objects and high-e orbits.
