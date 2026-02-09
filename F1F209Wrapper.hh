// -*- C++ -*-

/* class F1F209Wrapper
 * Wrapper class to get f1, f2, and xs with P. Bosted fitting.
 * Unit of xs is ub/MeV-sr.
 * Valid for all W<5 GeV and all Q2<11 GeV2.
 * Use -fno-leading-underscore -fno-second-underscore when compiling F1F209.f
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   may 2014, Comments modified By Jixie Zhang

#ifndef F1F209WRAPPER_H
#define F1F209WRAPPER_H

class F1F209Wrapper {
public:
  F1F209Wrapper();
  ~F1F209Wrapper();

  void GetF1F2IN09(double Z, double A, double Q2, double W2, double &F1,
                   double &F2, double &rc);
  void Getpind(double W2, double Q2, double &F1, double &R, double &sigt,
               double &sigl);
  void GetF1F2QE09(double Z, double A, double Q2, double W2, double &F1,
                   double &F2);
  void GetMEC2009(double Z, double A, double nu, double Q2, double W2,
                  double &F1, double &F2);

  void Getelas(float eb, float theta, float &crsect);
  void Getelasrad(float es, float theta_d, float t, float wcut,
                  float &crsect_elast_rad, float &corrfact);

  double GetXS(double Z, double A, double Ei, double Ef, double theta);
  double GetXS_QE(double Z, double A, double Ei, double Ef, double theta);
  double GetXS_MEC(double Z, double A, double Ei, double Ef, double theta);
  double GetXS_Ferm_p(double W, double Ei, double theta);

  float GetXSelas(float eb, float theta);
  float GetRadCorrFact(float es, float theta_d, float t, float wcut);
  float GetXSelas_rad(float es, float theta_d, float t, float wcut);
};

#endif
