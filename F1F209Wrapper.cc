// -*- C++ -*-

/* class F1F209Wrapper
 * Wrapper class to get f1, f2, and xs with P. Bosted fitting.
 * Unit of xs is ub/MeV-sr.
 * Valid for all W<5 GeV and all Q2<11 GeV2.
 * Use -fno-leading-underscore -fno-second-underscore when compiling F1F209.f
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   May 2014, Comments modified By Jixie Zhang

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "F1F209Wrapper.hh"

using namespace std;

extern "C" {
void f1f2in09_(double *Z, double *A, double *Q2, double *W2, double *F1,
               double *F2, double *rc);
void f1f2qe09_(double *Z, double *A, double *Q2, double *W2, double *F1,
               double *F2);
void mec2009_(double *Z, double *A, double *Q2, double *W2, double *F1);
double fitemc_(double *x, double *A, bool *goodfit);
void pind_(double *W2, double *Q2, double *F1, double *R, double *sigt,
           double *sigl);
void elas_(float *eb, float *theta, float *crsect_elast);
void elasrad_(float *es, float *theta_d, float *t, float *wcut,
              float *crsect_elast_rad, float *corrfact);
}

F1F209Wrapper::F1F209Wrapper() {
  // Nothing to do
}

F1F209Wrapper::~F1F209Wrapper() {
  // Nothing to do
}

void F1F209Wrapper::GetF1F2IN09(double Z, double A, double Q2, double W2,
                                double &F1, double &F2, double &rc) {
  double F1o, F2o, rco;
  f1f2in09_(&Z, &A, &Q2, &W2, &F1o, &F2o, &rco);
  F1 = F1o;
  F2 = F2o;
  rc = rco;
}

void F1F209Wrapper::Getpind(double W2, double Q2, double &F1, double &R,
                            double &sigt, double &sigl) {
  double F1o, Ro, sigto, siglo;
  pind_(&W2, &Q2, &F1o, &Ro, &sigto, &siglo);
  F1 = F1o;
  R = Ro;
  sigt = sigto;
  sigl = siglo;
}

void F1F209Wrapper::GetF1F2QE09(double Z, double A, double Q2, double W2,
                                double &F1, double &F2) {
  double F1o, F2o, rco;
  f1f2qe09_(&Z, &A, &Q2, &W2, &F1o, &F2o);
  F1 = F1o;
  F2 = F2o;
}

void F1F209Wrapper::GetMEC2009(double Z, double A, double nu, double Q2,
                               double W2, double &F1, double &F2) {
  const double M = 0.5 * (0.938272 + 0.939565);
  double x = Q2 / 2. / M / nu;
  double w1o, w2o, f1o, R, sigt, sigl;
  double emc_factor;
  bool goodfit = false;

  // printf("%9.6lf\n", Ef);

  pind_(&W2, &Q2, &f1o, &R, &sigt, &sigl);
  const double P6 = 1.10950e-1;
  const double P23 = 4.928e-3;
  R = R * (1. + P6 + P23 * A);

  /*
  if (x > 0.7) {
    F1 = 0;
    F2 = 0;
    return;
  }
  */

  mec2009_(&Z, &A, &Q2, &W2, &w1o);
  w2o = w1o * (1. + R) / (1. + nu * nu / Q2);

  // emc_factor = fitemc_(&x, &A, &goodfit);
  //------------- Manually computing EMC factor-------------- //
  double C[3];
  double ALPHA[9];
  double xp = std::min(0.7, x);

  ALPHA[0] = -6.98871401e-02;
  ALPHA[1] = 2.18888887e+00;
  ALPHA[2] = -2.46673765e+01;
  ALPHA[3] = 1.45290967e+02;
  ALPHA[4] = -4.97236711e+02;
  ALPHA[5] = 1.01312929e+03;
  ALPHA[6] = -1.20839250e+03;
  ALPHA[7] = 7.75766802e+02;
  ALPHA[8] = -2.05872410e+02;

  C[0] = 1.69029097e-02;
  C[1] = 1.80889367e-02;
  C[2] = 5.04268396e-03;

  double c = std::exp(C[0] + C[1] * xp + C[2] * xp * xp);
  double exponent = ALPHA[0];
  for (int i = 1; i <= 8; i++) {
    exponent += ALPHA[i] * std::pow(xp, (double)i);
  }
  emc_factor = c * std::pow(A, exponent);
  //---------------------------------------------------------//

  // std::cout << emc_factor << "\n";
  F1 = M * w1o * emc_factor;
  F2 = nu * w2o * emc_factor;
}

void F1F209Wrapper::Getelas(float eb, float theta, float &crsect) {
  float crsecto;
  elas_(&eb, &theta, &crsecto);
  crsect = crsecto;
}

void F1F209Wrapper::Getelasrad(float es, float theta_d, float t, float wcut,
                               float &crsect_elast_rad, float &corrfact) {
  float crsecto, factor;
  elasrad_(&es, &theta_d, &t, &wcut, &crsecto, &factor);
  crsect_elast_rad = crsecto;
  corrfact = factor;

  //    printf("%9.6lf\n",factor);
}

/*double F1F209Wrapper::GetXS(double Z, double A, double Ei, double Ef, double
theta)
{
    const double Mp = 0.938272;
    const double Mn = 0.939565;
    const double M = 0.5*(Mp+Mn);
  //  double Ef;
    double pi = 3.14159265358979323846;
    theta=theta*3.14159265358979323846/180.;
 //   Ef= (W*W - 2.*M*Ei - M*M)/(2.*Ei*(cos(theta)-1)-2.*M);

    double nu = Ei-Ef;
    double Q2 = 4.*Ei*Ef*(sin(abs(theta)/2.))*(sin(abs(theta)/2.));

   double w2 = M*M + 2.*M*nu-Q2;

    double F1, F2, r;
    double xs1, xs2;

    double x_bjork,F3,F4;

 //   printf("%9.6lf %9.6lf  \n", Ef, sqrt(w2));


    double Ge_p, q, denom1, Gm_p, denom2;
    double Ge_n, gd, temp, Gm_n, denom3;
    double tau_p, tau_n,tau, R, eps;
    double frac, Mott,xsec, Jacobian, x_bojrk;
    double TH, L, b1,b2,b3,R_par;

    tau_p = Q2/4./Mp/Mp;
    tau_n = Q2/4./Mn/Mn;
    tau = Q2/4./M/M;
    q = sqrt(Q2);



   denom1 = 1.0 + 0.62*q + 0.68*pow(q,2) + 2.80*pow(q,3) + 0.83*pow(q,4);
   Ge_p = 1/denom1;


   denom2 = 1.0 + 0.35*q + 2.44*pow(q,2) + 0.50*pow(q,3) + 1.04*pow(q,4) +
0.34*pow(q,5); Gm_p = 2.793/denom2;


  temp = 1+Q2/0.71;
  gd = pow(temp,-2);
  Ge_n = (-1.0*1.25*-1.913*tau*gd)/(1+18.3*tau);


  denom3 = 1.0 - 1.74*q + 9.29*pow(q,2) - 7.63*pow(q,3) + 4.63*pow(q,4);
  Gm_n = -1.913/denom3;


   // R = (Ge_p+Ge_n)*(Ge_p+Ge_n)/(tau*(Gm_p+Gm_n)*(Gm_p+Gm_n));

   R =  (Ge_p*Ge_p + Ge_n*Ge_n)/(tau*(Gm_p*Gm_p+Gm_n*Gm_n));

   x_bojrk = Q2/2./M/nu;
   TH = 1. + 12.*Q2/(Q2+1)*0.125*0.125/(0.125*0.125+x_bojrk*x_bojrk);

  L = 0.2;
  b1 = 0.635;
  b2 = 0.5747;
  b3 = -0.3534;

   R_par = b1/log(Q2/L/L)*TH + b2/Q2 + b3/(Q2*Q2 + 0.3*0.3);



    eps = 1./(1. +2.*(nu*nu+Q2)/Q2*tan(abs(theta) / 2.)*tan(abs(theta) / 2.));

    frac = 1. + ((1.-eps)/eps)*(1./(1+R));

    Mott =  389.379*(2. / 137. * Ef / Q2 * cos(theta / 2.))*(2. / 137. * Ef / Q2
* cos(theta / 2.));


    Jacobian = x_bojrk*Ei*Ef/pi/nu;
 //   xsec = 1./nu*F2*Mott*frac;



 //   GetF1F2IN09(Z, A, Q2, w2, F1, F2, r);

 //   xs1 = (2. / 137. * Ef / Q2 * cos(theta / 2.))*(2. / 137. * Ef / Q2 *
cos(theta / 2.)); // mott
  //  xs1 = xs1 * (2. / M * F1 * tan(abs(theta) / 2.) * tan(abs(theta) / 2.) +
F2 / nu);
  //  xs1 = xs1 * 389.379;
//
//   GetF1F2QE09(Z, A, 0.475, 0.88, F1, F2);

//   printf("%9.6lf %9.6lf %9.6lf %9.6lf\n ",Q2,sqrt(w2),F1, F2/2);

    GetF1F2QE09(Z, A, Q2, w2, F1, F2);

     xsec = 1./nu*F2*Mott*frac;


  //  xs2 = (2. / 137. * Ef / Q2 * cos(theta / 2.))*(2. / 137. * Ef / Q2 *
cos(theta / 2.)); // mott
 //   xs2 = xs2 * (2. / M * F1 * tan(abs(theta) / 2.) * tan(abs(theta) / 2.) +
F2 / nu);
  //  xs2 = xs2 * 389.379;//hc squared

   //  xs2 = xs2;


     //*(M-2*Ei*sin(theta/2.)*sin(theta/2.))/W;
     //*(M+2*Ei*sin(theta/2.)*sin(theta/2.))/W
     //*(4.*M*Ei*Ef*sin(theta))
    return xsec; // ub/GeV-sr
}*/

/*double F1F209Wrapper::GetXS(double Z, double A, double Ei, double Ef, double
theta)
{
    const double M = 0.93825;
 theta=theta*3.14159265358979323846/180.;
    double nu = Ei - Ef;
    double Q2 = 4. * Ei * Ef * (sin(abs(theta) / 2.))*(sin(abs(theta) / 2.));
    double w2 = M * M + 2. * M * nu - Q2;

    double F1, F2, r,tau,eps, Mott;
    double xs1, xs2;

 //   GetF1F2IN09(Z, A, Q2, w2, F1, F2, r);

 //   xs1 = (2. / 137. * Ef / Q2 * cos(theta / 2.))*(2. / 137. * Ef / Q2 *
cos(theta / 2.)); // mott
//    xs1 = xs1 * (2. / M * F1 * tan(abs(theta) / 2.) * tan(abs(theta) / 2.) +
F2 / nu);
//    xs1 = xs1 * 389.379;

//    GetF1F2QE09(Z, A, Q2, w2, F1, F2);
//    Mott =  389.379*(2. / 137. * Ef / Q2 * cos(theta / 2.))*(2. / 137. * Ef /
Q2 * cos(theta / 2.)); // mott



//        tau = Q2/4./M/M;
//    eps = 1./(1. +2.*(nu*nu+Q2)/Q2*tan(abs(theta) / 2.)*tan(abs(theta) / 2.));

//   xs2 =   Mott*(F1 + tau/eps*F2)/(1. + tau);

//xs2=  Mott*((F1 +tau*(F2))/(1+tau) +2*tan(theta/2.)*tan(theta/2.)*tau*(F2));
// Xsec_p = Mott*(Ge_p*Ge_p + tau/eps*Gm_p*Gm_p)/(1. + tau);
//Xsec_n = Mott*(Ge_n*Ge_n + tau/eps*Gm_n*Gm_n)/(1. + tau);
  //  xs2 = xs2 * (2. / M * F1 * tan(abs(theta) / 2.) * tan(abs(theta) / 2.) +
F2 / nu);
 //   xs2 = xs2 * 389.379;

  //   xs2

//    xs2=  Mott*(F1+F2);
    GetF1F2QE09(Z, A, Q2, w2, F1, F2);
    xs2 = (2. / 137. * Ef / Q2 * cos(theta / 2.))*(2. / 137. * Ef / Q2 *
cos(theta / 2.)); // mott
   // xs2 = xs2 * (2. / M * F1 * tan(abs(theta) / 2.) * tan(abs(theta) / 2.) +
F2 / nu); xs2 = xs2 * (F1+F2); xs2 = xs2 * 389.379;

    return xs2; // ub/GeV-sr
}*/

double F1F209Wrapper::GetXS_QE(double Z, double A, double Ei, double Ef,
                               double theta) {
  const double M = 0.5 * (0.938272 + 0.939565);
  theta = theta * 3.14159265358979323846 / 180.;
  double nu = Ei - Ef;
  double Q2 = 4. * Ei * Ef * (sin(abs(theta) / 2.)) * (sin(abs(theta) / 2.));
  double w2 = M * M + 2. * M * nu - Q2;

  double F1, F2, r;
  double xs1, xs2;

  GetF1F2QE09(Z, A, Q2, w2, F1, F2);
  xs2 = (2. / 137. * Ef / Q2 * cos(theta / 2.)) *
        (2. / 137. * Ef / Q2 * cos(theta / 2.)); // mott
  xs2 = xs2 *
        (2. / M * F1 * tan(abs(theta) / 2.) * tan(abs(theta) / 2.) + F2 / nu);
  xs2 = xs2 * 389.379;

  return xs2; //    ub/GeV sr
}

double F1F209Wrapper::GetXS_MEC(double Z, double A, double Ei, double Ef,
                                double theta) {
  const double M = 0.5 * (0.938272 + 0.939565);
  theta = theta * 3.14159265358979323846 / 180.;
  double nu = Ei - Ef;
  double Q2 = 4. * Ei * Ef * (sin(abs(theta) / 2.)) * (sin(abs(theta) / 2.));
  double w2 = M * M + 2. * M * nu - Q2;

  double F1, F2;
  double xs;

  if ((w2 <= 0) || (A < 2.5))
    return 0.0;

  GetMEC2009(Z, A, nu, Q2, w2, F1, F2);

  xs = (2. / 137. * Ef / Q2 * cos(theta / 2.)) *
       (2. / 137. * Ef / Q2 * cos(theta / 2.)); // mott;

  xs = xs *
       (2. / M * F1 * tan(abs(theta) / 2.) * tan(abs(theta) / 2.) + F2 / nu);
  xs = xs * 389.379;
  return xs;
}

double F1F209Wrapper::GetXS(double Z, double A, double Ei, double Ef,
                            double theta) {
  const double M = 0.5 * (0.938272 + 0.939565);
  theta = theta * 3.14159265358979323846 / 180.;
  double nu = Ei - Ef;
  double Q2 = 4. * Ei * Ef * (sin(abs(theta) / 2.)) * (sin(abs(theta) / 2.));
  double w2 = M * M + 2. * M * nu - Q2;

  double F1, F2, r;
  double xs1, xs2;

  GetF1F2IN09(Z, A, Q2, w2, F1, F2, r);

  xs1 = (2. / 137. * Ef / Q2 * cos(theta / 2.)) *
        (2. / 137. * Ef / Q2 * cos(theta / 2.)); // mott
  xs1 = xs1 *
        (2. / M * F1 * tan(abs(theta) / 2.) * tan(abs(theta) / 2.) + F2 / nu);
  xs1 = xs1 * 389.379;

  GetF1F2QE09(Z, A, Q2, w2, F1, F2);
  xs2 = (2. / 137. * Ef / Q2 * cos(theta / 2.)) *
        (2. / 137. * Ef / Q2 * cos(theta / 2.)); // mott
  xs2 = xs2 *
        (2. / M * F1 * tan(abs(theta) / 2.) * tan(abs(theta) / 2.) + F2 / nu);
  xs2 = xs2 * 389.379;

  return xs1 + xs2; // ub/GeV-sr
}

double F1F209Wrapper::GetXS_Ferm_p(double W, double Ei, double theta) {
  const double M = 0.93825;

  double Ef;
  double pi = 3.14159265358979323846;
  double alpha = 1. / 137.;
  theta = theta * pi / 180.;
  Ef = (W * W - 2. * M * Ei - M * M) / (2. * Ei * (cos(theta) - 1) - 2. * M);

  double nu = Ei - Ef;
  double Q2 = 4. * Ei * Ef * (sin(abs(theta) / 2.)) * (sin(abs(theta) / 2.));
  //  double w2 = M * M + 2. * M * nu - Q2;
  double w2 = W * W;
  double F1, R, sigt, sigl, epsilon, SIGMA, Gamma;

  // printf("%9.6lf\n", Ef);

  Getpind(w2, Q2, F1, R, sigt, sigl);
  epsilon = 1. / (1. + 2. * (1. + nu * nu / (Q2)) * tan(theta / 2.) *
                           tan(theta / 2.));

  Gamma = (alpha * (w2 - M * M) * Ef) /
          (2 * pi * pi * Q2 * 2 * M * Ei * (1. - epsilon));
  // Gamma = (alpha*W*(w2-M*M))/(4*pi*Ei*Ei*M*M*Q2*(1-epsilon));
  //  Gamma = (alpha*W*(w2-M*M))/(4*pi*Ei*Ei*M*M*Q2*(1-epsilon));
  SIGMA = (sigt + epsilon * sigl) * Gamma * W / M;
  //    printf("%9.6lf\n", epsilon);
  //   printf("%9.6lf\n", Q2);
  // GetF1F2IN09(Z, A, Q2, w2, F1, F2, r);
  /*
      xs1 = (2. / 137. * Ef / Q2 * cos(theta / 2.))*(2. / 137. * Ef / Q2 *
     cos(theta / 2.)); // mott xs1 = xs1 * (2. / M * F1 * tan(abs(theta) / 2.) *
     tan(abs(theta) / 2.) + F2 / nu); xs1 = xs1 * 389.379;

      GetF1F2QE09(Z, A, Q2, w2, F1, F2);
      xs2 = (2. / 137. * Ef / Q2 * cos(theta / 2.))*(2. / 137. * Ef / Q2 *
     cos(theta / 2.)); // mott xs2 = xs2 * (2. / M * F1 * tan(abs(theta) / 2.) *
     tan(abs(theta) / 2.) + F2 / nu); xs2 = xs2 * 389.379;

   */
  return SIGMA; // ub/GeV-sr
}

float F1F209Wrapper::GetXSelas(float eb, float theta) {
  float crsect;
  Getelas(eb, theta, crsect);
  return crsect;
};

float F1F209Wrapper::GetRadCorrFact(float es, float theta_d, float t,
                                    float wcut) {
  float factor, crsect;
  // double esr,theta_dr,tr,wcutr,factr;
  Getelasrad(es, theta_d, t, wcut, crsect, factor);
  // printf("%9.6lf\n",factor);
  // factor = factr;
  return factor;
};

float F1F209Wrapper::GetXSelas_rad(float es, float theta_d, float t,
                                   float wcut) {
  float factor, crsect;
  // double esr,theta_dr,tr,wcutr,factr;
  Getelasrad(es, theta_d, t, wcut, crsect, factor);
  // printf("%9.6lf\n",factr);
  // factor = factr;
  return crsect;
};
