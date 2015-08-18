#include "nufd.h"

// quick example of the non-uniform finite difference algorithm
// to evaluate an analytical functions
// need only fdcoef and fd functions.
int main() {

  // set number of points
  int ngrid(25);
  vector<double> dx(ngrid, 0.0), xgrid(ngrid, 0.0), f(ngrid, 0.0), df(ngrid, 0.0), ddf(ngrid, 0.0), dddf(ngrid, 0.0);
  vector<double> fp(ngrid, 0.0), fpp(ngrid, 0.0), fppp(ngrid, 0.0);
  int nb_points(0), diff_order(0);

  // just to add a bit of random in the 1D grid so the size isn't uniform
  std::mt19937_64 rng;
  uint64_t timeSeed;
  timeSeed = chrono::high_resolution_clock::now().time_since_epoch().count();
  seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
  rng.seed(ss);
  uniform_real_distribution<double> unif(0, 0.25);

  for (int i(1); i<ngrid; i++) {
    dx[i] = 1.5707963/2.0/double(ngrid-1) + unif(rng);
    xgrid[i] = xgrid[i-1] + dx[i];
  }

  for (int i(0); i<ngrid; i++) {
    f[i]   = sin(xgrid[i]);
    df[i]  = cos(xgrid[i]);
    ddf[i] = -f[i]; // - sin
    dddf[i] = -df[i]; // -cos
    // sin
    // cos
    // ...
  }

  // compute the numerical first derivative
  nb_points = 10;
  diff_order = 1;
  fp = fd(diff_order+1, nb_points, xgrid, f);

  // compute the numerical second derivative
  nb_points = 10;
  diff_order = 2;
  fpp = fd(diff_order+1, nb_points, xgrid, f);

  // compute the numerical third derivative
  nb_points = 10;
  diff_order = 3;
  fppp = fd(diff_order+1, nb_points, xgrid, f);

  // compare exact and numerical values
  double diff1(0.0), diff2(0.0), diff3(0.0);
  cout <<  "i" << setw(15)  << "x" << setw(15) << "dx" << setw(15) << "f";
  cout << setw(15) << "df/dx exact" << setw(15) <<  "df/dx numer" << setw(15) <<  "error";
  cout << setw(15) << "2nd exact" << setw(15) << "2nd numer" << setw(15) << "error";
  cout << setw(15) << "3rd exact" << setw(15) << "3rd numer" << setw(15) << "error" << endl;
  for (int i(0); i<ngrid; i++) {
    diff1 = 0.0;
    if (df[i] != 0.0)
      diff1 = (df[i] - fp[i])/df[i];

    diff2 = 0.0;
    if (ddf[i] != 0.0)
      diff2 = (ddf[i] - fpp[i])/ddf[i];

    diff3 = 0.0;
    if (dddf[i] != 0.0)
      diff3 = (dddf[i] - fppp[i])/dddf[i];

    cout << i << setw(15) << xgrid[i] << setw(15) << dx[i] << setw(15) << f[i];
    cout << setw(15) << df[i] << setw(15)  << fp[i] << setw(15)  << diff1;
    cout << setw(15)  << ddf[i] << setw(15)  << fpp[i] << setw(15)  << diff2;
    cout << setw(15) << dddf[i] << setw(15)  << fppp[i] << setw(15)  << diff3 << endl;
  }
  return 0;
}