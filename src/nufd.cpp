#include "nufd.h"

vector<double> fdcoef(unsigned int mord, unsigned int nord, double x0, const vector<double>::const_iterator grid) {
  // this routine implements simple recursions for calculating the weights
  // of finite difference formulas for any order of derivative and any order
  // of accuracy on one-dimensional grids with arbitrary spacing.

  // from Bengt Fornberg's article
  // generation of finite difference formulas on arbitrary spaced grids.
  // math. comp., 51(184):699-706, 1988.

  // input:
  // mord       = the order of the derivative
  // nord       = order of accuracy n
  // x0         = point at which to evaluate the coefficients
  // grid[nord] = array containing the grid starting at the lowest bound
  //              use during finite difference scheme

  // output:
  // coef[nord] = coefficients of the finite difference formula
  vector<double> coef(nord, 0.0);

  // local variables
  int nmmin(min(nord, mord));
  double c1, c2, c3, c4, alpha;

  // more precision for weight calculations results
  // in a smaller error on output coefficients
  long double weight[nmmin + 1][nord][nord];
  for (int i(0); i < nmmin + 1; i++)
    for (int j(0); j < nord; j++)
      for (int k(0); k < nord; k++)
        weight[i][j][k] = 0.0;

  // recursive algorithm implementation
  weight[0][0][0] = 1.0;
  c1 = 1.0;
  for (int nn(1); nn < nord; nn++) {
    c2 = 1.0;
    for (int nu(0); nu < nn; nu++) {
      c3 = grid[nn] - grid[nu];
      c2 = c2 * c3;
      c4 = 1.0 / c3;
      alpha = grid[nn] - x0;
      weight[0][nn][nu] = c4 * (alpha * weight[0][nn - 1][nu]);

      for (int mm(1); mm < nmmin + 1; mm++)
        weight[mm][nn][nu] = c4 * (alpha * weight[mm][nn - 1][nu] - mm * weight[mm - 1][nn - 1][nu]);
    }
    alpha = grid[nn - 1] - x0;
    weight[0][nn][nn] = c1 / c2 * (-alpha * weight[0][nn - 1][nn - 1]);
    c4 = c1 / c2;

    for (int mm(1); mm < nmmin + 1; mm++)
      weight[mm][nn][nn] = c4 * (mm * weight[mm - 1][nn - 1][nn - 1] - alpha * weight[mm][nn - 1][nn - 1]);

    c1 = c2;
  }

  // load the coefficients
  for (int nu(0); nu < nord; nu++)
    coef[nu] = double(weight[mord - 1][nord - 1][nu]);

  return coef;
}

vector<double> fd(unsigned int m, unsigned int n, const vector<double> &grid, const vector<double> &u) {
  // this routine computes the order m derivatives
  // using n points on an arbitrary grid

  // input:
  // m           1=value, 2=1st diff, 3=2nd diff, 4=3rd diff, ...
  // n           = number of points use in fd schemes
  // grid[ngrid] = array of independent values
  // u[ngrid]    = function values at the grid points

  // output:
  // du[ngrid]   = first derivative values at the grid points
  size_t ngrid(grid.size());
  vector<double> du(ngrid, 0.0);
  vector<double> coef(n, 0.0);

  // validate the size of the grid and number of points
  // used in the finite difference scheme
  assert(n <= ngrid);

  // use to point at the first element used
  // in the finite diff scheme to pass to fdcoef
  auto begin = grid.begin();
  auto end = grid.end();

  // number of forward and backward points
  int fb((n - 1) / 2);

  // beginning of the grid (forward differences)
  for (int i(0); i < fb; i++) {
    coef = fdcoef(m, n, grid[i], begin);
    for (int j(0); j < n; j++)
      du[i] = du[i] + coef[j] * u[j];
  }

  // middle of the grid (central differences)
  for (int i(fb); i < ngrid - fb; i++) {
    coef = fdcoef(m, n, grid[i], begin + i - fb);
    for (int j(0); j < n; j++)
      du[i] = du[i] + coef[j] * u[i - fb + j];
  }

  // end of grid (backward differences)
  for (size_t i(ngrid - fb); i < ngrid; i++) {
    coef = fdcoef(m, n, grid[i], end - n);
    for (int j(0); j < n; j++)
      du[i] = du[i] + coef[j] * u[ngrid - n + j];
  }

  return du;
}
