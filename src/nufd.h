#ifndef _nufd_
#define _nufd_

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <random>
#include <chrono>
#include <assert.h>

using namespace std;

vector<double> fdcoef(unsigned int mord, unsigned int nord, double x0, const vector<double>::const_iterator grid);
vector<double> fd(unsigned int m, unsigned int n, const vector<double>& grid, const vector<double>& u);

#endif //_nufd_
