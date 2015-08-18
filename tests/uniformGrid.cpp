#include "gtest/gtest-spi.h"
#include "gtest/gtest.h"
#include "../nufd.h"

class uniformGrid : public ::testing::Test {
  protected:

  // fixture used in all test cases
  virtual void SetUp() {
    ngrid = 9;
    vector<double> dx(ngrid, 0.0);
    xgrid.resize((unsigned long) ngrid);

    // build an uniform grid form 0-1
    for (int i(1); i < ngrid; i++) {
      dx[i] = 1.0 / double(ngrid - 1);
      xgrid[i] = xgrid[i - 1] + dx[i];
    }

    // fix size
    grid_size = dx[1]-dx[0];
  }
  virtual void TeadDown() {}

  unsigned int ngrid;
  vector<double> xgrid;
  double grid_size;
};

// validate all coefficient for regular grid
// https://en.wikipedia.org/wiki/Finite_difference_coefficient

// Central difference first derivative
TEST_F(uniformGrid, CentralSchemeD1O2) {
  int nb_points(3);
  int diff_order(1);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*grid_size,-1.0/2.0);
  EXPECT_DOUBLE_EQ(coef[1]*grid_size, 0.0);
  EXPECT_DOUBLE_EQ(coef[2]*grid_size, 1.0/2.0);
}

TEST_F(uniformGrid, CentralSchemeD1O4) {
  int nb_points(5);
  int diff_order(1);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*grid_size, 1.0/12.0);
  EXPECT_DOUBLE_EQ(coef[1]*grid_size,-2.0/3.0);
  EXPECT_DOUBLE_EQ(coef[2]*grid_size, 0.0);
  EXPECT_DOUBLE_EQ(coef[3]*grid_size, 2.0/3.0);
  EXPECT_DOUBLE_EQ(coef[4]*grid_size,-1.0/12.0);
}

TEST_F(uniformGrid, CentralSchemeD1O6) {
  int nb_points(7);
  int diff_order(1);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*grid_size,-1.0/60.0);
  EXPECT_DOUBLE_EQ(coef[1]*grid_size, 3.0/20.0);
  EXPECT_DOUBLE_EQ(coef[2]*grid_size,-3.0/4.0);
  EXPECT_DOUBLE_EQ(coef[3]*grid_size, 0.0);
  EXPECT_DOUBLE_EQ(coef[4]*grid_size, 3.0/4.0);
  EXPECT_DOUBLE_EQ(coef[5]*grid_size,-3.0/20.0);
  EXPECT_DOUBLE_EQ(coef[6]*grid_size, 1.0/60.0);
}

TEST_F(uniformGrid, CentralSchemeD1O8) {
  int nb_points(9);
  int diff_order(1);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*grid_size, 1.0/280.0);
  EXPECT_DOUBLE_EQ(coef[1]*grid_size,-4.0/105.0);
  EXPECT_DOUBLE_EQ(coef[2]*grid_size, 1.0/5.0);
  EXPECT_DOUBLE_EQ(coef[3]*grid_size,-4.0/5.0);
  EXPECT_DOUBLE_EQ(coef[4]*grid_size, 0.0);
  EXPECT_DOUBLE_EQ(coef[5]*grid_size, 4.0/5.0);
  EXPECT_DOUBLE_EQ(coef[6]*grid_size,-1.0/5.0);
  EXPECT_DOUBLE_EQ(coef[7]*grid_size, 4.0/105.0);
  EXPECT_DOUBLE_EQ(coef[8]*grid_size,-1.0/280.0);
}

// Central difference second derivative
TEST_F(uniformGrid, CentralSchemeD2O2) {
  int nb_points(3);
  int diff_order(2);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 2.0), 1.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 2.0),-2.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 2.0), 1.0);
}

TEST_F(uniformGrid, CentralSchemeD2O4) {
  int nb_points(5);
  int diff_order(2);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 2.0),-1.0/12.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 2.0), 4.0/3.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 2.0),-5.0/2.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 2.0), 4.0/3.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 2.0),-1.0/12.0);
}

TEST_F(uniformGrid, CentralSchemeD2O6) {
  int nb_points(7);
  int diff_order(2);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 2.0), 1.0/90.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 2.0),-3.0/20.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 2.0), 3.0/2.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 2.0),-49.0/18.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 2.0), 3.0/2.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 2.0),-3.0/20.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 2.0), 1.0/90.0);
}

TEST_F(uniformGrid, CentralSchemeD2O8) {
  int nb_points(9);
  int diff_order(2);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 2.0),-1.0/560.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 2.0), 8.0/315.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 2.0),-1.0/5.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 2.0), 8.0/5.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 2.0),-205.0/72.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 2.0), 8.0/5.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 2.0),-1.0/5.0);
  EXPECT_DOUBLE_EQ(coef[7]*pow(grid_size, 2.0), 8.0/315.0);
  EXPECT_DOUBLE_EQ(coef[8]*pow(grid_size, 2.0),-1.0/560.0);
}

// Central difference third derivative
TEST_F(uniformGrid, CentralSchemeD3O2) {
  int nb_points(5);
  int diff_order(3);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 3.0),-1.0/2.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 3.0), 1.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 3.0), 0.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 3.0),-1.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 3.0), 1.0/2.0);
}

TEST_F(uniformGrid, CentralSchemeD3O4) {
  int nb_points(7);
  int diff_order(3);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 3.0), 1.0/8.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 3.0),-1.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 3.0), 13.0/8.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 3.0), 0.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 3.0),-13.0/8.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 3.0), 1.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 3.0),-1.0/8.0);
}

TEST_F(uniformGrid, CentralSchemeD3O6) {
  int nb_points(9);
  int diff_order(3);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 3.0),-7.0/240.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 3.0), 3.0/10.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 3.0),-169.0/120.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 3.0), 61.0/30.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 3.0), 0.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 3.0),-61.0/30.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 3.0), 169.0/120.0);
  EXPECT_DOUBLE_EQ(coef[7]*pow(grid_size, 3.0),-3.0/10.0);
  EXPECT_DOUBLE_EQ(coef[8]*pow(grid_size, 3.0), 7.0/240.0);
}

// Central difference fourth derivative
TEST_F(uniformGrid, CentralSchemeD4O2) {
  int nb_points(5);
  int diff_order(4);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 4.0), 1.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 4.0),-4.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 4.0), 6.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 4.0),-4.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 4.0), 1.0);
}

TEST_F(uniformGrid, CentralSchemeD4O4) {
  int nb_points(7);
  int diff_order(4);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 4.0),-1.0/6.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 4.0), 2.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 4.0),-13.0/2.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 4.0), 28.0/3.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 4.0),-13.0/2.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 4.0), 2.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 4.0),-1.0/6.0);
}

TEST_F(uniformGrid, CentralSchemeD4O6) {
  int nb_points(9);
  int diff_order(4);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 4.0), 7.0/240.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 4.0),-2.0/5.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 4.0), 169.0/60.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 4.0),-122.0/15.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 4.0), 91.0/8.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 4.0),-122.0/15.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 4.0), 169.0/60.0);
  EXPECT_DOUBLE_EQ(coef[7]*pow(grid_size, 4.0),-2.0/5.0);
  EXPECT_DOUBLE_EQ(coef[8]*pow(grid_size, 4.0), 7.0/240.0);
}

// Central difference fifth derivative
TEST_F(uniformGrid, CentralSchemeD5O2) {
  int nb_points(7);
  int diff_order(5);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 5.0),-1.0/2.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 5.0), 2.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 5.0),-5.0/2.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 5.0), 0.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 5.0), 5.0/2.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 5.0),-2.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 5.0), 1.0/2.0);
}

// Central difference sixth derivative
TEST_F(uniformGrid, CentralSchemeD6O2) {
  int nb_points(7);
  int diff_order(6);
  auto mid = xgrid.begin()+(xgrid.size()-1)/2;
  vector<double> coef = fdcoef(diff_order+1, nb_points, *mid, mid-(nb_points-1)/2);
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 6.0), 1.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 6.0),-6.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 6.0), 15.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 6.0),-20.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 6.0), 15.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 6.0),-6.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 6.0), 1.0);
}

// Forward difference first derivative
TEST_F(uniformGrid, ForwardSchemeD1O1) {
  int nb_points(2);
  int diff_order(1);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*grid_size,-1.0);
  EXPECT_DOUBLE_EQ(coef[1]*grid_size, 1.0);
}

TEST_F(uniformGrid, ForwardSchemeD1O2) {
  int nb_points(3);
  int diff_order(1);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*grid_size,-3.0/2.0);
  EXPECT_DOUBLE_EQ(coef[1]*grid_size, 2.0);
  EXPECT_DOUBLE_EQ(coef[2]*grid_size,-1.0/2.0);
}

TEST_F(uniformGrid, ForwardSchemeD1O3) {
  int nb_points(4);
  int diff_order(1);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*grid_size,-11.0/6.0);
  EXPECT_DOUBLE_EQ(coef[1]*grid_size, 3.0);
  EXPECT_DOUBLE_EQ(coef[2]*grid_size,-3.0/2.0);
  EXPECT_DOUBLE_EQ(coef[3]*grid_size, 1.0/3.0);
}

TEST_F(uniformGrid, ForwardSchemeD1O4) {
  int nb_points(5);
  int diff_order(1);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";

  EXPECT_DOUBLE_EQ(coef[0]*grid_size,-25.0/12.0);
  EXPECT_DOUBLE_EQ(coef[1]*grid_size, 4.0);
  EXPECT_DOUBLE_EQ(coef[2]*grid_size,-3.0);
  EXPECT_DOUBLE_EQ(coef[3]*grid_size, 4.0/3.0);
  EXPECT_DOUBLE_EQ(coef[4]*grid_size,-1.0/4.0);
}

TEST_F(uniformGrid, ForwardSchemeD1O5) {
  int nb_points(6);
  int diff_order(1);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*grid_size,-137.0/60.0);
  EXPECT_DOUBLE_EQ(coef[1]*grid_size, 5.0);
  EXPECT_DOUBLE_EQ(coef[2]*grid_size,-5.0);
  EXPECT_DOUBLE_EQ(coef[3]*grid_size, 10.0/3.0);
  EXPECT_DOUBLE_EQ(coef[4]*grid_size,-5.0/4.0);
  EXPECT_DOUBLE_EQ(coef[5]*grid_size, 1.0/5.0);
}

TEST_F(uniformGrid, ForwardSchemeD1O6) {
  int nb_points(7);
  int diff_order(1);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*grid_size,-49.0/20.0);
  EXPECT_DOUBLE_EQ(coef[1]*grid_size, 6.0);
  EXPECT_DOUBLE_EQ(coef[2]*grid_size,-15.0/2.0);
  EXPECT_DOUBLE_EQ(coef[3]*grid_size, 20.0/3.0);
  EXPECT_DOUBLE_EQ(coef[4]*grid_size,-15.0/4.0);
  EXPECT_DOUBLE_EQ(coef[5]*grid_size, 6.0/5.0);
  EXPECT_DOUBLE_EQ(coef[6]*grid_size,-1.0/6.0);
}

// Forward difference second derivative
TEST_F(uniformGrid, ForwardSchemeD2O1) {
  int nb_points(3);
  int diff_order(2);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 2.0), 1.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 2.0),-2.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 2.0), 1.0);
}

TEST_F(uniformGrid, ForwardSchemeD2O2) {
  int nb_points(4);
  int diff_order(2);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 2.0), 2.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 2.0),-5.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 2.0), 4.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 2.0),-1.0);
}

TEST_F(uniformGrid, ForwardSchemeD2O3) {
  int nb_points(5);
  int diff_order(2);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 2.0), 35.0/12.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 2.0),-26.0/3.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 2.0), 19.0/2.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 2.0),-14.0/3.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 2.0), 11.0/12.0);
}

TEST_F(uniformGrid, ForwardSchemeD2O4) {
  int nb_points(6);
  int diff_order(2);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 2.0), 15.0/4.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 2.0),-77.0/6.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 2.0), 107.0/6.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 2.0),-13.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 2.0), 61.0/12.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 2.0),-5.0/6.0);
}

TEST_F(uniformGrid, ForwardSchemeD2O5) {
  int nb_points(7);
  int diff_order(2);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 2.0), 203.0/45.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 2.0),-87.0/5.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 2.0), 117.0/4.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 2.0),-254.0/9.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 2.0), 33.0/2.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 2.0),-27.0/5.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 2.0), 137.0/180.0);
}

TEST_F(uniformGrid, ForwardSchemeD2O6) {
  int nb_points(8);
  int diff_order(2);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 2.0), 469.0/90.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 2.0),-223.0/10.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 2.0), 879.0/20.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 2.0),-949.0/18.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 2.0), 41.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 2.0),-201.0/10.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 2.0), 1019.0/180.0);
  EXPECT_DOUBLE_EQ(coef[7]*pow(grid_size, 2.0),-7.0/10.0);
}

// Forward difference third derivative
TEST_F(uniformGrid, ForwardSchemeD3O1) {
  int nb_points(4);
  int diff_order(3);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 3.0),-1.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 3.0), 3.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 3.0),-3.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 3.0), 1.0);
}

TEST_F(uniformGrid, ForwardSchemeD3O2) {
  int nb_points(5);
  int diff_order(3);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 3.0),-5.0/2.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 3.0), 9.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 3.0),-12.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 3.0), 7.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 3.0),-3.0/2.0);
}

TEST_F(uniformGrid, ForwardSchemeD3O3) {
  int nb_points(6);
  int diff_order(3);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 3.0),-17.0/4.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 3.0), 71.0/4.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 3.0),-59.0/2.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 3.0), 49.0/2.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 3.0),-41.0/4.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 3.0), 7.0/4.0);
}

TEST_F(uniformGrid, ForwardSchemeD3O4) {
  int nb_points(7);
  int diff_order(3);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 3.0),-49.0/8.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 3.0), 29.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 3.0),-461.0/8.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 3.0), 62.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 3.0),-307.0/8.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 3.0), 13.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 3.0),-15.0/8.0);
}

TEST_F(uniformGrid, ForwardSchemeD3O5) {
  int nb_points(8);
  int diff_order(3);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 3.0),-967.0/120.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 3.0), 638.0/15.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 3.0),-3929.0/40.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 3.0), 389.0/3.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 3.0),-2545.0/24.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 3.0), 268.0/5.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 3.0),-1849.0/120.0);
  EXPECT_DOUBLE_EQ(coef[7]*pow(grid_size, 3.0), 29.0/15.0);
}

TEST_F(uniformGrid, ForwardSchemeD3O6) {
  int nb_points(9);
  int diff_order(3);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 3.0),-801.0/80.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 3.0), 349.0/6.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 3.0),-18353.0/120.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 3.0), 2391.0/10.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 3.0),-1457.0/6.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 3.0), 4891.0/30.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 3.0),-561.0/8.0);
  EXPECT_DOUBLE_EQ(coef[7]*pow(grid_size, 3.0), 527.0/30.0);
  EXPECT_DOUBLE_EQ(coef[8]*pow(grid_size, 3.0),-469.0/240.0);
}

// Forward difference fourth derivative
TEST_F(uniformGrid, ForwardSchemeD4O1) {
  int nb_points(5);
  int diff_order(4);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 4.0), 1.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 4.0),-4.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 4.0), 6.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 4.0),-4.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 4.0), 1.0);
}

TEST_F(uniformGrid, ForwardSchemeD4O2) {
  int nb_points(6);
  int diff_order(4);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 4.0), 3.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 4.0),-14.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 4.0), 26.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 4.0),-24.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 4.0), 11.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 4.0),-2.0);
}

TEST_F(uniformGrid, ForwardSchemeD4O3) {
  int nb_points(7);
  int diff_order(4);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 4.0), 35.0/6.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 4.0),-31.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 4.0), 137.0/2.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 4.0),-242.0/3.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 4.0), 107.0/2.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 4.0),-19.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 4.0), 17.0/6.0);
}

TEST_F(uniformGrid, ForwardSchemeD4O4) {
  int nb_points(8);
  int diff_order(4);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 4.0), 28.0/3.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 4.0),-111.0/2.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 4.0), 142.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 4.0),-1219/6.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 4.0), 176.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 4.0),-185.0/2.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 4.0), 82.0/3.0);
  EXPECT_DOUBLE_EQ(coef[7]*pow(grid_size, 4.0), -7.0/2.0);
}

TEST_F(uniformGrid, ForwardSchemeD4O5) {
  int nb_points(9);
  int diff_order(4);
  vector<double> coef = fdcoef(diff_order+1, nb_points, 0.0, xgrid.begin());
  ASSERT_EQ(coef.size(), nb_points) << "There must be one coefficient by point.";
  EXPECT_DOUBLE_EQ(coef[0]*pow(grid_size, 4.0), 1069.0/80.0);
  EXPECT_DOUBLE_EQ(coef[1]*pow(grid_size, 4.0),-1316.0/15.0);
  EXPECT_DOUBLE_EQ(coef[2]*pow(grid_size, 4.0), 15289.0/60.0);
  EXPECT_DOUBLE_EQ(coef[3]*pow(grid_size, 4.0),-2144.0/5.0);
  EXPECT_DOUBLE_EQ(coef[4]*pow(grid_size, 4.0), 10993.0/24.0);
  EXPECT_DOUBLE_EQ(coef[5]*pow(grid_size, 4.0),-4772.0/15.0);
  EXPECT_DOUBLE_EQ(coef[6]*pow(grid_size, 4.0), 2803.0/20.0);
  EXPECT_DOUBLE_EQ(coef[7]*pow(grid_size, 4.0), -536.0/15.0);
  EXPECT_DOUBLE_EQ(coef[8]*pow(grid_size, 4.0), 967.0/240.0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
};