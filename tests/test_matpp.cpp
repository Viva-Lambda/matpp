// test file for quaternion
#include "../mat.hpp"
#include <chrono>
#include <ctest.h>
#include <ratio>

/*! @{
 */

typedef float real;
using namespace matpp;

/*! @{ Test constructors of the matrix
 */
CTEST(suite, test_empty_constructor) {
  MatN<real, 0, 0> m;
  unsigned int vsize = 0;
  auto r = m.get_size(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(0));
  ASSERT_EQUAL(r.status, SUCCESS);
}
CTEST(suite, test_double_vec_constructor) {

  std::vector<std::vector<real>> mv;
  std::vector<real> r1 = {0, 1, 2};
  std::vector<real> r2 = {0, 2, 4};
  mv.push_back(r1);
  mv.push_back(r2);

  MatN<real, 2, 3> m(mv);
  unsigned int vsize = 0;

  auto r = m.get_size(vsize);
  ASSERT_EQUAL(r.status, SUCCESS);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(6));

  r = m.get_col_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(3));
  ASSERT_EQUAL(r.status, SUCCESS);

  r = m.get_row_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(2));
  ASSERT_EQUAL(r.status, SUCCESS);

  std::array<real, 6> arr;
  r = m.get(arr);
  ASSERT_EQUAL(r.status, SUCCESS);
  ASSERT_EQUAL(arr[0], static_cast<real>(0));
  ASSERT_EQUAL(arr[1], static_cast<real>(1));
  ASSERT_EQUAL(arr[2], static_cast<real>(2));
  ASSERT_EQUAL(arr[3], static_cast<real>(0));
  ASSERT_EQUAL(arr[4], static_cast<real>(2));
  ASSERT_EQUAL(arr[5], static_cast<real>(4));
}
CTEST(suite, test_double_vec_constructor_bigger_row) {

  std::vector<std::vector<real>> mv;
  std::vector<real> r1 = {0, 1, 2};
  std::vector<real> r2 = {0, 2, 4};
  std::vector<real> r3 = {1, 2, 4};
  mv.push_back(r1);
  mv.push_back(r2);
  mv.push_back(r3);

  MatN<real, 2, 3> m(mv);
  unsigned int vsize = 0;

  auto r = m.get_size(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(6));
  ASSERT_EQUAL(r.status, SUCCESS);

  r = m.get_col_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(3));
  ASSERT_EQUAL(r.status, SUCCESS);

  r = m.get_row_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(2));
  ASSERT_EQUAL(r.status, SUCCESS);

  std::array<real, 6> arr;
  r = m.get(arr);
  ASSERT_EQUAL(r.status, SUCCESS);
  ASSERT_EQUAL(arr[0], static_cast<real>(0));
  ASSERT_EQUAL(arr[1], static_cast<real>(1));
  ASSERT_EQUAL(arr[2], static_cast<real>(2));
  ASSERT_EQUAL(arr[3], static_cast<real>(0));
  ASSERT_EQUAL(arr[4], static_cast<real>(2));
  ASSERT_EQUAL(arr[5], static_cast<real>(4));
}
CTEST(suite, test_single_vec_constructor) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  unsigned int vsize = 0;

  auto r = m.get_size(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(6));
  ASSERT_EQUAL(r.status, SUCCESS);

  r = m.get_col_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(3));
  ASSERT_EQUAL(r.status, SUCCESS);

  r = m.get_row_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(2));
  ASSERT_EQUAL(r.status, SUCCESS);

  std::array<real, 6> arr;
  r = m.get(arr);
  ASSERT_EQUAL(r.status, SUCCESS);
  ASSERT_EQUAL(arr[0], static_cast<real>(0));
  ASSERT_EQUAL(arr[1], static_cast<real>(1));
  ASSERT_EQUAL(arr[2], static_cast<real>(2));
  ASSERT_EQUAL(arr[3], static_cast<real>(0));
  ASSERT_EQUAL(arr[4], static_cast<real>(2));
  ASSERT_EQUAL(arr[5], static_cast<real>(4));
}
CTEST(suite, test_single_vec_constructor_bigger) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4, 56, 35};

  MatN<real, 2, 3> m(mv);
  unsigned int vsize = 0;

  auto r = m.get_size(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(6));
  ASSERT_EQUAL(r.status, SUCCESS);

  r = m.get_col_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(3));
  ASSERT_EQUAL(r.status, SUCCESS);

  r = m.get_row_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(2));
  ASSERT_EQUAL(r.status, SUCCESS);

  std::array<real, 6> arr;
  r = m.get(arr);
  ASSERT_EQUAL(r.status, SUCCESS);
  ASSERT_EQUAL(arr[0], static_cast<real>(0));
  ASSERT_EQUAL(arr[1], static_cast<real>(1));
  ASSERT_EQUAL(arr[2], static_cast<real>(2));
  ASSERT_EQUAL(arr[3], static_cast<real>(0));
  ASSERT_EQUAL(arr[4], static_cast<real>(2));
  ASSERT_EQUAL(arr[5], static_cast<real>(4));
}
CTEST(suite, test_single_vec_constructor_smaller) {

  std::vector<real> mv = {0, 1, 2};

  MatN<real, 2, 3> m(mv);
  unsigned int vsize = 0;

  auto r = m.get_size(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(6));
  ASSERT_EQUAL(r.status, SUCCESS);

  r = m.get_col_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(3));
  ASSERT_EQUAL(r.status, SUCCESS);

  r = m.get_row_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(2));
  ASSERT_EQUAL(r.status, SUCCESS);

  std::array<real, 6> arr;
  r = m.get(arr);
  ASSERT_EQUAL(r.status, SUCCESS);
  ASSERT_EQUAL(arr[0], static_cast<real>(0));
  ASSERT_EQUAL(arr[1], static_cast<real>(1));
  ASSERT_EQUAL(arr[2], static_cast<real>(2));
  ASSERT_EQUAL(arr[3], static_cast<real>(0));
  ASSERT_EQUAL(arr[4], static_cast<real>(0));
  ASSERT_EQUAL(arr[5], static_cast<real>(0));
}
/*! @}
 */

CTEST(suite, test_static_from_row_cols) {

  MatN<real, 2, 3> m;
  auto r = MatN<real>::from_row_cols(m);
  ASSERT_EQUAL(r.status, SUCCESS);

  unsigned int vsize = 0;

  r = m.get_size(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(6));
  ASSERT_EQUAL(r.status, SUCCESS);

  r = m.get_col_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(3));
  ASSERT_EQUAL(r.status, SUCCESS);

  r = m.get_row_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(2));
  ASSERT_EQUAL(r.status, SUCCESS);

  std::array<real, 6> arr;
  r = m.get(arr);
  ASSERT_EQUAL(r.status, SUCCESS);
  ASSERT_EQUAL(arr[0], static_cast<real>(0));
  ASSERT_EQUAL(arr[1], static_cast<real>(0));
  ASSERT_EQUAL(arr[2], static_cast<real>(0));
  ASSERT_EQUAL(arr[3], static_cast<real>(0));
  ASSERT_EQUAL(arr[4], static_cast<real>(0));
  ASSERT_EQUAL(arr[5], static_cast<real>(0));
}
CTEST(suite, test_static_from_row_cols_fill) {

  MatN<real, 2, 3> m;
  MatN<real, 10, 10>::from_row_cols(20, m);
  unsigned int vsize = 0;

  m.get_size(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(6));

  m.get_col_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(3));

  m.get_row_nb(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(2));

  real out = 0;
  auto r = m.get(0, out);
  ASSERT_EQUAL(out, static_cast<real>(20));
  ASSERT_EQUAL(r.status, SUCCESS);
  std::array<real, 6> arr;
  r = m.get(arr);
  ASSERT_EQUAL(r.status, SUCCESS);
  ASSERT_EQUAL(arr[0], static_cast<real>(20));
  ASSERT_EQUAL(arr[1], static_cast<real>(20));
  ASSERT_EQUAL(arr[2], static_cast<real>(20));
  ASSERT_EQUAL(arr[3], static_cast<real>(20));
  ASSERT_EQUAL(arr[4], static_cast<real>(20));
  ASSERT_EQUAL(arr[5], static_cast<real>(20));
}
CTEST(suite, test_fill) {

  real vsize = static_cast<real>(2);
  MatN<real, 2, 3> m(vsize);

  real out = 0;
  auto res = m.get(0, out);
  ASSERT_EQUAL(res.status, SUCCESS);
  ASSERT_EQUAL(out, static_cast<real>(2));
  res = m.get(3, out);
  ASSERT_EQUAL(out, static_cast<real>(2));
  ASSERT_EQUAL(res.status, SUCCESS);
}
CTEST(suite, test_transpose) {

  MatN<real, 2, 3> m;
  MatN<real, 10>::from_row_cols(m);
  MatN<real, 3, 2> out;
  MatN<real, 10>::from_row_cols(out);

  auto res = m.transpose(out);
  ASSERT_EQUAL(res.status, SUCCESS);
  unsigned int nb = 0;
  //
  res = out.get_col_nb(nb);
  ASSERT_EQUAL(res.status, SUCCESS);
  ASSERT_EQUAL(nb, 2);
  //
  res = out.get_row_nb(nb);
  ASSERT_EQUAL(res.status, SUCCESS);
  ASSERT_EQUAL(nb, 3);
  //
}
CTEST(suite, test_col_nb) {

  MatN<real, 2, 3> m;
  MatN<real, 1>::from_row_cols(m);
  unsigned int cnb = 0;
  auto r = m.get_col_nb(cnb);
  ASSERT_EQUAL(r.status, SUCCESS);
  ASSERT_EQUAL(cnb, static_cast<unsigned int>(3));
  //
}
CTEST(suite, test_row_nb) {
  MatN<real, 2, 3> m;
  MatN<real, 1>::from_row_cols(m);
  unsigned int rnb = 0;
  auto r = m.get_row_nb(rnb);
  ASSERT_EQUAL(r.status, SUCCESS);
  ASSERT_EQUAL(rnb, static_cast<unsigned int>(2));
  //
}
CTEST(suite, test_get_size) {
  MatN<real, 2, 3> m;
  MatN<real, 1>::from_row_cols(m);
  unsigned int snb = 0;
  auto r = m.get_size(snb);
  ASSERT_EQUAL(r.status, SUCCESS);
  ASSERT_EQUAL(snb, static_cast<unsigned int>(6));
  //
}
CTEST(suite, test_add_scalar_value) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  auto r = m.add(2, out);
  ASSERT_EQUAL(r.status, SUCCESS);
  real tval = static_cast<real>(0);
  ASSERT_EQUAL(SUCCESS, (out.get(0, tval)).status);
  ASSERT_EQUAL(tval, 2);
  ASSERT_EQUAL(SUCCESS, (out.get(1, tval)).status);
  ASSERT_EQUAL(tval, 3);
  ASSERT_EQUAL(SUCCESS, out.get(2, tval).status);
  ASSERT_EQUAL(tval, 4);
  ASSERT_EQUAL(SUCCESS, out.get(3, tval).status);
  ASSERT_EQUAL(tval, 2);
  ASSERT_EQUAL(SUCCESS, out.get(4, tval).status);
  ASSERT_EQUAL(tval, 4);
  ASSERT_EQUAL(SUCCESS, out.get(5, tval).status);
  ASSERT_EQUAL(tval, 6);
  //
}
CTEST(suite, test_add_mat_value) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  MatN<real, 2, 3> tout(static_cast<real>(2));
  // filled matrix
  auto r = m.add(tout, out);
  ASSERT_EQUAL(r.status, SUCCESS);
  real tval = static_cast<real>(0);
  ASSERT_EQUAL(SUCCESS, out.get(0, tval).status);
  ASSERT_EQUAL(tval, 2);
  ASSERT_EQUAL(SUCCESS, out.get(1, tval).status);
  ASSERT_EQUAL(tval, 3);
  ASSERT_EQUAL(SUCCESS, out.get(2, tval).status);
  ASSERT_EQUAL(tval, 4);
  ASSERT_EQUAL(SUCCESS, out.get(3, tval).status);
  ASSERT_EQUAL(tval, 2);
  ASSERT_EQUAL(SUCCESS, out.get(4, tval).status);
  ASSERT_EQUAL(tval, 4);
  ASSERT_EQUAL(SUCCESS, out.get(5, tval).status);
  ASSERT_EQUAL(tval, 6);
  //
}
CTEST(suite, test_subtract_scalar_value) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  auto r = m.subtract(2, out);
  ASSERT_EQUAL(r.status, SUCCESS);
  real tval = static_cast<real>(0);
  ASSERT_EQUAL(SUCCESS, out.get(0, tval).status);
  ASSERT_EQUAL(tval, -2);
  ASSERT_EQUAL(SUCCESS, out.get(1, tval).status);
  ASSERT_EQUAL(tval, -1);
  ASSERT_EQUAL(SUCCESS, out.get(2, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(3, tval).status);
  ASSERT_EQUAL(tval, -2);
  ASSERT_EQUAL(SUCCESS, out.get(4, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(5, tval).status);
  ASSERT_EQUAL(tval, 2);
  //
}
CTEST(suite, test_subtract_mat_value) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;

  std::vector<real> mv2 = {0, 1, 2, 0, 2, 4};
  MatN<real, 2, 3> tout(mv2);
  // filled matrix
  auto r = m.subtract(tout, out);
  ASSERT_EQUAL(r.status, SUCCESS);
  real tval = static_cast<real>(0);
  ASSERT_EQUAL(SUCCESS, out.get(0, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(1, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(2, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(3, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(4, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(5, tval).status);
  ASSERT_EQUAL(tval, 0);
  //
}
CTEST(suite, test_hadamard_product_scalar_value) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  auto r = m.hadamard_product(2, out);
  ASSERT_EQUAL(r.status, SUCCESS);
  real tval = static_cast<real>(0);
  ASSERT_EQUAL(SUCCESS, out.get(0, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(1, tval).status);
  ASSERT_EQUAL(tval, 2);
  ASSERT_EQUAL(SUCCESS, out.get(2, tval).status);
  ASSERT_EQUAL(tval, 4);
  ASSERT_EQUAL(SUCCESS, out.get(3, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(4, tval).status);
  ASSERT_EQUAL(tval, 4);
  ASSERT_EQUAL(SUCCESS, out.get(5, tval).status);
  ASSERT_EQUAL(tval, 8);
  //
}
CTEST(suite, test_hadamard_product_mat_value) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;

  std::vector<real> mv2 = {0, 1, 2, 0, 2, 4};
  MatN<real, 2, 3> tout(mv2);
  // filled matrix
  auto r = m.hadamard_product(tout, out);
  ASSERT_EQUAL(r.status, SUCCESS);
  real tval = static_cast<real>(0);
  ASSERT_EQUAL(SUCCESS, out.get(0, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(1, tval).status);
  ASSERT_EQUAL(tval, 1);
  ASSERT_EQUAL(SUCCESS, out.get(2, tval).status);
  ASSERT_EQUAL(tval, 4);
  ASSERT_EQUAL(SUCCESS, out.get(3, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(4, tval).status);
  ASSERT_EQUAL(tval, 4);
  ASSERT_EQUAL(SUCCESS, out.get(5, tval).status);
  ASSERT_EQUAL(tval, 16);
  //
}
CTEST(suite, test_divide_scalar_value_false) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  auto r = m.divide(0, out);
  ASSERT_EQUAL(r.status, ARG_ERROR);
  //
}
CTEST(suite, test_divide_scalar_value_false_check_macro) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  auto r = CHECK(m.divide(0, out));
  ASSERT_EQUAL(r, false);
  //
}
CTEST(suite, test_divide_mat_false) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;

  std::vector<real> mv2 = {0, 1, 2, 0, 2, 4};
  MatN<real, 2, 3> tout(mv2);
  // filled matrix
  auto r = m.divide(tout, out);
  ASSERT_EQUAL(r.status, ARG_ERROR);
  //
}
CTEST(suite, test_divide_scalar_value_true) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;
  auto r = m.divide(2, out);
  ASSERT_EQUAL(r.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  ASSERT_EQUAL(SUCCESS, out.get(0, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(1, tval).status);
  ASSERT_EQUAL(tval, static_cast<real>(0.5));
  ASSERT_EQUAL(SUCCESS, out.get(2, tval).status);
  ASSERT_EQUAL(tval, 1);
  ASSERT_EQUAL(SUCCESS, out.get(3, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(4, tval).status);
  ASSERT_EQUAL(tval, 1);
  ASSERT_EQUAL(SUCCESS, out.get(5, tval).status);
  ASSERT_EQUAL(tval, 2);
}
CTEST(suite, test_divide_mat_true) {

  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  MatN<real, 2, 3> m(mv);
  MatN<real, 2, 3> out;

  std::vector<real> mv2 = {1, 1, 2, 1, 2, 4};
  MatN<real, 2, 3> tout(mv2);
  // filled matrix
  auto r = m.divide(tout, out);
  ASSERT_EQUAL(r.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  ASSERT_EQUAL(SUCCESS, out.get(0, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(1, tval).status);
  ASSERT_EQUAL(tval, 1);
  ASSERT_EQUAL(SUCCESS, out.get(2, tval).status);
  ASSERT_EQUAL(tval, 1);
  ASSERT_EQUAL(SUCCESS, out.get(3, tval).status);
  ASSERT_EQUAL(tval, 0);
  ASSERT_EQUAL(SUCCESS, out.get(4, tval).status);
  ASSERT_EQUAL(tval, 1);
  ASSERT_EQUAL(SUCCESS, out.get(5, tval).status);
  ASSERT_EQUAL(tval, 1);
}
CTEST(suite, test_multiply_mat) {

  // values from S. Lang, Introduction to Linear Algebra,
  // 1986, p. 48

  std::vector<real> Amat_values = {2, 1, 5, 1, 3, 2};
  std::vector<real> Bmat_values = {3, 4, -1, 2, 2, 1};

  MatN<real, 2, 3> A(Amat_values);
  MatN<real, 3, 2> B(Bmat_values);
  MatN<real, 2, 2> out;
  //

  auto start = std::chrono::steady_clock::now();
  auto res = A.multiply(B, out);
  auto stop = std::chrono::steady_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  CTEST_LOG("duration in microseconds %i", (int)duration.count());

  ASSERT_EQUAL(res.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  ASSERT_EQUAL(SUCCESS, out.get(0, tval).status);
  ASSERT_EQUAL(tval, 15);
  ASSERT_EQUAL(SUCCESS, out.get(1, tval).status);
  ASSERT_EQUAL(tval, 15);
  ASSERT_EQUAL(SUCCESS, out.get(2, tval).status);
  ASSERT_EQUAL(tval, 4);
  ASSERT_EQUAL(SUCCESS, out.get(3, tval).status);
  ASSERT_EQUAL(tval, 12);
}
CTEST(suite, test_gaxpy) {

  // values from Golub, Van Loan 2013, p.5

  std::array<real, 6> Amat_values = {1, 2, 3, 4, 5, 6};
  std::array<real, 2> x = {7, 8};
  std::array<real, 3> y = {0.0, 0.0, 0.0};

  MatN<real, 3, 2> A(Amat_values);
  auto start = std::chrono::high_resolution_clock::now();
  //
  auto res = A.gaxpy(x, y);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  CTEST_LOG("duration in microseconds %i", (int)duration.count());

  ASSERT_EQUAL(res.status, SUCCESS);

  //
  ASSERT_EQUAL(y[0], 23);
  ASSERT_EQUAL(y[1], 53);
  ASSERT_EQUAL(y[2], 83);
}
CTEST(suite, test_outer_product) {

  // values from Golub, Van Loan 2013, p.7

  std::array<real, 3> x = {1, 2, 3};
  std::array<real, 2> y = {4, 5};
  MatN<real, 3, 2> out(0);
  MatN<real, 3, 2> m(0);
  //
  auto start = std::chrono::steady_clock::now();
  auto res = m.outer_product<3, 2>(x, y, out);
  auto stop = std::chrono::steady_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  CTEST_LOG("duration in microseconds %i", (int)duration.count());

  ASSERT_EQUAL(res.status, SUCCESS);

  std::array<real, 6> d;
  out.get(d);

  //
  ASSERT_EQUAL(d[0], 4);
  ASSERT_EQUAL(d[1], 5);
  ASSERT_EQUAL(d[2], 8);
  ASSERT_EQUAL(d[3], 10);
  ASSERT_EQUAL(d[4], 12);
  ASSERT_EQUAL(d[5], 15);
}
CTEST(suite, test_multiply_scalar) {

  // values from S. Lang, Introduction to Linear Algebra,
  // 1986, p. 48

  std::vector<real> Amat_values = {2, 1, 5, 1, 3, 2};
  real B = 2;

  MatN<real, 2, 3> A(Amat_values);
  MatN<real, 2, 2> out;
  //
  auto res = A.multiply(B, out);
  //
  ASSERT_EQUAL(res.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  ASSERT_EQUAL(SUCCESS, out.get(0, tval).status);
  ASSERT_EQUAL(tval, 16);
  ASSERT_EQUAL(SUCCESS, out.get(1, tval).status);
  ASSERT_EQUAL(tval, 16);
  ASSERT_EQUAL(SUCCESS, out.get(2, tval).status);
  ASSERT_EQUAL(tval, 12);
  ASSERT_EQUAL(SUCCESS, out.get(3, tval).status);
  ASSERT_EQUAL(tval, 12);
}
CTEST(suite, test_dot_mat) {

  // values from S. Lang, Introduction to Linear Algebra,
  // 1986, p. 48

  std::vector<real> Amat_values = {2, 1, 5, 1, 3, 2};
  std::vector<real> Bmat_values = {3, 4, -1, 2, 2, 1};

  MatN<real, 2, 3> A(Amat_values);
  MatN<real, 3, 2> B(Bmat_values);
  MatN<real, 2, 2> out;
  //
  auto res = A.dot(B, out);
  ASSERT_EQUAL(res.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  ASSERT_EQUAL(SUCCESS, out.get(0, tval).status);
  ASSERT_EQUAL(tval, 15);
  ASSERT_EQUAL(SUCCESS, out.get(1, tval).status);
  ASSERT_EQUAL(tval, 15);
  ASSERT_EQUAL(SUCCESS, out.get(2, tval).status);
  ASSERT_EQUAL(tval, 4);
  ASSERT_EQUAL(SUCCESS, out.get(3, tval).status);
  ASSERT_EQUAL(tval, 12);
}
CTEST(suite, test_dot_scalar) {

  // values from S. Lang, Introduction to Linear Algebra,
  // 1986, p. 48

  std::vector<real> Amat_values = {2, 1, 5, 1, 3, 2};
  real B = 2;

  MatN<real, 2, 3> A(Amat_values);
  MatN<real, 2, 2> out;
  //
  auto res = A.dot(B, out);
  //
  ASSERT_EQUAL(res.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  ASSERT_EQUAL(SUCCESS, out.get(0, tval).status);
  ASSERT_EQUAL(tval, 16);
  ASSERT_EQUAL(SUCCESS, out.get(1, tval).status);
  ASSERT_EQUAL(tval, 16);
  ASSERT_EQUAL(SUCCESS, out.get(2, tval).status);
  ASSERT_EQUAL(tval, 12);
  ASSERT_EQUAL(SUCCESS, out.get(3, tval).status);
  ASSERT_EQUAL(tval, 12);
}
CTEST(suite, test_dot_vector) {
  // values from
  // https://people.math.umass.edu/~havens/m235Lectures/Lecture05.pdf

  std::vector<real> Amat_values = {1, -4, 7, -2, 5, -8, 3, -6, 9};
  // std::array<real, 3> B = {2, 1, -1};
  std::array<real, 3> Bv = {2, 1, -1};

  MatN<real, 3, 3> A(Amat_values);
  MatN<real, 3, 1> B(Bv);
  MatN<real, 3, 1> out;
  //
  auto res = A.dot(B, out);
  //
  ASSERT_EQUAL(res.status, SUCCESS);
  //
  real tval = static_cast<real>(0);
  ASSERT_EQUAL(SUCCESS, out.get(0, tval).status);
  ASSERT_EQUAL(tval, -9);
  ASSERT_EQUAL(SUCCESS, out.get(1, tval).status);
  ASSERT_EQUAL(tval, 9);
  ASSERT_EQUAL(SUCCESS, out.get(2, tval).status);
  ASSERT_EQUAL(tval, -9);
}
CTEST(suite, test_add_rows) {

  std::array<real, 9> new_rows = {1, -4, 7, -2, 5, -8, 3, -6, 9};
  // std::array<real, 3> B = {2, 1, -1};
  std::array<real, 3> Avals = {2, 1, -1};

  MatN<real, 1, 3> A(Avals);
  MatN<real, 4, 3> out;
  //
  auto res = A.add_rows<new_rows.size()>(new_rows, out);
  //
  ASSERT_EQUAL(res.status, SUCCESS);
  std::array<real, 12> check;
  out.get(check);
  //
  ASSERT_EQUAL(check[0], 2);
  ASSERT_EQUAL(check[1], 1);
  ASSERT_EQUAL(check[2], -1);
  ASSERT_EQUAL(check[3], 1);
  ASSERT_EQUAL(check[4], -4);
  ASSERT_EQUAL(check[5], 7);
  ASSERT_EQUAL(check[6], -2);
  ASSERT_EQUAL(check[7], 5);
  ASSERT_EQUAL(check[8], -8);
  ASSERT_EQUAL(check[9], 3);
  ASSERT_EQUAL(check[10], -6);
  ASSERT_EQUAL(check[11], 9);
}
CTEST(suite, test_add_columns) {

  std::array<real, 9> new_columns = {1, -4, 7, -2, 5, -8, 3, -6, 9};
  // std::array<real, 3> B = {2, 1, -1};
  std::array<real, 3> Avals = {2, 1, -1};

  MatN<real, 3, 1> A(Avals);
  MatN<real, 3, 4> out;
  //
  auto res = A.add_columns<new_columns.size()>(new_columns, out);
  //
  ASSERT_EQUAL(res.status, SUCCESS);
  std::array<real, 12> check;
  out.get(check);
  //

  ASSERT_EQUAL(check[0], 2);
  ASSERT_EQUAL(check[1], 1);
  ASSERT_EQUAL(check[2], -2);
  ASSERT_EQUAL(check[3], 3);
  ASSERT_EQUAL(check[4], 1);
  ASSERT_EQUAL(check[5], -4);
  ASSERT_EQUAL(check[6], 5);
  ASSERT_EQUAL(check[7], -6);
  ASSERT_EQUAL(check[8], -1);
  ASSERT_EQUAL(check[9], 7);
  ASSERT_EQUAL(check[10], -8);
  ASSERT_EQUAL(check[11], 9);

}
