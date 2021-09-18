# matpp
Simple, tested, header only matrix template for C++11

It should be fairly easy to use and integrate into existing projects. Once
created, all matrices are immutable. Internally we are using fixed sized
arrays to hold matrix data. So all size related functionality is checked by
the compiler.


As you might have guessed, we never throw exceptions. So `throw` does not
occur anywhere in code. No explicit use of `new` and `delete` either.


# Usage Examples

Constructing the matrix

```cpp
#include <ostream>
#include <array>
#include <vector>
#include "matpp.hpp"

typedef float real;

int main(){
  return 0;
  // the empty constructor
  MatN<real, 0, 0> me;
  unsigned int vsize = 0;
  auto r = me.get_size(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(0));
  ASSERT_EQUAL(r.status, SUCCESS);

  // constructing with a double vector (vector of vectors)
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

  // constructing it with an array

  std::array<real> r1 = {0, 1, 2, 0, 2, 4};
  m = MatN<real, 2, 3>(r1);
  real vval = 0.0;

  auto r = m.get(0, 2, vval);
  ASSERT_EQUAL(r.status, SUCCESS);
  ASSERT_EQUAL(vval, static_cast<real>(2));

  // constructing it with a single vector where the size of the vector
  // is bigger than the matrix size
  std::vector<real> mv1 = {0, 1, 2, 0, 2, 4, 56, 35};

  m = MatN<real, 2, 3>(mv1);
  // basically it does not really matter. If there are more elements we
  // ignore them

  unsigned int vsize = 0;

  auto r = m.get_size(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(6));
  ASSERT_EQUAL(r.status, SUCCESS);

  // constructing it with a single vector where the size of the vector
  // is smaller than the matrix size
  std::vector<real> mv2 = {0, 1, 2};

  m = MatN<real, 2, 3>(mv2);
  // again it does not matter, we just fill the space with zeros

  vsize = 0;

  auto r = m.get_size(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(6));
  ASSERT_EQUAL(r.status, SUCCESS);

  // constructing it with a scalar value
  m = MatN<real, 2, 3>(static_cast<real>(56.0));

  vsize = 0;

  auto r = m.get_size(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(6));
  ASSERT_EQUAL(r.status, SUCCESS);
  real vval = 0.0;

  auto r = m.get(0, 2, vval);
  ASSERT_EQUAL(r.status, SUCCESS);
  ASSERT_EQUAL(vval, static_cast<real>(56));

  // create a matrix based on another matrix
  MatN<real, 2, 3> m;

  // the function uses the size from the outer matrix
  // and fills the created matrix with zeros casted to the provided type.
  // There is also an overloaded version where you can change the fill value
  auto r = MatN<real>::from_row_cols(m);
  ASSERT_EQUAL(r.status, SUCCESS);

  unsigned int vsize = 0;

  r = m.get_size(vsize);
  ASSERT_EQUAL(vsize, static_cast<unsigned int>(6));
  ASSERT_EQUAL(r.status, SUCCESS);
  //
}

```

Filling a matrix with a certain value:

```cpp
#include <ostream>
#include "matpp.hpp"

typedef float real;

int main(){
  // create the fill value
  real vsize = static_cast<real>(2);

  // create the matrix. Notice they are both the same type
  MatN<real, 2, 3> m(vsize);

  real out = 0;
  auto res = m.get(0, out);
  ASSERT_EQUAL(res.status, SUCCESS);
  ASSERT_EQUAL(out, static_cast<real>(2));
  res = m.get(3, out);
  ASSERT_EQUAL(out, static_cast<real>(2));
  ASSERT_EQUAL(res.status, SUCCESS);

  // if you want to create a value filled with different type just
  // use the from_row_cols method
}

```

Transposing a matrix:

```cpp
#include <ostream>
#include "matpp.hpp"

typedef float real;

int main(){
  // some matrix
  MatN<real, 2, 3> m;
  MatN<real, 10>::from_row_cols(m);

  // output matrix that would hold the transposed data
  MatN<real, 3, 2> out;
  MatN<real, 10>::from_row_cols(out);

  // transposing method. Notice that since all our methods are qualified with
  // const, we are creating a new matrix for holding the transposed data
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

}

```

Doing element wise arithmetic operations


```cpp
#include <ostream>
#include "matpp.hpp"

typedef float real;

int main(){
  std::vector<real> mv = {0, 1, 2, 0, 2, 4};

  // some matrix
  MatN<real, 2, 3> m(mv);

  // some other matrix of the same size
  MatN<real, 2, 3> tout(static_cast<real>(2));

  // out holds the added matrix
  MatN<real, 2, 3> out;

  // MResult holds the state of the operations
  // whether it is successful etc
  MResult r;

  // scalar addition
  r = m.add(2, out);

  // matrix addition
  r = m.add(tout, out);

  // scalar subtraction
  r = m.subtract(2, out);

  // matrix subtraction
  r = m.subtract(tout, out);

  // scalar elementwise product (scalar hadamard product)
  r = m.hadamard_product(2, out);

  // matrix elementwise product (hadamard product)
  r = m.hadamard_product(tout, out);

  // elementwise division
  // checks for zero divisions so be sure to check your result
  // object ensuring the success of the operation

  // scalar division with zero
  r = m.divide(0, out);
  ASSERT_EQUAL(r.status, ARG_ERROR); // true

  // scalar division
  r = m.divide(2, out);
  ASSERT_EQUAL(r.status, SUCCESS); // true

  // elementwise matrix division
  mv2 = {1, 1, 2, 1, 2, 4};
  tout = MatN<real, 2, 3>(mv2);
  r = m.divide(tout, out);

}

```

Gaxpy (General Scale by x Plus by y) y = ax + y


```cpp
#include <ostream>
#include "matpp.hpp"

typedef float real;

int main(){
  // values from Golub, Van Loan 2013, p.5

  std::array<real, 6> Amat_values = {1, 2, 3, 4, 5, 6};
  std::array<real, 2> x = {7, 8};
  std::array<real, 3> y = {0.0, 0.0, 0.0};

  MatN<real, 3, 2> A(Amat_values);
  //
  A.gaxpy(x, y);
  ASSERT_EQUAL(y[0], 23);
  ASSERT_EQUAL(y[1], 53);
  ASSERT_EQUAL(y[2], 83);
}

```

Outer product is implemented as a series of saxpy (scale by x plus by y)
operations

```cpp
#include <ostream>
#include "matpp.hpp"

typedef float real;

int main()
{
  // values from Golub, Van Loan 2013, p.5
  std::array<real, 3> x = {1, 2, 3};
  std::array<real, 2> y = {4, 5};
  MatN<real, 3, 2> out(0);
  MatN<real, 3, 2> m(0);
  //
  auto res = m.outer_product<3, 2>(x, y, out);

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

```

Matrix multiplication is implemented as a serie of outer products

```cpp
#include <ostream>
#include "matpp.hpp"

typedef float real;

int main()
{
  std::array<real> Amat_values = {2, 1, 5, 1, 3, 2};
  std::array<real> Bmat_values = {3, 4, -1, 2, 2, 1};

  MatN<real, 2, 3> A(Amat_values);
  MatN<real, 3, 2> B(Bmat_values);

  // holds the output matrix. Notice that you need to give the size to matrix
  // yourself using the rule A^{NxM} \cdot B^{MxP} = C^{NxP}
  MatN<real, 2, 2> out;
  //

  auto res = A.multiply(B, out);

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

```


