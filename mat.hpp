/*
MIT License

Copyright (c) 2021 Viva Lambda email
<76657254+Viva-Lambda@users.noreply.github.com>

Permission is hereby granted, free of charge, to any person
obtaining a copy
of this software and associated documentation files (the
"Software"), to deal
in the Software without restriction, including without
limitation the rights
to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall
be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef MATPP_HPP
#define MATPP_HPP
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <math.h>
#include <ostream>
#include <stdio.h>
#include <vector>

namespace matpp {

enum status_t : std::uint_least8_t {
  SUCCESS = 1,
  INDEX_ERROR = 2,
  SIZE_ERROR = 3,
  ARG_ERROR = 4,
  LU_ERROR = 5,
  NOT_IMPLEMENTED = 6
};

std::ostream &operator<<(std::ostream &out, status_t flag) {
  switch (flag) {
  case SUCCESS: {
    out << "SUCCESS" << std::endl;
    break;
  }
  case SIZE_ERROR: {
    out << "SIZE_ERROR "
        << " :: " << __FILE__ << " :: " << __LINE__;
    break;
  }
  case INDEX_ERROR: {
    out << "INDEX_ERROR "
        << " :: " << __FILE__ << " :: " << __LINE__;
    break;
  }
  case ARG_ERROR: {
    out << "ARG_ERROR "
        << " :: " << __FILE__ << " :: " << __LINE__;
    break;
  }
  case NOT_IMPLEMENTED: {
    out << "NOT_IMPLEMENTED "
        << " :: " << __FILE__ << " :: " << __LINE__;
    break;
  }
  case LU_ERROR: {
    out << "LU_DECOMPOSITION_ERROR "
        << " :: " << __FILE__ << " :: " << __LINE__;
    break;
  }
  }
  return out;
}

struct MResult {
  //
  unsigned int line_info = 0;
  std::string file_name = "";
  std::string fn_name = "";
  std::string call_name = "";
  std::string duration_info = "";

  status_t status;
  bool success = false;

  MResult() {}
  MResult(unsigned int line, const std::string &fname,
          const std::string &funcname, const std::string &cname, status_t op)
      : line_info(line), file_name(fname), fn_name(funcname), call_name(cname),
        status(op), success(op == SUCCESS) {}
  MResult(unsigned int line, const char *fname, const char *funcname,
          const char *cname, status_t op)
      : line_info(line), file_name(fname), fn_name(funcname), call_name(cname),
        status(op), success(op == SUCCESS) {}
};

template <class T = float, unsigned int RowNb = 1, unsigned int ColNb = RowNb>
class MatN {
  /** holds the vector data*/
  std::array<T, ColNb * RowNb> data;
  static const unsigned int nb_rows = RowNb;
  static const unsigned int nb_cols = ColNb;
  static const unsigned int size = ColNb * RowNb;

public:
  MatN() {
    for (unsigned int i = 0; i < RowNb * ColNb; i++) {
      data[i] = static_cast<T>(0);
    }
    // lu_decomposition = LUdcmp<T, RowNb, ColNb>(data);
  }
  MatN(const std::array<T, RowNb * ColNb> &vd) : data(vd) {}
  MatN(const std::vector<std::vector<T>> &vdata) {

    int rdiff = vdata.size() - RowNb;
    int cdiff = vdata[0].size() - ColNb;
    unsigned int i = 0;
    unsigned int j = 0;

    if (rdiff >= 0 && cdiff >= 0) {
      // row and col size is bigger than current matrix
      for (i = 0; i < RowNb; i++) {
        for (j = 0; j < ColNb; j++) {
          set(i, j, vdata[i][j]);
        }
      }
    } else if (rdiff < 0 && cdiff >= 0) {
      // row size smaller col size is bigger than current
      // matrix
      for (i = 0; i < RowNb; i++) {
        for (j = 0; j < ColNb; j++) {
          T val = static_cast<T>(0);
          if (i < vdata.size()) {
            val = vdata[i][j];
          }
          set(i, j, val);
        }
      }
    } else if (rdiff < 0 && cdiff < 0) {
      // row and col sizes are smaller than current matrix
      for (i = 0; i < RowNb; i++) {
        for (j = 0; j < ColNb; j++) {
          T val = static_cast<T>(0);
          if (i < vdata.size() && j < vdata[0].size()) {
            val = vdata[i][j];
          }
          set(i, j, val);
        }
      }
    } else {
      // row size is bigger and col size is smaller than
      // current matrix
      for (i = 0; i < RowNb; i++) {
        for (j = 0; j < ColNb; j++) {
          T val = static_cast<T>(0);
          if (j < vdata[0].size()) {
            val = vdata[i][j];
          }
          set(i, j, val);
        }
      }
    }
    //
    // lu_decomposition = LUdcmp<T, RowNb, ColNb>(data);
  }
  MatN(const std::vector<T> &vdata) {
    int size_nb = static_cast<unsigned int>(vdata.size()) - RowNb * ColNb;
    unsigned int i = 0;
    if (size_nb < 0) {
      // vector size is smaller than matrix size
      for (i = 0; i < size; i++) {
        if (i < vdata.size()) {
          set(i, vdata[i]);
        } else {
          set(i, static_cast<T>(0));
        }
      }
    } else {
      // vector size is bigger than matrix size
      for (i = 0; i < size; i++) {
        set(i, vdata[i]);
      }
    }

    // lu_decomposition = LUdcmp<T, RowNb, ColNb>(data);
  }
  MatN(T fill_value) {
    for (unsigned int i = 0; i < size; i++) {
      set(i, fill_value);
    }
  }
  template <class K, unsigned int R, unsigned int C>
  friend std::ostream &operator<<(std::ostream &out, MatN<K, R, C> m);

  /**\brief Create matrix based on argument matrix*/
  template <unsigned int OutRowNb = RowNb, unsigned int OutColNb = ColNb>
  static MResult from_row_cols(MatN<T, OutRowNb, OutColNb> &out) {
    out = MatN<T, OutRowNb, OutColNb>(static_cast<T>(0));
    return MResult(__LINE__, __FILE__, __FUNCTION__, "from_row_cols", SUCCESS);
  }
  template <unsigned int OutRowNb = RowNb, unsigned int OutColNb = ColNb>
  static MResult from_row_cols(T v, MatN<T, OutRowNb, OutColNb> &out) {
    MatN<T, OutRowNb, OutColNb> mat(v);
    out = mat;
    return MResult(__LINE__, __FILE__, __FUNCTION__, "from_row_cols", SUCCESS);
  }
  template <unsigned int OutRowNb = RowNb, unsigned int OutColNb = ColNb>
  static MResult identity(unsigned int nb, MatN<T, OutRowNb, OutColNb> &out) {
    MatN<T, OutRowNb, OutColNb> mat;
    auto r = from_row_cols<OutRowNb, OutColNb>(mat);
    if (r.status != SUCCESS)
      return r;
    for (unsigned int i = 0; i < nb; i++) {
      mat.set(i, i, static_cast<T>(1));
    }
    out = mat;
    return MResult(__LINE__, __FILE__, __FUNCTION__, "identity", SUCCESS);
  }
  MResult apply(const MatN<T, RowNb, ColNb> &vmat,
                const std::function<T(T, T)> &fn,
                MatN<T, RowNb, ColNb> &out) const {
    for (unsigned int i = 0; i < data.size(); i++) {
      T tout = static_cast<T>(0);
      vmat.get(i, tout);
      T val = fn(data[i], tout);
      auto r = out.set(i, val);
      if (r.status != SUCCESS)
        return r;
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "apply", SUCCESS);
  }
  MResult apply(const T &v, const std::function<T(T, T)> &fn,
                MatN<T, RowNb, ColNb> &out) const {
    for (unsigned int i = 0; i < data.size(); i++) {
      T val = fn(data[i], v);
      auto r = out.set(i, val);
      if (r.status != SUCCESS)
        return r;
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "apply", SUCCESS);
  }
  // tested
  template <unsigned int OutRowNb = RowNb, unsigned int OutColNb = ColNb>
  MResult fill(T v, MatN<T, OutRowNb, OutColNb> &out) const {
    unsigned int s = 0;
    out.get_size(s);
    for (unsigned int i = 0; i < s; i++) {
      out.set(i, v);
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "fill", SUCCESS);
  }
  // tested
  MResult transpose(MatN<T, ColNb, RowNb> &out) const {

    for (unsigned int i = 0; i < nb_rows; i++) {
      for (unsigned int j = 0; j < nb_cols; j++) {
        T tout = static_cast<T>(0);
        get(i, j, tout);
        out.set(j, i, tout);
      }
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "transpose", SUCCESS);
  }
  // tested
  MResult get_col_nb(unsigned int &v) const {
    v = nb_cols;
    return MResult(__LINE__, __FILE__, __FUNCTION__, "get_col_nb", SUCCESS);
  }
  // tested
  MResult get_row_nb(unsigned int &v) const {
    v = nb_rows;
    return MResult(__LINE__, __FILE__, __FUNCTION__, "get_row_nb", SUCCESS);
  }
  // tested
  MResult get_size(unsigned int &out) const {
    out = size;
    return MResult(__LINE__, __FILE__, __FUNCTION__, "get_size", SUCCESS);
  }
  MResult get(unsigned int row, unsigned int col, T &out) const {
    unsigned int index = row * nb_cols + col;
    MResult r = get(index, out);
    if (r.status != SUCCESS)
      return r;
    return MResult(__LINE__, __FILE__, __FUNCTION__, "get", SUCCESS);
  }
  MResult get(unsigned int index, T &out) const {
    if (index >= data.size())
      return MResult(__LINE__, __FILE__, __FUNCTION__, "get", INDEX_ERROR);
    out = data[index];
    return MResult(__LINE__, __FILE__, __FUNCTION__, "get", SUCCESS);
  }
  MResult get(std::array<T, RowNb * ColNb> &out) const {
    out = data;
    return MResult(__LINE__, __FILE__, __FUNCTION__, "get", SUCCESS);
  }
  MResult set(unsigned int row, unsigned int col, T el) {
    unsigned int index = row * nb_cols + col;
    auto r = set(index, el);
    if (r.status != SUCCESS)
      return r;
    return MResult(__LINE__, __FILE__, __FUNCTION__, "set", SUCCESS);
  }
  MResult set(unsigned int index, T el) {
    if (index >= data.size())
      return MResult(__LINE__, __FILE__, __FUNCTION__, "set", INDEX_ERROR);

    data[index] = el;
    return MResult(__LINE__, __FILE__, __FUNCTION__, "set", SUCCESS);
  }
  MResult get_column(unsigned int index, std::array<T, RowNb> &out) const {
    if (index >= ColNb) {
      return MResult(__LINE__, __FILE__, __FUNCTION__, "get_column",
                     INDEX_ERROR);
    }
    for (unsigned int i = 0; i < RowNb; i++) {
      out[i] = data[i * ColNb + index];
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "get_column", SUCCESS);
  }
  MResult set_column(unsigned int index, const std::array<T, RowNb> &idata) {
    if (index >= ColNb) {
      return MResult(__LINE__, __FILE__, __FUNCTION__, "set_column",
                     INDEX_ERROR);
    }
    for (unsigned int i = 0; i < RowNb; i++) {
      data[i * ColNb + index] = idata[i];
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "set_column", SUCCESS);
  }
  MResult get_row(unsigned int index, std::array<T, ColNb> &out) const {
    if (index >= RowNb) {
      return MResult(__LINE__, __FILE__, __FUNCTION__, "get_row", INDEX_ERROR);
    }
    for (unsigned int i = 0; i < ColNb; i++) {
      out[i] = data[index * ColNb + i];
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "get_row", SUCCESS);
  }
  MResult set_row(unsigned int index, const std::array<T, ColNb> &idata) {
    if (index >= RowNb) {
      return MResult(__LINE__, __FILE__, __FUNCTION__, "set_row", INDEX_ERROR);
    }
    for (unsigned int i = 0; i < ColNb; i++) {
      data[index * ColNb + i] = idata[i];
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "set_row", SUCCESS);
  }

  /**Obtain submatrix TODO*/
  MResult submat(unsigned int row_start, unsigned int col_start,
                 MatN<T, RowNb, ColNb> &out) const {
    unsigned int row_size = nb_rows - row_start;
    unsigned int col_size = nb_cols - col_start;
    return MResult(__LINE__, __FILE__, __FUNCTION__, "submat", NOT_IMPLEMENTED);
  }
  MResult add(const MatN<T, RowNb, ColNb> &v,
              MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv + val; };
    return apply(v, fn, out);
  }
  MResult add(T v, MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv + val; };
    return apply(v, fn, out);
  }
  MResult subtract(const MatN<T, RowNb, ColNb> &v,
                   MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv - val; };
    return apply(v, fn, out);
  }
  MResult subtract(T v, MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv - val; };
    return apply(v, fn, out);
  }
  MResult hadamard_product(const MatN<T, RowNb, ColNb> &v,
                           MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv * val; };
    return apply(v, fn, out);
  }
  MResult hadamard_product(T v, MatN<T, RowNb, ColNb> &out) const {
    auto fn = [](T matv, T val) { return matv * val; };
    return apply(v, fn, out);
  }
  MResult divide(const MatN<T, RowNb, ColNb> &v,
                 MatN<T, RowNb, ColNb> &out) const {
    unsigned int osize = 0;
    v.get_size(osize);
    for (unsigned int i = 0; i < osize; i++) {
      T tout = static_cast<T>(0);
      v.get(i, tout);
      if (tout == static_cast<T>(0)) {
        // zero division risk
        return MResult(__LINE__, __FILE__, __FUNCTION__, "divide", ARG_ERROR);
      }
    }
    auto fn = [](T matv, T val) { return matv / val; };
    return apply(v, fn, out);
  }
  MResult divide(T v, MatN<T, RowNb, ColNb> &out) const {
    if (v == static_cast<T>(0)) {
      return MResult(__LINE__, __FILE__, __FUNCTION__, "divide", ARG_ERROR);
    }
    auto fn = [](T matv, T val) { return matv / val; };
    return apply(v, fn, out);
  }
  /**Declares inner vector product*/
  template <unsigned int N = RowNb>
  MResult vdot(const std::array<T, N> &x, const std::array<T, N> &y, T &out) {
    if (N == 0) {
      return MResult(__LINE__, __FILE__, __FUNCTION__, "vdot", SIZE_ERROR);
    }

    out = static_cast<T>(0);
    for (unsigned int i = 0; i < N; i++) {
      out += x[i] * y[i];
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "vdot", SUCCESS);
  }

  /**Declares inner vector product with scalars*/
  template <unsigned int N = RowNb>
  MResult vdot_s(const std::array<T, N> &x, const T &a,
                 const std::array<T, N> &out) {

    if (N == 0) {
      return MResult(__LINE__, __FILE__, __FUNCTION__, "vdot_s", SIZE_ERROR);
    }
    for (unsigned int i = 0; i < N; i++) {
      out[i] = x[i] * a;
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "vdot_s", SUCCESS);
  }
  /**Implements saxpy algorithm from Golub, Van Loan 2013, p. 4 alg.1.1.2*/
  template <unsigned int N = RowNb>
  MResult saxpy(const T &a, const std::array<T, N> &x,
                std::array<T, N> &y) const {
    if (N == 0) {
      return MResult(__LINE__, __FILE__, __FUNCTION__, "saxpy", SIZE_ERROR);
    }
    for (unsigned int i = 0; i < N; i++) {
      y[i] += x[i] * a; //
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "saxpy", SUCCESS);
  }
  /**
    Implements gaxpy algorithm from Golub, Van Loan 2013, p. 4 alg.1.1.3

    as specified in p. 6-7
   */
  MResult gaxpy(const std::array<T, ColNb> &x, std::array<T, RowNb> &y) const {
    for (unsigned int j = 0; j < ColNb; j++) {
      std::array<T, RowNb> c_j;
      get_column(j, c_j);
      saxpy<RowNb>(x[j], c_j, y);
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "gaxpy", SUCCESS);
  }
  /**
     Implements outer product update from Golub, Van Loan 2013, p. 7
     as a series of saxpy operations
    */
  template <unsigned int Rn, unsigned int Cn>
  MResult outer_product(const std::array<T, Rn> &x, const std::array<T, Cn> &y,
                        MatN<T, Rn, Cn> &out) const {
    for (unsigned int i = 0; i < Rn; i++) {
      std::array<T, Cn> A_i;
      out.get_row(i, A_i);
      saxpy<Cn>(x[i], y, A_i);
      out.set_row(i, A_i);
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "outer_product", SUCCESS);
  }
  template <unsigned int OutColNb = RowNb>
  MResult multiply(T v, MatN<T, RowNb, OutColNb> &out) const {
    // m x n \cdot  vmat (n x l) = out (m x l)
    // RowNb x ColNb \codt (n x l) = out (OutRowNb x
    // OutColNb)
    MatN<T, ColNb, OutColNb> vmat(v);

    auto r = multiply(vmat, out);
    if (r.status != SUCCESS)
      return r;
    return MResult(__LINE__, __FILE__, __FUNCTION__, "multiply", SUCCESS);
  }
  /*matrix to matrix multiplication*/
  template <unsigned int OutColNb = RowNb>
  MResult dot(const MatN<T, ColNb, OutColNb> &v,
              MatN<T, RowNb, OutColNb> &out) const {
    return multiply(v, out);
  }
  /*matrix to scalar multiplication*/
  template <unsigned int OutColNb = RowNb>
  MResult dot(T v, MatN<T, RowNb, OutColNb> &out) const {
    return multiply(v, out);
  }
  /*matrix to vector multiplication*/
  MResult dot(const std::array<T, ColNb> &v, MatN<T, RowNb, 1> &out) const {
    MatN<T, ColNb, 1> vmat(v);
    auto r = multiply<1>(vmat, out);
    if (r.status != SUCCESS)
      return r;
    return MResult(__LINE__, __FILE__, __FUNCTION__, "dot", SUCCESS);
  }

  /**
    m x n \cdot  vmat (n x l) = out (m x l)
    RowNb x ColNb \codt (OutRowNb x OutColNb) = out (RowNb x
    OutColNb)

    We are using the kij (row outer product) variant from Golub, van Loan
    2013, p. 11 alg. 1.1.8 due to implementing this algorithm in C++. For
    fortran etc one should use jki since it access matrices by column.  For a
    comparison of algorithms see table 1.1.1 in p. 9

    tested
   */
  template <unsigned int OutColNb = RowNb>
  MResult multiply(const MatN<T, ColNb, OutColNb> &B,
                   MatN<T, RowNb, OutColNb> &out) const {

    // fill out matrix with zero
    out = MatN<T, RowNb, OutColNb>(static_cast<T>(0));
    for (unsigned int k = 0; k < ColNb; k++) {
      // x vector
      std::array<T, RowNb> A_k;
      get_column(k, A_k);

      // y vector
      std::array<T, OutColNb> B_k;
      B.get_row(k, B_k);

      // compute their outer product
      outer_product<RowNb, OutColNb>(A_k, B_k, out);
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "multiply", SUCCESS);
  }
  /**
    add row
   */
  MResult add_row(const std::array<T, ColNb> &r_data,
                  MatN<T, RowNb + 1, ColNb> &out) const {
    return add_rows<ColNb>(r_data, out);
  }
  /**
    add rows if the incoming data has a size of multiple of
    number of columns
    of this array
  */
  template <unsigned int InRow = ColNb>
  MResult add_rows(const std::array<T, InRow> &r_data,
                   MatN<T, RowNb + (InRow / ColNb), ColNb> &out) const {
    if ((InRow % ColNb) != 0) {
      return MResult(__LINE__, __FILE__, __FUNCTION__, "add_rows", SIZE_ERROR);
    }
    // fill output matrix with zeros
    from_row_cols(out);

    // fill with the output matrix with current matrix
    // elements
    unsigned int i = 0;
    unsigned int j = 0;
    for (i = 0; i < RowNb; i++) {
      for (j = 0; j < ColNb; j++) {
        T value = static_cast<T>(0);
        get(i, j, value);
        out.set(i, j, value);
      }
    }

    // fill from r_data the remaining values
    unsigned int nb_of_rows_to_add = static_cast<unsigned int>(InRow / ColNb);
    for (i = 0; i <= nb_of_rows_to_add; i++) {
      unsigned int row = RowNb + i;
      for (unsigned int j = 0; j < ColNb; j++) {
        T row_val = r_data[i * ColNb + j];
        out.set(row, j, row_val);
      }
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "add_rows", SUCCESS);
  }
  /**
    add column
  */
  MResult add_column(const std::array<T, RowNb> &r_data,
                     MatN<T, RowNb, ColNb + 1> &out) const {
    return add_columns<RowNb>(r_data, out);
  }
  template <unsigned int InCol = RowNb>
  MResult add_columns(const std::array<T, InCol> &c_data,
                      MatN<T, RowNb, ColNb + (InCol / RowNb)> &out) const {
    if ((InCol % RowNb) != 0) {
      return MResult(__LINE__, __FILE__, __FUNCTION__, "add_columns",
                     SIZE_ERROR);
    }
    // fill output matrix with zeros
    from_row_cols(out);

    // fill with the output matrix with current matrix
    // elements
    unsigned int i = 0;
    unsigned int j = 0;
    for (i = 0; i < RowNb; i++) {
      for (j = 0; j < ColNb; j++) {
        T value = static_cast<T>(0);
        get(i, j, value);
        out.set(i, j, value);
      }
    }
    // fill from c_data the remaining values
    unsigned int nb_of_cols_to_add = static_cast<unsigned int>(InCol / RowNb);

    // even if there are zero columns to add the output should be one
    for (i = 0; i < nb_of_cols_to_add; i++) {
      unsigned int col = ColNb + i;
      for (j = 0; j < RowNb; j++) {
        T col_val = c_data[i * RowNb + j];
        out.set(j, col, col_val);
      }
    }
    return MResult(__LINE__, __FILE__, __FUNCTION__, "add_columns", SUCCESS);
  }
  MResult to_double_vec(std::vector<std::vector<T>> &ovec) const {
    //
    std::vector<std::vector<T>> out(nb_rows, std::vector<T>(nb_cols));
    for (unsigned int i = 0; i < nb_rows; i++) {
      for (unsigned int j = 0; j < nb_cols; j++) {
        get(i, j, out[i][j]);
      }
    }
    ovec = out;
    return MResult(__LINE__, __FILE__, __FUNCTION__, "to_double_vec", SUCCESS);
  }
};

template <typename T, unsigned int R, unsigned int C = R>
std::ostream &operator<<(std::ostream &out, MatN<T, R, C> m) {
  std::array<T, R * C> arr;
  m.get(arr);
  for (unsigned int i = 0; i < arr.size(); i++) {
    if (i % C == 0) {
      out << std::endl;
    }
    if (arr[i] >= 0) {
      out << " " << arr[i] << " ";
    } else {
      out << arr[i] << " ";
    }
  }
  return out;
}

bool CHECK(MResult r) { return r.status == SUCCESS; }

/**Checks if operation was successful*/
#define CHECK_M(call, res)                                                     \
  do {                                                                         \
    res = call;                                                                \
    res.call_name = #call;                                                     \
    res.line_info = __LINE__;                                                  \
    res.file_name = __FILE__;                                                  \
  } while (0)

MResult INFO(const MResult res) {
  std::cerr << res.status << " :: " << res.file_name << " :: " << res.line_info
            << " :: " << res.fn_name << " :: " << res.call_name << std::endl;
  return res;
}

/**Prints information about the operation*/
#define INFO_M(call, res)                                                      \
  do {                                                                         \
    res = call;                                                                \
    res.call_name = #call;                                                     \
    res.line_info = __LINE__;                                                  \
    res.file_name = __FILE__;                                                  \
    if (res.status != SUCCESS) {                                               \
      res = INFO(res);                                                         \
    }                                                                          \
  } while (0)

MResult INFO_VERBOSE(MResult res) {
  if (res.status == SUCCESS) {
    std::cerr << "SUCCESS "
              << " :: " << res.file_name << " :: " << res.fn_name
              << " :: " << res.call_name << " :: " << res.duration_info
              << " microseconds" << std::endl;
    return res;
  }
  return INFO(res);
}

/**Prints everything about operation result*/
#define INFO_VERBOSE_M(call, res)                                              \
  do {                                                                         \
    auto start = std::chrono::steady_clock::now();                             \
    res = call;                                                                \
    auto stop = std::chrono::steady_clock::now();                              \
    auto duration =                                                            \
        std::chrono::duration_cast<std::chrono::microseconds>(stop - start);   \
    res.duration_info = to_string(static_cast<int>(duration.count()));         \
    res.call_name = #call;                                                     \
    res.line_info = __LINE__;                                                  \
    res.file_name = __FILE__;                                                  \
    INFO_VERBOSE(res);                                                         \
  } while (0)

} // namespace matpp

#endif
