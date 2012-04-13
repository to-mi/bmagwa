/* BMAGWA software v2.0
 *
 * symmmatrix.hpp
 *
 * http://becs.aalto.fi/en/research/bayes/bmagwa/
 * Copyright 2012 Tomi Peltola <tomi.peltola@aalto.fi>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef SYMMMATRIX_HPP_
#define SYMMMATRIX_HPP_

#include <assert.h>
#include <stdlib.h>
#include <stdexcept>
#include <cmath>
#include <cblas.h>

// LAPACK / LINPACK
extern "C" void dpotrf_(char *uplo, int *n, double *a, int *lda,
                           int *info); // cholesky
extern "C" void dpotri_(char *uplo, int *n, double *a, int *lda,
                           int *info); // inverse matrix from cholesky factor
extern "C" void dchex_(double* r,int* ldr,int* p,int* k,int* l,double* z,
                       int* ldz,int* nz,double* c,double* s,int* job);

namespace bmagwa {

// forward declarations
class VectorView;
class MatrixView;

//! Symmetric or upper triangular matrix which does not own its elements.
/*!
 *  SymmMatrixView provides most of the functionality for symmetric and upper
 *  triangular matrices, but does not own or allocate memory for its elements.
 *  The implementation generally relies on only the upper triangle (incl.
 *  diagonal) of the matrix being accessed, although the implementation
 *  assumes memory layout for full square matrix (and may access/modify the
 *  lower triangle also).
 *
 *  Does not do bounds checking!
 *  Mostly a simple wrapper to BLAS and LAPACK operations.
 *
 *  Subclass SymmMatrix provides memory allocation.
 */
class SymmMatrixView
{
  public:
    // constructors and destructors
    SymmMatrixView(double *data, const size_t length_mem)
    : data_(data), length_(length_mem), length_mem_(length_mem),
      stride_(length_mem)
    {}

    SymmMatrixView(double *data, const size_t length_mem, const size_t length)
    : data_(data), length_(length), length_mem_(length_mem), stride_(length_mem)
    {}

    SymmMatrixView(const SymmMatrixView &m)
    : data_(m.data_), length_(m.length_), length_mem_(m.length_mem_),
      stride_(m.stride())
    {}

    SymmMatrixView(const SymmMatrixView &m, const size_t index,
                   const size_t length)
    : data_(m.data_ + m.stride() * index + index),
      length_(length),
      length_mem_(m.length_mem() - index),
      stride_(m.stride())
    {
      assert(length_mem_ >= 0);
      assert(length_ <= length_mem_);
    }

    virtual ~SymmMatrixView() {}

    // operators
    //! Does a copy of the values (requires equal sized matrices).
    SymmMatrixView& operator=(const SymmMatrixView &m);

    double& operator()(const size_t row, const size_t col)
    {
      return data_[col * stride() + row];
    }

    const double operator()(const size_t row, const size_t col) const
    {
      return data_[col * stride() + row];
    }

    SymmMatrixView& operator=(const double val)
    {
      set_to(val);
      return *this;
    }
    SymmMatrixView& operator*=(const double val);

    // getters and setters
    double* data() { return data_; }
    double* data() const { return data_; }
    const size_t length() const { return length_; }
    const size_t length_mem() const { return length_mem_; }
    const size_t rows() const { return length(); }
    const size_t cols() const { return length(); }
    const size_t stride() const { return stride_; }
    VectorView column(const size_t col);
    VectorView column(const size_t col, const size_t length);
    void resize(const size_t length)
    {
      if (length <= length_mem_) {
        length_ = length;
      }
      else
        throw std::runtime_error(
             "Cannot resize symmetric matrix to larger than allocated memory.");
    }
    void set_to(const double val);
    void set_mem_to(const double val);
    void set_to_identity();
    void set_to_innerproduct(const MatrixView &m);

    // methods

    //! Computes Cholesky decomposition and returns true if matrix is pos.def.
    bool cholesky()
    {
      char uplo = 'U';
      int info = -1, len = length(), offset = stride();

      // xx is upper triangular, cholesky l will also be
      dpotrf_(&uplo, &len, data_, &offset, &info);
      if (info != 0){
        assert(info > 0);
        return false;
        /*if (info < 0){
            throw std::logic_error("Invalid argument to chol");
          } else {
            return false;
          }*/
      }
      return true;
    }

  protected:
    double* data_;
    size_t length_, length_mem_, stride_;
};

//! Symmetric or upper triangular matrix which owns its elements.
/*!
 *  SymmMatrix extends SymmMatrixView to provide memory allocation and
 *  freeing.
 */
class SymmMatrix : public SymmMatrixView
{
  public:
    explicit SymmMatrix(const size_t length_mem)
    : SymmMatrixView(NULL, length_mem)
    {
      init();
    }

    SymmMatrix(const size_t length_mem, const size_t length)
    : SymmMatrixView(NULL, length_mem, length)
    {
      init();
    }

    SymmMatrix(const SymmMatrix &m)
    : SymmMatrixView(NULL, m.length_mem_, m.length_)
    {
      init();
      for (size_t c = 0; c < length_; ++c) {
        size_t offset = c * stride();
        for (size_t r = 0; r <= c; ++r)
          data_[r + offset] = m.data_[r + offset];
      }
    }


    SymmMatrix(const SymmMatrixView &m)
    : SymmMatrixView(NULL, m.length_mem(), m.length())
    {
      init();
      for (size_t c = 0; c < length_; ++c) {
        size_t offset = c * stride();
        for (size_t r = 0; r <= c; ++r)
          data_[r + offset] = m.data()[r + offset];
      }
    }

    ~SymmMatrix()
    {
#ifdef N_ALIGN
      free(data_);
#else
      delete[] data_;
#endif
    }

    // assignment operator is not inherited
    SymmMatrix& operator=(const double val)
    {
      set_to(val);
      return *this;
    }

    //! Does a copy of the values (requires equal sized matrices).
    SymmMatrix& operator=(const SymmMatrix& m);
    //! Does a copy of the values (requires equal sized matrices).
    SymmMatrix& operator=(const SymmMatrixView& m);

    const double* data() const { return data_; }

    //! Removes a column and row from the matrix
    void remove_colrow(const size_t col_rem);

    //! Updates Cholesky on expansion with single row/column.
    /*!
     *  On the expansion of pos.def. matrix XX with single row and column to the
     *  ends, this updates the associated Cholesky decomposition.
     *  \param newcol New column of XX (= transpose of new row of XX).
     *  \param inv_tau2_value If XX is really XX + Tau2^-1, then inv_tau2_value
     *                        should be last element of Tau2^-1 (Tau2 should be
     *                        diagonal). If XX is XX, then give 0.0.
     */
    bool cholesky_update(const VectorView& newcol, const double inv_tau2_value);
    //! Downdates Cholesky on removing a row and column.
    /*!
     *  Note: proper Cholesky decomposition is assumed as a starting point.
     */
    void cholesky_downdate(const size_t col_rem, double* c_tmp, double* s_tmp);

    void swap_colrow(const size_t colrow);

    //! Swaps adjacent rows/cols in Cholesky factorization
    /*!
     *  Swap variables with indices col and col+1.
     *
     *  Optionally update a vector v also (e.g., when v = U'^-1 X y, where U
     *  is the upper triangular Cholesky factor which is updated here).
     */
    void cholesky_swapadj(const size_t col, VectorView* v = NULL);


  private:
    void init()
    {
#ifdef N_ALIGN
      length_mem_ = (size_t)ceil(((double)length_mem_) / N_ALIGN) * N_ALIGN;
      stride_ = length_mem_;
      int res = posix_memalign((void **)&data_, N_ALIGN,
                                length_mem_ * length_mem_ * sizeof(double));
      if (res != 0) throw std::bad_alloc();
#else
      data_ = new double[length_mem_ * length_mem_];
#endif
    }
};

} // namespace bmagwa

#endif /* SYMMMATRIX_HPP_ */
