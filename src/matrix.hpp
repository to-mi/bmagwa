/* BMAGWA software v2.0
 *
 * matrix.hpp
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


#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <assert.h>
#include <stdlib.h>
#include <stdexcept>
#include <cmath>
#include "vector.hpp"

namespace bmagwa {

// forward declarations
class SymmMatrixView;

//! A matrix which does not own its elements.
/*!
 *  MatrixView implements most of the functionality of general matrices, but
 *  does not own or allocate memory for its elements. MatrixView is base class
 *  of Matrix, which include memory allocation. MatrixView can also be mapped to
 *  plain double arrays.
 *
 *  Does not do bounds checking!
 *  Mostly a simple wrapper to BLAS and LAPACK operations.
 *
 *  Subclass Matrix provides memory allocation.
 */
class MatrixView
{
  public:
    // constructors and destructors
    MatrixView(double *data, const size_t rows_mem, const size_t cols_mem)
    : data_(data), rows_(rows_mem), cols_(cols_mem),
      rows_mem_(rows_mem), cols_mem_(cols_mem), stride_(rows_mem)
    {}

    MatrixView(double *data, const size_t rows_mem, const size_t cols_mem,
                 const size_t rows, const size_t cols)
    : data_(data), rows_(rows), cols_(cols),
      rows_mem_(rows_mem), cols_mem_(cols_mem), stride_(rows_mem)
    {}

    MatrixView(const MatrixView &m)
    : data_(m.data_), rows_(m.rows_), cols_(m.cols_),
      rows_mem_(m.rows_mem_), cols_mem_(m.cols_mem_), stride_(m.stride())
    {}

    //! Takes a block above diagonal of a symmetric or upper triangular matrix.
    MatrixView(const SymmMatrixView &m, size_t col);

    virtual ~MatrixView() {}

    // operators
    //! Does a copy of the values (requires equal sized matrices).
    MatrixView& operator=(const MatrixView &m);

    double& operator()(const size_t row, const size_t col)
    {
      return data_[col * stride() + row];
    }
    const double operator()(const size_t row, const size_t col) const
    {
      return data_[col * stride() + row];
    }
    MatrixView& operator=(const double val)
    {
      set_to(val);
      return *this;
    }
    MatrixView& operator*=(const double val);

    // getters and setters
    double* data() { return data_; }
    double* data() const { return data_; }
    const size_t rows() const { return rows_; }
    const size_t cols() const { return cols_; }
    const size_t rows_mem() const { return rows_mem_; }
    const size_t cols_mem() const { return cols_mem_; }
    const size_t stride() const { return stride_; }
    VectorView column(const size_t col)
    {
      assert(col < cols());
      return VectorView(data_ + col * stride(), rows_mem_, rows_);
    }
    const VectorView column(const size_t col) const
    {
      assert(col < cols());
      return VectorView(data_ + col * stride(), rows_mem_, rows_);
    }
    MatrixView columnblock(const size_t begin_col, const size_t ncols)
    {
      assert(begin_col < cols() && begin_col + ncols <= cols());
      return MatrixView(data_ + begin_col * stride(), rows_mem_,
                        cols_mem_ - begin_col, rows_, ncols);
    }
    const MatrixView columnblock(const size_t begin_col, const size_t ncols) const
    {
      assert(begin_col < cols() && begin_col + ncols <= cols());
      return MatrixView(data_ + begin_col * stride(), rows_mem_,
                        cols_mem_ - begin_col, rows_, ncols);
    }
    void resize(const size_t rows, const size_t cols)
    {
      if (rows <= rows_mem_ && cols <= cols_mem_) {
        rows_ = rows;
        cols_ = cols;
      }
      else
        throw std::runtime_error(
                       "Cannot resize matrix to larger than allocated memory.");
    }
    void set_to(const double val);
    void set_mem_to(const double val);

  protected:
    double* data_;
    size_t rows_, cols_, rows_mem_, cols_mem_, stride_;
};

//! A matrix which owns its elements.
/*!
 *  Matrix extends MatrixView to provide memory allocation and freeing.
 *
 *  TODO: reallocation? This would invalidate any mapped MatrixViews.
 */
class Matrix : public MatrixView
{
  public:
    Matrix(const size_t rows_mem, const size_t cols_mem)
    : MatrixView(NULL, rows_mem, cols_mem)
    {
      init();
    }

    Matrix(const size_t rows_mem, const size_t cols_mem,
            const size_t rows, const size_t cols)
    : MatrixView(NULL, rows_mem, cols_mem, rows, cols)
    {
      init();
    }

    Matrix(const Matrix &m)
    : MatrixView(NULL, m.rows_mem(), m.cols_mem(), m.rows(), m.cols())
    {
      init();
      for (size_t c = 0; c < cols_; ++c) {
        size_t offset = c * stride();
        for (size_t r = 0; r < rows_; ++r)
          data_[r + offset] = m.data()[r + offset];
      }
    }

    Matrix(const MatrixView &m)
    : MatrixView(NULL, m.rows_mem(), m.cols_mem(), m.rows(), m.cols())
    {
      init();
      for (size_t c = 0; c < cols_; ++c) {
        size_t offset = c * stride();
        for (size_t r = 0; r < rows_; ++r)
          data_[r + offset] = m.data()[r + offset];
      }
    }

    ~Matrix()
    {
#ifdef N_ALIGN
      free(data_);
#else
      delete[] data_;
#endif
    }

    // assignment operator is not inherited
    Matrix& operator=(const double val)
    {
      set_to(val);
      return *this;
    }

    //! Does a copy of the values (requires equal sized matrices).
    Matrix& operator=(const Matrix& m);
    //! Does a copy of the values (requires equal sized matrices).
    Matrix& operator=(const MatrixView& m);

    const double* data() const { return data_; }

    void remove_column(const size_t col)
    {
    	assert(col < cols());
    	memmove(data_ + col * stride(), data_ + (col + 1) * stride(),
    			(cols() - col - 1) * stride() * sizeof(double));
    	resize(rows(), cols() - 1);
    }

  private:
    void init()
    {
#ifdef N_ALIGN
      rows_mem_ = (size_t)ceil(((double)rows_mem_) / N_ALIGN) * N_ALIGN;
      cols_mem_ = (size_t)ceil(((double)cols_mem_) / N_ALIGN) * N_ALIGN;
      stride_ = rows_mem_;
      int res = posix_memalign((void **)&data_, N_ALIGN,
                                rows_mem_ * cols_mem_ * sizeof(double));
      if (res != 0) throw std::bad_alloc();
#else
      data_ = new double[rows_mem_ * cols_mem_];
#endif
    }
};

} // namespace bmagwa

#endif /* MATRIX_HPP_ */
