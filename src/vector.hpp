/* BMAGWA software v2.0
 *
 * vector.hpp
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


#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include <assert.h>
#include <stdlib.h>
#include <stdexcept>
#include <cmath>
#include <cblas.h>
#include <cstring>
#include <algorithm>

namespace bmagwa {

// forward declarations
class Vector;
class MatrixView;
class SymmMatrixView;

//! A column vector which does not own its elements.
/*!
 *  VectorView implements most of the functionality of vectors, but does not
 *  own or allocate memory for its elements. VectorView is base class of
 *  Vector, which include memory allocation. VectorView can also be mapped to
 *  plain double arrays.
 *
 *  Does not do bounds checking!
 *  Mostly a simple wrapper to BLAS and LAPACK operations.
 *
 *  Subclass Vector provides memory allocation.
 */
class VectorView
{
  public:
    // constructors and destructors
    VectorView()
    : data_(NULL), length_(0), length_mem_(0)
    {}

    VectorView(double *data, const size_t length_mem)
    : data_(data), length_(length_mem), length_mem_(length_mem)
    {}

    VectorView(double *data, const size_t length_mem, const size_t length)
    : data_(data), length_(length), length_mem_(length_mem)
    {}

    VectorView(const VectorView &v, const size_t index, const size_t length)
    : data_(v.data_ + index), length_(length),
      length_mem_(v.length_mem_ - index)
    {}

    VectorView(const VectorView &v)
    : data_(v.data_), length_(v.length_), length_mem_(v.length_mem_)
    {}

    virtual ~VectorView() {}

    // operators
    double& operator()(const size_t index)
    {
      return data_[index];
    }
    const double operator()(const size_t index) const
    {
      return data_[index];
    }
    //! Does a copy of the values (requires equal sized vectors).
    VectorView& operator=(const VectorView &v);
    VectorView& operator=(const double val)
    {
      set_to(val);
      return *this;
    }
    Vector operator-(const double val);
    VectorView& operator-=(const VectorView& v);
    VectorView& operator-=(const double val);
    VectorView& operator/=(const double val);

    // getters and setters
    double* data() { return data_; }
    double* data() const { return data_; }
    const size_t length() const { return length_; }
    const size_t length_mem() const
    {
      return length_mem_;
    }
    VectorView block(const size_t begin, const size_t length)
    {
      assert(begin + length <= length_);
      return VectorView(data_ + begin, length_mem_ - begin, length);
    }
    const VectorView block(const size_t begin, const size_t length) const
    {
      assert(begin + length <= length_);
      return VectorView(data_ + begin, length_mem_ - begin, length);
    }
    void resize(const size_t length)
    {
      if (length <= length_mem_)
        length_ = length;
      else
        throw std::runtime_error(
                       "Cannot resize vector to larger than allocated memory.");
    }
    void set_to(const double val);
    void set_mem_to(const double val);
    void set_to_product(const MatrixView& m, const VectorView& v,
                        const bool transpose);
    void set_to_vector_diff(const VectorView& v1, const VectorView& v2);

    // reductions
    double sum() const;
    inline double mean() const
    {
      return sum() / length_;
    }
    double var() const;
    double max() const;
    double min() const;

    // methods
    VectorView& square();
    void multiply_by_triangularmatrix(const SymmMatrixView& m,
                                      const bool transpose);
    void multiply_by_invtriangularmatrix(const SymmMatrixView &m,
                                         const bool transpose);
    void add_vector(const VectorView& v);
    void add_vector(const VectorView& v, const double scale);

    // static methods
    static double dotproduct(const VectorView &v1,
                                 const VectorView &v2)
    {
      assert(v1.length() == v2.length());
      return cblas_ddot(v1.length(), v1.data_, 1, v2.data_, 1);
    }

  protected:
    double* data_;
    size_t length_, length_mem_;

  /*private:
    friend class Vector;*/
};

//! A column vector which owns its elements.
/*!
 *  Vector extends VectorView to provide memory allocation and freeing.
 *
 *  TODO: reallocation? This would invalidate any mapped VectorViews.
 */
class Vector : public VectorView
{
  public:
    // constructors and destructors
    explicit Vector(size_t length_mem) : VectorView(NULL, length_mem, length_mem)
    {
      init();
    }

    Vector(size_t length_mem, size_t length)
    : VectorView(NULL, length_mem, length)
    {
      init();
    }

    Vector(const Vector &v) : VectorView(NULL, v.length_mem_, v.length_)
    {
      init();
      for (size_t i = 0; i < length_; ++i)
        data_[i] = v.data_[i];
    }

    Vector(const VectorView &v) : VectorView(NULL, v.length_mem(), v.length())
    {
      init();
      for (size_t i = 0; i < length_; ++i)
        data_[i] = v.data()[i];
    }

    ~Vector()
    {
#ifdef N_ALIGN
      free(data_);
#else
      delete[] data_;
#endif
    }

    // operators (assignment operators are not inherited)
    //! Does a copy of the values (requires equal sized vectors).
    Vector& operator=(const Vector &v);
    //! Does a copy of the values (requires equal sized vectors).
    Vector& operator=(const VectorView &v);
    Vector& operator=(const double val)
    {
      set_to(val);
      return *this;
    }

    const double* data() const { return data_; }

    void remove(const size_t ind)
    {
      assert(ind < length());
      memmove(data_ + ind, data_ + ind + 1,
              (length() - ind - 1) * sizeof(double));
      resize(length() - 1);
    }


  private:
    void init()
    {
#ifdef N_ALIGN
      length_mem_ = (size_t)ceil(((double)length_mem_) / N_ALIGN) * N_ALIGN;
      int res = posix_memalign((void **)&data_, N_ALIGN,
          length_mem_ * sizeof(double));
      if (res != 0) throw std::bad_alloc();
#else
      data_ = new double[length_mem_];
#endif
    }
};

} // namespace bmagwa

#endif /* VECTOR_HPP_ */
