/* BMAGWA software v1.0
 *
 * vector.cpp
 *
 * http://www.lce.hut.fi/research/mm/bmagwa/
 * Copyright 2011 Tomi Peltola <tomi.peltola@aalto.fi>
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


#include "vector.hpp"
#include "matrix.hpp"
#include "symmmatrix.hpp"

namespace bmagwa {

//
// VectorView
//

// operators
VectorView& VectorView::operator=(const VectorView& v)
{
  if (this != &v){
    assert(length_ == v.length_);
    for (size_t i = 0; i < length_; ++i)
      data_[i] = v.data_[i];
  }
  return *this;
}

Vector VectorView::operator-(const double val)
{
  Vector v(*this);
  for (size_t i = 0; i < v.length(); ++i) v.data_[i] -= val;
  return v;
}

VectorView& VectorView::operator-=(const VectorView& v)
{
  assert(length_ == v.length_);
  for (size_t i = 0; i < length_; ++i) data_[i] -= v(i);
  return *this;
}

VectorView& VectorView::operator-=(const double val)
{
  for (size_t i = 0; i < length_; ++i) data_[i] -= val;
  return *this;
}

VectorView& VectorView::operator/=(const double val)
{
  for (size_t i = 0; i < length_; ++i) data_[i] /= val;
  return *this;
}

// getters and setters
void VectorView::set_to(const double val)
{
  for (size_t i = 0; i < length_; ++i) data_[i] = val;
}

void VectorView::set_mem_to(const double val)
{
  for (size_t i = 0; i < length_mem_; ++i) data_[i] = val;
}

void VectorView::set_to_product(const MatrixView& m, const VectorView& v,
                                const bool transpose)
{
  CBLAS_TRANSPOSE tr;
  if (transpose){
    tr = CblasTrans;
    assert(m.cols() == length());
    assert(m.rows() == v.length());
  } else {
    tr = CblasNoTrans;
    assert(m.cols() == v.length());
    assert(m.rows() == length());
  }

  cblas_dgemv(CblasColMajor, tr, m.rows(), m.cols(), 1.0, m.data_,
              m.stride(), v.data_, 1, 0.0, data_, 1);
}

void VectorView::set_to_vector_diff(const VectorView& v1, const VectorView& v2)
{
  assert(length_ == v1.length_ && length_ == v2.length_);
  for (size_t i = 0; i < length_; ++i) data_[i] = v1(i) - v2(i);
}

// reductions
double VectorView::sum() const
{
  double sum = 0;
  for (size_t i = 0; i < length_; ++i) sum += data_[i];
  return sum;
}

double VectorView::var() const
{
  double sqr_sum = 0;
  double sum = 0;

  for (size_t i = 0; i < length_; ++i){
    sqr_sum += data_[i] * data_[i];
    sum += data_[i];
  }

  return (sqr_sum - sum * sum / length_) / (length_ - 1);
}

double VectorView::max() const
{
  if (length_ < 1) return NAN;
  double max_ = data_[0];
  for (size_t i = 1; i < length_; ++i)
    max_ = std::max(max_, data_[i]);
  return max_;
}

double VectorView::min() const
{
  if (length_ < 1) return NAN;
  double min_ = data_[0];
  for (size_t i = 1; i < length_; ++i)
    min_ = std::min(min_, data_[i]);
  return min_;
}

// methods
VectorView& VectorView::square(){
  for (size_t i = 0; i < length_; ++i) data_[i] = data_[i] * data_[i];
  return *this;
}

void VectorView::multiply_by_triangularmatrix(const SymmMatrixView& m,
                                              const bool transpose)
{
  CBLAS_TRANSPOSE tr;
  if (transpose){
    tr = CblasTrans;
  } else {
    tr = CblasNoTrans;
  }
  assert(m.length() == length());

  cblas_dtrmv(CblasColMajor, CblasUpper, tr, CblasNonUnit,
                m.length(), m.data_, m.stride(), data_, 1);
}

void VectorView::multiply_by_invtriangularmatrix(const SymmMatrixView &m,
                                                 const bool transpose)
{
  CBLAS_TRANSPOSE tr;
  if (transpose){
    tr = CblasTrans;
  } else {
    tr = CblasNoTrans;
  }
  assert(m.length() == length());

  cblas_dtrsv(CblasColMajor, CblasUpper, tr, CblasNonUnit,
                m.length(), m.data_, m.stride(), data_, 1);
}

void VectorView::add_vector(const VectorView& v)
{
  assert(length_ == v.length_);
  for (size_t i = 0; i < length_; ++i) data_[i] += v(i);
}

void VectorView::add_vector(const VectorView& v, const double scale)
{
  assert(length_ == v.length_);
  for (size_t i = 0; i < length_; ++i) data_[i] += scale * v(i);
}


//
// Vector
//

// operators
Vector& Vector::operator=(const Vector &v)
{
  if (this != &v){
    assert(length_ == v.length_);
    for (size_t i = 0; i < length_; ++i)
      data_[i] = v.data_[i];
  }
  return *this;
}

Vector& Vector::operator=(const VectorView &v)
{
  if (this != &v){
    assert(length_ == v.length_);
    for (size_t i = 0; i < length_; ++i)
      data_[i] = v.data_[i];
  }
  return *this;
}


} // namespace bmagwa
