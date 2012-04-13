/* BMAGWA software v2.0
 *
 * matrix.cpp
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


#include "matrix.hpp"
#include "symmmatrix.hpp"

namespace bmagwa {

MatrixView::MatrixView(const SymmMatrixView &m, size_t col)
: data_(m.data() + m.stride() * col),
  rows_(col),
  cols_(m.length() - col),
  rows_mem_(m.length_mem()),
  cols_mem_(m.length_mem() - col),
  stride_(m.stride())
{
  assert(col > 0 && cols_mem_ > 0);
}

MatrixView& MatrixView::operator=(const MatrixView &m)
{
  if (this != &m){
    assert(rows_ == m.rows_ && cols_ == m.cols_);
    for (size_t c = 0; c < cols_; ++c) {
      size_t offset = c * stride();
      for (size_t r = 0; r < rows_; ++r)
        data_[r + offset] = m.data_[r + offset];
    }
  }
  return *this;
}

MatrixView& MatrixView::operator*=(const double val)
{
  for (size_t c = 0; c < cols_; ++c) {
    size_t offset = c * stride();
    for (size_t r = 0; r < rows_; ++r)
      data_[r + offset] *= val;
  }
  return *this;
}

void MatrixView::set_to(const double val)
{
  for (size_t c = 0; c < cols_; ++c) {
    size_t offset = c * stride();
    for (size_t r = 0; r < rows_; ++r)
      data_[r + offset] = val;
  }
}

void MatrixView::set_mem_to(const double val)
{
  for (size_t c = 0; c < cols_mem_; ++c) {
    size_t offset = c * stride();
    for (size_t r = 0; r < rows_mem_; ++r)
      data_[r + offset] = val;
  }
}

Matrix& Matrix::operator=(const Matrix& m)
{
  if (this != &m){
    assert(rows_ == m.rows_ && cols_ == m.cols_);
    for (size_t c = 0; c < cols_; ++c) {
      size_t offset = c * stride();
      for (size_t r = 0; r < rows_; ++r)
        data_[r + offset] = m.data_[r + offset];
    }
  }
  return *this;
}

Matrix& Matrix::operator=(const MatrixView& m)
{
  if (this != &m){
    assert(rows_ == m.rows() && cols_ == m.cols());
    for (size_t c = 0; c < cols_; ++c) {
      size_t offset = c * stride();
      for (size_t r = 0; r < rows_; ++r)
        data_[r + offset] = m.data()[r + offset];
    }
  }
  return *this;
}

} // namespace bmagwa
