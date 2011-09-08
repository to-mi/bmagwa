/* BMAGWA software v1.0
 *
 * matrix.cpp
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


#include "matrix.hpp"

namespace bmagwa {

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
  for (size_t i = 0; i < cols_mem_ * rows_mem_; ++i) {
    data_[i] = val;
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
    assert(rows_ == m.rows_ && cols_ == m.cols_);
    for (size_t c = 0; c < cols_; ++c) {
      size_t offset = c * stride();
      for (size_t r = 0; r < rows_; ++r)
        data_[r + offset] = m.data_[r + offset];
    }
  }
  return *this;
}

} // namespace bmagwa
