/* BMAGWA software v1.0
 *
 * symmmatrix.cpp
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


#include "symmmatrix.hpp"
#include "vector.hpp"
#include "matrix.hpp"

namespace bmagwa {

SymmMatrixView& SymmMatrixView::operator=(const SymmMatrixView &m)
{
  if (this != &m){
    assert(length_ == m.length_);
    for (size_t c = 0; c < length_; ++c) {
      size_t offset = c * stride();
      for (size_t r = 0; r <= c; ++r)
        data_[r + offset] = m.data_[r + offset];
    }
  }
  return *this;
}

VectorView SymmMatrixView::column(const size_t col)
{
  assert(col < length());
  return VectorView(data_ + col * stride(), length_mem_, length_);
}


void SymmMatrixView::set_to(const double val)
{
  for (size_t c = 0; c < length_; ++c) {
    size_t offset = c * stride();
    for (size_t r = 0; r < length_; ++r)
      data_[r + offset] = val;
  }
}

void SymmMatrixView::set_mem_to(const double val)
{
  for (size_t i = 0; i < length_mem_ * length_mem_; ++i) {
    data_[i] = val;
  }
}

void SymmMatrixView::set_to_identity()
{
  for (size_t c = 0; c < length_; ++c) {
    size_t offset = c * stride();
    for (size_t r = 0; r < length_; ++r)
      data_[r + offset] = (r == c) ? 1.0 : 0.0;
  }
}

void SymmMatrixView::set_to_innerproduct(const MatrixView &m)
{
  assert(length() == m.cols());
  cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, length(), m.rows(),
              1.0, m.data_, m.stride(), 0.0, data_, stride());
}


SymmMatrix& SymmMatrix::operator=(const SymmMatrix& m)
{
  if (this != &m){
    assert(length_ == m.length_);
    for (size_t c = 0; c < length_; ++c) {
      size_t offset = c * stride();
      for (size_t r = 0; r <= c; ++r)
        data_[r + offset] = m.data_[r + offset];
    }
  }
  return *this;
}

SymmMatrix& SymmMatrix::operator=(const SymmMatrixView& m)
{
  if (this != &m){
    assert(length_ == m.length_);
    for (size_t c = 0; c < length_; ++c) {
      size_t offset = c * stride();
      for (size_t r = 0; r <= c; ++r)
        data_[r + offset] = m.data_[r + offset];
    }
  }
  return *this;
}

void SymmMatrix::remove_colrow(const size_t col_rem)
{
  assert(col_rem < length());

  for (size_t i = col_rem; i < (length() - 1); ++i){
    for(size_t j = 0; j <= i; ++j){
      if (j < col_rem){
        // copy from one col right
        (*this)(j, i) = (*this)(j, i + 1);
      } else {
        // copy from one col right and one row down
        (*this)(j, i) = (*this)(j + 1, i + 1);
      }
    }
  }
  resize(length() - 1);
}

bool SymmMatrix::cholesky_update(const VectorView& newcol,
                                     const double inv_tau2_value)
{
  size_t col = length();
  assert(col + 1 == newcol.length());

  resize(col + 1);
  VectorView lcol(column(col));
  lcol = newcol;
  lcol.resize(col); // drop last temporarily
  resize(col);    // drop last temporarily
  lcol.multiply_by_invtriangularmatrix(*this, true);

  double lmm2 = newcol(col) + inv_tau2_value
                - Vector::dotproduct(lcol, lcol);

  resize(col + 1);
  if (lmm2 <= 0) return false; // not pos.def.
  (*this)(col, col) = sqrt(lmm2);
  return true;
}

void SymmMatrix::cholesky_downdate(const size_t col_rem,
                                       double* c_tmp, double* s_tmp)
{
  // downdate Cholesky
  // l should be of size m x m on entry
  int col_rem_f = col_rem + 1; // fortran indexing starts from 1
  int offset_f = (int)stride();
  int m = (int)length();
  int tmp = 0, job = 2;

  // dchex_ moves a single row and column to a new place in the Cholesky
  // decomposition, which is here used to move the removed variable to last...
  dchex_(data_, &offset_f, &m, &col_rem_f, &m, 0, &tmp, &tmp, c_tmp, s_tmp,
         &job);

  // ... and then the last is just dropped
  resize(m - 1);

  // fix to the dchex bug
  // (ref. M.Seeger tech. report from 2008:
  //  Low Rank Updates for the Cholesky Decomposition)
  // don't know if this happens but I dare not remove this...
  for (size_t i = col_rem; i < length(); ++i){
    if ((*this)(i, i) < 0){
      for (size_t j = i; j < length(); ++j){
        (*this)(i, j) *= -1;
      }
    }
  }
}

} // namespace bmagwa
