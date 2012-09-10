/* BMAGWA software v2.0
 *
 * symmmatrix_tests.hpp
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

#ifndef SYMMMATRIX_TESTS_HPP_
#define SYMMMATRIX_TESTS_HPP_

#include "gtest/gtest.h"
#include "../linalg.hpp"
#include "../rand.hpp"

namespace {

  TEST(SymmMatrixTests, SetTo) {
    bmagwa::SymmMatrix m(1000);
    double val = 1.1;

    m.set_to(val);
    for (size_t c = 0; c < m.cols(); ++c)
      for (size_t r = 0; r < m.rows(); ++r)
        ASSERT_EQ(val, m(r, c));
  }

  TEST(SymmMatrixTests, Cholesky) {
    double threshold = 1e-3;
    bmagwa::SymmMatrix m(10);

    m.set_to_identity();

    ASSERT_TRUE(m.cholesky());
    for (size_t c = 0; c < m.cols(); ++c)
      for (size_t r = 0; r <= c; ++r)
        ASSERT_EQ((c == r) ? 1.0 : 0.0, m(r, c));

    m.resize(2);
    m(0, 0) = 1; m(0, 1) = 0.9; m(1, 1) = 1;
    ASSERT_TRUE(m.cholesky());
    EXPECT_NEAR(1, m(0, 0), threshold);
    EXPECT_NEAR(0.9, m(0, 1), threshold);
    EXPECT_NEAR(0.4359, m(1, 1), threshold);

    m(0, 0) = 1; m(0, 1) = 1; m(1, 1) = 1;
    ASSERT_FALSE(m.cholesky());

    m(0, 0) = 1; m(0, 1) = 1.1; m(1, 1) = 1;
    ASSERT_FALSE(m.cholesky());
  }

  TEST(SymmMatrixTests, CholeskyUpdate) {
    double threshold = 1e-10;
    size_t nreplications = 100;
    size_t n = 100; size_t m = 10;

    bmagwa::Rand rng(1234, 100);

    for (size_t rep = 0; rep < nreplications; ++rep){
      bmagwa::Matrix x(n, m);
      for (size_t r = 0; r < n; ++r)
        for (size_t c = 0; c < m; ++c)
          x(r, c) = 10 * rng.rand_normal() + 1;

      bmagwa::SymmMatrix xx(m);
      bmagwa::SymmMatrix uchol(m);
      bmagwa::SymmMatrix uchol_upd(m);


      xx.set_to_innerproduct(x);
      uchol = xx;
      if (!uchol.cholesky())
        continue;

      xx.resize(m-1);
      uchol_upd.resize(m-1);
      uchol_upd = xx;
      uchol_upd.cholesky();

      xx.resize(m);
      bmagwa::VectorView xxcol(xx.column(m-1));
      ASSERT_TRUE(uchol_upd.cholesky_update(xxcol, 0.0));

      for (size_t c = 0; c < uchol_upd.length(); ++c)
        for (size_t r = 0; r <= c; ++r){
          //EXPECT_DOUBLE_EQ(uchol(r, c), uchol_upd(r, c));
          EXPECT_NEAR(uchol(r, c), uchol_upd(r, c), threshold);
        }
    }
  }


  TEST(SymmMatrixTests, CholeskyDowndate) {
    double threshold = 1e-10;
    size_t nreplications = 100;
    size_t n = 100; size_t m = 10;
    double c_tmp[m], s_tmp[m];

    bmagwa::Rand rng(1234, 100);

    for (size_t rep = 0; rep < nreplications; ++rep){
      bmagwa::Matrix x(n, m);
      for (size_t r = 0; r < n; ++r)
        for (size_t c = 0; c < m; ++c)
          x(r, c) = 10 * rng.rand_normal() + 1;

      bmagwa::SymmMatrix xx(m);
      bmagwa::SymmMatrix uchol(m - 1);
      bmagwa::SymmMatrix uchol_dow(m);

      size_t col_rem = floor(m * rng.rand_01());

      xx.set_to_innerproduct(x);
      uchol_dow = xx;
      if (!uchol_dow.cholesky())
        continue;

      uchol_dow.cholesky_downdate(col_rem, c_tmp, s_tmp);

      xx.remove_colrow(col_rem);
      uchol = xx;
      ASSERT_TRUE(uchol.cholesky());

      for (size_t c = 0; c < uchol_dow.length(); ++c)
        for (size_t r = 0; r <= c; ++r){
          //EXPECT_DOUBLE_EQ(uchol(r, c), uchol_dow(r, c));
          EXPECT_NEAR(uchol(r, c), uchol_dow(r, c), threshold);
        }
    }

  }

  TEST(SymmMatrixTests, CholeskySwapadj) {
    double threshold = 1e-10;
    size_t nreplications = 100;
    size_t n = 100; size_t m = 10;

    bmagwa::Rand rng(1234, 100);

    for (size_t rep = 0; rep < nreplications; ++rep){
      bmagwa::Matrix x(n, m);
      for (size_t r = 0; r < n; ++r)
        for (size_t c = 0; c < m; ++c)
          x(r, c) = 10 * rng.rand_normal() + 1;

      bmagwa::SymmMatrix xx(m);
      bmagwa::SymmMatrix uchol(m);
      bmagwa::SymmMatrix uchol2(m);

      size_t col_swap = round((m - 2) * rng.rand_01());

      xx.set_to_innerproduct(x);
      uchol = xx;
      if (!uchol.cholesky())
        continue;

      uchol.cholesky_swapadj(col_swap);

      // swap cols/rows in xx
      xx.swap_colrow(col_swap);

      uchol2 = xx;
      ASSERT_TRUE(uchol2.cholesky());

      for (size_t c = 0; c < uchol.length(); ++c)
        for (size_t r = 0; r <= c; ++r){
          //EXPECT_DOUBLE_EQ(uchol(r, c), uchol_dow(r, c));
          //std::cout << r << "," << c << " " << uchol(r, c) << " " << uchol2(r, c) << std::endl;
          EXPECT_NEAR(uchol(r, c), uchol2(r, c), threshold) << r << ", " << c;
        }
    }

    {
      // this gives negative r in chol_swapadj Givens rotation
      bmagwa::SymmMatrix xx(2);
      bmagwa::SymmMatrix chol(2);

      xx(0, 0) = 100.0;
      xx(0, 1) = -50.0;
      xx(1, 1) = 26.0;

      chol = xx;

      ASSERT_TRUE(chol.cholesky());
      chol.cholesky_swapadj(0);

      ASSERT_TRUE(chol(0,0) > 0);
      ASSERT_TRUE(chol(1,1) > 0);

      xx.swap_colrow(0);
      xx.cholesky();
      EXPECT_NEAR(xx(0, 0), chol(0, 0), threshold);
      EXPECT_NEAR(xx(0, 1), chol(0, 1), threshold);
      EXPECT_NEAR(xx(1, 1), chol(1, 1), threshold);
    }

    {
      // this gives positive r in chol_swapadj Givens rotation
      bmagwa::SymmMatrix xx(2);
      bmagwa::SymmMatrix chol(2);

      xx(0, 0) = 100.0;
      xx(0, 1) = 50.0;
      xx(1, 1) = 80.0;

      chol = xx;

      ASSERT_TRUE(chol.cholesky());
      chol.cholesky_swapadj(0);

      ASSERT_TRUE(chol(0,0) > 0);
      ASSERT_TRUE(chol(1,1) > 0);

      xx.swap_colrow(0);
      xx.cholesky();
      EXPECT_NEAR(xx(0, 0), chol(0, 0), threshold);
      EXPECT_NEAR(xx(0, 1), chol(0, 1), threshold);
      EXPECT_NEAR(xx(1, 1), chol(1, 1), threshold);
    }

  }

}

#endif /* SYMMMATRIX_TESTS_HPP_ */
