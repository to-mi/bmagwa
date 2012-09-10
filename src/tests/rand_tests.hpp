/* BMAGWA software v2.0
 *
 * rand_tests.hpp
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

#ifndef RAND_TESTS_HPP_
#define RAND_TESTS_HPP_

#include "gtest/gtest.h"
#include "../rand.hpp"
#include "../linalg.hpp"

namespace {
  TEST(RandTests, SInvChi2) {
    size_t nsamples = 100000000;
    double threshold_var = 1e-3;
    double threshold_mu = 1e-5;
    //double nu = 10, s2 = 5.5;
    double nu = 5, s2 = 0.0136;
    double mean = nu * s2 / (nu - 2);
    double var = 2 * nu * nu * s2 * s2 / ((nu - 2) * (nu - 2) * (nu - 4));
    bmagwa::Rand rand(1234, nu);

    {
      bmagwa::Vector v(nsamples);
      for (size_t i = 0; i < nsamples; ++i)
        v(i) = rand.rand_sinvchi2(s2);
      EXPECT_NEAR(mean, v.mean(), threshold_mu);
      EXPECT_NEAR(v.var(), var, threshold_var);
    }

    {
      bmagwa::Vector v(nsamples);
      for (size_t i = 0; i < nsamples; ++i)
        v(i) = rand.rand_sinvchi2(nu, s2);
      EXPECT_NEAR(mean, v.mean(), threshold_mu);
      EXPECT_NEAR(v.var(), var, threshold_var);
    }
  }

  TEST(RandTests, MVNormalWithIdentitySigma) {
    double threshold = 1e-2;
    double mu_val = 1.1; double scale = 0.1;
    bmagwa::Rand rand(1234, 100);

    bmagwa::SymmMatrix m(1000);
    m.set_to_identity();

    bmagwa::Vector v(1000);
    bmagwa::Vector mu(1000);
    mu = mu_val;

    rand.rand_mvnormal_chol(mu, m, scale, v);

    EXPECT_NEAR(mu_val, v.mean(), threshold);
    EXPECT_NEAR(scale * scale, v.var(), threshold);

    v.set_to(5.5);
    rand.rand_mvnormal_invchol(mu, m, scale, v);
    EXPECT_NEAR(mu_val, v.mean(), threshold);
    EXPECT_NEAR(scale * scale, v.var(), threshold);
  }

  TEST(RandTests, MVNormal2D) {
    double threshold = 1e-2; size_t n = 1000000;
    bmagwa::Rand rand(1234, 100);

    // m is cholesky of covariance matrix
    bmagwa::SymmMatrix m(2);
    m(0,0) = 1; m(0,1) = 0.9; m(1, 1) = 0.4359;

    bmagwa::Matrix v(n, 2);
    bmagwa::Vector mu(2);
    mu(0) = 1.1; mu(1) = 2.1;

    // test with chol(Sigma)
    for (size_t i = 0; i < n; ++i) {
      bmagwa::Vector tmp(2);
      rand.rand_mvnormal_chol(mu, m, 1.0, tmp);
      v(i, 0) = tmp(0); v(i, 1) = tmp(1);
    }
    double covar = bmagwa::VectorView::dotproduct(
                                    v.column(0) - v.column(0).mean(),
                                    v.column(1) - v.column(1).mean()) / (n - 1);
    EXPECT_NEAR(0.9, covar, threshold);

    // test with chol(Sigma^-1)
    v.set_to(0.0);
    m(0,0) = 2.2942; m(0,1) = -2.0647; m(1, 1) = 1.000;
    for (size_t i = 0; i < n; ++i) {
      bmagwa::Vector tmp(2);
      rand.rand_mvnormal_invchol(mu, m, 1.0, tmp);
      v(i, 0) = tmp(0); v(i, 1) = tmp(1);
    }
    covar = bmagwa::VectorView::dotproduct(v.column(0) - v.column(0).mean(),
                                           v.column(1) - v.column(1).mean()
                                           ) / (n - 1);
    EXPECT_NEAR(0.9, covar, threshold);
  }

  TEST(RandTests, MVNormal3D) {
    double threshold = 1e-2; size_t n = 1000000;
    bmagwa::Rand rand(1234, 100);

    bmagwa::SymmMatrix m(3);

    // m is cholesky of covariance matrix
    m(0,0) = 1; m(0, 1) = 0.9;    m(0, 2) = 0.5;
                m(1, 1) = 0.4359; m(1, 2) = -0.3441;
                                  m(2, 2) = 0.7947;


    bmagwa::Matrix v(n, 3);
    bmagwa::Vector mu(3);
    mu(0) = 1.1; mu(1) = 2.1; mu(2) = 0.5;
    double scale = 0.6;

    // test with chol(Sigma)
    for (size_t i = 0; i < n; ++i) {
      bmagwa::Vector tmp(3);
      rand.rand_mvnormal_chol(mu, m, scale, tmp);
      v(i, 0) = tmp(0); v(i, 1) = tmp(1); v(i, 2) = tmp(2);
    }
    double covar12 = bmagwa::VectorView::dotproduct(
                                    v.column(0) - v.column(0).mean(),
                                    v.column(1) - v.column(1).mean()) / (n - 1);
    double covar13 = bmagwa::VectorView::dotproduct(
                                    v.column(0) - v.column(0).mean(),
                                    v.column(2) - v.column(2).mean()) / (n - 1);
    double covar23 = bmagwa::VectorView::dotproduct(
                                    v.column(1) - v.column(1).mean(),
                                    v.column(2) - v.column(2).mean()) / (n - 1);

    EXPECT_NEAR(scale * scale * 0.9, covar12, threshold);
    EXPECT_NEAR(scale * scale * 0.5, covar13, threshold);
    EXPECT_NEAR(scale * scale * 0.3, covar23, threshold);

    // test with chol(Sigma^-1)
    v.set_to(0.0);
    m(0,0) = 2.7538; m(0, 1) = -2.2696; m(0, 2) = -0.6960;
                      m(1, 1) = 1.0483; m(1, 2) = -0.3145;
                                        m(2, 2) = 1.0000;
    for (size_t i = 0; i < n; ++i) {
      bmagwa::Vector tmp(3);
      rand.rand_mvnormal_invchol(mu, m, scale, tmp);
      v(i, 0) = tmp(0); v(i, 1) = tmp(1); v(i, 2) = tmp(2);
    }
    covar12 = bmagwa::VectorView::dotproduct(v.column(0) - v.column(0).mean(),
                                             v.column(1) - v.column(1).mean()
                                             ) / (n - 1);
    covar13 = bmagwa::VectorView::dotproduct(v.column(0) - v.column(0).mean(),
                                             v.column(2) - v.column(2).mean()
                                             ) / (n - 1);
    covar23 = bmagwa::VectorView::dotproduct(v.column(1) - v.column(1).mean(),
                                             v.column(2) - v.column(2).mean()
                                             ) / (n - 1);

    EXPECT_NEAR(scale * scale * 0.9, covar12, threshold);
    EXPECT_NEAR(scale * scale * 0.5, covar13, threshold);
    EXPECT_NEAR(scale * scale * 0.3, covar23, threshold);
  }
}

#endif /* RAND_TESTS_HPP_ */
