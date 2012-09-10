/* BMAGWA software v2.0
 *
 * utils_tests.hpp
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

#ifndef UTILS_TESTS_HPP_
#define UTILS_TESTS_HPP_

#include "gtest/gtest.h"
#include "../utils.hpp"
#include "../rand.hpp"

namespace {

  TEST(UtilsTests, SampleDiscreteNaive)
  {
    double p_moves_cumsum[4] = { 0.4, 0.8, 1, 1 };
    size_t nsamples = 1000000;
    double threshold = 0.001;

    size_t res[4] = { 0, 0, 0, 0 };
    double res_freq[4] = { 0, 0, 0, 0 };

    bmagwa::Rand rng(1234, 100);

    for (size_t i = 0; i < nsamples; ++i)
      ++res[bmagwa::Utils::sample_discrete_naive(p_moves_cumsum, 4, rng)];

    for (size_t i = 0; i < 4; i++)
      res_freq[i] = (double)res[i] / nsamples;

    EXPECT_NEAR(0.4, res_freq[0], threshold);
    EXPECT_NEAR(0.4, res_freq[1], threshold);
    EXPECT_NEAR(0.2, res_freq[2], threshold);
    EXPECT_EQ(0, res_freq[3]);

    ASSERT_EQ(nsamples, res[0] + res[1] + res[2] + res[3]);

  }

  TEST(UtilsTests, Gammaln)
  {
    EXPECT_DOUBLE_EQ(bmagwa::Utils::gammaln(4.1),
                     log(3.1) + bmagwa::Utils::gammaln(3.1));
  }

  TEST(UtilsTests, GeometricDist)
  {
    bmagwa::Rand rng(1234, 100);

    int maxsize = 10;
    double values[10];
    for (int i = 0; i < 1; ++i){
      double p = rng.rand_01();
      bmagwa::Utils::geometric_dist_cdf(maxsize, p, values);
      ASSERT_DOUBLE_EQ(p, values[0]);

      for (int j = 0; j < maxsize; ++j){
        double val = 0.0;
        for (int k = 0; k <= j; ++k) val += p * std::pow(1 - p, k);
        EXPECT_DOUBLE_EQ(val, values[j]) << "Problem at " << j << " when p = " << p;
      }
    }
  }

}


#endif /* UTILS_TESTS_HPP_ */
