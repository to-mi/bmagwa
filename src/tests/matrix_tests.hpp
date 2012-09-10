/* BMAGWA software v2.0
 *
 * matrix_tests.hpp
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

#ifndef MATRIX_TESTS_HPP_
#define MATRIX_TESTS_HPP_

#include "gtest/gtest.h"
#include "../linalg.hpp"

namespace {

  TEST(MatrixTests, SetTo) {
    bmagwa::Matrix m(1000, 2);
    double val = 1.1;

    m.set_to(val);
    for (size_t c = 0; c < m.cols(); ++c)
      for (size_t r = 0; r < m.rows(); ++r)
        ASSERT_EQ(val, m(r, c));
  }

}

#endif /* MATRIX_TESTS_HPP_ */
