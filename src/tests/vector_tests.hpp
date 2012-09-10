/* BMAGWA software v2.0
 *
 * vector_tests.hpp
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

#ifndef VECTOR_TESTS_HPP_
#define VECTOR_TESTS_HPP_

#include "gtest/gtest.h"
#include "../linalg.hpp"

namespace {

  TEST(VectorTests, Copy)
  {
    bmagwa::Vector v(10);
    v = 1.5;

    bmagwa::VectorView vv(v);
    bmagwa::Vector v2(v);
    bmagwa::Vector v3(vv);

    ASSERT_EQ(vv.length(), v.length());
    ASSERT_EQ(v2.length(), v.length());
    ASSERT_EQ(v3.length(), v.length());

    for (int i = 0; i < 10; ++i){
      ASSERT_EQ(vv(i), v(i));
      ASSERT_EQ(v2(i), v(i));
      ASSERT_EQ(v3(i), v(i));

      v(i) *= 2;
    }

    for (int i = 0; i < 10; ++i){
      ASSERT_EQ(vv(i), v(i));
      ASSERT_NE(v2(i), v(i));
      ASSERT_NE(v3(i), v(i));
    }
  }

  TEST(VectorTests, Block)
  {
    bmagwa::Vector v(4);
    v(0) = 1.1; v(1) = 2.2; v(2) = 3.3; v(3) = 4.4;

    bmagwa::VectorView vv(v.block(1, 2));
    ASSERT_EQ(2, vv.length());
    ASSERT_EQ(vv(0), v(1));
    ASSERT_EQ(vv(1), v(2));
  }

  TEST(VectorTests, SetTo)
  {
    bmagwa::Vector v(10);
    double val = 1.1;

    v.set_to(val);
    for (size_t i = 0; i < v.length(); ++i)
      ASSERT_EQ(val, v(i));

    val = 1.5;
    v.set_mem_to(val);
    for (size_t i = 0; i < v.length_mem(); ++i)
      ASSERT_EQ(val, v(i));
  }

  TEST(VectorTests, DotProduct)
  {
    bmagwa::Vector v(3);
    bmagwa::Vector v2(3);
    v(0) = 1.0; v(1) = 2.0; v(2) = 3.0;
    v2(0) = 4.0; v2(1) = 5.0; v2(2) = 6.0;

    ASSERT_EQ(bmagwa::VectorView::dotproduct(v, v2), 1.0*4.0+2.0*5.0+3.0*6.0);
  }
}


#endif /* VECTOR_TESTS_HPP_ */
