/* BMAGWA software v2.0
 *
 * options_tests.hpp
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

#ifndef OPTIONS_TESTS_HPP_
#define OPTIONS_TESTS_HPP_

#include <vector>
#include "gtest/gtest.h"
#include "../options.hpp"

namespace {
  TEST(OptionsTests, Initialization) {
    bmagwa::Options options("testdata/plinktest.ini");

    std::vector<bmagwa::Options*> opts;
    for (size_t i = 0; i < options.n_threads; ++i){
      opts.push_back(options.clone(i));
    }



    size_t n = opts.size();
    for (size_t i = 0; i < n; ++i){
      delete opts.back();
      opts.pop_back();
    }

  }
}

#endif /* OPTIONS_TESTS_HPP_ */
