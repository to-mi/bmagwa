/* BMAGWA software v2.0
 *
 * main.cpp
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

#include "gtest/gtest.h"

// tests
#include "rand_tests.hpp"
#include "vector_tests.hpp"
#include "matrix_tests.hpp"
#include "symmmatrix_tests.hpp"
#include "data_tests.hpp"
#include "data_model_tests.hpp"
#include "model_tests.hpp"
#include "options_tests.hpp"
#include "utils_tests.hpp"
#include "prior_tests.hpp"
#include "sampler_tests.hpp"
#include "discrete_distribution_tests.hpp"
#include "precomputed_snp_covariances_tests.hpp"
#include "raoblackwellizer_tests.hpp"


int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
