 /* BMAGWA software v2.0
 *
 * data_tests.hpp
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

#ifndef DATA_TESTS_HPP_
#define DATA_TESTS_HPP_

#include "gtest/gtest.h"
#include "../data.hpp"

namespace {

  TEST(DataTests, LoadAndHandlePlinkDataNoMissing)
  {
    const int n = 30;
    const int m = 7;
    bmagwa::Data* data = new bmagwa::Data(n, m, 0,
                            "testdata/small_modelspace.fam",
                            "testdata/small_modelspace.bed",
                            false,
                            "testdata/small_modelspace.e");

    const char genotypes[7][30] =
    {{2,1,0,2,1,2,2,1,2,2,1,0,2,1,1,0,2,2,1,2,2,2,1,2,2,2,2,2,2,1},
     {1,2,2,2,2,1,2,2,1,2,2,2,1,2,2,2,2,2,0,2,2,2,2,2,2,2,2,1,1,2},
     {1,1,1,1,1,1,1,2,1,2,1,1,1,1,1,0,0,0,0,1,1,0,1,2,2,2,1,2,2,2},
     {2,2,2,2,2,2,2,2,1,2,2,2,0,2,2,1,2,2,2,1,2,2,2,1,1,2,2,2,2,2},
     {2,0,2,0,1,2,2,1,1,1,1,1,1,0,1,1,0,2,1,1,1,2,2,0,2,1,2,1,1,2},
     {2,1,2,2,2,1,0,1,0,2,2,2,0,2,2,2,1,2,2,2,2,1,1,0,2,2,2,2,2,1},
     {0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,2,1,2,1,2,2,1,1,1,0,2,1}};

    const double y[30] = {
        0.290913,
        1.72566,
        -1.14075,
         -0.775,
        0.283488,
        -0.0173902,
        0.359494,
        0.91743,
        1.17233,
        0.0284589,
        -0.3621,
        -0.63587,
        -0.518061,
        -0.175326,
        -0.0953311,
        1.88926,
        -0.311815,
        -0.935637,
        0.342461,
        0.401009,
        0.348721,
        1.65077,
        1.53112,
        2.58911,
        -0.162152,
        -0.196501,
        -0.719778,
        0.384602,
        1.43458,
        0.207898};


    for (int i = 0; i < n; ++i){
      for (int j = 0; j < m; ++j){
        ASSERT_EQ(data->get_genotype(i, j), genotypes[j][i]);
      }
    }

    for (int j = 0; j < m; ++j){
      bmagwa::Vector snp(n);
      data->get_genotypes_additive(j, snp);
      for (int i = 0; i < n; ++i){
        ASSERT_EQ(snp(i), genotypes[j][i]);
      }

      data->get_genotypes_recessive(j, snp);
      for (int i = 0; i < n; ++i){
        ASSERT_EQ(snp(i), genotypes[j][i] == 2);
      }

      data->get_genotypes_dominant(j, snp);
      for (int i = 0; i < n; ++i){
        ASSERT_EQ(snp(i), genotypes[j][i] > 0);
      }

      data->get_genotypes_heterozygous(j, snp);
      for (int i = 0; i < n; ++i){
        ASSERT_EQ(snp(i), genotypes[j][i] == 1);
      }
    }

    for (int i = 0; i < n; ++i){
      ASSERT_EQ(data->y()(i), y[i]);
    }


    delete data;
  }


  TEST(DataTests, LoadAndHandlePlinkDataWithMissing)
  {
    const int n = 5;
    const int m = 10;
    bmagwa::Data* data = new bmagwa::Data(n, m, 0,
                                "testdata/plinktest.fam",
                                "testdata/plinktest.bed",
                                false,
                                "testdata/plinktest.e");

    const char genotypes[10][5] =
    {{1,2,2,-1,0},
     {1,1,1,2,0},
     {1,1,2,1,1},
     {0,-1,0,1,2},
     {0,2,2,2,2},
     {1,2,1,1,0},
     {2,1,0,2,1},
     {1,0,1,2,1},
     {1,0,1,1,0},
     {0,1,1,1,1}};

    for (int i = 0; i < n; ++i){
      for (int j = 0; j < m; ++j){
        ASSERT_EQ(data->get_genotype(i, j), genotypes[j][i]);
      }
    }

    for (int j = 0; j < m; ++j){
      bmagwa::Vector snp(n);
      data->get_genotypes_additive(j, snp);
      for (int i = 0; i < n; ++i){
        ASSERT_EQ(snp(i), genotypes[j][i]);
      }

      data->get_genotypes_recessive(j, snp);
      for (int i = 0; i < n; ++i){
        ASSERT_EQ(snp(i), (genotypes[j][i] == -1) ? -1 : (genotypes[j][i] == 2));
      }

      data->get_genotypes_dominant(j, snp);
      for (int i = 0; i < n; ++i){
        ASSERT_EQ(snp(i), (genotypes[j][i] == -1) ? -1 : (genotypes[j][i] > 0));
      }

      data->get_genotypes_heterozygous(j, snp);
      for (int i = 0; i < n; ++i){
        ASSERT_EQ(snp(i), (genotypes[j][i] == -1) ? -1 : (genotypes[j][i] == 1));
      }
    }

    // check missing handling
    for (int j = 0; j < 7; ++j){
      if (j == 0 || j == 3){
        ASSERT_EQ(data->miss_loc()[j][0], 1); // single missing value
        if (j == 0){
          ASSERT_EQ(data->miss_loc()[j][1], 3);       // at position 3
          ASSERT_EQ(data->miss_prior()[j][0], 1);     // cumulative counts
          ASSERT_EQ(data->miss_prior()[j][1], 1+1);
          ASSERT_EQ(data->miss_prior()[j][2], 1+1+2);
        }
        if (j == 3){
          ASSERT_EQ(data->miss_loc()[j][1], 1); // at position 1
          ASSERT_EQ(data->miss_prior()[j][0], 2);     // cumulative counts
          ASSERT_EQ(data->miss_prior()[j][1], 2+1);
          ASSERT_EQ(data->miss_prior()[j][2], 2+1+1);
        }
      } else {
        ASSERT_TRUE(data->miss_loc()[j] == NULL);
        ASSERT_TRUE(data->miss_prior()[j] == NULL);
      }
    }

    delete data;
  }


}

#endif /* DATA_TESTS_HPP_ */
