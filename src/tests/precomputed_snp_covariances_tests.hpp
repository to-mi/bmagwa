/* BMAGWA software v2.0
 *
 * precomputed_snp_covariances_tests.hpp
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

#ifndef PRECOMPUTED_SNP_COVARIANCES_TESTS_HPP_
#define PRECOMPUTED_SNP_COVARIANCES_TESTS_HPP_

#include "gtest/gtest.h"
#include "../precomputed_snp_covariances.hpp"

namespace {

  TEST(PrecomputedSNPCovariancesTests, ComputationA)
  {
    const int n = 5;
    const int m = 10;
    bmagwa::Data data(n, m, 0, "testdata/plinktest.fam", "testdata/plinktest.bed",
                      false, "testdata/test.e");
    std::vector<bmagwa::DataModel::ef_t> types;
    types.push_back(bmagwa::DataModel::A);

    bmagwa::Rand rng(1234, 100);
    bmagwa::DataModel data_model(&data, types);

    bmagwa::PrecomputedSNPCovariances pre_xxcov(&data_model);

    for (int j = 0; j < m; ++j){
      bmagwa::Vector snp(n);
      data_model.get_genotypes_additive(j, snp);

      double* xx = pre_xxcov.xx + j * pre_xxcov.offset;

      ASSERT_DOUBLE_EQ(*xx, snp.mean() * n);
      ++xx;
      ASSERT_DOUBLE_EQ(*xx, snp.var() * (n - 1));
    }
  }

  TEST(PrecomputedSNPCovariancesTests, ComputationAAH)
  {
    const double threshold = 1e-10;
    const int n = 5;
    const int m = 10;
    bmagwa::Data data(n, m, 0, "testdata/plinktest.fam", "testdata/plinktest.bed",
                      false, "testdata/test.e");
    std::vector<bmagwa::DataModel::ef_t> types;
    types.push_back(bmagwa::DataModel::A);
    types.push_back(bmagwa::DataModel::AH);

    bmagwa::Rand rng(1234, 100);
    bmagwa::DataModel data_model(&data, types);

    bmagwa::PrecomputedSNPCovariances pre_xxcov(&data_model);

    for (int j = 0; j < m; ++j){
      bmagwa::Vector snp_a(n);
      bmagwa::Vector snp_h(n);
      data_model.get_genotypes_additive(j, snp_a);
      data_model.get_genotypes_heterozygous(j, snp_h);

      double* xx = pre_xxcov.xx + j * pre_xxcov.offset;


      ASSERT_DOUBLE_EQ(*xx, snp_a.mean() * n);
      ++xx;
      ASSERT_DOUBLE_EQ(*xx, snp_a.var() * (n - 1));
      ++xx;

      ASSERT_DOUBLE_EQ(*xx, snp_h.mean() * n);
      ++xx;
      ASSERT_DOUBLE_EQ(*xx, snp_h.var() * (n - 1));
      ++xx;

      EXPECT_NEAR(*xx, bmagwa::VectorView::dotproduct(snp_a, snp_h) - n * snp_a.mean() * snp_h.mean(), threshold);
    }
  }

  TEST(PrecomputedSNPCovariancesTests, ComputationOnlyAll)
  {
    const double threshold = 1e-10;
    const int n = 5;
    const int m = 10;
    bmagwa::Data data(n, m, 0, "testdata/plinktest.fam", "testdata/plinktest.bed",
                      false, "testdata/test.e");
    std::vector<bmagwa::DataModel::ef_t> types;
    types.push_back(bmagwa::DataModel::A);
    types.push_back(bmagwa::DataModel::H);
    types.push_back(bmagwa::DataModel::D);
    types.push_back(bmagwa::DataModel::R);
    types.push_back(bmagwa::DataModel::AH);

    bmagwa::Rand rng(1234, 100);
    bmagwa::DataModel data_model(&data, types);

    bmagwa::PrecomputedSNPCovariances pre_xxcov(&data_model);

    for (int j = 0; j < m; ++j){
      bmagwa::Vector snp_a(n);
      bmagwa::Vector snp_h(n);
      bmagwa::Vector snp_d(n);
      bmagwa::Vector snp_r(n);
      data_model.get_genotypes_additive(j, snp_a);
      data_model.get_genotypes_heterozygous(j, snp_h);
      data_model.get_genotypes_dominant(j, snp_d);
      data_model.get_genotypes_recessive(j, snp_r);


      double* xx = pre_xxcov.xx + j * pre_xxcov.offset;


      ASSERT_DOUBLE_EQ(*xx, snp_a.mean() * n);
      ++xx;
      ASSERT_DOUBLE_EQ(*xx, snp_a.var() * (n - 1));
      ++xx;

      ASSERT_DOUBLE_EQ(*xx, snp_h.mean() * n);
      ++xx;
      ASSERT_DOUBLE_EQ(*xx, snp_h.var() * (n - 1));
      ++xx;

      ASSERT_DOUBLE_EQ(*xx, snp_d.mean() * n);
      ++xx;
      ASSERT_DOUBLE_EQ(*xx, snp_d.var() * (n - 1));
      ++xx;

      ASSERT_DOUBLE_EQ(*xx, snp_r.mean() * n);
      ++xx;
      ASSERT_DOUBLE_EQ(*xx, snp_r.var() * (n - 1));
      ++xx;

      EXPECT_NEAR(*xx, bmagwa::VectorView::dotproduct(snp_a, snp_h) - n * snp_a.mean() * snp_h.mean(), threshold);
    }
  }
}


#endif /* PRECOMPUTED_SNP_COVARIANCES_HPP_ */
