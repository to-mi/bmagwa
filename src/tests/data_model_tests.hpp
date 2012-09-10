/* BMAGWA software v2.0
 *
 * data_model_tests.hpp
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

#ifndef DATA_MODEL_TESTS_HPP_
#define DATA_MODEL_TESTS_HPP_

#include <cstring>
#include "gtest/gtest.h"
#include "../data_model.hpp"
#include "../precomputed_snp_covariances.hpp"

namespace {

  TEST(DataModelTests, Initialization) {
    bmagwa::Data data(5, 10, 0, "testdata/plinktest.fam", "testdata/plinktest.bed",
                      false, "testdata/test.e");
    std::vector<bmagwa::DataModel::ef_t> types;
    types.push_back(bmagwa::DataModel::A);
    types.push_back(bmagwa::DataModel::AH);

    bmagwa::DataModel data_model(&data, types);

    ASSERT_EQ(data_model.types.size(), types.size());
    for (size_t i = 0; i < types.size(); ++i)
      EXPECT_EQ(data_model.types[i], types[i]);

    ASSERT_TRUE(data_model.allow_terms(bmagwa::DataModel::A));
    ASSERT_TRUE(data_model.allow_terms(bmagwa::DataModel::H));
    ASSERT_FALSE(data_model.allow_terms(bmagwa::DataModel::D));
    ASSERT_FALSE(data_model.allow_terms(bmagwa::DataModel::R));

    ASSERT_TRUE(data_model.allow_types(bmagwa::DataModel::A));
    ASSERT_FALSE(data_model.allow_types(bmagwa::DataModel::H));
    ASSERT_FALSE(data_model.allow_types(bmagwa::DataModel::D));
    ASSERT_FALSE(data_model.allow_types(bmagwa::DataModel::R));
    ASSERT_TRUE(data_model.allow_types(bmagwa::DataModel::AH));
  }


  TEST(DataModelTests, GetGenotypes) {
    const int n = 5;
    const int m = 10;

    bmagwa::Data data(n, m, 0, "testdata/plinktest.fam", "testdata/plinktest.bed",
                      false, "testdata/test.e");
    std::vector<bmagwa::DataModel::ef_t> types;
    types.push_back(bmagwa::DataModel::A);
    types.push_back(bmagwa::DataModel::AH);

    bmagwa::DataModel data_model(&data, types);

    // genotypes with zeroed missing values
    const char genotypes[10][5] =
    {{1,2,2,0,0},
     {1,1,1,2,0},
     {1,1,2,1,1},
     {0,0,0,1,2},
     {0,2,2,2,2},
     {1,2,1,1,0},
     {2,1,0,2,1},
     {1,0,1,2,1},
     {1,0,1,1,0},
     {0,1,1,1,1}};

    for (int j = 0; j < m; ++j){
      bmagwa::Vector snp(n);
      data_model.get_genotypes_additive(j, snp);
      for (int i = 0; i < n; ++i){
        ASSERT_EQ(snp(i), genotypes[j][i]);
      }
    }

    // additive
    for (int j = 0; j < m; ++j){
      bmagwa::Vector snp(n);
      (data_model.*data_model.bmagwa::DataModel::get_genotypes[bmagwa::DataModel::A])(j, snp);
      for (int i = 0; i < n; ++i){
        ASSERT_EQ(snp(i), genotypes[j][i]);
      }
    }

    // heterozygous
    for (int j = 0; j < m; ++j){
      bmagwa::Vector snp(n);
      (data_model.*data_model.bmagwa::DataModel::get_genotypes[bmagwa::DataModel::H])(j, snp);
      for (int i = 0; i < n; ++i){
        ASSERT_EQ(snp(i), genotypes[j][i] == 1);
      }
    }

    // dominant
    for (int j = 0; j < m; ++j){
      bmagwa::Vector snp(n);
      (data_model.*data_model.bmagwa::DataModel::get_genotypes[bmagwa::DataModel::D])(j, snp);
      for (int i = 0; i < n; ++i){
        ASSERT_EQ(snp(i), genotypes[j][i] > 0);
      }
    }

    // recessive
    for (int j = 0; j < m; ++j){
      bmagwa::Vector snp(n);
      (data_model.*data_model.bmagwa::DataModel::get_genotypes[bmagwa::DataModel::R])(j, snp);
      for (int i = 0; i < n; ++i){
        ASSERT_EQ(snp(i), genotypes[j][i] == 2);
      }
    }
  }


  TEST(DataModelTests, MissingSampling) {
    size_t replications = 100000;
    double threshold = 1e-2;

    bmagwa::Data data(5, 10, 0, "testdata/plinktest.fam", "testdata/plinktest.bed",
                      false, "testdata/test.e");
    std::vector<bmagwa::DataModel::ef_t> types;
    types.push_back(bmagwa::DataModel::A);
    types.push_back(bmagwa::DataModel::AH);

    bmagwa::Rand rng(1234, 100);
    bmagwa::DataModel data_model(&data, types);

    bmagwa::Vector freqs(3);

    bmagwa::Vector v(5);
    data.get_genotypes_additive(0, v);
    ASSERT_EQ(-1.0, v(3));
    data_model.get_genotypes_additive(0, v);
    ASSERT_EQ(0.0, v(3));

    // test sample_missing_single
    freqs.set_to(0);
    for (size_t r = 0; r < replications; ++r){
      data_model.sample_missing_single(0, rng);
      data_model.get_genotypes_additive(0, v);
      ASSERT_LE(0, v(3));
      freqs(v(3)) += 1;
    }

    double sum = freqs.sum();
    freqs(0) = freqs(0) / sum;
    freqs(1) = freqs(1) / sum;
    freqs(2) = freqs(2) / sum;

    EXPECT_NEAR(0.25, freqs(0), threshold);
    EXPECT_NEAR(0.25, freqs(1), threshold);
    EXPECT_NEAR(0.5, freqs(2), threshold);

    // test sample_missing
    std::vector<int32_t> model_inds(10, -1);
    v.set_to(-1.0);
    freqs.set_to(0);
    for (size_t r = 0; r < replications; ++r){
      data_model.sample_missing(model_inds, rng);
      data_model.get_genotypes_additive(0, v);
      ASSERT_LE(0, v(3));
      freqs(v(3)) += 1;
    }

    sum = freqs.sum();
    freqs(0) = freqs(0) / sum;
    freqs(1) = freqs(1) / sum;
    freqs(2) = freqs(2) / sum;

    EXPECT_NEAR(0.25, freqs(0), threshold);
    EXPECT_NEAR(0.25, freqs(1), threshold);
    EXPECT_NEAR(0.5, freqs(2), threshold);
  }


  TEST(DataModelTests, XXUpdatesA)
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

    for (int r = 0; r < 100; ++r){
      data_model.sample_missing_single(0, rng);

      for (int j = 0; j < 2; ++j){
        int jj = (j == 0) ? 0 : 3;

        bmagwa::Vector snp(n);
        double xx[2];
        xx[0] = pre_xxcov.xx[pre_xxcov.offset * jj]; // n*mean of first snp
        xx[1] = pre_xxcov.xx[pre_xxcov.offset * jj + 1]; // (n-1)*var of first snp

        data_model.update_prexx_cov(jj, xx);
        data_model.get_genotypes_additive(jj, snp);

        ASSERT_DOUBLE_EQ(xx[0], n * snp.mean());
        ASSERT_DOUBLE_EQ(xx[1], (n-1) * snp.var());
      }
    }
  }


  TEST(DataModelTests, XXUpdatesAAH)
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

    for (int r = 0; r < 100; ++r){
      data_model.sample_missing_single(0, rng);

      for (int j = 0; j < 2; ++j){
        int jj = (j == 0) ? 0 : 3;

        bmagwa::Vector snp_a(n);
        bmagwa::Vector snp_h(n);
        double xx[5];
        xx[0] = pre_xxcov.xx[pre_xxcov.offset * jj];
        xx[1] = pre_xxcov.xx[pre_xxcov.offset * jj + 1];
        xx[2] = pre_xxcov.xx[pre_xxcov.offset * jj + 2];
        xx[3] = pre_xxcov.xx[pre_xxcov.offset * jj + 3];
        xx[4] = pre_xxcov.xx[pre_xxcov.offset * jj + 4];

        data_model.update_prexx_cov(jj, xx);
        data_model.get_genotypes_additive(jj, snp_a);
        data_model.get_genotypes_heterozygous(jj, snp_h);

        ASSERT_DOUBLE_EQ(xx[0], n * snp_a.mean());
        ASSERT_DOUBLE_EQ(xx[1], (n-1) * snp_a.var());
        ASSERT_DOUBLE_EQ(xx[2], n * snp_h.mean());
        ASSERT_DOUBLE_EQ(xx[3], (n-1) * snp_h.var());
        EXPECT_NEAR(xx[4], bmagwa::VectorView::dotproduct(snp_a, snp_h) - n * snp_a.mean() * snp_h.mean(), threshold);
      }
    }
  }


  TEST(DataModelTests, XXUpdatesAll)
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

    for (int r = 0; r < 100; ++r){
      data_model.sample_missing_single(0, rng);

      for (int j = 0; j < 2; ++j){
        int jj = (j == 0) ? 0 : 3;

        bmagwa::Vector snp_a(n);
        bmagwa::Vector snp_h(n);
        bmagwa::Vector snp_d(n);
        bmagwa::Vector snp_r(n);
        double xx[9];
        for (int i = 0; i < 9; ++i)
          xx[i] = pre_xxcov.xx[pre_xxcov.offset * jj + i];

        data_model.update_prexx_cov(jj, xx);
        data_model.get_genotypes_additive(jj, snp_a);
        data_model.get_genotypes_heterozygous(jj, snp_h);
        data_model.get_genotypes_dominant(jj, snp_d);
        data_model.get_genotypes_recessive(jj, snp_r);

        ASSERT_DOUBLE_EQ(xx[0], n * snp_a.mean());
        ASSERT_DOUBLE_EQ(xx[1], (n-1) * snp_a.var());
        ASSERT_DOUBLE_EQ(xx[2], n * snp_h.mean());
        ASSERT_DOUBLE_EQ(xx[3], (n-1) * snp_h.var());
        ASSERT_DOUBLE_EQ(xx[4], n * snp_d.mean());
        ASSERT_DOUBLE_EQ(xx[5], (n-1) * snp_d.var());
        ASSERT_DOUBLE_EQ(xx[6], n * snp_r.mean());
        ASSERT_DOUBLE_EQ(xx[7], (n-1) * snp_r.var());
        EXPECT_NEAR(xx[8], bmagwa::VectorView::dotproduct(snp_a, snp_h) - n * snp_a.mean() * snp_h.mean(), threshold);
      }
    }
  }

}

#endif /* DATA_MODEL_TESTS_HPP_ */
