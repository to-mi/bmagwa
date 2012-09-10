/* BMAGWA software v2.0
 *
 * model_tests.hpp
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

#ifndef MODEL_TESTS_HPP_
#define MODEL_TESTS_HPP_

#include "gtest/gtest.h"
#include <vector>
#include "../data.hpp"
#include "../data_model.hpp"
#include "../prior.hpp"
#include "../model.hpp"
#include "../rand.hpp"
#include "../sampler.hpp"
#include "../samplerstats.hpp"

namespace {
  TEST(ModelTests, LikelihoodComputation) {
    double threshold = 0.00000001;

    size_t n = 5;

    bmagwa::Data data(n, 10, 0, "testdata/plinktest.fam",
        "testdata/plinktest.bed", false, "testdata/test.e");
    std::vector<bmagwa::DataModel::ef_t> types;
    types.push_back(bmagwa::DataModel::A);
    types.push_back(bmagwa::DataModel::AH);

    bmagwa::DataModel data_model(&data, types);
    double types_prior[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
    bmagwa::Vector inv_tau2_e(1); inv_tau2_e = 0;
    double nu_sigma = 4, s2_sigma = 1;
    double nu_tau2[] = {4, 4, 4, 4};
    double s2_tau2[] = {1, 1, 1, 1};
    double tau2 = nu_tau2[0] * s2_tau2[0] / (nu_tau2[0] - 2);
    double alpha = 2.0;
    bmagwa::Prior prior(&data, &data_model, types_prior, 2, 4, inv_tau2_e,
        nu_sigma, s2_sigma, nu_tau2, s2_tau2, alpha, false);
    bmagwa::SamplerStats samplerstats;

    bmagwa::Model model(&data, &data_model, &prior, &samplerstats);

    double inv_tau2_alpha2s[2];
    inv_tau2_alpha2s[0] = 1.0 / (alpha * alpha * tau2);
    inv_tau2_alpha2s[1] = 1.0 / (alpha * alpha * tau2);

    model.add_term(1, 0, inv_tau2_alpha2s);
    //model.compute_log_likelihood();

    {
      // hand/matlab computed value (alpha is computed into the X matrix here
      // and does not appear in log(tau2) and log(4.0 * ...))
      double Syy = 0.58 - 0.51670588235294123702;
      double log_likelihood = -0.5 * log(tau2)
                                -0.5 * log(4.0 * 10.0 + 5.0 / tau2)
                                -0.5 * (nu_sigma + n) * log(nu_sigma * s2_sigma
                                                            + Syy);
      EXPECT_NEAR(log_likelihood, model.log_likelihood(), threshold);
    }

    model.add_term(5, 1, inv_tau2_alpha2s);
    {
      // hand/matlab computed value
      double Syy = 0.58 - 0.56537856592135671274;
      double log_likelihood = -0.5 * 3.0 * log(tau2)
                              -0.5 * log(810.62499999999954525265)
                              -0.5 * (nu_sigma + n) * log(nu_sigma * s2_sigma
                                                          + Syy);

      EXPECT_NEAR(log_likelihood, model.log_likelihood(), threshold);
    }
  }

  TEST(ModelTests, ExhaustiveComputation) {
    double threshold = 1e-10;

    size_t n = 5;

    bmagwa::Data data(n, 10, 0, "testdata/plinktest.fam",
        "testdata/plinktest.bed", false, "testdata/test.e");
    std::vector<bmagwa::DataModel::ef_t> types;
    types.push_back(bmagwa::DataModel::A);
    types.push_back(bmagwa::DataModel::AH);

    bmagwa::DataModel data_model(&data, types);
    double types_prior[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
    bmagwa::Vector inv_tau2_e(1); inv_tau2_e = 0;
    double nu_sigma = 4, s2_sigma = 1;
    double nu_tau2[] = {4, 4, 4, 4};
    double s2_tau2[] = {1, 1, 1, 1};
    double alpha = 2.0;
    bmagwa::Prior prior(&data, &data_model, types_prior, 2, 4, inv_tau2_e,
        nu_sigma, s2_sigma, nu_tau2, s2_tau2, alpha, false);
    bmagwa::SamplerStats samplerstats;

    bmagwa::Model model(&data, &data_model, &prior, &samplerstats);

    double tau2 = nu_tau2[0] * s2_tau2[0] / (nu_tau2[0] - 2);
    double inv_tau2_alpha2s[2];
    inv_tau2_alpha2s[0] = 1.0 / (alpha * alpha * tau2);
    inv_tau2_alpha2s[1] = 1.0 / (alpha * alpha * tau2);

    model.add_term(0, 0, inv_tau2_alpha2s);
    model.add_term(1, 1, inv_tau2_alpha2s);
    model.add_term(2, 0, inv_tau2_alpha2s);
    model.add_term(3, 0, inv_tau2_alpha2s);
    model.add_term(4, 1, inv_tau2_alpha2s);
    assert(std::isfinite(model.log_likelihood()));
    char addtypes[] = {0, 1, 0, 0, 1};

    {
      bmagwa::ExhModel exh_model(&data, &data_model, &prior, &samplerstats);
      exh_model.update_to_model(model, 0);
      exh_model.update_on_add();
      exh_model.update_on_add();
      exh_model.update_on_add();
      exh_model.update_on_add();
      exh_model.update_on_add();
      ASSERT_NEAR(model.log_likelihood(), exh_model.log_likelihood(), threshold);

      exh_model.update_on_moveleft();

      bmagwa::Model model2(&data, &data_model, &prior, &samplerstats);
      model2.add_term(0, 0, inv_tau2_alpha2s);
      model2.add_term(1, 1, inv_tau2_alpha2s);
      model2.add_term(2, 0, inv_tau2_alpha2s);
      model2.add_term(4, 1, inv_tau2_alpha2s);
      ASSERT_NEAR(model2.log_likelihood(), exh_model.log_likelihood(), threshold);
    }

    bmagwa::ExhModel exh_model(&data, &data_model, &prior, &samplerstats);
    exh_model.update_to_model(model, 0);
    double* model_probs = new double[(size_t)std::pow(2, 5)];

    double log_max_model;
    bmagwa::compute_exhaustive_modelset(5, &exh_model, model_probs, log_max_model);

    for (size_t i = 0; i < std::pow(2, 5); ++i){
      bmagwa::Model model2(&data, &data_model, &prior, &samplerstats);
      double log_model_prior = 0;
      int Ns[] = {0, 0, 0, 0, 0};
      int nloci = 0;
      for (size_t j = 0; j < 5; ++j){
        if (((i >> j) & 0x1) == 0x1){
          model2.add_term(j, addtypes[j], inv_tau2_alpha2s);
          log_model_prior += prior.compute_log_change_on_add(Ns, nloci, data_model.types[addtypes[j]]);
          ++Ns[data_model.types[addtypes[j]]];
          ++nloci;
        }
      }
      ASSERT_NEAR(log_model_prior + model2.log_likelihood(),
                  model_probs[i], threshold);
    }

    delete[] model_probs;
  }

  TEST(ModelTests, ExhaustiveComputationWithConstLoci) {
    double threshold = 1e-10;

    size_t n = 5;

    bmagwa::Data data(n, 10, 0, "testdata/plinktest.fam",
        "testdata/plinktest.bed", false, "testdata/test.e");
    std::vector<bmagwa::DataModel::ef_t> types;
    types.push_back(bmagwa::DataModel::A);
    types.push_back(bmagwa::DataModel::AH);

    bmagwa::DataModel data_model(&data, types);
    double types_prior[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
    bmagwa::Vector inv_tau2_e(1); inv_tau2_e = 0;
    double nu_sigma = 4, s2_sigma = 1;
    double nu_tau2[] = {4, 4, 4, 4};
    double s2_tau2[] = {1, 1, 1, 1};
    double alpha = 2.0;
    bmagwa::Prior prior(&data, &data_model, types_prior, 2, 4, inv_tau2_e,
        nu_sigma, s2_sigma, nu_tau2, s2_tau2, alpha, false);
    bmagwa::SamplerStats samplerstats;

    bmagwa::Model model(&data, &data_model, &prior, &samplerstats);

    double tau2 = nu_tau2[0] * s2_tau2[0] / (nu_tau2[0] - 2);
    double inv_tau2_alpha2s[2];
    inv_tau2_alpha2s[0] = 1.0 / (alpha * alpha * tau2);
    inv_tau2_alpha2s[1] = 1.0 / (alpha * alpha * tau2);

    model.add_term(0, 1, inv_tau2_alpha2s);
    model.add_term(1, 0, inv_tau2_alpha2s);
    model.add_term(2, 0, inv_tau2_alpha2s);
    model.add_term(3, 1, inv_tau2_alpha2s);
    model.add_term(4, 0, inv_tau2_alpha2s);
    model.add_term(5, 0, inv_tau2_alpha2s);
    model.add_term(6, 0, inv_tau2_alpha2s);
    model.add_term(7, 0, inv_tau2_alpha2s);
    assert(std::isfinite(model.log_likelihood()));
    char addtypes[] = {1, 0, 0, 1, 0, 0, 0, 0};


    bmagwa::ExhModel exh_model(&data, &data_model, &prior, &samplerstats);
    exh_model.update_to_model(model, 2);
    double* model_probs = new double[(size_t)std::pow(2, 6)];

    double max_log_model;
    bmagwa::compute_exhaustive_modelset(6, &exh_model, model_probs, max_log_model);

    for (size_t i = 0; i < std::pow(2, 6); ++i){
      bmagwa::Model model2(&data, &data_model, &prior, &samplerstats);
      model2.add_term(0, 1, inv_tau2_alpha2s);
      model2.add_term(1, 0, inv_tau2_alpha2s);
      double log_model_prior = 0;
      int Ns[] = {1, 0, 0, 0, 1};
      int nloci = 2;
      for (size_t j = 0; j < 6; ++j){
        if (((i >> j) & 0x1) == 0x1){
          model2.add_term(j+2, addtypes[j+2], inv_tau2_alpha2s);
          log_model_prior += prior.compute_log_change_on_add(Ns, nloci, data_model.types[addtypes[j+2]]);
          ++Ns[data_model.types[addtypes[j+2]]];
          ++nloci;
        }
      }
      ASSERT_NEAR(log_model_prior + model2.log_likelihood(),
                  model_probs[i], threshold);
    }

    delete[] model_probs;
  }

  TEST(ModelTests, SimpleAddRemove) {
    double threshold = 0.000001;

    bmagwa::Data data(5, 10, 0, "testdata/plinktest.fam",
                      "testdata/plinktest.bed", false, "testdata/test.e");
    std::vector<bmagwa::DataModel::ef_t> types;
    types.push_back(bmagwa::DataModel::A);
    types.push_back(bmagwa::DataModel::AH);

    bmagwa::DataModel data_model(&data, types);
    double nu_tau2[] = {4, 4, 4, 4};
    double s2_tau2[] = {1, 1, 1, 1};
    double types_prior[2] = { 1.0, 1.0 };
    double alpha = 2.0;
    bmagwa::Vector inv_tau2_e(1); inv_tau2_e = 0;
    bmagwa::Prior prior(&data, &data_model, types_prior, 2, 4, inv_tau2_e,
    		            4, 1, nu_tau2, s2_tau2, alpha, false);
    bmagwa::SamplerStats samplerstats;

    bmagwa::Model model(&data, &data_model, &prior, &samplerstats);

    double tau2 = nu_tau2[0] * s2_tau2[0] / (nu_tau2[0] - 2);
    double inv_tau2_alpha2s[2];
    inv_tau2_alpha2s[0] = 1.0 / (alpha * alpha * tau2);
    inv_tau2_alpha2s[1] = 1.0 / (alpha * alpha * tau2);

    ASSERT_EQ((size_t)0, model.size());

    {
      double likelihood0_1 = model.log_likelihood();
      model.add_term(1, 0, inv_tau2_alpha2s);
      ASSERT_EQ((size_t)1, model.size());
      double likelihood1_1 = model.log_likelihood();
      model.compute_log_likelihood();
      double likelihood1_2 = model.log_likelihood();
      model.remove_term(0);
      ASSERT_EQ((size_t)0, model.size());
      double likelihood0_2 = model.log_likelihood();

      EXPECT_NEAR(likelihood0_1, likelihood0_2, threshold);
      EXPECT_NEAR(likelihood1_1, likelihood1_2, threshold);

      EXPECT_TRUE(fabs(likelihood0_1 - likelihood1_1) > threshold);
    }

    {
      double likelihood0_1 = model.log_likelihood();
      model.add_term(1, 0, inv_tau2_alpha2s);
      ASSERT_EQ((size_t)1, model.size());
      model.add_term(7, 1, inv_tau2_alpha2s);
      ASSERT_EQ((size_t)2, model.size());
      double likelihood1_1 = model.log_likelihood();
      model.compute_log_likelihood();
      double likelihood1_2 = model.log_likelihood();
      model.remove_term(0);
      ASSERT_EQ((size_t)1, model.size());
      model.remove_term(0);
      ASSERT_EQ((size_t)0, model.size());
      double likelihood0_2 = model.log_likelihood();

      EXPECT_NEAR(likelihood0_1, likelihood0_2, threshold);
      EXPECT_NEAR(likelihood1_1, likelihood1_2, threshold);

      EXPECT_TRUE(fabs(likelihood0_1 - likelihood1_1) > threshold);
    }

  }

  TEST(ModelTests, Sigma2BetaSampling) {
    double threshold = 0.0001;

    size_t n = 5, nsamples = 10000000;

    bmagwa::Data data(n, 10, 0, "testdata/plinktest.fam",
        "testdata/plinktest.bed", false, "testdata/test.e");
    std::vector<bmagwa::DataModel::ef_t> types;
    types.push_back(bmagwa::DataModel::A);
    types.push_back(bmagwa::DataModel::AH);

    bmagwa::DataModel data_model(&data, types);
    double nu_tau2[] = {4, 4, 4, 4};
    double s2_tau2[] = {1, 1, 1, 1};
    double types_prior[1] = { 1.0 };
    bmagwa::Vector inv_tau2_e(1); inv_tau2_e = 0;
    double nu_sigma = 0, s2_sigma = 1;
    double alpha = 1.0;
    bmagwa::Prior prior(&data, &data_model, types_prior, 2, 4, inv_tau2_e,
        nu_sigma, s2_sigma, nu_tau2, s2_tau2, alpha, false);
    bmagwa::Rand rand(1234, n + nu_sigma);
    bmagwa::SamplerStats samplerstats;

    bmagwa::Model model(&data, &data_model, &prior, &samplerstats);

    double tau2 = nu_tau2[0] * s2_tau2[0] / (nu_tau2[0] - 2);
    double inv_tau2_alpha2s[2];
    inv_tau2_alpha2s[0] = 1.0 / (alpha * alpha * tau2);
    inv_tau2_alpha2s[1] = 1.0 / (alpha * alpha * tau2);

    bmagwa::Vector s2(nsamples); s2 = 0;
    bmagwa::Vector b1(nsamples); b1 = 0;
    bmagwa::Vector b2(nsamples); b2 = 0;

    model.compute_log_likelihood();
    //model.add_term(1, 0);
    {
      double Syy = 0.0680;
      // true values
      double mu_s2 = (nu_sigma * s2_sigma + Syy) / (nu_sigma + n - 2);
      double var_s2 = mu_s2 * mu_s2 * 2.0 / (nu_sigma + n - 4.0);
      double mu_b1 = 0.32;
      double var_b1 = (nu_sigma * s2_sigma + Syy) / (nu_sigma + n) / n * (nu_sigma + n) / (nu_sigma + n - 2);

      // sampling
      for (size_t i = 0; i != nsamples; ++i){
        model.sample_beta_sigma2(rand);
        s2(i) = model.get_sigma2_value();
        b1(i) = model.get_beta()(0);
      }

      EXPECT_NEAR(mu_s2, s2.mean(), threshold);
      EXPECT_NEAR(var_s2, s2.var(), threshold);
      EXPECT_NEAR(mu_b1, b1.mean(), threshold);
      EXPECT_NEAR(var_b1, b1.var(), threshold);
    }


    model.add_term(5, 1, inv_tau2_alpha2s);
    {
      double Syy = 0.02988235294117647101;
      // true values (assume nu_sigma = 0)
      double mu_s2 = (nu_sigma * s2_sigma + Syy) / (nu_sigma + n - 2);
      double var_s2 = mu_s2 * mu_s2 * 2.0 / (nu_sigma + n - 4.0);
      double mu_b1 = -0.11999999999999999556;
      double mu_b2 = -0.03529411764705892018;
      double var_b1 = 0.4000 * Syy / (n - 2.0); // is this correct?!
      double var_b2 = 0.58823529411764718944 * Syy / (n - 2.0); // is this correct?!

      // sampling
      for (size_t i = 0; i != nsamples; ++i){
        model.sample_beta_sigma2(rand);
        s2(i) = model.get_sigma2_value();
        b1(i) = model.get_beta()(1);
        b2(i) = model.get_beta()(2);
      }

      EXPECT_NEAR(mu_s2, s2.mean(), threshold);
      EXPECT_NEAR(var_s2, s2.var(), threshold);
      EXPECT_NEAR(mu_b1, b1.mean(), threshold);
      EXPECT_NEAR(mu_b2, b2.mean(), threshold);
      EXPECT_NEAR(var_b1, b1.var(), threshold);
      EXPECT_NEAR(var_b2, b2.var(), threshold);
    }

  }
}

#endif /* MODEL_TESTS_HPP_ */
