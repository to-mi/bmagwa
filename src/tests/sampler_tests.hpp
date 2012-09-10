/* BMAGWA software v2.0
 *
 * sampler_tests.hpp
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

#ifndef SAMPLER_TESTS_HPP_
#define SAMPLER_TESTS_HPP_

#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include <boost/math/distributions/non_central_f.hpp>
#include "../options.hpp"
#include "../data.hpp"
#include "../sampler.hpp"
#include "../data_model.hpp"
#include "../prior.hpp"
#include "../model.hpp"
#include "../precomputed_snp_covariances.hpp"
#include "../utils.hpp"
#include "../Test/Cuba-3.0/cuba.h"

namespace {

  class SamplerTestsF : public ::testing::Test {
    protected:
      SamplerTestsF()
      : opt(NULL), data(NULL), data_model(NULL),
        pre_xxcov(NULL), sampler(NULL),
        rao_threshold(0.00000000000001), // these thresholds cause the tests to fail
        freq_threshold(0.00000000000001)
      {}

      ~SamplerTestsF()
      {
        delete sampler;
        delete pre_xxcov;
        delete data_model;
        delete data;
        delete opt;
      }

      virtual void SetUp()
      {}

      void SetUpSampler(const std::string options_file)
      {
        opt = new bmagwa::Options(options_file);
        data = new bmagwa::Data(opt->n, opt->m_g, opt->m_e,
                                opt->file_fam, opt->file_g,
                                opt->recode_g_to_minor_allele_count,
                                opt->file_e, opt->file_y);
        data_model = new bmagwa::DataModel(data, opt->types);
        pre_xxcov = new bmagwa::PrecomputedSNPCovariances(data_model);
        sampler = new bmagwa::Sampler(*opt, data, pre_xxcov);

        sampler->initialize_p_proposal_flat();
      }

      void compare_to_rao(const std::string rao_file,
                            const bmagwa::VectorView pips) const
      {
        std::ifstream rao(rao_file.c_str(), std::ios::in | std::ios::binary);
        for (size_t i = 0; i < opt->m_g; ++i){
          double rao_val;
          rao.read(reinterpret_cast<char*>(&rao_val), sizeof(rao_val));
          EXPECT_NEAR(rao_val, pips(i), rao_threshold);
        }
        rao.close();
      }

      void compare_to_freq(const std::string freq_file_basename,
                              const bmagwa::VectorView pips,
                              const size_t burnin, const size_t niter) const
      {
        // read model size and loci data and compute the frequencies
        std::string ms_filename = freq_file_basename + std::string("_modelsize.dat");
        std::string loci_filename = freq_file_basename + std::string("_loci.dat");

        std::ifstream ms_file(ms_filename.c_str(), std::ios::in | std::ios::binary);
        std::ifstream loci_file(loci_filename.c_str(), std::ios::in | std::ios::binary);

        bmagwa::Vector freqs(pips.length());
        freqs = 0;

        uint32_t ms = 0;
        uint32_t loci[7];
        size_t nsamples = 0;
        for (size_t i = 0; i < niter; ++i){
          ms_file.read(reinterpret_cast<char*>(&ms), sizeof(ms));
          loci_file.read(reinterpret_cast<char*>(&loci[0]), ms * sizeof(loci[0]));
          if (i < burnin){
            for (unsigned int j = 0; j < ms; ++j){
              ++freqs(loci[j]);
            }
            ++nsamples;
          }
        }
        freqs /= nsamples;

        ms_file.close();
        loci_file.close();

        // check
        for (size_t i = 0; i < pips.length(); ++i){
             EXPECT_NEAR(freqs(i), pips(i), freq_threshold);
        }
      }

      // virtual void TearDown() {}

      bmagwa::SamplerStats samplerstats;
      bmagwa::Options* opt;
      bmagwa::Data* data;
      bmagwa::DataModel* data_model;
      bmagwa::PrecomputedSNPCovariances* pre_xxcov;
      bmagwa::Sampler* sampler;
      double rao_threshold;
      double freq_threshold;
  };

  struct idata {
    bmagwa::Model* model;
    bmagwa::Prior* prior;
    double nu_tau2;
    double s2_tau2;
    double mu_alpha;
    double empty_model_log_likelihood;
    int alphaidx;
  };
  static int IntegrandShared(const int *ndim, const double xx[],
      const int *ncomp, double ff[], void *userdata)
  {
    idata* idt = (idata*)userdata;
    double t2scale = 1; // larger than in IntegrandIndividual (i.e. covers more mass and is slower for the same accuracy...)
    double alphascale = 8;

    // model size is ok as all effect are A!
    double lt = 0.0;

    // assumes mu_alpha = 1 (for having equal area on both sides of it)
    double alpha = alphascale * xx[idt->alphaidx] - 0.5 * alphascale + idt->mu_alpha;
    if (alpha != 0){
      idt->prior->set_alpha(alpha, idt->model); // raises exception if alpha = 0
    }
    double tmp = alpha - idt->mu_alpha;
    lt += -0.5 * (tmp * tmp);

    // prior is zero at zero
    if (xx[0] <= 0.0){
      ff[0] = 0.0;
      return 0;
    }

    // skip constant term
    double xval = xx[0] * t2scale;
    for (unsigned int i = 0; i < idt->model->size(); ++i){
      idt->model->set_inv_tau2_alpha2_value(i+1, 1.0 / (xval * alpha * alpha));
    }
    lt += bmagwa::Utils::sinvchi2_lpdf(xval, idt->nu_tau2, idt->s2_tau2);

    double m = 0.0;
    if (alpha != 0.0){ // zero alpha means empty model
      idt->model->compute_log_likelihood();
      m = idt->model->log_likelihood() - idt->empty_model_log_likelihood;
    }

    // 0.5 from scaling tau2 val to [0,0.5] and 8.0 from scaling alpha
    ff[0] = exp(m + lt) * alphascale * t2scale;

    return 0;
  }
  static int IntegrandIndividual(const int *ndim, const double xx[],
                       const int *ncomp, double ff[], void *userdata)
  {
    idata* idt = (idata*)userdata;

    // model size is ok as all effect are A!
    double lt = 0.0;

    // assumes mu_alpha = 1 (for having equal area on both sides of it)
    double alpha = 8.0 * xx[idt->alphaidx] - 3.0;
    idt->prior->set_alpha(alpha, idt->model);
    double tmp = alpha - idt->mu_alpha;
    lt += -0.5 * (tmp * tmp);

    for (unsigned int i = 0; i < idt->model->size(); ++i){
      // prior is zero at zero
      if (xx[i] <= 0.0){
        ff[0] = 0.0;
        return 0;
      }

      // skip constant term
      double xval = xx[i] * 0.5;
      idt->model->set_inv_tau2_alpha2_value(i+1, 1.0 / (xval * alpha * alpha));
      lt += bmagwa::Utils::sinvchi2_lpdf(xval, idt->nu_tau2, idt->s2_tau2);

      //lt += log(boost::math::pdf(*(idt->ncf), xx[i] / idt->s2_tau2)) - log(idt->s2_tau2);
    }

    double m = 0.0;
    if (alpha != 0.0){ // zero alpha means empty model
      idt->model->compute_log_likelihood();
      m = idt->model->log_likelihood() - idt->empty_model_log_likelihood;
    }

    // pow(0.5, ms) from scaling tau2 vals to [0,0.5] and 8.0 from scaling alpha
    ff[0] = exp(m + lt) * 8.0 * std::pow(0.5, idt->model->size());

    return 0;
  }

  TEST_F(SamplerTestsF, OneSNPNoTau2) {
      SetUpSampler("testdata/onesnp_notau2.ini");

      ASSERT_FALSE(opt->use_individual_tau2);
      ASSERT_TRUE(opt->m_g == 1);


      // sample (don't sample missing or tau2; although tau2 is sampled once at
      //         the beginning)
      sampler->sample();


      // compute comparison value
      bmagwa::Prior* prior_ = sampler->get_prior();
      bmagwa::DataModel* data_model_ = sampler->get_data_model();

      bmagwa::Vector posteriors(2);
      posteriors = 0;

      double tmp = 1.0;
      for(int i = 0; i < 2; ++i){
        bmagwa::Model model(data, data_model_, prior_, &samplerstats);
        if (i == 1) {
          model.add_term(0, 0, &tmp);
          model.set_inv_tau2_alpha2_value(1, prior_->get_inv_tau2_alpha2_value(bmagwa::DataModel::A));
          model.compute_log_likelihood();
        }

        posteriors(i) = model.log_likelihood()+prior_->compute_log_model(model.get_Ns());
      }

      bmagwa::Vector pips(1);
      pips(0) = exp(posteriors(1)-posteriors(0));
      pips(0) = pips(0) / (pips(0) + 1.0);

      // compare to sampled values
      compare_to_rao("/tmp/chain_onesnp_notau2_rao.dat", pips);
      size_t burnin = opt->n_rao_burnin * opt->n_rao / opt->thin;
      compare_to_freq("/tmp/chain_onesnp_notau2",
                      pips, burnin, opt->do_n_iter / opt->thin);
    }


  TEST_F(SamplerTestsF, SmallModelSpaceWithOneSNP) {
    SetUpSampler("testdata/onesnp_tau2.ini");

    ASSERT_FALSE(opt->use_individual_tau2);
    ASSERT_TRUE(opt->m_g == 1);

    // sample (don't sample missing or tau2; although tau2 is sampled once at
    //         the beginning)
    sampler->sample();

    // compute posterior over whole model space
    double s2_tau2[] = {0,0,0,0};

    s2_tau2[0] = opt->s2_tau2[0];
    bmagwa::Prior* prior_ = sampler->get_prior();

    double empty_model_log_likelihood;
    {
      bmagwa::Model model(data, data_model, prior_, &samplerstats);
      model.compute_log_likelihood();
      empty_model_log_likelihood = model.log_likelihood();
    }

    bmagwa::Vector posteriors(2);
    posteriors = 0;

    double tmp = 1.0;
    for(int i = 0; i < 2; ++i){
      bmagwa::Model model(data, data_model, prior_, &samplerstats);
      if (i == 1) model.add_term(0, 0, &tmp);

      // integrate over tau2 and alpha
      const double EPSREL = 1e-10;
      const double EPSABS = 1e-28;
      const int MINEVAL = 100;
      const int MAXEVAL = 30000000;
      const int NNEW = 3000;
      const double FLATNESS = 25.;
      const int KEY = 0;

      int nregions, neval, fail;
      double integral[1], error[1], prob[1];
      idata idt;
      idt.model = &model;
      idt.prior = prior_;
      idt.mu_alpha = opt->mu_alpha;
      idt.nu_tau2 = opt->nu_tau2[0];
      idt.s2_tau2 = s2_tau2[0];
      idt.empty_model_log_likelihood = empty_model_log_likelihood;
      idt.alphaidx = 1;

      Cuhre(2, 1, IntegrandShared, &idt,
        EPSREL, EPSABS, 0,
        MINEVAL, MAXEVAL, KEY,
        &nregions, &neval, &fail, integral, error, prob);
      /*Suave(2, 1, IntegrandShared, &idt,
        EPSREL, EPSABS, 0, 67534,
        MINEVAL, MAXEVAL, NNEW, FLATNESS,
        &nregions, &neval, &fail, integral, error, prob);*/

      posteriors(i) = log(integral[0])
                                     +prior_->compute_log_model(model.get_Ns());
    }

    bmagwa::Vector pips(1);
    pips(0) = exp(posteriors(1)-posteriors(0));
    pips(0) = pips(0) / (pips(0) + 1.0);

    // compare to sampled values
    compare_to_rao("/tmp/chain_onesnp_rao.dat", pips);
    size_t burnin = opt->n_rao_burnin * opt->n_rao / opt->thin;
    compare_to_freq("/tmp/chain_onesnp",
                    pips, burnin, opt->do_n_iter / opt->thin);
  }

  TEST_F(SamplerTestsF, SmallModelSpaceWithSharedTau2) {
    SetUpSampler("testdata/small_modelspace_tau2.ini");

    ASSERT_FALSE(opt->use_individual_tau2);

    // sample (don't sample missing or tau2; although tau2 is sampled once at
    //         the beginning)
    sampler->sample();

    // compute posterior over whole model space
    double s2_tau2[] = {0,0,0,0};

    s2_tau2[0] = opt->s2_tau2[0];
    bmagwa::Prior* prior_ = sampler->get_prior();


    bmagwa::Vector likelihoods(std::pow(2, opt->m_g));
    bmagwa::Vector priors(std::pow(2, opt->m_g));
    likelihoods = 0;
    priors = 0;

    double empty_model_log_likelihood;
    {
      bmagwa::Model model(data, data_model, prior_, &samplerstats);
      model.compute_log_likelihood();
      empty_model_log_likelihood = model.log_likelihood();
    }

    double tmp = 1.0;
    for (size_t i = 0; i < std::pow(2, opt->m_g); ++i){
      bmagwa::Model model(data, data_model, prior_, &samplerstats);
      for (size_t j = 0; j < opt->m_g; ++j){
        if (((i >> j) & 0x1) == 0x1){
          model.add_term(j, 0, &tmp);
        }
      }

      idata idt;
      idt.model = &model;
      idt.prior = prior_;
      idt.mu_alpha = opt->mu_alpha;
      idt.nu_tau2 = opt->nu_tau2[0];
      idt.s2_tau2 = s2_tau2[0];
      idt.empty_model_log_likelihood = empty_model_log_likelihood;
      idt.alphaidx = 1;

      int nregions, neval, fail;
      double integral[1], error[1], prob[1];

      const double EPSREL = 1e-10;
      const double EPSABS = 1e-28;
      const int MINEVAL = 100;
      const int MAXEVAL = 3000000;
      const int NNEW = 3000;
      const double FLATNESS = 25.;
      const int KEY = 0;

      /*Cuhre(model.size()+1, 1, Integrand, &idt,
        EPSREL, EPSABS, 0,
        MINEVAL, MAXEVAL, KEY,
        &nregions, &neval, &fail, integral, error, prob);*/
      Suave(2, 1, IntegrandShared, &idt,
        EPSREL, EPSABS, 0, 67534,
        MINEVAL, MAXEVAL, NNEW, FLATNESS,
        &nregions, &neval, &fail, integral, error, prob);

      likelihoods(i) = log(integral[0]);

      priors(i) = prior_->compute_log_model(model.get_Ns());
    }

    double max_l = likelihoods.max();
    double max_p = priors.max();
    likelihoods -= max_l;
    priors -= max_p;

    // compute posterior inclusion probs etc.
    double totalprob = 0;
    bmagwa::Vector pips(opt->m_g);
    pips = 0;

    for (size_t i = 0; i < std::pow(2, opt->m_g); ++i){
      totalprob += exp(likelihoods(i) + priors(i));

      for (size_t j = 0; j < opt->m_g; ++j){
        if (((i >> j) & 0x1) == 0x1){
          pips(j) += exp(likelihoods(i) + priors(i));
        }
      }
    }

    pips /= totalprob;

    // compare to sampled values
    compare_to_rao("/tmp/chain_rao.dat", pips);
    size_t burnin = opt->n_rao_burnin * opt->n_rao / opt->thin;
    compare_to_freq("/tmp/chain",
                    pips, burnin, opt->do_n_iter / opt->thin);
  }

  TEST_F(SamplerTestsF, SmallModelSpaceWithIndividualTau2) {
      SetUpSampler("testdata/small_modelspace_tau2_individual.ini");

      ASSERT_TRUE(opt->use_individual_tau2);

      // sample (don't sample missing or tau2; although tau2 is sampled once at
      //         the beginning)
      sampler->sample();

      // compute posterior over whole model space
      double s2_tau2[] = {0,0,0,0};

      s2_tau2[0] = opt->s2_tau2[0];
      bmagwa::Prior* prior_ = sampler->get_prior();

      bmagwa::Vector likelihoods(std::pow(2, opt->m_g));
      bmagwa::Vector priors(std::pow(2, opt->m_g));
      likelihoods = 0;
      priors = 0;

      double empty_model_log_likelihood;
      {
        bmagwa::Model model(data, data_model, prior_, &samplerstats);
        model.compute_log_likelihood();
        empty_model_log_likelihood = model.log_likelihood();
      }

      double tmp = 1.0;
      for (size_t i = 0; i < std::pow(2, opt->m_g); ++i){
        bmagwa::Model model(data, data_model, prior_, &samplerstats);
        for (size_t j = 0; j < opt->m_g; ++j){
          if (((i >> j) & 0x1) == 0x1){
            model.add_term(j, 0, &tmp);
          }
        }

        idata idt;
        idt.model = &model;
        idt.prior = prior_;
        idt.mu_alpha = opt->mu_alpha;
        idt.nu_tau2 = opt->nu_tau2[0];
        idt.s2_tau2 = s2_tau2[0];
        idt.empty_model_log_likelihood = empty_model_log_likelihood;
        idt.alphaidx = model.size();

        int nregions, neval, fail;
        double integral[1], error[1], prob[1];

        const double EPSREL = 1e-10;
        const double EPSABS = 1e-28;
        const int MINEVAL = 200;
        const int MAXEVAL = 6000000;
        const int NNEW = 3000;
        const double FLATNESS = 25.;
        const int KEY = 0;

        /*Cuhre(model.size()+1, 1, Integrand, &idt,
          EPSREL, EPSABS, 0,
          MINEVAL, MAXEVAL, KEY,
          &nregions, &neval, &fail, integral, error, prob);*/
        Suave(model.size()+1, 1, IntegrandIndividual, &idt,
          EPSREL, EPSABS, 0, 67534,
          MINEVAL, MAXEVAL, NNEW, FLATNESS,
          &nregions, &neval, &fail, integral, error, prob);

        likelihoods(i) = log(integral[0]);

        priors(i) = prior_->compute_log_model(model.get_Ns());
      }

      double max_l = likelihoods.max();
      double max_p = priors.max();
      likelihoods -= max_l;
      priors -= max_p;

      // compute posterior inclusion probs etc.
      double totalprob = 0;
      bmagwa::Vector pips(opt->m_g);
      pips = 0;

      for (size_t i = 0; i < std::pow(2, opt->m_g); ++i){
        totalprob += exp(likelihoods(i) + priors(i));
        for (size_t j = 0; j < opt->m_g; ++j){
          if (((i >> j) & 0x1) == 0x1){
            pips(j) += exp(likelihoods(i) + priors(i));
          }
        }
      }

      pips /= totalprob;

      // check against rao results
      compare_to_rao("/tmp/chain_individual_rao.dat", pips);
      size_t burnin = opt->n_rao_burnin * opt->n_rao / opt->thin;
      compare_to_freq("/tmp/chain_individual",
                      pips, burnin, opt->do_n_iter / opt->thin);
    }


  TEST(SamplerTests, MissingSampling)
  {
    double threshold = 1e-10;

    // ensure that X' X and X' y are updated correctly

    // setup
    bmagwa::Options opt("testdata/plinktest.ini");
    bmagwa::Data data(opt.n, opt.m_g, opt.m_e,
                      opt.file_fam, opt.file_g,
                      opt.recode_g_to_minor_allele_count,
                      opt.file_e, opt.file_y);
    bmagwa::DataModel data_model_tmp(&data, opt.types);
    bmagwa::PrecomputedSNPCovariances pre_xxcov(&data_model_tmp);
    bmagwa::Sampler sampler(opt, &data, &pre_xxcov);
    bmagwa::Rand rng(1234, opt.n + opt.nu_sigma2);

    bmagwa::Prior* prior_ = sampler.get_prior();
    bmagwa::Model* cmodel = sampler.get_current_model();
    bmagwa::Model* nmodel = sampler.get_new_model();
    bmagwa::DataModel* data_model = sampler.get_data_model();

    // 0 and 3 contain missing values
    cmodel->compute_log_likelihood();
    double inv_val[2];
    inv_val[0] = prior_->get_inv_tau2_alpha2_value(bmagwa::DataModel::A);
    inv_val[1] = prior_->get_inv_tau2_alpha2_value(bmagwa::DataModel::H);
    cmodel->add_term(0, 0, inv_val);
    cmodel->add_term(1, 0, inv_val);
    cmodel->add_term(3, 1, inv_val);

    // do some repetitions of the test
    for (int r = 0; r < 1000; ++r){
      cmodel->sample_beta_sigma2(rng);
      *nmodel = *cmodel;

      sampler.sample_missing();

      // check that cmodel is correctly updated

      // compute xx, xy
      bmagwa::SymmMatrix xx(5);
      bmagwa::Vector xy(5);
      bmagwa::Matrix x(opt.n, 5);
      x = 1;

      data_model->get_genotypes_additive(0, x.column(1));
      data_model->get_genotypes_additive(1, x.column(2));
      data_model->get_genotypes_additive(3, x.column(3));
      data_model->get_genotypes_heterozygous(3, x.column(4));

      xx.set_to_innerproduct(x);
      xy.set_to_product(x, data.y(), true);

      // check xy and the upper triangle of xx
      for (int i = 0; i < 5; ++i){
        EXPECT_NEAR(cmodel->get_xy()(i), xy(i), threshold);
        for (int j = i; j < 5; ++j){
          EXPECT_NEAR(cmodel->get_xx()(i, j), xx(i, j), threshold);
        }
      }
    }
  }


  TEST(SamplerTests, ExhProposals) {
    double threshold = 1e-10;

    // assume situation with snp 1 and 2 are in the model
    // and we have proposed {+0,-2,+3}
    // snp 1 is not involved
    int movesize = 3;
    double log_prop_probs[8] = {0,0,0,0,0,0,0,0}; // 2**3
    double q_add[] = {0.3, 0.4, 0.5, 0.6};
    double q_rem[] = {0.7, 0.6, 0.5, 0.4};
    double z_add = 0.3 + 0.5 + 0.6; // "empty" model (only snp1 in it)
    double z_rem = 0.6;             // "empty" model
    size_t const_loci = 1;
    size_t m_g = 4;
    double* log_q_add_types = NULL; // assume no types
    unsigned char bit_to_normalized_order[3];

    // this is how bit_to_normalized_order is filled:
    //for i_move = 0 to (movesize-1)
    //  ind = move_inds_map[i_move]; // {+0,+2,+3}
    //  model_ind = new_model->model_inds[ind]; // this the full model (with removals added to the end)
    //  //model_inds in the full model are {1,0,3,2} (snp 1 is always there, snp 0 is added the, snp 3 then and finally snp 2 (as it was removal)
    //  dr_bit_to_normalized_order[model_ind - const_loci] = i_move;
    //  dr_q_add[i_move] = q_add_[ind];
    //  dr_q_rem[i_move] = q_rem_[ind];
    bit_to_normalized_order[0] = 1 - 1;
    bit_to_normalized_order[1] = 3 - 1;
    bit_to_normalized_order[2] = 2 - 1;
    double dr_q_add[] = {0.3, 0.5, 0.6};
    double dr_q_rem[] = {0.7, 0.5, 0.4};


    bmagwa::compute_proposal_probs_for_exh_modelset(false,
                                                    movesize,
                                                    bit_to_normalized_order,
                                                    dr_q_add,
                                                    dr_q_rem,
                                                    z_add,
                                                    z_rem,
                                                    const_loci,
                                                    m_g,
                                                    log_q_add_types,
                                                    log_prop_probs);

    // log_prop_probs come out in bit order

    // right most = snp 0, second right most = snp 3, left most = snp 2
    // 0 000: {+0,+2,+3}
    EXPECT_NEAR(log_prop_probs[0], +log(q_add[0]/(q_add[0] + q_add[2] + q_add[3]))+LOG_HALF
                                   +log(q_add[2]/(q_add[2] + q_add[3]))+LOG_HALF
                                   +log(q_add[3]/q_add[3])+LOG_HALF
                                 , threshold);
    // 1 001: {-0,+2,+3}
    EXPECT_NEAR(log_prop_probs[1], +log(q_rem[0]/(q_rem[0]+q_rem[1]))+LOG_HALF
                                   +log(q_add[2]/(q_add[2] + q_add[3]))+LOG_HALF
                                   +log(q_add[3]/q_add[3])+LOG_HALF
                                 , threshold);
    // 2 010: {+0,+2,-3}
    EXPECT_NEAR(log_prop_probs[2], +log(q_add[0]/(q_add[0] + q_add[2]))+LOG_HALF
                                   +log(q_add[2]/q_add[2])+LOG_HALF
                                   +log(q_rem[3]/(q_rem[1]+q_rem[3]))
                                 , threshold);
    // 3 011: {-3,+2,-0}
    EXPECT_NEAR(log_prop_probs[3], +log(q_rem[3]/(q_rem[0]+q_rem[1]+q_rem[3]))+LOG_HALF
                                   +log(q_add[2]/q_add[2])+LOG_HALF
                                   +log(q_rem[0]/(q_rem[1]+q_rem[0]))
                                 , threshold);
    // 4 100: {+0,-2,+3}
    EXPECT_NEAR(log_prop_probs[4], +log(q_add[0]/(q_add[0] + q_add[3]))+LOG_HALF
                                   +log(q_rem[2]/(q_rem[1]+q_rem[2]))+LOG_HALF
                                   +log(q_add[3]/q_add[3])+LOG_HALF
                                 , threshold);
    // 5 101: {-2,-0,+3}
    EXPECT_NEAR(log_prop_probs[5], +log(q_rem[2]/(q_rem[0]+q_rem[1]+q_rem[2]))+LOG_HALF
                                   +log(q_rem[0]/(q_rem[1]+q_rem[0]))+LOG_HALF
                                   +log(q_add[3]/q_add[3])+LOG_HALF
                                 , threshold);
    // 6 110: {+0,-3,-2}
    EXPECT_NEAR(log_prop_probs[6], +log(q_add[0]/q_add[0])+LOG_HALF
                                   +log(q_rem[3]/(q_rem[3]+q_rem[1]+q_rem[2]))
                                   +log(q_rem[2]/(q_rem[1]+q_rem[2]))
                                 , threshold);
    // 7 111: {-3,-2,-0}
    EXPECT_NEAR(log_prop_probs[7], +log(q_rem[3]/(q_rem[0]+q_rem[3]+q_rem[1]+q_rem[2]))
                                   +log(q_rem[2]/(q_rem[0]+q_rem[1]+q_rem[2]))
                                   +log(q_rem[0]/(q_rem[0]+q_rem[1]))
                                 , threshold);

  }


}

#endif /* SAMPLER_TESTS_HPP_ */
