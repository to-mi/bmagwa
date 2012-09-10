/* BMAGWA software v2.0
 *
 * raoblackwellizer_tests.hpp
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

#ifndef RAOBLACKWELLIZER_TESTS_HPP_
#define RAOBLACKWELLIZER_TESTS_HPP_

#include "gtest/gtest.h"
#include "../sampler.hpp"


namespace {

  TEST(RaoBlackwellizerTests, EmptyCurrentModel)
  {
    double threshold = 1e-10;

    // setup (note: shared tau2)
    bmagwa::Options opt("testdata/small_modelspace_tau2.ini");
    bmagwa::Data data(opt.n, opt.m_g, opt.m_e,
                      opt.file_fam, opt.file_g,
                      opt.recode_g_to_minor_allele_count,
                      opt.file_e, opt.file_y);
    bmagwa::DataModel data_model(&data, opt.types);
    bmagwa::PrecomputedSNPCovariances pre_xxcov(&data_model);
    bmagwa::Sampler sampler(opt, &data, &pre_xxcov);
    bmagwa::Rand rng(1234, opt.n + opt.nu_sigma2);

    bmagwa::RaoBlackwellizer rb(&data_model, &pre_xxcov);

    bmagwa::Prior* prior_ = sampler.get_prior();
    bmagwa::Model* cmodel = sampler.get_current_model();
    cmodel->compute_log_likelihood();
    cmodel->sample_beta_sigma2(rng);
    double pves[3];
    cmodel->compute_pve(pves);
    cmodel->compute_log_likelihood();

    double empty_model_log_likelihood = cmodel->log_likelihood();

    double p_r[7]; // m_g = 7
    double** p_r_types = NULL;

    bmagwa::Vector y_hat(opt.n);
    y_hat = 0; // residual is y

    // compute RB
    ASSERT_EQ(cmodel->size(), 0);
    rb.p_raoblackwell(&sampler, y_hat, p_r, p_r_types);
    ASSERT_EQ(cmodel->size(), 0);

    // compare to "hand-computed"
    double sigma2 = cmodel->get_sigma2_value();
    double inv_tau2_alpha2 = prior_->get_inv_tau2_alpha2_value(bmagwa::DataModel::A);
    int Ns[] = {0, 0, 0, 0, 0};
    for (int j = 0; j < 7; ++j){
      bmagwa::Vector x(opt.n);
      data_model.get_genotypes_additive(j, x);

      double xx = bmagwa::VectorView::dotproduct(x, x);
      double yx = bmagwa::VectorView::dotproduct(data_model.y(), x);
      double xsum = x.sum();
      double rsum = data_model.y().sum();

      double det = opt.n * (xx + inv_tau2_alpha2) - xsum * xsum;
      double rxxr = (rsum * rsum * (xx + inv_tau2_alpha2) - 2.0 * yx * rsum * xsum + yx * yx * opt.n) / det;

      double pip = exp(0.5 * (log(opt.n) - log(det)) + 0.5 * log(inv_tau2_alpha2)
                       -0.5 * (rsum * rsum / opt.n - rxxr) / sigma2
                       + prior_->compute_log_change_on_add(Ns, 0, bmagwa::DataModel::A));

      pip = pip / (1.0 + pip);

      EXPECT_NEAR(p_r[j], pip, threshold);
    }
  }

}


#endif /* RAOBLACKWELLIZER_TESTS_HPP_ */
