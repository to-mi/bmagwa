/* BMAGWA software v2.0
 *
 * prior.cpp
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


#include "prior.hpp"
#include "model.hpp"
#include "rand.hpp"
#include "linalg.hpp"

namespace bmagwa {

double Prior::sample_alpha(Model* model, const VectorView& y, Rand& rng)
{
  // then use full conditional of alpha to update it
  size_t nterms = (model->beta).length() - m_e;

  double a = alpha_;

  do { // disallow exact zero value
    if (nterms > 0){
      // propose switching sign of alpha and eta (sign of eta has no effect as long
      // as we actually sample beta; note that eta has symmetric zero mean
      // distribution); this is metropolis move with deterministic
      // proposal; no effect is mu_alpha is zero (might as well skip...)
      if ((a <= 0.0) || log(rng.rand_01()) <= -2.0 * mu_alpha * a){
        a = -a;
      }

      Vector xb(n);
      Vector eb(n);

      xb.set_to_product((model->x).columnblock(m_e, nterms),
          (model->beta).block(m_e, nterms), false);
      eb.set_to_product((model->x).columnblock(0, m_e),
          (model->beta).block(0, m_e), false);
      eb -= y; // note: this is E * b - y, hence there is minus below in mu formul.

      double alpha_sigma2 = a * model->sigma2;
      double var = 1.0 / (1.0 + VectorView::dotproduct(xb, xb) / (a * alpha_sigma2));
      double mu = var * (mu_alpha - VectorView::dotproduct(eb, xb) / alpha_sigma2);

      a = sqrt(var) * rng.rand_normal() + mu;
    } else {
      // sample from prior
      a = rng.rand_normal() + mu_alpha;
    }
  } while (a == 0 || a * a == 0);
  model->mu_beta_computed = false;

  return a;
}

void Prior::sample_tau2(Model* model, Rand& rng)
{
  const std::vector<DataModel::ef_t>& x_types = model->x_types;
  const VectorView& beta = model->beta;

  double alpha2_sigma2 = alpha2_ * model->sigma2;

  if (use_individual_tau2){
    for (size_t i = m_e; i < beta.length(); ++i){
      double nu = nu_tau2[x_types[i]] + 1.0;
      double s2 = (nus2_tau2[x_types[i]] + beta(i) * beta(i) / alpha2_sigma2) / nu;

      model->inv_tau2_alpha2(i) = 1.0 / (alpha2_ * rng.rand_sinvchi2(nu, s2));
      assert(model->inv_tau2_alpha2(i) > 0);
    }
  } else {
    size_t m = beta.length();

    double beta2sum[] = {0.0, 0.0, 0.0, 0.0};
    size_t m_types[] = {0, 0, 0, 0};

    for (size_t i = m_e; i < m; ++i){
      beta2sum[x_types[i]] += beta(i) * beta(i);
      ++m_types[x_types[i]];
    }

    for (size_t t = 0; t < 4; ++t){
      if (allow_terms[t]){
        double nu = nu_tau2[t] + m_types[t];
        double s2 = (nus2_tau2[t] + beta2sum[t] / alpha2_sigma2) / nu;
        assert(std::isfinite(s2));

        inv_tau2_alpha2[t] = 1.0 / (alpha2_ * rng.rand_sinvchi2(nu, s2));
        assert(inv_tau2_alpha2[t] > 0);
      }
    }

    for (size_t i = m_e; i < m; ++i){
      model->inv_tau2_alpha2(i) = inv_tau2_alpha2[x_types[i]];
    }
  }

  model->mu_beta_computed = false;
}

void Prior::set_alpha(const double val, Model* model)
{
  assert(val != 0 && val * val != 0);
  double old_alpha2 = alpha2_;
  alpha_ = val;
  alpha2_ = alpha_ * alpha_;

  // need to update inv_tau2_alpha2 values
  if (use_individual_tau2){
    for (size_t i = m_e; i < model->beta.length(); ++i){
      model->inv_tau2_alpha2(i) = model->inv_tau2_alpha2(i) * old_alpha2 / alpha2_;
    }
  } else {
    const std::vector<DataModel::ef_t>& x_types = model->x_types;

    for (size_t t = 0; t < 4; ++t){
      if (allow_terms[t]){
        inv_tau2_alpha2[t] = inv_tau2_alpha2[t] * old_alpha2 / alpha2_;
      }
    }

    for (size_t i = m_e; i < model->beta.length(); ++i){
      model->inv_tau2_alpha2(i) = inv_tau2_alpha2[x_types[i]];
    }
  }
}

} // namespace bmagwa
