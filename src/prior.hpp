/* BMAGWA software v2.0
 *
 * prior.hpp
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


#ifndef PRIOR_HPP_
#define PRIOR_HPP_

#include <fstream>
#include "data.hpp"
#include "data_model.hpp"
#include "linalg.hpp"
#include "utils.hpp"

namespace bmagwa {

// forward declaration
class Model;

//! Implements the prior structure and sampling procedures.
class Prior
{
  public:
    const double n_plus_nu, minus_n_plus_nu_2, nus2_plus_yy;
    const double mu_alpha;
    const bool use_individual_tau2;

    Prior(const Data* data, const DataModel* data_model,
          const double* _types_prior,
          const double e_qg, const double var_qg,
          const Vector& _inv_tau2_e,
          const double _nu_sigma2, const double _s2_sigma2,
          const double* _nu_tau2, const double* _s2_tau2,
          const double _mu_alpha, const bool _use_individual_tau2)
    : n_plus_nu(data->n + _nu_sigma2), minus_n_plus_nu_2(-0.5 * n_plus_nu),
      nus2_plus_yy(_nu_sigma2 * _s2_sigma2 + data->yy()),
      mu_alpha(_mu_alpha),
      use_individual_tau2(_use_individual_tau2),
      n(data->n), m_g(data->m_g), m_e(data->m_e), n_types(data_model->n_types),
      log_n_types(log(n_types)), g_a(NAN), g_b(NAN), inv_tau2_e(_inv_tau2_e),
      nu_sigma2(_nu_sigma2), s2_sigma2(_s2_sigma2),
      nus2_sigma2(nu_sigma2 * s2_sigma2),
      alpha_(std::max(0.5, mu_alpha)),
      alpha2_(alpha_ * alpha_),
      types_prior_sum(0),
      pve_pow(0.5 * (0.35 - 1.0))
    {
      // check sigma2 params
      if (s2_sigma2 <= 0){
        throw std::runtime_error(
                   "Prior specification error: s2_sigma2 is not positive "
                   "(Probably option s2_sigma2 is not specified and "
                   "using R2mode_sigma2 leads to the invalid value).");
      }
      if (nu_sigma2 < 0){
        throw std::runtime_error(
                   "Prior specification error: nu_sigma2 is negative.");
      }

      // solve model prior parameters
      solve_betadist_params(e_qg, var_qg);

      // initialize types
      for (int i = 0; i < 5; ++i){
        allow_types[i] = data_model->allow_types((DataModel::ef_t)i);
        if (allow_types[i]){
          types_prior[i] = _types_prior[i];
          types_prior_sum += _types_prior[i];
        } else {
          types_prior[i] = NAN;
        }
      }
      for (int i = 0; i < 4; ++i){
        allow_terms[i] = data_model->allow_terms((DataModel::ef_t)i);
        if (allow_terms[i]){
          nu_tau2[i] = _nu_tau2[i];
          s2_tau2[i] = _s2_tau2[i];

          // check params
          if (s2_tau2[i] <= 0){
            std::string msg("Prior specification error: s2_tau2 for ");
            msg += data_model->ef_name[i];
            msg += " is not positive.";
            throw std::runtime_error(msg.c_str());
          }
          if (nu_tau2[i] <= 0){
            std::string msg("Prior specification error: nu_tau2 for ");
            msg += data_model->ef_name[i];
            msg += " is not positive.";
            throw std::runtime_error(msg.c_str());
          }

          nus2_tau2[i] = nu_tau2[i] * s2_tau2[i];
          // init to expectation (or smaller if nu_tau2 < 2.1)
          double tau2_init = nus2_tau2[i] / std::max(0.1, (nu_tau2[i] - 2));
          inv_tau2_alpha2[i] = 1.0 / (alpha2_ * tau2_init);
        } else {
          nu_tau2[i] = NAN;
          s2_tau2[i] = NAN;
          nus2_tau2[i] = NAN;
          inv_tau2_alpha2[i] = NAN;
        }
      }
      allow_terms[4] = false; // never allow AH (instead it's split to A and H)
    }

    ~Prior() {}

    const double& get_inv_tau2_alpha2_value(const DataModel::ef_t type) const
    {
      if (use_individual_tau2)
        throw std::runtime_error("Cannot call get_inv_tau2_alpha2_value if use_individual_tau2 is set.");

      return inv_tau2_alpha2[type];
    }

    const double& alpha() const
    {
      return alpha_;
    }
    const double& alpha2() const
    {
      return alpha2_;
    }
    void set_alpha(const double val, Model* model);

    double compute_log_change_on_add(const int* Ns_old, const int n_loci_old,
                                     const DataModel::ef_t type) const
    {
      return (log(g_a + n_loci_old) - log(g_b + m_g - n_loci_old - 1)) +
             (log(types_prior[type] + Ns_old[type])
              - log(types_prior_sum + n_loci_old));
    }

    double compute_log_change_on_rem(const int* Ns_old, const int n_loci_old,
                                     const DataModel::ef_t type) const
    {
      return (log(g_b + m_g - n_loci_old) - log(g_a + n_loci_old - 1)) +
             (log(types_prior_sum + n_loci_old - 1)
              - log(types_prior[type] + Ns_old[type] - 1));
    }

    double compute_log_change_on_swi(const int* Ns_old,
                                     const DataModel::ef_t type_add,
                                     const DataModel::ef_t type_rem) const
    {
      return (type_add == type_rem) ? 0 :
                            log(types_prior[type_add] + Ns_old[type_add])
                            - log(types_prior[type_rem] + Ns_old[type_rem] - 1);
    }

    double compute_log_model(const int *Ns) const
    {
      double lp = 0;
      int n_g = 0;
      for (size_t i = 0; i < 5; ++i){
        if (allow_types[i]){
          lp += Utils::gammaln(types_prior[i] + (double)Ns[i]);
          n_g += Ns[i];
        }
      }

      lp += Utils::gammaln(g_a + (double)n_g) + Utils::gammaln(g_b + (double)m_g - (double)n_g)
            - Utils::gammaln(types_prior_sum + (double)n_g);
      return lp;
    }

    void set_inv_tau2_e(VectorView& inv_tau2) const
    {
      // no alpha for covariates
      assert(inv_tau2.length() >= m_e);
      inv_tau2.block(0, m_e) = inv_tau2_e;
    }

    double sample_single_inv_tau2_from_prior_times_alpha2(const DataModel::ef_t type, Rand& rng) const
    {
      assert(type != DataModel::AH); // type should be term type
      if (use_individual_tau2)
        return 1.0 / (alpha2_ * rng.rand_sinvchi2(nu_tau2[type], s2_tau2[type]));
      else
        return inv_tau2_alpha2[type];
    }

    void sample_alpha_and_tau2(Model* model, const VectorView& y, Rand& rng)
    {
      /* note: must sample in this order!
       * This is because alpha_ is used in tau2 sampling to scale beta to eta,
       * but tau2 sampling should actually depend only on eta (i.e., must use
       * the same alpha in scaling which was in use when beta was sampled).
       */
      sample_tau2(model, rng);
      set_alpha(sample_alpha(model, y, rng), model);
    }

    double log_det_invQ_e() const
    {
      double log_sum_2 = 0.0;
      for (size_t i = 1; i < m_e; ++i){
        log_sum_2 += log(inv_tau2_e(i));
      }
      return log_sum_2;
    }

    void print(std::string filename)
    {
      std::ofstream* f_p = Utils::openfile(filename, false);

      *f_p << "a_omega = " << g_a << std::endl
           << "b_omega = " << g_b << std::endl
           << "nu_sigma2 = " << nu_sigma2 << std::endl
           << "s2_sigma2 = " << s2_sigma2 << std::endl
           << "mu_alpha = " << mu_alpha << std::endl;
      for (size_t t = 0; t < 5; ++t){
        *f_p << "type_" << DataModel::ef_name[t] << " = " << types_prior[t] << " (allowed? "
             << allow_types[t] << ")" << std::endl;
      }
      for (size_t t = 0; t < 4; ++t){
        *f_p << "term_" << DataModel::ef_name[t] << " allowed = " << allow_terms[t] << std::endl
             << "nu_tau2_" << DataModel::ef_name[t] << " = " << nu_tau2[t] << std::endl
             << "s2_tau2_" << DataModel::ef_name[t] << " = " << s2_tau2[t] << std::endl;
      }


      f_p->close();
      delete f_p;
    }

    double e_g() const
    {
      return g_a / (g_a + g_b) * m_g;
    }
  private:
    const size_t n, m_g, m_e;
    const int n_types;
    const double log_n_types;
    double g_a, g_b;
    Vector inv_tau2_e;
    double nu_sigma2, s2_sigma2, nus2_sigma2;
    // note: only four tau2s, since AH splits to A and H...
    double nu_tau2[4], s2_tau2[4], nus2_tau2[4];
    double alpha_, alpha2_;
    double inv_tau2_alpha2[4];
    double types_prior[5];
    double types_prior_sum;
    bool allow_types[5];
    bool allow_terms[5];
    double pve_pow;

    // disallow copy constructor and assignment
    Prior(const Prior& prior);
    Prior& operator=(const Prior& prior);

    double sample_alpha(Model* model, const VectorView& y, Rand& rng);
    void sample_tau2(Model* model, Rand& rng);

    void solve_betadist_params(const double e_q, const double var_q)
    {
      if (m_g > 1){
        double z = (var_q - e_q * (1 - e_q)) / ((m_g - 1) * e_q);

        g_a = (z - 1) / (1 - m_g * z / e_q);
        g_b = (m_g / e_q - 1) * g_a;
      } else {
        double z = e_q / (1 - e_q);

        g_b = z / var_q / pow(1 + z, 3) - 1 / (1 + z);
        g_a = z * g_b;
      }

      if (g_a <= 0 || g_b <= 0 || std::isnan(g_a) || std::isnan(g_b)){
        throw std::runtime_error("Invalid a or b.\n");
      }
    }
};

} // namespace bmagwa

#endif /* PRIOR_HPP_ */
