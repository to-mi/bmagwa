/* BMAGWA software v1.0
 *
 * prior.hpp
 *
 * http://www.lce.hut.fi/research/mm/bmagwa/
 * Copyright 2011 Tomi Peltola <tomi.peltola@aalto.fi>
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

class Prior
{
  public:
    double n_plus_nu, minus_n_plus_nu_2, nus2_plus_yy;

    Prior(const Data* data, const DataModel* data_model,
          const double* _types_prior,
          const double e_qg, const double var_qg,
          const Vector& _inv_tau2_e,
          const double _nu_sigma2, const double _s2_sigma2,
          const double* _nu_tau2, const double* _s2_tau2)
    : n_plus_nu(data->n + _nu_sigma2), minus_n_plus_nu_2(-0.5 * n_plus_nu),
      nus2_plus_yy(_nu_sigma2 * _s2_sigma2 + data->yy()),
      m_g(data->m_g), m_e(data->m_e), n_types(data_model->n_types),
      log_n_types(log(n_types)), g_a(NAN), g_b(NAN), inv_tau2_e(_inv_tau2_e),
      nu_sigma2(_nu_sigma2), s2_sigma2(_s2_sigma2),
      nus2_sigma2(nu_sigma2 * s2_sigma2),
      types_prior_sum(0)
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
          t2_tau2[i] = nus2_tau2[i] / std::max(0.1, (nu_tau2[i] - 2));
        } else {
          nu_tau2[i] = NAN;
          s2_tau2[i] = NAN;
          nus2_tau2[i] = NAN;
          t2_tau2[i] = NAN;
        }
      }
      allow_terms[4] = false; // never allow AH (instead it's split to A and H)
    }

    ~Prior() {}

    const double& get_tau2_value(const DataModel::ef_t type) const
    {
      return t2_tau2[type];
    }
    void set_tau2_value(const DataModel::ef_t type, const double val)
    {
      t2_tau2[type] = val;
    }

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

    void sample_tau2(const std::vector<DataModel::ef_t>& x_types,
                     const VectorView& beta, const double sigma2, Rand& rng)
    {
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
          double s2 = (nus2_tau2[t] + beta2sum[t] / sigma2) / nu;

          t2_tau2[t] = rng.rand_sinvchi2(nu, s2);
        }
      }
    }

    double do_invQ(SymmMatrixView& xx,
                   const std::vector<DataModel::ef_t>& x_types) const
    {
      assert(xx.length() >= m_e);

      double log_sum_2 = 0;

      // handle E
      // assumption: first term is constant
      xx(0, 0) += inv_tau2_e(0);
      for (size_t i = 1; i < m_e; ++i){
        xx(i, i) += inv_tau2_e(i);
        log_sum_2 += log(inv_tau2_e(i));
      }

      // G effects
      double inv_val[] = {1.0 / t2_tau2[0],
                          1.0 / t2_tau2[1],
                          1.0 / t2_tau2[2],
                          1.0 / t2_tau2[3]};
      for (size_t i = m_e; i < xx.length(); ++i){
        xx(i, i) += inv_val[x_types[i]];
        log_sum_2 += log(inv_val[x_types[i]]);
      }
      //log_sum_2 += log(inv_val) * (xx.length() - m_e);

      return 0.5 * log_sum_2;
    }

    double log_pdf_tau2(const DataModel::ef_t type)
    {
      return Utils::sinvchi2_lpdf(t2_tau2[type], nu_tau2[type], s2_tau2[type]);
    }

    void print(std::string filename)
    {
      std::ofstream* f_p = Utils::openfile(filename, false);

      *f_p << "a_omega = " << g_a << std::endl
           << "b_omega = " << g_b << std::endl
           << "nu_sigma2 = " << nu_sigma2 << std::endl
           << "s2_sigma2 = " << s2_sigma2 << std::endl;
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

  private:
    const size_t m_g, m_e;
    const int n_types;
    const double log_n_types;
    double g_a, g_b;
    Vector inv_tau2_e;
    double nu_sigma2, s2_sigma2, nus2_sigma2;
    // note: only four tau2s, since AH splits to A and H...
    double nu_tau2[4], s2_tau2[4], nus2_tau2[4], t2_tau2[4];
    double types_prior[5];
    double types_prior_sum;
    bool allow_types[5];
    bool allow_terms[5];

    // disallow copy constructor and assignment
    Prior(const Prior& prior);
    Prior& operator=(const Prior& prior);

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
