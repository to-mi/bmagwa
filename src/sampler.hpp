/* BMAGWA software v1.0
 *
 * sampler.hpp
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


#ifndef SAMPLER_HPP_
#define SAMPLER_HPP_

#include <string>
#include <cmath>
#include <fstream>
#include "linalg.hpp"
#include "rand.hpp"
#include "data_model.hpp"
#include "prior.hpp"
#include "model.hpp"
#include "options.hpp"
#include "precomputed_snp_covariances.hpp"

namespace bmagwa {

// forward declaration
class Sampler;

class RaoBlackwellizer
{
  public:
    RaoBlackwellizer(const DataModel* data_model,
                     const PrecomputedSNPCovariances* pre_xxcov)
    : pre_xxcov(pre_xxcov),
      m_g(data_model->m_g),
      // 2 = constant term + one genetic term
      v(2 + data_model->allow_types(DataModel::AH)),
      r_cm(data_model->n),
      r_om(data_model->n),
      x(data_model->n, 2 + data_model->allow_types(DataModel::AH)),
      xxcov(2 + data_model->allow_types(DataModel::AH))
    {
      x = 1;
    }

    ~RaoBlackwellizer() {}

    void p_raoblackwell(const Sampler* sampler, const Vector& y_hat,
                        double *p_r, double **p_r_types);

  private:
    // xx contains pre-computed covariance matrices for all SNPs in all types
    // xx_types contains pointers to xx in a different order
    const PrecomputedSNPCovariances* pre_xxcov;
    /*double **xx;         // in order of data_model->types
    double *xx_types[5]; // in order of ef_t enum
    size_t xx_offset[5]; // in order of ef_t enum*/
    const size_t m_g;

    Vector v, r_cm, r_om;
    Matrix x;
    SymmMatrix xxcov;

    // disallow copy constructor and assignment
    RaoBlackwellizer(const RaoBlackwellizer& rb);
    RaoBlackwellizer& operator=(const RaoBlackwellizer& rb);
};

class Sampler
{
  public:
    // TODO: group variables into manageable classes or something
    Sampler(const Options& opts, const Data* data,
            const PrecomputedSNPCovariances* pre_xxcov)
    : n_iter(0),
      do_n_iter(opts.do_n_iter),
      n_accepted(0),
      n_rao(opts.n_rao),
      n_rao_burnin(opts.n_rao_burnin),
      n_sample_tau2_and_missing(opts.n_sample_tau2_and_missing),
      adaptation(opts.adaptation),
      save_beta(opts.save_beta),
      t_proposal(opts.t_proposal),
      thin(opts.thin),
      verbosity(opts.verbosity),
      basename(opts.basename),
      level(std::max(0, (int)round((log(10) - log(data->m_g)) / log(0.5)))),
      m_g(data->m_g),
      n_types(opts.types.size()),
      use_types(n_types > 1),
      seed(opts.seed),
      rng(seed, data->n + opts.nu_sigma2),
      p_proposal(NULL),
      p_proposal_types(NULL),
      p_proposal_n(1),
      p_rao(NULL),
      p_rao_types(NULL),
      p_rao_n(0),
      q_add_(NULL),
      q_p_add_(NULL),
      new_q_add_(NULL),
      new_q_p_add_(NULL),
      p_r_(NULL),
      p_r_types_(NULL),
      q_add_types_(NULL),
      q_p_add_types_(NULL)
    {
      data_model = new DataModel(data, opts.types);

      // solve prior variance scale params
      double s2_sigma2 = 0.0;
      if (opts.s2_sigma2 > 0){
        s2_sigma2 = opts.s2_sigma2;
      } else {
        s2_sigma2 = data->var_y() * (1 - opts.R2mode_sigma2)
                           * (opts.nu_sigma2 + 2) / opts.nu_sigma2;
      }

      double s2_tau2[] = {0.0, 0.0, 0.0, 0.0};
      for (size_t t = 0; t < 4; ++t){
        if (opts.s2_tau2[t] > 0){
          s2_tau2[t] = opts.s2_tau2[t];
        } else {
          s2_tau2[t] = opts.eh_tau2[t] * (opts.nu_tau2[t] - 2)
                         / (opts.nu_tau2[t]
                         * (data->var_x() + data->mean_x() * data->mean_x())
                         * (1 - opts.R2mode_sigma2));
        }
      }

      assert(opts.inv_tau2_e != NULL);
      prior = new Prior(data, data_model, opts.types_prior, opts.e_qg,
                        opts.var_qg, *(opts.inv_tau2_e), opts.nu_sigma2,
                        s2_sigma2, opts.nu_tau2, s2_tau2);

      current_model = new Model(data, data_model, prior);
      new_model = new Model(data, data_model, prior);

      // allocate proposal distribution (initialization is done outside of ctor)
      p_proposal = new double[m_g];
      for (size_t i = 0; i < m_g; ++i) p_proposal[i] = 0.0;
      if (data_model->n_types > 1){
        p_proposal_types = new double*[m_g];
        for (size_t i = 0; i < m_g; ++i){
          p_proposal_types[i] = new double[data_model->n_types];
          for (size_t j = 0; j < data_model->n_types; ++j)
            p_proposal_types[i][j] = 0;
        }
      }

      // allocate p_rao/_types
      p_rao = new double[m_g];
      for (size_t i = 0; i < m_g; ++i) p_rao[i] = 0;
      if (data_model->n_types > 1){
        p_rao_types = new double*[m_g];
        for (size_t i = 0; i < m_g; ++i){
          p_rao_types[i] = new double[data_model->n_types];
          for (size_t j = 0; j < data_model->n_types; ++j)
            p_rao_types[i][j] = 0;
        }
      }

      raob = new RaoBlackwellizer(data_model, pre_xxcov);

      q_add_ = new double[m_g];
      q_p_add_ = new double[m_g];
      new_q_add_ = new double[m_g];
      new_q_p_add_ = new double[m_g];
      p_r_ = new double[m_g];

      if (use_types){
        p_r_types_ = new double*[m_g];
        q_add_types_ = new double*[m_g];
        q_p_add_types_ = new double*[m_g];
        for (size_t i = 0; i < m_g; ++i){
          p_r_types_[i] = new double[n_types];
          q_add_types_[i] = new double[n_types];
          q_p_add_types_[i] = new double[n_types];
        }
      }
    }

    ~Sampler()
    {
      delete[] q_add_;
      delete[] q_p_add_;
      delete[] new_q_add_;
      delete[] new_q_p_add_;
      delete[] p_r_;
      if (use_types){
        for (size_t i = 0; i < m_g; ++i){
          delete[] p_r_types_[i];
          delete[] q_add_types_[i];
          delete[] q_p_add_types_[i];
        }
        delete[] p_r_types_;
        delete[] q_add_types_;
        delete[] q_p_add_types_;
      }

      delete[] p_rao;
      if (p_rao_types != NULL){
        for (size_t i = 0; i < m_g; ++i)
          delete[] p_rao_types[i];
        delete[] p_rao_types;
      }
      delete raob;
      delete[] p_proposal;
      if (p_proposal_types != NULL){
        for (size_t i = 0; i < m_g; ++i){
          delete[] p_proposal_types[i];
        }
        delete[] p_proposal_types;
      }
      delete current_model;
      delete new_model;
      delete prior;
      delete data_model;
    }

    void sample();
    void initialize_p_proposal()
    {
      data_model->sample_missing(current_model->model_inds, rng);
      current_model->sample_beta_sigma2(rng);

      Vector y_hat(data_model->n);
      y_hat.set_to_product(current_model->x, current_model->beta, false);

      raob->p_raoblackwell(this, y_hat, p_proposal,p_proposal_types);
    }
    void save_p_proposal()
    {
      std::ofstream* f_p = Utils::openfile(basename + "_singlesnp_post.dat",
                                                                          true);
      f_p->write(reinterpret_cast<const char*>(p_proposal),
                 m_g * sizeof(p_proposal[0]));
      f_p->close();
      delete f_p;

      if (use_types){
        std::ofstream* f_pt = Utils::openfile(basename
                                              + "_singlesnp_post_types.dat",
                                              true);
        for (size_t i = 0; i < m_g; ++i){
          f_pt->write(reinterpret_cast<const char*>(p_proposal_types[i]),
                      n_types * sizeof(p_proposal_types[i][0]));
        }
        f_pt->close();
        delete f_pt;
      }
    }
    void print_prior()
    {
      std::string filename(
               basename.substr(0, basename.length() - 1) + "_prior.dat", false);

      prior->print(filename);
    }

    const double get_tau2_value(const DataModel::ef_t type) const {
      return prior->get_tau2_value(type);
    }

  private:
    size_t n_iter, do_n_iter, n_accepted, n_rao, n_rao_burnin;
    size_t n_sample_tau2_and_missing;
    bool adaptation, save_beta;
    double t_proposal;
    size_t thin, verbosity;
    std::string basename;

    int level;
    const size_t m_g;
    const size_t n_types;
    const bool use_types;
    const uint32_t seed;

    Rand rng;
    //Data* data;
    Prior* prior;
    DataModel* data_model;
    Model* current_model;
    Model* new_model;

    double *p_proposal;
    double **p_proposal_types;
    size_t p_proposal_n;

    RaoBlackwellizer *raob;
    double *p_rao;
    double **p_rao_types;
    int p_rao_n;

    //
    double p_moves_[4], p_moves_cumsum_[4];
    double new_p_moves_[4], new_p_moves_cumsum_[4];

    double* q_add_;
    double* q_p_add_;
    double* new_q_add_;
    double* new_q_p_add_;
    double* p_r_;

    double** p_r_types_;
    double** q_add_types_;
    double** q_p_add_types_;
    double q_p_add_types_tmp_[5];


    // methods
    void sample_missing();
    void compute_p_moves(const size_t n_loci,
                         double* p_moves, double* p_moves_cumsum) const;
    void compute_proposal_dist_t(const Model* model,
                                 const double *p_proposal,
                                 double* q_add, double* q_p_add) const;
    void remdate_proposal_dist_t(const size_t ind, double p_proposal,
                                 double* q_add, double* q_p_add) const;
    void adddate_proposal_dist_t(const size_t ind, double p_proposal,
                                 double* q_add, double* q_p_add) const;
    void switchdate_proposal_dist_t(const size_t ind_add,
                                    double p_proposal_add,
                                    const size_t ind_rem,
                                    double p_proposal_rem,
                                    double* q_add, double* q_p_add) const;
    void compute_proposal_types_dist_t(const double*const* p_proposal_types,
                                       double** q_add_types,
                                       double** q_p_add_types) const;
    size_t sample_discrete(const double* cumsum, const size_t m);

    char do_addition()
    {
      size_t ind = sample_discrete(q_p_add_, m_g);
      int t_ind = 0;
      if (use_types) t_ind = Utils::sample_discrete_naive(q_p_add_types_[ind],
                                                          n_types, rng);
      DataModel::ef_t type = data_model->types[t_ind];

      // new model
      // sample missing values
      data_model->sample_missing_single(ind, rng);

      // add term
      new_model->add_term(ind, t_ind);

      compute_p_moves(new_model->size(), new_p_moves_, new_p_moves_cumsum_);

      double r = prior->compute_log_change_on_add(current_model->Ns,
                                                  current_model->size(),
                                                  type)
                 + new_model->log_likelihood_ - current_model->log_likelihood_
                 - log(new_model->size()) + log(new_p_moves_[1])
                 + log(q_p_add_[m_g - 1]) - log(q_add_[ind]) - log(p_moves_[0]);
      if (use_types)
        r -= log(q_add_types_[ind][t_ind] / q_p_add_types_[ind][n_types-1]);

      if (r >= 0 || log(rng.rand_01()) <= r){
        ++n_accepted;
        *current_model = *new_model;

        memcpy(p_moves_, new_p_moves_, 4 * sizeof(double));
        memcpy(p_moves_cumsum_, new_p_moves_cumsum_, 4 * sizeof(double));

        adddate_proposal_dist_t(ind, p_proposal[ind], q_add_, q_p_add_);
        memcpy(new_q_add_, q_add_, m_g  * sizeof(double));
        memcpy(new_q_p_add_, q_p_add_, m_g * sizeof(double));

        return 1;
      } else {
        *new_model = *current_model;
        return 0;
      }
    }

    char do_removal()
    {
      // find model ind
      int model_ind = (int)round(
                              (current_model->size() - 1) * rng.rand_01());

      uint32_t ind = current_model->loci[model_ind];
      int t_ind = current_model->type_inds[model_ind];
      DataModel::ef_t type = data_model->types[t_ind];

      new_model->remove_term(model_ind);

      compute_p_moves(new_model->size(), new_p_moves_, new_p_moves_cumsum_);
      remdate_proposal_dist_t(ind, p_proposal[ind], new_q_add_, new_q_p_add_);

      double r = prior->compute_log_change_on_rem(current_model->Ns,
                                                  current_model->size(),
                                                  type)
                 + new_model->log_likelihood_ - current_model->log_likelihood_
                 - log(new_q_p_add_[m_g - 1]) + log(new_q_add_[ind])
                 + log(new_p_moves_[0])
                 + log(current_model->size()) - log(p_moves_[1]);
      if (use_types)
        r += log(q_add_types_[ind][t_ind] / q_p_add_types_[ind][n_types-1]);

      if (r >= 0 || log(rng.rand_01()) <= r){
        ++n_accepted;
        *current_model = *new_model;

        memcpy(p_moves_, new_p_moves_, 4 * sizeof(double));
        memcpy(p_moves_cumsum_, new_p_moves_cumsum_, 4 * sizeof(double));

        memcpy(q_add_, new_q_add_, m_g * sizeof(double));
        memcpy(q_p_add_, new_q_p_add_, m_g * sizeof(double));

        return 1;
      } else {
        *new_model = *current_model;

        memcpy(new_q_add_, q_add_, m_g * sizeof(double));
        memcpy(new_q_p_add_, q_p_add_, m_g * sizeof(double));

        return 0;
      }
    }

    char do_switch_of_snps()
    {
      size_t ind = sample_discrete(q_p_add_, m_g);
      int t_ind = 0;
      if (use_types) t_ind = Utils::sample_discrete_naive(q_p_add_types_[ind], n_types, rng);
      DataModel::ef_t type_add = data_model->types[t_ind];

      // sample missing values
      data_model->sample_missing_single(ind, rng);
      int model_ind = (int)round(
                              (current_model->size() - 1) * rng.rand_01());
      size_t ind_rem = current_model->loci[model_ind];
      int t_ind_rem = current_model->type_inds[model_ind];
      DataModel::ef_t type_rem = data_model->types[t_ind_rem];

      // remove and add
      new_model->remove_term(model_ind);
      new_model->add_term(ind, t_ind);

      compute_p_moves(new_model->size(), new_p_moves_, new_p_moves_cumsum_);
      switchdate_proposal_dist_t(ind, p_proposal[ind],
                                 ind_rem, p_proposal[ind_rem],
                                 new_q_add_, new_q_p_add_);

      double r = prior->compute_log_change_on_swi(current_model->Ns, type_add,
                                                  type_rem)
                 + new_model->log_likelihood_ - current_model->log_likelihood_
                 - log(new_q_p_add_[m_g - 1]) + log(new_q_add_[ind_rem])
                 - log(new_model->size()) + log(new_p_moves_[2])
                 + log(q_p_add_[m_g - 1]) - log(q_add_[ind])
                 + log(current_model->size()) - log(p_moves_[2]);
      if (use_types)
        r += log(q_add_types_[ind_rem][t_ind_rem] / q_p_add_types_[ind_rem][n_types-1])
             - log(q_add_types_[ind][t_ind] / q_p_add_types_[ind][n_types-1]);

      if (r >= 0 || log(rng.rand_01()) <= r){
        ++n_accepted;
        *current_model = *new_model;

        memcpy(p_moves_, new_p_moves_, 4 * sizeof(double));
        memcpy(p_moves_cumsum_, new_p_moves_cumsum_, 4 * sizeof(double));

        memcpy(q_add_, new_q_add_, m_g * sizeof(double));
        memcpy(q_p_add_, new_q_p_add_, m_g * sizeof(double));

        return 1;
      } else {
        *new_model = *current_model;

        memcpy(new_q_add_, q_add_, m_g * sizeof(double));
        memcpy(new_q_p_add_, q_p_add_, m_g * sizeof(double));

        return 0;
      }
    }

    char do_switch_of_types()
    {
      int model_ind = (int)round((current_model->size() - 1) * rng.rand_01());
      size_t ind = current_model->loci[model_ind];

      size_t t_ind_rem = current_model->type_inds[model_ind];

      // cumulative sum without the current type
      q_p_add_types_tmp_[0] = (t_ind_rem == 0) ? 0 : q_add_types_[ind][0];
      for (size_t i = 1; i < n_types; ++i){
        q_p_add_types_tmp_[i] = q_p_add_types_tmp_[i-1]
                               + ((t_ind_rem == i) ? 0 : q_add_types_[ind][i]);
      }

      int t_ind = Utils::sample_discrete_naive(q_p_add_types_tmp_, n_types, rng);
      DataModel::ef_t type_add = data_model->types[t_ind];
      DataModel::ef_t type_rem = data_model->types[t_ind_rem];

      // remove and add
      new_model->remove_term(model_ind);
      new_model->add_term(ind, t_ind);

      // Delta posterior - proposal prob + reverse prob
      double r = prior->compute_log_change_on_swi(current_model->Ns, type_add,
                                                  type_rem)
                 + new_model->log_likelihood_ - current_model->log_likelihood_
                 - log(q_add_types_[ind][t_ind] / q_p_add_types_tmp_[n_types-1]) \
                 + log(q_add_types_[ind][t_ind_rem] / (q_p_add_types_[ind][n_types-1] - q_add_types_[ind][t_ind]));

      if (r >= 0 || log(rng.rand_01()) <= r){
        ++n_accepted;
        *current_model = *new_model;

        return 1;
      } else {
        *new_model = *current_model;

        return 0;
      }
    }

    // disallow copy constructor and assignment
    Sampler(const Sampler& sampler);
    Sampler& operator=(const Sampler& sampler);

    friend class RaoBlackwellizer;
};

} // namespace bmagwa

#endif /* SAMPLER_HPP_ */
