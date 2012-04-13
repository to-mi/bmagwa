/* BMAGWA software v2.0
 *
 * sampler.hpp
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
#include "discrete_distribution.hpp"

namespace bmagwa {

#define LOG_HALF -0.69314718055994528622676398299518041312694549560546875

// forward declaration
class Sampler;

//! Implements Rao-Blackwellization procedure for computing marginal SNP posteriors.
/*!
 *  Computes simple linear regression for all SNPs against the residual of the
 *  current sampled linear model (except when the SNP in question is in the
 *  current model; then the residual is recomputed without the SNP) and converts
 *  to marginal probability estimate.
 *
 *  Original reference:
 *  Y Guan and M Stephens, Bayesian Variable Selection Regression for
 *  Genome-wide Association Studies, and other Large-Scale Problems. Annals of
 *  Applied Statistics. 5(3): 1780-1815, September 2011.
 */
class RaoBlackwellizer
{
  public:
    RaoBlackwellizer(const DataModel* data_model,
                     const PrecomputedSNPCovariances* pre_xxcov)
    : pre_xxcov(pre_xxcov),
      m_g(data_model->m_g),
      r_cm(data_model->n),
      r_om(data_model->n),
      x(data_model->n)
    {
      x = 1;
    }

    ~RaoBlackwellizer() {}

    void p_raoblackwell(Sampler* sampler, const Vector& y_hat,
                        double *p_r, double **p_r_types);

  private:
    const PrecomputedSNPCovariances* pre_xxcov;
    const size_t m_g;

    Vector r_cm, r_om, x;

    // disallow copy constructor and assignment
    RaoBlackwellizer(const RaoBlackwellizer& rb);
    RaoBlackwellizer& operator=(const RaoBlackwellizer& rb);
};


void compute_exhaustive_modelset(const size_t n_inds, ExhModel* model,
                                 double* model_probabilities, double& max_log_model);
void compute_proposal_probs_for_exh_modelset(const bool& use_types,
                                             const int& n_inds,
                                             const unsigned char* bit_to_normalized_order,
                                             const double* q_add,
                                             const double* q_rem,
                                             const double& z_add,
                                             const double& z_rem,
                                             const size_t& const_loci,
                                             const size_t& m_g,
                                             const double* log_q_add_types,
                                             double* log_prop_probs);

//! Implements the main sampling routines and holds the state of the sampler.
class Sampler
{
  public:
    // TODO: group variables into manageable classes or something
    Sampler(const Options& opts, const Data* data,
            const PrecomputedSNPCovariances* pre_xxcov)
    : sampler_type(opts.sampler_type),
      n_iter(0),
      do_n_iter(opts.do_n_iter),
      n_accepted(0),
      n_rao(opts.n_rao),
      n_rao_burnin(opts.n_rao_burnin),
      n_sample_tau2_and_missing(opts.n_sample_tau2_and_missing),
      flat_proposal_dist(opts.flat_proposal_dist),
      adaptation(opts.adaptation),
      save_beta(opts.save_beta),
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
      q_add_min_(1.0 / (m_g - opts.e_qg)),
      q_rem_(NULL),
      q_rem_min_(1.0 / opts.e_qg),
      dd_add(NULL),
      dd_rem(NULL),
      p_r_(NULL),
      p_r_types_(NULL),
      q_add_types_(NULL),
      q_p_add_types_(NULL),
      max_move_size_(std::min(std::min(m_g, (size_t)255), (size_t)opts.max_move_size)),
      max_SNP_neighborhood_size_(opts.max_SNP_neighborhood_size),
      p_move_size_(opts.p_move_size),
      p_move_size_nbs_(opts.p_move_size_nbs), // near-by-SNP-switch
      p_move_size_nbc_(opts.p_move_size_nbc), // near-by-SNP-state change
      adapt_p_move_size(opts.adapt_p_move_size && max_move_size_ > 1),
      q_p_move_size_(NULL),
      q_p_move_size_nbs_(NULL), // near-by-SNP-switch
      q_p_move_size_nbc_(NULL), // near-by-SNP-state change
      r_move_size_sum_(NULL),
      r_move_size_n_(NULL),
      q_p0_move_size_(NULL),
      acpt_move_size_goal_(opts.p_move_size_acpt_goal),
      move_inds_add(NULL),
      move_inds_rem(NULL),
      move_inds(NULL),
      move_inds_map(NULL),
      move_isadd(NULL),
      movesize(0),
      delay_rejection(std::min((size_t)max_move_size_, opts.delay_rejection)),
      dr_model_probabilities(NULL),
      dr_q_add(NULL),
      dr_q_rem(NULL),
      dr_log_q_add_types(NULL),
      dr_log_p_moves(NULL),
      dr_bit_to_normalized_order(NULL)
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
        if (!data_model->allow_terms((DataModel::ef_t)t)) continue;
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
                        s2_sigma2, opts.nu_tau2, s2_tau2, opts.mu_alpha,
                        opts.use_individual_tau2);

      current_model = new Model(data, data_model, prior, &samplerstats);
      new_model = new Model(data, data_model, prior, &samplerstats);
      exh_model = new ExhModel(data, data_model, prior, &samplerstats);

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
      q_rem_ = new double[m_g];
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

      move_inds_add = new size_t[max_move_size_];
      move_inds_rem = new size_t[max_move_size_];
      move_inds = new size_t[max_move_size_];
      move_inds_map = new int64_t[max_move_size_];
      move_isadd = new bool[max_move_size_];

      // set "prior" values (to avoid possible zeroes case)
      n_acpt_moves_[0] = 1; n_acpt_moves_[1] = 1; n_acpt_moves_[2] = 1;
      n_moves_[0] = 2; n_moves_[1] = 2; n_moves_[2] = 2;
      q_p_move_size_ = new double[max_move_size_];
      r_move_size_sum_ = new double[max_move_size_];
      r_move_size_n_ = new double[max_move_size_];
      q_p0_move_size_ = new double[max_move_size_];
      size_t max_move_size2_ = (max_move_size_ + 1)/2;
      q_p_move_size_nbs_ = new double[max_move_size2_];
      q_p_move_size_nbc_ = new double[max_move_size_];

      Utils::geometric_dist_cdf(max_move_size_, p_move_size_, q_p_move_size_);
      Utils::geometric_dist_cdf(max_move_size2_, p_move_size_nbs_, q_p_move_size_nbs_);
      Utils::geometric_dist_cdf(max_move_size_, p_move_size_nbc_, q_p_move_size_nbc_);
      for (unsigned int i = 0; i < max_move_size_; ++i){
        r_move_size_sum_[i] = 0.0;
        r_move_size_n_[i] = 0.0;
        q_p0_move_size_[i] = 0.0;
      }

      if (delay_rejection > 0){
        dr_model_probabilities = new double[(size_t)std::pow(2, delay_rejection)];
        dr_q_add = new double[delay_rejection];
        dr_q_rem = new double[delay_rejection];
        dr_log_q_add_types = new double[delay_rejection];
        for (int i = 0; i < delay_rejection; ++i)
          dr_log_q_add_types[i] = 0.0;
        dr_log_p_moves = new double[delay_rejection + 1];
        dr_bit_to_normalized_order = new unsigned char[delay_rejection];
      }
    }

    ~Sampler()
    {
      delete[] dr_bit_to_normalized_order;
      delete[] dr_log_p_moves;
      delete[] dr_log_q_add_types;
      delete[] dr_q_add;
      delete[] dr_q_rem;
      delete[] dr_model_probabilities;
      delete[] q_p_move_size_;
      delete[] q_p_move_size_nbs_;
      delete[] q_p_move_size_nbc_;
      delete[] r_move_size_sum_;
      delete[] r_move_size_n_;
      delete[] q_p0_move_size_;
      delete[] move_inds_add;
      delete[] move_inds_rem;
      delete[] move_inds;
      delete[] move_inds_map;
      delete[] move_isadd;

      delete[] q_add_;
      delete[] q_rem_;
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
      delete exh_model;
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
    void initialize_p_proposal_flat()
    {
      double val = prior->e_g() / (double)m_g;
      double val_types = 1.0 / n_types;
      for (size_t i = 0; i < m_g; ++i){
        p_proposal[i] = val;
        if (use_types){
          for (size_t t = 0; t < n_types; ++t)
            p_proposal_types[i][t] = val_types;
        }
      }
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
               basename.substr(0, basename.length() - 1) + "_prior.txt", false);

      prior->print(filename);
    }

    const double get_alpha_value() const {
      return prior->alpha();
    }

    void compute_p_moves(double* p_moves, double* p_moves_cumsum) const;

  private:
#define NMOVES 7
    char sampler_type;
    size_t n_iter, do_n_iter, n_accepted, n_rao, n_rao_burnin;
    size_t n_sample_tau2_and_missing;
    bool flat_proposal_dist;
    bool adaptation, save_beta;
    size_t thin, verbosity;
    std::string basename;

    int level;
    const size_t m_g;
    const size_t n_types;
    const bool use_types;
    const uint32_t seed;

    SamplerStats samplerstats;
    Rand rng;
    Prior* prior;
    DataModel* data_model;
    Model* current_model;
    Model* new_model;
    ExhModel* exh_model;

    double* p_proposal;
    double** p_proposal_types;
    size_t p_proposal_n;

    RaoBlackwellizer *raob;
    double* p_rao;
    double** p_rao_types;
    int p_rao_n;

    //
    double p_moves_[NMOVES], p_moves_cumsum_[NMOVES];

    double* q_add_;
    double q_add_min_;
    double* q_rem_;
    double q_rem_min_;
    DiscreteDistribution* dd_add;
    DiscreteDistribution* dd_rem;
    double* p_r_;

    double** p_r_types_;
    double** q_add_types_;
    double** q_p_add_types_;
    double q_p_add_types_tmp_[5];

    unsigned char max_move_size_;
    size_t max_SNP_neighborhood_size_;
    double p_move_size_, p_move_size_nbs_, p_move_size_nbc_;
    bool adapt_p_move_size;
    double* q_p_move_size_;
    double* q_p_move_size_nbs_;
    double* q_p_move_size_nbc_;

    double* r_move_size_sum_;
    double* r_move_size_n_;
    double* q_p0_move_size_;
    size_t n_acpt_moves_[NMOVES];
    size_t n_moves_[NMOVES];
    double acpt_move_size_goal_;

    size_t* move_inds_add;
    size_t* move_inds_rem;
    size_t* move_inds;
    int64_t* move_inds_map;
    bool* move_isadd;
    unsigned char movesize;

    unsigned char delay_rejection;
    double* dr_model_probabilities;
    double* dr_q_add;
    double* dr_q_rem;
    double* dr_log_q_add_types;
    double* dr_log_p_moves;
    unsigned char* dr_bit_to_normalized_order;

    // methods
    void sample_missing();
    void compute_proposal_dist(const double* p_proposal, const double q_add_min,
                               double* q_add) const;
    void compute_proposal_types_dist(const double*const* p_proposal_types,
                                     const double q_add_min,
                                     double** q_add_types,
                                     double** q_p_add_types) const;

    // preparations for addition of a new term to the model
    // separate from add_term as add_term may be called even if we don't want
    // to sample missing values or tau2
    void prepare_add_new_term_without_missing_sampling(const size_t ind, const DataModel::ef_t type,
                                                       double* inv_tau2_alpha2, Rand& rng)
    {
      if (type == DataModel::AH){
        inv_tau2_alpha2[0] = prior->sample_single_inv_tau2_from_prior_times_alpha2(DataModel::A, rng);
        inv_tau2_alpha2[1] = prior->sample_single_inv_tau2_from_prior_times_alpha2(DataModel::H, rng);
      } else {
        inv_tau2_alpha2[0] = prior->sample_single_inv_tau2_from_prior_times_alpha2(type, rng);
      }
    }
    void prepare_add_new_term(const size_t ind, const DataModel::ef_t type,
                              double* inv_tau2_alpha2, Rand& rng)
    {
      data_model->sample_missing_single(ind, rng);
      prepare_add_new_term_without_missing_sampling(ind, type, inv_tau2_alpha2, rng);
    }



    char do_block_nk()
    {
      movesize = Utils::sample_discrete_naive(q_p_move_size_, max_move_size_, rng) + 1;

      for (unsigned char i = 0; i < movesize; ++i){
        bool unique;
        size_t ind;
        // sample a unique ind
        do {
          unique = true;
          ind = std::floor(m_g * rng.rand_01());
          assert(ind < m_g);
          // check that ind is unique
          for (unsigned char j = 0; j < i; ++j){
            if (ind == move_inds[j]){
              unique = false;
              break;
            }
          }
        } while(!unique);
        move_inds[i] = ind;
      }

      // NOTE: NO support for effect types other than additive
      unsigned char moved = 0;
      double log_mpc = 0.0, log_q_forward = 0.0, log_q_backward = 0.0;
      double inv_tau2_alpha2_tmp[2];

      for (unsigned char i = 0; i < movesize; ++i){
        size_t ind = move_inds[i];
        double gamma_one_prob = std::min(0.99, std::max(0.01, p_proposal[ind]));
        // sample \gamma_i = 0 or 1 according to prior for each
        if (rng.rand_01() < gamma_one_prob){
          // set gamma_i = 1
          log_q_forward += log(gamma_one_prob);

          if (new_model->model_inds[ind] < 0){
            log_mpc += prior->compute_log_change_on_add(new_model->Ns, new_model->size(), DataModel::A);

            prepare_add_new_term(ind, data_model->types[0], inv_tau2_alpha2_tmp, rng);
            new_model->add_term(ind, 0, inv_tau2_alpha2_tmp);
            ++moved;

            log_q_backward += log(1 - gamma_one_prob); // change to 0
          } else {
            log_q_backward += log(gamma_one_prob); // keep at 1
          }
        } else {
          // set gamma_i = 0
          log_q_forward += log(1 - gamma_one_prob);

          int32_t model_ind = new_model->model_inds[ind];
          if (model_ind >= 0){
            log_mpc += prior->compute_log_change_on_rem(new_model->Ns, new_model->size(), DataModel::A);
            new_model->remove_term(model_ind);
            ++moved;

            log_q_backward += log(gamma_one_prob); // change to 1
          } else {
            log_q_backward += log(1 - gamma_one_prob); // keep at 0
          }
        }
      }

      if (moved == 0){
        r_move_size_sum_[movesize - 1] += 0.0;
        r_move_size_n_[movesize - 1] += 1.0;
        return 0;
      }

      assert(std::isfinite(log_q_forward));
      assert(std::isfinite(log_q_backward));

      double log_r = log_q_backward - log_q_forward;
      log_r += log_mpc + new_model->log_likelihood_ - current_model->log_likelihood_;
      assert(!std::isnan(log_r));

      r_move_size_sum_[movesize - 1] += (double)moved / (double)movesize * (log_r >= 0 ? 1.0 : exp(log_r));
      r_move_size_n_[movesize - 1] += 1.0;

      if (log_r >= 0 || log(rng.rand_01()) <= log_r){
        // accepted move
        *current_model = *new_model;

        return moved;
      } else {
        *new_model = *current_model;

        return 0;
      }

    }


    char do_block_ksc()
    {
      movesize = Utils::sample_discrete_naive(q_p_move_size_, max_move_size_, rng) + 1;

      for (unsigned char i = 0; i < movesize; ++i){
        bool unique;
        size_t ind;
        // sample a unique ind
        do {
          unique = true;
          ind = std::floor(m_g * rng.rand_01());
          assert(ind < m_g);
          // check that ind is unique
          for (unsigned char j = 0; j < i; ++j){
            if (ind == move_inds[j]){
              unique = false;
              break;
            }
          }
        } while(!unique);
        move_inds[i] = ind;
      }

      // NOTE: NO support for effect types other than additive
      unsigned char moved = 0, n_gamma1_forward = 0, n_gamma1_backward = 0;
      double log_mpc = 0.0;
      double inv_tau2_alpha2_tmp[2];

      double gamma_one_prob = exp(prior->compute_log_change_on_add(current_model->Ns, current_model->size(), DataModel::A));
      gamma_one_prob = gamma_one_prob / (gamma_one_prob + 1.0);

      for (unsigned char i = 0; i < movesize; ++i){
        size_t ind = move_inds[i];
        // sample \gamma_i = 0 or 1 according to prior for each
        if (rng.rand_01() < gamma_one_prob){
          // set gamma_i = 1
          ++n_gamma1_forward;

          if (new_model->model_inds[ind] < 0){
            log_mpc += prior->compute_log_change_on_add(new_model->Ns, new_model->size(), DataModel::A);

            prepare_add_new_term(ind, data_model->types[0], inv_tau2_alpha2_tmp, rng);
            new_model->add_term(ind, 0, inv_tau2_alpha2_tmp);

            ++moved;
          } else {
            ++n_gamma1_backward; // keep at 1
          }
        } else {
          // set gamma_i = 0
          int32_t model_ind = new_model->model_inds[ind];
          if (model_ind >= 0){
            log_mpc += prior->compute_log_change_on_rem(new_model->Ns, new_model->size(), DataModel::A);
            new_model->remove_term(model_ind);
            ++moved;
            ++n_gamma1_backward; // change to 1
          }
        }
      }

      if (moved == 0){
        r_move_size_sum_[movesize - 1] += 0.0;
        r_move_size_n_[movesize - 1] += 1.0;
        return 0;
      }

      double log_q_forward = n_gamma1_forward * log(gamma_one_prob) + (movesize - n_gamma1_forward) * log(1.0 - gamma_one_prob);
      assert(std::isfinite(log_q_forward));

      double bw_gamma_one_prob = exp(prior->compute_log_change_on_add(new_model->Ns, new_model->size(), DataModel::A));
      bw_gamma_one_prob = bw_gamma_one_prob / (bw_gamma_one_prob + 1.0);

      double log_q_backward = n_gamma1_backward * log(bw_gamma_one_prob) + (movesize - n_gamma1_backward) * log(1.0 - bw_gamma_one_prob);
      assert(std::isfinite(log_q_backward));

      double log_r = log_q_backward - log_q_forward;
      log_r += log_mpc + new_model->log_likelihood_ - current_model->log_likelihood_;
      assert(!std::isnan(log_r));

      r_move_size_sum_[movesize - 1] += (double)moved / (double)movesize * (log_r >= 0 ? 1.0 : exp(log_r));
      r_move_size_n_[movesize - 1] += 1.0;

      if (log_r >= 0 || log(rng.rand_01()) <= log_r){
        // accepted move
        *current_model = *new_model;

        return moved;
      } else {
        *new_model = *current_model;

        return 0;
      }

    }


    char do_block_gibbs()
    {
      assert(max_move_size_ == delay_rejection);
      // maximum move size = cm->size() removals + (m_g - cm->size()) additions
      //                   = m_g
      // maxmovesize parameter is enforced to be equal or less than m_g
      movesize = Utils::sample_discrete_naive(q_p_move_size_, max_move_size_, rng) + 1;

      double inv_tau2_alpha2_tmp[2];

      // special case for movesize == 1
      if (movesize == 1){
        dr_model_probabilities[0] = 1.0;
        dr_model_probabilities[1] = 0.0;

        size_t ind = std::floor(m_g * rng.rand_01());
        assert(ind < m_g);

        int32_t model_ind = new_model->model_inds[ind];
        if (model_ind >= 0){
          new_model->remove_term(model_ind);

          // NOTE: NO support for effect types other than additive!
          dr_model_probabilities[1] += prior->compute_log_change_on_rem(
                        current_model->Ns, current_model->size(), DataModel::A);
        } else {
          // NOTE: NO support for effect types other than additive!
          prepare_add_new_term(ind, data_model->types[0], inv_tau2_alpha2_tmp, rng);
          new_model->add_term(ind, 0, inv_tau2_alpha2_tmp);

          dr_model_probabilities[1] += prior->compute_log_change_on_add(
                        current_model->Ns, current_model->size(), DataModel::A);
        }
        dr_model_probabilities[1] += new_model->log_likelihood() - current_model->log_likelihood();
        // cumulative sum, where dr_model_probabilities[0] = 1.0
        dr_model_probabilities[1] = 1.0 + exp(dr_model_probabilities[1]);

        int sampled_model = Utils::sample_discrete_naive(dr_model_probabilities, 2, rng);

        if (sampled_model == 0){
          // remained in current_model
          *new_model = *current_model;

          // record the proportion of variables involved which changed status
          r_move_size_sum_[movesize - 1] += 0.0;
          r_move_size_n_[movesize - 1] += 1.0;

          return 0;
        } else {
          // moved to new_model
          *current_model = *new_model;

          // record the proportion of variables involved which changed status
          r_move_size_sum_[movesize - 1] += 1.0;
          r_move_size_n_[movesize - 1] += 1.0;

          return 1;
        }
      }

      // following is only executed on movesize > 1
      unsigned long currentmodel_binary = 0;
      unsigned char m_inmodel = 0;
      for (unsigned char i = 0; i < movesize; ++i){
        bool unique;
        size_t ind;
        // sample a unique ind
        do {
          unique = true;
          ind = std::floor(m_g * rng.rand_01());
          assert(ind < m_g);
          // check that ind is unique
          for (unsigned char j = 0; j < i; ++j){
            if (ind == move_inds[j]){
              unique = false;
              break;
            }
          }
        } while(!unique);
        move_inds[i] = ind;

        int32_t model_ind = new_model->model_inds[ind];
        if (model_ind >= 0){
          ++m_inmodel;
          // wasteful to remove and then add back, but exhaustive comp requires
          // that the variables involved are at the back (instead of remove, we
          // could just move the term to be the last one)
          new_model->remove_term(model_ind);

          currentmodel_binary |= (1 << i);
        }
      }

      // add terms to get the full model
      for (unsigned char i = 0; i < movesize; ++i){
        // NOTE: NO support for effect types other than additive!
        prepare_add_new_term(move_inds[i], data_model->types[0], inv_tau2_alpha2_tmp, rng);
        new_model->add_term(move_inds[i], 0, inv_tau2_alpha2_tmp);

        assert(((currentmodel_binary >> i) & 1) == (current_model->model_inds[move_inds[i]] >= 0));
      }
      assert(new_model->size() == current_model->size() - m_inmodel + movesize);

      // do the exhaustive computation
      size_t const_loci = current_model->size() - m_inmodel;
      exh_model->update_to_model(*new_model, const_loci);

      double max_log_model;
      compute_exhaustive_modelset(movesize, exh_model,
                                  dr_model_probabilities, max_log_model);

      // compute cumulative sum
      unsigned long nmodels = std::pow(2, movesize);
      dr_model_probabilities[0] = exp(dr_model_probabilities[0] - max_log_model);
      for (unsigned long i = 1; i < nmodels; ++i){
        dr_model_probabilities[i] = dr_model_probabilities[i-1] +
                                 exp(dr_model_probabilities[i] - max_log_model);
      }

      // sample a model
      unsigned long sampled_model = Utils::sample_discrete(dr_model_probabilities, nmodels, (int)movesize/2 - 3, rng);
      assert(sampled_model < nmodels);

      unsigned char moved = 0;
      if (sampled_model == currentmodel_binary){
        *new_model = *current_model;
      } else {
        // remove terms from the full model to arrive at the sampled model
        for (unsigned char i = 0; i < movesize; ++i){
          if (0 == ((sampled_model >> i) & 1)){
            new_model->remove_term(new_model->model_inds[move_inds[i]]);

            if (1 == ((currentmodel_binary >> i) & 1)) ++moved;
          } else {
            if (0 == ((currentmodel_binary >> i) & 1)) ++moved;
          }
        }

        // update current model
        *current_model = *new_model;
      }

      // record the proportion of variables involved which changed status
      r_move_size_sum_[movesize - 1] += (double)moved / (double)movesize;
      r_move_size_n_[movesize - 1] += 1.0;

      return moved;
    }

    char do_multistep_additions_and_removals()
    {
      // maximum move size = cm->size() removals + (m_g - cm->size()) additions
      //                   = m_g
      // maxmovesize parameter is enforced to be equal or less than m_g
      movesize = Utils::sample_discrete_naive(q_p_move_size_, max_move_size_, rng) + 1;

      double log_q_forward = 0.0, log_q_backward = 0.0, log_mpc = 0.0;
      unsigned char ms_rem = 0;
      prepare_addrem(log_q_forward, ms_rem, movesize);
      do_addrem(log_q_forward, log_q_backward, log_mpc, movesize);
      backward_prepare_addrem(log_q_backward, movesize);

      double log_r = log_q_backward - log_q_forward;
      log_r += log_mpc + new_model->log_likelihood_ - current_model->log_likelihood_;
      assert(!std::isnan(log_r));

      if (delay_rejection == 0)
        r_move_size_sum_[movesize - 1] += (log_r >= 0 ? 1.0 : exp(log_r));
      r_move_size_n_[movesize - 1] += 1.0;

      if (log_r >= 0 || log(rng.rand_01()) <= log_r){
        // accepted move
        *current_model = *new_model;

        if (delay_rejection != 0)
          r_move_size_sum_[movesize - 1] += 1.0;

        return movesize;
      } else {
        // the delayed rejection algorithm can only propose current_model
        // if movesize == 1. Skip DR in this case.
        if (movesize <= delay_rejection && movesize > 1){
          // the number of variables in the model, which are not involved here
          size_t const_loci = current_model->size() - ms_rem;
          unsigned char ms_add = movesize - ms_rem;

          // binary representation of the current model
          double z_add = dd_add->total_w(); // these are for the new_model
          double z_rem = dd_rem->total_w(); // but will be for the "empty" model
          unsigned long newmodel_binary = ((1 << ms_add) - 1);

          // add removed back into new_model
          for (unsigned char i_move = 0; i_move < movesize; ++i_move){
            size_t ind = move_inds[i_move];
            if (move_isadd[i_move]){
              // remove added from zs to make them for the "empty" model
              z_add += q_add_[ind];
              z_rem -= q_rem_[ind];
            } else {
              char t_ind = current_model->type_inds[current_model->model_inds[ind]];

              DataModel::ef_t type = data_model->types[t_ind];
              if (type == DataModel::AH){
                double inv_tau2_alpha2_val[2];
                inv_tau2_alpha2_val[0] = current_model->inv_tau2_alpha2(current_model->x_ind1[current_model->model_inds[ind]]);
                inv_tau2_alpha2_val[1] = current_model->inv_tau2_alpha2(current_model->x_ind2[current_model->model_inds[ind]]);
                new_model->add_term(ind, t_ind, inv_tau2_alpha2_val);
              } else {
                double inv_tau2_alpha2_val = current_model->inv_tau2_alpha2(current_model->x_ind1[current_model->model_inds[ind]]);
                new_model->add_term(ind, t_ind, &inv_tau2_alpha2_val);
              }


            }
          }

          // move_inds_map contains the inds in the order +1,+2,+3,... move
          for (unsigned char i_move = 0; i_move < movesize; ++i_move){
            size_t ind = move_inds_map[i_move];
            size_t model_ind = new_model->model_inds[ind];
            assert((int)model_ind - (int)const_loci >= 0 && (int)model_ind - (int)const_loci < movesize);
            dr_bit_to_normalized_order[model_ind - const_loci] = i_move;
            dr_q_add[i_move] = q_add_[ind];
            dr_q_rem[i_move] = q_rem_[ind];

            if (n_types > 1){
              unsigned char t_ind = new_model->type_inds[new_model->model_inds[ind]];
              dr_log_q_add_types[i_move] = log(q_add_types_[ind][t_ind] / q_p_add_types_[ind][n_types-1]);
            }
          }

          exh_model->update_to_model(*new_model, const_loci);

          double max_log_model;
          compute_exhaustive_modelset(movesize, exh_model,
                                      dr_model_probabilities, max_log_model);
          assert(std::abs(log_r - (log_q_backward - log_q_forward)
              - (dr_model_probabilities[newmodel_binary] - dr_model_probabilities[(~newmodel_binary) & ((1 << movesize) - 1)])) < 1e-10);

          compute_proposal_probs_for_exh_modelset(use_types,
                                                  movesize,
                                                  dr_bit_to_normalized_order,
                                                  dr_q_add,
                                                  dr_q_rem,
                                                  z_add,
                                                  z_rem,
                                                  const_loci,
                                                  m_g,
                                                  dr_log_q_add_types,
                                                  dr_model_probabilities);

          // compute cumulative of proposal distribution and sample
          unsigned long nmodels = std::pow(2, movesize);
          unsigned long mask = nmodels - 1;
          assert((mask >> (movesize - 1) & 1) == 1);
          assert((mask >> movesize & 1) == 0);
          assert(std::abs(dr_model_probabilities[newmodel_binary] - dr_model_probabilities[(~newmodel_binary) & mask] - log_r) < 1e-10);

          double sum = 0.0;
          for (unsigned long i = 0; i < nmodels/2; ++i){
            unsigned int j = (~i) & mask;
            double a = dr_model_probabilities[j] - dr_model_probabilities[i];

            if (a > 0){
              dr_model_probabilities[i] = exp(dr_model_probabilities[j] - max_log_model) * (1 - exp(-a));
              dr_model_probabilities[j] = -1.0; // j is nonzero
            } else {
              dr_model_probabilities[i] = exp(dr_model_probabilities[i] - max_log_model) * (1 - exp(a));
              dr_model_probabilities[j] = 1.0; // i is nonzero
            }
            // cumsum
            dr_model_probabilities[i] += sum;
            sum = dr_model_probabilities[i];
          }

          unsigned long sampled_model = Utils::sample_discrete(dr_model_probabilities, nmodels/2, (int)movesize/2 - 3, rng);
          if (dr_model_probabilities[(~sampled_model) & mask] < 0){
            sampled_model = (~sampled_model) & mask;
          }

          // update models, dd_add/rem, new_dd_add/rem
          assert(sampled_model != newmodel_binary);
          unsigned char moved = 0;
          if (sampled_model == ((~newmodel_binary) & mask)){
            // not moving anywhere

            // no need to add zero
            // r_move_size_sum_[movesize - 1] += 0.0;

            *new_model = *current_model;

            for (unsigned char i_move = 0; i_move < movesize; ++i_move){
              size_t ind = move_inds[i_move];
              if (move_isadd[i_move]){
                dd_add->remdate(ind);
                dd_rem->adddate(ind);
              } else {
                dd_add->adddate(ind);
                dd_rem->remdate(ind);
              }

              assert(dd_add->zeroed(ind) == (current_model->model_inds[ind] >= 0));
              assert(dd_rem->zeroed(ind) == (current_model->model_inds[ind] < 0));
            }
          } else {
            // update dd_add/rem by examining difference to xmodel_binary
            // update new_dd_add/rem by examining difference to (~xmodel_binary & mask)
            for (unsigned char i_move = 0; i_move < movesize; ++i_move){
              if (1 == ((sampled_model >> i_move) & 1)){
                if (0 == ((newmodel_binary >> i_move) & 1)){
                  // sampled model has this one, but new model does not
                  size_t ind = move_inds_map[dr_bit_to_normalized_order[i_move]];
                  dd_add->adddate(ind);
                  dd_rem->remdate(ind);
                } else {
                  // sampled model and new model has this, so it's different to current model
                  ++moved;
                }
              } else {
                if (1 == ((newmodel_binary >> i_move) & 1)){
                  // new model has this one, but sampled model does not
                  size_t ind = move_inds_map[dr_bit_to_normalized_order[i_move]];
                  dd_add->remdate(ind);
                  dd_rem->adddate(ind);
                } else {
                  // sampled model and new model does not have this, but current model does
                  ++moved;
                }
              }
            }

            // check for differences to full_model and update models
            unsigned long removals = mask & sampled_model;
            for (unsigned char i = 0; i < movesize; ++i){
              if ((removals & 1) == 0){
                new_model->remove_term(const_loci + i);
                // removed one so const_loci + i for next removal must account
                // for this. let's do this by decreasing const_loci
                --const_loci;
              }
              removals >>= 1;
            }

            r_move_size_sum_[movesize - 1] += (double)moved / (double)movesize;

            *current_model = *new_model;
          }

          return moved;
        } else {
          *new_model = *current_model;

          // no need to add zero
          //if (delay_rejection != 0)
          //  r_move_size_sum_[movesize - 1] += 0.0;

          for (unsigned char i_move = 0; i_move < movesize; ++i_move){
            size_t ind = move_inds[i_move];
            if (move_isadd[i_move]){
              dd_add->remdate(ind);
              dd_rem->adddate(ind);
            } else {
              dd_add->adddate(ind);
              dd_rem->remdate(ind);
            }

            assert(dd_add->zeroed(ind) == (current_model->model_inds[ind] >= 0));
            assert(dd_rem->zeroed(ind) == (current_model->model_inds[ind] < 0));
          }

          return 0;
        }
      }
    }

    // prepares for the accept/reject step
    void prepare_addrem(double& log_q_forward, unsigned char& n_removals,
                        const unsigned char ms)
    {
      size_t max_adds = m_g - current_model->size();
      size_t max_rems = current_model->size();
      n_removals = 0;

      for (unsigned char i_move = 0; i_move < ms; ++i_move){
        bool sample_add; // add or rem

        if (max_rems == 0) sample_add = true;
        else if (max_adds == 0) sample_add = false;
        else {
          log_q_forward += LOG_HALF;
          sample_add = (rng.rand_01() < 0.5);
        }

        if (sample_add){
          --max_adds;
          move_isadd[i_move] = true;

          // addition
          size_t ind = dd_add->sample();
          move_inds[i_move] = ind;
          move_inds_map[i_move] = ind;
          log_q_forward += log(q_add_[ind]) - log(dd_add->total_w());
          dd_add->adddate(ind);
        } else {
          ++n_removals;
          --max_rems;
          move_isadd[i_move] = false;

          size_t ind = dd_rem->sample();
          move_inds[i_move] = ind;
          move_inds_map[i_move] = -1;
          log_q_forward += log(q_rem_[ind]) - log(dd_rem->total_w());
          dd_rem->adddate(ind);
        }
      }
    }

    // prepares for the accept/reject step
    void backward_prepare_addrem(double& log_q_backward, const unsigned char ms)
    {
      size_t max_adds = m_g - new_model->size();
      size_t max_rems = new_model->size();

      double w_add = dd_add->total_w();
      double w_rem = dd_rem->total_w();
      unsigned char last_rem_pos = ms - 1;
      size_t ind;

      for (unsigned char i_move = 0; i_move < ms; ++i_move){
        if (max_adds > 0 && max_rems > 0) log_q_backward += LOG_HALF;

        // move_isadd in reverse
        if (move_isadd[i_move]){
          // this is removal
          --max_rems;

          while (!move_isadd[last_rem_pos]) --last_rem_pos;
          ind = move_inds_map[last_rem_pos];
          --last_rem_pos;

          log_q_backward += log(q_rem_[ind]) - log(w_rem);
          w_rem -= q_rem_[ind];
        } else {
          // this is addition
          --max_adds;

          ind = move_inds_map[i_move];

          log_q_backward += log(q_add_[ind]) - log(w_add);
          w_add -= q_add_[ind];
        }
      }
    }

    void do_addrem(double& log_q_forward, double& log_q_backward, double& log_mpc, const unsigned char ms)
    {

      size_t ind;
      unsigned char last_map_pos = ms - 1;
      // do removals first
      for (unsigned char i_move = 0; i_move < ms; ++i_move){
        if (move_isadd[i_move]) continue;

        ind = move_inds[i_move];
        dd_add->remdate(ind);

        // complete move_inds_map
        while (move_inds_map[last_map_pos] >= 0) --last_map_pos;
        move_inds_map[last_map_pos] = ind;

        size_t model_ind = new_model->model_inds[ind];
        int t_ind = new_model->type_inds[model_ind];
        DataModel::ef_t type = data_model->types[t_ind];

        // accumulate log prior change (must do before updating the model)
        // (could be done more efficiently in one go?)
        log_mpc += prior->compute_log_change_on_rem(new_model->Ns,
            new_model->size(),
            type);

        // accumulate type proposal p (for backward move)
        if (use_types)
          log_q_backward += log(q_add_types_[ind][t_ind] / q_p_add_types_[ind][n_types-1]);

        // remove term
        new_model->remove_term(model_ind);
      }

      // then additions
      double inv_tau2_alpha2_tmp[2];
      for (unsigned char i_move = 0; i_move < ms; ++i_move){
        if (!move_isadd[i_move]) continue;

        ind = move_inds[i_move];
        dd_rem->remdate(ind);

        int t_ind = 0;
        if (use_types) t_ind = Utils::sample_discrete_naive(q_p_add_types_[ind],
            n_types, rng);
        DataModel::ef_t type = data_model->types[t_ind];

        // accumulate log prior change (must do before updating the model)
        // (could be done more efficiently in one go?)
        log_mpc += prior->compute_log_change_on_add(new_model->Ns,
            new_model->size(),
            type);

        // accumulate type proposal p
        if (use_types)
          log_q_forward += log(q_add_types_[ind][t_ind] / q_p_add_types_[ind][n_types-1]);

        // new model
        // sample missing values
        prepare_add_new_term(ind, type, inv_tau2_alpha2_tmp, rng);
        new_model->add_term(ind, t_ind, inv_tau2_alpha2_tmp);

        //new_model->add_term(ind, t_ind, rng);
      }
    }


    char do_switch_of_nearby_snps()
    {
      size_t full_model_move = m_g - current_model->size();
      size_t max_move_size2_ = (max_move_size_ + 1) / 2;
      unsigned char maxmovesize = std::min(max_move_size2_, (full_model_move+1)/2);
      maxmovesize = std::min((size_t)maxmovesize, current_model->size());
      assert(maxmovesize > 0);
      assert(q_p_move_size_nbs_[maxmovesize-1] > 0.0);
      movesize = Utils::sample_discrete_naive(q_p_move_size_nbs_, maxmovesize, rng)
                 + 1;

      // size of neighborhood (to one direction from SNP to be removed)
      size_t neighborhood_size = std::max((size_t)1, (full_model_move - movesize) / 2);
      neighborhood_size = std::min(max_SNP_neighborhood_size_, neighborhood_size);
      assert(neighborhood_size > 0);

      size_t ind_rem, ind_add, model_ind_rem;
      // find proposals first, then do removals and additions
      for (unsigned char i_move = 0; i_move < movesize; ++i_move){

        // choose model_ind to remove
        model_ind_rem = std::floor(rng.rand_01() * new_model->size());
        ind_rem = new_model->loci[model_ind_rem];
        move_inds_rem[i_move] = ind_rem;

        // find replacement
        // NOTE: we abuse current_model->model_inds to have values >= 0 for those that will be added!
        // This is fine, because no other changes are done to current_model and
        // the model_inds are reverted immediately after this loop
        int j = 0;
        int move = (int)std::floor(rng.rand_01() * neighborhood_size) + 1;
        if (rng.rand_01() < 0.5){
          // go up
          ind_add = ind_rem;
          j = 0;
          while (j < move){
            ++ind_add;
            if (ind_add < m_g){
              if (current_model->model_inds[ind_add] < 0){
                ++j;
              }
            } else {
              break;
            }
          }
          if (ind_add == m_g){
            // hit upper border, must go down instead (assume we cannot hit lower border)
            ind_add = ind_rem;
            move = move + move - j;
            j = 0;
            while (j < move){
              --ind_add;
              while (current_model->model_inds[ind_add] >= 0){
                --ind_add;
              }
              ++j;
            }
          }
        } else {
          // go down
          ind_add = ind_rem;
          j = 0;
          if (ind_add > 0){
            while (j < move){
              --ind_add;
              if (ind_add > 0){
                if (current_model->model_inds[ind_add] < 0){
                  ++j;
                }
              } else {
                break;
              }
            }
          }
          if (j < move){
            // about to hit lower border
            assert(ind_add == 0);
            if (current_model->model_inds[ind_add] < 0){
              ++j;
            }
            if (j < move){
              // hit lower border, must go up
              ind_add = ind_rem;
              move = move + move - j;
              j = 0;
              while(j < move){
                ++ind_add;
                while(current_model->model_inds[ind_add] >= 0){
                  ++ind_add;
                }
                ++j;
              }
            }
          }

        }
        assert(current_model->model_inds[ind_add] < 0 && new_model->model_inds[ind_add] < 0);
        assert(ind_add != ind_rem);
        move_inds_add[i_move] = ind_add;

        // remove term
        char t_ind = new_model->type_inds[model_ind_rem];
        new_model->remove_term(model_ind_rem);
        // add term (cannot add yet)
        //data_model->sample_missing_single(ind_add, rng);
        //new_model->add_term(ind_add, t_ind);
        current_model->model_inds[ind_add] = t_ind; // abuse of model_inds in cmodel!
      }

      double inv_tau2_alpha2_tmp[2];
      for (unsigned char i_move = 0; i_move < movesize; ++i_move){
        // add term
        ind_add = move_inds_add[i_move];

        char t_ind = (char)current_model->model_inds[ind_add];
        DataModel::ef_t type = data_model->types[t_ind];

        prepare_add_new_term(ind_add, type, inv_tau2_alpha2_tmp, rng);
        new_model->add_term(ind_add, t_ind, inv_tau2_alpha2_tmp);

        current_model->model_inds[ind_add] = -1; // revert!
      }

      double log_r = new_model->log_likelihood_ - current_model->log_likelihood_;

      if (log_r >= 0 || log(rng.rand_01()) <= log_r){
        *current_model = *new_model;

        for (unsigned char i_move = 0; i_move < movesize; ++i_move){
          dd_add->remdate(move_inds_rem[i_move]);
          dd_rem->adddate(move_inds_rem[i_move]);

          dd_add->adddate(move_inds_add[i_move]);
          dd_rem->remdate(move_inds_add[i_move]);

          assert(dd_add->zeroed(move_inds_add[i_move]) == (current_model->model_inds[move_inds_add[i_move]] >= 0));
          assert(dd_rem->zeroed(move_inds_add[i_move]) == (current_model->model_inds[move_inds_add[i_move]] < 0));
          assert(dd_add->zeroed(move_inds_rem[i_move]) == (current_model->model_inds[move_inds_rem[i_move]] >= 0));
          assert(dd_rem->zeroed(move_inds_rem[i_move]) == (current_model->model_inds[move_inds_rem[i_move]] < 0));
        }

        return 2 * movesize;
      } else {
        *new_model = *current_model;

        return 0;
      }
    }


    char do_statechange_of_nearby_snps()
    {
      // TODO: does not account for multiple effect types at this time

      // select one SNP from the model and propose changes to near-by SNPs
      // it is assumed that current model size > 0 and m_g > 2!

      // maximum number of proposed changes (m_g - 1 since one is chose as the
      // "center" piece and is not changed)
      size_t maxmovesize = std::min(max_SNP_neighborhood_size_ * 2,
                                    std::min(m_g - 1, (size_t)max_move_size_));

      assert(maxmovesize > 0);
      assert(q_p_move_size_nbc_[maxmovesize-1] > 0.0);

      movesize = Utils::sample_discrete_naive(q_p_move_size_nbc_, maxmovesize, rng)
                 + 1;

      // size of neighborhood (to one direction from SNP to be removed)
      // we'll have the neighborhood wrapping at edges
      size_t neighborhood_size = std::min(max_SNP_neighborhood_size_, m_g/2);
      assert(neighborhood_size > 0);

      // first find the central SNP
      size_t model_ind_c = std::floor(rng.rand_01() * new_model->size());
      size_t ind_c = new_model->loci[model_ind_c];

      size_t* nbpermutation = new size_t[2 * neighborhood_size];

      for (size_t i = 0; i < 2 * neighborhood_size; ++i){
        nbpermutation[i] = i;
      }

      // adjust the area around ind_c to wrap at correct index
      // (i.e., if we count up from ind_c we'll wrap back to values < ind_c when
      // we go over wrap)
      size_t wrap = ind_c + neighborhood_size;
      if (ind_c < neighborhood_size){
        wrap += (neighborhood_size - ind_c);
        assert(wrap < m_g);
      } else {
        if (wrap > m_g - 1) wrap = m_g - 1;
      }

      // then start proposing changes to near-by SNPs
      // find proposals first, then do removals and additions
      for (unsigned char i_move = 0; i_move < movesize; ++i_move){
        // get ind to change
        size_t ind = std::floor(rng.rand_01() * (neighborhood_size * 2 - i_move)) + i_move;
        std::swap(nbpermutation[i_move], nbpermutation[ind]);
      }

      // move to ind_c coords and wrap and divide into additions and removals
      unsigned char ms_rem = 0, ms_add = 0;
      for (unsigned char i_move = 0; i_move < movesize; ++i_move){
        size_t ind = nbpermutation[i_move] + ind_c + 1;
        if (ind > wrap) ind = ind_c - (ind - wrap);

        if (new_model->model_inds[ind] < 0){
          move_inds_add[ms_add++] = ind;
        } else {
          move_inds_rem[ms_rem++] = ind;
        }

      }

      delete[] nbpermutation;


      double log_mpc = 0.0;
      char t_ind;
      DataModel::ef_t type;

      // do removals
      for (unsigned char i_move = 0; i_move < ms_rem; ++i_move){
        int model_ind_rem = new_model->model_inds[move_inds_rem[i_move]];
        type = data_model->types[new_model->type_inds[model_ind_rem]];

        log_mpc += prior->compute_log_change_on_rem(new_model->Ns, new_model->size(), type);

        new_model->remove_term(model_ind_rem);
      }

      double inv_tau2_alpha2_tmp[2];
      // do additions (using the effect type of ind_c for all!)
      t_ind = new_model->type_inds[model_ind_c];
      type = data_model->types[t_ind];
      for (unsigned char i_move = 0; i_move < ms_add; ++i_move){
        size_t ind_add = move_inds_add[i_move];

        log_mpc += prior->compute_log_change_on_add(new_model->Ns, new_model->size(), type);

        prepare_add_new_term(ind_add, type, inv_tau2_alpha2_tmp, rng);
        new_model->add_term(ind_add, t_ind, inv_tau2_alpha2_tmp);

        //new_model->add_term(ind_add, t_ind, rng);
      }

      assert(new_model->model_inds[ind_c] >= 0);

      // accept/reject
      double log_r = log_mpc + new_model->log_likelihood_ - current_model->log_likelihood_;
      log_r += log(current_model->size()) - log(new_model->size());

      if (log_r >= 0 || log(rng.rand_01()) <= log_r){
        *current_model = *new_model;

        for (unsigned char i_move = 0; i_move < ms_add; ++i_move){
          dd_add->adddate(move_inds_add[i_move]);
          dd_rem->remdate(move_inds_add[i_move]);

          assert(dd_add->zeroed(move_inds_add[i_move]));
          assert(!dd_rem->zeroed(move_inds_add[i_move]));
        }
        for (unsigned char i_move = 0; i_move < ms_rem; ++i_move){
          dd_add->remdate(move_inds_rem[i_move]);
          dd_rem->adddate(move_inds_rem[i_move]);

          assert(!dd_add->zeroed(move_inds_rem[i_move]));
          assert(dd_rem->zeroed(move_inds_rem[i_move]));
        }

        return movesize;
      } else {
        // the delayed rejection algorithm can only propose current_model
        // if movesize == 1. Skip DR in this case.
        if (movesize <= delay_rejection && movesize > 1){
          // the number of variables in the model, which are not involved here
          size_t const_loci = current_model->size() - ms_rem;
          assert(const_loci > 0);

          // binary representation of the current model
          unsigned long newmodel_binary = ((1 << ms_add) - 1);

          // add removed back into new_model
          for (unsigned char i_move = 0; i_move < ms_rem; ++i_move){
            size_t ind = move_inds_rem[i_move];
            t_ind = current_model->type_inds[current_model->model_inds[ind]];
            DataModel::ef_t type = data_model->types[t_ind];

            if (type == DataModel::AH){
              double inv_tau2_alpha2_val[2];
              inv_tau2_alpha2_val[0] = current_model->inv_tau2_alpha2(current_model->x_ind1[current_model->model_inds[ind]]);
              inv_tau2_alpha2_val[1] = current_model->inv_tau2_alpha2(current_model->x_ind2[current_model->model_inds[ind]]);
              new_model->add_term(ind, t_ind, inv_tau2_alpha2_val);
            } else {
              double inv_tau2_alpha2_val = current_model->inv_tau2_alpha2(current_model->x_ind1[current_model->model_inds[ind]]);
              new_model->add_term(ind, t_ind, &inv_tau2_alpha2_val);
            }

          }

          exh_model->update_to_model(*new_model, const_loci);

          double max_log_model;
          compute_exhaustive_modelset(movesize, exh_model,
              dr_model_probabilities, max_log_model);
          assert(std::abs(log_r - (log(current_model->size()) - log(new_model->size()-ms_rem))
              - (dr_model_probabilities[newmodel_binary] - dr_model_probabilities[(~newmodel_binary) & ((1 << movesize) - 1)])) < 1e-10);

          unsigned long nmodels = std::pow(2, movesize);
          unsigned long mask = nmodels - 1;
          assert((mask >> (movesize - 1) & 1) == 1);
          assert((mask >> movesize & 1) == 0);

          //add ind_c selection probabilities
          for (unsigned int i = 0; i < nmodels; ++i){
            dr_model_probabilities[i] -= log(const_loci + __builtin_popcount(i));
          }

          // compute cumulative of proposal distribution and sample
          double sum = 0.0;
          for (unsigned long i = 0; i < nmodels/2; ++i){
            unsigned int j = (~i) & mask;
            double a = dr_model_probabilities[j] - dr_model_probabilities[i];

            if (a > 0){
              dr_model_probabilities[i] = exp(dr_model_probabilities[j] - max_log_model) * (1 - exp(-a));
              dr_model_probabilities[j] = -1.0; // j is nonzero
            } else {
              dr_model_probabilities[i] = exp(dr_model_probabilities[i] - max_log_model) * (1 - exp(a));
              dr_model_probabilities[j] = 1.0; // i is nonzero
            }
            // cumsum
            dr_model_probabilities[i] += sum;
            sum = dr_model_probabilities[i];
          }

          unsigned long sampled_model = Utils::sample_discrete(dr_model_probabilities, nmodels/2, (int)movesize/2 - 3, rng);
          if (dr_model_probabilities[(~sampled_model) & mask] < 0){
            sampled_model = (~sampled_model) & mask;
          }

          // update models, dd_add/rem, new_dd_add/rem
          assert(sampled_model != newmodel_binary);
          unsigned char moved = 0;
          if (sampled_model == ((~newmodel_binary) & mask)){
            // not moving anywhere
            *new_model = *current_model;
          } else {
            // update dd_add/rem by examining difference to newmodel_binary
            for (unsigned char i_move = 0; i_move < movesize; ++i_move){
              if (1 == ((sampled_model >> i_move) & 1)){
                if (1 == ((newmodel_binary >> i_move) & 1)){
                  // sampled model and new model has this, so it's different to current model
                  ++moved;

                  // this must have been an addition
                  size_t ind = move_inds_add[i_move];
                  dd_add->adddate(ind);
                  dd_rem->remdate(ind);
                }
              } else {
                if (0 == ((newmodel_binary >> i_move) & 1)){
                  // sampled model and new model does not have this, but current model does
                  ++moved;

                  // this must have been a removal
                  size_t ind = move_inds_rem[i_move - ms_add];
                  dd_add->remdate(ind);
                  dd_rem->adddate(ind);
                }
              }
            }

            // check for differences to full_model and update models
            unsigned long removals = mask & sampled_model;
            for (unsigned char i = 0; i < movesize; ++i){
              if ((removals & 1) == 0){
                new_model->remove_term(const_loci + i);
                // removed one so const_loci + i for next removal must account
                // for this. let's do this by decreasing const_loci
                --const_loci;
              }
              removals >>= 1;
            }

            *current_model = *new_model;
          }

          return moved;

        } else {
          *new_model = *current_model;

          return 0;
        }
      }
    }


    char do_switch_of_types()
    {
      movesize = 1;
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
      double inv_tau2_alpha2_tmp[2];
      new_model->remove_term(model_ind);
      prepare_add_new_term_without_missing_sampling(ind, type_add, inv_tau2_alpha2_tmp, rng);
      new_model->add_term(ind, t_ind, inv_tau2_alpha2_tmp);

      // Delta posterior - proposal prob + reverse prob
      double r = prior->compute_log_change_on_swi(current_model->Ns, type_add,
                                                  type_rem)
                 + new_model->log_likelihood_ - current_model->log_likelihood_
                 - log(q_add_types_[ind][t_ind] / q_p_add_types_tmp_[n_types-1]) \
                 + log(q_add_types_[ind][t_ind_rem] / (q_p_add_types_[ind][n_types-1] - q_add_types_[ind][t_ind]));

      if (r >= 0 || log(rng.rand_01()) <= r){
        *current_model = *new_model;

        return 1;
      } else {
        *new_model = *current_model;

        return 0;
      }
    }

    void adapt_p_move_size_acptrate(double& p_ms, size_t& n_acpt_moves,
                                    size_t& n_moves, const double t)
    {
      // acceptance probability coercion
      double acpt_rate = (double)n_acpt_moves / (double)n_moves;
      p_ms = p_ms + (acpt_move_size_goal_ - acpt_rate) / t;
      p_ms = std::max(std::min(p_ms, 0.99), 0.01);

      // set "prior" values (to avoid possible zeroes case)
      n_acpt_moves = 1; n_moves = 2;
    }


    void adapt_p_move_size_jd_mb(double &p_ms, const double* r_ms,
                                 const double* r_n, const unsigned char max_ms)
    {
      double q = 1 - p_ms;
      double zp = 1 - std::pow(q, max_ms); // normalization due to truncation
      double q_p = p_ms / zp;
      for (unsigned int i = 0; i < max_ms; ++i){
        q_p0_move_size_[i] += q_p;
        q_p *= q;
      }

      double p_best = 0.01, p = 0.01;
      double h_best = -INFINITY, h;
      double step = 0.98 / 49.0;
      double ratio = 0.0;
      // Evaluate at 50 grid points from 0.01 to 0.99 and take the best one.
      // This also ensures that p_ms cannot go outside of [0.01,0.99].
      for (int j = 0; j < 50; ++j){
        // no need for normalization of truncation here
        h = 0.0;
        zp = 0.0; // normalization of the weights
        q = 1.0 - p;
        q_p = p;

        for (unsigned int i = 0; i < max_ms; ++i){
          ratio = q_p / q_p0_move_size_[i];
          h += (i+1) * r_ms[i] * ratio;
          zp += r_n[i] * ratio;
          q_p *= q;
        }
        h /= zp;

        if (h > h_best){
          h_best = h;
          p_best = p;
        }

        p += step;
      }

      p_ms = p_best;
    }

    // disallow copy constructor and assignment
    Sampler(const Sampler& sampler);
    Sampler& operator=(const Sampler& sampler);

    friend class RaoBlackwellizer;
};

} // namespace bmagwa

#endif /* SAMPLER_HPP_ */
