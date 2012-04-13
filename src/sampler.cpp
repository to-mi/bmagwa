/* BMAGWA software v2.0
 *
 * sampler.cpp
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


#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include "sampler.hpp"
#include "utils.hpp"

namespace bmagwa {

void RaoBlackwellizer::p_raoblackwell(Sampler* sampler,
                                      const Vector& y_hat,
                                      double *p_r, double **p_r_types)
{
  const Model* cmodel = sampler->current_model;
  const DataModel* data_model = sampler->data_model;
  const Prior* prior = sampler->prior;
  Rand& rng = sampler->rng;

  const double sigma2_times_2 = 2 * cmodel->sigma2;

  // set up types
  const size_t n_types = data_model->n_types;
  const bool use_types = (n_types > 1);

  // residual (r = y - y_hat)
  r_cm.set_to_vector_diff(data_model->y(), y_hat);
  double mr_cm = r_cm.mean();

  // log model prior
  int Ns[5];
  memcpy(Ns, cmodel->Ns, 5 * sizeof(int));
  size_t n_loci = cmodel->size();

  double lmp_g_add[5], lmp_g_rem[5][5];
  for (int t = 0; t < 5; ++t){
    if (data_model->allow_types((DataModel::ef_t)t)){
      if ((size_t)(Ns[t]+1) <= m_g){
        lmp_g_add[t] = sampler->prior->compute_log_change_on_add(Ns, n_loci,
                                               static_cast<DataModel::ef_t>(t));
      }
      Ns[t]--; // remove one
      n_loci--;
      if (Ns[t] >= 0){
        for (int j = 0; j < 5; ++j){
          if (data_model->allow_types((DataModel::ef_t)j)){
            lmp_g_rem[t][j] = sampler->prior->compute_log_change_on_add(Ns,
                                       n_loci, static_cast<DataModel::ef_t>(j));
          }
        }
      }
      Ns[t]++; // back to normal
      n_loci++;
    }
  }

  double inv_tau2_alpha2_vals[4];
  if (!prior->use_individual_tau2){
    for (int t = 0; t < 4; ++t){
      DataModel::ef_t te = (DataModel::ef_t)t;
      if (data_model->allow_terms(te)){
        inv_tau2_alpha2_vals[t] = prior->get_inv_tau2_alpha2_value(te);
      }
    }
  }

  double* pre_xx_tmp = new double[pre_xxcov->offset];
  // main loop
  for (size_t i = 0; i < m_g; ++i){
    Vector* residual;
    double mr, det, exp_term;
    double* lmp;
    double sum_log_Q;

    int32_t model_ind = cmodel->model_inds[i];

    // sample tau2 values for different term types
    if (prior->use_individual_tau2){
      for (int t = 0; t < 4; ++t){
        DataModel::ef_t te = (DataModel::ef_t)t;
        if (data_model->allow_terms(te)){
          inv_tau2_alpha2_vals[t] = prior->sample_single_inv_tau2_from_prior_times_alpha2(te, rng);
        }
      }
    }

    // prepare residual and choose log prior change from tabulated values
    if (model_ind < 0){
      // not in model

      residual = &r_cm;
      mr = mr_cm;
      lmp = lmp_g_add;
    } else {
      // in model

      unsigned int x_ind = cmodel->x_ind1[model_ind];
      DataModel::ef_t type = data_model->types[cmodel->type_inds[model_ind]];

      // override sampled tau2 for this type
      if (prior->use_individual_tau2){
        if (type == DataModel::AH){
          inv_tau2_alpha2_vals[DataModel::A] = cmodel->inv_tau2_alpha2(x_ind);
          inv_tau2_alpha2_vals[DataModel::H] = cmodel->inv_tau2_alpha2(cmodel->x_ind2[model_ind]);
        } else {
          inv_tau2_alpha2_vals[type] = cmodel->inv_tau2_alpha2(x_ind);
        }
      }

      /* Need to compute residual without SNP i:
       * 1) r_cm = y - y_hat
       * 2) y_hat' = y_hat - x_rem * beta_rem
       * 3) r_om = y - y_hat' = y - y_hat + x_rem * beta_rem
       * => r_om = r_cm + x_rem * beta_rem
       */
      r_om = r_cm;
      r_om.add_vector(cmodel->x.column(x_ind), cmodel->beta(x_ind));
      if (type == DataModel::AH){
        // two terms
        x_ind = cmodel->x_ind2[model_ind];
        r_om.add_vector(cmodel->x.column(x_ind), cmodel->beta(x_ind));
      }

      residual = &r_om;
      mr = r_om.mean();

      lmp = lmp_g_rem[type];
    }

    double max_types = -INFINITY;
    // compute rao-blackwellization for each type
    std::memcpy(pre_xx_tmp, pre_xxcov->xx + i * pre_xxcov->offset,
                pre_xxcov->offset * sizeof(double));
    data_model->update_prexx_cov(i, pre_xx_tmp);

    for (size_t t = 0; t < n_types; ++t){
      DataModel::ef_t type = data_model->types[t];

      if (type == DataModel::AH){
        double sa = pre_xx_tmp[0];
        double va = pre_xx_tmp[1] + inv_tau2_alpha2_vals[DataModel::A];
        double sh = pre_xx_tmp[2];
        double vh = pre_xx_tmp[3] + inv_tau2_alpha2_vals[DataModel::H];
        double vah = pre_xx_tmp[pre_xxcov->offset - 1];

        sum_log_Q = -log(inv_tau2_alpha2_vals[DataModel::A])
                    - log(inv_tau2_alpha2_vals[DataModel::H]);

        data_model->get_genotypes_additive(i, x);
        double rxa = Vector::dotproduct(x, *residual) - sa * mr;
        data_model->get_genotypes_heterozygous(i, x);
        double rxh = Vector::dotproduct(x, *residual) - sh * mr;

        det = va * vh - vah * vah;

        exp_term = rxa * rxa * vh - 2.0 * rxa * rxh * vah + rxh * rxh * va;
        exp_term /= det;
      } else {
        size_t offset = pre_xxcov->offset_type[type];

        double s = pre_xx_tmp[offset];
        det = pre_xx_tmp[offset+1] + inv_tau2_alpha2_vals[type];

        sum_log_Q = -log(inv_tau2_alpha2_vals[type]);

        (data_model->*data_model->bmagwa::DataModel::get_genotypes[type])(i, x);
        double rx = Vector::dotproduct(x, *residual) - s * mr;

        exp_term = (rx * rx) / det;
      }

      double p = (exp_term) / sigma2_times_2
                 - 0.5 * (log(det) + sum_log_Q)
                 + *(lmp + type);

      if (use_types){
        if (p > max_types) max_types = p;
        p_r_types[i][t] = p;
      } else {
        p = exp(p);
        if (std::isfinite(p)){
          p_r[i] = p / (1 + p);
        } else {
          p_r[i] = 1;
        }
      }
    } // end of loop over types (t)
    if (use_types){
      // handle infinite case of max_types
      if (!std::isfinite(max_types)){
        std::cout << "Warning: raob max_types is not finite." << std::endl;
        if (std::isnan(max_types)){
          throw std::runtime_error("raob max_types is nan");
        }
        if (max_types > 0){
          //infinite
          double sum_types = 0;
          for (size_t t = 0; t < n_types; ++t){
            p_r_types[i][t] = (!std::isfinite(p_r_types[i][t]) &&
                               (p_r_types[i][t] > 0)) ? 1 : 0;
            sum_types += p_r_types[i][t];
          }
          // normalize types distribution
          for (size_t t = 0; t < n_types; ++t){
            p_r_types[i][t] /= sum_types;
          }
          p_r[i] = 1;
        } else {
          //-infinite
          for (size_t t = 0; t < n_types; ++t){
            p_r_types[i][t] = 1.0/n_types;
          }
          p_r[i] = 0;
        }

        continue;
      }

      // finite case of max_types
      double sum_types = 0;
      // normalize by max_types to avoid numerical issues
      for (size_t t = 0; t < n_types; ++t){
        p_r_types[i][t] = exp(p_r_types[i][t] - max_types);
        sum_types += p_r_types[i][t];
      }
      // normalize types distribution
      for (size_t t = 0; t < n_types; ++t){
        p_r_types[i][t] /= sum_types;
      }
      // compute prob for variant being associated
      sum_types = exp(log(sum_types) + max_types);
      if (std::isfinite(sum_types)){
        p_r[i] = sum_types / (1 + sum_types);
      } else {
        p_r[i] = 1;
      }
    } // end of use_types if
  } // end of main loop (i)
  delete[] pre_xx_tmp;
} // end of p_raoblackwell method


void Sampler::sample_missing()
{
  assert(current_model->size() == new_model->size());

  const size_t n = data_model->n;
  const double sigma2_times_2 = current_model->sigma2 * 2;
  double lprior[3], likelihood[3];

  // compute y_hat and residual
  Vector y_hat(n);
  const Vector& _beta = current_model->beta;
  y_hat.set_to_product(current_model->x, _beta, false);

  Vector residual(n);
  residual.set_to_vector_diff(data_model->y(), y_hat);
  double r2 = 0;
  for (size_t i = 0; i < residual.length(); ++i)
    r2 += residual(i) * residual(i);

  // sample new values
  for (size_t i = 0; i < current_model->size(); ++i){
    uint32_t ind = current_model->loci[i];

    // skip if no missing data
    if (data_model->miss_loc()[ind] == NULL) continue;

    DataModel::ef_t type = data_model->types[current_model->type_inds[i]];
    unsigned int x_ind1 = current_model->x_ind1[i];
    unsigned int x_ind2 = (type == DataModel::AH) ? current_model->x_ind2[i] : 0;

    VectorView tmp_x1(current_model->x.column(x_ind1));
    VectorView tmp_x2(current_model->x.column(x_ind2));

    // compute log prior
    lprior[0] = log(data_model->miss_prior()[ind][0]);
    lprior[1] = log(data_model->miss_prior()[ind][1]
                    - data_model->miss_prior()[ind][0]);
    lprior[2] = log(data_model->miss_prior()[ind][2]
                    - data_model->miss_prior()[ind][1]);

    for (size_t j = 1; j <= data_model->miss_loc()[ind][0]; ++j){
      size_t i_miss = data_model->miss_loc()[ind][j];

      // note the possible change in AH term also
      double old_term2 = residual(i_miss) * residual(i_miss);
      double new_term;

      switch (type){
        case DataModel::A:
          for (int k = 0; k < 3; ++k){
            new_term = residual(i_miss)
                   - _beta(x_ind1) * ((double)k - tmp_x1(i_miss));
            likelihood[k] = lprior[k] - (r2 + new_term * new_term - old_term2) / sigma2_times_2;
          }
          break;
        case DataModel::H:
          for (int k = 0; k < 3; ++k){
            new_term = residual(i_miss) - _beta(x_ind1) * ((double)(k == 1) - tmp_x1(i_miss));
            likelihood[k] = lprior[k] - (r2 + new_term * new_term - old_term2) / sigma2_times_2;
          }
          break;
        case DataModel::R:
          for (int k = 0; k < 3; ++k){
            new_term = residual(i_miss) - _beta(x_ind1) * ((double)(k == 2) - tmp_x1(i_miss));
            likelihood[k] = lprior[k] - (r2 + new_term * new_term - old_term2) / sigma2_times_2;
          }
          break;
        case DataModel::D:
          for (int k = 0; k < 3; ++k){
            new_term = residual(i_miss) - _beta(x_ind1) * ((double)(k > 0) - tmp_x1(i_miss));
            likelihood[k] = lprior[k] - (r2 + new_term * new_term - old_term2) / sigma2_times_2;
          }
          break;
        case DataModel::AH:
          for (int k = 0; k < 3; ++k){
            new_term = residual(i_miss) - _beta(x_ind1) * ((double)k - tmp_x1(i_miss));
            new_term -= _beta(x_ind2) * ((double)(k == 1) - tmp_x2(i_miss));
            likelihood[k] = lprior[k] - (r2 + new_term * new_term - old_term2) / sigma2_times_2;
          }
          break;
      }


      // cumsum for sampling (normalize by likelihood[0] for numerical stability)
      likelihood[1] -= likelihood[0];
      likelihood[2] -= likelihood[0];
      likelihood[0] = 1; // exp(0) = 1
      likelihood[1] = exp(likelihood[1]) + likelihood[0];
      likelihood[2] = exp(likelihood[2]) + likelihood[1];

      // need to update AD term also here! (and y_hat and residual, r2)
      int k = Utils::sample_discrete_naive(likelihood, 3, rng);
      data_model->miss_val()[ind][j] = (char)k;

      switch (type){
        case DataModel::A:
          y_hat(i_miss) += _beta(x_ind1) * ((double)k - tmp_x1(i_miss));
          tmp_x1(i_miss) = (double)k;
          break;
        case DataModel::H:
          k = (k == 1);
          y_hat(i_miss) += _beta(x_ind1) * ((double)k - tmp_x1(i_miss));
          tmp_x1(i_miss) = (double)k;
          break;
        case DataModel::R:
          k = (k == 2);
          y_hat(i_miss) += _beta(x_ind1) * ((double)k - tmp_x1(i_miss));
          tmp_x1(i_miss) = (double)k;
          break;
        case DataModel::D:
          k = (k > 0);
          y_hat(i_miss) += _beta(x_ind1) * ((double)k - tmp_x1(i_miss));
          tmp_x1(i_miss) = (double)k;
          break;
        case DataModel::AH:
          y_hat(i_miss) += _beta(x_ind1) * ((double)k - tmp_x1(i_miss));
          tmp_x1(i_miss) = (double)k;
          k = (k == 1);
          y_hat(i_miss) += _beta(x_ind2) * ((double)k - tmp_x2(i_miss));
          tmp_x2(i_miss) = (double)k;
          break;
      }

      residual(i_miss) = data_model->y()(i_miss) - y_hat(i_miss);
      r2 += residual(i_miss) * residual(i_miss) - old_term2;
    } // end of miss_loc loop (j)
  } // end of sample new values loop (i)

  // update xx, xy
  for (size_t i = 0; i < current_model->size(); ++i){
    uint32_t ind = current_model->loci[i];
    DataModel::ef_t type =data_model->types[current_model->type_inds[i]];

    // skip if no missing data
    if (data_model->miss_loc()[ind] == NULL) continue;

    // term of non-AH or first term of AH
    {
      unsigned int x_ind = current_model->x_ind1[i];
      VectorView xnew(current_model->x.column(x_ind));
      VectorView xold(new_model->x.column(x_ind));

      for (size_t j = 1; j <= data_model->miss_loc()[ind][0]; ++j){
        size_t i_miss = data_model->miss_loc()[ind][j];

        // x' * y
        current_model->xy(x_ind) += data_model->y()(i_miss) * (xnew(i_miss) - xold(i_miss));

        // x' * x (upper part is only accessed)
        // if x_ij and x_kj have both changed, need to be wary of updating (xx)_ij twice
        unsigned int k = 0;
        for (; k < x_ind; ++k){
          //if (k < data_model->m_e || data->g[cmodel->x_ind_to_g_ind[k] * data->n + i_miss] >= 0){
          if (k < data_model->m_e || !data_model->genotype_missing(i_miss, current_model->x_ind_to_g_ind[k])){
            current_model->xx(k, x_ind) += xnew(i_miss) * current_model->x(i_miss, k) - xold(i_miss) * new_model->x(i_miss, k);
          }
        }
        for (; k < current_model->xx.length(); ++k){
          current_model->xx(x_ind, k) += xnew(i_miss) * current_model->x(i_miss, k) - xold(i_miss) * new_model->x(i_miss, k);
        }
      }
    }
    // second term of AH
    if (type == DataModel::AH){
      unsigned int x_ind = current_model->x_ind2[i];
      VectorView xnew(current_model->x.column(x_ind));
      VectorView xold(new_model->x.column(x_ind));

      for (size_t j = 1; j <= data_model->miss_loc()[ind][0]; ++j){
        size_t i_miss = data_model->miss_loc()[ind][j];

        // x' * y
        current_model->xy(x_ind) += data_model->y()(i_miss) * (xnew(i_miss) - xold(i_miss));

        // x' * x (upper part is only accessed)
        unsigned int k = 0;
        for (; k < x_ind; ++k){
          if (k < data_model->m_e || !data_model->genotype_missing(i_miss, current_model->x_ind_to_g_ind[k])){
            current_model->xx(k, x_ind) += xnew(i_miss) * current_model->x(i_miss, k) - xold(i_miss) * new_model->x(i_miss, k);
          }
        }
        for (; k < current_model->xx.length(); ++k){
          current_model->xx(x_ind, k) += xnew(i_miss) * current_model->x(i_miss, k) - xold(i_miss) * new_model->x(i_miss, k);
        }
      }
    }
  }

  current_model->mu_beta_computed = false;
} // sample_missing

void Sampler::compute_p_moves(double* p_moves, double* p_moves_cumsum) const
{
  // TODO: should make this into an user defined option...
  switch (sampler_type){
    case 0: // PMV
      if (use_types){
        p_moves[0] = 0.8;
        p_moves[1] = 0.1;
        p_moves[2] = 0.0;
        p_moves[3] = 0.1;
        p_moves[4] = 0.0;
        p_moves[5] = 0.0;
        p_moves[6] = 0.0;
      } else {
        p_moves[0] = 0.7;
        p_moves[1] = 0.15;
        p_moves[2] = 0.15;
        p_moves[3] = 0.0;
        p_moves[4] = 0.0;
        p_moves[5] = 0.0;
        p_moves[6] = 0.0;
      }
      break;
    case 1: // NK
      p_moves[0] = 0.0;
      p_moves[1] = 0.15;
      p_moves[2] = 0.15;
      p_moves[3] = 0.0;
      p_moves[4] = 0.7;
      p_moves[5] = 0.0;
      p_moves[6] = 0.0;
      break;
    case 2: // KSC
      p_moves[0] = 0.0;
      p_moves[1] = 0.15;
      p_moves[2] = 0.15;
      p_moves[3] = 0.0;
      p_moves[4] = 0.0;
      p_moves[5] = 0.7;
      p_moves[6] = 0.0;
      break;
    case 3: // Gibbs
      p_moves[0] = 0.0;
      p_moves[1] = 0.15;
      p_moves[2] = 0.15;
      p_moves[3] = 0.0;
      p_moves[4] = 0.0;
      p_moves[5] = 0.0;
      p_moves[6] = 0.7;
      break;
    default:
      // this should not happen ever if options works...
      throw std::runtime_error("Error! Sampler type is invalid.");
  }

  // cumulative sum
  p_moves_cumsum[0] = p_moves[0];
  for (int i = 1; i < 7; ++i){
    p_moves_cumsum[i] = p_moves_cumsum[i-1] + p_moves[i];
  }
}

void Sampler::compute_proposal_dist(const double* p_proposal,
                                    const double q_add_min,
                                    double* q_add) const
{
  for (size_t i = 0; i < m_g; ++i){
    q_add[i] = std::max(p_proposal[i], q_add_min);
  }
}

void Sampler::compute_proposal_types_dist(const double*const* p_proposal_types,
                                          const double q_add_min,
                                          double** q_add_types,
                                          double** q_p_add_types) const
{
  double val;

  for (size_t i = 0; i < m_g; ++i){
    val = std::max(p_proposal_types[i][0], q_add_min);
    q_add_types[i][0] = val;
    q_p_add_types[i][0] = val;

    for (size_t t = 1; t < n_types; ++t){
      val = std::max(p_proposal_types[i][t], q_add_min);
      q_add_types[i][t] = val;
      q_p_add_types[i][t] = q_p_add_types[i][t - 1] + val;
    }
  }
}


#define OPENFILE(item) std::ofstream* f_ ## item = Utils::openfile(basename + "_" + #item + ".txt", false)
#define OPENBFILE(item) std::ofstream* f_ ## item = Utils::openfile(basename + "_" + #item + ".dat", true)
#define CLOSEFILE(item) {(f_ ## item)->close(); delete f_ ## item;}

void Sampler::sample()
{
  if (!std::isfinite(current_model->log_likelihood_))
    throw std::runtime_error("Sampler cannot start from non-finite likelihood.");

  double pves[3]; pves[0] = 0.0; pves[1] = 0.0; pves[2] = 0.0;

  double k_move_size = 1000.0 / n_rao_burnin;

  // open files
  OPENFILE(log);
  *f_log << std::setprecision(3) << std::fixed;
  OPENBFILE(loci);
  OPENBFILE(modelsize);
  OPENBFILE(jumpdistance);
  OPENBFILE(log_likelihood);
  OPENBFILE(log_prior);
  OPENBFILE(move_type);
  OPENBFILE(move_size);
  OPENBFILE(pve);
  OPENBFILE(alpha);
  /*std::ofstream* f_tau2[] = {NULL, NULL, NULL, NULL};
  for (int i = 0; i < 4; ++i){
    if (data_model->allow_terms((DataModel::ef_t)i)){
      f_tau2[i] =  Utils::openfile(basename + "_tau2_" + DataModel::ef_name[i]
                                                                + ".dat", true);
    }
  }*/
  //OPENBFILE(tau2);
  OPENBFILE(sigma2);
  std::ofstream* f_types = NULL;
  if (use_types) f_types = Utils::openfile(basename + "_types.dat", true);
  std::ofstream* f_beta = NULL;
  if (save_beta) f_beta = Utils::openfile(basename + "_beta.dat", true);

  *f_log << "output " << basename << std::endl
         << "seed " << seed << std::endl
         << "type (PMV = 0, NK = 1, KSC = 2, G = 3) " << (int)sampler_type << std::endl
         << "----------------------------------------------------" << std::endl;

  // models need to be initialized to same
  *new_model = *current_model;

  // compute proposal dist
  compute_proposal_dist(p_proposal, q_add_min_, q_add_);
  for (size_t i = 0; i < m_g; ++i){
    q_rem_[i] = std::max(1 - p_proposal[i], q_rem_min_);
  }
  dd_add = new DiscreteDistribution(true, q_add_, m_g, rng);
  dd_rem = new DiscreteDistribution(true, q_rem_, m_g, rng);

  // zero all rems
  for (size_t i = 0; i < m_g; ++i) {
    dd_rem->adddate(i);
  }
  if (use_types){
    compute_proposal_types_dist(p_proposal_types, q_add_min_,
                                q_add_types_, q_p_add_types_);
  }

  compute_p_moves(p_moves_, p_moves_cumsum_);

  // sample beta and sigma2 so that tau2 can be sampled
  current_model->sample_beta_sigma2(rng);
  sample_missing(); // this is not really needed here?
  prior->sample_alpha_and_tau2(current_model, data_model->y(), rng);
  // compute likelihood after missing and tau2 sampling
  current_model->compute_log_likelihood();
  // update new_model to be equal to current
  *new_model = *current_model;

  // start timing
  samplerstats.reset_timer();

  size_t _n_iter = n_iter + do_n_iter;
  for (size_t iter = n_iter; iter < _n_iter; ++iter){
    // sample missing genotypes (requires current and new models to be the same)
    if ((iter + 1) % n_sample_tau2_and_missing == 0){
      sample_missing();

      prior->sample_alpha_and_tau2(current_model, data_model->y(), rng);
      // recompute likelihood and update new model to current
      current_model->compute_log_likelihood();
      *new_model = *current_model;
    }

    // select move
    unsigned char move = (unsigned char)Utils::sample_discrete_naive(p_moves_cumsum_, NMOVES, rng);
    unsigned char jumpdistance = 0;

    assert(std::isfinite(current_model->log_likelihood_));
    samplerstats.mhmove_tic();
    switch (move){
      case 0: // additions or removals
        jumpdistance = do_multistep_additions_and_removals();
        break;
      case 1: // SNP switch
        // require at least two free SNPs
        if (current_model->size() > 0 && current_model->size() < (m_g - 1)){
          jumpdistance = do_switch_of_nearby_snps();
        }
        break;
      case 2: //
        if (current_model->size() > 0 && m_g > 2){
          jumpdistance = do_statechange_of_nearby_snps();
        }
        break;
      case 3: // type switch
        if (current_model->size() > 0){
          jumpdistance = do_switch_of_types();
        }
        break;
      case 4: // Nott & Kohn
        jumpdistance = do_block_nk();
        break;
      case 5: // Kohn & Smith & Chan
        jumpdistance = do_block_ksc();
        break;
      case 6: // block Gibbs
        jumpdistance = do_block_gibbs();
        break;
      default:
        throw std::logic_error("Move is not in 0...6");
    }
    samplerstats.mhmove_toc();
    ++n_moves_[move];
    n_acpt_moves_[move] += jumpdistance > 0;
    n_accepted += jumpdistance > 0;
    assert(std::isfinite(current_model->log_likelihood_));

    if ((iter + 1) % thin == 0 || (iter + 2) % n_sample_tau2_and_missing == 0){
      // note: needs to be done before Model_pve and before sampling tau2
      current_model->sample_beta_sigma2(rng);
    }

    // collect sample
    f_jumpdistance->write(reinterpret_cast<const char*>(&jumpdistance), sizeof(jumpdistance));
    f_move_type->write(reinterpret_cast<const char*>(&move), sizeof(move));
    f_move_size->write(reinterpret_cast<const char*>(&movesize),
                       sizeof(movesize));
    if ((iter + 1) % thin == 0){
      uint32_t n_loci = current_model->size();
      f_modelsize->write(reinterpret_cast<const char*>(&n_loci),
                         sizeof(n_loci));
      f_loci->write(reinterpret_cast<const char*>(&current_model->loci[0]),
                    n_loci * sizeof(current_model->loci[0]));
      if (use_types){
        f_types->write(
                    reinterpret_cast<const char*>(&current_model->type_inds[0]),
                    n_loci * sizeof(current_model->type_inds[0]));
      }
      f_log_likelihood->write(
                  reinterpret_cast<const char*>(&current_model->log_likelihood_),
                  sizeof(current_model->log_likelihood_));
      double log_prior = prior->compute_log_model(current_model->Ns);
      f_log_prior->write(
                  reinterpret_cast<const char*>(&log_prior), sizeof(log_prior));
      f_sigma2->write(reinterpret_cast<const char*>(&current_model->sigma2),
                      sizeof(current_model->sigma2));
      if (save_beta){
        f_beta->write(reinterpret_cast<const char*>(current_model->beta.data()),
          current_model->beta.length() * sizeof(current_model->beta.data()[0]));
      }

      // note sample_beta_sigma2 needs to be done before this!
      // also: this call should leave current_model->y_hat to be correct for
      //       Rao-Blackwellization below
      current_model->compute_pve(pves);
      f_pve->write(reinterpret_cast<const char*>(pves), 3 * sizeof(pves[0]));
      f_alpha->write(reinterpret_cast<const char*>(&prior->alpha()),
                     sizeof(prior->alpha()));
      /*for (int i = 0; i < 4; ++i){
        if (data_model->allow_terms((DataModel::ef_t)i)){
          (f_tau2[i])->write(reinterpret_cast<const char*>(&prior->get_tau2_value((DataModel::ef_t)i)),
              sizeof(prior->get_tau2_value((DataModel::ef_t)i)));
        }
      }*/
    }

    // do rao
    if ((iter + 1) % n_rao == 0){
      if (!flat_proposal_dist || n_rao_burnin <= 0){
        data_model->sample_missing(current_model->model_inds, rng);
        samplerstats.rao_tic();
        raob->p_raoblackwell(this, current_model->y_hat, p_r_, p_r_types_);
        samplerstats.rao_toc();
      }

      if (n_rao_burnin <= 0){
        double z1 = (double)(p_rao_n + 1);
        double z2 = (double)p_rao_n / z1;
        for (size_t i = 0; i < m_g; ++i){
          p_rao[i] = z2 * p_rao[i] + p_r_[i] / z1;
        }

        if (use_types){
          for (size_t i = 0; i < m_g; ++i){
            for (size_t t = 0; t < n_types; ++t){
              p_rao_types[i][t] = z2 * p_rao_types[i][t] + p_r_types_[i][t] / z1;
            }
          }
        }

        ++p_rao_n;
      } else {
        n_rao_burnin--;
      }

      if (adaptation || n_rao_burnin > 0){
        if (!flat_proposal_dist){
          // compute new p_proposal
          double z1 = (double)(p_proposal_n + 1);
          double z2 = (double)p_proposal_n / z1;
          for (size_t i = 0; i < m_g; ++i){
            p_proposal[i] = z2 * p_proposal[i] + p_r_[i] / z1;
          }

          if (use_types){
            for (size_t i = 0; i < m_g; ++i){
              for (size_t t = 0; t < n_types; ++t){
                p_proposal_types[i][t] = z2 * p_proposal_types[i][t]
                                              + p_r_types_[i][t] / z1;
              }
            }
          }
        }

        // update move size distribution
        if (adapt_p_move_size){
          double t = std::max(1.0, k_move_size * (double)p_proposal_n);
          //unsigned char max_move_size2_ = (max_move_size_ + 1) / 2;
          if (acpt_move_size_goal_ > 0){
            adapt_p_move_size_acptrate(p_move_size_, n_acpt_moves_[0],
                      n_moves_[0], 2.0 + t);
          } else {
            adapt_p_move_size_jd_mb(p_move_size_, r_move_size_sum_, r_move_size_n_,
                                 max_move_size_);
          }
          Utils::geometric_dist_cdf(max_move_size_, p_move_size_, q_p_move_size_);
        }


        // update count
        ++p_proposal_n;

        if (!flat_proposal_dist){
          // update q_add_ and q_rem_ etc.
          compute_proposal_dist(p_proposal, q_add_min_, q_add_);
          for (size_t i = 0; i < m_g; ++i){
            q_rem_[i] = std::max(1 - p_proposal[i], q_rem_min_);
          }
          dd_add->update_weights(q_add_);
          dd_rem->update_weights(q_rem_);

          if (use_types) compute_proposal_types_dist(p_proposal_types,
              q_add_min_,
              q_add_types_,
              q_p_add_types_);
        }
      }
    } // end of rao block

    if (verbosity > 0 && (iter + 1) % verbosity == 0){
      *f_log << "(" << (iter + 1) << ")"
                 << " acc.rate " << (double)n_accepted / (iter + 1)
                 << " acc.rate2 " << (double)n_acpt_moves_[1] / n_moves_[1]
                 << " acc.rate3 " << (double)n_acpt_moves_[2] / n_moves_[2]
                 << " p_move_size " << p_move_size_
                 << " model size " << current_model->size()
                 << " log p " << current_model->log_likelihood_
                                 + prior->compute_log_model(current_model->Ns)
                 << " PVE " << pves[0] << "/" << pves[1] << "/" << pves[2]
                 << " sigma2 " << current_model->sigma2;
      /*for (int t = 0; t < 4; ++t){
        if (data_model->allow_terms((DataModel::ef_t)t)){
          *f_log << " a2*tau2" << DataModel::ef_name[t] << " "
                 << prior->get_tau2_value((DataModel::ef_t)t)
                        * prior->alpha() * prior->alpha();
        }
      }*/
      *f_log << std::endl;
      f_log->flush(); // write immediately
    } // end of verbosity block
  } // main loop (iter)
  n_iter = _n_iter;

  // record sampling time and other stats
  samplerstats.print(basename + "_samplerstats.txt", p_move_size_);
  // put sampling time also in sampling transcript
  *f_log << "Elapsed time: " << samplerstats.elapsed_time() << " seconds."
             << std::endl;

  delete dd_add;
  delete dd_rem;

  // write rao data to disk
  OPENBFILE(rao);
  f_rao->write(reinterpret_cast<const char*>(p_rao), m_g * sizeof(p_rao[0]));
  CLOSEFILE(rao);

  if (use_types){
    OPENBFILE(rao_types);
    for (size_t i = 0; i < m_g; ++i){
      f_rao_types->write(reinterpret_cast<const char*>(p_rao_types[i]),
                         n_types * sizeof(p_rao_types[i][0]));
    }
    CLOSEFILE(rao_types);
  }

  // close files
  CLOSEFILE(log);
  CLOSEFILE(loci);
  CLOSEFILE(modelsize);
  CLOSEFILE(jumpdistance);
  CLOSEFILE(log_likelihood);
  CLOSEFILE(log_prior);
  CLOSEFILE(move_type);
  CLOSEFILE(move_size);
  CLOSEFILE(pve);
  CLOSEFILE(alpha);
  //CLOSEFILE(tau2);
  /*for (int i = 0; i < 4; ++i){
    if (f_tau2[i] != NULL){
      CLOSEFILE(tau2[i]);
    }
  }*/
  CLOSEFILE(sigma2);
  if (use_types) CLOSEFILE(types);
  if (save_beta) CLOSEFILE(beta);
}

void compute_exhaustive_modelset(const size_t n_inds, ExhModel* exhmodel,
                                 double* log_model_probabilities,
                                 double& max_log_model)
{
  size_t nmodels = 0, tmp;
  size_t cmodel_binary = 0;
  int model_size = 0;

  size_t *cmodel_inds = new size_t[n_inds];

  for (size_t i = 0; i < n_inds; ++i)
    cmodel_inds[i] = i;

  log_model_probabilities[cmodel_binary] = exhmodel->log_prob();
  max_log_model = log_model_probabilities[cmodel_binary];

  for (size_t i = 0; i < n_inds; ++i){
    // comp full model with i+1 "new" variables
    ++model_size;
    if (i > 1){
      ++model_size;
      exhmodel->update_on_add();
      assert(std::isfinite(exhmodel->log_likelihood()));
    }
    cmodel_binary = std::pow(2, model_size) - 1;
    log_model_probabilities[cmodel_binary] = exhmodel->update_on_add();
    if (log_model_probabilities[cmodel_binary] > max_log_model)
      max_log_model = log_model_probabilities[cmodel_binary];
    assert(std::isfinite(exhmodel->log_likelihood()));
    assert(std::isfinite(log_model_probabilities[cmodel_binary]));

    // move new variable to left-most location
    for (size_t j = 0; j < i; ++j){
      --model_size;
      tmp = cmodel_inds[model_size - 1];
      cmodel_inds[model_size - 1] = cmodel_inds[model_size];
      cmodel_inds[model_size] = tmp;

      cmodel_binary &= ~(1 << tmp);
      log_model_probabilities[cmodel_binary] = exhmodel->update_on_moveleft();
      if (log_model_probabilities[cmodel_binary] > max_log_model)
        max_log_model = log_model_probabilities[cmodel_binary];
      assert(std::isfinite(exhmodel->log_likelihood()));
      assert(std::isfinite(log_model_probabilities[cmodel_binary]));
    }

    // compute rest of the models
    nmodels = std::pow(2, i) - i - 1;
    size_t j = 0;
    size_t K, nK = 0;
    char do_lefts = 0;
    while (j < nmodels){
      if (do_lefts < 2){
        // two new + swap
        ++model_size;
        tmp = cmodel_inds[model_size];
        cmodel_inds[model_size] = cmodel_inds[model_size - 1];
        cmodel_inds[model_size - 1] = tmp;

        cmodel_binary |= (1 << tmp);
        log_model_probabilities[cmodel_binary] = exhmodel->update_on_twonewswap();
        if (log_model_probabilities[cmodel_binary] > max_log_model)
          max_log_model = log_model_probabilities[cmodel_binary];
        assert(std::isfinite(exhmodel->log_likelihood()));
        assert(std::isfinite(log_model_probabilities[cmodel_binary]));

        ++j;
        ++do_lefts;
      } else {
        // K+1 swaps to left
        ++nK;
        K = 0;
        // values of K follow sequence 0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,...
        while (((nK >> K) & 1) == 0) K++;
        // do the swaps
        for (size_t k = 0; k <= K; ++k){
          --model_size;
          tmp = cmodel_inds[model_size - 1];
          cmodel_inds[model_size - 1] = cmodel_inds[model_size];
          cmodel_inds[model_size] = tmp;

          cmodel_binary &= ~(1 << tmp);
          log_model_probabilities[cmodel_binary] = exhmodel->update_on_moveleft();
          if (log_model_probabilities[cmodel_binary] > max_log_model)
            max_log_model = log_model_probabilities[cmodel_binary];
          assert(std::isfinite(exhmodel->log_likelihood()));
          assert(std::isfinite(log_model_probabilities[cmodel_binary]));

          ++j;
        }
        do_lefts = 0;
      }
    }

  }

  delete[] cmodel_inds;

}

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
                                             double* log_prop_probs)
{
  bool* isadd = new bool[n_inds];

  unsigned long nmodels = std::pow(2, n_inds);
  for (unsigned long i = 0; i < nmodels; ++i){
    size_t max_adds = m_g - const_loci;
    size_t max_rems = const_loci;
    size_t model_size = 0;
    double z_a = z_add; // for the "empty" model
    double z_r = z_rem;

    // find out the "normalized" add/rem pattern
    for (int j = 0; j < n_inds; ++j){
      unsigned char bit = (i >> j) & 1;
      unsigned char nind = bit_to_normalized_order[j];
      if (bit){
        // in model, move is removal
        isadd[nind] = false;
        z_a -= q_add[nind];
        z_r += q_rem[nind];
        ++max_rems; --max_adds;
        ++model_size;
      } else {
        // not in model, move is addition
        isadd[nind] = true;
      }
    }

    // compute the proposal probability (loops in the "normalized" order)
    int last_rem_pos = n_inds - 1;
    for (int j = 0; j < n_inds; ++j){
      if (max_adds > 0 && max_rems > 0) log_prop_probs[i] += LOG_HALF;

      if (isadd[j]){
        // adds are performed in order
        assert(max_adds > 0);
        --max_adds;
        log_prop_probs[i] += log(q_add[j]) - log(z_a);
        z_a -= q_add[j];

        if (use_types) log_prop_probs[i] += log_q_add_types[j];
      } else {
        // rems are performed in reverse order
        assert(max_rems > 0);
        --max_rems;
        while (isadd[last_rem_pos]) --last_rem_pos;

        log_prop_probs[i] += log(q_rem[last_rem_pos]) - log(z_r);
        z_r -= q_rem[last_rem_pos];
        --last_rem_pos;
      }
    }
    assert(std::isfinite(log_prop_probs[i]));
  }

  delete[] isadd;
}

} // namespace bmagwa
