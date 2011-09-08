/* BMAGWA software v1.0
 *
 * sampler.cpp
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


#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include "sampler.hpp"
#include "utils.hpp"

namespace bmagwa {

void RaoBlackwellizer::p_raoblackwell(const Sampler* sampler,
                                      const Vector& y_hat,
                                      double *p_r, double **p_r_types)
{
  const Model* cmodel = sampler->current_model;
  const DataModel* data_model = sampler->data_model;
  const Prior* prior = sampler->prior;

  const size_t n = data_model->n;
  //const double sum_log_Q = log(prior->get_tau2_value());
  const double log_n_2 = 0.5 * log(n);
  const double sigma2_times_2 = 2 * cmodel->sigma2;

  // set up types
  const size_t n_types = data_model->n_types;
  const bool use_types = (n_types > 1);

  // residual (r = y - y_hat)
  r_cm.set_to_vector_diff(data_model->y(), y_hat);
  double r2_cm = r_cm.mean();
  r2_cm = r2_cm * r2_cm * n;

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

  // main loop
  for (size_t i = 0; i < m_g; ++i){
    Vector* residual;
    double r2;
    double* lmp;
    double sum_log_Q;

    int32_t model_ind = cmodel->model_inds[i];

    // prepare residual and choose log prior change from tabulated values
    if (model_ind < 0){
      // not in model
      residual = &r_cm;
      r2 = r2_cm;
      lmp = lmp_g_add;
    } else {
      // in model

      unsigned int x_ind = cmodel->x_ind1[model_ind];
      DataModel::ef_t type = data_model->types[cmodel->type_inds[model_ind]];
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
        unsigned int x_ind = cmodel->x_ind2[model_ind];
        r_om.add_vector(cmodel->x.column(x_ind), cmodel->beta(x_ind));
      }

      residual = &r_om;
      r2 = r_om.mean();
      r2 = r2 * r2 * n;

      lmp = lmp_g_rem[type];
    }

    double max_types = -INFINITY;
    // compute rao-blackwellization for each type
    for (size_t t = 0; t < n_types; ++t){
      DataModel::ef_t type = data_model->types[t];
      size_t offset = pre_xxcov->xx_offset[type];

      double *xx_tmp = pre_xxcov->xx_types[type] + offset * offset * i;

      if (type == DataModel::AH){
        xxcov.resize(3);
        x.resize(n, 3);
        v.resize(3);

        // copy upper triangular part
        xxcov(0, 0) = xx_tmp[0]; // 1,1
        xxcov(0, 1) = xx_tmp[3]; // 1,2
        xxcov(1, 1) = xx_tmp[4]; // 2,2
        xxcov(0, 2) = xx_tmp[6]; // 1,3
        xxcov(1, 2) = xx_tmp[7]; // 2,3
        xxcov(2, 2) = xx_tmp[8]; // 3,3

        data_model->update_xx_ah(i, xxcov);

        xxcov(1, 1) += 1.0 / prior->get_tau2_value(DataModel::A);
        xxcov(2, 2) += 1.0 / prior->get_tau2_value(DataModel::H);

        sum_log_Q = log(prior->get_tau2_value(DataModel::A))
                    + log(prior->get_tau2_value(DataModel::H));

        data_model->get_genotypes_additive(i, x.column(1));
        data_model->get_genotypes_heterozygous(i, x.column(2));
      } else {
        xxcov.resize(2);
        x.resize(n, 2);
        v.resize(2);

        // copy upper triangular part
        xxcov(0, 0) = xx_tmp[0];          // 1,1
        xxcov(0, 1) = xx_tmp[offset];     // 1,2
        xxcov(1, 1) = xx_tmp[offset + 1]; // 2,2
        offset = 2;

        (data_model->*data_model->bmagwa::DataModel::update_xx[type])(i, xxcov);

        xxcov(1, 1) += 1.0 / prior->get_tau2_value(type);

        sum_log_Q = log(prior->get_tau2_value(type));

        (data_model->*data_model->bmagwa::DataModel::get_genotypes[type])(i,
                                                                   x.column(1));
      }

      if (!xxcov.cholesky()){
        // xx is not pos def.
        std::cout << "Warning, chol: xx is not positive definitive. rao QTL:"
                  << i << std::endl;
        if (use_types) p_r_types[i][t] = -INFINITY;
        else p_r[i] = 0;
        continue;
      }

      double sum_log_diag_l = 0;
      for (size_t j = 0; j < offset; ++j){
        sum_log_diag_l += log(xxcov(j, j));
      }

      // compute v = x' * r
      v.set_to_product(x, *residual, true);

      // v = L' \ v
      v.multiply_by_invtriangularmatrix(xxcov, true);

      double p = (VectorView::dotproduct(v, v) - r2) / sigma2_times_2
                 - sum_log_diag_l - 0.5 * sum_log_Q
                 + *(lmp + type) + log_n_2;

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
            p_r_types[i][t] = 1/n_types;
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
} // end of p_raoblackwell method


void Sampler::sample_missing()
{
  assert(current_model->size() == new_model->size());

  const size_t n = data_model->n;
  const double sigma2_times_2 = current_model->sigma2 * 2;
  double lprior[3], likelihood[3];

  // compute y_hat and residual
  Vector y_hat(n);
  y_hat.set_to_product(current_model->x, current_model->beta, false);

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
                   - current_model->beta(x_ind1) * ((double)k - tmp_x1(i_miss));
            likelihood[k] = lprior[k] - (r2 + new_term * new_term - old_term2) / sigma2_times_2;
          }
          break;
        case DataModel::H:
          for (int k = 0; k < 3; ++k){
            new_term = residual(i_miss) - current_model->beta(x_ind1) * ((double)(k == 1) - tmp_x1(i_miss));
            likelihood[k] = lprior[k] - (r2 + new_term * new_term - old_term2) / sigma2_times_2;
          }
          break;
        case DataModel::R:
          for (int k = 0; k < 3; ++k){
            new_term = residual(i_miss) - current_model->beta(x_ind1) * ((double)(k == 2) - tmp_x1(i_miss));
            likelihood[k] = lprior[k] - (r2 + new_term * new_term - old_term2) / sigma2_times_2;
          }
          break;
        case DataModel::D:
          for (int k = 0; k < 3; ++k){
            new_term = residual(i_miss) - current_model->beta(x_ind1) * ((double)(k > 0) - tmp_x1(i_miss));
            likelihood[k] = lprior[k] - (r2 + new_term * new_term - old_term2) / sigma2_times_2;
          }
          break;
        case DataModel::AH:
          for (int k = 0; k < 3; ++k){
            new_term = residual(i_miss) - current_model->beta(x_ind1) * ((double)k - tmp_x1(i_miss));
            new_term -= current_model->beta(x_ind2) * ((double)(k == 1) - tmp_x2(i_miss));
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
          y_hat(i_miss) += current_model->beta(x_ind1) * ((double)k - tmp_x1(i_miss));
          tmp_x1(i_miss) = (double)k;
          break;
        case DataModel::H:
          k = (k == 1);
          y_hat(i_miss) += current_model->beta(x_ind1) * ((double)k - tmp_x1(i_miss));
          tmp_x1(i_miss) = (double)k;
          break;
        case DataModel::R:
          k = (k == 2);
          y_hat(i_miss) += current_model->beta(x_ind1) * ((double)k - tmp_x1(i_miss));
          tmp_x1(i_miss) = (double)k;
          break;
        case DataModel::D:
          k = (k > 0);
          y_hat(i_miss) += current_model->beta(x_ind1) * ((double)k - tmp_x1(i_miss));
          tmp_x1(i_miss) = (double)k;
          break;
        case DataModel::AH:
          y_hat(i_miss) += current_model->beta(x_ind1) * ((double)k - tmp_x1(i_miss));
          tmp_x1(i_miss) = (double)k;
          k = (k == 1);
          y_hat(i_miss) += current_model->beta(x_ind2) * ((double)k - tmp_x2(i_miss));
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

  current_model->beta_sigma2_sampled = false;
} // sample_missing

void Sampler::compute_p_moves(const size_t n_loci,
                              double* p_moves, double* p_moves_cumsum) const
{
  if (n_loci > 0){
    if (use_types){
      p_moves[0] = 0.4;
      p_moves[1] = 0.4;
      p_moves[2] = 0.1;
      p_moves[3] = 0.1;
      p_moves_cumsum[0] = 0.4;
      p_moves_cumsum[1] = 0.8;
      p_moves_cumsum[2] = 0.9;
      p_moves_cumsum[3] = 1;
    } else {
      p_moves[0] = 0.4;
      p_moves[1] = 0.4;
      p_moves[2] = 0.2;
      p_moves[3] = 0;
      p_moves_cumsum[0] = 0.4;
      p_moves_cumsum[1] = 0.8;
      p_moves_cumsum[2] = 1;
      p_moves_cumsum[3] = 1;
    }

  } else {
    p_moves[0] = 1;
    p_moves[1] = 0;
    p_moves[2] = 0;
    p_moves[3] = 0;
    p_moves_cumsum[0] = 1;
    p_moves_cumsum[1] = 1;
    p_moves_cumsum[2] = 1;
    p_moves_cumsum[3] = 1;
  }
}

void Sampler::compute_proposal_dist_t(const Model* model,
                                      const double* p_proposal,
                                      double* q_add, double* q_p_add) const
{
  double val = (model->model_inds[0] < 0) ? pow(p_proposal[0], t_proposal) : 0;

  q_add[0] = val;
  q_p_add[0] = val;
  for (size_t i = 1; i < m_g; ++i){
    val = (model->model_inds[i] < 0) ? pow(p_proposal[i], t_proposal) : 0;
    q_add[i] = val;
    q_p_add[i] = q_p_add[i - 1] + val;
  }
}

void Sampler::remdate_proposal_dist_t(const size_t ind, double p_proposal,
                                      double* q_add, double* q_p_add) const
{
  p_proposal = pow(p_proposal, t_proposal);
  q_add[ind] = p_proposal;
  for (size_t i = ind; i < m_g; ++i){
    q_p_add[i] += p_proposal;
  }
}


void Sampler::adddate_proposal_dist_t(const size_t ind, double p_proposal,
                                      double* q_add, double* q_p_add) const
{
  p_proposal = pow(p_proposal, t_proposal);
  q_add[ind] = 0;
  for (size_t i = ind; i < m_g; ++i){
    q_p_add[i] -= p_proposal;
  }
}


void Sampler::switchdate_proposal_dist_t(const size_t ind_add,
                                         double p_proposal_add,
                                         const size_t ind_rem,
                                         double p_proposal_rem,
                                         double* q_add, double* q_p_add) const
{
  p_proposal_add = pow(p_proposal_add, t_proposal);
  p_proposal_rem = pow(p_proposal_rem, t_proposal);

  q_add[ind_add] = 0;
  q_add[ind_rem] = p_proposal_rem;

  double diff = (p_proposal_rem - p_proposal_add);

  if (ind_add < ind_rem){
    size_t i = ind_add;
    for (; i < ind_rem; ++i){
      q_p_add[i] -= p_proposal_add;
    }
    for (; i < m_g; ++i){
      q_p_add[i] += diff;
    }
  } else {
    size_t i = ind_rem;
    for (; i < ind_add; ++i){
      q_p_add[i] += p_proposal_rem;
    }
    for (; i < m_g; ++i){
      q_p_add[i] += diff;
    }
  }
}

void Sampler::compute_proposal_types_dist_t(const double*const* p_proposal_types,
                                            double** q_add_types,
                                            double** q_p_add_types) const
{
  double val;

  for (size_t i = 0; i < m_g; ++i){
    val = pow(p_proposal_types[i][0], t_proposal);
    q_add_types[i][0] = val;
    q_p_add_types[i][0] = val;

    for (size_t t = 1; t < n_types; ++t){
      val = pow(p_proposal_types[i][t], t_proposal);
      q_add_types[i][t] = val;
      q_p_add_types[i][t] = q_p_add_types[i][t - 1] + val;
    }
  }
}

size_t Sampler::sample_discrete(const double* cumsum, const size_t m){
  size_t a = 0;
  size_t b = m - 1;
  size_t c;
  int level_ = level;

  double z = cumsum[b];
  double r = rng.rand_01() * z;

  while (level_ > 0){
    c = (b - a) / 2;
    if (r <= cumsum[c]){
      b = c;
    } else {
      a = c + 1;
    }
    level_--;
  }

  for (size_t i = a; i <= b; ++i){
    if (r <= cumsum[i]){
      return i;
    }
  }

  throw std::logic_error("sample_discrete reached end");
}


#define OPENFILE(item) std::ofstream* f_ ## item = Utils::openfile(basename + "_" + #item + ".dat", false)
#define OPENBFILE(item) std::ofstream* f_ ## item = Utils::openfile(basename + "_" + #item + ".dat", true)
#define CLOSEFILE(item) {(f_ ## item)->close(); delete f_ ## item;}

void Sampler::sample()
{
  if (!std::isfinite(current_model->log_likelihood_))
    throw std::runtime_error("Sampler cannot start from non-finite likelihood.");

  time_t t_elapsed = time(NULL);
  double pves[3]; pves[0] = 0.0; pves[1] = 0.0; pves[2] = 0.0;

  // open files
  OPENFILE(verbose);
  *f_verbose << std::setprecision(3) << std::fixed;
  OPENBFILE(loci);
  OPENBFILE(modelsize);
  OPENBFILE(accepted);
  OPENBFILE(log_likelihood);
  OPENBFILE(log_prior);
  OPENBFILE(move_type);
  OPENBFILE(pve);
  std::ofstream* f_tau2[] = {NULL, NULL, NULL, NULL};
  for (int i = 0; i < 4; ++i){
    if (data_model->allow_terms((DataModel::ef_t)i)){
      f_tau2[i] =  Utils::openfile(basename + "_tau2_" + DataModel::ef_name[i]
                                                                + ".dat", true);
    }
  }
  //OPENBFILE(tau2);
  OPENBFILE(sigma2);
  std::ofstream* f_types = NULL;
  if (use_types) f_types = Utils::openfile(basename + "_types.dat", true);
  std::ofstream* f_beta = NULL;
  if (save_beta) f_beta = Utils::openfile(basename + "_beta.dat", true);

  *f_verbose << "Sampler name is " << basename << " and seed is " << seed
             << std::endl;

  // models need to be initialized to same
  *new_model = *current_model;

  // compute proposal dist
  compute_proposal_dist_t(current_model, p_proposal, q_add_, q_p_add_);
  memcpy(new_q_add_, q_add_, m_g * sizeof(double));
  memcpy(new_q_p_add_, q_p_add_, m_g * sizeof(double));
  if (use_types)
    compute_proposal_types_dist_t(p_proposal_types, q_add_types_, q_p_add_types_);

  compute_p_moves(current_model->size(), p_moves_, p_moves_cumsum_);

  // sample beta and sigma2 so that tau2 can be sampled
  current_model->sample_beta_sigma2(rng);
  sample_missing(); // this is not really needed here?
  prior->sample_tau2(current_model->x_types, current_model->beta,
                     current_model->sigma2, rng);
  // compute likelihood after missing and tau2 sampling
  current_model->compute_log_likelihood();
  // update new_model to be equal to current
  *new_model = *current_model;

  size_t _n_iter = n_iter + do_n_iter;
  for (size_t iter = n_iter; iter < _n_iter; ++iter){
    // sample missing genotypes (requires current and new models to be the same)
    if ((iter + 1) % n_sample_tau2_and_missing == 0){
      sample_missing();
      prior->sample_tau2(current_model->x_types, current_model->beta,
                         current_model->sigma2, rng);
      // recompute likelihood and update new model to current
      current_model->compute_log_likelihood();
      *new_model = *current_model;
    }

    // select move
    assert(current_model->size() == 0 || p_moves_cumsum_[0] == 0.4);
    assert(current_model->size() == 0 || p_moves_cumsum_[1] == 0.8);
    assert(current_model->size() == 0 || p_moves_cumsum_[2] == 1 || p_moves_cumsum_[2] == 0.9);

    char move = (char)Utils::sample_discrete_naive(p_moves_cumsum_, 4, rng);
    char accept = 0;

    assert(std::isfinite(current_model->log_likelihood_));
    switch (move){
      case 0: // addition
        accept = do_addition();
        break;
      case 1: // removal
        accept = do_removal();
        break;
      case 2: // SNP switch
        accept = do_switch_of_snps();
        break;
      case 3: // type switch
        accept = do_switch_of_types();
        break;
      default:
        throw std::logic_error("Move is not in 0...3");
    }
    assert(std::isfinite(current_model->log_likelihood_));

    if ((iter + 1) % thin == 0 || (iter + 2) % n_sample_tau2_and_missing == 0){
      // note: needs to be done before Model_pve and before sampling tau2
      current_model->sample_beta_sigma2(rng);
    }

    // collect sample
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
      f_accepted->write(reinterpret_cast<const char*>(&accept), sizeof(accept));
      f_log_likelihood->write(
                  reinterpret_cast<const char*>(&current_model->log_likelihood_),
                  sizeof(current_model->log_likelihood_));
      double log_prior = prior->compute_log_model(current_model->Ns);
      f_log_prior->write(
                  reinterpret_cast<const char*>(&log_prior), sizeof(log_prior));
      f_move_type->write(reinterpret_cast<const char*>(&move), sizeof(move));
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

      for (int i = 0; i < 4; ++i){
        if (data_model->allow_terms((DataModel::ef_t)i)){
          (f_tau2[i])->write(reinterpret_cast<const char*>(&prior->get_tau2_value((DataModel::ef_t)i)),
              sizeof(prior->get_tau2_value((DataModel::ef_t)i)));
        }
      }
    }

    // do rao
    if ((iter + 1) % n_rao == 0){
      data_model->sample_missing(current_model->model_inds, rng);
      raob->p_raoblackwell(this, current_model->y_hat, p_r_, p_r_types_);

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

        ++p_proposal_n;

        compute_proposal_dist_t(current_model, p_proposal, q_add_, q_p_add_);
        memcpy(new_q_add_, q_add_, m_g * sizeof(double));
        memcpy(new_q_p_add_, q_p_add_, m_g * sizeof(double));

        if (use_types) compute_proposal_types_dist_t(p_proposal_types,
                                                     q_add_types_,
                                                     q_p_add_types_);
      }
    } // end of rao block

    if (verbosity > 0 && (iter + 1) % verbosity == 0){
      *f_verbose << "(" << (iter + 1) << ")"
                 << " acc.rate " << (double)n_accepted / (iter + 1)
                 << " model size " << current_model->size()
                 << " log p " << current_model->log_likelihood_
                                 + prior->compute_log_model(current_model->Ns)
                 << " PVE " << pves[0] << "/" << pves[1] << "/" << pves[2]
                 << " sigma2 " << current_model->sigma2;
      for (int t = 0; t < 4; ++t){
        if (data_model->allow_terms((DataModel::ef_t)t)){
          *f_verbose << " tau2" << DataModel::ef_name[t] << " "
                     << prior->get_tau2_value((DataModel::ef_t)t);
        }
      }
      *f_verbose << std::endl;
      f_verbose->flush(); // write immediately
    } // end of verbosity block
  } // main loop (iter)
  n_iter = _n_iter;

  // write rao data to disk
  OPENFILE(rao);
  f_rao->write(reinterpret_cast<const char*>(p_rao), m_g * sizeof(p_rao[0]));
  CLOSEFILE(rao);

  if (use_types){
    OPENFILE(rao_types);
    for (size_t i = 0; i < m_g; ++i){
      f_rao_types->write(reinterpret_cast<const char*>(p_rao_types[i]),
                         n_types * sizeof(p_rao_types[i][0]));
    }
    CLOSEFILE(rao_types);
  }

  // get elapsed time
  *f_verbose << "Elapsed time: " << (time(NULL) - t_elapsed) << " seconds."
             << std::endl;

  // close files
  CLOSEFILE(verbose);
  CLOSEFILE(loci);
  CLOSEFILE(modelsize);
  CLOSEFILE(accepted);
  CLOSEFILE(log_likelihood);
  CLOSEFILE(log_prior);
  CLOSEFILE(move_type);
  CLOSEFILE(pve);
  //CLOSEFILE(tau2);
  for (int i = 0; i < 4; ++i){
    if (f_tau2[i] != NULL){
      CLOSEFILE(tau2[i]);
    }
  }
  CLOSEFILE(sigma2);
  if (use_types) CLOSEFILE(types);
  if (save_beta) CLOSEFILE(beta);
}

} // namespace bmagwa
