/* BMAGWA software v2.0
 *
 * model.hpp
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


#ifndef MODEL_HPP_
#define MODEL_HPP_

#ifndef MAXMODELSIZE
#define MAXMODELSIZE 1024
#endif

#include <vector>
#include <stdexcept>
#include <stdlib.h>
#include <iostream>
#include "linalg.hpp"
#include "data.hpp"
#include "data_model.hpp"
#include "prior.hpp"
#include "samplerstats.hpp"

namespace bmagwa {

// forward declaration
class RaoBlackwellizer;

//! Holds the state of the linear model and implements likelihood computation etc.
class Model
{
  public:
    Model(const Data* data, const DataModel* data_model, const Prior* prior,
          SamplerStats* samplerstats)
    : data_(data), data_model_(data_model), prior_(prior),
      samplerstats_(samplerstats),
      n_model_inds(data->m_g + data->m_e),
      x(data->n, MAXMODELSIZE, data->n, data->m_e),
      xx(MAXMODELSIZE, data->m_e),
      l(MAXMODELSIZE, data->m_e),
      xy(MAXMODELSIZE, data->m_e),
      v(MAXMODELSIZE, data->m_e),
      inv_tau2_alpha2(MAXMODELSIZE, data->m_e),
      y_hat(data->n),
      y_hat_tmp_e(data->n),
      c_tmp(NULL),
      s_tmp(NULL),
      beta(MAXMODELSIZE, data->m_e),
      mu_beta(MAXMODELSIZE, data->m_e),
      sigma2(0),
      syx_plus_vs2(0),
      log_det_invQ(0),
      log_det_invQ_plus_xx(0),
      mu_beta_computed(false),
      log_likelihood_(NAN),
      loci(0),
      x_ind_to_g_ind(data->m_e, 0),
      model_inds(n_model_inds, -1),
      type_inds(0),
      x_types(0),
      x_ind1(0),
      x_ind2(0)
    {
      for (size_t i = 0; i < 5; i++){
        Ns[i] = 0;
      }

      // reserve memory to avoid reallocs
      loci.reserve(MAXMODELSIZE);
      x_ind_to_g_ind.reserve(MAXMODELSIZE);
      type_inds.reserve(MAXMODELSIZE);
      x_ind1.reserve(MAXMODELSIZE);
      x_ind2.reserve(MAXMODELSIZE);
      x_types.reserve(MAXMODELSIZE);

      // add E terms to model
      if (data_->m_e <= MAXMODELSIZE){
        x = data_->e(); // copy
        xx.set_to_innerproduct(x);   // x' * x
        xy.set_to_product(x, data->y(), true); // x' * y

        // use AH type for fixed covariates. AH should never be in x_types
        // elsewhere (since AH splits to A and H)
        for (size_t i = 0; i < data_->m_e; ++i)
          x_types.push_back(DataModel::AH);

        prior->set_inv_tau2_e(inv_tau2_alpha2);

        compute_log_likelihood();
      } else {
        throw std::runtime_error(
                      "Not enough memory for initializing environmental data.");
      }

      c_tmp = new double[MAXMODELSIZE];
      s_tmp = new double[MAXMODELSIZE];
    }

    Model& operator=(const Model& src)
    {
      if (this != &src){
        if (n_model_inds != src.n_model_inds)
          throw std::logic_error(
              "Models must have same number of model_inds on copy");

        data_ = src.data_;
        data_model_ = src.data_model_;
        prior_ = src.prior_;

        x.resize(src.x.rows(), src.x.cols());
        x = src.x;
        xx.resize(src.xx.length());
        xx = src.xx;
        l.resize(src.l.length());
        l = src.l;
        xy.resize(src.xy.length());
        xy = src.xy;
        v.resize(src.v.length());
        v = src.v;
        inv_tau2_alpha2.resize(src.inv_tau2_alpha2.length());
        inv_tau2_alpha2 = src.inv_tau2_alpha2;
        //c = src.c;
        //s = src.s;
        beta.resize(src.beta.length());
        beta = src.beta;
        mu_beta.resize(src.mu_beta.length());
        mu_beta = src.mu_beta;
        sigma2 = src.sigma2;
        syx_plus_vs2 = src.syx_plus_vs2;
        log_det_invQ = src.log_det_invQ;
        log_det_invQ_plus_xx = src.log_det_invQ_plus_xx;
        memcpy(Ns, src.Ns, 5 * sizeof(int));
        mu_beta_computed = src.mu_beta_computed;
        log_likelihood_ = src.log_likelihood_;

        loci = src.loci;
        x_ind_to_g_ind = src.x_ind_to_g_ind;
        model_inds = src.model_inds;
        type_inds = src.type_inds;
        x_types = src.x_types;
        x_ind1 = src.x_ind1;
        x_ind2 = src.x_ind2;
      }

      return *this;
    }

    virtual ~Model() {
      delete[] s_tmp;
      delete[] c_tmp;
    }

    const size_t size() const {
      return loci.size();
    }
    const double log_likelihood() const {
      return log_likelihood_;
    }
    const double get_sigma2_value() const {
        return sigma2;
    }
    const VectorView get_beta() const {
      return beta;
    }
    // for testing purposes mostly:
    const int* get_Ns() const {
      return Ns;
    }

    void set_inv_tau2_alpha2_value(const unsigned int x_ind, double value)
    {
      if (x_ind < inv_tau2_alpha2.length()){
        inv_tau2_alpha2(x_ind) = value;
      }
    }

    void compute_log_likelihood()
    {
      if (x.cols() < 1){
        log_likelihood_ = NAN;
        return;
      }

      samplerstats_->likelihood_computation();

      log_det_invQ = prepare_l_chol_and_get_det_invQ();

      if (!l.cholesky()){
        /*std::cout << "Warning, chol: xx is not positive definitive. QTL:";
        for (size_t i = 0; i < loci.size(); ++i){
          std::cout << " " << loci[i];
        }
        std::cout << std::endl;*/
        log_likelihood_ = -INFINITY;
        syx_plus_vs2 = INFINITY;
        mu_beta_computed = false;
        return;
      }

      // update v to x' * y and compute its real value (i.e. l' \ (x' * y))
      v = xy;
      v.multiply_by_invtriangularmatrix(l, true);

      // v s2 + y y - v' * v
      syx_plus_vs2 = prior_->nus2_plus_yy - Vector::dotproduct(v, v);

      log_det_invQ_plus_xx = 0;
      for (size_t i = 0; i < l.length(); ++i){
        log_det_invQ_plus_xx += log(l(i, i));
      }

      log_likelihood_ = log_det_invQ - log_det_invQ_plus_xx
                       + prior_->minus_n_plus_nu_2 * log(syx_plus_vs2);
      mu_beta_computed = false;
    }

    void add_term(const uint32_t ind, const char t_ind, const double* inv_tau2_alpha2_val,
                  const bool update_likelihood = true)
    {
      DataModel::ef_t type = data_model_->types[t_ind];
      bool type_is_AH = (type == DataModel::AH);

      size_t model_ind = loci.size();

      if (x.cols() + 1 + (size_t)type_is_AH <= MAXMODELSIZE){
        loci.push_back(ind);
        ++Ns[type];
        type_inds.push_back(t_ind);
        model_inds[ind] = model_ind;
        mu_beta_computed = false;

        if (type_is_AH){
          update_likelihood_on_add(ind, model_ind, DataModel::A, false, inv_tau2_alpha2_val[0], update_likelihood);
          if (std::isfinite(syx_plus_vs2))
            update_likelihood_on_add(ind, model_ind, DataModel::H, true, inv_tau2_alpha2_val[1], update_likelihood);
        } else {
          update_likelihood_on_add(ind, model_ind, type, false, inv_tau2_alpha2_val[0], update_likelihood);
        }
      } else {
        throw std::runtime_error("Not enough memory for adding term.");
      }

      assert(x.cols() == v.length());
      assert(x.cols() == xy.length());
      assert(x.cols() == l.length());
    }

    void remove_term(const int model_ind, const bool update_likelihood = true)
    {
      uint32_t ind = loci[model_ind];
      DataModel::ef_t type = data_model_->types[type_inds[model_ind]];
      unsigned int xi1 = x_ind1[model_ind];
      unsigned int xi2 = x_ind2[model_ind];

      model_inds[ind] = -1;
      --Ns[type];

      loci.erase(loci.begin() + model_ind);
      type_inds.erase(type_inds.begin() + model_ind);
      x_ind1.erase(x_ind1.begin() + model_ind);
      x_ind2.erase(x_ind2.begin() + model_ind);

      int xmove = (type == DataModel::AH) ? 2 : 1;
      for (size_t i = model_ind; i < loci.size(); ++i){
        // assumes that AH were added to x consecutively!
        // and that if model_ind > model_ind_rem, then x_ind > x_ind_rem!
        x_ind1[i] -= xmove;
        x_ind2[i] -= xmove;
        model_inds[loci[i]]--;
      }

      if (type == DataModel::AH){
        // remove the one with larger index first
        if (xi1 > xi2){
          update_likelihood_on_remove(xi1, update_likelihood);
          update_likelihood_on_remove(xi2, update_likelihood);
        } else {
          update_likelihood_on_remove(xi2, update_likelihood);
          update_likelihood_on_remove(xi1, update_likelihood);
        }
      } else {
        update_likelihood_on_remove(xi1, update_likelihood);
      }

      mu_beta_computed = false;

      assert(x.cols() == v.length());
      assert(x.cols() == xy.length());
      assert(x.cols() == l.length());
    }

    void compute_mu_beta()
    {
      if (!mu_beta_computed){
        // update mu
        mu_beta.resize(v.length());
        mu_beta = v;
        mu_beta.multiply_by_invtriangularmatrix(l, false);

        mu_beta_computed = true;
      }
    }

    void sample_beta_sigma2(Rand& rng)
    {
      // sample beta and sigma2
      sigma2 = rng.rand_sinvchi2(syx_plus_vs2 / prior_->n_plus_nu);

      assert(v.length() == x.cols());

      // find mean of beta
      compute_mu_beta();

      // sample beta
      beta.resize(mu_beta.length());
      rng.rand_mvnormal_invchol(mu_beta, l, sqrt(sigma2), beta);

      assert(beta.length() == x.cols());
      assert(std::isfinite(beta(0)));
    }

    // need to leave y_hat to be correct
    void compute_pve(double *pves)
    {
      size_t m_e = data_model_->m_e;

      bool compute_e = (m_e > 1);
      if (compute_e){
        y_hat_tmp_e.set_to_product(x.columnblock(0, m_e), beta.block(0, m_e),
                                   false);
        pves[2] = y_hat_tmp_e.var();
      } else {
        pves[2] = 0.0;
        y_hat_tmp_e = beta(0); // constant term
      }

      bool compute_g = (loci.size() > 0);
      if (compute_g){
        y_hat.set_to_product(x.columnblock(m_e, x.cols() - m_e),
                             beta.block(m_e, x.cols() - m_e),
                             false);
        pves[1] = y_hat.var();
      } else {
        pves[1] = 0.0;
      }

      if (compute_e || compute_g){
        for (size_t i = 0; i < y_hat.length(); ++i){
          y_hat(i) += y_hat_tmp_e(i);
        }
        if (!compute_e){ // if only g is active
          pves[0] = pves[1];
        } else if (!compute_g){ // if only e is active
          pves[0] = pves[2];
        } else { // if both are active
          pves[0] = y_hat.var();
        }
      } else {
        y_hat = y_hat_tmp_e; // only constant term
        pves[0] = 0;
        pves[1] = 0;
        pves[2] = 0;
        return;
      }

      double z = (pves[0] + sigma2);
      pves[0] /= z;
      pves[1] /= z;
      pves[2] /= z;
    }

  protected:
    const Data* data_;
    const DataModel* data_model_;
    const Prior* prior_;
    SamplerStats* samplerstats_;

    //size_t n;
    const uint32_t n_model_inds;

    Matrix x;
    SymmMatrix xx, l;
    Vector xy, v;
    Vector inv_tau2_alpha2;
    Vector y_hat, y_hat_tmp_e; // pre-allocated tmps for pve computations
    double* c_tmp;
    double* s_tmp;
    Vector beta, mu_beta;
    double sigma2;
    double syx_plus_vs2, log_det_invQ, log_det_invQ_plus_xx;
    int Ns[5];
    bool mu_beta_computed;
    double log_likelihood_;

    std::vector<uint32_t> loci;          // length nloci
    std::vector<uint32_t> x_ind_to_g_ind;// length x.length()
    std::vector<int32_t> model_inds;     // length m_g
    std::vector<unsigned char> type_inds;         // length nloci
    std::vector<DataModel::ef_t> x_types; // length x.length()
    std::vector<unsigned int> x_ind1;    // length nloci
    std::vector<unsigned int> x_ind2;    // length nloci

    Model(const Model& m)
    : data_(m.data_), data_model_(m.data_model_), prior_(m.prior_),
      n_model_inds(0), x(0, 0), xx(0), l(0), xy(0), v(0), inv_tau2_alpha2(0),
      y_hat(0), y_hat_tmp_e(0), c_tmp(NULL),
      s_tmp(NULL), beta(0), mu_beta(0)
    {}

    void update_likelihood_on_add(const uint32_t ind,
                                  const size_t model_ind,
                                  const DataModel::ef_t type,
                                  const bool is_x_ind2,
                                  const double inv_tau2_alpha2_val,
                                  const bool update_likelihood)
    {
      if (type == DataModel::AH)
        throw std::logic_error("Do not call model update with type AH."
                              "Instead do update twice, first with A, then H.");

      size_t col = x.cols();
      if (is_x_ind2){
        x_ind2[model_ind] = col;
      } else {
        x_ind1.push_back(col);
        x_ind2.push_back(0); // don't care or set later to correct
      }
      x_ind_to_g_ind.push_back(ind);
      x_types.push_back(type);

      // get genotype data
      x.resize(x.rows(), col + 1);
      VectorView xcol(x.column(col));
      (data_model_->*data_model_->bmagwa::DataModel::get_genotypes[type])
                                                                    (ind, xcol);

      // update xy
      xy.resize(col + 1);
      xy(col) = Vector::dotproduct(xcol, data_->y());

      // update inv_tau2
      inv_tau2_alpha2.resize(col + 1);
      inv_tau2_alpha2(col) = inv_tau2_alpha2_val;

      // update xx
      xx.resize(col + 1);
      VectorView xxcol(xx.column(col));
      xxcol.set_to_product(x, xcol, true);

      if (update_likelihood){
        samplerstats_->likelihood_update_on_add();

        // update cholesky
        if (!l.cholesky_update(xxcol, inv_tau2_alpha2_val)){
          //std::cout << "Warning, chol update: xx is not positive definitive."
          //          << " QTL: " << loci[model_ind] << std::endl;
          log_likelihood_ = -INFINITY;
          syx_plus_vs2 = INFINITY;
          mu_beta_computed = false;
          v.resize(col + 1);
          return;
        }

        // update log posterior
        VectorView lcol(l.column(col));
        lcol.resize(col); // drop last
        double vm2 = xy(col) - Vector::dotproduct(lcol, v);
        lcol.resize(col + 1);
        v.resize(col + 1);
        v(col) = vm2 / lcol(col); // this is the sqrt of line below
        vm2 = v(col) * v(col); //vm2 = vm2 * vm2 / lmm2;

        syx_plus_vs2 = syx_plus_vs2 - vm2;

        log_det_invQ += 0.5 * log(inv_tau2_alpha2_val);
        log_det_invQ_plus_xx += log(lcol(col));

        log_likelihood_ = log_det_invQ - log_det_invQ_plus_xx
            + prior_->minus_n_plus_nu_2 * log(syx_plus_vs2);
      } else {
        l.resize(col + 1);
        v.resize(col + 1);
      }
    }

    void update_likelihood_on_remove(const unsigned int x_ind,
                                     const bool update_likelihood)
    {
      x.remove_column(x_ind);
      x_ind_to_g_ind.erase(x_ind_to_g_ind.begin() + x_ind);
      xy.remove(x_ind);
      xx.remove_colrow(x_ind);
      x_types.erase(x_types.begin() + x_ind);
      double inv_val = inv_tau2_alpha2(x_ind);
      inv_tau2_alpha2.remove(x_ind);

      if (update_likelihood){
        samplerstats_->likelihood_update_on_rem();

        if (x_ind == x.cols()){
          // last variable, i.e. chol downdate == drop last col/row
          log_det_invQ_plus_xx -= log(l(x_ind, x_ind));
          syx_plus_vs2 += v(x_ind) * v(x_ind);


          l.resize(xx.length());
          v.resize(xy.length());
        } else {
          // update chol
          l.cholesky_downdate(x_ind, c_tmp, s_tmp);

          // update v
          v.resize(xy.length());
          v = xy;

          v.multiply_by_invtriangularmatrix(l, true);

          log_det_invQ_plus_xx = 0;
          for (size_t i = 0; i < l.length(); i++){
            log_det_invQ_plus_xx += log(l(i, i));
          }

          syx_plus_vs2 = prior_->nus2_plus_yy - VectorView::dotproduct(v, v);
        }

        log_det_invQ -= 0.5 * log(inv_val);

        log_likelihood_ = log_det_invQ - log_det_invQ_plus_xx
            + prior_->minus_n_plus_nu_2 * log(syx_plus_vs2);
      } else {
        l.resize(xx.length());
        v.resize(xy.length());
      }
    }

    double prepare_l_chol_and_get_det_invQ()
    {
      // prepare l for cholesky
      assert(l.length() == xx.length());
      l = xx;

      // do invQ
      double log_sum_2 = 0;

      // handle E & G
      // assumption: first term is constant
      l(0, 0) += inv_tau2_alpha2(0); // may be zero, do not take log
      for (size_t i = 1; i < l.length(); ++i){
        l(i, i) += inv_tau2_alpha2(i);
        log_sum_2 += log(inv_tau2_alpha2(i));
      }

      return 0.5 * log_sum_2;
    }

    friend class RaoBlackwellizer;
    friend class Sampler;
    friend class Prior; // should use accessors...
    friend class ExhModel;
};

//! Utility class for computing an exhaustive set of models.
class ExhModel : public Model
{
  public:
    ExhModel(const Data* data, const DataModel* data_model, const Prior* prior,
             SamplerStats* samplerstats)
    : Model(data, data_model, prior, samplerstats),
      model_size_(0),
      const_loci_(0),
      v_size_(0),
      log_model_prior_(0.0)
    {}

    double log_prob() const
    {
      return log_likelihood_ + log_model_prior_;
    }

    double update_to_model(const Model& src, const size_t const_loci)
    {
      // reduced version of model copy
      {
        if (n_model_inds != src.n_model_inds)
          throw std::logic_error(
              "Models must have same number of model_inds on copy");

        data_ = src.data_;
        data_model_ = src.data_model_;
        prior_ = src.prior_;
        samplerstats_ = src.samplerstats_;

        l.resize(src.l.length());
        l = src.l;
        xy.resize(src.xy.length());
        xy = src.xy;
        v.resize(src.v.length());
        v = src.v;
        inv_tau2_alpha2.resize(src.inv_tau2_alpha2.length());
        inv_tau2_alpha2 = src.inv_tau2_alpha2;

        type_inds = src.type_inds;
        x_types = src.x_types;
        x_ind1 = src.x_ind1;
        x_ind2 = src.x_ind2;
      }

      const_loci_ = const_loci;
      model_size_ = const_loci_;
      if (model_size_ > 0){
        if (DataModel::AH == data_model_->types[type_inds[model_size_ - 1]])
          v_size_ = x_ind2[model_size_ - 1] + 1;
        else
          v_size_ = x_ind1[model_size_ - 1] + 1;
      } else {
        v_size_ = data_model_->m_e;
      }

      // update log likelihood
      VectorView vpart(v, 0, v_size_);
      syx_plus_vs2 = prior_->nus2_plus_yy - VectorView::dotproduct(vpart, vpart);
      assert(syx_plus_vs2 > 0);

      log_det_invQ_plus_xx = 0;
      for (size_t i = 0; i < v_size_; i++){
        log_det_invQ_plus_xx += log(l(i, i));
      }

      log_det_invQ = prior_->log_det_invQ_e();

      for (size_t i = data_model_->m_e; i < v_size_; ++i){
        log_det_invQ += log(inv_tau2_alpha2(i));
      }
      log_det_invQ *= 0.5;

      log_likelihood_ = log_det_invQ - log_det_invQ_plus_xx
                  + prior_->minus_n_plus_nu_2 * log(syx_plus_vs2);

      // update model prior
      log_model_prior_ = 0.0;

      for (size_t i = 0; i < 5; ++i){
        Ns[i] = 0;
      }
      for (size_t i = 0; i < model_size_; ++i)
        ++Ns[data_model_->types[type_inds[i]]];

      return log_likelihood_ + log_model_prior_;
    }

    double update_on_add()
    {
      DataModel::ef_t type = data_model_->types[type_inds[model_size_]];

      syx_plus_vs2 -= v(v_size_) * v(v_size_);
      log_det_invQ += 0.5 * log(inv_tau2_alpha2(v_size_));
      log_det_invQ_plus_xx += log(l(v_size_, v_size_));
      ++v_size_;

      if (DataModel::AH == type){
        syx_plus_vs2 -= v(v_size_) * v(v_size_);
        log_det_invQ += 0.5 * log(inv_tau2_alpha2(v_size_));
        log_det_invQ_plus_xx += log(l(v_size_, v_size_));
        ++v_size_;
      }
      assert(syx_plus_vs2 > 0);

      // update log_likelihood
      log_likelihood_ = log_det_invQ - log_det_invQ_plus_xx
                  + prior_->minus_n_plus_nu_2 * log(syx_plus_vs2);

      // update model prior
      log_model_prior_ += prior_->compute_log_change_on_add(Ns, model_size_, type);

      ++Ns[type];
      ++model_size_;

      return log_likelihood_ + log_model_prior_;
    }

    double update_on_moveleft()
    {
      --model_size_;
      size_t ind_keep = model_size_;
      size_t ind_rem = model_size_ - 1;
      DataModel::ef_t type_keep = data_model_->types[type_inds[ind_keep]];
      DataModel::ef_t type_rem = data_model_->types[type_inds[ind_rem ]];

      // remove old values of keep var
      --v_size_;
      syx_plus_vs2 += v(v_size_) * v(v_size_);
      log_det_invQ_plus_xx -= log(l(v_size_, v_size_));

      if (DataModel::AH == type_keep){
        --v_size_;
        syx_plus_vs2 += v(v_size_) * v(v_size_);
        log_det_invQ_plus_xx -= log(l(v_size_, v_size_));
      }

      // remove old values of rem var
      --v_size_;
      syx_plus_vs2 += v(v_size_) * v(v_size_);
      log_det_invQ -= 0.5 * log(inv_tau2_alpha2(v_size_));
      log_det_invQ_plus_xx -= log(l(v_size_, v_size_));

      if (DataModel::AH == type_rem){
        --v_size_;
        syx_plus_vs2 += v(v_size_) * v(v_size_);
        log_det_invQ -= 0.5 * log(inv_tau2_alpha2(v_size_));
        log_det_invQ_plus_xx -= log(l(v_size_, v_size_));
      }

      // update l, xy, v, type_inds, inv_tau2, x_ind1, x_ind2, x_types
      std::swap(type_inds[ind_keep], type_inds[ind_rem]);

      if (DataModel::AH == type_rem){
        if (DataModel::AH == type_keep){
          // both AH
          size_t swap_ind = x_ind2[ind_rem];
          l.cholesky_swapadj(swap_ind, &v);
          ++swap_ind;
          l.cholesky_swapadj(swap_ind, &v);
          swap_ind -= 2;
          l.cholesky_swapadj(swap_ind, &v);
          ++swap_ind;
          l.cholesky_swapadj(swap_ind, &v);

          std::swap(xy(x_ind1[ind_keep]), xy(x_ind2[ind_rem]));
          std::swap(xy(x_ind1[ind_keep]), xy(x_ind2[ind_keep]));
          std::swap(xy(x_ind1[ind_rem]), xy(x_ind2[ind_rem]));
          std::swap(xy(x_ind1[ind_keep]), xy(x_ind2[ind_rem]));

          std::swap(inv_tau2_alpha2(x_ind1[ind_keep]), inv_tau2_alpha2(x_ind2[ind_rem]));
          std::swap(inv_tau2_alpha2(x_ind1[ind_keep]), inv_tau2_alpha2(x_ind2[ind_keep]));
          std::swap(inv_tau2_alpha2(x_ind1[ind_rem]), inv_tau2_alpha2(x_ind2[ind_rem]));
          std::swap(inv_tau2_alpha2(x_ind1[ind_keep]), inv_tau2_alpha2(x_ind2[ind_rem]));


          // x_types and x_ind1/2 should be correct without any changes
//          std::swap(x_types[x_ind1[ind_keep]], x_types[x_ind2[ind_rem]]);
//          std::swap(x_types[x_ind1[ind_keep]], x_types[x_ind2[ind_keep]]);
//          std::swap(x_types[x_ind1[ind_rem]], x_types[x_ind2[ind_rem]]);
//          std::swap(x_types[x_ind1[ind_keep]], x_types[x_ind2[ind_rem]]);
        } else {
          // type_rem is AH
          size_t swap_ind = x_ind2[ind_rem];
          l.cholesky_swapadj(swap_ind, &v);
          --swap_ind;
          l.cholesky_swapadj(swap_ind, &v);

          std::swap(xy(x_ind2[ind_rem]), xy(x_ind1[ind_keep]));
          std::swap(xy(x_ind1[ind_rem]), xy(x_ind2[ind_rem]));
          std::swap(inv_tau2_alpha2(x_ind2[ind_rem]), inv_tau2_alpha2(x_ind1[ind_keep]));
          std::swap(inv_tau2_alpha2(x_ind1[ind_rem]), inv_tau2_alpha2(x_ind2[ind_rem]));
          std::swap(x_types[x_ind2[ind_rem]], x_types[x_ind1[ind_keep]]);
          std::swap(x_types[x_ind1[ind_rem]], x_types[x_ind2[ind_rem]]);

          x_ind2[ind_keep] = x_ind1[ind_keep];
          --x_ind1[ind_keep];
        }
      } else {
        if (DataModel::AH == type_keep){
          // type_keep is AH
          size_t swap_ind = x_ind1[ind_rem];
          l.cholesky_swapadj(swap_ind, &v);
          ++swap_ind;
          l.cholesky_swapadj(swap_ind, &v);

          std::swap(xy(x_ind1[ind_rem]), xy(x_ind1[ind_keep]));
          std::swap(xy(x_ind1[ind_keep]), xy(x_ind2[ind_keep]));
          std::swap(inv_tau2_alpha2(x_ind1[ind_rem]), inv_tau2_alpha2(x_ind1[ind_keep]));
          std::swap(inv_tau2_alpha2(x_ind1[ind_keep]), inv_tau2_alpha2(x_ind2[ind_keep]));
          std::swap(x_types[x_ind1[ind_rem]], x_types[x_ind1[ind_keep]]);
          std::swap(x_types[x_ind1[ind_keep]], x_types[x_ind2[ind_keep]]);

          x_ind2[ind_rem] = x_ind1[ind_keep];
          x_ind1[ind_keep] = x_ind2[ind_keep];
        } else {
          // both not AH
          size_t swap_ind = x_ind1[ind_rem];
          l.cholesky_swapadj(swap_ind, &v);

          std::swap(xy(x_ind1[ind_keep]), xy(x_ind1[ind_rem]));
          std::swap(inv_tau2_alpha2(x_ind1[ind_keep]), inv_tau2_alpha2(x_ind1[ind_rem]));
          std::swap(x_types[x_ind1[ind_keep]], x_types[x_ind1[ind_rem]]);
          // x_ind1 should be correct without any changes
        }
      }

      // update
      syx_plus_vs2 -= v(v_size_) * v(v_size_);
      log_det_invQ_plus_xx += log(l(v_size_, v_size_));
      ++v_size_;

      if (DataModel::AH == type_keep){
        syx_plus_vs2 -= v(v_size_) * v(v_size_);
        log_det_invQ_plus_xx += log(l(v_size_, v_size_));
        ++v_size_;
      }
      assert(syx_plus_vs2 > 0);

      // update log_likelihood
      log_likelihood_ = log_det_invQ - log_det_invQ_plus_xx
                  + prior_->minus_n_plus_nu_2 * log(syx_plus_vs2);

      log_model_prior_ += prior_->compute_log_change_on_rem(Ns, model_size_+1, type_rem);
      --Ns[type_rem];

      return log_likelihood_ + log_model_prior_;
    }

    double update_on_twonewswap()
    {
      update_on_add();
      update_on_add();
      return update_on_moveleft();
    }

  private:
    size_t model_size_, const_loci_, v_size_;
    double log_model_prior_;
};

} // namespace bmagwa

#endif /* MODEL_HPP_ */
