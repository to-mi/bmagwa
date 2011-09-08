/* BMAGWA software v1.0
 *
 * model.hpp
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

namespace bmagwa {

// forward declaration
class RaoBlackwellizer;

class Model
{
  public:
    Model(const Data* data, const DataModel* data_model, const Prior* prior)
    : data_(data), data_model_(data_model), prior_(prior),
      n_model_inds(data->m_g + data->m_e),
      x(data->n, MAXMODELSIZE, data->n, data->m_e), xx(MAXMODELSIZE, data->m_e),
      l(MAXMODELSIZE, data->m_e),
      xy(MAXMODELSIZE, data->m_e), v(MAXMODELSIZE, data->m_e),
      y_hat(data->n), y_hat_tmp_e(data->n),
      c_tmp(NULL), s_tmp(NULL), beta(MAXMODELSIZE, data->m_e),
      mu_beta(MAXMODELSIZE, data->m_e), sigma2(0), syx_plus_vs2(0),
      log_det_invQ(0),
      log_det_invQ_plus_xx(0), beta_sigma2_sampled(false), log_likelihood_(NAN),
      loci(0), x_ind_to_g_ind(data->m_e, 0), model_inds(n_model_inds, -1),
      type_inds(0), x_types(0), x_ind1(0), x_ind2(0)
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
        beta_sigma2_sampled = src.beta_sigma2_sampled;
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

    ~Model() {
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
      if (beta_sigma2_sampled)
        return sigma2;
      else
        return NAN;
    }
    const VectorView get_beta() const {
      return beta;
    }
    // for testing purposes mostly:
    const int* get_Ns() const {
      return Ns;
    }


    void compute_log_likelihood()
    {
      if (x.cols() < 1){
        log_likelihood_ = NAN;
        return;
      }

      l = xx;
      log_det_invQ = prior_->do_invQ(l, x_types);

      if (!l.cholesky()){
        /*std::cout << "Warning, chol: xx is not positive definitive. QTL:";
        for (size_t i = 0; i < loci.size(); ++i){
          std::cout << " " << loci[i];
        }
        std::cout << std::endl;*/
        log_likelihood_ = -INFINITY;
        syx_plus_vs2 = INFINITY;
        beta_sigma2_sampled = false;
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
      beta_sigma2_sampled = false;
    }

    void add_term(const uint32_t ind, const char t_ind)
    {
      DataModel::ef_t type = data_model_->types[t_ind];
      bool type_is_AH = (type == DataModel::AH);

      size_t model_ind = loci.size();

      if (x.cols() + 1 + (size_t)type_is_AH <= MAXMODELSIZE){
        loci.push_back(ind);
        ++Ns[type];
        type_inds.push_back(t_ind);
        model_inds[ind] = model_ind;
        beta_sigma2_sampled = false;

        if (type_is_AH){
          update_likelihood_on_add(ind, model_ind, DataModel::A, false);
          if (std::isfinite(syx_plus_vs2))
            update_likelihood_on_add(ind, model_ind, DataModel::H, true);
        } else {
          update_likelihood_on_add(ind, model_ind, type, false);
        }
      } else {
        throw std::runtime_error("Not enough memory for adding term.");
      }

      assert(x.cols() == v.length());
      assert(x.cols() == xy.length());
      assert(x.cols() == l.length());
    }

    void remove_term(const int model_ind)
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
          update_likelihood_on_remove(xi1);
          update_likelihood_on_remove(xi2);
        } else {
          update_likelihood_on_remove(xi2);
          update_likelihood_on_remove(xi1);
        }
      } else {
        update_likelihood_on_remove(xi1);
      }

      beta_sigma2_sampled = false;

      assert(x.cols() == v.length());
      assert(x.cols() == xy.length());
      assert(x.cols() == l.length());
    }

    void sample_beta_sigma2(Rand& rng)
    {
      // sample beta and sigma2
      sigma2 = rng.rand_sinvchi2(syx_plus_vs2 / prior_->n_plus_nu);

      assert(v.length() == x.cols());

      // find mean of beta
      if (!beta_sigma2_sampled){
        // update mu
        mu_beta.resize(v.length());
        mu_beta = v;
        mu_beta.multiply_by_invtriangularmatrix(l, false);
      }

      // sample beta
      beta.resize(mu_beta.length());
      rng.rand_mvnormal_invchol(mu_beta, l, sqrt(sigma2), beta);

      beta_sigma2_sampled = true;

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

  private:
    const Data* data_;
    const DataModel* data_model_;
    const Prior* prior_;

    //size_t n;
    const uint32_t n_model_inds;

    Matrix x;
    SymmMatrix xx, l;
    Vector xy, v;
    Vector y_hat, y_hat_tmp_e; // pre-allocated tmps for pve computations
    double* c_tmp;
    double* s_tmp;
    Vector beta, mu_beta;
    double sigma2;
    double syx_plus_vs2, log_det_invQ, log_det_invQ_plus_xx;
    int Ns[5];
    bool beta_sigma2_sampled;
    double log_likelihood_;

    std::vector<uint32_t> loci;          // length nloci
    std::vector<uint32_t> x_ind_to_g_ind;// length x.length()
    std::vector<int32_t> model_inds;     // length m_g
    std::vector<char> type_inds;         // length nloci
    std::vector<DataModel::ef_t> x_types; // length x.length()
    std::vector<unsigned int> x_ind1;    // length nloci
    std::vector<unsigned int> x_ind2;    // length nloci

    Model(const Model& m)
    : data_(m.data_), data_model_(m.data_model_), prior_(m.prior_),
      n_model_inds(0), x(0, 0), xx(0), l(0), xy(0), v(0),
      y_hat(0), y_hat_tmp_e(0), c_tmp(NULL),
      s_tmp(NULL), beta(0), mu_beta(0)
    {}

    void update_likelihood_on_add(const uint32_t ind,
                                  const size_t model_ind,
                                  const DataModel::ef_t type,
                                  const bool is_x_ind2)
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

      // update xx
      xx.resize(col + 1);
      VectorView xxcol(xx.column(col));
      xxcol.set_to_product(x, xcol, true);

      // update cholesky
      if (!l.cholesky_update(xxcol, 1 / prior_->get_tau2_value(type))){
        //std::cout << "Warning, chol update: xx is not positive definitive."
        //          << " QTL: " << loci[model_ind] << std::endl;
        log_likelihood_ = -INFINITY;
        syx_plus_vs2 = INFINITY;
        beta_sigma2_sampled = false;
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

      log_det_invQ -= 0.5 * log(prior_->get_tau2_value(type));
      log_det_invQ_plus_xx += log(lcol(col));

      log_likelihood_ = log_det_invQ - log_det_invQ_plus_xx
                       + prior_->minus_n_plus_nu_2 * log(syx_plus_vs2);
    }

    void update_likelihood_on_remove(const unsigned int x_ind)
    {
      x.remove_column(x_ind);
      x_ind_to_g_ind.erase(x_ind_to_g_ind.begin() + x_ind);
      xy.remove(x_ind);
      xx.remove_colrow(x_ind);
      DataModel::ef_t type = x_types[x_ind];
      x_types.erase(x_types.begin() + x_ind);

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

      log_det_invQ += 0.5 * log(prior_->get_tau2_value(type));

      log_likelihood_ = log_det_invQ - log_det_invQ_plus_xx
                       + prior_->minus_n_plus_nu_2 * log(syx_plus_vs2);
    }

    friend class RaoBlackwellizer;
    friend class Sampler;

};

} // namespace bmagwa

#endif /* MODEL_HPP_ */
