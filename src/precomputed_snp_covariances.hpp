/* BMAGWA software v2.0
 *
 * precomputed_snp_covariances.hpp
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


#ifndef PRECOMPUTED_SNP_COVARIANCES_HPP_
#define PRECOMPUTED_SNP_COVARIANCES_HPP_

#include "data_model.hpp"
#include "linalg.hpp"

namespace bmagwa {

// forward declaration
class RaoBlackwellizer;

//! Precomputes and holds SNP (co)variances for use in RB regression.
/*!
 *  Computes and holds the (co)variances required in RB regression. That is,
 *  no cross-covariances between different SNPs are computed. Note, the pre-
 *  computation needs to be done at a stage where all missing data is set to
 *  zeros, so that DataModel.update_prexx_cov can be used to account for the
 *  current sampled values for the missing data.
 */
class PrecomputedSNPCovariances
{
  public:
    //! Length of the elements for single SNP
    const int offset;
    //! Offsets for effect term types in the order of DataModel::ef_t (only allowed terms)
    int offset_type[4];
    //! The precomputed values
    /*!
     *  Elements for jth SNP begin at j * offset. Offset for tth type is added
     *  to this (e.g., j * offset + offset_type[t]). For A,H,D,R values include
     *  n * sample mean and unnormalized sample variance. Unnormalized
     *  covariance between A and H is included if AH is allowed.
     */
    double* xx;         // in order of data_model->types

    PrecomputedSNPCovariances(const DataModel* data_model)
    : offset(2 * data_model->n_types +
             // if AH, then space for the crossterm and possibly A or H is needed
             //        if they are not allowed individually. note that AH is counted
             //        in n_types!
             (data_model->allow_types(DataModel::AH) ? 1 : 0)     // if AH
             * (1 +                                               // cross term
                (data_model->allow_types(DataModel::A) &&         // if both -2
                 data_model->allow_types(DataModel::H)     ? -2 : 0) +
                (data_model->allow_types(DataModel::A) ||         // if neither +2
                 data_model->allow_types(DataModel::H)     ? 0 : 2))),
      xx(NULL)
    {
      const size_t m_g = data_model->m_g;

      xx = new double[m_g * offset];

      int i = 0;
      for (int j = 0; j < 4; ++j){
        if (data_model->allow_terms((DataModel::ef_t)j)){
          offset_type[j] = i;
          ++i;
        } else {
          offset_type[j] = -1;
        }
      }

      // precompute
      precompute(data_model);
    }


    ~PrecomputedSNPCovariances()
    {
      delete[] xx;
    }
  private:
    void precompute(const DataModel* data_model)
    {
      size_t n = data_model->n;

      Vector x(n);
      x = 0;
      Vector xAH(n);
      xAH = 0;

      double* xxloc = xx;

      for (size_t i = 0; i < data_model->m_g; ++i){
        assert(xxloc == xx + i * offset);
        for (int t = 0; t < 4; ++t){
          DataModel::ef_t type = (DataModel::ef_t)t;
          if (!data_model->allow_terms(type))
            continue;

          (data_model->*data_model->bmagwa::DataModel::get_genotypes[type])(i, x);
          double xsum = x.sum();
          *xxloc = xsum;
          ++xxloc;
          *xxloc = Vector::dotproduct(x, x) - xsum * xsum / (double)n;
          ++xxloc;
        }
        if (data_model->allow_types(DataModel::AH)){
          (data_model->*data_model->bmagwa::DataModel::get_genotypes[DataModel::A])(i, x);
          (data_model->*data_model->bmagwa::DataModel::get_genotypes[DataModel::H])(i, xAH);

          double asum = x.sum();
          double hsum = xAH.sum();

          *xxloc = Vector::dotproduct(x, xAH) - asum * hsum / (double)n;
          ++xxloc;
        }
      }
    }

    // disallow copy constructor and assignment
    PrecomputedSNPCovariances(const PrecomputedSNPCovariances& precov);
    PrecomputedSNPCovariances& operator=(const PrecomputedSNPCovariances& precov);

    friend class RaoBlackwellizer;
};

} // namespace bmagwa

#endif /* PRECOMPUTED_SNP_COVARIANCES_HPP_ */
