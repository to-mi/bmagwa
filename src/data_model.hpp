/* BMAGWA software v1.0
 *
 * data_model.hpp
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


#ifndef DATA_MODEL_HPP_
#define DATA_MODEL_HPP_

#include <vector>
#include "data.hpp"
#include "rand.hpp"

namespace bmagwa {

class DataModel
{
  public:
    // enums and typedefs
    enum ef_t {A = 0, H = 1, D = 2, R = 3, AH = 4}; // effect type
    static const char* ef_name[];

    typedef void (DataModel::*get_genotypes_fp)(const size_t snp, VectorView v)
    		                                                                  const;
    typedef void (DataModel::*update_xx_fp)(const size_t snp,
                                            SymmMatrixView& xx) const;

    // member variables
    const size_t n, m_g, m_e, n_types;
    const std::vector<ef_t> types;
    get_genotypes_fp get_genotypes[4];
    update_xx_fp update_xx[5];

    // constructors and destructors
    DataModel(const Data* data, const std::vector<ef_t> _types)
    : n(data->n), m_g(data->m_g), m_e(data->m_e), n_types(_types.size()),
      types(_types), data_(data), miss_val_(NULL)
    {
      get_genotypes[A] = &DataModel::get_genotypes_additive;
      get_genotypes[H] = &DataModel::get_genotypes_heterozygous;
      get_genotypes[D] = &DataModel::get_genotypes_dominant;
      get_genotypes[R] = &DataModel::get_genotypes_recessive;

      update_xx[A] = &DataModel::update_xx_a;
      update_xx[H] = &DataModel::update_xx_h;
      update_xx[D] = &DataModel::update_xx_d;
      update_xx[R] = &DataModel::update_xx_r;
      update_xx[AH] = &DataModel::update_xx_ah;

      for (size_t i = 0; i < 5; ++i){
        allow_types_[i] = false;
        allow_terms_[i] = false;
      }
      for (size_t i = 0; i < n_types; ++i){
        allow_types_[types[i]] = true;
        if (types[i] == AH){
          allow_terms_[A] = true;
          allow_terms_[H] = true;
        } else {
          allow_terms_[types[i]] = true;
        }
      }

      // initialize miss_val
      miss_val_ = new char*[m_g];

      for (size_t i = 0; i < m_g; ++i){
        if (data_->miss_loc()[i] != NULL){
          miss_val_[i] = new char[data_->miss_loc()[i][0] + 1];
          miss_val_[i][0] = 0; // skip first for same indexing with miss_loc
          for (size_t j = 1; j <= data_->miss_loc()[i][0]; ++j){
            miss_val_[i][j] = 0;
          }
        } else {
          miss_val_[i] = NULL;
        }
      }
    }

    ~DataModel()
    {
      if (miss_val_ != NULL){
        for (size_t i = 0; i < m_g; ++i){
          delete[] miss_val_[i];
        }
        delete[] miss_val_;
      }
    }

    char** miss_val() { return miss_val_; }
    const size_t*const* miss_loc() const { return data_->miss_loc(); }
    const double*const* miss_prior() const { return data_->miss_prior(); }
    bool genotype_missing(const size_t individual, const size_t snp) const {
      return (data_->get_genotype(individual, snp) < 0);
    }


    const bool allow_terms(DataModel::ef_t type) const
    {
      return allow_terms_[type];
    }
    const bool allow_types(DataModel::ef_t type) const
    {
      return allow_types_[type];
    }

    /*
     * These fetch genotypes. Missing values are filled from data_model.miss_val.
     */
    void get_genotypes_additive(const size_t snp, VectorView v) const;
    void get_genotypes_heterozygous(const size_t snp, VectorView v) const;
    void get_genotypes_dominant(const size_t snp, VectorView v) const;
    void get_genotypes_recessive(const size_t snp, VectorView v) const;

    /*
     * These are for updating the precomputed snp-wise covariance matrices.
     * Note: Only a local copy of xx should be updated, not the actual
     *       precomputed xx.
     *       This is because the updates assume that missing values were set to
     *       0 (and if they are not, the updates produce nonsense).
     */
    void update_xx_a(const size_t snp, SymmMatrixView& xx) const;
    void update_xx_h(const size_t snp, SymmMatrixView& xx) const;
    void update_xx_d(const size_t snp, SymmMatrixView& xx) const;
    void update_xx_r(const size_t snp, SymmMatrixView& xx) const;
    void update_xx_ah(const size_t snp, SymmMatrixView& xx) const;

    /*
     * These update miss_val (either for all SNPs with model_inds[snp] < 0 or
     * single SNP) from prior.
     */
    void sample_missing(const std::vector<int32_t>& model_inds, Rand& rng);
    void sample_missing_single(const uint32_t snp, Rand& rng);

    const Vector& y() const { return data_->y(); }

  private:
    const Data* data_;
    bool allow_types_[5];
    bool allow_terms_[5];
    char** miss_val_; // note: if miss_val[i] != NULL, miss_val[0] is the number
                      //       of missing values (same as in data->miss_loc)

    // disallow copy constructor and assignment
    DataModel(const DataModel& data_model);
    DataModel& operator=(const DataModel& data_model);
};

} // namespace bmagwa

#endif /* DATA_MODEL_HPP_ */
