/* BMAGWA software v1.0
 *
 * data_model.cpp
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


#include "data_model.hpp"
#include "utils.hpp"

namespace bmagwa {

const char* DataModel::ef_name[] = {"A", "H", "D", "R", "AH"};

void DataModel::get_genotypes_additive(const size_t snp, VectorView v) const
{
  data_->get_genotypes_additive(snp, v);

  if (data_->miss_loc()[snp] == NULL) return;

  for (size_t i = 1; i <= data_->miss_loc()[snp][0]; ++i){
    v(data_->miss_loc()[snp][i]) = (double)miss_val_[snp][i];
  }
}

void DataModel::get_genotypes_heterozygous(const size_t snp, VectorView v) const
{
  data_->get_genotypes_heterozygous(snp, v);

  if (data_->miss_loc()[snp] == NULL) return;

  for (size_t i = 1; i <= data_->miss_loc()[snp][0]; ++i){
    v(data_->miss_loc()[snp][i]) = (double)(miss_val_[snp][i] == 1);
  }
}

void DataModel::get_genotypes_dominant(const size_t snp, VectorView v) const
{
  data_->get_genotypes_dominant(snp, v);

  if (data_->miss_loc()[snp] == NULL) return;

  for (size_t i = 1; i <= data_->miss_loc()[snp][0]; ++i){
    v(data_->miss_loc()[snp][i]) = (double)(miss_val_[snp][i] > 0);
  }
}

void DataModel::get_genotypes_recessive(const size_t snp, VectorView v) const
{
  data_->get_genotypes_recessive(snp, v);

  if (data_->miss_loc()[snp] == NULL) return;

  for (size_t i = 1; i <= data_->miss_loc()[snp][0]; ++i){
    v(data_->miss_loc()[snp][i]) = (double)(miss_val_[snp][i] == 2);
  }
}

/*
 * Samples missing data from prior for all SNPs except those having
 * model_ind >= 0 (i.e. is in model).
 */
void DataModel::sample_missing(const std::vector<int32_t>& model_inds,
                               Rand& rng)
{
  for (size_t snp = 0; snp < m_g; ++snp){
    // skip if no missing data or if locus in model
    if (miss_val_[snp] == NULL || model_inds[snp] >= 0) continue;

    for (size_t i = 1; i <= data_->miss_loc()[snp][0]; ++i){
      miss_val_[snp][i] = Utils::sample_discrete_naive(data_->miss_prior()[snp], 3,
                                                      rng);
    }
  }
}

/*
 * Samples missing data from prior for single SNPs.
 */
void DataModel::sample_missing_single(const uint32_t snp, Rand& rng)
{
  if (miss_val_[snp] == NULL) return;

  for (size_t i = 1; i <= data_->miss_loc()[snp][0]; ++i){
    miss_val_[snp][i] = Utils::sample_discrete_naive(data_->miss_prior()[snp], 3,
                                                    rng);
  }
}

void DataModel::update_xx_a(const size_t snp, SymmMatrixView& xx) const
{
  if (miss_val_[snp] == NULL) return;

  // update xx
  // x' * x (upper part is only accessed)
  for (size_t i = 1; i <= data_->miss_loc()[snp][0]; ++i){
    double val = (double)miss_val_[snp][i];

    // constant 1 * x_new (at 0,1) (old value was 0 for missing x_ij!)
    xx(0, 1) += val;
    // x_new * x_new (at 1,1)
    xx(1, 1) += (val * val);
  }
}

void DataModel::update_xx_h(const size_t snp, SymmMatrixView& xx) const
{
  if (miss_val_[snp] == NULL) return;

  // update xx
  // x' * x (upper part is only accessed)
  for (size_t i = 1; i <= data_->miss_loc()[snp][0]; ++i){
    double val = (double)(miss_val_[snp][i] == 1);

    // constant 1 * x_new (at 0,1) (old value was 0 for missing x_ij!)
    xx(0, 1) += val;
    // x_new * x_new (at 1,1)
    xx(1, 1) += (val * val);
  }
}

void DataModel::update_xx_d(const size_t snp, SymmMatrixView& xx) const
{
  if (miss_val_[snp] == NULL) return;

  // update xx
  // x' * x (upper part is only accessed)
  for (size_t i = 1; i <= data_->miss_loc()[snp][0]; ++i){
    double val = (double)(miss_val_[snp][i] > 0);

    // constant 1 * x_new (at 0,1) (old value was 0 for missing x_ij!)
    xx(0, 1) += val;
    // x_new * x_new (at 1,1)
    xx(1, 1) += (val * val);
  }
}

void DataModel::update_xx_r(const size_t snp, SymmMatrixView& xx) const
{
  if (miss_val_[snp] == NULL) return;

  // update xx
  // x' * x (upper part is only accessed)
  for (size_t i = 1; i <= data_->miss_loc()[snp][0]; ++i){
    double val = (double)(miss_val_[snp][i] == 2);

    // constant 1 * x_new (at 0,1) (old value was 0 for missing x_ij!)
    xx(0, 1) += val;
    // x_new * x_new (at 1,1)
    xx(1, 1) += (val * val);
  }
}

void DataModel::update_xx_ah(const size_t snp, SymmMatrixView& xx) const
{
  if (miss_val_[snp] == NULL) return;

  // update xx
  // x' * x (upper part is only accessed)
  for (size_t i = 1; i <= data_->miss_loc()[snp][0]; ++i){
    double val_a = (double)miss_val_[snp][i];
    double val_h = (double)(val_a == 1.0);

    // constant * additive at 1,2
    xx(0, 1) += val_a;
    // additive * additive at 2,2
    xx(1, 1) += (val_a * val_a);
    // constant * heterozygosity at 1,3
    xx(0, 2) += val_h;
    // additive * heterozygosity at 2,3
    xx(1, 2) += (val_a * val_h);
    // heterozygosity * heterozygosity at 3,3
    xx(2, 2) += (val_h * val_h);
  }
}

} // namespace bmagwa
