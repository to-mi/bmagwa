/* BMAGWA software v2.0
 *
 * data_model.cpp
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

void DataModel::update_prexx_cov(const size_t snp, double* pre_xx) const
{
  if (miss_val_[snp] == NULL) return;

  double* xxloc = pre_xx;
  char* miss_val_snp = miss_val_[snp];
  char val = 0;
  double val_g = 0.0, sum_val_g, sum_val_g2, sum_a_times_sum_h = 0.0;

  if (allow_types(AH)){
    sum_a_times_sum_h = pre_xx[0] * pre_xx[2];
  }

  for (int t = 0; t < 4; ++t){
    if (!allow_terms((ef_t)t))
      continue;

    sum_val_g = 0.0;
    sum_val_g2 = 0.0;
    for (size_t i = 1; i <= data_->miss_loc()[snp][0]; ++i){
      val = miss_val_snp[i];
      switch (types[t]){
        case A:
          val_g = (double)val;
          break;
        case H:
          val_g = (double)(val == 1);
          break;
        case D:
          val_g = (double)(val > 0);
          break;
        case R:
          val_g = (double)(val == 2);
          break;
        default:
          throw std::logic_error("type not in A,H,D,R");
      }
      sum_val_g += val_g;
      sum_val_g2 += val_g * val_g;
    }
    // update sum (contribution was 0)
    *xxloc += sum_val_g;
    xxloc++;

    // update var
    //*xxloc += val_g * val_g - (2 * (*xxloc) * val_g + val_g * val_g)/n;
    *xxloc += sum_val_g2 - sum_val_g * ((*xxloc) * 2.0 + sum_val_g) / (double)n;
    xxloc++;
  }
  if (allow_types(AH)){
    *xxloc += (sum_a_times_sum_h - pre_xx[0] * pre_xx[2]) / (double)n;
    for (size_t i = 1; i <= data_->miss_loc()[snp][0]; ++i){
      val = miss_val_snp[i];

      // addition: val * (val == 1)
      if (val == 1){
        *xxloc += (double)val;
      }
    }
    // this is last, no need to xxloc++
  }
}

} // namespace bmagwa
