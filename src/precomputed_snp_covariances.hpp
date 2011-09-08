/* BMAGWA software v1.0
 *
 * precomputed_snp_covariances.hpp
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


#ifndef PRECOMPUTED_SNP_COVARIANCES_HPP_
#define PRECOMPUTED_SNP_COVARIANCES_HPP_

#include "data_model.hpp"
#include "linalg.hpp"

namespace bmagwa {

// forward declaration
class RaoBlackwellizer;

class PrecomputedSNPCovariances
{
  public:
	PrecomputedSNPCovariances(const DataModel* data_model)
	{
	  const size_t m_g = data_model->m_g;

	  // allocate xx and xx_types
	  for (size_t i = 0; i < 5; ++i) xx_types[i] = NULL;
	  xx = new double*[data_model->n_types];

	  if (data_model->allow_types(DataModel::AH)){
	    size_t j = 0;
	    xx[j] = new double[m_g * 3 * 3];
	    xx_offset[DataModel::AH] = 3;
	    xx_types[DataModel::AH] = xx[j];

	    if (data_model->allow_types(DataModel::A)){
	      // allocated xx for A and AH
	      xx_offset[DataModel::A] = 3;
	      xx_types[DataModel::A] = xx[j];
	    }

	    // allocate for others (skipping first and last; these should be A and AH, respectively as !)
	    j = 1;
	    for (size_t i = 1; i < 4; ++i){
	      if (!(data_model->allow_types((DataModel::ef_t)i))) continue;
	      xx[j] = new double[m_g * 2 * 2];
	      xx_offset[i] = 2;
	      xx_types[i] = xx[j];
	      ++j;
	    }
	  } else {
	    size_t j = 0;
	    for (size_t i = 0; i < 4; ++i){
	      if (!(data_model->allow_types((DataModel::ef_t)i))) continue;
	      xx[j] = new double[m_g * 2 * 2];
	      xx_offset[i] = 2;
	      xx_types[i] = xx[j];
	      ++j;
	    }
	  }

	  // precompute
	  precompute(data_model);
	}


    ~PrecomputedSNPCovariances()
    {
      if (xx != NULL){
        for (size_t i = 0; i < 5; ++i){
          // use last xx_types pointer only to free memory
          for (size_t j = i + 1; j < 5; ++j){
            if (xx_types[i] == xx_types[j]) xx_types[i] = NULL;
          }
          delete[] xx_types[i];
        }
        delete[] xx;
      }
    }
  private:
    double **xx;         // in order of data_model->types
    double *xx_types[5]; // in order of ef_t enum
    size_t xx_offset[5]; // in order of ef_t enum

    void precompute(const DataModel* data_model)
    {
      size_t n = data_model->n;

      Matrix x(n, 2 + data_model->allow_types(DataModel::AH));
      x = 1;

      if (data_model->allow_types(DataModel::AH)){
        for (size_t i = 0; i < data_model->m_g; ++i){
          size_t j = 0;
          {
            // AH first
            x.resize(n, 3);

            data_model->get_genotypes_additive(i, x.column(1));
            data_model->get_genotypes_heterozygous(i, x.column(2));


            SymmMatrixView tmp_xx(xx[j] + i * 9, 3);
            tmp_xx.set_to_innerproduct(x);
          }

          // then others (skipping A)
          x.resize(n, 2);
          j = 1;
          for (size_t k = 1; k < 4; ++k){
            if (!(data_model->allow_types((DataModel::ef_t)k))) continue;

            (data_model->*data_model->bmagwa::DataModel::get_genotypes[k])(i,
            		                                               x.column(1));

            SymmMatrixView tmp_xx(xx[j] + i * 4, 2);
            tmp_xx.set_to_innerproduct(x);

            ++j;
          }
        }
      } else {
        for (size_t i = 0; i < data_model->m_g; ++i){
          size_t j = 0;
          for (size_t k = 0; k < 4; ++k){
            if (!data_model->allow_types((DataModel::ef_t)k)) continue;
            (data_model->*data_model->bmagwa::DataModel::get_genotypes[k])(i,
                                                                   x.column(1));

            SymmMatrixView tmp_xx(xx[j] + i * 4, 2);
            tmp_xx.set_to_innerproduct(x);

            ++j;
          }
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
