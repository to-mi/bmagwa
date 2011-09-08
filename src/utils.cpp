/* BMAGWA software v1.0
 *
 * utils.cpp
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


#include <stdexcept>
#include <fstream>
#include <boost/math/special_functions/gamma.hpp>
#include "utils.hpp"
#include "rand.hpp"
#include "vector.hpp"

namespace bmagwa {
namespace Utils {

double compute_allele_freq(const VectorView snp_genotypes)
{
  double macount = 0;
  size_t ngenos = 0;
  for (size_t i = 0; i < snp_genotypes.length(); ++i){
    if (snp_genotypes(i) >= 0){
      macount += snp_genotypes(i);
      ++ngenos;
    }
  }
  return macount / (2.0 * ngenos);
}

int sample_discrete_naive(const double *cumsum, const int m, Rand &rng){
  double z = cumsum[m - 1];
  double r = rng.rand_01() * z;

  for (int i = 0; i < m; ++i){
    if (r <= cumsum[i]){
      return i;
    }
  }

  throw std::logic_error("sample_discrete_naive reached end, exiting");
}

double gammaln(const double val)
{
  return boost::math::lgamma(val);
}

double sinvchi2_lpdf(const double x, const double nu, const double s2){
  double nu_2 = 0.5 * nu;

  return log(nu_2) * nu_2 - gammaln(nu_2) + log(s2) * nu_2
		  - log(x) * (nu_2 + 1) - nu_2 * s2 / x;
}

std::ofstream* openfile(const std::string filename, const bool binary)
{
  std::ofstream* tmp = new std::ofstream;
  if (binary)
    tmp->open(filename.c_str(), std::ios::binary);
  else
    tmp->open(filename.c_str());
  if (!tmp->is_open())
    throw std::runtime_error("Failed to open file: " + filename);
  return tmp;
}


} // namespace Utils
} // namespace bmagwa
