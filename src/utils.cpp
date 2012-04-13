/* BMAGWA software v2.0
 *
 * utils.cpp
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

int sample_discrete_naive_ft(const double *cumsum, const int from, const int to, Rand &rng){
  double d = ((from > 0) ? cumsum[from-1] : 0.0);
  double z = cumsum[to] - d;
  double r = rng.rand_01() * z + d;

  for (int i = from; i <= to; ++i){
    if (r <= cumsum[i]){
      return i;
    }
  }

  throw std::logic_error("sample_discrete_naive_ft reached end, exiting");
}

size_t sample_discrete(const double* cumsum, const size_t m, int level, Rand& rng)
{
  size_t a = 0;
  size_t b = m - 1;
  size_t c;

  double z = cumsum[b];
  double r = rng.rand_01() * z;

  // bisection search
  while (level > 0){
    c = (a + b) / 2;
    if (r <= cumsum[c]){
      b = c;
    } else {
      a = c + 1;
    }
    --level;
  }

  // linear search
  for (size_t i = a; i <= b; ++i){
    if (r <= cumsum[i]){
      return i;
    }
  }

  throw std::logic_error("sample_discrete reached end");
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

double lognchoosek(size_t n, size_t k)
{
  return gammaln(n + 1) - gammaln(n - k + 1) - gammaln(k + 1);
}

// p = prob of success
// values[i] = prob of first occurence of success on i + 1 trial
// support {1,...,maxsize}, i.e. truncated at maxsize (but not normalized)
void geometric_dist_pdf(const int maxsize, const double p, double *values)
{
  double logp = log(p);
  double log1mp = log1p(-p);
  for (int i = 0; i < maxsize; ++i){
    values[i] = exp(i * log1mp + logp);
  }
}

void geometric_dist_cdf(const int maxsize, const double p, double *values)
{
  double q = (1 - p);
  double qp = q;
  for (int i = 0; i < maxsize; ++i){
    values[i] = 1 - q;
    q *= qp;
  }
}

void binomial_cdf_phalf(const int trials, double *values)
{
  // normalized by mode (and leave the cdf to sum to more than one)
  double norm = lognchoosek(trials, trials/2);
  values[0] = std::exp(lognchoosek(trials, 0) - norm);
  for(int i = 1; i <= trials; ++i){
    values[i] = values[i-1] + std::exp(lognchoosek(trials, i) - norm);
  }
}

} // namespace Utils
} // namespace bmagwa
