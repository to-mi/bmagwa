/* BMAGWA software v2.0
 *
 * utils.hpp
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


#ifndef UTILS_HPP_
#define UTILS_HPP_

namespace bmagwa {

// forward declarations
class Rand;
class VectorView;

namespace Utils {

double compute_allele_freq(const VectorView snp_genotypes);
int sample_discrete_naive(const double* cumsum, const int m, Rand& rng);
int sample_discrete_naive_ft(const double* cumsum, const int from, const int to, Rand& rng);
size_t sample_discrete(const double* cumsum, const size_t m, int level, Rand& rng);
double gammaln(const double val);
double sinvchi2_lpdf(const double x, const double nu, const double s2);
std::ofstream* openfile(const std::string filename, const bool binary);
double lognchoosek(size_t n, size_t k);
void geometric_dist_pdf(const int maxsize, const double p, double *values);
void geometric_dist_cdf(const int maxsize, const double p, double *values);
void binomial_cdf_phalf(const int trials, double *values);

} // namespace Utils
} // namespace bmagwa

#endif /* UTILS_HPP_ */
