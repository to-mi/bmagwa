/* BMAGWA software v1.0
 *
 * utils.hpp
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


#ifndef UTILS_HPP_
#define UTILS_HPP_

namespace bmagwa {

// forward declarations
class Rand;
class VectorView;

namespace Utils {

double compute_allele_freq(const VectorView snp_genotypes);
int sample_discrete_naive(const double *cumsum, const int m, Rand &rng);
double gammaln(const double val);
double sinvchi2_lpdf(const double x, const double nu, const double s2);
std::ofstream* openfile(const std::string filename, const bool binary);

} // namespace Utils
} // namespace bmagwa

#endif /* UTILS_HPP_ */
