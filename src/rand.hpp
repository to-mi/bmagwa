/* BMAGWA software v1.0
 *
 * rand.hpp
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


#ifndef RAND_HPP_
#define RAND_HPP_

#include <boost/random.hpp>
#include "linalg.hpp"

namespace bmagwa {

//! Random number generation facilities.
/*!
 *  Random number generation using Mersenne twister for uniform on 0...1,
 *  standard normal, scaled inverse-Chi^2 and multivariate normal.
 */
class Rand
{
  public:
    //! Sets up the generators.
    /*!
     *  \param seed Seed for the base generator.
     *  \param sinvchi2_nu Degrees of freedom for the scaled inv-Chi^2
     *         generator. Scale is given as an argument on generator call.
     */
    Rand(const uint32_t seed, const double sinvchi2_nu)
    : _sinvchi2_nu(sinvchi2_nu), _gamma_alpha(0.5 * sinvchi2_nu),
      rng(seed),
      norm_dist(0, 1), unif_dist(), gamma_dist(_gamma_alpha),
      normal_sampler(rng, norm_dist), unif_sampler(rng, unif_dist),
      gamma_sampler(rng, gamma_dist)
    {}

    //! Random value from standard normal.
    inline double rand_normal() { return normal_sampler(); }
    //! Random value from uniform on 0...1.
    inline double rand_01() { return unif_sampler(); }
    //! Random number from scaled inverse-Chi^2.
    /*!
     *  \param s2 Scale.
     */
    double rand_sinvchi2(const double s2)
    {
      /*
       * Sampling from scaled-Inv-Chi2(nu, s2):
       *  1) sample X from Chi2(nu)
       *   1.1) Chi2(nu) is Gamma(alpha = nu/2, beta = 0.5)
       *                     = 2 * Gamma(nu/2, beta = 1)
       *  2) return nu * s2 / X
       * ref. Gelman et al. 1995
       */
      return _sinvchi2_nu * s2 / (2.0 * gamma_sampler());
    }
    //! Random number from scaled inverse-Chi^2.
    /*!
     *  \param nu Degrees of freedom.
     *  \param s2 Scale.
     */
    double rand_sinvchi2(const double nu, const double s2)
    {
      boost::gamma_distribution<double> gamma_dist_(0.5 * nu);
      boost::variate_generator<boost::mt19937&,
          boost::gamma_distribution<double> > gamma_sampler_(rng, gamma_dist_);

      return nu * s2 / (2.0 * gamma_sampler_());
    }

    //! Generates a draw from multivariate normal distribution.
    /*!
     *  \param mu Mean vector.
     *  \param U Upper cholesky factor of Sigma_unscaled (i.e.
     *         U' * U = Sigma / scale^2).
     *  \param scale (i.e. Sigma = scale^2 * Sigma_unscaled)
     *  \param x Vector(View) of length d (dimension), where the random number
     *           will be put into.
     */
    void rand_mvnormal_chol(const VectorView& mu,
        const SymmMatrixView& U, const double scale,
        VectorView& x)
    {
      /*
       * Sampling multivariate normal;
       *  1) find U' U = Sigma
       *  2) sample z (independent standard normal)
       *  3) return U' z + mu
       *  ref. Gelman et al. 1995
       */

      // independent standard normal variates
      for (size_t i = 0; i < x.length(); ++i){
        x(i) = rand_normal();
      }

      // multiply by lower cholesky, though we have upper,
      // so need to transpose it
      x.multiply_by_triangularmatrix(U, true);

      // scale and add mean
      for (size_t i = 0; i < x.length(); ++i){
        x(i) = scale * x(i) + mu(i);
      }
    }

    //! Generates a draw from multivariate normal distribution.
    /*!
     *  \param mu Mean vector.
     *  \param U Upper Cholesky factor of Sigma_unscaled^-1 (i.e.
     *         U' * U = Sigma^-1 * scale^2).
     *  \param scale (i.e. Sigma = scale^2 * Sigma_unscaled)
     *  \param x Vector(View) of length d (dimension), where the random number
     *           will be put into.
     */
    void rand_mvnormal_invchol(const VectorView& mu,
        const SymmMatrixView& Uinv,
        const double scale,
        VectorView& x)
    {
      /*
       * Sampling multivariate normal;
       *  1) find (U' U)^-1 = U^-1 U'^-1 = Sigma
       *  2) sample z (independent standard normal)
       *  3) return U^-1 z + mu
       *  ref. Gelman et al. 1995
       */

      // independent standard normal variates
      for (size_t i = 0; i < x.length(); ++i){
        x(i) = rand_normal();
      }

      x.multiply_by_invtriangularmatrix(Uinv, false);

      // scale and add mean
      for (size_t i = 0; i < x.length(); ++i){
        x(i) = scale * x(i) + mu(i);
      }
    }


  private:
    const double _sinvchi2_nu, _gamma_alpha;

    boost::mt19937 rng;
    boost::normal_distribution<double> norm_dist;
    boost::uniform_01<double> unif_dist;
    boost::gamma_distribution<double> gamma_dist;
    boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<double> >
    normal_sampler;
    boost::variate_generator<boost::mt19937&,
    boost::uniform_01<double> >
    unif_sampler;
    boost::variate_generator<boost::mt19937&,
    boost::gamma_distribution<double> >
    gamma_sampler;

    // disallow copy constructor and assignment
    Rand(const Rand& rng);
    Rand& operator=(const Rand& rng);
};

} // namespace bmagwa

#endif /* RAND_HPP_ */
