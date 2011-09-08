/* BMAGWA software v1.0
 *
 * data.hpp
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


#ifndef DATA_HPP_
#define DATA_HPP_

#include <cmath>
#include <map>
#include "linalg.hpp"

namespace bmagwa {

// Data needs to be thread-safe
// (achieved via being read-only after construction)
class Data
{
  public:
    const size_t n, m_g, m_e;

    // constructor
    Data(size_t n_, size_t m_g_, size_t m_e_, const std::string file_fam,
        const std::string file_g, const bool recode_g_to_ma_count,
        const std::string file_e, const std::string file_y = "")
    : n(n_), m_g(m_g_), m_e(m_e_ + 1),
      _y(n), _e(n, m_e), _g(NULL), miss_loc_(NULL), miss_prior_(NULL),
      _snp_offset((size_t)ceil((double)n_ / 4.0)),
      _n_4((size_t)floor((double)n / 4.0))
    {
      _e = 1; // note m_e += 1 in init list to include a column of ones!

      read_fam(file_fam);
      if (!file_y.empty())
        read_y(file_y);
      if (m_e > 1)
        read_e(file_e);
      read_g(file_g);

      if (recode_g_to_ma_count) recode_g_to_minor_allele_count();

      handle_missing_g();

      // precompute
      _yy = VectorView::dotproduct(_y, _y);
      //_var_y = (_y - _y.mean()).square().sum() / (n - 1);
      _var_y = _y.var();
      //_var_x = compute_g_var();
      compute_g_var_and_mean(_mean_x, _var_x);
    }

    // methods
    inline const Vector& y() const { return _y; }
    inline const Matrix& e() const { return _e; }
    inline const double var_x() const { return _var_x; }
    inline const double mean_x() const { return _mean_x; }
    inline const double var_y() const { return _var_y; }
    inline const double yy() const { return _yy; }
    const size_t*const* miss_loc() const { return miss_loc_; }
    const double*const* miss_prior() const { return miss_prior_; }
    const double get_genotype(const size_t individual, const size_t snp) const;

    void get_genotypes_additive(const size_t snp, VectorView v) const;
    void get_genotypes_heterozygous(const size_t snp, VectorView v) const;
    void get_genotypes_dominant(const size_t snp, VectorView v) const;
    void get_genotypes_recessive(const size_t snp, VectorView v) const;

    // destructor
    ~Data()
    {
      if (miss_loc_ != NULL) {
        for (size_t snp = 0; snp < m_g; ++snp)
          delete[] miss_loc_[snp];
        delete[] miss_loc_;
      }
      if (miss_prior_ != NULL) {
        for (size_t snp = 0; snp < m_g; ++snp)
          delete[] miss_prior_[snp];
        delete[] miss_prior_;
      }
      delete[] _g;
    }

  private:
    Vector _y;
    Matrix _e;
    unsigned char *_g;

    size_t** miss_loc_;
    double** miss_prior_;

    double _var_x, _mean_x;
    double _var_y;
    double _yy; // y' * y

    size_t _snp_offset; // ceil(n / 4)
    size_t _n_4; // floor(n / 4)
    // maps from FID, IID to fam-position (not really required to keep this in
    // memory after constructor...)
    std::map<std::pair<std::string, std::string>, size_t> fam_mapping;

    // genotype look-up tables
    static const double _SINGLE_GENOTYPE_TABLE[4]; // additive for single snp
    static const double _ADDITIVE_GENOTYPE_TABLE[256][4];
    static const double _HETEROZYGOUS_GENOTYPE_TABLE[256][4];
    static const double _DOMINANT_GENOTYPE_TABLE[256][4];
    static const double _RECESSIVE_GENOTYPE_TABLE[256][4];

    // disallow copy constructor and assignment
    Data(const Data& data);
    Data& operator=(const Data& data);

    // methods
    void read_fam(const std::string& file_fam);
    void read_y(const std::string& file_y);
    void read_e(const std::string& file_e);
    void read_g(const std::string& file_g);
    void set_genotype(const size_t individual, const size_t snp,
                      const double genotype);
    void recode_g_to_minor_allele_count();
    void handle_missing_g();
    double compute_g_var() const;
    void compute_g_var_and_mean(double&  mean, double& var) const;
};

} // namespace bmagwa

#endif /* DATA_HPP_ */
