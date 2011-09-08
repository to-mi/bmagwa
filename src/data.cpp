/* BMAGWA software v1.0
 *
 * data.cpp
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


#include <vector>
#include <fstream>
#include <utility>
#include <stdexcept>
#include "data.hpp"
#include "utils.hpp"

/*
 TODO: read fam, y, e by line to string stream and then to values.(more robust)
 */

namespace bmagwa {

const double Data::_SINGLE_GENOTYPE_TABLE[4] = { 0.0, -1.0, 1.0, 2.0 };
#include "genotype_tables.hpp"

// public:
const double Data::get_genotype(const size_t individual, const size_t snp) const
{
  size_t offset = _snp_offset * snp + floor(individual / 4);
  // get the correct byte, shift by 2 * (individual % 4) and take the two bits
  unsigned char byte = (_g[offset] >> (2 * (individual % 4))) & 0x03;

  /*switch (byte){
     case 0x01: return -1.0;
     case 0x00: return 0.0;
     case 0x02: return 1.0;
     case 0x03: return 2.0;
     default: throw std::runtime_error("get_genotype failed");
     }*/
  return _SINGLE_GENOTYPE_TABLE[byte];
}

void Data::get_genotypes_additive(const size_t snp, VectorView v) const
{
  size_t offset = _snp_offset * snp;
  size_t i = 0;
  for (size_t ind4 = 0; ind4 < _n_4; ++ind4) {
    const double *tbl = _ADDITIVE_GENOTYPE_TABLE[_g[offset++]];
    v(i++) = tbl[0];
    v(i++) = tbl[1];
    v(i++) = tbl[2];
    v(i++) = tbl[3];
  }
  // rest of the individuals
  if (i < n) {
    size_t nr = n - i;
    const double *tbl = _ADDITIVE_GENOTYPE_TABLE[_g[offset]];
    for (size_t ind = 0; ind < nr; ++ind) {
      v(i++) = tbl[ind];
    }
  }
}

void Data::get_genotypes_heterozygous(const size_t snp, VectorView v) const
{
  size_t offset = _snp_offset * snp;
  size_t i = 0;
  for (size_t ind4 = 0; ind4 < _n_4; ++ind4) {
    const double *tbl = _HETEROZYGOUS_GENOTYPE_TABLE[_g[offset++]];
    v(i++) = tbl[0];
    v(i++) = tbl[1];
    v(i++) = tbl[2];
    v(i++) = tbl[3];
  }
  // rest of the individuals
  if (i < n) {
    size_t nr = n - i;
    const double *tbl = _HETEROZYGOUS_GENOTYPE_TABLE[_g[offset]];
    for (size_t ind = 0; ind < nr; ++ind) {
      v(i++) = tbl[ind];
    }
  }
}

void Data::get_genotypes_dominant(const size_t snp, VectorView v) const
{
  size_t offset = _snp_offset * snp;
  size_t i = 0;
  for (size_t ind4 = 0; ind4 < _n_4; ++ind4) {
    const double *tbl = _DOMINANT_GENOTYPE_TABLE[_g[offset++]];
    v(i++) = tbl[0];
    v(i++) = tbl[1];
    v(i++) = tbl[2];
    v(i++) = tbl[3];
  }
  // rest of the individuals
  if (i < n) {
    size_t nr = n - i;
    const double *tbl = _DOMINANT_GENOTYPE_TABLE[_g[offset]];
    for (size_t ind = 0; ind < nr; ++ind) {
      v(i++) = tbl[ind];
    }
  }
}

void Data::get_genotypes_recessive(const size_t snp, VectorView v) const
{
  size_t offset = _snp_offset * snp;
  size_t i = 0;
  for (size_t ind4 = 0; ind4 < _n_4; ++ind4) {
    const double *tbl = _RECESSIVE_GENOTYPE_TABLE[_g[offset++]];
    v(i++) = tbl[0];
    v(i++) = tbl[1];
    v(i++) = tbl[2];
    v(i++) = tbl[3];
  }
  // rest of the individuals
  if (i < n) {
    size_t nr = n - i;
    const double *tbl = _RECESSIVE_GENOTYPE_TABLE[_g[offset]];
    for (size_t ind = 0; ind < nr; ++ind) {
      v(i++) = tbl[ind];
    }
  }
}

// private:
void Data::read_fam(const std::string& file_fam)
{
  size_t i = 0;
  // create a map of FID IID orders for y and e files
  std::ifstream file(file_fam.c_str());

  while (file.good()) {
    std::string fid, iid, tmp;
    double tmp_y = NAN;

    file >> fid; // family id
    file >> iid; // individual id
    file >> tmp; // paternal id
    file >> tmp; // maternal id
    file >> tmp; // sex
    file >> tmp_y; // phenotype

    if (file.good() && !fid.empty() && !iid.empty() && !std::isnan(tmp_y)) {
      fam_mapping[make_pair(fid, iid)] = i;
      _y(i) = tmp_y;
      ++i;
    }
  }

  file.close();

  if (i != n || fam_mapping.size() != n)
    throw std::runtime_error("FAM file size does not match given n");
}

void Data::read_y(const std::string& file_y)
{
  size_t i = 0;

  std::ifstream file(file_y.c_str());

  while (file.good()) {
    std::string fid, iid;
    std::map<std::pair<std::string, std::string>, size_t>::const_iterator it;
    double tmp_y = NAN;

    file >> fid; // family id
    file >> iid; // individual id
    file >> tmp_y; // phenotype

    if (file.good() && !fid.empty() && !iid.empty() && !std::isnan(tmp_y)) {
      it = fam_mapping.find(make_pair(fid, iid));
      if (it != fam_mapping.end()) {
        _y(it->second) = tmp_y;
        ++i;
      }
    }
  }

  file.close();

  if (i != n)
    throw std::runtime_error(
        "Alternate phenotype file does not contain all phenotypes");
}

void Data::read_e(const std::string& file_e)
{
  size_t i = 0;

  std::ifstream file(file_e.c_str());

  while (file.good()) {
    std::string fid, iid;
    std::map<std::pair<std::string, std::string>, size_t>::const_iterator it;
    std::vector<double> tmp_y(m_e - 1, NAN); // could use Eigen::VectorXd

    file >> fid; // family id
    file >> iid; // individual id
    for (size_t j = 0; j < (m_e - 1); ++j)
      file >> tmp_y[j]; // covariates

    if (file.good() && !fid.empty() && !iid.empty()) {
      // check for nans (is this neccessary?)
      bool nans = 0;
      for (size_t j = 0; j < (m_e - 1); ++j) {
        if (std::isnan(tmp_y[j])) {
          nans = 1;
          break;
        }
      }
      if (nans)
        continue;

      it = fam_mapping.find(make_pair(fid, iid));
      if (it != fam_mapping.end()) {
        for (size_t j = 0; j < (m_e - 1); ++j)
          _e(it->second, j + 1) = tmp_y[j];
        ++i;
      }
    }
  }

  file.close();

  if (i != n)
    throw std::runtime_error("Covariate file does not contain all covariates");
}

void Data::read_g(const std::string& file_g)
{
  std::ifstream file(file_g.c_str(), std::ios::in | std::ios::binary);
  if (!file.good())
    throw std::runtime_error("BED file could not be opened");

  // magic number and snp-major: 01101100 00011011 00000001
  char byte = 0;
  file.read(&byte, 1);
  if (byte != 0x6C)
    throw std::runtime_error(
        "BED file not recognised (magic number does not match)");
  file.read(&byte, 1);
  if (byte != 0x1B)
    throw std::runtime_error(
        "BED file not recognised (magic number does not match)");
  file.read(&byte, 1);
  if (byte != 0x01)
    throw std::runtime_error("BED file not in snp-major format");

  // rest of the data
  size_t len = _snp_offset * m_g;
  _g = new unsigned char[len];
  file.read((char *) _g, len);

  if (!file.good() || (size_t)file.gcount() != len)
    throw std::runtime_error("Reading the BED file failed");
  file.close();
}

void Data::set_genotype(const size_t individual, const size_t snp,
                 const double genotype)
{
  if (!(genotype == -1.0 || genotype == 0.0 ||
        genotype == 1.0 || genotype == 2.0))
    throw std::domain_error("Genotype must be -1, 0, 1 or 2.");

  const int geno = static_cast<const int>(genotype);

  size_t offset = _snp_offset * snp + floor(individual / 4);
  unsigned char* byte = _g + offset;

  //(_g[offset] >> (2 * (individual % 4))) & 0x03

  unsigned char null_mask = 0x03;
  unsigned char mask = 0;
  switch (geno){
    case -1:
      mask = 0x01;
      break;
    case 0:
      mask = 0x00;
      break;
    case 1:
      mask = 0x02;
      break;
    case 2:
      mask = 0x03;
      break;
    default:
      throw std::domain_error("Genotype must be -1, 0, 1 or 2.");
  }

  null_mask = ~(null_mask << (2 * (individual % 4)));
  mask = mask << (2 * (individual % 4));

  // null_mask clears the bits and mask sets them to correct genotype
  *byte = (*byte & null_mask) | mask;
//
//  if (!(genotype == -1.0 || genotype == 0.0 ||
//        genotype == 1.0 || genotype == 2.0))
//    throw std::domain_error("Genotype must be -1, 0, 1 or 2.");
  /*case 0x01: return -1.0;
  case 0x00: return 0.0;
  case 0x02: return 1.0;
  case 0x03: return 2.0;*/

}

void Data::recode_g_to_minor_allele_count()
{
  Vector snp_genotypes(n);
  for (size_t snp = 0; snp < m_g; ++snp){
    get_genotypes_additive(snp, snp_genotypes);

    // need to recode if allele freq > 0.5, i.e. swap 0s to 2s and 2s to 0s
    if (Utils::compute_allele_freq(snp_genotypes) > 0.5){
      for (size_t ind = 0; ind < n; ++ind){
        double genotype = get_genotype(ind, snp);
        if (genotype == 0.0) set_genotype(ind, snp, 2.0);
        else if (genotype == 2.0) set_genotype(ind, snp, 0.0);
      }
    }
  }
}

void Data::handle_missing_g()
{
  miss_loc_ = new size_t*[m_g];
  miss_prior_ = new double*[m_g];

  for (size_t snp = 0; snp < m_g; ++snp) {
    miss_loc_[snp] = NULL;
    miss_prior_[snp] = NULL;

    int n_miss = 0;
    for (size_t ind = 0; ind < n; ++ind) {
      if (get_genotype(ind, snp) < 0)
        ++n_miss;
    }
    if (n_miss > 0) {
      miss_loc_[snp] = new size_t[n_miss + 1];
      miss_prior_[snp] = new double[3];
      miss_prior_[snp][0] = 0;
      miss_prior_[snp][1] = 0;
      miss_prior_[snp][2] = 0;

      miss_loc_[snp][0] = n_miss;
      n_miss = 0;
      for (size_t ind = 0; ind < n; ++ind) {
        double geno = get_genotype(ind, snp);
        if (geno < 0)
          miss_loc_[snp][++n_miss] = ind;
        else
          ++miss_prior_[snp][(int) geno];
      }
      // cumulative sum
      miss_prior_[snp][1] += miss_prior_[snp][0];
      miss_prior_[snp][2] += miss_prior_[snp][1];
    }
  }
}

double Data::compute_g_var() const
{
  double total_var = 0;
  size_t nsnps = 0;

  for (size_t snp = 0; snp < m_g; ++snp) {
    double squaresum_snp = 0;
    double sum_snp = 0;
    size_t ngenos = 0;
    for (size_t ind = 0; ind < n; ++ind) {
      double geno = get_genotype(ind, snp);
      if (geno >= 0) {
        squaresum_snp += geno * geno;
        sum_snp += geno;
        ++ngenos;
      }
    }
    if (ngenos > 1){
      total_var += (squaresum_snp - sum_snp * sum_snp / ngenos) / (ngenos - 1);
      ++nsnps;
    }
  }
  return total_var / nsnps;
}

void Data::compute_g_var_and_mean(double&  mean, double& var) const
{
  double total_var = 0;
  double total_mean = 0;
  size_t nsnps_mean = 0;
  size_t nsnps_var = 0;

  for (size_t snp = 0; snp < m_g; ++snp) {
    double squaresum_snp = 0;
    double sum_snp = 0;
    size_t ngenos = 0;
    for (size_t ind = 0; ind < n; ++ind) {
      double geno = get_genotype(ind, snp);
      if (geno >= 0) {
        squaresum_snp += geno * geno;
        sum_snp += geno;
        ++ngenos;
      }
    }
    if (ngenos > 1){
      total_mean += sum_snp / ngenos;
      total_var += (squaresum_snp - sum_snp * sum_snp / ngenos) / (ngenos - 1);
      ++nsnps_mean;
      ++nsnps_var;
    } else if (ngenos == 1){
      total_mean += sum_snp / ngenos;
      ++nsnps_mean;
    }
  }
  var = total_var / nsnps_var;
  mean = total_mean / nsnps_mean;
}

} // namespace bmagwa
