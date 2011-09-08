/* BMAGWA software v1.0
 *
 * options.hpp
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


#ifndef OPTIONS_HPP_
#define OPTIONS_HPP_

#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
//#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "inih/cpp/INIReader.h"
#include "data_model.hpp"
#include "linalg.hpp"

namespace bmagwa {

class Options
{
  public:
    std::string file_fam, file_g, file_e, file_y;
    size_t n, m_g, m_e;
    bool recode_g_to_minor_allele_count;

    std::string types_str;
    std::vector<DataModel::ef_t> types;

    size_t do_n_iter, n_rao, n_rao_burnin, thin, n_sample_tau2_and_missing;
    double t_proposal;
    bool adaptation, save_beta;
    size_t verbosity;

    std::string basename;
    size_t n_threads;
    std::vector<uint32_t> seeds;
    uint32_t seed;

    double e_qg, var_qg;
    double nu_sigma2, s2_sigma2, R2mode_sigma2;
    double nu_tau2[4], s2_tau2[4], eh_tau2[4];
    double inv_tau2_e_const_val, inv_tau2_e_val;
    Vector* inv_tau2_e;
    double types_prior[5];


    // constructors and destructors
    Options(const std::string options_file)
    : inv_tau2_e(NULL)
    {
      INIReader reader(options_file);

      if (reader.ParseError() != 0)
        throw std::runtime_error("Cannot load/parse configuration file.");

      parse_options(reader);
    }

    ~Options() {
      delete inv_tau2_e;
    }

    Options* clone(size_t index) const
    {
      Options* tmp = new Options(*this, index);
      return tmp;
    }

  private:
    // disallow assignment
    Options& operator=(const Options& options);
    Options(const Options& opt, size_t index)
    : file_fam(opt.file_fam), file_g(opt.file_g), file_e(opt.file_e),
      file_y(opt.file_y), n(opt.n), m_g(opt.m_g), m_e(opt.m_e),
      types(opt.types), do_n_iter(opt.do_n_iter), n_rao(opt.n_rao),
      n_rao_burnin(opt.n_rao_burnin), thin(opt.thin),
      n_sample_tau2_and_missing(opt.n_sample_tau2_and_missing),
      t_proposal(opt.t_proposal), adaptation(opt.adaptation),
      save_beta(opt.save_beta), verbosity(opt.verbosity),
      basename(opt.basename + boost::lexical_cast<std::string>(index)),
      n_threads(opt.n_threads), seeds(opt.seeds),
      seed(seeds[index]),
      e_qg(opt.e_qg), var_qg(opt.var_qg), nu_sigma2(opt.nu_sigma2),
      s2_sigma2(opt.s2_sigma2),
      R2mode_sigma2(opt.R2mode_sigma2),
      inv_tau2_e_const_val(opt.inv_tau2_e_const_val),
      inv_tau2_e_val(opt.inv_tau2_e_val)
    {
      inv_tau2_e = new Vector(*(opt.inv_tau2_e));
      memcpy(types_prior, opt.types_prior, 5 * sizeof(double));
      memcpy(nu_tau2, opt.nu_tau2, 4 * sizeof(double));
      memcpy(eh_tau2, opt.eh_tau2, 4 * sizeof(double));
      memcpy(s2_tau2, opt.s2_tau2, 4 * sizeof(double));
    }

    void parse_options(INIReader& r)
    {
      n = parse_integer_pos(r, "sizes", "n", false, 0);
      m_g = parse_integer_pos(r, "sizes", "m_g", false, 0);
      m_e = parse_integer_nonneg(r, "sizes", "m_e", true, 0);

      file_fam = parse_string(r, "datafiles", "file_fam", false, "");
      file_g = parse_string(r, "datafiles", "file_g", false, "");
      if (m_e > 0)
        file_e = parse_string(r, "datafiles", "file_e", false, "");
      else
    	file_e = parse_string(r, "datafiles", "file_e", true, "");
      file_y = parse_string(r, "datafiles", "file_y", true, "");
      recode_g_to_minor_allele_count = parse_bool(r, "datafiles",
          "recode_g_to_minor_allele_count", true, 0);

      types_str = parse_string(r, "model", "types", false, "");

      do_n_iter = parse_integer_pos(r, "sampler", "do_n_iter", false, 0);
      n_rao = parse_integer_nonneg(r, "sampler", "n_rao", false, 0);
      n_rao_burnin =
    		       parse_integer_nonneg(r, "sampler", "n_rao_burnin", false, 0);
      thin = parse_integer_pos(r, "sampler", "thin", false, 0);
      n_sample_tau2_and_missing =
    	 parse_integer_pos(r, "sampler", "n_sample_tau2_and_missing", false, 0);
      t_proposal = parse_double_01(r, "sampler", "t_proposal", false, 0.0);
      adaptation = parse_bool(r, "sampler", "adaptation", true, 0);
      save_beta = parse_bool(r, "sampler", "save_beta", true, 0);
      verbosity = parse_integer_nonneg(r, "sampler", "verbosity", false, 0);

      basename = parse_string(r, "thread", "basename", true, "chain");
      n_threads = parse_integer_pos(r, "thread", "n_threads", true, 1);
      std::string seeds_str = parse_string(r, "thread", "seeds", false, "");

      e_qg = parse_double_pos(r, "prior", "e_qg", false, 0.0);
      var_qg = parse_double_pos(r, "prior", "var_qg", false, 0.0);
      nu_sigma2 = parse_double_pos(r, "prior", "nu_sigma2", false, 0.0);
      R2mode_sigma2 = parse_double_01(r, "prior", "R2mode_sigma2", true, 0.0);
      s2_sigma2 = parse_double_nonneg(r, "prior", "s2_sigma2", true, 0.0);
      if (s2_sigma2 <= 0 && R2mode_sigma2 <= 0)
        throw std::runtime_error(
            "Config error: s2_sigma2 or R2mode_sigma2 must be > 0.");


      inv_tau2_e_const_val = parse_double_nonneg(r, "prior",
    		                                     "inv_tau2_e_const_val", false,
    		                                     0.0);
      if (m_e > 0)
        inv_tau2_e_val = parse_double_nonneg(r, "prior", "inv_tau2_e_val", false,
        		                             0.0);
      else
    	inv_tau2_e_val = parse_double_nonneg(r, "prior", "inv_tau2_e_val", true,
    			                             0.0);

      types_prior[DataModel::A] = parse_double(r, "prior", "type_A", true, 1.0);
      types_prior[DataModel::H] = parse_double(r, "prior", "type_H", true, 1.0);
      types_prior[DataModel::R] = parse_double(r, "prior", "type_R", true, 1.0);
      types_prior[DataModel::D] = parse_double(r, "prior", "type_D", true, 1.0);
      types_prior[DataModel::AH] = parse_double(r, "prior", "type_AH", true,
    		                                    1.0);

      // additional checks
      if (n_rao % thin != 0)
        throw std::runtime_error(
        		            "Config error: n_rao should be a multiple of thin");

      // parse types
      std::vector<std::string> types_sv;
      boost::split(types_sv, types_str, boost::is_any_of(","));
      for (std::vector<std::string>::iterator it = types_sv.begin();
           it != types_sv.end(); ++it){
        boost::trim(*it);
        if (*it == "A")
          types.push_back(DataModel::A);
        else if (*it == "H")
          types.push_back(DataModel::H);
        else if (*it == "D")
          types.push_back(DataModel::D);
        else if (*it == "R")
          types.push_back(DataModel::R);
        else if (*it == "AH")
          types.push_back(DataModel::AH);
        else {
          throw std::runtime_error("Config error: Unknown model type");
        }
      }
      std::sort(types.begin(), types.end());
      if (types.empty())
        throw std::runtime_error("Config error: Types cannot be empty");

      if (std::find(types.begin(), types.end(), DataModel::A) != types.end() ||
          std::find(types.begin(), types.end(), DataModel::AH) != types.end()){
        nu_tau2[DataModel::A] = parse_double_pos(r, "prior", "nu_tau2_A", false, 0.0);
        eh_tau2[DataModel::A] = parse_double_01(r, "prior", "eh_tau2_A", true, 0.0);
        s2_tau2[DataModel::A] = parse_double_nonneg(r, "prior", "s2_tau2_A", true, 0.0);
        if (s2_tau2[DataModel::A] <= 0 && eh_tau2[DataModel::A] <= 0)
          throw std::runtime_error(
              "Config error: s2_tau2_A or eh_tau2_A must be > 0.");
      }
      if (std::find(types.begin(), types.end(), DataModel::H) != types.end() ||
          std::find(types.begin(), types.end(), DataModel::AH) != types.end()){
        nu_tau2[DataModel::H] = parse_double_pos(r, "prior", "nu_tau2_H", false, 0.0);
        eh_tau2[DataModel::H] = parse_double_01(r, "prior", "eh_tau2_H", true, 0.0);
        s2_tau2[DataModel::H] = parse_double_nonneg(r, "prior", "s2_tau2_H", true, 0.0);
        if (s2_tau2[DataModel::H] <= 0 && eh_tau2[DataModel::H] <= 0)
          throw std::runtime_error(
              "Config error: s2_tau2_H or eh_tau2_H must be > 0.");
      }
      if (std::find(types.begin(), types.end(), DataModel::D) != types.end()){
        nu_tau2[DataModel::D] = parse_double_pos(r, "prior", "nu_tau2_D", false, 0.0);
        eh_tau2[DataModel::D] = parse_double_01(r, "prior", "eh_tau2_D", true, 0.0);
        s2_tau2[DataModel::D] = parse_double_nonneg(r, "prior", "s2_tau2_D", true, 0.0);
        if (s2_tau2[DataModel::D] <= 0 && eh_tau2[DataModel::D] <= 0)
          throw std::runtime_error(
              "Config error: s2_tau2_D or eh_tau2_D must be > 0.");
      }
      if (std::find(types.begin(), types.end(), DataModel::R) != types.end()){
        nu_tau2[DataModel::R] = parse_double_pos(r, "prior", "nu_tau2_R", false, 0.0);
        eh_tau2[DataModel::R] = parse_double_01(r, "prior", "eh_tau2_R", true, 0.0);
        s2_tau2[DataModel::R] = parse_double_nonneg(r, "prior", "s2_tau2_R", true, 0.0);
        if (s2_tau2[DataModel::R] <= 0 && eh_tau2[DataModel::R] <= 0)
          throw std::runtime_error(
              "Config error: s2_tau2_R or eh_tau2_R must be > 0.");
      }

      // parse seeds
      std::vector<std::string> seeds_sv;
      boost::split(seeds_sv, seeds_str, boost::is_any_of(","));
      for (std::vector<std::string>::iterator it = seeds_sv.begin();
           it != seeds_sv.end(); ++it){
        boost::trim(*it);
        seeds.push_back(str_to<uint32_t>(*it));
      }
      if (seeds.size() != n_threads){
        throw std::runtime_error(
                          "Config error: Number of seeds must equal n_threads");
      }
      for (size_t i = 0; i < (n_threads - 1); ++i){
        for (size_t j = (i + 1); j < n_threads; ++j){
          if (seeds[i] == seeds[j])
            throw std::runtime_error("Seeds are not unique");
        }
      }

      seed = seeds[0];

      // do inv_tau2_e
      inv_tau2_e = new Vector(m_e + 1);
      *inv_tau2_e = inv_tau2_e_val;
      (*inv_tau2_e)(0) = inv_tau2_e_const_val;
    }


    std::string parse_string(INIReader& reader, const std::string section,
    		                 const std::string name, const bool allow_default,
    		                 const std::string default_value)
    {
      std::string res(reader.Get(section, name, "DEFAULT_VALUE"));
      bool no_value = (res.empty() || res == "DEFAULT_VALUE");
      if (no_value && !allow_default)
    	throw std::runtime_error("Cannot leave option empty.");
      if (no_value){
    	  std::cout << "Warning: using default value of " << default_value
    	            << " for option " << section << "." << name << std::endl;
    	  return default_value;
      }
      return res;
    }

    double parse_double(INIReader& reader, const std::string section,
                        const std::string name, const bool allow_default,
                        const double default_value)
    {
      std::string res(reader.Get(section, name, "DEFAULT_VALUE"));
      bool no_value = (res.empty() || res == "DEFAULT_VALUE");
      if (no_value && !allow_default)
        throw std::runtime_error("Cannot leave option empty.");
      if (no_value){
        std::cout << "Warning: using default value of " << default_value
                  << " for option " << section << "." << name << std::endl;
        return default_value;
      }
      return str_to<double>(res);
    }

    double parse_double_01(INIReader& reader, const std::string section,
                           const std::string name, const bool allow_default,
                           const double default_value)
    {
    	double res = parse_double(reader, section, name, allow_default,
    			                  default_value);
    	if (res < 0 || res > 1)
          throw std::runtime_error("Value in 0..1 required");
    	return res;
    }

    double parse_double_pos(INIReader& reader, const std::string section,
                            const std::string name, const bool allow_default,
                            const double default_value)
    {
      double res = parse_double(reader, section, name, allow_default,
                            default_value);
      if (res <= 0)
        throw std::runtime_error("Positive value required for " + section + "."
        		                 + name);
      return res;
    }

    double parse_double_nonneg(INIReader& reader, const std::string section,
                               const std::string name, const bool allow_default,
                               const double default_value)
    {
      double res = parse_double(reader, section, name, allow_default,
                            default_value);
      if (res < 0)
        throw std::runtime_error("Non-negative value required for " + section +
        		                 "." + name);
      return res;
    }

    int parse_integer(INIReader& reader, const std::string section,
                      const std::string name, const bool allow_default,
                      const int default_value)
    {
      std::string res(reader.Get(section, name, "DEFAULT_VALUE"));
      bool no_value = (res.empty() || res == "DEFAULT_VALUE");
      if (no_value && !allow_default)
        throw std::runtime_error("Cannot leave option empty.");
      if (no_value){
        std::cout << "Warning: using default value of " << default_value
                  << " for option " << section << "." << name << std::endl;
        return default_value;
      }
      return str_to<int>(res);
    }

    int parse_integer_nonneg(INIReader& reader, const std::string section,
                             const std::string name, const bool allow_default,
                             const int default_value)
    {
      int res = parse_integer(reader, section, name, allow_default,
    		                  default_value);
      if (res < 0)
    	  throw std::runtime_error("Nonnegative integer required");
      return res;
    }

    int parse_integer_pos(INIReader& reader, const std::string section,
                          const std::string name, const bool allow_default,
                          const int default_value)
    {
      int res = parse_integer(reader, section, name, allow_default,
                          default_value);
      if (res <= 0)
        throw std::runtime_error("Positive integer required");
      return res;
    }

    bool parse_bool(INIReader& reader, const std::string section,
                    const std::string name, const bool allow_default,
                    const int default_value)
    {
      int res = parse_integer(reader, section, name, allow_default,
    		                  default_value);
      if (res != 0 && res != 1)
    	throw std::runtime_error("Boolean value required (0 or 1)");
      return (res == 1);
    }

    template<typename T>
    T str_to(const std::string& s, bool disallow_leftover_chars = true)
    {
      T x;
      convert(s, x, disallow_leftover_chars);
      return x;
    }

    template<typename T>
    void convert(const std::string& s, T& x,
    		     bool disallow_leftover_chars = true)
    {
      std::istringstream i(s);
      char c;
      if (!(i >> x) || (disallow_leftover_chars && i.get(c)))
        throw std::runtime_error("Invalid conversion.");
    }
};

} // namespace bmagwa

#endif /* OPTIONS_HPP_ */
