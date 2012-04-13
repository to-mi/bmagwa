/* BMAGWA software v2.0
 *
 * samplerstats.hpp
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


#ifndef SAMPLERSTATS_HPP_
#define SAMPLERSTATS_HPP_

#include <ctime>
#include <fstream>
#include "utils.hpp"

namespace bmagwa {

//! Book-keeping
class SamplerStats
{
  public:
    SamplerStats()
    : n_likelihood_updates_on_add_(0),
      n_likelihood_updates_on_rem_(0),
      n_likelihood_computations_(0),
      sampling_start_time_(time(NULL)),
      mhmove_elapsed(0),
      rao_elapsed(0)
    {}

    void reset_timer()
    {
      time(&sampling_start_time_);
    }

    double elapsed_time() const
    {
      return difftime(time(NULL), sampling_start_time_);
    }

    void mhmove_tic()
    {
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &mhmove_state);
    }

    void mhmove_toc()
    {
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &mhmove_state2);
      mhmove_elapsed += (mhmove_state2.tv_sec - mhmove_state.tv_sec) + (double)(mhmove_state2.tv_nsec - mhmove_state.tv_nsec) / 1000000000.0;
    }

    void rao_tic()
    {
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &rao_state);
    }

    void rao_toc()
    {
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &rao_state2);
      rao_elapsed += (rao_state2.tv_sec - rao_state.tv_sec) + (double)(rao_state2.tv_nsec - rao_state.tv_nsec) / 1000000000.0;
    }

    void likelihood_update_on_add()
    {
      ++n_likelihood_updates_on_add_;
    }

    void likelihood_update_on_rem()
    {
      ++n_likelihood_updates_on_rem_;
    }

    void likelihood_computation()
    {
      ++n_likelihood_computations_;
    }

    void print(std::string filename, const double p_ms) const
    {
      struct timespec clock_res;
      clock_getres(CLOCK_THREAD_CPUTIME_ID, &clock_res);
      std::ofstream* f_p = Utils::openfile(filename, false);

      *f_p << "n_likelihood_updates_on_add " << n_likelihood_updates_on_add_ << std::endl
           << "n_likelihood_updates_on_rem " << n_likelihood_updates_on_rem_ << std::endl
           << "n_likelihood_computations " << n_likelihood_computations_ << std::endl
           << "sampling_time_in_seconds " << elapsed_time() << std::endl
           << "mhmove_time_in_seconds " << mhmove_elapsed << std::endl
           << "rao_time_in_seconds " << rao_elapsed << std::endl
           << "clock_resolution_in_seconds " << clock_res.tv_sec + (double)clock_res.tv_nsec / 1000000000.0 << std::endl
           << "p_movesize " << p_ms << std::endl;


      f_p->close();
      delete f_p;
    }
  private:
    unsigned long n_likelihood_updates_on_add_;
    unsigned long n_likelihood_updates_on_rem_;
    unsigned long n_likelihood_computations_;
    time_t sampling_start_time_;
    struct timespec mhmove_state;
    struct timespec mhmove_state2;
    double mhmove_elapsed;
    struct timespec rao_state;
    struct timespec rao_state2;
    double rao_elapsed;

    // disallow copy constructor and assignment
    SamplerStats(const SamplerStats& ss);
    SamplerStats& operator=(const SamplerStats& ss);

};

} // namespace bmagwa

#endif /* SAMPLERSTATS_HPP_ */
