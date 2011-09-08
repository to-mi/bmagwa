/* BMAGWA software v1.0
 *
 * main.cpp
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


// these should be replaced automatically by release script
#define VERSION_MAJOR 1
#define VERSION_MINOR 0

#include <iostream>
#include <pthread.h>
#include "options.hpp"
#include "data.hpp"
#include "prior.hpp"
#include "sampler.hpp"

void* Run_sampler(void* threadarg);

int main(int argc, char* argv[])
{
  printf("-------------------------------------------------------------\n"
         "BMAGWA software version %d.%d\n"
         "http://www.lce.hut.fi/research/mm/bmagwa/\n"
         "-------------------------------------------------------------\n\n",
         VERSION_MAJOR, VERSION_MINOR);
  if (argc != 2){
    printf("Usage: %s INIFILE\n", argv[0]);
    return 0;
  }
  std::string options_file(argv[1]);
  bmagwa::Options options1(options_file);

  bmagwa::Options** options = new bmagwa::Options*[options1.n_threads];
  bmagwa::Sampler** samplers = new bmagwa::Sampler*[options1.n_threads];
  pthread_t* threads = new pthread_t[options1.n_threads];

  const bmagwa::Data* data = new bmagwa::Data(options1.n, options1.m_g,
                                              options1.m_e,
                                              options1.file_fam, options1.file_g,
                                              options1.recode_g_to_minor_allele_count,
                                              options1.file_e, options1.file_y);

  std::cout << "var y = " << data->var_y() << std::endl
            << "var x = " << data->var_x() << std::endl
            << "mean x = " << data->mean_x() << std::endl;

  std::cout << "Precomputing SNP covariances" << std::endl;
  bmagwa::DataModel* data_model = new bmagwa::DataModel(data, options1.types);
  const bmagwa::PrecomputedSNPCovariances* pre_xxcov =
                              new bmagwa::PrecomputedSNPCovariances(data_model);
  delete data_model;

  for (size_t t = 0; t < options1.n_threads; ++t){
    options[t] = options1.clone(t);

    std::cout << "Initializing sampler " << t << std::endl;

    samplers[t] = new bmagwa::Sampler(*options[t], data, pre_xxcov);
  }
  samplers[0]->print_prior();

  for(size_t t = 0; t < options1.n_threads; ++t){
    std::cout << "Creating thread " << t << std::endl;
    int rc = pthread_create(&threads[t], NULL, &Run_sampler, (void *)(samplers[t]));
    if (rc){
      throw std::runtime_error("Creating thread failed");
    }
  }

  for(size_t t = 0; t < options1.n_threads; ++t) {
    void* status;
    int rc = pthread_join(threads[t], &status);
    if (rc) {
      throw std::runtime_error("Thread join failed");
    }
    std::cout << "Completed join with thread " << t << " having a status of "
              << (long)status << std::endl;
  }

  for (size_t t = 0; t < options1.n_threads; ++t){
    delete samplers[t];
    delete options[t];
  }
  delete data;
  delete[] samplers;
  delete[] threads;
  delete[] options;
  delete pre_xxcov;

  return 0;
}


void* Run_sampler(void* threadarg){
  bmagwa::Sampler *sampler = reinterpret_cast<bmagwa::Sampler*>(threadarg);

  sampler->initialize_p_proposal();
  sampler->save_p_proposal();

  // run sampler
  sampler->sample();

  return NULL;
}
