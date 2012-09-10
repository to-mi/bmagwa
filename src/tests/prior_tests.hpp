/* BMAGWA software v2.0
 *
 * prior_tests.hpp
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

#ifndef PRIOR_TESTS_HPP_
#define PRIOR_TESTS_HPP_

#include "gtest/gtest.h"
#include "../prior.hpp"


namespace {

  TEST(PriorTests, SimpleAddRemove)
  {
    double threshold = 0.000001;

    bmagwa::Data data(5, 10, 0, "testdata/plinktest.fam",
                      "testdata/plinktest.bed", false, "testdata/test.e");
    std::vector<bmagwa::DataModel::ef_t> types;
    types.push_back(bmagwa::DataModel::A);
    types.push_back(bmagwa::DataModel::AH);
    double types_prior[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
    double nu_tau2[4] = { 3.0, 3.0, 3.0, 3.0 };
    double s2_tau2[4] = { 1.0, 1.0, 1.0, 1.0 };

    bmagwa::DataModel data_model(&data, types);
    bmagwa::Vector inv_tau2_e(1); inv_tau2_e = 0;

    bmagwa::Prior prior(&data, &data_model,
          types_prior,
          3, 9,
          inv_tau2_e,
          1, 1,
          nu_tau2, s2_tau2, 1.0, false);

    //prior.print("testi.prior");

    // A type
    {
      int Ns[5] = { 0, 0, 0, 0, 0 };
      int Ns_new[5] = { 1, 0, 0, 0, 0 };
      EXPECT_NEAR(prior.compute_log_model(Ns_new) - prior.compute_log_model(Ns),
                  prior.compute_log_change_on_add(Ns, 0, bmagwa::DataModel::A),
                  threshold);
      EXPECT_NEAR(prior.compute_log_model(Ns) - prior.compute_log_model(Ns_new),
                  prior.compute_log_change_on_rem(Ns_new, 1, bmagwa::DataModel::A),
                  threshold);
    }

    {
      int Ns[5] = { 2, 0, 0, 0, 3 };
      int Ns_new[5] = { 3, 0, 0, 0, 3 };
      EXPECT_NEAR(prior.compute_log_model(Ns_new) - prior.compute_log_model(Ns),
                  prior.compute_log_change_on_add(Ns, 5, bmagwa::DataModel::A),
                  threshold);
      EXPECT_NEAR(prior.compute_log_model(Ns) - prior.compute_log_model(Ns_new),
                  prior.compute_log_change_on_rem(Ns_new, 6, bmagwa::DataModel::A),
                  threshold);
    }

    // AH type
    {
      int Ns[5] = { 0, 0, 0, 0, 0 };
      int Ns_new[5] = { 0, 0, 0, 0, 1 };
      EXPECT_NEAR(prior.compute_log_model(Ns_new) - prior.compute_log_model(Ns),
                  prior.compute_log_change_on_add(Ns, 0, bmagwa::DataModel::AH),
                  threshold);
      EXPECT_NEAR(prior.compute_log_model(Ns) - prior.compute_log_model(Ns_new),
                  prior.compute_log_change_on_rem(Ns_new, 1, bmagwa::DataModel::AH),
                  threshold);
    }

    {
      int Ns[5] = { 3, 0, 0, 0, 2 };
      int Ns_new[5] = { 3, 0, 0, 0, 3 };
      EXPECT_NEAR(prior.compute_log_model(Ns_new) - prior.compute_log_model(Ns),
                  prior.compute_log_change_on_add(Ns, 5, bmagwa::DataModel::AH),
                  threshold);
      EXPECT_NEAR(prior.compute_log_model(Ns) - prior.compute_log_model(Ns_new),
                  prior.compute_log_change_on_rem(Ns_new, 6, bmagwa::DataModel::AH),
                  threshold);
    }

    // switches
    {
      int Ns[5] = { 3, 0, 0, 0, 3 };
      int Ns_new[5] = { 3, 0, 0, 0, 3 };
      EXPECT_NEAR(prior.compute_log_model(Ns_new) - prior.compute_log_model(Ns),
                  prior.compute_log_change_on_swi(Ns, bmagwa::DataModel::AH,
                                                  bmagwa::DataModel::AH),
                  threshold);
    }

    {
      int Ns[5] = { 3, 0, 0, 0, 2 };
      int Ns_new[5] = { 2, 0, 0, 0, 3 };
      EXPECT_NEAR(prior.compute_log_model(Ns_new) - prior.compute_log_model(Ns),
                  prior.compute_log_change_on_swi(Ns, bmagwa::DataModel::AH,
                                                  bmagwa::DataModel::A),
                  threshold);
      EXPECT_NEAR(prior.compute_log_model(Ns) - prior.compute_log_model(Ns_new),
                  prior.compute_log_change_on_swi(Ns_new, bmagwa::DataModel::A,
                                                  bmagwa::DataModel::AH),
                  threshold);
    }
  }

}


#endif /* PRIOR_TESTS_HPP_ */
