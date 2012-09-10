/* BMAGWA software v2.0
 *
 * discrete_distribution_tests.hpp
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

#ifndef DISCRETE_DISTRIBUTION_TESTS_HPP_
#define DISCRETE_DISTRIBUTION_TESTS_HPP_

#include "gtest/gtest.h"
#include <iostream>
#include "../discrete_distribution.hpp"


namespace {
  TEST(DiscreteDistributionTests, Construction) {
    bmagwa::Rand rng(1234, 10);
    size_t m_g = 5;
    double p[] = { 0.1, 0.2, 0.3, 0.4, 0.5 };

    bmagwa::DiscreteDistribution dd(true, p, m_g, rng);

    bmagwa::DDNode* nodes = dd.nodes();
    for (size_t i = 0; i < 5; ++i)
      ASSERT_EQ(i, nodes[i].index);

    ASSERT_TRUE(nodes[0].parent == NULL);
    ASSERT_TRUE(nodes[1].parent == nodes);
    ASSERT_TRUE(nodes[2].parent == nodes);
    ASSERT_TRUE(nodes[3].parent == nodes + 1);
    ASSERT_TRUE(nodes[4].parent == nodes + 1);

    ASSERT_TRUE(nodes[0].left_child == nodes + 1);
    ASSERT_TRUE(nodes[0].right_child == nodes + 2);
    ASSERT_TRUE(nodes[1].left_child == nodes + 3);
    ASSERT_TRUE(nodes[1].right_child == nodes + 4);
    for (size_t i = 2; i < 5; ++i){
      ASSERT_TRUE(nodes[i].left_child == NULL);
      ASSERT_TRUE(nodes[i].right_child == NULL);
    }

    ASSERT_TRUE(nodes[0].upright_parent == NULL);
    ASSERT_TRUE(nodes[1].upright_parent == nodes);
    ASSERT_TRUE(nodes[2].upright_parent == NULL);
    ASSERT_TRUE(nodes[3].upright_parent == nodes + 1);
    ASSERT_TRUE(nodes[4].upright_parent == nodes);

    ASSERT_TRUE(dd.rightmost_node() == nodes + 2);

    ASSERT_EQ(0.2 + 0.4, nodes[1].value);
    ASSERT_EQ(0.4, nodes[3].value);
    ASSERT_EQ(0.5, nodes[4].value);
    ASSERT_EQ(0.3, nodes[2].value);
    ASSERT_EQ(0.1 + 0.2 + 0.4 + 0.5, nodes[0].value);

    ASSERT_EQ(0.1 + 0.2 + 0.4 + 0.5 + 0.3, dd.total_w());
  }

  TEST(DiscreteDistributionTests, Updates) {
    double threshold = 1e-15;
    bmagwa::Rand rng(1234, 10);
    size_t m_g = 5;
    double p[] = { 0.1, 0.2, 0.3, 0.4, 0.5 };

    bmagwa::DiscreteDistribution dd(true, p, m_g, rng);
    bmagwa::DDNode* nodes = dd.nodes();

    dd.adddate(0);
    ASSERT_TRUE(dd.zeroed(0));

    ASSERT_NEAR(0.2 + 0.4, nodes[1].value, threshold);
    ASSERT_NEAR(0.4, nodes[3].value, threshold);
    ASSERT_NEAR(0.5, nodes[4].value, threshold);
    ASSERT_NEAR(0.3, nodes[2].value, threshold);
    ASSERT_NEAR(0.2 + 0.4 + 0.5, nodes[0].value, threshold);
    ASSERT_NEAR(0.2 + 0.4 + 0.5 + 0.3, dd.total_w(), threshold);

    dd.adddate(4);
    ASSERT_TRUE(dd.zeroed(4));

    ASSERT_NEAR(0.2 + 0.4, nodes[1].value, threshold);
    ASSERT_NEAR(0.4, nodes[3].value, threshold);
    ASSERT_NEAR(0.0, nodes[4].value, threshold);
    ASSERT_NEAR(0.3, nodes[2].value, threshold);
    ASSERT_NEAR(0.2 + 0.4, nodes[0].value, threshold);
    ASSERT_NEAR(0.2 + 0.4 + 0.3, dd.total_w(), threshold);

    dd.adddate(2);
    ASSERT_TRUE(dd.zeroed(2));

    ASSERT_NEAR(0.2 + 0.4, nodes[1].value, threshold);
    ASSERT_NEAR(0.4, nodes[3].value, threshold);
    ASSERT_NEAR(0.0, nodes[4].value, threshold);
    ASSERT_NEAR(0.0, nodes[2].value, threshold);
    ASSERT_NEAR(0.2 + 0.4, nodes[0].value, threshold);
    ASSERT_NEAR(0.2 + 0.4, dd.total_w(), threshold);

    dd.remdate(0);
    ASSERT_FALSE(dd.zeroed(0));

    ASSERT_NEAR(0.2 + 0.4, nodes[1].value, threshold);
    ASSERT_NEAR(0.4, nodes[3].value, threshold);
    ASSERT_NEAR(0.0, nodes[4].value, threshold);
    ASSERT_NEAR(0.0, nodes[2].value, threshold);
    ASSERT_NEAR(0.1 + 0.2 + 0.4, nodes[0].value, threshold);
    ASSERT_NEAR(0.1 + 0.2 + 0.4, dd.total_w(), threshold);

    double tw = dd.total_w();
    dd.switchdate(3, 2);
    dd.switchdate(2, 3);
    ASSERT_DOUBLE_EQ(tw, dd.total_w());
    dd.switchdate(3, 2);

    ASSERT_TRUE(dd.zeroed(3));
    ASSERT_FALSE(dd.zeroed(2));

    ASSERT_NEAR(0.2, nodes[1].value, threshold);
    ASSERT_NEAR(0.0, nodes[3].value, threshold);
    ASSERT_NEAR(0.0, nodes[4].value, threshold);
    ASSERT_NEAR(0.3, nodes[2].value, threshold);
    ASSERT_NEAR(0.1 + 0.2, nodes[0].value, threshold);
    ASSERT_NEAR(0.1 + 0.2 + 0.3, dd.total_w(), threshold);

    dd.remdate(4);
    dd.remdate(3);
    ASSERT_FALSE(dd.zeroed(4));
    ASSERT_FALSE(dd.zeroed(3));

    ASSERT_NEAR(0.2 + 0.4, nodes[1].value, threshold);
    ASSERT_NEAR(0.4, nodes[3].value, threshold);
    ASSERT_NEAR(0.5, nodes[4].value, threshold);
    ASSERT_NEAR(0.3, nodes[2].value, threshold);
    ASSERT_NEAR(0.1 + 0.2 + 0.4 + 0.5, nodes[0].value, threshold);
    ASSERT_NEAR(0.1 + 0.2 + 0.4 + 0.5 + 0.3, dd.total_w(), threshold);

    dd.adddate(0);
    dd.adddate(1);
    dd.adddate(2);
    dd.adddate(3);
    dd.adddate(4);
    ASSERT_TRUE(dd.zeroed(0));
    ASSERT_TRUE(dd.zeroed(1));
    ASSERT_TRUE(dd.zeroed(2));
    ASSERT_TRUE(dd.zeroed(3));
    ASSERT_TRUE(dd.zeroed(4));

    ASSERT_NEAR(0.0, nodes[1].value, threshold);
    ASSERT_NEAR(0.0, nodes[3].value, threshold);
    ASSERT_NEAR(0.0, nodes[4].value, threshold);
    ASSERT_NEAR(0.0, nodes[2].value, threshold);
    ASSERT_NEAR(0.0, nodes[0].value, threshold);
    ASSERT_NEAR(0.0, dd.total_w(), threshold);

    dd.remdate(0);
    dd.remdate(1);
    dd.remdate(2);
    dd.remdate(3);
    dd.remdate(4);
    ASSERT_FALSE(dd.zeroed(0));
    ASSERT_FALSE(dd.zeroed(1));
    ASSERT_FALSE(dd.zeroed(2));
    ASSERT_FALSE(dd.zeroed(3));
    ASSERT_FALSE(dd.zeroed(4));

    for (size_t i = 0; i < 1000000; ++i){
      size_t ind = (size_t)round(4 * rng.rand_01());
      dd.adddate(ind);
      ASSERT_TRUE(dd.total_w() >= 0);
      dd.adddate((ind + 2) % 5);
      ASSERT_TRUE(dd.total_w() >= 0);
      dd.switchdate((ind + 1) % 5, ind);
      ASSERT_TRUE(dd.total_w() >= 0);
      dd.remdate((ind + 2) % 5);
      ASSERT_TRUE(dd.total_w() >= 0);
      dd.remdate((ind + 1) % 5);
      ASSERT_TRUE(dd.total_w() >= 0);
    }

    ASSERT_NEAR(0.2 + 0.4, nodes[1].value, threshold);
    ASSERT_NEAR(0.4, nodes[3].value, threshold);
    ASSERT_NEAR(0.5, nodes[4].value, threshold);
    ASSERT_NEAR(0.3, nodes[2].value, threshold);
    ASSERT_NEAR(0.1 + 0.2 + 0.4 + 0.5, nodes[0].value, threshold);
    ASSERT_NEAR(0.1 + 0.2 + 0.4 + 0.5 + 0.3, dd.total_w(), threshold);
  }

  TEST(DiscreteDistributionTests, Sample) {
    size_t nsamples = 1000000;
    double threshold = 0.001;

    bmagwa::Rand rng(1234, 10);
    size_t m_g = 5;
    double p[] = { 0.1, 0.001, 0.3, 0.099, 0.5 };
    double vals[] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

    bmagwa::DiscreteDistribution dd(true, p, m_g, rng);

    for (size_t i = 0; i < nsamples; ++i)
      vals[dd.sample()] += 1.0;

    for (size_t i = 0; i < 5; ++i){
      vals[i] /= nsamples;
      EXPECT_NEAR(p[i], vals[i], threshold);
    }
  }

  TEST(DiscreteDistributionTests, Sample2) {
    size_t nsamples = 5000000;
    double threshold = 0.001;

    bmagwa::Rand rng(1234, 10);
    size_t m_g = 10;
    double p[10];
    double vals[10];

    double s = 0.0;
    for (int i = 0; i < 10; ++i){
      p[i] = rng.rand_01();
      s += p[i];
      vals[i] = 0.0;
    }
    for (int i = 0; i < 10; ++i){
      p[i] = p[i] / s;
    }

    bmagwa::DiscreteDistribution dd(true, p, m_g, rng);
    bmagwa::DDNode* nodes = dd.nodes();

    ASSERT_DOUBLE_EQ(p[1]+p[3]+p[7]+p[8],nodes[1].value);
    ASSERT_DOUBLE_EQ(p[3]+p[7],nodes[3].value);
    ASSERT_DOUBLE_EQ(p[7],nodes[7].value);
    ASSERT_DOUBLE_EQ(p[8],nodes[8].value);
    ASSERT_DOUBLE_EQ(p[9],nodes[9].value);
    ASSERT_DOUBLE_EQ(p[4]+p[9],nodes[4].value);
    ASSERT_DOUBLE_EQ(p[0]+p[1]+p[3]+p[4]+p[7]+p[8]+p[9],nodes[0].value);
    ASSERT_DOUBLE_EQ(1.0, dd.total_w());

    for (size_t i = 0; i < nsamples; ++i)
      vals[dd.sample()] += 1.0;

    for (size_t i = 0; i < 10; ++i){
      vals[i] /= nsamples;
      EXPECT_NEAR(p[i], vals[i], threshold);
    }
  }

  TEST(DiscreteDistributionTests, Sample3) {
    size_t nsamples = 5000000;
    double threshold = 0.001;

    bmagwa::Rand rng(1234, 10);
    size_t m_g = 14;
    double p[14];
    double vals[14];

    double s = 0.0;
    for (int i = 0; i < 14; ++i){
      p[i] = rng.rand_01();
      s += p[i];
      vals[i] = 0.0;
    }
    for (int i = 0; i < 14; ++i){
      p[i] = p[i] / s;
    }

    bmagwa::DiscreteDistribution dd(true, p, m_g, rng);

    ASSERT_DOUBLE_EQ(1.0, dd.total_w());

    for (size_t i = 0; i < nsamples; ++i)
      vals[dd.sample()] += 1.0;

    for (size_t i = 0; i < 14; ++i){
      vals[i] /= nsamples;
      EXPECT_NEAR(p[i], vals[i], threshold);
    }
  }

  TEST(DiscreteDistributionTests, SampleManySizes) {
    size_t nsamples = 5000000;
    double threshold = 0.001;

    bmagwa::Rand rng(1234, 10);

    double p[32];
    double vals[32];

    for (size_t size = 1; size < 32; ++size){
      size_t m_g = size;

      double s = 0.0;
      for (size_t i = 0; i < size; ++i){
        p[i] = rng.rand_01();
        s += p[i];
        vals[i] = 0.0;
      }
      for (size_t i = 0; i < size; ++i){
        p[i] = p[i] / s;
      }

      bmagwa::DiscreteDistribution dd(true, p, m_g, rng);

      ASSERT_DOUBLE_EQ(1.0, dd.total_w());

      for (size_t i = 0; i < nsamples; ++i)
        vals[dd.sample()] += 1.0;

      for (size_t i = 0; i < size; ++i){
        vals[i] /= nsamples;
        EXPECT_NEAR(p[i], vals[i], threshold);
      }
    }
  }

  TEST(DiscreteDistributionTests, SampleSomeZeroed) {
    size_t nsamples = 5000000;
    double threshold = 0.001;

    bmagwa::Rand rng(1234, 10);

    double p[16];
    double vals[16];

    for (size_t size = 1; size < 16; ++size){
      size_t m_g = size;

      double s = 0.0;
      for (size_t i = 0; i < size; ++i){
        p[i] = rng.rand_01();
        s += p[i];
        vals[i] = 0.0;
      }
      for (size_t i = 0; i < size; ++i){
        p[i] = p[i] / s;
      }

      bmagwa::DiscreteDistribution dd(true, p, m_g, rng);

      for (size_t i = 0; i < (size / 2); ++i){
        size_t ind;
        do {
          ind = (size_t)round((size-1) * rng.rand_01());
        } while(dd.zeroed(ind));
        dd.adddate(ind);
      }

      // renormalize...
      s = 0.0;
      for (size_t i = 0; i < size; ++i){
        if (dd.zeroed(i))
          p[i] = 0.0;
        else
          s += p[i];
      }
      for (size_t i = 0; i < size; ++i){
        p[i] = p[i] / s;
      }
      dd.update_weights(p);

      ASSERT_DOUBLE_EQ(1.0, dd.total_w());

      for (size_t i = 0; i < nsamples; ++i)
        vals[dd.sample()] += 1.0;

      for (size_t i = 0; i < size; ++i){
        vals[i] /= nsamples;
        EXPECT_NEAR(p[i], vals[i], threshold);
      }
    }
  }

}

#endif /* DISCRETE_DISTRIBUTION_TESTS_HPP_ */
