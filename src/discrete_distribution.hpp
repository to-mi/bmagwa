/* BMAGWA software v2.0
 *
 * discrete_distribution.hpp
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


#ifndef DISCRETE_DISTRIBUTION_HPP_
#define DISCRETE_DISTRIBUTION_HPP_

#include "rand.hpp"

namespace bmagwa {

//! Node (item) in DiscreteDistribution.
class DDNode
{
  public:
    DDNode()
    : parent(NULL),
      left_child(NULL),
      right_child(NULL),
      upright_parent(NULL),
      value(0.0),
      zeroed(false),
      index(0)
    {}

    ~DDNode() {}

    DDNode* parent;
    DDNode* left_child;
    DDNode* right_child;
    // go up towards the root in the tree; the first time you go up and right
    // you arrive at upright_parent
    DDNode* upright_parent;
    double value;
    bool zeroed;
    // snp index for the node
    size_t index;
};

//! Sampling from discrete distribution with given weights per item.
/*!
 *  Threaded binary tree implementation for sampling from a discrete
 *  distribution with sampling probabilities proportional to given weights.
 *  Provides means to zero and non-zero the weights of individual items.
 */
class DiscreteDistribution
{
  public:
    DiscreteDistribution(const bool use_external_weights,
                         double* weights, const size_t size, Rand& rng)
    : use_external_weights_(use_external_weights),
      w_(NULL),
      size_(size),
      total_w_(0.0),
      rng_(rng)
    {
      if (use_external_weights_){
        w_ = weights;
      } else {
        w_ = new double[size_];
        memcpy(w_, weights, size_ * sizeof(double));
      }
      nodes_ = new DDNode[size_];
      construct_tree();
    }

    ~DiscreteDistribution()
    {
      if (!use_external_weights_) delete[] w_;
      delete[] nodes_;
    }

    bool use_external_weights() const
    {
      return use_external_weights_;
    }

    void update_weights(double* weights){
      if (use_external_weights()){
        w_ = weights;
      } else {
        memcpy(w_, weights, size_ * sizeof(double));
      }
      update_tree_weights();
    }

    double weight(size_t ind) const {
      return w_[ind];
    }

    bool zeroed(size_t ind) const {
      return nodes_[ind].zeroed;
    }

    size_t size() const { return size_; }

    double* weights() { return w_; }

    DDNode* root() { return root_; }

    DDNode* rightmost_node() { return rightmost_node_; }

    DDNode* nodes() { return nodes_; }

    double total_w() const { return total_w_; }

    size_t sample() const
    {
      double r = rng_.rand_01() * total_w_;
      // avoid minus by carrying left subtree weight on selecting right subtree
      double carry = 0.0;

      DDNode* node = root_;
      while (node->left_child != NULL){
        if (r < (node->value + carry)){
          node = node->left_child;
        } else {
          carry += node->value;
          if (node->right_child != NULL){
            node = node->right_child;
          } else {
            assert(node->upright_parent->zeroed == false);
            return node->upright_parent->index;
          }
        }
      }

      if (r < (node->value + carry)){
        assert(node->zeroed == false);
        return node->index;
      } else {
        assert(node->upright_parent->zeroed == false);
        return node->upright_parent->index;
      }
    }

    // zero some value
    void adddate(const size_t ind)
    {
      DDNode* node = nodes_ + ind;
      node->zeroed = true;
      if (node->left_child != NULL){
        node->value = node->left_child->value;
        DDNode* tmpnode = node->left_child->right_child;
        while (tmpnode != NULL){
          node->value += tmpnode->value;
          tmpnode = tmpnode->right_child;
        }
      } else {
        node->value = 0.0;
      }

      while (node->upright_parent != NULL){
        node = node->upright_parent;
        node->value = (node->zeroed ? 0.0 : w_[node->index]) +
                      node->left_child->value;
        DDNode* tmpnode = node->left_child->right_child;
        while (tmpnode != NULL){
          node->value += tmpnode->value;
          tmpnode = tmpnode->right_child;
        }
       }

      compute_total_w();
    }

    // non-zero some value
    void remdate(const size_t ind)
    {
      double value = w_[ind];
      DDNode* node = nodes_ + ind;
      node->zeroed = false;
      node->value += value;

      while (node->upright_parent != NULL){
        node = node->upright_parent;
        node->value += value;
      }

      total_w_ += value;
    }

    void switchdate(const size_t ind_add, const size_t ind_rem){
      adddate(ind_add);
      remdate(ind_rem);
    }

  private:
    bool use_external_weights_;
    double* w_;
    size_t size_;
    double total_w_;
    Rand& rng_;

    DDNode* root_;
    DDNode* rightmost_node_;
    DDNode* nodes_;

    void construct_tree()
    {
      DDNode* node = NULL;

      // sort indices? would need to check all methods if changing nodes-order

      // build tree
      root_ = nodes_;
      root_->index = 0;
      root_->value = w_[root_->index];
      root_->zeroed = false;

      for (size_t ind = 1; ind < size_; ++ind){
        node = nodes_ + ind;
        node->index = ind;
        node->zeroed = false;

        // find parent ind and connect
        size_t pind = (ind - 1) / 2;
        node->parent = nodes_ + pind;
        if (ind % 2 == 0){
          // right child
          node->parent->right_child = node;
          // find upright index (if any)
          DDNode* tmp = node->parent;
          while (tmp->parent != NULL){
            if (tmp->parent->left_child == tmp){
              node->upright_parent = tmp->parent;
              break;
            }
            tmp = tmp->parent;
          }
        } else {
          // left child
          node->parent->left_child = node;
          node->upright_parent = node->parent;
        }
      }

      rightmost_node_ = root_;
      while (rightmost_node_->right_child != NULL)
        rightmost_node_ = rightmost_node_->right_child;

      update_tree_weights();
    }

    void update_tree_weights()
    {
      DDNode* node = root_;
      DDNode* tmp_node = NULL;
      bool went_upright = false;

      // depth first order
      while (node != rightmost_node_){
        if (went_upright){
          went_upright = false;
          // can't go left as we came up from there, so
          // go to right child OR upright
          if (node->right_child != NULL){
            node = node->right_child;
          } else {
            // going upright (here upright_parent == parent)
            tmp_node = node;
            node = node->upright_parent;
            while (tmp_node != node){
              node->value += tmp_node->value;
              tmp_node = tmp_node->parent;
            }
            went_upright = true;
          }
        } else {
          // visiting new node
          node->value = node->zeroed ? 0.0 : w_[node->index];
          // go to left child OR upright
          if (node->left_child != NULL){
            node = node->left_child;
          } else {
            // going upright (here upright_parent may or may not be == parent)
            tmp_node = node;
            node = node->upright_parent;
            while (tmp_node != node){
              node->value += tmp_node->value;
              tmp_node = tmp_node->parent;
            }
            went_upright = true;
          }
        }
      }
      // update rightmost node and its possible left child
      node->value = node->zeroed ? 0.0 : w_[node->index];
      if (node->left_child != NULL){
        node->left_child->value = node->left_child->zeroed ? 0.0 :
                                                   w_[node->left_child->index];
        node->value += node->left_child->value;
      }

      compute_total_w();
    }

    void compute_total_w()
    {
      DDNode* node = root_;
      total_w_ = node->value;
      while(node->right_child != NULL)
      {
        node = node->right_child;
        total_w_ += node->value;
      }
    }

    // disallow copy constructor and assignment
    DiscreteDistribution(const DiscreteDistribution& dd);
    DiscreteDistribution& operator=(const DiscreteDistribution& dd);
};

} // namespace bmagwa

#endif /* DISCRETE_DISTRIBUTION_HPP_ */
