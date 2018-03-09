#include <vector>
#include <algorithm>
#include<memory>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <utility>
#include <map>
#include <tuple>
#include <cmath>
#include <set>
#include <array>
#include <cassert>
#include <random>

using namespace std;

static int total_colors = 3;
static int iterations = 1000;
static int burn_in = 5000;
static vector<double> color_weights = {1, 2, 3};

using Edge_t = pair<int, int>;
using EdgeCompare_t = struct edge_compare {
 bool operator()(const Edge_t &lhs, const Edge_t &rhs) const {
  if (lhs.first == rhs.first) {
   return lhs.second < rhs.second;
  } else if (lhs.first == rhs.second) {
   return (lhs.second < rhs.first);
  }
  return lhs.first < rhs.first;
 }
};
using Probabilities_t = map<Edge_t, vector<double>, EdgeCompare_t>;
using Graph_t = set<Edge_t, EdgeCompare_t>;

static auto graph = make_shared<Graph_t>();

void read_input(string path) {
 ifstream ifile(path);
 if (!ifile.is_open()) {
  cout << "file opening problem! " << path << "\n";
  exit(1);
 }
 string line;
 int current_vertex = 0;
 while (getline(ifile, line)) {
  current_vertex++;
  stringstream strm(line);
  string val;
  int inner_vertex = 0;
  while (getline(strm, val, ',')) {
   inner_vertex++;
   if (stoi(val) == 1) {
    graph->insert(Edge_t{current_vertex, inner_vertex});
   }
  }
 }
 ifile.close();
}

shared_ptr<map<Edge_t, int, EdgeCompare_t> > find_suitable_sample() {
 vector<int> all_colors(total_colors, 0);
 for_each(all_colors.begin(), all_colors.end(),
          [](int &val) {
           static int i = 0;
           val = ++i;
          });
 bool perm_found = false;
 auto edge_coloring = make_shared<map<Edge_t, int, EdgeCompare_t> >();
 auto edge_itr = graph->begin();
 while (!perm_found) {
  //find this nodes neighbours
  vector<Graph_t::iterator> neighbours;
  for (Graph_t::iterator inner_edge_itr = graph->begin(); inner_edge_itr != graph->end(); ++inner_edge_itr) {
   if (!(inner_edge_itr->first == edge_itr->first && inner_edge_itr->second == edge_itr->second)) {
    if (edge_itr->first == inner_edge_itr->first || edge_itr->first == inner_edge_itr->second
        || edge_itr->second == inner_edge_itr->first || edge_itr->second == inner_edge_itr->second) {
     neighbours.push_back(inner_edge_itr);
    }
   }
  }
  vector<int> temp_colors(all_colors.size());
  copy(all_colors.begin(), all_colors.end(), temp_colors.begin());
  //this node and neighbours should be assigned unique values
  for (const auto &nbr_itr:neighbours) {
   Edge_t temp_edge{nbr_itr->first, nbr_itr->second};
   if (edge_coloring->find(temp_edge) != edge_coloring->end()) {
    temp_colors.erase(remove(temp_colors.begin(), temp_colors.end(), edge_coloring->at(temp_edge)), temp_colors.end());
   }
  }
  (*edge_coloring)[*edge_itr] = temp_colors.back();
  if (edge_coloring->size() == graph->size()) {
   perm_found = true;
   break;
  }
  temp_colors.pop_back();
  ++edge_itr;
 }
 if (perm_found) {
  cout << "First sample found:\n";
  for (const auto &edge:*edge_coloring) {
   cout << edge.first.first << " to " << edge.first.second << " is: " << edge.second << "\n";
  }
 }
 return edge_coloring;
}

void gibbs(shared_ptr<map<Edge_t, int, EdgeCompare_t> > first_sample) {

 std::default_random_engine generator;
 std::uniform_real_distribution<double> distribution(0.0,1.0);
 int total_rounds = burn_in + iterations;
 while (total_rounds > 0) {
  for (auto edge_val_itr = first_sample->begin(); edge_val_itr != first_sample->end();
       ++edge_val_itr) {
   //vector<map<Edge_t, int, EdgeCompare_t>::iterator> neghibours;
   int sum_weight = 0;
   for (auto inner_edge_val_itr = first_sample->begin(); inner_edge_val_itr != first_sample->end();
        ++inner_edge_val_itr) {
    if (inner_edge_val_itr != edge_val_itr) {
     if(edge_val_itr->first.first == inner_edge_val_itr->first.first
        || edge_val_itr->first.first == inner_edge_val_itr->first.second
        || edge_val_itr->first.second == inner_edge_val_itr->first.first
        || edge_val_itr->first.second == inner_edge_val_itr->first.second){
      //neghibours.push_back(inner_edge_val_itr);
      sum_weight += color_weights[inner_edge_val_itr->second - 1];
     }
    }
    vector<double> probs(total_colors);
    for(int color=0; color< total_colors;++color) {
     probs[color] = exp(sum_weight+color_weights[color]);
    }
    //Normalize
    double prob_sum = 0.0;
    for(const auto &val:probs){
     prob_sum+=val;
    }
    map<double_t ,int> rev_prob;
    double cummulative_sum = 0.0;
    for(int color = 0;color<total_colors;++total_colors){
     probs[color] = (double)probs[color]/prob_sum;
     cummulative_sum+=probs[color];
     rev_prob[cummulative_sum] = color;
    }
    //it's already sorted
    //divide the range 0-1
    double slot = distribution(generator);
    double prev_range = 0.0;
    for(const auto &val:rev_prob){
     if ( slot < val.first && slot >= prev_range){
      cout << "updated the edge from " << edge_val_itr->first.first << " to " << edge_val_itr->first.second
           << " from value: " <<  edge_val_itr->second << " to value: " << val.second << "\n";
      edge_val_itr->second = val.second;
      break;
     }
     else if (next(rev_prob.find(val.first),1) ==  rev_prob.end()){
      cout << "updated the edge from " << edge_val_itr->first.first << " to " << edge_val_itr->first.second
           << " from value: " <<  edge_val_itr->second << " to value: " << val.second << "\n";
      edge_val_itr->second = val.second;
      break;
     }
     prev_range = val.first;
    }
   }
  }
  if(total_rounds <= iterations) {

  }
  --total_rounds;
 }
}

int main(int argc, char **argv) {
 read_input("../graph2.txt");
 shared_ptr<map<Edge_t, int, EdgeCompare_t> > first_sample = find_suitable_sample();
 gibbs(first_sample);
 return 0;
}

