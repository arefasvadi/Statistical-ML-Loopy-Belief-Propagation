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

static int total_colors = 0;
static int iterations = 0;
static int burn_in = 0;
static vector<double> color_weights;
static int vertex_size = 0;

using Edge_t = pair<int, int>;
using EdgeCompare_t = struct edge_compare {
 bool operator()(const Edge_t &lhs, const Edge_t &rhs) const {
  if (lhs.first == rhs.first) {
   return lhs.second < rhs.second;
  } else if (lhs.first == rhs.second) {
   return (lhs.second < rhs.first);
  } else if (lhs.second == rhs.first){
   return lhs.first < rhs.second;
  } else if (lhs.second == rhs.second){
   return lhs.first < rhs.first;
  }
  return lhs.first < rhs.first;
 }
};
using Probabilities_t = map<Edge_t, vector<double>, EdgeCompare_t>;
using Graph_t = set<Edge_t, EdgeCompare_t>;
using Sample_t = map<Edge_t, int, EdgeCompare_t>;
using Final_Table_t = vector<vector<vector<double >>>;
static auto graph = make_shared<Graph_t>();
static auto final_table = make_shared<Final_Table_t>();
static shared_ptr<vector<Sample_t> > samples;

void read_input(string graph_path,string weight_path) {
 ifstream graph_ifile(graph_path);
 if (!graph_ifile.is_open()) {
  cout << "file opening problem! " << graph_path << "\n";
  exit(1);
 }
 string line;
 int current_vertex = 0;
 while (getline(graph_ifile, line)) {
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
  vertex_size = inner_vertex;
 }
 graph_ifile.close();

 ifstream weights_ifile(weight_path);
 if (!weights_ifile.is_open()) {
  cout << "file opening problem! " << weight_path << "\n";
  exit(1);
 }
 getline(weights_ifile, line);
 stringstream strm(line);
 while (getline(strm, line, ',')) {
  total_colors++;
  color_weights.push_back(stoi(line));
 }
 weights_ifile.close();
 
 final_table->resize(vertex_size);
 for (int i = 0; i < vertex_size; ++i) {
  (*final_table)[i] = vector<vector<double >>(vertex_size);
  for (int j = 0; j < vertex_size; ++j) {
   (*final_table)[i][j] = vector<double>(total_colors, 0.0);
  }
 }
 
}

shared_ptr<Sample_t> find_suitable_sample() {
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

shared_ptr<vector<Sample_t> > gibbs(shared_ptr<Sample_t> first_sample) {
 
 samples = make_shared<vector<Sample_t> >();
 std::random_device r;
 std::default_random_engine generator{r()};
 std::uniform_real_distribution<double> distribution(0.0, 1.0);
 int total_rounds = burn_in + iterations;
 while (total_rounds > 0) {
  cout << "Iteration: " << (burn_in + iterations - total_rounds + 1) << "\n";
  for (auto edge_itr = first_sample->begin(); edge_itr != first_sample->end() && total_rounds > 0; ++edge_itr) {
   set<int> neghibour_colors;
   int sum_weight = 0;
   for (auto inner_edge_itr = first_sample->begin(); inner_edge_itr != first_sample->end(); ++inner_edge_itr) {
    if (inner_edge_itr != edge_itr) {
     if (edge_itr->first.first == inner_edge_itr->first.first
         || edge_itr->first.first == inner_edge_itr->first.second
         || edge_itr->first.second == inner_edge_itr->first.first
         || edge_itr->first.second == inner_edge_itr->first.second) {
      neghibour_colors.insert(inner_edge_itr->second);
      sum_weight += color_weights[inner_edge_itr->second - 1];
     }
    }
   }
   vector<double> probs(total_colors);
   for (int color = 0; color < total_colors; ++color) {
    if (neghibour_colors.find(color + 1) != neghibour_colors.end()) {
     probs[color] = 0.0;
    } else {
     probs[color] = exp(sum_weight + color_weights[color]);
    }
   }
   //Normalize
   double prob_sum = 0.0;
   for (const auto &val:probs) {
    prob_sum += val;
   }
   map<double_t, int> rev_prob;
   double cummulative_sum = 0.0;
   for (int color = 0; color < total_colors; ++color) {
    if (probs[color] == 0.0) {
     continue;
    }
    probs[color] = (double) probs[color] / prob_sum;
    cummulative_sum += probs[color];
    rev_prob[cummulative_sum] = color + 1;
   }
   //it's already sorted
   //divide the range 0 - 1
   double slot = distribution(generator);
   double prev_range = 0.0;
   for (const auto &val:rev_prob) {
    if (slot < val.first && slot >= prev_range) {
     cout << "\tupdated the edge from " << edge_itr->first.first << " to " << edge_itr->first.second
          << " from value: " << edge_itr->second << " to value: " << val.second << "\n";
     edge_itr->second = val.second;
     break;
    } else if (next(rev_prob.find(val.first), 1) == rev_prob.end()) {
     cout << "\tupdated the edge from " << edge_itr->first.first << " to " << edge_itr->first.second
          << " from value: " << edge_itr->second << " to value: " << val.second << "\n";
     edge_itr->second = val.second;
     break;
    }
    prev_range = val.first;
   }
   if (total_rounds <= iterations) {
    samples->push_back(*first_sample);
    //increment in table!
    for(const auto &edge:*first_sample) {
     (*final_table)[edge.first.first-1][edge.first.second-1][edge.second-1] += 1.0;
     (*final_table)[edge.first.second-1][edge.first.first-1][edge.second-1] += 1.0;
    }
   }
   --total_rounds;
  }
 }
 cout << samples->size() << " samples generated!\n";
 return samples;
}

void fill_print_table() {
 cout << "\n****\n";
 for (auto edge_itr = graph->begin(); edge_itr != graph->end(); ++edge_itr) {
//  Edge_t current_edge = *edge_itr;
  double prob_sum = 0.0;
  for(int color=0;color<total_colors;++color) {
   prob_sum += (*final_table)[edge_itr->first-1][edge_itr->second-1][color];
  }
  for(int color=0;color<total_colors;++color) {
   (*final_table)[edge_itr->first-1][edge_itr->second-1][color] =
    (*final_table)[edge_itr->first-1][edge_itr->second-1][color]/prob_sum;
   (*final_table)[edge_itr->second-1][edge_itr->first-1][color] =
    (*final_table)[edge_itr->first-1][edge_itr->second-1][color];
   cout << "\tedge between " << edge_itr->first << " and " << edge_itr->second
        << " with value " << color+1 << " has probability of: "
        << (*final_table)[edge_itr->first-1][edge_itr->second-1][color] << "\n";
  }
 }
}

int main(int argc, char **argv) {
 if (argc != 5){
  cout << "Correct usage requires 4 parameters:\n"
   "\t./gibbs-sampler \"graph_file\" \"weights_file\" \"num_burnins\" \"num_iterations\"";
  exit(1);
 }
 burn_in = stoi(string(argv[3]));
 iterations = stoi(string(argv[4]));
 read_input(string(argv[1]),string(argv[2]));
 auto first_sample = find_suitable_sample();
 samples = gibbs(first_sample);
 fill_print_table();
 return 0;
}