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

using namespace std;

static constexpr int total_colors = 3;
static constexpr int iterations = 1000;
static constexpr int burn_in = 5000;
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
/* do {
  auto edge_itr = graph->begin();
  vector<int> temp_colors(all_colors.size());
  copy(all_colors.begin(), all_colors.end(), temp_colors.begin());
  while (edge_itr != graph->end()) {
   int current_color = all_colors.back();
   bool neighbour_has_it = false;
   auto edge_coloring_itr = edge_coloring->begin();
   while (edge_coloring_itr != edge_coloring->end()) {
    if (edge_coloring_itr->first.first != edge_itr->first
        && edge_coloring_itr->first.second != edge_itr->second) { // check to see it's not the same
     if (edge_itr->first == edge_coloring_itr->first.first
         || edge_itr->first == edge_coloring_itr->first.second
         || edge_itr->second == edge_coloring_itr->first.first
         || edge_itr->second == edge_coloring_itr->first.second) {
      //its neighbour
      if (edge_coloring_itr->second == current_color) {
       neighbour_has_it = true;
       break;
      }
     }
    }
    ++edge_coloring_itr;
   }
   if (!neighbour_has_it) {
    (*edge_coloring)[Edge_t{edge_itr->first,edge_itr->second}] = current_color;
    temp_colors.pop_back();
    if (edge_coloring->size() == graph->size()) {
     perm_found = true;
     break;
    }
   }
   ++edge_itr;
  }
 } while (next_permutation(all_colors.begin(), all_colors.end()) && !perm_found);
 if(perm_found){
  cout << "First sample found:\n";
  for(const auto &edge:*edge_coloring){
   cout << edge.first.first << "to " << edge.first.second << "is: " << edge.second <<"\n";
  }
 }*/

 //test
/* (*edge_coloring)[Edge_t{1,2}] = 3;
 if(edge_coloring->find(Edge_t{3,1}) != edge_coloring->end()){
  cout << "found!\n";
 }*/
 //test
 auto edge_itr = graph->begin();
 while(!perm_found){
/*  if(edge_coloring->find(*edge_itr) != edge_coloring->end()){
   continue;
  }*/
  //find this nodes neighbours
  vector<Graph_t ::iterator> neighbours;
  for(Graph_t::iterator inner_edge_itr = graph->begin();inner_edge_itr!=graph->end();++inner_edge_itr){
   if (!(inner_edge_itr->first == edge_itr->first && inner_edge_itr->second == edge_itr->second)){
    if (edge_itr->first == inner_edge_itr->first || edge_itr->first == inner_edge_itr->second
        || edge_itr->second == inner_edge_itr->first || edge_itr->second == inner_edge_itr->second){
     neighbours.push_back(inner_edge_itr);
    }
   }
  }
  vector<int> temp_colors(all_colors.size());
  copy(all_colors.begin(), all_colors.end(), temp_colors.begin());
  //this node and neighbours should be assigned unique values
  for(const auto & nbr_itr:neighbours){
   Edge_t temp_edge{nbr_itr->first,nbr_itr->second};
   if (edge_coloring->find(temp_edge) != edge_coloring->end()){
    temp_colors.erase(remove(temp_colors.begin(),temp_colors.end(),edge_coloring->at(temp_edge)),temp_colors.end());
   }
  }
  (*edge_coloring)[*edge_itr] = temp_colors.back();
  if(edge_coloring->size() == graph->size()){
   perm_found = true;
   break;
  }
  temp_colors.pop_back();
  ++edge_itr;
 }
 if(perm_found) {
  cout << "First sample found:\n";
  for (const auto &edge:*edge_coloring) {
   cout << edge.first.first << " to " << edge.first.second << " is: " << edge.second << "\n";
  }
 }
 return edge_coloring;
}

int main(int argc, char **argv) {

 read_input("../graph2.txt");
 find_suitable_sample();
 //assert(graph->size() == 3);
 return 0;
}

