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

using namespace std;

static int total_colors = 3;
static int iterations = 100;
static vector<double> color_weights = {1, 2, 3};
const static double norm_factor = 0.01;

using ADJACENCY_T = vector<vector<int> >;
using VARS_T = int;
union table;
using TABLE_T = map<int, table>;
using MESSAGE_T = map<int, double>;
union table {
	shared_ptr<TABLE_T> table;
	double value;
};
struct Singleton;

struct Cluster {
	int index;
	map<int, shared_ptr<Singleton> > single_neighb;
	map<int, shared_ptr<MESSAGE_T> > incoming_messages;

	void initialize_messages() {
		//set outgoing messages to 1 for every outgoing message!
		auto itr = single_neighb.begin();
		while(itr != single_neighb.end()){
			incoming_messages[itr->first] = make_shared<MESSAGE_T>();
			auto &message = incoming_messages[itr->first];
			for(int i=1;i<= total_colors; ++i) {
				//setting every messages beliefs to 1 in the beginning
				message->at(i) = 1.0;
			}
			++itr;
		}
	}

	void send_messages(){
		auto itr = single_neighb.begin();
		while(itr != single_neighb.end()){
			auto &message = itr->second->incoming_messages.at(itr->first);
			for (int i = 1; i < total_colors; ++i) {
				//Now loop through incoming messages
				double temp = 1.0;
				auto incom_itr = incoming_messages.begin();
				while(incom_itr != incoming_messages.end()){
					if(incom_itr->first != itr->first){
						temp *= incom_itr->second->at(i);
					}
					++incom_itr;
				}
				message->at(i) = norm_factor * belief[i] * temp;
			}
			++itr;
		}
	}
};

struct Singleton {
	int index;
	//vector<shared_ptr<Cluster> > clust_neighb;
	map<int,shared_ptr<Cluster> > clust_neighb;
	map<int, double> belief;
	map<int, shared_ptr<MESSAGE_T> > incoming_messages;

	void initializeBeliefs() {
		for (int i = 1; i <= total_colors; ++i) {
			belief[i] = exp(color_weights[i]);
		}
	}

	void initialize_messages() {
		//set outgoing messages to 1 for every outgoing message!
		auto itr = clust_neighb.begin();
		while(itr != clust_neighb.end()){
			incoming_messages[itr->first] = make_shared<MESSAGE_T>();
			auto &message = incoming_messages[itr->first];
			for(int i=1;i<= total_colors; ++i) {
				//setting every messages beliefs to 1 in the beginning
				message->at(i) = 1.0;
			}
			++itr;
		}
	}

	void send_messages(){
		auto itr = clust_neighb.begin();
		while(itr != clust_neighb.end()){
			auto &message = itr->second->incoming_messages.at(itr->first);
			for (int i = 1; i < total_colors; ++i) {
				//Now loop through incoming messages
				double temp = 1.0;
				auto incom_itr = incoming_messages.begin();
				while(incom_itr != incoming_messages.end()){
					if(incom_itr->first != itr->first){
						temp *= incom_itr->second->at(i);
					}
					++incom_itr;
				}
				message->at(i) = norm_factor * belief[i] * temp;
			}
			++itr;
		}
	}
};

using ALL_BELIEFS = struct AllBeliefs {
	vector<shared_ptr<Cluster> > clusters;
	vector<shared_ptr<Singleton> > singletons;
	int max_degree;

	AllBeliefs(const vector<shared_ptr<Cluster> > &clusters, const vector<shared_ptr<Singleton> > &singletons,
						 int max_degree) : clusters(clusters), singletons(singletons), max_degree(max_degree) {}
};

ADJACENCY_T read_input(string path) {
	ifstream ifile(path);
	vector<vector<int> > adjacency;
	if (!ifile.is_open()) {
		cout << "file opening problem! " << path << "\n";
		exit(1);
	}
	string line;
	while (getline(ifile, line)) {
		vector<int> vertex;
		stringstream strm(line);
		string val;
		while (getline(strm, val, ',')) {
			vertex.push_back(stoi(val));
		}
		adjacency.push_back(vertex);
	}
	ifile.close();
	return adjacency;
}

ALL_BELIEFS prepare_beliefs(ADJACENCY_T &adjacency) {
	int total_edges = 0;
	int max_degree = 0;
	vector<shared_ptr<Cluster> > clusters;
	vector<shared_ptr<Singleton> > singletons;
	vector<vector<shared_ptr<Singleton> > > edge_map(adjacency.size());
	for (int i = 0; i < adjacency.size(); i++) {
		vector<int> &vertex = adjacency[i];
		auto vert_clust = make_shared<Cluster>();
		clusters.push_back(vert_clust);
		vert_clust->index = ++i;
		edge_map[i] = vector<shared_ptr<Singleton> >(vertex.size());

		int current_max_degree = 0;
		for (int j = 0; j < vertex.size(); ++j) {
			if (vertex[j] == 1) {
				current_max_degree++;
				if (i < j) {
					auto vert_single = make_shared<Singleton>();
					singletons.push_back(vert_single);
					vert_single->index = ++total_edges;
					vert_single->clust_neighb[vert_clust->index] =  vert_clust;
					edge_map[i][j] = vert_single;
					vert_clust->single_neighb[vert_single->index] = vert_single;

				} else {
					edge_map[i][j] = edge_map[j][i];
					edge_map[i][j]->clust_neighb[vert_clust->index] = vert_clust;
					vert_clust->single_neighb[edge_map[i][j]->index] = edge_map[i][j];
				}
			}
		}
		if (current_max_degree > max_degree) {
			max_degree = current_max_degree;
		}
	}
	cout << "Total Number of Edges are: " << total_edges << "\n";
	return AllBeliefs(clusters, singletons, max_degree);
}

double compute_sum_product(ALL_BELIEFS &all_beliefs, vector<double> &color_weights,
													 int total_colors, int interations) {


	return 0.0;
}

int main(int argc, char **argv) {

	ADJACENCY_T adjacency = read_input("../graph.txt");
	ALL_BELIEFS beliefs = prepare_beliefs(adjacency);

	return 0;
}
