//
// Project: Top-k Reliable Edge Colors in Uncertain Graphs
//
// Author:	Andreas Nufer
// E-Mail:	anufer@student.ethz.ch
// Date:    2015-07-01
// Summary: Contains the definitions of custom types, global constants and global variables
//

#include <vector>
#include <LEDA/core/dictionary.h>
#include <LEDA/core/h_array.h>
#include <LEDA/core/set.h>
#include <LEDA/core/random_source.h>
#include <LEDA/core/map2.h>
#include <LEDA/graph/graph.h>
#include <LEDA/graph/node_array.h>
#include <LEDA/core/p_queue.h>
#include <LEDA/core/slist.h>

using namespace std;
using namespace leda;

// Node, Color and Path types
typedef long Node;
typedef set<Node> NodeSet;
typedef h_array <leda_node, long> LedaNodeSet; // when using set<leda_node> GCC cannot compile

typedef long Color;
typedef set<Color> ColorSet;
typedef list<Color> ColorList;

typedef slist<leda_edge> Path;
typedef h_array <long, Path> PathSet;

// Helper types for top-k colors
typedef h_array <long, leda_node> AddressMap;

#ifndef COLREL_GRAPH_TYPES
#define COLREL_GRAPH_TYPES

// Graph: node and edge attribute definitons
typedef struct EdgeAttr_t { 
	double capacity;
	int color; 
	EdgeAttr_t() {}
	EdgeAttr_t(double capacity, int color) : capacity(capacity), color(color) {}
} EdgeAttr;

typedef struct NodeAttr_t { 
	Node id;
	int cost;
	int runId;
	leda_edge edgeToParent;
	NodeAttr_t() : id(0), cost(-1), runId(0), edgeToParent(NULL) {}
	NodeAttr_t(Node id) : id(id), cost(-1), runId(0), edgeToParent(NULL) {}
} NodeAttr;

// Graph definition
typedef GRAPH <NodeAttr, EdgeAttr> RelationGraph;

typedef struct ColoredGraph_t { 
	RelationGraph* Graph;
	map<long, leda_node> Nodes;
	set<long> Colors;
	bool directed;
	bool logCapacity; //edge capacity in this graph was transformed
	int runId;

	ColoredGraph_t() {}
	ColoredGraph_t(RelationGraph* Graph) : Graph(Graph) {}
} ColoredGraph;

//global variables
extern int BASELINE1_LOOP; //default: 1000
extern int BASELINE2_LOOP; //default: 100
extern bool USE_EPPSTEIN; //default: false
extern clock_t start;

#endif

//global contants
const Color NEUTRAL_COLOR = -1; //Neutral color gets special treatment during sampling
const int EVAL_LOOP = 50;
const int MONTECARLO_LOOP = 1000;
const int DIJKSTRA_LIMIT = 1000000;

//
// operators in main.cpp
//

// operators needed to use EdgeAttr and NodeAttr within GRAPH.
std::ostream &operator<<(std::ostream &os, EdgeAttr const &m);
std::istream &operator>>(std::istream &is, EdgeAttr const &m);

std::ostream &operator<<(std::ostream &os, NodeAttr const &m);
std::istream &operator>>(std::istream &is, NodeAttr const &m);