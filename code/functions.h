//
// Project: Top-k Reliable Edge Colors in Uncertain Graphs
//
// Author:	Andreas Nufer
// E-Mail:	anufer@student.ethz.ch
// Date:    2015-07-01
// Summary: Contains the definitions of the projects global accessible functions
//

#include <LEDA/core/dictionary.h>
#include <LEDA/core/h_array.h>
#include <LEDA/core/set.h>
#include <LEDA/core/random_source.h>

#include "types.h"

//
// functions in main.cpp
//
void decrementRunId(ColoredGraph& G);

//
// functions in single.cpp
//
ColorList* single(ColoredGraph& G, Node s, Node t, int k, int r, double* p);
double single_evaluateMetric(ColoredGraph& G, Node s, Node t, const h_array <long, int>& selectedColors);

//common used functions/algorithms in single.cpp
double monteCarlo(leda_node s1, leda_node t1, random_source& rand, ColoredGraph& Gc1);
leda_node addNode(RelationGraph& G1, Node u, AddressMap& M1);
void topRShortestPaths(ColoredGraph& G, Node s, Node t, int r, PathSet& topr_paths);
ColorList* const topKColors(PathSet& P, ColoredGraph& G, Node s, Node t, int k, double* p) ;

//
// functions in eppstein.cpp
//
void topRShortestPaths_eppstein(ColoredGraph& G, Node s, Node t, int r, PathSet& topr_paths);

//
// functions in baseline.cpp
//
ColorList* baseline1(ColoredGraph& G, Node s, Node t, int k, double *p);
ColorList* baseline2(ColoredGraph& G, Node s, Node t, int k, double *p);

//
// functions in max_sum.cpp
//
ColorList* maxSum(ColoredGraph& G, NodeSet& S, NodeSet& T, int k, int r);
ColorList* maxSum_baseline1(ColoredGraph& G, NodeSet& S, NodeSet& T, int k);
ColorList* maxSum_baseline2(ColoredGraph& G, NodeSet& S, NodeSet& T, int k);
double maxsum_evaluateMetric(ColoredGraph& G, NodeSet S, NodeSet T, const h_array <long, int>& selectedColors);

//
// functions in max_avg.cpp
//
double maxavg_evaluateMetric(ColoredGraph& G, NodeSet S, NodeSet T, const h_array <long, int>& selectedColors);
ColorList* maxAvg(ColoredGraph& G, NodeSet& S, NodeSet& T, int k, int r);
ColorList* maxAvg_baseline1(ColoredGraph& G, NodeSet& S, NodeSet& T, int k);
ColorList* maxAvg_baseline2(ColoredGraph& G, NodeSet& S, NodeSet& T, int k);

//
// functions in max_max.cpp
//
ColorList* maxMax(ColoredGraph& G, NodeSet& S, NodeSet& T, int k, int r);
ColorList* maxMax_baseline1(ColoredGraph& G, NodeSet& S, NodeSet& T, int k);
ColorList* maxMax_baseline2(ColoredGraph& G, NodeSet& S, NodeSet& T, int k);
double maxmax_evaluateMetric(ColoredGraph& G, NodeSet S, NodeSet T, const h_array <long, int>& selectedColors);

//
// functions in max_min.cpp
//
ColorList* maxMin(ColoredGraph& G, NodeSet& S, NodeSet& T, int k, int r);
ColorList* maxMin_baseline1(ColoredGraph& G, NodeSet& S, NodeSet& T, int k);
ColorList* maxMin_baseline2(ColoredGraph& G, NodeSet& S, NodeSet& T, int k);
double maxMin_evaluateMetric(ColoredGraph& G, NodeSet S, NodeSet T, const h_array <long, int>& selectedColors);

//
// functions in max_conn.cpp
//
ColorList* maxConn(ColoredGraph& G, LedaNodeSet& S, int k, int r);
ColorList* maxConn_baseline1(ColoredGraph& G, LedaNodeSet& S, int k);
ColorList* maxConn_baseline2(ColoredGraph& G, LedaNodeSet& S, int topk);
double maxConn_evaluateMetric(ColoredGraph& G, LedaNodeSet& S, const h_array <long, int>& selectedColors);

//
// functions in max_avg.cpp
//
double maxRel_evaluateMetric(ColoredGraph& G, NodeSet S, NodeSet T, const h_array <long, int>& selectedColors);
ColorList* maxRel(ColoredGraph& G, NodeSet& S, NodeSet& T, int k, int r);

//
// functions in debug.cpp
//
std::string printColorSet(ColorSet& cs);
std::string printPath(Path& p, ColoredGraph& G1);
std::string printPathSet(PathSet& ps);
std::string printGraph(RelationGraph& G);
std::string printHashtable(h_array<leda_node, leda_edge>& hmap);
std::string printEdge(const ColoredGraph& G, leda_edge edge);
std::string printEdges(slist<leda_edge> edges, const ColoredGraph& G1);
