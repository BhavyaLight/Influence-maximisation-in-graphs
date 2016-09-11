//
// Project: Top-k Reliable Edge Colors in Uncertain Graphs
//
// Author:	Andreas Nufer
// E-Mail:	anufer@student.ethz.ch
// Date:    2015-07-01
// Summary: Contains the heuristic and baseline implementations for the MAX-sum optimization. 
//

#include "functions.h"

// private constants
const long NEW_S = -1;
const long NEW_T = -2;

//
// Gets the edge between u and v from graph G.
// Accessibility: private
//
static inline leda_edge getEdge(Node u, Node v, const ColoredGraph& G)
{
	leda_node u1 = G.Nodes[u];
	leda_node v1 = G.Nodes[v];
	
	leda_edge e;
	forall_adj_edges(e, u1)
	{
		leda_node t = G.Graph->target(e);

		if(t==v1)
		{
			return e;
		}
	}

	cout << "Error: getEdge: Cannot find edge between u and v! (u,v)=(" << u << "," << v << ")" << endl;
	exit(-1);
	//throw std::exception("Error: getEdge: Cannot find edge between u and v!");
}

//
// Same functionality as topRShortestPaths in single.cpp but with a set of destination nodes. This
// will reduce the search depth of the shortest paths by 1, wich may lead to better results under
// certain circumstances. The paths will be added to the parameter topr_paths.
// Accessibility: private
//
static void topRShortestPathsTSet(const ColoredGraph& G, Node s, NodeSet T, int r, PathSet& topr_paths)
{
	struct EdgeWithIndex {
		leda_edge edge;
		long prevIndex;
		EdgeWithIndex() {}
		EdgeWithIndex(leda_edge edge, long prevIndex) : edge(edge), prevIndex(prevIndex) {}
	};
	h_array <long, EdgeWithIndex> path_index; 

	p_queue <double, long> all_paths;
	h_array <long, long> path_end, path_end_prev;

	path_index.clear();
	all_paths.clear();
	path_end.clear();
	path_end_prev.clear();
	topr_paths.clear();

	int i = 0;
	path_index[i] = EdgeWithIndex(NULL, -1);
	all_paths.insert(0, i);
	path_end[i] = s;
	path_end_prev[i] = -1;

	int pathsFound = 0;
	if(T.member(s)) pathsFound++;

	const NodeSet end_nodes = T;
	
	while(all_paths.size() != 0 && pathsFound < r && i < DIJKSTRA_LIMIT) 
	{
		pq_item min_item = all_paths.find_min();
        long j = all_paths.inf(min_item);
        double cost = all_paths.prio(min_item);
        all_paths.del_min();

        Node s1 = path_end[j];
		leda_node u = G.Nodes[s1];
        if(!end_nodes.member(s1)) 
		{
			leda_edge e;
			forall_adj_edges(e, u) 
			{
				leda_node v = G.Graph->target(e);
				Node t1 = G.Graph->inf(v).id;
				if(t1 != path_end_prev[j]) 
				{
					i++;
					path_index[i] = EdgeWithIndex(e, j);

					double capacity = -log(G.Graph->inf(e).capacity);
					all_paths.insert(cost + capacity, i);
					path_end[i] = t1;
					path_end_prev[i] = s1;

					if(end_nodes.member(t1)) 
					{
						pathsFound++;

						Path& pathList = topr_paths[i];
						pathList.clear();
						pathList.push(getEdge(t1, NEW_T, G));
						
						int z = i;
						while (z > 0)
						{
							EdgeWithIndex& elem = path_index[z];
							z = elem.prevIndex;
							pathList.push(elem.edge);
						}
					}
				}
			}
        }  
    }

	path_index.clear();
	all_paths.clear();
	path_end.clear();
	path_end_prev.clear();
}

//
// Adds a new source (destination) point that are connected to all nodes in S (T) with color NEUTRAL_COLOR and 
// capacity of newCapacity.
// Accessibility: private
//
std::vector<leda_edge> extendGraph(ColoredGraph* G1,  NodeSet& S, NodeSet& T, double newCapacity)
{
	//add new s1 and t1 to graph
	RelationGraph& Gr = *(G1->Graph);
	// Add new source node with id -1
	leda_node s1 = Gr.new_node(NodeAttr(NEW_S));
	// Add node-1 = s1
	G1->Nodes[NEW_S] = s1;
	//  Add new target node
	leda_node t1 = Gr.new_node(NodeAttr(NEW_T));
	// Add node in node list with id -2
	G1->Nodes[NEW_T] = t1;

	std::vector<leda_edge> newEdges;

	Node s;
	forall(s, S) 
	{
		EdgeAttr edgeAttr = EdgeAttr(newCapacity, NEUTRAL_COLOR);
		newEdges.push_back(Gr.new_edge(s1, G1->Nodes[s], edgeAttr));
	}

	Node t;
	forall(t, T) 
	{
		EdgeAttr edgeAttr = EdgeAttr(newCapacity, NEUTRAL_COLOR);
		newEdges.push_back(Gr.new_edge(G1->Nodes[t], t1, edgeAttr));
	}

	return newEdges;
}

//
// Removes the edges and nodes that were added in extendGraph
// Accessibility: private
//
static void reduceGraph(ColoredGraph& G1, std::vector<leda_edge> newEdges)
{
	//revert changes in graph
	for(leda_edge e : newEdges)
	{
		G1.Graph->del_edge(e);
	}
	leda_node s1 = G1.Nodes[NEW_S];
	leda_node t1 = G1.Nodes[NEW_T];

	G1.Graph->del_node(s1);
	G1.Graph->del_node(t1);
	G1.Nodes[NEW_S] = NULL;
	G1.Nodes[NEW_T] = NULL;
}

//
// Implementation of heuristic for MAX sum.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    k: the number of colors to return
//    r: the number of top-reliable paths to use.
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
//
ColorList* maxSum(ColoredGraph& G, NodeSet& S, NodeSet& T, int k, int r) 
{
	//add new s1 and t1 to graph
	double newCapacity = G.logCapacity ? 0 : 1;
	std::vector<leda_edge> newEdges = extendGraph(&G, S, T, newCapacity);

	//get top k-colors
	double p = 0;
	PathSet P;
	if(USE_EPPSTEIN)
		topRShortestPaths_eppstein(G, NEW_S, NEW_T, r, P);
	else
		topRShortestPathsTSet(G, NEW_S, T, r, P);

	ColorList* L_top = topKColors(P, G, NEW_S, NEW_T, k, &p);

	// clean up path set
	long i;
	forall_defined(i, P) 
		P[i].clear();

	P.clear();
	
	//revert changes in graph
	reduceGraph(G, newEdges);
	return L_top;
}

//
// Implementation of baseline1 for MAX sum.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    k: the number of colors to return
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
//
ColorList* maxSum_baseline1(ColoredGraph& G, NodeSet& S, NodeSet& T, int k) 
{
	std::vector<leda_edge> newEdges = extendGraph(&G, S, T, 1);

	//get top k-colors
	double p = 0;
	ColorList* L_top = baseline1(G, NEW_S, NEW_T, k, &p);

	reduceGraph(G, newEdges);
	return L_top;
}

//
// Implementation of baseline2 for MAX sum.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    k: the number of colors to return
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected.
//
ColorList* maxSum_baseline2(ColoredGraph& G, NodeSet& S, NodeSet& T, int k) 
{
	std::vector<leda_edge> newEdges = extendGraph(&G, S, T, 1);

	//get top k-colors
	double p = 0;
	ColorList* L_top = baseline2(G, NEW_S, NEW_T, k, &p);

	reduceGraph(G, newEdges);
	return L_top;
}

//
// Implementation of the result evaluation for MAX sum. Uses sampling to estimate the result, 
// hence the result value is not deterministic.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    selectedColors: the set of selected colors to evaluate
// Returns: The evaluation metric for the selected colors. Value is between 0 and EVAL_LOOP (see types.h).
//
double maxsum_evaluateMetric(ColoredGraph& G, NodeSet S, NodeSet T, const h_array <long, int>& selectedColors) 
{
	std::vector<leda_edge> newEdges = extendGraph(&G, S, T, 1);

	h_array <long, int> selectedColorsPlusNeutral(selectedColors);
	selectedColorsPlusNeutral[NEUTRAL_COLOR] = 1;

	double metric = single_evaluateMetric(G, NEW_S, NEW_T, selectedColorsPlusNeutral);

	reduceGraph(G, newEdges);
	return metric;
}