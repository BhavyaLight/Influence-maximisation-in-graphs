//
// Project: Top-k Reliable Edge Colors in Uncertain Graphs
//
// Author:	Andreas Nufer
// E-Mail:	anufer@student.ethz.ch
// Date:    2015-07-01
// Summary: Contains the eppstein implementation for the single source/destination case. 
//

#include "functions.h"
#include <LEDA/core/sortseq.h> 

// type definitions
typedef sortseq <double, leda_edge> hout_heap;
typedef sortseq <double, leda_node> hg_heap;
typedef struct 
{ 
	bool isDirected; 
	bool isInit; 
	leda_node t; //root node of tree
	p_queue <int, leda_node> front_nodes; 
	int nodes_done; 
} ShortestPathTreeState;

// static member
static ShortestPathTreeState* last_state;

//
// debug functions
//
std::string printHout(hout_heap hmap, const ColoredGraph& G)
{
	std::stringstream ss;

	hout_heap::item it;
	forall_items(it, hmap)
	{
		leda_edge e = hmap.inf(it);
		double key = hmap.key(it);
		auto ea = G.Graph->inf(e);
		ss  << key << ": (" << ea.capacity <<"," <<ea.color << " -> (" << G.Graph->inf(G.Graph->source(e)).id <<", " << G.Graph->inf(G.Graph->target(e)).id << ") " << endl;
	}
	return ss.str();
}

std::string printPQueue(hg_heap hmap,const ColoredGraph& G)
{
	std::stringstream ss;

	hg_heap::item it;
	forall_items(it, hmap)
	{
		leda_node e = hmap.inf(it);
		auto na = G.Graph->inf(e);
		ss << hmap.key(it) << ": (" << na.id << "," << na.cost << ") " << endl;
	}
	return ss.str();
}
std::string printHashtable(h_array<leda_edge, int>& hmap,const ColoredGraph& G)
{
	std::stringstream ss;

	h_array<leda_edge, int>::item it;
	forall_items(it, hmap)
	{
		leda_edge e = hmap.key(it);
		auto ea = G.Graph->inf(e);
		ss << "(" << ea.capacity <<"," <<ea.color << " -> (" << hmap.inf(it) << ") " << endl;
	}
	return ss.str();
}
std::string printNode(const ColoredGraph& G, leda_node node)
{
	std::stringstream ss;

	NodeAttr na = G.Graph->inf(node);
	ss << "(" << na.id <<"," << na.cost << ") " << endl;
	return ss.str();
}

leda_node inf(hg_heap::item it)
{
	hg_heap h;
	return h.inf(it);
}
leda_edge inf2(hout_heap::item it)
{
	hout_heap h;
	return h.inf(it);
}

//
// helper functions
//

//
// Calcuates the shortest path tree of every node in the graph G to the target node (state.t). If G is undirected the tree building will be interrupted as
// soon as the shortest path from stop_node to state.t has been found. The current state of the calculation is stored in state and can be resumed any time later.
// The resulting tree is stored on the nodes and edges of G.
// Accessibility: private
// Parameters:
//	  G: the graph.
//	  state: stores the datastructures needed to calcuate the shortest path tree and allows to resume the algorithm after is has been interrupted.
//    stop_node: Interrupts the tree building as soon a shortest path to stop_node has been found.
//
static void buildShortestPathTree(ColoredGraph& G, ShortestPathTreeState* state, leda_node stop_node)
{
	// nodes can have 2 states: the states are flagged using NodeAttr.runId
	// node.runId == Graph.runId: node cost and edgetoparent are initialized (=has a Path To T). They may be changed if there is a better path
	// node.runId == Graph.runId + 1: node has been selected as min item in the priority queue. the shortest path to this node is now complete and wont be changed.

	//T = shortest path tree 
	//E = G - T
	int& nodesDone = state->nodes_done;
	leda_node t1 = state->t;
	p_queue <int, leda_node>& front_nodes = state->front_nodes;
	
	if(!state->isInit)
	{
		state->isInit = true;
		decrementRunId(G);
		decrementRunId(G);
		
		front_nodes.insert(0, t1);
		NodeAttr na_t1 = G.Graph->inf(t1);
		na_t1.cost = 0;
		na_t1.edgeToParent = NULL;
		G.Graph->assign(t1, na_t1); //update node cost
	}

	int runId = G.runId;
	while(!front_nodes.empty()) 
	{
		pq_item min_item = front_nodes.find_min();
        leda_node u = front_nodes.inf(min_item);
        double cost = front_nodes.prio(min_item);
        front_nodes.del_min();

		NodeAttr na_u = G.Graph->inf(u);
		if(na_u.runId == G.runId + 1) continue; //node was already done
		na_u.runId = G.runId + 1;
		G.Graph->assign(u, na_u); //update node cost
		++nodesDone;

		leda_edge e;
		forall_in_edges(e, u) 
		{
			leda_node v = G.Graph->source(e);
			int v_runId = G.Graph->inf(v).runId;
			bool visited = v_runId == runId + 1;
			
			if(!visited) //not yet visited
			{
				bool hasPathToT = v_runId == runId;
				const EdgeAttr& ea = G.Graph->inf(e);
				const NodeAttr& na = G.Graph->inf(v);
				int totalCost = cost + (int)ea.capacity;
				if(!hasPathToT || totalCost < na.cost)
				{
					NodeAttr na2;
					na2.cost = totalCost;
					na2.edgeToParent = e;
					na2.runId = G.runId;
					na2.id = na.id;
					G.Graph->assign(v, na2); //update node cost
					front_nodes.insert(na2.cost, v);
				}
			}
		}

		//cannot abort early in a directed graph. there would be no way to know if a node is unreachable or hasn't just found the shortest path yet.
		if(u == stop_node && !G.directed) return; 
    }
}

//
// Calcuates the penalty of the given sidetrack compared to the optimal shortest path from source(e). It uses state to resume the shortest path tree building,
// if a shortest path for target(e) hasn't been found yet.
// Accessibility: private
// Parameters:
//	  G: the graph.
//    e: edge of sidetrack to get the penalty for.
//	  state: current state of the shortest path tree building.
// Returns: the penalty of using e. -1 if e cannot be used in valid a shortest path.
//
static inline int getPenalty(ColoredGraph& G, leda_edge e, ShortestPathTreeState* state)
{
	if(e == NULL) return -1;
	leda_node u = G.Graph->source(e);
	leda_node v = G.Graph->target(e);
	if(G.Graph->inf(u).edgeToParent == e) return -1;

	int v_runId = G.Graph->inf(v).runId;
	bool vHasShortestPath = v_runId == G.runId + 1; //v is connected to t
	if(!vHasShortestPath) {
		if(state != NULL && !state->isDirected){
			buildShortestPathTree(G, state, v);
			return getPenalty(G, e, NULL); //state=NULL to avoid infinte loop
		}
		return -1;
	}

	int edgeCost = (int)G.Graph->inf(e).capacity;
	int vCost = G.Graph->inf(v).cost;
	int uCost = G.Graph->inf(u).cost;
	int penalty = (vCost + edgeCost) - uCost;
	return penalty;
}

//
// Builds a heap of the out edges of node v with their corresponding penalty as key. Will return cached version if heap has been bulit before.
// Accessibility: private
// Parameters:
//	  G: the graph.
//    v: node
//    all_out: cache of all built heaps.
//	  state: current state of the shortest path tree building.
//
static hout_heap& buildHout(ColoredGraph& G, leda_node v, h_array<leda_node, hout_heap>& all_hout, ShortestPathTreeState* state)
{
	if(all_hout.defined(v)) return all_hout[v];

	hout_heap& hout = all_hout[v];
	leda_edge e;
	forall_out_edges(e, v)
	{
		int p = getPenalty(G, e, state);
		if(p >= 0) //p=-1 means e is in T (shortest path tree) or cannot be used for other reasons.
		{
			hout.insert(p, e);
		}
	}
	return hout;
}

//
// Builds a heap of all nodes on the shortest path from v to t their corresponding min penalty (see buildHout) as key. Will return cached version if heap has been bulit before.
// Accessibility: private
// Parameters:
//	  G: the graph.
//    v: the node to build the heap for
//    t: the target node of the shortest path
//    all_hg: cache of all built hg-heaps.
//    all_hout: cache of all built hout-heaps.
//	  state: current state of the shortest path tree building.
// Returns: h(v)
//
static hg_heap& buildH_G(ColoredGraph& G, leda_node v, leda_node t, h_array<leda_node, hg_heap>& all_hg, h_array<leda_node, hout_heap>& all_hout, 
				  ShortestPathTreeState* state)
{
	//test if hg_v has already been built
	if(all_hg.defined(v))
	{
		return all_hg[v];
	}

	hg_heap& hgThis = all_hg[v];

	if(v == t)
	{
		//if(penalty >= 0) //we are not interested in paths that loop in t as it will give us no benefit to find better colors.
		//	hgThis.insert(penalty, v); //-1=no sidepaths for this node
		return hgThis;
	}
	else
	{
		//build h_out
		hout_heap& h_out_v = buildHout(G, v, all_hout, state);
		int penalty = h_out_v.empty() ? -1 : h_out_v.key(h_out_v.find_min());

		leda_node nextT = G.Graph->opposite( G.Graph->inf(v).edgeToParent, v);
		hg_heap& hgNextT = buildH_G(G, nextT, t, all_hg, all_hout, state);

		hg_heap::item it;
		forall_items(it, hgNextT)
		{
			hgThis.insert(hgNextT.key(it), hgNextT.inf(it));
		}
		if(penalty >= 0) 
			hgThis.insert(penalty, v); //-1=no sidepaths for this node	
	}

	return hgThis;
}

//
// Implementation of the Eppstein algorithm to find the top-r shortest paths between s and t.
// Accessibility: public
// Parameters:
//	  G: the graph.
//	  s: the id of the source node. Must be element of G.
//    t: the id of the destination node. Must be element of G.
//    r: the number of paths to return.
//    topr_paths: Return parameter. A set to store the paths in.
//
void topRShortestPaths_eppstein(ColoredGraph& G, Node s, Node t, int r, PathSet& topr_paths)
{
	struct QueueItem {
		hout_heap::item sidetrack;
		hout_heap* hout_heap_p;
		Path parentPath;
		hg_heap::item hg_heap_item;
		hg_heap* hg_heap_p;
		QueueItem() {}
		QueueItem(hout_heap* hout_heap_p, hout_heap::item sidetrack, Path parentPath, hg_heap* hg_heap_p, hg_heap::item hg_heap_item) 
			: hout_heap_p(hout_heap_p), sidetrack(sidetrack), parentPath(parentPath), hg_heap_item(hg_heap_item), hg_heap_p(hg_heap_p) 
		{}
	};

	if(G.logCapacity == false) {
		cout << "edges must be log converted" << endl;
		exit(-1);
	}
	leda_node s1 = G.Nodes[s];
	leda_node t1 = G.Nodes[t];

	if(last_state == NULL)
	{
		last_state = new ShortestPathTreeState();
		last_state->t = NULL;
	}

	ShortestPathTreeState& state = *last_state;
	int sss = state.front_nodes.size();
	if(state.t != t1)
	{
		state.front_nodes.clear();
		state.nodes_done = 0;
		state.t = t1;
		state.isDirected = G.directed;
		state.isInit = false;
	}
	else {
		state.nodes_done = 0;
	}

	buildShortestPathTree(G, &state, s1);

	bool s1_visited = G.Graph->inf(s1).runId == G.runId + 1;
	if(!s1_visited) return; //s and t are not connected
	
	// add shortest path to topr_paths
	{
		Path& ps = topr_paths[0];
		leda_node n = s1;
		while(n != t1)
		{
			leda_edge e = G.Graph->inf(n).edgeToParent;
			ps.append(e);
			n = G.Graph->target(e);
		}
		cout << printPath(ps, G) << endl;
	}

	h_array<leda_node, hout_heap> all_hout;
	h_array<leda_node, hg_heap> all_hg;	
	hg_heap& hg_s = buildH_G(G, s1, t1, all_hg, all_hout, &state); 
	
	p_queue <double, QueueItem> Q;
	hg_heap::item hgMinS1 = hg_s.find_min();
	hout_heap emptyHeap;
	hout_heap& hout_s1 = hgMinS1 == NULL ? emptyHeap : all_hout[hg_s.inf(hgMinS1)];
	hout_heap::item fs1 = hout_s1.empty() ? NULL : hout_s1.find_min();
	Q.insert(0, QueueItem(&hout_s1, fs1, Path(), &hg_s, hgMinS1));

	int pathsFound = 1;
	while(!Q.empty() && pathsFound < r)
	{
		pq_item min_item = Q.find_min();
		QueueItem item = Q.inf(min_item);
		Path S_k = item.parentPath;
        double cost_currentSidePath = Q.prio(min_item);//cost_currentSidePath = cost of S_k + item.sidetrack
		Q.del_min();

		if(item.sidetrack == NULL) continue;
		leda_edge lastSidetrack = item.hout_heap_p->inf(item.sidetrack);
		//add to topr_paths
		{
			Path S_k1(S_k);
			S_k1.append(lastSidetrack);

			Path& ps = topr_paths[pathsFound++]; //create new path in topr_paths
			leda_node n = s1;
			Path::item it = S_k1.first_item();
			while(n!=t1)
			{
				if(it != NULL && G.Graph->source(S_k1.inf(it)) == n)
				{
					ps.append(S_k1.inf(it));
					n = G.Graph->target(S_k1.inf(it));
					it = S_k1.next_item(it);
				}
				else
				{
					leda_edge e = G.Graph->inf(n).edgeToParent;
					ps.append(e);
					n = G.Graph->target(e);
				}
			}
			cout << cost_currentSidePath <<"> " << printPath(ps, G) << endl;
		}
		
		leda_node v = G.Graph->target(lastSidetrack);
		
		//get next best along new path reachable throgh latest sidetrack
		hg_heap& hg_head_e = buildH_G(G, v,  t1, all_hg, all_hout, &state);

		hg_heap::item hgMin = hg_head_e.find_min();
		if(hgMin != NULL)
		{
			hout_heap::item f = all_hout[hg_head_e.inf(hgMin)].find_min();
			if(f != NULL)
			{
				// insert S_k \dot f
				Path S_k1(S_k);
				S_k1.append(lastSidetrack);

				double score = cost_currentSidePath + getPenalty(G, all_hout[hg_head_e.inf(hgMin)].inf(f), NULL);
				Q.insert(score, QueueItem(&all_hout[hg_head_e.inf(hgMin)], f, S_k1, &hg_head_e, hgMin));
			}
		}

		//add paths from D(G)
		
		//get next best in the same node as the current sidetrack
		//get egde from queue from this
		hout_heap::item nextEdge = item.hout_heap_p->next_item(item.sidetrack);
		if(nextEdge != NULL)
		{
			double score = cost_currentSidePath - getPenalty(G, lastSidetrack, NULL) + getPenalty(G, item.hout_heap_p->inf(nextEdge), NULL);
			Q.insert(score, QueueItem(item.hout_heap_p, nextEdge, S_k, NULL, NULL)); // no dot further continue within hg_heap
		}
		
		//get next best along the old path as the current sidetrack
		hg_heap::item nextHgItem = item.hg_heap_p != NULL ? item.hg_heap_p->next_item(item.hg_heap_item) : NULL;
		if(nextHgItem != NULL)
		{
			leda_node nextNode = item.hg_heap_p->inf(nextHgItem);
			hout_heap& nextHout = all_hout[nextNode];
			hout_heap::item nextEdge2 = nextHout.first_item();
			if(nextEdge2 != NULL)
			{
				double score = cost_currentSidePath - getPenalty(G, lastSidetrack, NULL) + getPenalty(G, nextHout.inf(nextEdge2), NULL);
				Q.insert(score, QueueItem(&nextHout, nextEdge2, S_k, item.hg_heap_p, nextHgItem));
			}
		}
	}
	cout << "Nodes collected: " << state.nodes_done << endl;
	Q.clear();
	all_hg.clear();
}