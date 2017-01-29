//
// Project: Top-k Reliable Edge Colors in Uncertain Graphs
//
// Author:	Andreas Nufer
// E-Mail:	anufer@student.ethz.ch
// Date:    2015-07-01
// Summary: Contains the heuristic implementation for the single source/destination case. 
//

#include "functions.h"

//
// Returns the colors of path P wich is a subgraph of G.
// Accessibility: private
//
static inline ColorSet getColors(const Path& P, ColoredGraph& G)
{
	ColorSet L_p;
	L_p.clear();

	leda_edge e;
	forall(e, P)
	{
		if(G.Graph->inf(e).color != NEUTRAL_COLOR) //NEUTRAL_COLOR is used to get max sum
			L_p.insert(G.Graph->inf(e).color);
	}
	return L_p;
}

//
// Implementation of the Dijkstra algorithm to find the top-r shortest paths between s and t.
// Accessibility: public
// Parameters:
//	  G: the graph.
//	  s: the id of the source node. Must be element of G.
//    t: the id of the destination node. Must be element of G.
//    r: the number of paths to return.
//    topr_paths: Return parameter. A set to store the paths in.
//
void topRShortestPaths(ColoredGraph& G, Node s, Node t, int r, PathSet& topr_paths)
{
	if(USE_EPPSTEIN)
	{
		topRShortestPaths_eppstein(G, s, t, r, topr_paths);
		return;
	}

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
	if(s == t) pathsFound++;

	const Node end_node = t;
	
	while(all_paths.size() != 0 && pathsFound < r && i < DIJKSTRA_LIMIT) 
	{
		pq_item min_item = all_paths.find_min();
        long j = all_paths.inf(min_item);
        double cost = all_paths.prio(min_item);
        all_paths.del_min();

        Node s1 = path_end[j];
		leda_node u = G.Nodes[s1];
        if(s1 != end_node) 
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

					if(t1 == end_node)
					{
						pathsFound++;

						Path& pathList = topr_paths[i];
						pathList.clear();

						int z = i;
						while (z > 0)
						{
							EdgeWithIndex& elem = path_index[z];
							z = elem.prevIndex;
							pathList.push(elem.edge);
						}

						cout << i <<"> " << printPath(pathList, G) << endl;
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
// Monte Carlo sampling
// Accessibility: public
// Parameters:
//	  s1:   the source node. Must be element of G.
//    t1:   the destination node. Must be element of G.
//    rand: a random number source
//    Gc1:  the graph.
// Returns: A reachability estimation between s1 and t1.
//
inline double monteCarlo(leda_node s1, leda_node t1, random_source& rand, ColoredGraph& Gc1) 
{
	RelationGraph& G1 = *Gc1.Graph;

	slist<leda_node> nodes1, nodes2;
	slist<leda_node> &old_n = nodes1, &new_n = nodes2;

	// actual MC Simulation
	double count = 0;
	for(int k = 1; k < MONTECARLO_LOOP; k++) 
	{  
		old_n.clear();
		decrementRunId(Gc1);

		old_n.append(s1);
		G1.assign(s1, NodeAttr(Gc1.runId));

		bool terminate = false;
		while(terminate != true && old_n.size() != 0) 
		{
			new_n.clear();

			leda_node u;
			forall(u, old_n) 
			{
				leda_edge e;
				forall_adj_edges(e, u) 
				{
					leda_node v = G1.target(e);
					int lastVisited = G1.inf(v).id; //absuse node id as last visited flag. node id is not needed for baseline calculation
					if(lastVisited != Gc1.runId)
					{
						double num;
						rand >> num;
						if(num <= G1.inf(e).capacity)
						{                      
							new_n.append(v);
							G1.assign(v, NodeAttr(Gc1.runId));
							if(v == t1) 
							{
								terminate = true;
								count++;
								goto nextwhile;
							}
						}
					}
				}
			}

			//swap
			slist<leda_node> &tmp = old_n;
			old_n = new_n;
			new_n = tmp;
		}
		nextwhile:;
	}
	return count;
}

//
// Adds node u to G1 and registers reverse lookup in M1
// Accessibility: public
//
// inline leda_node addNode(RelationGraph& G1, Node u, AddressMap& M1)
// {
// 	leda_node u1;
// 	if(!M1.defined(u)) {
// 		u1 = G1.new_node(NodeAttr(u));
// 		M1[u]= u1;
// 	}
// 	else
// 		u1 = M1[u];

// 	return u1;
// }

//
// Implementation of top-K colors algorithm.
// Accessibility: public
// Parameters:
//    P: a set of pre selected paths
//	  G: the graph.
//	  s: the id of the source node. Must be element of G.
//    t: the id of the destination node. Must be element of G.
//    k: the number of colors to return.
//    p: Return type. The sum of the sampling counts of the selected colors.
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
//
ColorList* const topKColors(PathSet& P, ColoredGraph& G, Node s, Node t, int k, double* p)
{
	*p = 0;
	bool pNeedsRefresh = false;
	random_source S(1,10);

	//init top-r paths from P
	list<long> topr_paths; //= P;
	topr_paths.clear();
	int elem;
	// Copy all paths into another list
	forall_defined(elem, P)
	{
		topr_paths.push(elem);
	}

	RelationGraph G1;
	ColoredGraph Gc1;
	Gc1.Graph = &G1;

	AddressMap M1;
	// G1 graph containes the present selected paths
	leda_node s_G1 = addNode(G1, s, M1);  //Create a new node S in G, M[s]=new node created in G

	int monteCarloCount = 0;

	ColorList* L1 = new ColorList(); //selected colors to return
	ColorSet selected_colors; //for fast lookup
	
	bool terminate = false;

    while(topr_paths.size() != 0 && terminate != true) 
	{
		list<long>::item it;
		//For each path 'it' in the list of topr_paths
		forall_items(it, topr_paths) //remove all path with the same colors as already selected
		{
			//Get the path
			int i = topr_paths.contents(it);
			Path& path = P[i];
			//Add all colors
			ColorSet current_colors = getColors(path, G).join(selected_colors);
			if(current_colors.size() - selected_colors.size() == 0) // current path doesn't add any new colors -> add to selected paths immediately
			{
				topr_paths.del_item(it);
                
				// Construct the path from source node s_G1 through all edges in the path till the target edge
				leda_node u_G1 = s_G1;
				leda_edge e;
				forall(e, path)
				{
					EdgeAttr e_G = G.Graph->inf(e);
					if(G.logCapacity) e_G.capacity = exp(-e_G.capacity / 1000);
					Node v_G = G.Graph->inf(G.Graph->target(e)).id;	//Return id of the target node of edge e
					leda_node v_G1 = addNode(G1, v_G, M1); //add node v_G to G, M1[v_G]= the newly created node in G

					G1.new_edge(u_G1, v_G1, e_G);     
					if(!G.directed) G1.new_edge(v_G1, u_G1, e_G);
					u_G1 = v_G1;
				}
				pNeedsRefresh = true;
			}
			// if it is >k, simply delete the path 
			else if(current_colors.size() > k)
			{
				topr_paths.del_item(it);
			}
		}

		list<long>::item best_path_item = NULL;
		double best_path_count = -1;//any negative value will work

		// At the present stage, finds the best path that will be added to get max reliability
		//colors of any path in topr_paths combined with selected colors is guaranteed to be lower than k.
        forall_items(it, topr_paths) 
		{

			int i = topr_paths.contents(it);
			Path& path = P[i];
			
			leda_node u_G1 = s_G1;

			leda::list<leda_edge> edgeList;
			leda_edge e;
			// Add each path in graph G1
			forall(e, P[i])
			{
				EdgeAttr e_G = G.Graph->inf(e);
				if(G.logCapacity) e_G.capacity = exp(-e_G.capacity / 1000);
				Node v_G = G.Graph->inf(G.Graph->target(e)).id;	
				leda_node v_G1 = addNode(G1, v_G, M1); //add node to G1

				edgeList.push(G1.new_edge(u_G1, v_G1, e_G));
				if(!G.directed) edgeList.push(G1.new_edge(v_G1, u_G1, e_G));

				u_G1 = v_G1;
			}

			monteCarloCount++;
			double count = monteCarlo(M1[s], M1[t], S, Gc1); //monte carlo sampling cost estimation
			if(count > best_path_count)
			{
				best_path_count = count;
				best_path_item = it;
			}

			*p = std::max(*p, count);

			// remove edges from current path from G1
			G1.del_edges(edgeList);
			edgeList.clear();
		}

		if(best_path_item != NULL) 
		{
			int i = topr_paths.contents(best_path_item);
			topr_paths.del_item(best_path_item);
			
			Path& P_best = P[i];
			
			leda_node u_G1 = s_G1;
			leda_edge e;
			forall(e, P_best)
			{
				EdgeAttr e_G = G.Graph->inf(e);
				Node v_G = G.Graph->inf(G.Graph->target(e)).id;	
				leda_node v_G1 = M1[v_G]; //node is guaranteed to be in M1

				G1.new_edge(u_G1, v_G1, e_G);
				if(!G.directed) G1.new_edge(v_G1, u_G1, e_G);
				u_G1 = v_G1;

				if(e_G.color != NEUTRAL_COLOR && !selected_colors.member(e_G.color)) //NEUTRAL_COLOR is used to get max sum
				{
					selected_colors.insert(e_G.color);
					L1->append(e_G.color);
				}
			}
			pNeedsRefresh = false;
		}
		else
			terminate = true;

    } //end while
	
	if(pNeedsRefresh) //refresh p when there were path selected without doing a monteCarlo estimation.
	{
		monteCarloCount++;
		double cost = monteCarlo(M1[s], M1[t], S, Gc1); //monte carlo sampling cost estimation
		std::max(*p, cost);
	}
	return L1;
}


//
// Implementation of heuristic for single source/destination.
// Accessibility: public
// Parameters:
//	  G1: the graph data.
//	  s: the id of the source node. Must be element of G.
//    t: the id of the destination node. Must be element of G.
//    k: the number of colors to return.
//    r: the number of top-reliable paths to use.
//    p: Return type. The sum of the sampling counts of the selected colors.
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
//
ColorList* single(ColoredGraph& G1, Node s, Node t, int k, int r, double *p) 
{
	PathSet P;
	topRShortestPaths(G1, s, t, r, P);
	ColorList* L = topKColors(P, G1, s, t, k, p);

	// clean up path set
	long i;
	forall_defined(i, P) P[i].clear();
	P.clear();

	return L;
}

//
// Implementation of the result evaluation for single source/destination. Uses sampling to estimate the result, 
// hence the result value is not deterministic.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  s: the id of the source node. Must be element of G
//    t: the id of the destination node. Must be element of G
//    selectedColors: the set of selected colors to evaluate
// Returns: The evaluation metric for the selected colors. Value is between 0 and EVAL_LOOP (see types.h).
//
double single_evaluateMetric(ColoredGraph& G, Node s, Node t, const h_array <long, int>& selectedColors) 
{
	random_source S(1,10);

	slist<leda_node> nodes1, nodes2;
	slist<leda_node> &old_n = nodes1, &new_n = nodes2;

	const leda_node s1 = G.Nodes[s];
	const leda_node t1 = G.Nodes[t];

	int count = 0;
	for(int i = 1; i <= EVAL_LOOP; i++) 
	{
		decrementRunId(G);
		old_n.clear();

		old_n.append(s1);
		G.Graph->assign(s1, NodeAttr(G.runId));

		bool terminate = false;
		while(old_n.size() > 0 && !terminate) 
		{
			new_n.clear();

			leda_node u;
			forall(u, old_n) 
			{
				leda_edge e;
				forall_adj_edges(e, u) 
				{
					leda_node v = G.Graph->target(e);
					const EdgeAttr& ea = G.Graph->inf(e);
					int lastVisited = G.Graph->inf(v).id; //absuse node id as last visited flag. node id is not needed for baseline calculation
					if(lastVisited != G.runId && selectedColors.defined(ea.color))
					{
						long num;
						S >> num;
						if(num <= ea.capacity * 10)
						{
							new_n.append(v);
							G.Graph->assign(v, NodeAttr(G.runId));

							if(v == t1) 
							{
								terminate = true;
								count++;
								goto exitwhile;
							}
						}
					}
				}
			}

			slist<leda_node> &tmp = old_n;
			old_n = new_n;
			new_n = tmp;
		}
exitwhile:;
	}

	return count;
}