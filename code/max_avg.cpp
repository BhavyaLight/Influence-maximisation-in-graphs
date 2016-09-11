//
// Project: Top-k Reliable Edge Colors in Uncertain Graphs
//
// Author:	Andreas Nufer
// E-Mail:	anufer@student.ethz.ch
// Date:    2015-07-01
// Summary: Contains the specific implementations for the MAX-average optimization. Heuristic and baselines are the same as for MAX-sum. 
//          See max_sum.cpp for the remaining implementations.
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
// Adds a new source (destination) point that are connected to all nodes in S (T) with color NEUTRAL_COLOR and 
// capacity of newCapacity.
// Accessibility: private
//
static std::vector<leda_edge> extendGraph(ColoredGraph* G1,  NodeSet& S, NodeSet& T, double newCapacity)
{
	//add new s1 and t1 to graph
	RelationGraph& Gr = *(G1->Graph);
	leda_node s1 = Gr.new_node(NodeAttr(NEW_S));
	G1->Nodes[NEW_S] = s1;
	leda_node t1 = Gr.new_node(NodeAttr(NEW_T));
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
// Implementation of heuristic for MAX maximum.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    k: the number of colors to return
//    r: the number of top-reliable paths to use.
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
//
ColorList* maxAvg(ColoredGraph& G, NodeSet& S, NodeSet& T, int k, int r) 
{
	//add new s1 and t1 to graph
	double newCapacity = G.logCapacity ? 0 : 1;
	std::vector<leda_edge> newEdges = extendGraph(&G, S, T, newCapacity);

	PathSet combinedPaths;

	int i = 1;
	Node t;
	forall(t, T) 
	{
		Node s;
		forall(s, S)
		{
			if(s == t) continue;
			PathSet P;
			topRShortestPaths(G, s, t, r, P);			

			Path path;
			forall(path, P) 
			{
				Path newPath;
				newPath.append(getEdge(NEW_S, s, G));
				leda_edge e; forall(e, path) newPath.append(e);
				newPath.append(getEdge(t, NEW_T, G));

				i++;
				combinedPaths[i] = newPath;
			}
			P.clear();
		}
	}

	//get top k-colors
	double p = 0;
	ColorList* L_top = topKColors(combinedPaths, G, NEW_S, NEW_T, k, &p);

	// clean up path set
	long j;
	forall_defined(j, combinedPaths) 
		combinedPaths[j].clear();

	combinedPaths.clear();
	
	//revert changes in graph
	reduceGraph(G, newEdges);
	return L_top;
}

//
// Implementation of baseline1 for MAX maximum.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    k: the number of colors to return
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contain less than k colors.
//
ColorList* maxAvg_baseline1(ColoredGraph& G, NodeSet& S, NodeSet& T, int k) 
{
	clock_t timeLimit = start + CLOCKS_PER_SEC * 1000000; //set time limit in sec
	bool abort = false;
	random_source rand(1,10);
	p_queue <double, long> sort_colors;
	sort_colors.clear();

	slist<leda_node> nodes1, nodes2;
	slist<leda_node> &old_n = nodes1, &new_n = nodes2;

	long cl;
	forall(cl, G.Colors) 
	{
		int count = 0;
		for(long i = 1; i <= BASELINE1_LOOP && !abort; i++) 
		{
			long s;
			forall(s, S) 
			{
				long t;
				forall(t, T) 
				{
					if(s == t) continue;

					leda_node s1 = G.Nodes[s];
					leda_node t1 = G.Nodes[t];
					old_n.clear();
					new_n.clear();
					decrementRunId(G);

					old_n.append(s1);
					G.Graph->assign(s1, NodeAttr(G.runId));
					abort = clock() > timeLimit;
					bool terminate = false;
					while(old_n.size() > 0 && !terminate && !abort) 
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
								if(lastVisited != G.runId && (ea.color == cl || ea.color == NEUTRAL_COLOR)) 
								{
									double num;
									rand >> num;
									if(num <= ea.capacity) 
									{
										new_n.append(v);
										G.Graph->assign(v, NodeAttr(G.runId));

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
						abort = clock() > timeLimit;
					}
				nextwhile:;
				}
			}
		}
		sort_colors.insert(-count, cl);
	}

	ColorList* L = new ColorList();
	for(int i = 1; i <= k; i++) 
	{
		pq_item min_color = sort_colors.find_min();
		double prio = sort_colors.prio(min_color);
		if(-prio > 0)
		{
			L->append(sort_colors.inf(min_color));
			cout << sort_colors.inf(min_color) << " p:" << -prio << endl;
		}
		sort_colors.del_min();
	}
	if(abort) L->append(-1);
	return L;
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
ColorList* maxAvg_baseline2(ColoredGraph& G, NodeSet& S, NodeSet& T, int topk) 
{
	clock_t timeLimit = start + CLOCKS_PER_SEC * 1000000; //set time limit in sec
	bool abort = false;
	ColorList* L = new ColorList();
	random_source rand(1,10);

	p_queue <double, long> sort_colors;
	sort_colors.clear();

	h_array <long, int> colors_del;
	colors_del.clear();
	set<Color> colors_add;

	slist<leda_node> nodes1, nodes2;
	slist<leda_node> &old_n = nodes1, &new_n = nodes2;

	{ // limit scope of variable cl
		Color cl;
		forall(cl, G.Colors)
		{
			colors_del[cl] = 1;
		}
	}
	

	for(int k = 1; k <= topk && !abort; k++) 
	{
		int maxn = 0;
		long appendCount = 0;
		sort_colors.clear();
		Color cl;
		forall_defined(cl, colors_del) 
		{
			int count = 0;
			for(int i = 1; i <= BASELINE2_LOOP && !abort; i++) 
			{
				long s;
				forall(s, S) 
				{
					long t;
					forall(t, T) 
					{
						if(s == t) continue;
						const leda_node s1 = G.Nodes[s];
						const leda_node t1 = G.Nodes[t];

						old_n.clear();
						new_n.clear();
						decrementRunId(G);

						old_n.append(s1);
						G.Graph->assign(s1, NodeAttr(G.runId));

						abort = clock() > timeLimit;
						bool terminate = false;
						while(old_n.size() > 0 && !terminate && !abort) 
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
									if(lastVisited != G.runId && (ea.color == cl || colors_add.member(ea.color) || ea.color == NEUTRAL_COLOR)) 
									{
										double num;
										rand >> num;
										if(num <= ea.capacity) 
										{
											new_n.append(v);
											G.Graph->assign(v, NodeAttr(G.runId));

											appendCount++;
											if(v == t1) {
												terminate = true;
												count++;
												goto nextwhile;
											}
										}
									}
								}
							}
							maxn = std::max(maxn, new_n.size());

							//swap
							slist<leda_node> &tmp = old_n;
							old_n = new_n;
							new_n = tmp;
							abort = clock() > timeLimit;
						}
						nextwhile:;
					}
				}
			}
			sort_colors.insert(-count, cl);
		}

		pq_item min_col = sort_colors.find_min();
		Color c = sort_colors.inf(min_col);
		colors_del.undefine(c);
		colors_add.insert(c);

		double prio = sort_colors.prio(min_col);

		L->append(c);
		cout << c << " p:" << -prio << " #" << appendCount << " mx:" << maxn << endl;
		//if(-prio == BASELINE2_LOOP) break;
	}
	if(abort) L->append(-1);
	return L;
}

//
// Implementation of the result evaluation for MAX average. Uses sampling to estimate the result, 
// hence the result value is not deterministic.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    selectedColors: the set of selected colors to evaluate
// Returns: The evaluation metric for the selected colors. Value is between 0 and EVAL_LOOP (see types.h).
//
double maxavg_evaluateMetric(ColoredGraph& G, NodeSet S, NodeSet T, const h_array <long, int>& selectedColors) 
{
	double metric_sum = 0;
	long s;
	forall(s, S) 
	{
		long t;
		forall(t, T) 
		{
			if(s == t) continue;
			double metric = single_evaluateMetric(G, s, t, selectedColors);
			metric_sum = metric_sum + metric;
		}
	}

	return metric_sum / (S.size() * T.size());
}
