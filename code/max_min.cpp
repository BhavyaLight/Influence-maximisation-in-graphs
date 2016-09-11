//
// Project: Top-k Reliable Edge Colors in Uncertain Graphs
//
// Author:	Andreas Nufer
// E-Mail:	anufer@student.ethz.ch
// Date:    2015-07-01
// Summary: Contains the heuristic and baseline implementations for the MAX-minimum optimization. 
//

#include "functions.h"
#include <LEDA/core/sortseq.h> 

typedef long PathId;

//
// Returns the colors of path P wich is a subgraph of G.
// Accessibility: private
//
static inline ColorSet getColors(Path& P, ColoredGraph& G)
{
	ColorSet L_p;
	L_p.clear();

	leda_edge e;
	forall(e, P)
	{
		L_p.insert(G.Graph->inf(e).color);
	}
	return L_p;
}

//
// Returns the source node of path P wich is a subgraph of G.
// Accessibility: private
//
static inline Node getStart(Path& p, ColoredGraph& G)
{
	leda_edge e = p.inf(p.first_item());
	return G.Graph->inf(G.Graph->source(e)).id;
}

//
// Returns the destination node of path P wich is a subgraph of G.
// Accessibility: private
//
static inline Node getEnd(Path& p, ColoredGraph& G)
{
	leda_edge e = p.inf(p.last_item());
	return G.Graph->inf(G.Graph->source(e)).id;
}

//
// Returns the colors set of minimal size such that each PathSet in G_pairs has at least
// one path that is fully covered by the selected colors. k denotes the maximum number of colors selected.
// The colors are added to the set passed in L_min.
// Accessibility: private
//
static void findMinimumColorSet(h_array<int, PathSet> G_pairs, ColoredGraph& G, int k, ColorSet& L_min)
{
	h_array<PathId, Path> all_paths;
	h_array<PathId, int> P_graph_map;
	sortseq<int, PathId> P_remain_sort;

	set<int> G_remain;
	set<int> G_done;

	// init G_remain and P_remain
	int j = 1;
	int i;
	forall_defined(i, G_pairs)
	{
		G_remain.insert(i);
		
		int k = 0;
		Path path;
		forall(path, G_pairs[i])
		{
			all_paths[j] = path;
			P_graph_map[j] = i;
			P_remain_sort.insert(k * 1000 + i, j); //make sure that the top paths at the beginning of a pathset are also at the beginning of the sequence
			j++;
			k++;
		}
	}

	while (!G_remain.empty() && !P_remain_sort.empty())
	{
		// find min color
		int colDiff = INT_MAX;

		PathId P = 0; // arg min...;

		seq_item it;
		forall_items(it, P_remain_sort) 
		{
			PathId pathId = P_remain_sort.inf(it);
			Path& path = all_paths[pathId];
			ColorSet col_path = getColors(path, G);
			int diffSize = L_min.diff(col_path).size();
			if(diffSize < colDiff)
			{
				colDiff = diffSize;
				P = pathId;
			}

			if(diffSize == 0)
				break;
		}
		
		L_min += getColors(all_paths[P], G);
		while(L_min.size() > k)
		{
			long c_remove = 0;
			long c;
			forall(c, L_min)
			{
				ColorSet L_c = L_min;
				L_c.del(c);
				bool connected = true;
				int g;
				forall(g, G_done)
				{
					bool innerCon = false;
					Path p;
					forall(p, G_pairs[g])
					{
						if(getColors(p, G).diff(L_c).empty())
						{
							innerCon = true;
							break;
						}
					}

					if(!innerCon)
					{
						connected = false;
						break;
					}
				}

				if(connected)
				{
					c_remove = c;
					break;
				}
			}

			if(c_remove == 0)
			{
				//abort, return any colors
				L_min.clear();//clear already selected colors
				return;
			}

			L_min.del(c_remove);
		}

		int g = P_graph_map[P]; // get graph for P and remove it from G_remain
		G_remain.del(g); 
		G_done.insert(g);

		//remove all paths from P_remain_sort that have the same path set as P
		seq_item it2;
		forall_items(it2, P_remain_sort) 
		{
			PathId path = P_remain_sort.inf(it2);
			if(P_graph_map[path] == g)
			{
				P_remain_sort.del_item(it2);
			}
		}
	}
}

//
// Estimates the probability to reach t from s with in the graph containing all edges in P with color in L.
// Accessibility: private
//
static double RCol_G(PathSet& P, Node s, Node t, ColorSet L, ColoredGraph& G)
{
	random_source rand(1,10);
	RelationGraph G1;
	ColoredGraph Gc1;
	Gc1.Graph = &G1;
	AddressMap M1;

	long i;
	forall_defined(i, P) 
	{
		Path& path = P[i];
		
		// add start node of path to G
		leda_edge e_first = path.inf(path.first_item());
		Node v_G = G.Graph->inf(G.Graph->source(e_first)).id;	

		leda_node u_G1 = addNode(G1, s, M1);
		leda_edge e;
		forall(e, path)
		{
			EdgeAttr e_G = G.Graph->inf(e);
			if(G.logCapacity) e_G.capacity = exp(-e_G.capacity / 1000);
			Node v_G = G.Graph->inf(G.Graph->target(e)).id;	
			leda_node v_G1 = addNode(G1, v_G, M1); //add node to G1

			if(L.member(e_G.color))
			{
				G1.new_edge(u_G1, v_G1, e_G);
				if(!G.directed) G1.new_edge(v_G1, u_G1, e_G);
			}
			u_G1 = v_G1;
		}
	}
	//s1: start node, t1: end node
	double cost = monteCarlo(M1[s], M1[t], rand, Gc1);
	G1.clear();
	M1.clear();

	return cost;
}

//
// Find additional colors to L1 to maximize the minimal probaility to reach any pair.
// Accessibility: private
//
static inline ColorSet& topKColorsOfMaxMin(h_array<int, PathSet> G_pairs, ColoredGraph& G, int k, ColorSet& L1)
{
	while(L1.size() < k)
	{
		set<long> selected_paths;
		selected_paths.clear();

		// select graph G_min with the lowest connectivity for L1
		double prob_min = MAXDOUBLE;
		int g_min = 0;
		int g;
		forall_defined(g, G_pairs)
		{
			PathSet& G_p = G_pairs[g];
			
			Path first_path = G_p.inf(G_p.first_item());
			double prob = RCol_G(G_p, getStart(first_path, G), getEnd(first_path, G), L1, G);

			if(prob < prob_min)
			{
				prob_min = prob;
				g_min = g;
			}
		}
		PathSet& G_min = G_pairs[g_min];

		double prob_opt = 0;
		PathSet::item it_opt = NULL;
		PathSet::item it;
		forall_items(it, G_min)
		{
			if(selected_paths.member(G_min.key(it) + g_min * 100000)) //key + g_min * 100000 = unique id of path in select_paths
				continue;

			Path p = G_min.inf(it);
			ColorSet L_p = getColors(p, G).join(L1);
			if(L_p.size() > k || L_p.size() - L1.size() == 0) // check for max color <= k constraint.
				continue;

			double prob = RCol_G(G_min, getStart(p, G), getEnd(p, G), L_p, G);

			if(prob > prob_opt)
			{
				prob_opt = prob;
				it_opt = it;
			}
		}

		if(it_opt != NULL)
		{
			long key = G_min.key(it_opt);
			Path P_opt= G_min.inf(it_opt);
			selected_paths.insert(key + g_min * 100000); //will work for r < 100000
			L1 += getColors(P_opt, G);
		}
		else
		{
			return L1;
		}
	}

	return L1;
}

//
// Implementation of heuristic for MAX minimum.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    k: the number of colors to return
//    r: the number of top-reliable paths to use.
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
//
ColorList* maxMin(ColoredGraph& G, NodeSet& S, NodeSet& T, int k, int r) 
{
	h_array<int, PathSet> topPathSets;

	int i = 1;
	Node t;
	forall(t, T) 
	{
		Node s;
		forall(s, S)
		{
			if(s == t) continue;
			PathSet& P = topPathSets[i];
			topRShortestPaths(G, s, t, r, P);
			i++;

			if(P.empty()) //if one s-t pair has no top-r paths, the max-min has only randomly chosen colors
			{
				return new ColorList();
			}
		}
	}

	ColorSet L_topk;
	findMinimumColorSet(topPathSets, G, k, L_topk);
	if(!L_topk.empty())
		topKColorsOfMaxMin(topPathSets, G, k, L_topk); 
	
	// clean up topPathSets
	PathSet P;
	forall(P, topPathSets)
		P.clear();
	topPathSets.clear();

	ColorList* CL1 = new ColorList();
	Color c1;
	forall(c1, L_topk) CL1->append(c1);
	return CL1;
}

//
// Sampling to estimate the probaility to reach t1 from s1 using the colors "cl" \union "base_colors"
// Accessibility: private
//
static long baseline_relP(ColoredGraph& G,
				   leda_node s1,
				   leda_node t1,
				   Color cl,
				   const h_array <long, int>* base_colors,
				   random_source& rand,
				   bool anyColor,
				   int loopCount) 
{
	clock_t timeLimit = start + CLOCKS_PER_SEC * 1000000; //set time limit in sec
	bool abort = false;

	slist<leda_node> nodes1, nodes2;
	slist<leda_node> &old_n = nodes1, &new_n = nodes2;

	int count = 0;
	for(int i = 1; i <= loopCount && !abort; i++) 
	{
		decrementRunId(G);
		
		old_n.clear();
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

					bool inBaseColors = anyColor || (base_colors != NULL && base_colors->defined(ea.color));

					int lastVisited = G.Graph->inf(v).id; //absuse node id as last visited flag. node id is not needed for baseline calculation
					if(lastVisited != G.runId && (ea.color == cl || inBaseColors || ea.color == NEUTRAL_COLOR))
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
								goto exitwhile;
							}
						}
					}
				}
			}

			slist<leda_node> &tmp = old_n;
			old_n = new_n;
			new_n = tmp;
			abort = clock() > timeLimit;
		} // end while
exitwhile:;
	}
	if(abort) count = -1;
	return count;
}

//
// Implementation of baseline1 for MAX minimum.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    k: the number of colors to return
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
//
ColorList* maxMin_baseline1(ColoredGraph& G, NodeSet& S, NodeSet& T, int k) 
{
	random_source rand(1,10);
	
	p_queue <double, long> sort_colors;

	sort_colors.clear();
	bool abort = false;
	long cl;
	forall(cl, G.Colors) 
	{
		long min_count = LONG_MAX;
		
		long t;
		forall(t, T) 
		{
			long s;
			forall(s, S) 	
			{
				if(s == t || min_count == 0 || abort) continue;

				leda_node s1 = G.Nodes[s];
				leda_node t1 = G.Nodes[t];
				
				long count = baseline_relP(G, s1,t1, cl, NULL, rand, false, BASELINE1_LOOP);
				min_count = std::min(min_count, count);
				abort = count == -1;
			}
		}
		sort_colors.insert(-min_count, cl);
	}

	ColorList* L = new ColorList();
	for(int i = 1; i <= k; i++) 
	{
		pq_item min_color = sort_colors.find_min();
		double prio = sort_colors.prio(min_color);
		if(-prio > 0)
		{
			L->append(sort_colors.inf(min_color));
		}
		sort_colors.del_min();
	}
	if(abort) L->append(-1);
	return L;
}

//
// Implementation of baseline2 for MAX minimium.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    topk: the number of colors to return
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected.
//
ColorList* maxMin_baseline2(ColoredGraph& G, NodeSet& S, NodeSet& T, int topk) 
{
	bool abort = false;
	ColorList* L = new ColorList();
	random_source rand(1,10);

	p_queue <double, long> sort_colors;
	sort_colors.clear();

	h_array <long, int> colors_del, colors_add;
	colors_del.clear();
	colors_add.clear();   

	h_array <leda_node, long> present, old_n, new_n; 

	{ // limit scope of variable cl
		Color cl;
		forall(cl, G.Colors)
		{
			colors_del[cl] = 1;
		}
	}
	
	for(int k = 1; k <= topk && !abort; k++) 
	{
		sort_colors.clear();
		Color cl;
		forall_defined(cl, colors_del) 
		{
			long min_count = LONG_MAX;
			long s;
			forall(s, S) 
			{
				long t;
				forall(t, T) 
				{
					if(s == t || min_count == 0 || abort) continue;

					leda_node s1 = G.Nodes[s];
					leda_node t1 = G.Nodes[t];
				
					long count = baseline_relP(G, s1, t1, cl, &colors_add, rand, false, BASELINE2_LOOP);
					min_count = std::min(min_count, count);
					abort = count == -1;
				}
			}
			sort_colors.insert(-min_count, cl);
		}

		pq_item min_col = sort_colors.find_min();
		Color c = sort_colors.inf(min_col);
		colors_del.undefine(c);
		colors_add[c]=1;
		
		double prio = sort_colors.prio(min_col);
		L->append(c);
		cout << c << " p:" << -prio << endl;
	}
	if(abort) L->append(-1);
	return L;
}

//
// Implementation of the result evaluation for MAX minimum. Uses sampling to estimate the result, 
// hence the result value is not deterministic.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    selectedColors: the set of selected colors to evaluate
// Returns: The evaluation metric for the selected colors. Value is between 0 and EVAL_LOOP (see types.h).
//
double maxMin_evaluateMetric(ColoredGraph& G, NodeSet S, NodeSet T, const h_array <long, int>& selectedColors) 
{
	double metric_min = DBL_MAX;
	long s;
	forall(s, S) 
	{
		long t;
		forall(t, T) 
		{
			if(s == t || metric_min == 0) continue; //don't need to test the other combinations if metric becomes 0
			double metric = single_evaluateMetric(G, s, t, selectedColors);
			metric_min = std::min(metric, metric_min);
		}
	}

	return metric_min;
}
