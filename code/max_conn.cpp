//
// Project: Top-k Reliable Edge Colors in Uncertain Graphs
//
// Author:	Andreas Nufer
// E-Mail:	anufer@student.ethz.ch
// Date:    2015-07-01
// Summary: Contains the heuristic and baseline implementations for the MAX-connectivity optimization. 
//

#include "functions.h"

using namespace std;

//
// datatype definitions
//
typedef struct Tree_t { 
	double cost; 
	long p_hash;
	leda_node v; 
	int looseEdges;
	pq_item pqueue_item; //only for use inside updateQueue()
	Tree_t* prev; 
	Tree_t* prev2;
	leda_edge edge;
	int nextAltRoute;
	bool initComplete;
} Tree; //T(v, p) = cost

typedef struct { int id; ColorSet colors; slist<leda_edge> edges; } SteinerTree;
typedef h_array<int, SteinerTree> SteinerTreeSet;
typedef h_array<long, Tree> TreeMap;

//
// functions
//

//
// Returns true if a steinter tree in sts has the same set of edges as in edges.
// Accessibility: private
//
static bool containsTree(SteinerTreeSet& sts, slist<leda_edge>& edges, RelationGraph &G)
{
	int i;
	forall_defined(i, sts)
	{
		slist<leda_edge>& edges2 = sts[i].edges;
		if(edges.size() == edges2.size()) 
		{
			bool contains = true;
			leda_edge e1;
			forall(e1, edges)
			{
				bool hasEdge = false;
				leda_edge e2;
				forall(e2, edges2)
				{
					if(e1 == e2)
						hasEdge = true;
				}
				if(!hasEdge) contains = false;
			}
			if(contains) return true;
		}
	}
	return false;
}

//
// Collects the edges of the steiner tree represented by the given Tree t and the set of alternative routes.
// Accessibility: private
//
static void collectEdges(Tree* t, set<int>* altRoutes, slist<leda_edge>& edges)
{
	if(t == NULL) return;
	
	//procceed with default route
	if(t->edge != NULL)
	{
		edges.append(t->edge);
	}

	collectEdges(t->prev, altRoutes, edges);
	collectEdges(t->prev2, altRoutes, edges);
}

//
// Update Tree t in queue with new priority of t.
// Accessibility: private
//
static inline void updateQueue(p_queue<double, Tree*>& queue, Tree& t)
{
	double prio = t.cost;

	if(t.pqueue_item == NULL)
	{
		t.pqueue_item = queue.insert(prio, &t);
		return;
	}

	double oldPrio = queue.prio(t.pqueue_item);
	if(prio < oldPrio) 
	{
		queue.decrease_p(t.pqueue_item, prio);
	}
	else if(prio > oldPrio) 
	{
		queue.del_item(t.pqueue_item);
		t.pqueue_item = queue.insert(prio, &t);
	}
}

//
// Lookup Tree in all_T. If tree doesn't exist it will be initizalied with v and p_hash and cost = MAXDLB
// Accessibility: private
// 
static inline Tree& getTree(h_array<leda_node, TreeMap>& all_T, leda_node v, long p_hash)
{
	TreeMap& v_trees = all_T[v];
	Tree& Tup = v_trees[p_hash];
	
	if(!Tup.initComplete)
	{
		Tup.p_hash = p_hash;
		Tup.v = v;
		Tup.cost = DBL_MAX;

		Tup.pqueue_item = NULL;
		Tup.prev = NULL;
		Tup.prev2 = NULL;
		Tup.nextAltRoute = 0;
		Tup.looseEdges = 0;
		Tup.edge = NULL;
	}

	return Tup;
}

//
// Implementation of the DPBF-r (get top-r SteinerTree) algorithm from Bolin Ding et al.
// Accessibility: private
// 
static inline void getTopRSteinerTrees(ColoredGraph& G, LedaNodeSet& S, int r, SteinerTreeSet& topR) //DPBF-r
{
	if(S.size() > 32) //not more than 32 elements in S supported
	{
		cout << "Error: not more than 32 elements in S supported!" << endl;
		exit(-1);
	}

	int altRouteIndex = 1;
	p_queue<double, Tree*> Q_T;
	h_array<leda_node, TreeMap> all_T; //using map<leda_node, TreeMap> will cause memory leaks

	long S_hash = (1 << S.size()) - 1;// the hash of the complete set has the last |S| bits set to 1

	// init trees for all nodes in S
	int nextHash = 1;
	leda_node v_init;
	forall_defined(v_init, S)
	{
		Tree& tvp = getTree(all_T, v_init, nextHash);
		tvp.cost = 0;
		tvp.initComplete = true;

		updateQueue(Q_T, tvp);
		nextHash = nextHash << 1; //update hash for next item
	}

	int found = 0;
	int i = 0;
	while (!Q_T.empty() && found < r && i < 1500000)//TODO determine i max
	{
		double minQ_TPrio = !Q_T.empty() ? Q_T.prio(Q_T.find_min()) : DBL_MAX;
		i++;

		if(Q_T.empty()) break; //leave while loop
		Tree& Tvp = *Q_T[Q_T.find_min()];
		Q_T.del_min();

		leda_node v = Tvp.v;
		TreeMap& v_trees = all_T[v];

		int p_hash = Tvp.p_hash;
		if(p_hash == S_hash /*&& Tvp.looseEdges==0*/)
		{
			slist<leda_edge> edges;
			collectEdges(&Tvp, NULL, edges);

			if(!containsTree(topR, edges, *G.Graph))  //tree was already added, expand tree will not result in new tree as all S are already captured
			{
				SteinerTree& st = topR[found + 1];
				st.edges = edges;
				st.id = found + 1;

				// get distinct colors of tree
				leda_edge e;
				forall(e, st.edges)
				{
					st.colors.insert(G.Graph->inf(e).color);
				}

				cout << "found " << found << ": " << printEdges(st.edges, G) << endl;
				found++;
				edges.clear();
			}
			continue; // no further expansion on this tree
		}

		//
		// expand on edges
		//
		leda_edge e;
		if(G.directed)
		{
			forall_adj_edges(e, v)
			{
				leda_node u = G.Graph->target(e);
				Tree& Tup = getTree(all_T, u, p_hash);

				double edgeCost = -log(G.Graph->inf(e).capacity);
				if(Tvp.cost + edgeCost < Tup.cost) //T(v,p) + (v,u) < T(u,p)
				{
					//update Tup
					Tup.cost = Tvp.cost + edgeCost;
					Tup.prev = &Tvp;
					Tup.prev2 = NULL;
					Tup.edge = e;
					Tup.initComplete = true;
					updateQueue(Q_T, Tup);
				}
			}
		}
		else
		{
			forall_inout_edges(e, v)
			{
				leda_node u = G.Graph->opposite(e, v);
				Tree& Tup = getTree(all_T, u, p_hash);

				double edgeCost = -log(G.Graph->inf(e).capacity);
				if(Tvp.cost + edgeCost < Tup.cost) //T(v,p) + (v,u) < T(u,p)
				{
					//update Tup
					Tup.cost = Tvp.cost + edgeCost;
					Tup.prev = &Tvp;
					Tup.prev2 = NULL;
					Tup.edge = e;
					Tup.initComplete = true;
					updateQueue(Q_T, Tup);
				}
			}
		}
		
		//
		// merge trees
		//
		int p1_hash = p_hash;

		slist<int> disjunctHashes;
		int p2_hash;
		forall_defined(p2_hash, v_trees) //need to copy hashes into a tree so we can modify the v_tree in the next loop
		{
			if((p1_hash & p2_hash) == 0) disjunctHashes.append(p2_hash);
		}

		forall(p2_hash, disjunctHashes)
		{
			Tree& Tvp2 = v_trees[p2_hash];

			int p12_hash = p1_hash + p2_hash;
			Tree& Tvp12 = getTree(all_T, v, p12_hash); //modifies v_trees
			if(Tvp.cost + Tvp2.cost < Tvp12.cost)
			{
				Tvp12.cost = Tvp.cost + Tvp2.cost;
				Tvp12.prev = &Tvp;
				Tvp12.prev2 = &Tvp2;
				Tvp12.initComplete = true;

				//Tvp12.looseEdges = Tvp.looseEdges + Tvp2.looseEdges;
				//if(!S.defined(v)) Tvp12.looseEdges -= 2;

				updateQueue(Q_T, Tvp12);
			}
		}
		disjunctHashes.clear();
	}
	
	//cleanup and release memory
	leda_node j; forall_defined(j, all_T) all_T[j].clear();
	all_T.clear();
	Q_T.clear();
}

//
// Copies node from sourceGraph to destinationGraph.
// Accessibility: private
//
static inline leda_node copyNode(leda_node sourceNode, ColoredGraph& sourceGraph, ColoredGraph& targetGraph)
{
	NodeAttr na = sourceGraph.Graph->inf(sourceNode);
	leda_node u1;
	if(!targetGraph.Nodes.defined(na.id)) {
		u1 = targetGraph.Graph->new_node(na);
		targetGraph.Nodes[na.id]= u1;
	}
	else
		u1 = targetGraph.Nodes[na.id];

	return u1;
}

//
// Copies edge from sourceGraph to destinationGraph.
// Accessibility: private
//
static inline void copyEdge(leda_edge sourceEdge, ColoredGraph& sourceGraph, ColoredGraph& targetGraph)
{
	leda_node u2 = copyNode(sourceGraph.Graph->source(sourceEdge),sourceGraph, targetGraph);
	leda_node v2 = copyNode(sourceGraph.Graph->target(sourceEdge),sourceGraph, targetGraph);

	Color c = sourceGraph.Graph->inf(sourceEdge).color;
	leda_edge e;
	forall_adj_edges(e, u2)
	{
		if(targetGraph.Graph->target(e) == v2 && targetGraph.Graph->inf(e).color == c) //edge already exists in target graph
			return;
	}

	EdgeAttr edgeAttr = sourceGraph.Graph->inf(sourceEdge);
	targetGraph.Graph->new_edge(u2, v2, edgeAttr);
	if(!targetGraph.directed) targetGraph.Graph->new_edge(v2, u2, edgeAttr);;
}

//
// Sampling to estime the connectivity between s1 and the remaining elements of S in G with the given colors.
// Accessibility: private
//
static inline long relG_col(ColoredGraph& G, LedaNodeSet& S, const leda_node s1,
				   Color cl, const h_array <long, int>* base_colors, random_source& rand,
				   bool anyColor, int loopCount) 
{
	h_array <leda_node, long> S1;
	slist<leda_node> nodes1, nodes2;
	slist<leda_node> &old_n = nodes1, &new_n = nodes2;

	int count = 0;
	for(int i = 1; i <= loopCount; i++) 
	{
		decrementRunId(G);
		
		S1.clear();
		leda_node n;
		forall_defined(n, S) 
		{
			S1[n] = 1;
		}
		S1.undefine(s1);

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

							if(S1.defined(v)) 
							{
								S1.undefine(v);
								if(S1.empty())
								{
									terminate = true;
									count++;
									goto exitwhile;
								}
							}
						}
					}
				}
			}

			slist<leda_node> &tmp = old_n;
			old_n = new_n;
			new_n = tmp;
		} // end while
exitwhile:;
	}

	return count;
}

//
// Sampling to estime the connectivity between the elements of S in G with the given colors.
// Accessibility: private
//
static inline long relG_col(ColoredGraph& G, LedaNodeSet& S,
				   Color cl, const h_array <long, int>* base_colors, random_source& rand,
				   bool anyColor, int loopCount) 
{
	if(!G.directed)
	{
		const leda_node s1 = S.key(S.first_item());
		return relG_col(G, S, s1, cl, base_colors, rand, anyColor, loopCount);
	}
	else
	{
		long max_count = 0;
		leda_node s1;
		forall_defined(s1, S)
		{
			long count = relG_col(G, S, s1, cl, base_colors, rand, anyColor, loopCount);
			max_count = std::max(count, max_count);
		}
		return max_count;
	}
}

//
// Implementation of top-K colors SteinerTree algorithm.
// Accessibility: private
// Parameters:
//	  G: the graph.
//    T: a set of pre selected Steiner trees
//	  S: the graph nodes where the connectivity should be maximized. Must be element of G.
//    k: the number of colors to return.
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
//
static inline ColorList* selectTopKColorSteinerTrees(ColoredGraph& G, SteinerTreeSet& T, LedaNodeSet& S, int k)
{
	random_source rand(1,10);

	// init selected trees graph for montecarlo
	ColoredGraph Gc1;
	RelationGraph G1;
	Gc1.Graph = &G1;
	Gc1.directed = G.directed;

	LedaNodeSet S1;
	leda_node nS; forall_defined(nS, S) 
		S1[G1.new_node(G.Graph->inf(nS))] = 1; //init ledaNodeSet S with nodes from G1

	h_array<long, long> T_remain;
	long i; forall_defined(i, T) T_remain[i] = 1;

	slist<int> T_selected;

	ColorSet L1; // colors for fast lookup
	ColorList* CL1 = new ColorList(); // colors in order to return

	bool terminate = false;
	L1.clear();
	while(L1.size() < k && !terminate)
	{
		long t_index;
		slist<int> remainingIndexes;
		forall_defined(t_index, T_remain) remainingIndexes.append(t_index); //copy index list so we may modify T_remain in the next loop
		
		forall(t_index, remainingIndexes)
		{
			SteinerTree& t = T[t_index];
			
			ColorSet L_t = t.colors.join(L1); //current colors
			if(L_t.size() - L1.size() == 0) //t doesn't add new colors to L1 -> add to selected
			{
				T_selected.push(t_index);
				T_remain.undefine(t_index);
			}
			else if(L_t.size() > k)
			{
				T_remain.undefine(t_index);
			}
			L_t.clear();
		}
		remainingIndexes.clear();

		double prob_opt = -1;
		SteinerTree* t_opt = NULL;
		forall_defined(t_index, T_remain)
		{
			SteinerTree& t = T[t_index];

			T_selected.push(t_index);

			//calc RCol_T start
			G1.del_all_edges();
			int i;
			forall(i, T_selected)
			{
				SteinerTree& t = T[i];
				leda_edge e;
				forall(e, t.edges)
				{
					copyEdge(e, G, Gc1);
				}
			}

			// actual MC Simulation
			double prob = relG_col(Gc1, S1, NEUTRAL_COLOR, NULL, rand, true, MONTECARLO_LOOP);
			//calc RCol_T end

			T_selected.pop();

			if(prob > prob_opt)
			{
				prob_opt = prob;
				t_opt = &t;
			}
		}

		if(t_opt != NULL)
		{
			T_selected.push(t_opt->id);
			T_remain.undefine(t_opt->id);

			Color c;
			forall(c, t_opt->colors)
			{
				if(!L1.member(c))
				{
					L1.insert(c);
					CL1->append(c);
				}
			}
		}
		else
		{
			terminate = true;
		}
	}

	//cleanup and release memory
	G1.clear();
	Gc1.Nodes.clear();

	S1.clear();
	L1.clear();
	T_selected.clear();
	T_remain.clear();
	return CL1;
}

//
// Implementation of heuristic for MAX connectivity.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    k: the number of colors to return
//    r: the number of top-reliable paths to use.
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
//
ColorList* maxConn(ColoredGraph& G, LedaNodeSet& S, int k, int r) 
{
	SteinerTreeSet T;
	getTopRSteinerTrees(G, S, r, T);
	cout << "trees: " << T.size() << endl;
	ColorList* L_topk = selectTopKColorSteinerTrees(G, T, S, k);

	//cleanup and release memory
	int i;
	forall_defined(i, T)
	{
		SteinerTree& t = T[i];
		t.colors.clear();
		t.edges.clear();
	}
	T.clear();
	return L_topk;
}

//
// Implementation of baseline1 for MAX connectivity.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    k: the number of colors to return
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
//
ColorList* maxConn_baseline1(ColoredGraph& G, LedaNodeSet& S, int k) 
{
	clock_t timeLimit = start + CLOCKS_PER_SEC * 1000000; //set time limit sec
	bool abort = false;

	random_source rand(1,10);
	p_queue <double, long> sort_colors;
	sort_colors.clear();

	long cl;
	forall(cl, G.Colors) 
	{
		abort = clock() > timeLimit;
		if(abort) break;
		long count = relG_col(G, S, cl, NULL, rand, false, BASELINE1_LOOP);
		sort_colors.insert(-count, cl);
	}

	ColorList* L = new ColorList();
	for(int i = 1; i <= k && !abort; i++) 
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
// Implementation of baseline2 for MAX connectivity.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    k: the number of colors to return
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected.
//
ColorList* maxConn_baseline2(ColoredGraph& G, LedaNodeSet& S, int topk) 
{
	clock_t timeLimit = start + CLOCKS_PER_SEC * 2000; //time limit 1000 sec
	bool abort = false;

	ColorList* L = new ColorList();
	random_source rand(1,10);

	p_queue <double, long> sort_colors;
	sort_colors.clear();

	h_array <long, int> colors_del, colors_add;
	colors_del.clear();
	colors_add.clear();   

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
			abort = clock() > timeLimit;
			if(abort) break;
			long count = relG_col(G, S, cl, &colors_add, rand, false, BASELINE2_LOOP);
			sort_colors.insert(-count, cl);
		}

		if(abort) continue;
		pq_item min_col = sort_colors.find_min();
		Color c = sort_colors.inf(min_col);
		colors_del.undefine(c);
		colors_add[c] = 1;

		double prio = sort_colors.prio(min_col);
		L->append(c);
		cout << c << " p:" << -prio << endl;
	}

	if(abort) L->append(-1);
	return L;
}

//
// Implementation of the result evaluation for MAX connectivity. Uses sampling to estimate the result, 
// hence the result value is not deterministic.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids. All ids must be element of G
//    selectedColors: the set of selected colors to evaluate
// Returns: The evaluation metric for the selected colors. Value is between 0 and EVAL_LOOP (see types.h).
//
double maxConn_evaluateMetric(ColoredGraph& G, LedaNodeSet& S, const h_array <long, int>& selectedColors) 
{
	random_source rand(1,10);
	long metric = relG_col(G, S, NEUTRAL_COLOR, &selectedColors, rand, false, EVAL_LOOP);
	return metric;
}
