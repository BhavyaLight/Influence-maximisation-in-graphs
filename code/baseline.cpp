//
// Project: Top-k Reliable Edge Colors in Uncertain Graphs
//
// Author:	Andreas Nufer
// E-Mail:	anufer@student.ethz.ch
// Date:    2015-07-01
// Summary: Contains the baseline implementations for single source/destination
//

#include "functions.h"

//
// Implementation of baseline1 for single source/destination.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  s: the id of the source node. Must be element of G
//    t: the id of the destination node. Must be element of G
//    k: the number of colors to return
//    p: Return type. The sum of the sampling counts of the selected colors
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
//
ColorList* baseline1(ColoredGraph& G, Node s, Node t, int k, double *p) 
{
	clock_t timeLimit = start + CLOCKS_PER_SEC * 1000000; //set time limit in sec
	bool abort = false;
	random_source S(1,10);
	p_queue <double, long> sort_colors;
	sort_colors.clear();

	slist<leda_node> nodes1, nodes2;
	slist<leda_node> &old_n = nodes1, &new_n = nodes2;

	leda_node s1 = G.Nodes[s];
	leda_node t1 = G.Nodes[t];

	long cl;
	forall(cl, G.Colors) 
	{
		int count = 0;
		for(long i = 1; i <= BASELINE1_LOOP && !abort; i++) 
		{
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
							S >> num;
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

		sort_colors.insert(-count, cl);
	}

	*p = 0;
	ColorList* L = new ColorList();
	for(int i = 1; i <= k; i++) 
	{
		pq_item min_color = sort_colors.find_min();
		double prio = sort_colors.prio(min_color);
		*p += -prio;
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
// Implementation of baseline2 for single source/destination.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  s: the id of the source node. Must be element of G
//    t: the id of the destination node. Must be element of G
//    topk: the number of colors to return
//    p: Return type. The sum of the sampling counts of the selected colors
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected.
//
ColorList* baseline2(ColoredGraph& G, Node s, Node t, int topk, double *p) 
{
	clock_t timeLimit = start + CLOCKS_PER_SEC * 1000000; //set time limit in sec
	bool abort = false;
	*p = 0;
	ColorList* L = new ColorList();
	random_source S(1,10);

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
	
	const leda_node s1 = G.Nodes[s];
	const leda_node t1 = G.Nodes[t];
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
								S >> num;
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
			sort_colors.insert(-count, cl);
		}

		pq_item min_col = sort_colors.find_min();
		Color c = sort_colors.inf(min_col);
		colors_del.undefine(c);
		colors_add.insert(c);

		double prio = sort_colors.prio(min_col);
		*p = std::max(*p, -prio);

		L->append(c);
		cout << c << " p:" << -prio << " #" << appendCount << " mx:" << maxn << endl;
		//if(-prio == BASELINE2_LOOP) break;
	}
	if(abort) L->append(-1);
	return L;
}