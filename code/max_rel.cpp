//
// Project: Top-k Reliable Edge Colors in Uncertain Graphs
//
// Author:	Bhavya Chandra
// E-Mail:	chan0727@e.ntu.edu.sg
// Date:    2016-10-11
// Summary: Contains the specific implementations for the MAX-REL optimization. Heuristic and baselines are the same as for MAX-sum. 
//          See max_sum.cpp for the remaining implementations.
//

#include "functions.h"

// private constants
const long NEW_S = -1;
const long NEW_T = -2;

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
// Adds node u to G1 and registers reverse lookup in M1
// Accessibility: public
//
inline leda_node addNode(RelationGraph& G1, Node u, AddressMap& M1)
{
	leda_node u1;
	if(!M1.defined(u)) {
		u1 = G1.new_node(NodeAttr(u));
		M1[u]= u1;
	}
	else
		u1 = M1[u];

	return u1;
}

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
	// leda_node t1 = Gr.new_node(NodeAttr(NEW_T));
	// G1->Nodes[NEW_T] = t1;

	std::vector<leda_edge> newEdges;

	Node s;
	forall(s, S) 
	{
		EdgeAttr edgeAttr = EdgeAttr(newCapacity, NEUTRAL_COLOR);
		newEdges.push_back(Gr.new_edge(s1, G1->Nodes[s], edgeAttr));
	}

	// Node t;
	// forall(t, T) 
	// {
	// 	EdgeAttr edgeAttr = EdgeAttr(newCapacity, NEUTRAL_COLOR);
	// 	newEdges.push_back(Gr.new_edge(G1->Nodes[t], t1, edgeAttr));
	// }

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
	// leda_node t1 = G1.Nodes[NEW_T];

	G1.Graph->del_node(s1);
	// G1.Graph->del_node(t1);
	G1.Nodes[NEW_S] = NULL;
	// G1.Nodes[NEW_T] = NULL;
}

//
// Implementation of the objective function to find sum of probabilities of each target.
// Accessibility: public
// Parameters: 
//			Colored Graph G1
//			NodeSet Target T1
// Returns: An estimation of the sum of probabilities of reaching each target node
double sumOfTargetNodeProbabilities(ColoredGraph& Gc1,NodeSet& Targets, AddressMap& M1){
	random_source S(1,10);
	double count=0,sum=0;
	cout<<"DOing sum monte carlo"<<endl;

    Node t;
	forall(t,Targets){
	count = monteCarlo(M1[NEW_S], M1[t], S, Gc1);
	sum+=count;
	}
	cout<<"END doing sum monte carlo. Summed value:"<<sum<<endl;
	return sum;
}

//
// Implementation of top-K colors algorithm for multiple targets.
// Accessibility: public
// Parameters:
//    P: a set of pre selected paths
//	  G: the graph.
//	  s: the id of the source node. Must be element of G.
//    t: the set of target nodes
//    k: the number of colors to return.
//    p: Return type. The sum of the sampling counts of the selected colors.
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
ColorList* const topKColorsOfMaxTarget(PathSet& P, ColoredGraph& G, Node s, NodeSet& T, int k, double* p)
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
	// G1 graph contains the present selected paths
	leda_node s_G1 = addNode(G1, s, M1);  //Create a new node S in G, M[s]=new node created in G

	int monteCarloCount = 0;

	ColorList* L1 = new ColorList(); //selected colors to return
	ColorSet selected_colors; //for fast lookup
	
	bool terminate = false;
	while(topr_paths.size()!=0 && terminate != true) 
	{
		cout<<"PATHS THAT REMAIN"<<endl;
		cout<<topr_paths<<endl;
		list<long>::item it;
		//For each path in the list of top r paths
		//Delete paths that have the colors that are already selected_colors
		//Delete the paths that exceed the color limit 
		forall_items(it, topr_paths) 
		{
			//Get the path
			int i = topr_paths.contents(it);
			Path& path = P[i];
			//Add all colors
			ColorSet current_colors = getColors(path, G).join(selected_colors);
			if(current_colors.size() - selected_colors.size() == 0) // current path doesn't add any new colors -> add to selected paths immediately
			{
				topr_paths.del_item(it);
                
				// Add the path to G1
				// Construct the path from source node s_G1 through all edges in the path till the target edge
				leda_node u_G1 = s_G1;
				leda_edge e;
				Node v_G;
				forall(e, path)
				{
					EdgeAttr e_G = G.Graph->inf(e);
					if(G.logCapacity) e_G.capacity = exp(-e_G.capacity / 1000);
					v_G = G.Graph->inf(G.Graph->target(e)).id;	//Return id of the target node of edge e
					leda_node v_G1 = addNode(G1, v_G, M1); //add node v_G to G, M1[v_G]= the newly created node in G

					G1.new_edge(u_G1, v_G1, e_G);     
					if(!G.directed) G1.new_edge(v_G1, u_G1, e_G);
					u_G1 = v_G1;
				}
				pNeedsRefresh = true;
			}
			// if it is >k, simply delete the path 
			if(current_colors.size() > k)
			{
				topr_paths.del_item(it);
			}
		}
		cout<<"MAIN AEAS ENTERED"<<endl;
		list<long>::item best_path_item = NULL;
		Node visited_target;
		double best_path_sum = -1;//any negative value will work

		// At the present stage, finds the best path that will be added to get max reliability
		//colors of any path in topr_paths combined with selected colors is guaranteed to be lower than k.
        forall_items(it, topr_paths) 
		{

			int i = topr_paths.contents(it);
			Path& path = P[i];
			
			leda_node u_G1 = s_G1;

			leda::list<leda_edge> edgeList;
			leda_edge e;
			Node v_G;
			// Add each edge of the current path in graph G1
			forall(e, P[i])
			{
				EdgeAttr e_G = G.Graph->inf(e);
				if(G.logCapacity) e_G.capacity = exp(-e_G.capacity / 1000);
				v_G = G.Graph->inf(G.Graph->target(e)).id;	
				leda_node v_G1 = addNode(G1, v_G, M1); //add node to G1

				edgeList.push(G1.new_edge(u_G1, v_G1, e_G));
				if(!G.directed) edgeList.push(G1.new_edge(v_G1, u_G1, e_G));

				u_G1 = v_G1;
			}

			monteCarloCount++;
			double sum = sumOfTargetNodeProbabilities(Gc1,T,M1);
			if(sum > best_path_sum)
			{
				best_path_sum = sum;
				best_path_item = it;

			}

			*p = std::max(*p, sum);

			// remove edges from current path from G1
			G1.del_edges(edgeList);
			edgeList.clear();
		}
		cout<<"ADDING PATHS TO GRAPH"<<endl;
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
			terminate = true; //paths finished

    } //end while
	
	// if(pNeedsRefresh) //refresh p when there were path selected without doing a monteCarlo estimation. //?? 
	// {
	// 	Node t;
	// 	forall(t,T){
	// 	monteCarloCount++;
	// 	double cost = monteCarlo(M1[NEW_S], M1[t], S, Gc1); //monte carlo sampling cost estimation
	// 	std::max(*p, cost);
	// 	}
	// }
	cout<<"RETURNING TIME"<<endl;
	return L1;
}

//
// Implementation of heuristic for MAX reliability in multiple source and targets.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    k: the number of colors to return
//    r: the number of top-reliable paths to use.
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
//
ColorList* maxRel(ColoredGraph& G, NodeSet& S, NodeSet& T, int k, int r) 
{
	//add new s1 and t1 to graph
	double newCapacity = G.logCapacity ? 0 : 1;
	std::vector<leda_edge> newEdges = extendGraph(&G, S, T, newCapacity);
    
	// h_array<Node,PathSet> topRPaths;
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
				// newPath.append(getEdge(t, NEW_T, G));

				i++;
				combinedPaths[i] = newPath;
				cout<<"Index:"<<i<<" PATH: "<<printPath(newPath,G)<<endl;
			}
			P.clear();
		}
		// topRPaths[t]=combinedPaths;
		// combinedPaths.clear();
	}

	//get top k-colors
	double p = 0;
	// find the top K colors for reaching the target set from source T
	ColorList* L_top = topKColorsOfMaxTarget(combinedPaths, G, NEW_S, T, k, &p);

	// clean up path set
	// PathSet P;
	// forall(P, topRPaths){
		long j;
		forall_defined(j,combinedPaths){
			combinedPaths[j].clear();
		}
	// }

	combinedPaths.clear();
	
	//revert changes in graph
	reduceGraph(G, newEdges);
	return L_top;
}


// Implementation of the result evaluation for MAX rel. Uses sampling to estimate the result, 
// hence the result value is not deterministic.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    selectedColors: the set of selected colors to evaluate
// Returns: The evaluation metric for the selected colors. Value is between 0 and EVAL_LOOP (see types.h).
//
double maxRel_evaluateMetric(ColoredGraph& G, NodeSet S, NodeSet T, const h_array <long, int>& selectedColors) 
{
	std::vector<leda_edge> newEdges = extendGraph(&G, S, T, 1);

	h_array <long, int> selectedColorsPlusNeutral(selectedColors);
	selectedColorsPlusNeutral[NEUTRAL_COLOR] = 1;

	double metric = single_evaluateMetric(G, NEW_S, NEW_T, selectedColorsPlusNeutral);

	reduceGraph(G, newEdges);
	return metric;
}