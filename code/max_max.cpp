//
// Project: Top-k Reliable Edge Colors in Uncertain Graphs
//
// Author:	Andreas Nufer
// E-Mail:	anufer@student.ethz.ch
// Date:    2015-07-01
// Summary: Contains the heuristic and baseline implementations for the MAX-maximum optimization. 
//

#include "functions.h"

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
ColorList* maxMax(ColoredGraph& G, NodeSet& S, NodeSet& T, int k, int r) 
{
	double p_top = 0;
	ColorList* L_top = new ColorList();

	long t;
	forall(t, T)  
	{
		long s;
		forall(s, S)
		{
			if(s == t) continue;

			PathSet P;
			topRShortestPaths(G, s, t, r, P);
			double p = 0;
			ColorList* L = topKColors(P, G, s, t, k, &p);
			
			if(p_top <= p) 
			{
				delete L_top;
				p_top = p;
				L_top = L;
			}
			else 
			{
				delete L;
			}
		}
	}

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
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected. The list might contains less than k colors.
//
ColorList* maxMax_baseline1(ColoredGraph& G, NodeSet& S, NodeSet& T, int k) 
{
	double p_top = 0;
	ColorList* L_top = new ColorList();

	long s;
	forall(s, S) 
	{
		long t;
		forall(t, T) 
		{
			if(s == t) continue;

			double p = 0;
			ColorList* L = baseline1(G, s, t, k, &p);
			
			if(p_top <= p) 
			{
				delete L_top;
				p_top = p;
				L_top = L;
			}
			else 
			{
				delete L;
			}
		}
	}
	return L_top;
}

//
// Implementation of baseline2 for MAX maximum.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    k: the number of colors to return
// Returns: An ordered list of colors for the top-k colors problem. The colors are ordered as they were selected.
//
ColorList* maxMax_baseline2(ColoredGraph& G, NodeSet& S, NodeSet& T, int k) 
{
	double p_top = 0;
	ColorList* L_top = new ColorList();

	long s;
	forall(s, S) 
	{
		long t;
		forall(t, T) 
		{
			if(s == t) continue;

			double p = 0;
			ColorList* L = baseline2(G, s, t, k, &p);
			
			if(p_top <= p) 
			{
				delete L_top;
				p_top = p;
				L_top = L;
			}
			else 
			{
				delete L;
			}
		}
	}
	return L_top;
}

//
// Implementation of the result evaluation for MAX maximum. Uses sampling to estimate the result, 
// hence the result value is not deterministic.
// Accessibility: public
// Parameters:
//	  G: the graph data
//	  S: the set of ids of the source nodes. All ids must be element of G
//    T: the set ids of the destination nodes. All ids must be element of G
//    selectedColors: the set of selected colors to evaluate
// Returns: The evaluation metric for the selected colors. Value is between 0 and EVAL_LOOP (see types.h).
//
double maxmax_evaluateMetric(ColoredGraph& G, NodeSet S, NodeSet T, const h_array <long, int>& selectedColors) 
{
	double metric_max = 0;
	long s;
	forall(s, S) 
	{
		long t;
		forall(t, T) 
		{
			if(s == t) continue;
			double metric = single_evaluateMetric(G, s, t, selectedColors);
			metric_max = std::max(metric, metric_max);
		}
	}

	return metric_max;
}
