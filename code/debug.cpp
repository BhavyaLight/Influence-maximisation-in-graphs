//
// Project: Top-k Reliable Edge Colors in Uncertain Graphs
//
// Author:	Andreas Nufer
// E-Mail:	anufer@student.ethz.ch
// Date:    2015-07-01
// Summary: Contains "toString()" functions for various data structures used for debugging
//

#include "types.h"

//
// Returns a human readable string representation of the given ColorSet object.
// Accessibility: public
//
std::string printColorSet(ColorSet& cs)
{
	stringstream ss;
	long l;
	forall(l, cs)
	{
		ss << " " << l ;
	}

	return ss.str();
}

//
// Returns a human readable string representation of the given Path object.
// Accessibility: public
//
std::string printPath(Path& p, ColoredGraph& G1)
{
	RelationGraph& G = *G1.Graph;
	std::stringstream ss;
	NodeAttr ns = G.inf(G.source(p.inf(p.first_item())));//(NodeAttr&)p.inf(p.first_item())->terminal(0)->data(0);
	ss << ns.id;
	
	double totalCost = 0;

	leda_edge e;
	forall(e, p)
	{
		double cost = 0;
		EdgeAttr ea = G.inf(e);
		if(G1.logCapacity)
			cost += ea.capacity / 1000;
		else
			cost += -log(ea.capacity);

		totalCost += cost;
		NodeAttr n = G.inf(G.target(e));// (NodeAttr&)e->terminal(1)->data(0);
		ss << "-(c:" << ea.color <<",cost:" << cost << ")->" << n.id;
	}

		ss << " path-cost: " << totalCost << " ";
	return ss.str();
}

//
// Returns a human readable string representation of the given PathSet object.
// Accessibility: public
//
std::string printPathSet(PathSet& ps)
{
	ColoredGraph G1;
	RelationGraph G;
	G1.Graph = &G;
	stringstream ss;
	Path p;
	forall(p, ps)
	{
		ss << "{" << printPath(p, G1) << "} " ;
	}

	return ss.str();
}

//
// Returns a human readable string representation of the given RelationGraph object.
// Accessibility: public
//
std::string printGraph(RelationGraph& G)
{
	std::stringstream ss;
	
	leda_edge e;
	forall_edges(e, G)
	{
		NodeAttr u = G.inf(G.source(e));
		NodeAttr v = G.inf(G.target(e));
		const EdgeAttr& ea = G.inf(e);
		ss << "(" << u << "," << v << " | c:" << ea.color << " p:" << ea.capacity << ") ";
	}

	leda_node n;
	forall_nodes(n,G)
	{
		NodeAttr na = G.inf(n);
		if(na.cost >= 0)
			ss << "[" << na.id << " | cost:" << na.cost << "] ";
	}
	return ss.str();
}

//
// Returns a human readable string representation of the given h_array object.
// Accessibility: public
//
std::string printHashtable(h_array<leda_node, leda_edge>& hmap)
{
	std::stringstream ss;

	h_array<leda_node, leda_edge>::item it;
	forall_items(it, hmap)
	{
		ss << "(" << hmap.key(it) << " -> (" << hmap.inf(it) << ") " << endl;
	}
	return ss.str();
}

//
// Returns a human readable string representation of the given edge object.
// Accessibility: public
//
std::string printEdge(const ColoredGraph& G, leda_edge edge)
{
	std::stringstream ss;

	NodeAttr u = G.Graph->inf(G.Graph->source(edge));
	NodeAttr v = G.Graph->inf(G.Graph->target(edge));
	const EdgeAttr ea = G.Graph->inf(edge);
	ss << "(" << u << "," << v << " | c:" << ea.color << " p:" << ea.capacity << ") ";
	return ss.str();
}

//
// Returns a human readable string representation of the given edge list object.
// Accessibility: public
//
std::string printEdges(slist<leda_edge> edges, const ColoredGraph& G1)
{
	RelationGraph& G = *G1.Graph;
	std::stringstream ss;

	leda_edge e;
	forall(e, edges)
	{
		NodeAttr u = G.inf(G.source(e));
		NodeAttr v = G.inf(G.target(e));
		const EdgeAttr& ea = G.inf(e);
		ss << "(" << u << "," << v << "|c:" << ea.color << " p:" << ea.capacity << ") ";
	}

	return ss.str();
}