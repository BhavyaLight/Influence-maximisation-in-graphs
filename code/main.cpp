//
// Project: Top-k Reliable Edge Colors in Uncertain Graphs
//
// Author:	Andreas Nufer
// E-Mail:	anufer@student.ethz.ch
// Date:    2015-07-01
// Summary: Contains the entry point of the application and some common global functions.
//

#include "functions.h"
#include <iostream>
#include <string>
#include <vector>
#include <time.h>
#include <LEDA/core/array.h> 
using namespace std;

//global variables
int BASELINE1_LOOP = 1000; //DBLP=1000, FLIXSTER=100
int BASELINE2_LOOP = 100; //DBLP=100, FLIXSTER=20
bool USE_EPPSTEIN = false;
clock_t start;

//
// Necessary operators to use EdgeAttr and NodeAttr within LEDA graph with no useful implementation.
// Accessibility: public
//
std::ostream &operator<<(std::ostream &os, EdgeAttr const &m) 
{ 
	return os << m.capacity;
}
std::istream &operator>>(std::istream &is, EdgeAttr const &m) 
{ 
	return is;
}
std::ostream &operator<<(std::ostream &os, NodeAttr const &m)
{ 
	return os << m.id;
}
std::istream &operator>>(std::istream &is, NodeAttr const &m)
{ 
	return is;
}

//
// Decrements the runId (used to mark visited nodes in sampling) handles underflows.
// Accessibility: public
//
void decrementRunId(ColoredGraph& G)
{
	//int runId = -10;//use run id as indicator of visited node. abuse node id for this. node id cannot be negative, special nodes in max_sum have ID of -1 and -2

	G.runId--;
	if(G.runId > 0)
	{ //after runId underflow
		leda_node v;
		forall_nodes(v, *G.Graph) G.Graph->assign(v, NodeAttr(0)); //reset node ids for baseline2()
		G.runId = -10;
	}
}

//
// Extracts space separated nodes from input line.
// Accessibility: private
//
static inline ColoredGraph const getGraph(std::string graphFile, bool directed, bool convertCapacity)
{
	clock_t start = clock();

	ColoredGraph Gc;
	RelationGraph* G = new RelationGraph();
	Gc.Graph = G;
	Gc.directed = directed;
	Gc.logCapacity = convertCapacity;
	Gc.Nodes.clear();
	G->make_directed(); //Make the graph directed by splitting Adj edge into in and out
	G->clear();  // Makes the graph empty
	Gc.runId = -10; // runId must be smaller than any possible node id

	std::ifstream input(graphFile);
    std::string line;
	map <long, leda_node>& nodes = Gc.Nodes ; //nodes reference to the graphs nodes

	G->clear();
    while(std::getline(input, line)) 
	{
		stringstream ss;
		ss << line;
		
		long u, v; // start/end node of egde		
		ss >> u >> v;// >> p >> c;

		leda_node u1, v1;
		//Create a new node if node doesn't exist'
		if (!nodes.defined(u)) 
		{
			// adds a new node v to G and returns it. 
			//v is inserted in front of (dir=leda::before) or behind (dir=leda::behind) node u
			u1 = G->new_node(NodeAttr(u));
			nodes[u] = u1;
		}
		else
			u1 = nodes[u];

		if (!nodes.defined(v))
		{
			v1 = G->new_node(NodeAttr(v));
			nodes[v] = v1;
		}
		else
			v1 = nodes[v];

		if(u != v) 
		{
			double p; // propability (capacity)
			long c; //color
			
			while(ss >> p >> c)
			{
				if(convertCapacity) p = -log(p) * 1000;
				EdgeAttr edgeAttr = EdgeAttr(p, c);
				G->new_edge(u1, v1, edgeAttr);
				if(!directed) G->new_edge(v1, u1, edgeAttr); //flixster graph is directed

				Gc.Colors.insert(c);
			} 
		}
    }

	clock_t end = clock();
	std::cout << "Reading Done in " << (double)(end-start)/CLOCKS_PER_SEC << "s" << std::endl;

	line.clear();
	input.close();
	return Gc;
}

//
// Helper function to extract space separated nodes from input line.
// Accessibility: private
//
static inline NodeSet extractNodeSet(std::string data)
{
	NodeSet nodes;

	std::stringstream ss(data);
	long node;
	while(ss >> node)
	{
		nodes.insert(node);
	}

	return nodes;
}

//
// Fills L1 with random colors from G until the size of L1 is k.
// Accessibility: private
//
static void addRandomColors(ColorList& L1, ColoredGraph& G, int k)
{
	random_source rand(1,10);
	leda::array<Color> colors(G.Colors.size());
	int index = 0;
	Color col;
	forall(col, G.Colors) 
	{
		colors[index] = col;
		index++;
	}
	colors.permute();

	forall(col, colors) 
	{
		//go through colors and select one if not already in selected colors
		//and not already k colors printed
		if(L1.size() < k && L1.rank(col) == 0) 
		{
			L1.append(col);
		}
	}
}

//
// Evaluates the given query-file and writes the results to the outFile.
// Accessibility: private
//
static int callEval(std::string algo, std::string dataFile, std::string inputFile, std::string outfile, bool directed) 
{
	// read graph data from <datafile>
	ColoredGraph G = getGraph(dataFile, directed, false);

	std::ofstream sout;
	sout.open(outfile);

	int lineCount = 0;
	std::ifstream input(inputFile);
    std::string line;
    while(std::getline(input, line)) 
	{
		// parse query line
		std::stringstream ss(line);
		std::string tokens;
		std::getline(ss, tokens, ';'); //start nodes
	
		NodeSet S = extractNodeSet(tokens);
		std::getline(ss, tokens, ';'); //target nodes
		NodeSet T = extractNodeSet(tokens);
		std::getline(ss, tokens, ';'); //colors
		h_array <long, int> selectedColors;

		std::stringstream ssc(tokens);
		long color;
		while(ssc >> color)
		{
			selectedColors[color] = 1;
		}

		Node n;
		LedaNodeSet S1;
		LedaNodeSet T1;
		bool hasError = false;
		forall(n, S) {
			if(!G.Nodes.defined(n)){ cout << "Skip line. Node is not element of graph: " << n << endl; hasError = true; }
			else S1[G.Nodes[n]] = 1;
		}
		
		forall(n, T) {
			if(!G.Nodes.defined(n)){ cout << "Skip line. Node is not element of graph: " << n << endl; hasError = true; }
			else T1[G.Nodes[n]] = 1;
		}
		if(hasError)
		{
			S.clear();
			T.clear();
			S1.clear();
			T1.clear();
			lineCount++;
			continue;
		}

		clock_t start = clock();

		// calculate top-k colors
		double metric = 0;
		if(algo == "s") 
		{
			Node s = S.inf(S.first_item());
			Node t = T.inf(T.first_item());
			metric = single_evaluateMetric(G, s, t, selectedColors);
		}
		else if(algo == "sum") metric = maxsum_evaluateMetric(G, S, T, selectedColors); 
		else if(algo == "avg") metric = maxavg_evaluateMetric(G, S, T, selectedColors); 
		else if(algo == "max") metric = maxmax_evaluateMetric(G, S, T, selectedColors); 
		else if(algo == "min") metric = maxMin_evaluateMetric(G, S, T, selectedColors);
		else if(algo == "conn") {
			LedaNodeSet ST; leda_node n1;
			forall_defined(n1, S1) ST[n1] = 1;
			forall_defined(n1, T1) ST[n1] = 1;
			metric = maxConn_evaluateMetric(G, ST, selectedColors);
		} else {
			cout << "Unknown <algo>. Abort!" << endl;
			return -1;
		}
		
		clock_t end = clock();
		double evalTime = (double)(end-start)/CLOCKS_PER_SEC;
		std::cout << "[" << lineCount << "]" << "Eval: " << metric << " " << selectedColors.size() << " in " << evalTime << std::endl;
		// print output to file
		sout << line << "; Eval: " << metric << " " << EVAL_LOOP << "; Evaltime: " << evalTime << std::endl;

		S.clear();
		T.clear();
		S1.clear();
		T1.clear();
		lineCount++;
	}
	sout.flush();
	sout.close();

	return 0;
}

//
// Helper function to find argument with specified prefix in command line parameters.
// Accessibility: private
//
static int findArg(std::string prefix, char* argv[], int argc, int defaultValue)
{
	for (int i = 0; i < argc; i++)
	{
		if(strncmp(prefix.c_str(), argv[i], strlen(prefix.c_str())) == 0)
		{
			std::string arg(argv[i]);
			return atol(arg.substr(prefix.length()).c_str());
		}
	}
	return defaultValue;
}

//
// Entry point. Evaluate command arguments and calls the corresponding algorithms.
//
int main(int argc, char* argv[])
{
	// if too few argurment then print program usage
	if(argc < 2 || (argc < 8 && strcmp(argv[1], "test") != 0))
	{
		std::cout << "Not enough arguments." << endl << endl;
		std::cout << "Usage: colrel <mode> <algo> <datafile> <queryfile> <outfile> <k> <r>" << endl;
		std::cout << "<mode>:      Program mode:" << endl;
		std::cout << "                h     Get colors using heuristic" << endl;
		std::cout << "                b1    Get colors using baseline 1" << endl;
		std::cout << "                b2    Get colors using baseline 2" << endl;
		std::cout << "                eval  Evaluate calculated colors" << endl;
		std::cout << "                test  Run unit tests. data and query file are ignored." << endl;
		std::cout << "<algo>:      Algorithm type:" << endl;
		std::cout << "                s     Single s-t" << endl;
		std::cout << "                sum   Multi s-t with max sum" << endl;
		std::cout << "                max   Multi s-t with max max" << endl;
		std::cout << "                min   Multi s-t with max min" << endl;
		std::cout << "                conn  Multi S with max connectivity" << endl;
		std::cout << "<datafile>:  File containing the graph data" << endl;
		std::cout << "<queryfile>: File containing the queries" << endl;
		std::cout << "<outfile>: File to store the results" << endl;
		std::cout << "<k>: Number of colors" << endl;
	    std::cout << "<r>: Number of top shortest paths" << endl;
		std::cout << "[<baselineloop>]: Number of top shortest paths" << endl;
		return -1;
	}

	// get cmd arguments
	std::string mode = argv[1];
	std::string algo = argv[2];
	std::string dataFile = argv[3];
	std::string queryFile = argv[4];
	std::string outfile = argv[5];
	
	// handle baseline loop count argument
	if(argc >= 9 && (mode == "b1" || mode == "b2"))
	{
		int baselineLoop = atol(argv[8]);
		if(mode == "b1") { BASELINE1_LOOP = baselineLoop; std::cout << "SET BASELINE 1 LOOP: " << BASELINE1_LOOP << endl; }
		if(mode == "b2") { BASELINE2_LOOP = baselineLoop; std::cout << "SET BASELINE 2 LOOP: " << BASELINE2_LOOP << endl; }
	} 
	if(mode == "h" && algo != "conn")
	{
		USE_EPPSTEIN = true;
		std::cout << "Use Eppstein heuristic" << endl;
	}
	if(mode == "d")
	{
		USE_EPPSTEIN = false;
		std::cout << "Use Dijkstra heuristic" << endl;
		mode = "h";
	}

    // Name should contain directed if it is a directed graph
	bool directedGraph = dataFile.find("directed") != std::string::npos;
	if(directedGraph) std::cout << "DIRECTED GRAPH" << endl;

	if(mode == "eval")
	{
		return callEval(algo, dataFile, queryFile, outfile, directedGraph);
	}

	int k = atol(argv[6]);
	int r = atol(argv[7]);

	// read graph data from <datafile>
	
	ColoredGraph G;
	if(algo == "conn" && mode == "h") 
	{
		G = getGraph(dataFile, true, false); // top-r steiner tree algorithm ignores direction on edges
		G.directed = directedGraph;
	}
	else
		G = getGraph(dataFile, directedGraph, USE_EPPSTEIN);

	if(G.Colors.size() < k)
	{
		cout << "Error: Graph contains less than k colors!" << endl;
		exit(-1);
		//throw std::exception("Error: Graph contains less than k colors!");
	}
	std::ofstream sout;
	sout.open(outfile);

	int skipLines = findArg("-skip:", argv, argc, 0);
	int takeLines = findArg("-take:", argv, argc, 1000000);
	int lineCount = -1;
	std::ifstream input(queryFile);
    std::string line;
    while(std::getline(input, line)) 
	{
		lineCount++;
		if(lineCount < skipLines || lineCount > skipLines + takeLines) continue;

		// parse query line
		std::stringstream ss(line);
		std::string nodes;
		std::getline(ss, nodes, ';');
		NodeSet S = extractNodeSet(nodes);
		std::getline(ss, nodes);
		NodeSet T = extractNodeSet(nodes);
		
		start = clock();

		// validate if S and T in graph
		Node n;
		LedaNodeSet S1;
		LedaNodeSet T1;
		bool hasError = false;
		forall(n, S) {
			if(!G.Nodes.defined(n)){ cout << "Skip line. Node is not element of graph: " << n << endl; hasError = true; }
			else S1[G.Nodes[n]] = 1;
		}
		
		forall(n, T) {
			if(!G.Nodes.defined(n)){ cout << "Skip line. Node is not element of graph: " << n << endl; hasError = true; }
			else T1[G.Nodes[n]] = 1;
		}
		if(hasError)
		{
			S.clear();
			T.clear();
			S1.clear();
			T1.clear();
			lineCount++;
			continue;
		}

		// calculate top-k colors
		ColorList* L1 = NULL;
		if(algo == "s") 
		{
			Node s = S.inf(S.first_item());
			Node t = T.inf(T.first_item());
			double p = 0;
			if(mode == "b1") L1 = baseline1(G, s, t, k, &p);
			else if(mode == "b2") L1 = baseline2(G, s, t, k, &p);
			else L1 = single(G, s, t, k, r, &p);
		}
		else if(algo == "sum" && mode == "b1") L1 = maxSum_baseline1(G, S, T, k);
		else if(algo == "sum" && mode == "b2") L1 = maxSum_baseline2(G, S, T, k); 
		else if(algo == "sum") L1 = maxSum(G, S, T, k, r);
		else if(algo == "avg" && mode == "b1") L1 = maxAvg_baseline1(G, S, T, k);
		else if(algo == "avg" && mode == "b2") L1 = maxAvg_baseline2(G, S, T, k); 
		else if(algo == "avg") L1 = maxAvg(G, S, T, k, r);					
		else if(algo == "max" && mode == "b1") L1 = maxMax_baseline1(G, S, T, k);
		else if(algo == "max" && mode == "b2") L1 = maxMax_baseline2(G, S, T, k); 
		else if(algo == "max") L1 = maxMax(G, S, T, k, r); 
		else if(algo == "min" && mode == "b1") L1 = maxMin_baseline1(G, S, T, k);
		else if(algo == "min" && mode == "b2") L1 = maxMin_baseline2(G, S, T, k); 
		else if(algo == "min") L1 = maxMin(G, S, T, k, r); 
		else if(algo == "conn") {
			LedaNodeSet ST; leda_node n1;
			forall_defined(n1, S1) ST[n1] = 1;
			forall_defined(n1, T1) ST[n1] = 1;
			if(mode == "b1") L1 = maxConn_baseline1(G, ST, k);
			else if(mode == "b2") L1 = maxConn_baseline2(G, ST, k); 
			else L1 = maxConn(G, ST, k, r); 
		} else {
			cout << "Unknown <algo>. Abort!" << endl;
			return -1;
		}

		clock_t end = clock();
		
		// print found colors to console
		Color c1;
		forall(c1, *L1)
		{
			std::cout << c1 << " ";
		}
		
		// select arbitrary colors up to k
		int foundColors = L1->size();
		addRandomColors(*L1, G, k);

		// print output to file
		sout << line << ";";
		long l;
		forall(l, *L1)
		{
			sout << " " << l ;
		}
		std::cout << "[" << lineCount << "] Time: " << (double)(end-start)/CLOCKS_PER_SEC << std::endl;

		int loops = 0;
		if(mode == "b1") loops = BASELINE1_LOOP;
		if(mode == "b2") loops = BASELINE2_LOOP;
		sout << "; " << foundColors << "; Loop: " << loops << "; Time: " << (double)(end-start)/CLOCKS_PER_SEC << std::endl;
		delete L1;

		S.clear();
		T.clear();
		S1.clear();
		T1.clear();
	}
	sout.flush();
	sout.close();

	return 0;
}
