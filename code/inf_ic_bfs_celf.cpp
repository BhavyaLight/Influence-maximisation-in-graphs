
// ----------------------------------------------------- //
//  Influence Maximization with IC Model+CELF [c++ stl]  //
//                   Arijit Khan                         //
//-------------------------------------------------------//

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <assert.h>
#include <fstream>
#include <algorithm>
#include <unordered_set>
#include <ctime>
#include <queue>
#include <set>

struct NodePrio {
int node;
double priority;
NodePrio() : node(-1), priority(0) {}
NodePrio(int node_, double priority_) : node(node_), priority(priority_) {}
};

class sortNodesByPrio {
public:
    bool operator () (const NodePrio &n1, const NodePrio  &n2)   const;
};

bool sortNodesByPrio::operator () (const NodePrio &n1, const NodePrio &n2) const {
  return n1.priority < n2.priority;
}

/* Input: sourcefile outputfile topK Max
*
*/
int main(int argc, char **argv){

  char sourcefile[256], outfile[256],targetsfile[256],colorsfile[256];
  strcpy(sourcefile, argv[1]);
  strcpy(outfile, argv[2]);
  strcpy(targetsfile, argv[3]);
  strcpy(colorsfile, argv[4]);
  int topk = atoi(argv[5]);  
  int max = atoi(argv[6]);  

  std::vector <std::vector<int> > graph(max + 1, std::vector<int>());
  std::vector <std::vector<int> > capacity(max + 1, std::vector<int>());
  std::set<int> targets;
  std::set<int> colors;
  std::vector <int> neighbors_size; 
  
  int Simulation = 1000;
 
  char line[50000];
  char *s, *s1,*target,*color;
  int wt, o1, o2, id, i;
  srand(time(NULL)); 

  neighbors_size.resize(max+1);   
  for(i=0; i<max; i++)
    neighbors_size[i] = 0;
  
  std::ifstream g_stream;
  std::ofstream sout;

  clock_t start, end, start1, end1;

/*   Read Color File   */
  g_stream.open(colorsfile);
   while(!g_stream.eof()) {
    g_stream.getline(line,500000, '\n');
    if (strlen(line) == 0)
      continue;
    else {
      color = strtok(line," ");
      while(color!=NULL){
      colors.insert(atol(color));
      //std::cout<<"Reading color "<<*color<<std::endl;
      color=strtok(NULL," ");
      }
      
    }
   }

 g_stream.close();

  /*  Read Uncertain Graph */
  /*  To initialise graph vectors with proper size  */
  g_stream.open(sourcefile);
  //g_stream.getline(line,500000, '\n');
  while(!g_stream.eof()) {
    g_stream.getline(line,500000, '\n');
    if (strlen(line) == 0)
      continue;
    else {
      o1 = atol(strtok(line," "));
      o2 = atol(strtok(NULL," "));
      wt = strtod(strtok(NULL, " "),NULL)*100;
      id = atol(strtok(NULL," "));
      if (!colors.empty() && colors.find(id)==colors.end())
        continue;
      neighbors_size[o1]++;
      //neighbors_size[o2]++;
    }
  }
  g_stream.close();

  for(i=0; i<max; i++){
   graph[i].resize(neighbors_size[i]+1);
   capacity[i].resize(neighbors_size[i]+1);
  }
 
  for(i=0; i<max; i++) 
    neighbors_size[i] = 0;

  g_stream.open(sourcefile);
  //g_stream.getline(line,500000, '\n');
  while(!g_stream.eof()) {
    g_stream.getline(line,500000, '\n');
    if (strlen(line) == 0)
      continue;
    else {
      o1 = atol(strtok(line," "));
      o2 = atol(strtok(NULL," "));
      wt = strtod(strtok(NULL, " "),NULL)*100;
      id = atol(strtok(NULL," "));
      
      if (!colors.empty() && colors.find(id)==colors.end())
        continue;
      i = neighbors_size[o1];
      graph[o1][i] = o2;
      capacity[o1][i] = wt;
      neighbors_size[o1]++;      

      //i = neighbors_size[o2];
      //graph[o2][i] = o1;
      //capacity[o2][i] = wt;
      //neighbors_size[o2]++;
    }
  }
  g_stream.close();
   
  // Reading Graph 
  // int j;
  // for(i=0; i<graph.size(); i++){
  //   for(j=0; j<graph[i].size(); j++){
  //     std::cout<<"Node "<<i<<"-- Node"<<graph[i][j]<<"---> Colour"<<j<<" Wt: "<<capacity[i][j]<<std::endl;
  //   }
  // }

//Targets
  g_stream.open(targetsfile);
   while(!g_stream.eof()) {
    g_stream.getline(line,500000, '\n');
    if (strlen(line) == 0)
      continue;
    else {
      target = strtok(line," ");
      while(target!=NULL){
      targets.insert(atol(target));
      //std::cout<<"Reading target "<<*target<<std::endl;
      target=strtok(NULL," ");
      }
      
    }
   }

 g_stream.close();



  // // typedef std::set<int>::iterator setIterator;
  // // setIterator it;

  // // for(it=targets.begin(); it != targets.end(); it++){
  // //   std::cout<<*it<<" "<<std::endl;
  // // }

  // // for(it=colors.begin(); it != colors.end(); it++){
  // //   std::cout<<*it<<" "<<std::endl;
  // // }

  /*  Find Best Seed Nodes */

  std::priority_queue <NodePrio, std::vector<NodePrio>, sortNodesByPrio> count;
  std::vector <int> curIds;
  std::vector <int> nextIds;
  std::vector <int> seeds;
  std::vector <int> flag;
  seeds.resize(topk+1);
  flag.resize (max+1);  

  for(i=0; i<topk+1; i++)
    seeds[i] = -1;

  int  j, k, n, m, p, o3, recompute, num;
  double total, cur;
  sout.open(outfile);

  start = clock();

  /* First Iteration to Find First Seed Node */
  
  i = 0;
  for(n=0; n<max; n++) {
    total = 0;
    // ------- MC sampling to get u.mg
    for(j=1; j<=Simulation; j++) {
      start1 = clock();
      curIds.clear();  
      nextIds.clear(); 
      std::vector <bool> visited(max + 1, false);

      curIds.push_back(n);
      visited[n] = true;
      if(targets.find(n)!=targets.end())
      total++;

      while (!curIds.empty()) {
        //std::cout << "depth" << std::endl;
        for (m=0; m<curIds.size(); m++) {
          o1 = curIds[m];
          //if (visited[o1]) continue;
          //visited[o1] = true;
          for (p=0; p<graph[o1].size(); p++) {
            o2 = graph[o1][p];
            if (visited[o2]) continue;
            num = (rand()/double(RAND_MAX))*1000;
            if(num < capacity[o1][p]) {
              nextIds.push_back(o2);
              visited[o2] = true;
              if(targets.find(o2)!=targets.end())
              total++;
            }
          }
        }
        curIds.clear();
        curIds.swap(nextIds);
      }
      end1 = clock();
      visited.clear();
      //std::cout << "Simulation " << (double) (end1-start1)/CLOCKS_PER_SEC << std::endl;          
    }
    //sout << "calculate " << i+1 << " " << n << " " << total/Simulation << std::endl;
    //------------ Update PQ and flag
    count.push(NodePrio(n, (double) total/Simulation));  
    flag[n] = i;
  }



  sout << i+1 << " " << (count.top()).node << " " << (count.top()).priority << std::endl;
  cur =  (count.top()).priority;
  seeds[i]=(count.top()).node;
  count.pop();

  i = 1;
  recompute = 0;
  while(i<topk) {
    o1 = (count.top()).node;
    total = (count.top()).priority;
    count.pop();
    if(flag[o1] == i) {
      sout << i+1 << " " << o1 << " " << total+cur << " " << recompute << std::endl;  
      seeds[i] = o1;
      cur = cur+total;
      i++;
      recompute = 0;
    }
    else {
      recompute++;
      //std::cout << "recompute " << i+1 << " " << o1 << " " << total << " " << flag[o1] << std::endl;
      total = 0;
      for(j=1; j<=Simulation; j++) {
        curIds.clear();
        nextIds.clear();
        std::vector <bool> visited(max + 1, false);

        for(k=0; k<i; k++) {
          o2 = seeds[k];
          curIds.push_back(o2);
          if(targets.find(o2)!=targets.end())
          total++;
          visited[o2] = true;
        }
        curIds.push_back(o1);
        visited[o1] = true;
        if(targets.find(o1)!=targets.end())
        total++;

        while (!curIds.empty()) {
          //std::cout << "depth" << std::endl;
          for (m=0; m<curIds.size(); m++) {
            o2 = curIds[m];
            for (p=0; p<graph[o2].size(); p++) {
              o3 = graph[o2][p];
              if (visited[o3]) continue;
              num = (rand()/double(RAND_MAX))*1000;
              if(num < capacity[o2][p]) {
                nextIds.push_back(o3);
                visited[o3] = true;
                if(targets.find(o3)!=targets.end())
                total++;
              }
            }
          }
          curIds.clear();
          curIds.swap(nextIds);
        }
        visited.clear();
      }
      //sout << "recompute_new " << i+1 << " " << o1 << " " << (double) total/Simulation << " " << flag[o1] << std::endl;
      if(cur < (double) total/Simulation)
        count.push(NodePrio(o1, (double) total/Simulation - cur));
      else
        count.push(NodePrio(o1, 0));
      flag[o1] = i;
    }
  }
    
  end = clock();
  sout << "Time " << (double) (end-start)/CLOCKS_PER_SEC << std::endl;  
  sout.close();      
} 




  