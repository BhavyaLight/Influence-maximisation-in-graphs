#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


main(){
  char line[50];
  char *node1, *node2;
  int max = 0;
    
  std::ifstream g_stream;
  std::ofstream sout;
  g_stream.open("DBLP_Graph_0.1.txt");
  while (!g_stream.eof()) {
    g_stream.getline(line,50,'\n');
    if (strlen(line) == 0)
      continue;
    else {
      node1 = strtok(line, "_");
      node2 = strtok(NULL, " ");
      if(max<atoi(node1))
        max = atoi(node1);
      if(max<atoi(node2))
        max = atoi(node2);
    }
  }
  std::cout << max << std::endl;
}     
