#include<stdio.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<limits.h>
#include<math.h>
#include<python2.7/Python.h>
#include<stdbool.h>
#include "functions64.h"

/*
struct result {
  int E;
  int N;
  int* NODES;
  int maxNODE;
};

struct state {
  int s2a;
  int p2a;
};

struct network {
  struct result EN;
  struct activation SIS;
  struct pointers POINT;
  int* NODES;
  int* LINKS;
  int* DEGREE;
  int* WEIGHT;
};

struct pointers {
  int* INI;
  int* FIN;
};

struct activation {
  int Nact;
  int Eact;
  int* NI;
  int** EA;
  int* pointerEA;
};

*/


bool inVECTOR(int k, int size, int* VECTOR){

  int j;
  bool found;

  j=0;
  found = false;
  while((j< size)&&(!found)){
    found = (k == VECTOR[j]);
    j += 1;
  }

  return found;

}

struct result get_E_N(char *filename){

  char str[100]; // 100 is an arbitrary number representing maximum line length to read
  //IMPORTANT: to be large enough to fit two large numbers (~10⁶) in one line
  char *line, *token; // line and token are pointers to strings of characters
  int i, node_i,node_j, count;
  struct result EN; // result = {E,N} with E = # of links, N = # of nodes
  bool found_i, found_j;
  float progress;

  FILE* fp;

  filename=strcat(filename,".txt");
  fp = fopen(filename,"r");

  if (fp == NULL){
    printf("ERROR!!! File path not found!!!\n");
    printf("\n");
    exit(0);
  }

  EN.E = 0;
  line=fgets(str,100,fp);

  printf("...>> getting E\n");
  while (line){ // Canviar line != NULL per line m'ha solucionat un error
    EN.E += 1;
    line = fgets(str,100,fp);
    printf(".\r");
    printf("..\r");
    printf("...\r");
  }

  fclose(fp);

  fp = fopen(filename,"r");

  EN.N = 0;
  line=fgets(str,100,fp);
  EN.NODES = malloc((EN.E+1)*sizeof(int)); //max possible number of nodes is E+1 if there is E links
  i = 0;
  EN.maxNODE = 0;
  printf("...>> getting N\n");
  count = 0;
  while (line != NULL){
    token = strtok(line," "); // Reading the first token
    node_i = atoi(token); // Parsing the char pointer to int
    found_i = inVECTOR(node_i,(EN.E+1),EN.NODES);
    if (!found_i){
      EN.NODES[i] = node_i;
      if (node_i >= EN.maxNODE){
        EN.maxNODE = node_i;
      }
      EN.N += 1;
      i += 1;
    }

    token = strtok(NULL," "); //reading following tokens
    node_j = atoi(token);  // Parsing the char pointer to int
    found_j = inVECTOR(node_j,(EN.E+1),EN.NODES);
    if (!found_j){
      EN.NODES[i] = node_j;
      if (node_j >= EN.maxNODE){
        EN.maxNODE = node_j;
      }
      EN.N += 1;
      i += 1;
    }
    count += 1;
    progress = count*100/(float)EN.E;
    printf("Progress %f %c\r", progress, '%');
    line = fgets(str,100,fp); //Read next line
  }
  fclose(fp);
  printf("\n");

  return EN;
}

void get_pointers(struct network *net, char* filename, int* DEGREE_0){

  // IMPORTANT: this function don't need the network txt to be ordered

  char str[100];
  char* line;
  char* token;
  int node_i, node_j, i, index;

  FILE* fp;
  fp = fopen(filename,"r");

  // Preparing degree for later filling
  for (i=0;i<net->EN.maxNODE;i++){
    DEGREE_0[i] = 0;
  }

  line=fgets(str,100,fp); // Setting first
  net->loops = 0;
  while (line){
    token = strtok(line," "); // Reading the first token
    node_i = atoi(token); // Parsing the char pointer to int
    token = strtok(NULL," "); //reading following tokens
    node_j = atoi(token);  // Parsing the char pointer to int

    if (node_i != node_j){ // This if works to avoid counting loops two times
      DEGREE_0[node_i-1] +=1;
      DEGREE_0[node_j-1] +=1;
    }
    else if(node_i == node_j){
      net->loops += 1;
    }
    line = fgets(str,100,fp); // Read next line
  }
  fclose(fp);


  net->POINT.INI[0] = 0;
  index=1;
  for (i=1;i<net->EN.maxNODE;i++){
    if (DEGREE_0[i-1] != 0){
      net->POINT.INI[index] = net->POINT.INI[index-1] + DEGREE_0[i-1]; //If there are missing nodes, we don't want to sum them
      index +=1;
    }
  }

  for (i=0;i<net->EN.N;i++){
    net->POINT.FIN[i] = net->POINT.INI[i]-1;
  }

}

int cmpfunc_c(const void * a, const void * b){
   return ( *(int*)a - *(int*)b );
 } // Function to be used for qsort()

int get_index(int ind, int* VECTOR, int len){

  int k, middle, first, last;

  first = 0;
  last = len - 1;
  middle = (first+last)/2;

  while (first <= last) {
    if (VECTOR[middle] < ind){
      first = middle + 1;
    }else if (VECTOR[middle] == ind){
      k = middle;
      break;
    }else{
      last = middle - 1;
    }
    middle = (first + last)/2;
  }

  if (first > last){
    printf("ERROR!!!! Not found! %d is not present in the list.\n", ind);
    k = INT_MAX;
    exit(1);
  }

  return k;
}

void last_link(struct network *net, char* filename){

  char str[100];
  char* line;
  char* token;
  int i, j, node_j, node_i, index_i, index_j;
  int count;
  float progress;
  bool found_i, found_j;

  FILE* fp;
  fp = fopen(filename,"r");

  // Counters to keep track of the index for filling LINKS
  for (i=0;i<net->EN.N;i++){
    net->DEGREE[i] = 0;
  }

  net->repetitions = 0;
  count = 0;
  line=fgets(str,100,fp); // Setting first
  while (line != NULL){
    token = strtok(line," "); // Reading the first token
    node_i = atoi(token); // Parsing the char pointer to int
    token = strtok(NULL," "); //reading following tokens
    node_j = atoi(token);  // Parsing the char pointer to int
    index_i = get_index(node_i, net->NODES, net->EN.N);
    if (node_i != node_j){
      index_j = get_index(node_j, net->NODES, net->EN.N);
      // INDEX i
      if (net->POINT.FIN[index_i]>=net->POINT.INI[index_i]){ // there is already something written
        i = net->POINT.INI[index_i];
        found_i = false;
        while((!found_i) && (i<=net->POINT.FIN[index_i])){
          found_i = (node_j == net->LINKS[i]); //look if the link is already in place
          i += 1;
        }
        if(!found_i){
          net->POINT.FIN[index_i] += 1;
          net->LINKS[net->POINT.INI[index_i]+net->DEGREE[index_i]] = node_j;
          net->DEGREE[index_i] +=1;
        }
        else{
          net->repetitions +=1;
        }
      }
      else{
        net->POINT.FIN[index_i] += 1;
        net->LINKS[net->POINT.INI[index_i]+net->DEGREE[index_i]] = node_j;
        net->DEGREE[index_i] +=1;
      }
      // INDEX j
      if (net->POINT.FIN[index_j]>=net->POINT.INI[index_j]){
        j = net->POINT.INI[index_j];
        found_j = false;
        while((!found_j) && (j<=net->POINT.FIN[index_j])){
          found_j = (node_i == net->LINKS[j]); //look if the link is already in place
          j += 1;
        }
        if(!found_j){
          net->POINT.FIN[index_j] += 1;
          net->LINKS[net->POINT.INI[index_j]+net->DEGREE[index_j]] = node_i;
          net->DEGREE[index_j] +=1;
        }
      }
      else{
        net->POINT.FIN[index_j] += 1;
        net->LINKS[net->POINT.INI[index_j]+net->DEGREE[index_j]] = node_i;
        net->DEGREE[index_j] +=1;
      }
    }
    else{
      //If node_i = node_j there is a loop, we just skip it
      //printf("ERROR!!!! Node i is the same as node j!!!!!!!\n");
    }

    count += 1;
    progress = count*100/(float)net->EN.E;
    printf("Progress %f %c\r", progress, '%');
    line = fgets(str,100,fp);
  }
  fclose(fp);
}

bool link_is_once(int* link, int** MATRIX, int len){

  int i, count;
  bool is_once;

  count = 0;
  for(i=0;i<len;i++){
    if(((MATRIX[0][i] == link[0]) && (MATRIX[1][i] == link[1]))
      || ((MATRIX[1][i] == link[0]) && (MATRIX[0][i] == link[1]))){
        count +=1;
      }
  }

  is_once = (count == 1);

  return is_once;
}

bool check_file(char* filename, struct network *net){

  int i, j, node_i, node_j, count;
  char str[100];
  char *line, *token;
  int** MATRIXfile;
  int* link;
  bool sf_sg;
  float progress;

  FILE* fr = fopen(filename, "r");

  line=fgets(str,100,fr); // Setting first

  MATRIXfile = malloc(  net->EN.E*sizeof(int*));
  for(i=0;i<net->EN.E;i++){
    MATRIXfile[i] = malloc(2*sizeof(int));
  }

  // 1) Writing a matrix with the links in the file
  printf("...storing the file in array\n");
  i=0;
  while (line != NULL){
    token = strtok(line," "); // Reading the first token
    node_i = atoi(token); // Parsing the char pointer to int
    token = strtok(NULL," "); //reading following tokens
    node_j = atoi(token);  // Parsing the char pointer to int
    MATRIXfile[0][i] = node_i;
    MATRIXfile[1][i] = node_j;
    line = fgets(str,100,fr);
    i+=1;
    progress = i*100/(float)net->EN.E;
    printf("Progress %f %c\r", progress, '%');
  }

  fclose(fr);

  // 2) Checking each link in the generated vectors are once and only once in the matrix
  printf("...checking links in stored array\n");
  count = 0;
  sf_sg = true;
  for(i=0;i<net->EN.N;i++){
    for(j=net->POINT.INI[i];j<net->POINT.FIN[i]+1;j++){
      link = malloc(2*sizeof(int));
      link[0] = net->NODES[i];
      link[1] = net->LINKS[j];
      sf_sg = link_is_once(link, MATRIXfile, net->EN.E);
      free(link);
      count += 1;
      progress = count*100/(float)(2*net->EN.E); // IMPORTANT: if loops, result won't add to 100
      printf(".\r");
      printf("..\r");
      printf("...\r");
    }
  }
  printf("\n");
  for(i=0;i<net->EN.E;i++){
    free(MATRIXfile[i]);
  }
  free(MATRIXfile);

  return sf_sg;
}

void OverwritePointer(struct network *net, int af,int lenLINKS,int tmp,int neighbor){

  int i,l, index_act, index_act1, index_act2;

  index_act = get_index(net->LINKS[neighbor],net->NODES,net->EN.N);
  for(l=net->POINT.INI[index_act];l<=net->POINT.FIN[index_act];l++){
    if(net->LINKS[l] == af){
      net->SIS.pointerEA[l] = 0;
    }
  }
  index_act1 = get_index(net->SIS.EA[0][tmp-1],net->NODES,net->EN.N);
  for(l=net->POINT.INI[index_act1];l<=net->POINT.FIN[index_act1];l++){
    if(net->LINKS[l] == net->SIS.EA[1][tmp-1]){
      net->SIS.pointerEA[l] = tmp;
    }
  }
  index_act2 = get_index(net->SIS.EA[1][tmp-1],net->NODES,net->EN.N);
  for(l=net->POINT.INI[index_act2];l<=net->POINT.FIN[index_act2];l++){
    if(net->LINKS[l] == net->SIS.EA[0][tmp-1]){
      net->SIS.pointerEA[l] = tmp;
    }
  }
  for(i=0;i<lenLINKS;i++){
    if(net->SIS.pointerEA[i] == net->SIS.Eact){
      net->SIS.pointerEA[i] = 0;
    }
  }
}

void Event_Happens(struct network *LAYERS, float inf, float* COEFFS){

  /********************************************************************************
  function intended to replace modify NI and EA depending on the event:

    - inf = 1 if activation
    - inf = 0 if recovery

  *********************************************************************************/

  int i,j;

  for(i = 0; i<LAYERS[0].SIS.Eact; i++){
    if(LAYERS[0].SIS.EA[0][i] == 0){
      printf("WARNING!!!!!!\n");
      printf("%f\n", inf);
      for(j=0;j<8;j++){
        printf("%f ", COEFFS[j]);
      }
      printf("\n");
      break;
    }
  }

    if((inf <= COEFFS[0]) && (LAYERS[0].SIS.Nact < LAYERS[0].EN.N)){ // viral reactivation

      Reactivation(&LAYERS[0]);

    }
    else if (((COEFFS[0] < inf) && (inf <= COEFFS[1])) && (LAYERS[0].SIS.Nact != 0)){ // deactivation

      Deactivation(&LAYERS[0]);

    }
    else if (((COEFFS[1] < inf)&&(inf <= COEFFS[2])) && ((LAYERS[0].SIS.Nsus != 0)&&(LAYERS[0].SIS.Nact<LAYERS[0].EN.N))){ // m.m. activation

      MMactivation(&LAYERS[0]);

    }
    else if (((COEFFS[2]<inf)&&(inf <= COEFFS[3])) && ((LAYERS[0].SIS.Nsus!=0) && (LAYERS[0].SIS.Nact<LAYERS[0].EN.N))){ // viral activation

      Activation(&LAYERS[0]);

    }
    else if(((COEFFS[3]<inf) && (inf <= COEFFS[4])) &&(LAYERS[1].SIS.Nact < LAYERS[1].EN.N)){ // viral reactivation

      Reactivation(&LAYERS[1]);

    }
    else if (((COEFFS[4] < inf) && (inf <= COEFFS[5])) && (LAYERS[1].SIS.Nact != 0)){ // deactivation

      Deactivation(&LAYERS[1]);

    }
    else if (((COEFFS[5] < inf)&&(inf <= COEFFS[6])) && ((LAYERS[1].SIS.Nsus != 0)&&(LAYERS[1].SIS.Nact<LAYERS[1].EN.N))){ // m.m. activation

      MMactivation(&LAYERS[1]);

    }
    else if (((COEFFS[6]<inf)&&(inf <= COEFFS[7])) && ((LAYERS[1].SIS.Nsus!=0) && (LAYERS[1].SIS.Nact<LAYERS[1].EN.N))){ // viral activation

      Activation(&LAYERS[1]);

    }


}

void initial_population(struct network *net, float fraction, float fitness){

    ///////////////////////////////////////////////////////////////
    // Initial infected population

  int i, j, k, l, index, index_act, index_act1,index_act2, tmp, count;
  int infected;
  bool found2;


  net->SIS.Nact = fraction*net->EN.N;
  net->SIS.Nsus = net->EN.N;
  net->SIS.Eact = 0; // Active Edges
  net->SIS.count.s2a = 0;
  net->SIS.count.p2a = 0;
  net->fitness = fitness;

  for(i=0;i<net->EN.N;i++){
    net->SIS.NI[i] = 0; // I do this bc I don't know what is there in the allocated memory
  }

  for(i=0;i<net->lenLINKS;i++){
    net->SIS.pointerEA[i] = 0; // I do this bc I don't know what is there in the allocated memory
  }
  for(i=0;i<net->EN.E;i++){
    net->SIS.EA[0][i] = 0;
    net->SIS.EA[1][i] = 0;
  }

  for(i=0;i<net->EN.N;i++){
    net->SIS.NS[i] = 0;
    net->SIS.pointNS[i] = 0;
  }


  for(k=0;k<net->SIS.Nact;k++){

    index = rand()%net->EN.N;
    while (inVECTOR(net->NODES[index],net->SIS.Nact,net->SIS.NI)){
      index = rand()%(net->EN.N);
    }

    infected = net->NODES[index];
    net->SIS.NS[index] = 1; // by definition of NS only nodes still non active can be chosen
    net->SIS.Nsus -= 1;
    net->SIS.NI[k] = infected;

    for(j=net->POINT.INI[index];j<=net->POINT.FIN[index];j++){

      index_act = get_index(net->LINKS[j],net->NODES,net->EN.N); // Where is the neighbor in NODES

      if ((net->SIS.pointerEA[j] == 0) && (!inVECTOR(net->LINKS[j],net->SIS.Nact,net->SIS.NI))){ // "Si algun veí de NODES[i] NO està infectat..."
        // Adding an active Edge
        net->SIS.EA[0][net->SIS.Eact] = infected;
        net->SIS.EA[1][net->SIS.Eact] = net->LINKS[j];
        net->SIS.Eact += 1;
        net->SIS.count.s2a += 1;

        net->SIS.pointerEA[j] = net->SIS.Eact;
        for(l=net->POINT.INI[index_act];l<=net->POINT.FIN[index_act];l++){
          if(net->LINKS[l] == infected){
            net->SIS.pointerEA[l] = net->SIS.Eact;
          }
        }
      }
      else if ((net->SIS.pointerEA[j] == 0)&&(inVECTOR(net->LINKS[j],net->SIS.Nact,net->SIS.NI))){
        printf("WARNING! #1 \n");
        continue;

      }
      else if((net->SIS.pointerEA[j] != 0)&&(inVECTOR(net->LINKS[j],net->SIS.Nact,net->SIS.NI))){ // "Si algun veí de NODES[i] està infectat..."

        tmp = net->SIS.pointerEA[j];
        net->SIS.EA[0][tmp-1] = net->SIS.EA[0][net->SIS.Eact-1];
        net->SIS.EA[1][tmp-1] = net->SIS.EA[1][net->SIS.Eact-1];
        net->SIS.pointerEA[j] = 0;

        for(l=net->POINT.INI[index_act];l<=net->POINT.FIN[index_act];l++){
          if(net->LINKS[l] == infected){
            net->SIS.pointerEA[l] = 0;
          }
        }
        index_act1 = get_index(net->SIS.EA[0][tmp-1],net->NODES,net->EN.N);
        for(l=net->POINT.INI[index_act1];l<=net->POINT.FIN[index_act1];l++){
          if(net->LINKS[l] == net->SIS.EA[1][tmp-1]){
            net->SIS.pointerEA[l] = tmp;
          }
        }
        index_act2 = get_index(net->SIS.EA[1][tmp-1],net->NODES,net->EN.N);
        for(l=net->POINT.INI[index_act2];l<=net->POINT.FIN[index_act2];l++){
          if(net->LINKS[l] == net->SIS.EA[0][tmp-1]){
            net->SIS.pointerEA[l] = tmp;
          }
        }
        net->SIS.EA[0][net->SIS.Eact-1] = 0;
        net->SIS.EA[1][net->SIS.Eact-1] = 0;
        for(i=0;i<net->lenLINKS;i++){
          if(net->SIS.pointerEA[i] == net->SIS.Eact){
            net->SIS.pointerEA[i] = 0;
          }
        }

        net->SIS.Eact -= 1;
        net->SIS.count.s2a -= 1;
      }
      else if ((net->SIS.pointerEA[j] != 0)&&(!inVECTOR(net->LINKS[j],net->SIS.Nact,net->SIS.NI))){
        printf("WARNING!!! #2 \n");
        continue;
      }
    }
  }


  found2 = true;
  i=1;
  while((i<=net->SIS.Eact)&&(found2)){
    count = 0;
    for(j=0;j<net->lenLINKS;j++){
      if(net->SIS.pointerEA[j] == i){
        count+=1;
      }
    }
    found2 = (count == 2);
    i+=1;
  }

  if (!found2){
    printf("\n");
    printf("\n");
    printf("ERROR!!!! problems with position %d in EA\n", i-1);
    printf("  [%d %d] is found %d times\n", net->SIS.EA[0][i-1], net->SIS.EA[1][i-1], count);
    printf("\n");
  }
}

void WriteNet_Act(struct network net, char* directory, int filecount){

  int i,j,node, index;
  char pathFile[1000];


  snprintf(pathFile,1000,"%s/net%d.txt",directory, filecount);
  FILE* f = fopen(pathFile,"w");

  int* WRITTEN = malloc(net.SIS.Nact*sizeof(int));
  for(i=0;i<net.SIS.Nact;i++){
    WRITTEN[i] = 0;
  }

  for(i=0;i<net.SIS.Nact;i++){
    node = net.SIS.NI[i];
    index = get_index(node,net.NODES,net.EN.N);
    for(j=net.POINT.INI[index];j<=net.POINT.FIN[index];j++){
      if(((!inVECTOR(net.LINKS[j],net.SIS.Nact,WRITTEN)) && (inVECTOR(net.LINKS[j],net.SIS.Nact,net.SIS.NI)))
      && ((!inVECTOR(node,net.SIS.Nact,WRITTEN)) && (inVECTOR(node,net.SIS.Nact,net.SIS.NI)))){
        //if LINKS[j] is in vector it means the link has been written previously
        // We only print the net of active nodes connected
        fprintf(f,"%d %d\n", node, net.LINKS[j]);
      }
    }
    WRITTEN[i] = node;
  }

  fclose(f);
  free(WRITTEN);
}

void WriteNet_Total(struct network net, char* directory, int filecount){

  int i,j,node, index, index_ngbr, total_nodes;
  char pathFile[1000];


  snprintf(pathFile,1000,"%s/net%d.txt",directory, filecount);
  FILE* f = fopen(pathFile,"w");

  total_nodes = net.EN.N-net.SIS.Nsus;

  int* WRITTEN = malloc(total_nodes*sizeof(int));
  for(i=0;i<total_nodes;i++){
    WRITTEN[i] = 0;
  }

  for(i=0;i<net.SIS.Nact;i++){
    node = net.SIS.NI[i];
    index = get_index(node,net.NODES,net.EN.N);
    for(j=net.POINT.INI[index];j<=net.POINT.FIN[index];j++){
      index_ngbr = get_index(net.LINKS[j],net.NODES,net.EN.N);
      if(((!inVECTOR(net.NODES[index_ngbr],total_nodes,WRITTEN)) && (net.SIS.pointNS[index_ngbr] != 0))
        && (((!inVECTOR(node,total_nodes,WRITTEN)) && (net.SIS.pointNS[index] != 0)))){
        //if LINKS[j] is in vector it means the link has been written previously
        // We only print the net of active nodes connected
        fprintf(f,"%d %d\n", node, net.NODES[index_ngbr]);
      }
    }
    WRITTEN[i] = node;
  }

  fclose(f);
  free(WRITTEN);
}

struct state ComputeEs2A(struct network *net){

  int i;
  int index, neighbor;
  int Es2a, Ep2a;
  struct state online;

  Es2a = 0;
  Ep2a = 0;
  for (i=0;i<net->SIS.Eact;i++){
    neighbor = net->SIS.EA[1][i];
    index = get_index(neighbor,net->NODES,net->EN.N);
    if(net->SIS.NS[index]==0){
      Es2a += 1;
    }
    else{
      Ep2a += 1;
    }
  }

  //using result structure only for convenience
  online.s2a = Es2a;
  online.p2a = Ep2a;

  return online;

}

void Update_pointNS(struct network *net, int cl){

  int i,j,k;

  for(i=0;i<cl;i++){
    if(net[i].SIS.Nsus != 0){
      j = 0;
      for(k=0;k<net[i].EN.N;k++){
        if(net[i].SIS.NS[k] == 0){
          net[i].SIS.pointNS[j] = k;
          j+=1;
        }
      }
    }
  }
}

void Update_pointEA_S(struct network *LAYERS, int cl){

  int i,j,ks,kp,index;

  for(i=0;i<cl;i++){

    ks = 0;
    kp = 0;

    for(j=0;j<LAYERS[i].SIS.Eact;j++){

      index = get_index(LAYERS[i].SIS.EA[1][j], LAYERS[i].NODES, LAYERS[i].EN.N);

      if(LAYERS[i].SIS.NS[index] == 0){

        LAYERS[i].SIS.pointEA_S[ks] = j;
        ks += 1;

      } else if(LAYERS[i].SIS.NS[index] == 2){

        LAYERS[i].SIS.pointEA_P[kp] = j;
        kp += 1;

      }
    }

    LAYERS[i].SIS.count.s2a = ks;
    LAYERS[i].SIS.count.p2a = kp;

    if(LAYERS[i].SIS.count.p2a + LAYERS[i].SIS.count.s2a != LAYERS[i].SIS.Eact){

      printf("\n");
      printf("------------------------------------------------\n");
      printf("WARNING!!!! Not counting active links properly\n");
      printf("\n");
      printf("NET %d \n",i+1);
      printf("\n");
      printf("  Es2a: %d \n", LAYERS[i].SIS.count.s2a);
      printf("  Ep2a: %d \n", LAYERS[i].SIS.count.p2a);
      printf(" _____________\n");
      printf("  Eact: %d \n", LAYERS[i].SIS.Eact);
      printf("------------------------------------------------\n");
      printf("\n");

      exit(1);
    }

  }

}

float ComputeMean(float* VECTOR,int len){

  int i;
  float mean;

  mean = 0;
  for(i=0;i<len;i++){
    mean += VECTOR[i];
  }

  return mean/(float)len;
}

void ResetVectors(struct network *LAYERS, struct network *LAYERS_0, int cl){

  int i,k;

  for (i=0;i<cl;i++){

    LAYERS[i].EN.N = LAYERS_0[i].EN.N;
    LAYERS[i].EN.E = LAYERS_0[i].EN.E;
    LAYERS[i].lenLINKS = LAYERS_0[i].lenLINKS;
    LAYERS[i].SIS.Nact = LAYERS_0[i].SIS.Nact;
    LAYERS[i].SIS.Eact = LAYERS_0[i].SIS.Eact;
    LAYERS[i].SIS.Nsus = LAYERS_0[i].SIS.Nsus;
    LAYERS[i].SIS.Nact = LAYERS_0[i].SIS.Nact;
    LAYERS[i].loops = LAYERS_0[i].loops;
    LAYERS[i].zeros = LAYERS_0[i].zeros;
    LAYERS[i].repetitions = LAYERS_0[i].repetitions;

    for(k=0;k<LAYERS_0[i].EN.N;k++){
      LAYERS[i].NODES[k] = LAYERS_0[i].NODES[k];
    }
    for(k=0;k<LAYERS_0[i].lenLINKS;k++){
      LAYERS[i].LINKS[k] = LAYERS_0[i].LINKS[k];
    }
    for(k=0;k<LAYERS_0[i].EN.N;k++){
      LAYERS[i].DEGREE[k] = LAYERS_0[i].DEGREE[k];
    }
    for(k=0;k<LAYERS_0[i].EN.N;k++){
      LAYERS[i].POINT.INI[k] = LAYERS_0[i].POINT.INI[k];
      LAYERS[i].POINT.FIN[k] = LAYERS_0[i].POINT.FIN[k];
    }
    for(k=0;k<LAYERS_0[i].EN.N;k++){
      LAYERS[i].SIS.NS[k] = LAYERS_0[i].SIS.NS[k];
    }
    for(k=0;k<LAYERS_0[i].EN.E;k++){
      LAYERS[i].SIS.EA[1][k] = LAYERS_0[i].SIS.EA[1][k];
      LAYERS[i].SIS.EA[0][k] = LAYERS_0[i].SIS.EA[0][k];
    }
    for(k=0;k<LAYERS_0[i].EN.N;k++){
      LAYERS[i].SIS.NI[k] = LAYERS_0[i].SIS.NI[k];
    }

    for(k=0;k<LAYERS_0[i].EN.N;k++){
      LAYERS[i].SIS.pointNS[k] = LAYERS_0[i].SIS.pointNS[k];
    }

    for(k=0;k<LAYERS_0[i].lenLINKS;k++){
      LAYERS[i].SIS.pointerEA[k] = LAYERS_0[i].SIS.pointerEA[k];
    }
    for(k=0;k<LAYERS_0[i].lenLINKS;k++){
      LAYERS[i].SIS.pointEA_S[k] = LAYERS_0[i].SIS.pointEA_S[k];
    }
    for(k=0;k<LAYERS_0[i].lenLINKS;k++){
      LAYERS[i].SIS.pointEA_P[k] = LAYERS_0[i].SIS.pointEA_P[k];
    }
  }
}

void InitializeNet(struct network *net, char* file, char* directory){

  int i, j, maxDEGREE, sumDEGREE, Kcount;
  float cummulateP;
  bool same_file;
  char pathFile[4000];
  char *answer;

  answer = malloc(10);

  FILE *fw;

  printf("...getting E & N\n");

  net->zeros = net->EN.maxNODE - net->EN.N; // Count the number of nodes with degree 0

  printf("\n");
  printf("      ---------------------------\n");
  printf("      Number of nodes (N) is %d\n", net->EN.N);
  printf("      Number of links (E) is %d\n", net->EN.E);
  printf("\n");
  printf("      (%d missing nodes)\n", net->zeros); //A missing node is a node with degree 0
  printf("      ---------------------------\n");
  printf("\n");


  /*
  For the last_link I need a vector, sorted from min to max index with the
  different nodes in the network, appearing only once each of them. These nodes
  are already in EN.NODES however EN.NODES is a larger vector, as it is defined
  before knowing the exact number of N. For this reason I define another vector
  of length N for accuracy.
  */

  for (i=0;i<net->EN.N;i++){
    net->NODES[i] = net->EN.NODES[i];
  }
  /*
  newDEGREE contains the same information that DEGREE but without the zeros,
  also it has the repetitions corrected
  */

  qsort(net->NODES, net->EN.N, sizeof(int),cmpfunc_c);

  net->repetitions = 0;
  printf("...filling LINKS\n");
  last_link(net,file);
  printf("\n");
  if (net->repetitions != 0){
    for(i=0;i<net->EN.N;i++){
      net->DEGREE[i] = net->POINT.FIN[i] - net->POINT.INI[i]+1;
    }
  }

  if (net->repetitions != 0){
    printf("\n");
    printf("      -------------------------------------------\n");
    printf("      %d repetitions have been found and cleaned\n", net->repetitions);
    printf("      -------------------------------------------\n");
    printf("\n");
  }

  maxDEGREE = net->DEGREE[0]; // Looking for the max degree
  for(i=1;i<net->EN.N;i++){
    if(net->DEGREE[i] > maxDEGREE){
      maxDEGREE = net->DEGREE[i];
    }
  }

  printf("...writing 'degree.txt'\n");
  net->DEG.P = malloc((maxDEGREE+1)*sizeof(float));
  net->DEG.Pc = malloc((maxDEGREE+1)*sizeof(float));
  net->DEG.P[0] = 0;

  for(i=1;i<=maxDEGREE;i++){
    Kcount = 0;
    for(j=0;j<net->EN.N;j++){
      if(net->DEGREE[j] == i){
        Kcount += 1; // Counting all nodes, that have a degree = i
      }
    }
    net->DEG.P[i] = Kcount/(float)net->EN.N;
  }

  sprintf(pathFile,"%s/degree.txt",directory);
  fw = fopen(pathFile, "w");
  if(fw == NULL){
    printf("Error opening file!\n");
    exit(1);
  }
  for(i=0;i<=maxDEGREE;i++){ // i is the degree
    cummulateP = 0.0;
    for(j=i;j<=maxDEGREE;j++){
      cummulateP += net->DEG.P[j];
    }
    if (net->DEG.P[i] != 0){ //Excluding the degrees with probability 0
      net->DEG.Pc[i] = cummulateP;
      fprintf(fw,"%d %f %f\n", i, net->DEG.P[i], net->DEG.Pc[i]);
    }
  }

  fclose(fw);


  sumDEGREE = 0;
  for (i=0;i<net->EN.N;i++){
    sumDEGREE += net->DEGREE[i];
  }

  printf("...cheking DEGREE matches LINKS dimension\n");
  if (sumDEGREE != (net->lenLINKS-2*net->repetitions)){
    printf("\n");
    printf("ERROR!!! DEGREE does not match size of LINKS !!\n");
    printf("sumDEGREE = %d || lenLINKS - 2*repetitions = %d \n", sumDEGREE, (net->lenLINKS-2*net->repetitions));

  }

  // CHECKING OUTPUT
  printf("\n");
  printf("      3) Do you want to compare current output with initial network file? [Y/N]: ");
  scanf("%s",answer);
  while((strcmp(answer, "Y") != 0)&&(strcmp(answer, "N") != 0)){
    printf("\n");
    printf("Can you repeat please? (remember the uppercase) :");
      scanf("%s",answer);
  }

  if(strcmp(answer, "Y") == 0){
    printf("\n");
    printf("sit back and relax, it can take a while...\n");
    printf(" \n");
    printf("...checking vector output matches initial file\n");
    same_file = check_file(file,net);
    if (same_file){
      printf("\n");
      printf("      CONGRATULATIONS!!! File read and created are the same!!\n");
    }
    else{
      printf(" \n");
      printf("      Too bad :( ... Better luck next time\n");
    }
  }
  else{
    printf(" \n");
  }


  free(answer);

}

void Compute_GllspCoeffs(struct network *LAYERS,int cl, float lambda, float sigma, float mu, float* COEFFS){

  int i;
  int A_0;

  for(i=0;i<4*cl;i++){
    COEFFS[i] = 0;
  }

  A_0 = 0;
  for(i=0;i<cl;i++){
    LAYERS[i].SIS.a_i.TOTAL = LAYERS[i].lambda*LAYERS[i].SIS.count.s2a
                              + LAYERS[i].lambda*LAYERS[i].SIS.count.p2a
                              + sigma*LAYERS[i].SIS.Nact
                              + mu*LAYERS[i].SIS.Nsus;

    A_0 += LAYERS[i].SIS.a_i.TOTAL;

  }


  for (i=0;i<cl;i++){

    LAYERS[i].SIS.a_i.react = LAYERS[i].lambda*LAYERS[i].SIS.count.p2a/(float)A_0;
    LAYERS[i].SIS.a_i.deact = sigma*LAYERS[i].SIS.Nact/(float)A_0;
    LAYERS[i].SIS.a_i.mmact = mu*LAYERS[i].SIS.Nsus/(float)A_0;
    LAYERS[i].SIS.a_i.act = LAYERS[i].lambda*LAYERS[i].SIS.count.s2a/(float)A_0;

    COEFFS[4*(cl-1)*i] = LAYERS[i].SIS.a_i.react;
    COEFFS[4*(cl-1)*i+1] = LAYERS[i].SIS.a_i.deact;
    COEFFS[4*(cl-1)*i+2] = LAYERS[i].SIS.a_i.mmact;
    COEFFS[4*(cl-1)*i+3] = LAYERS[i].SIS.a_i.act;

  }

  for(i=1;i<4*cl;i++){
    COEFFS[i] = COEFFS[i] + COEFFS[i-1];
  }


}

void Compute_Density(struct network *LAYERS, int cl){

  int i;

  for(i=0;i<cl;i++){
    LAYERS[i].rho.act = LAYERS[i].SIS.Nact/(float)LAYERS[i].EN.N;
    LAYERS[i].rho.susc = LAYERS[i].SIS.Nsus/(float)LAYERS[i].EN.N;
    LAYERS[i].rho.pas = 1 - LAYERS[i].rho.act - LAYERS[i].rho.susc;

    if((LAYERS[i].rho.act >1) || (LAYERS[i].rho.act < -0.00001)){
      printf("----------------------------------\n");
      printf("  ERROR in layer %d !!!!!\n", i+1);
      printf("    rho_act =  %f\n", LAYERS[i].rho.act);
      printf("    rho_pas =  %f\n", LAYERS[i].rho.pas);
      printf("    rho_susc =  %f\n", LAYERS[i].rho.susc);
      printf("----------------------------------\n");

      exit(1);
    }
    else if((LAYERS[i].rho.pas > 1) || (LAYERS[i].rho.pas < -0.00001)){
      printf("----------------------------------\n");
      printf("  ERROR in layer %d !!!!!\n", i+1);
      printf("    rho_act =  %f\n", LAYERS[i].rho.act);
      printf("    rho_pas =  %f\n", LAYERS[i].rho.pas);
      printf("    rho_susc =  %f\n", LAYERS[i].rho.susc);
      printf("----------------------------------\n");

      exit(1);
    }
  }

}

void Update_ActRate(struct network *LAYERS, int cl, float lambda){

  int i;

  for(i=0;i<cl;i++){
    LAYERS[i].lambda = lambda*LAYERS[i].weight;

  }
}

void Update_Weight(struct network *LAYERS, int cl, float epsilon, float coupling){

  int i;
  float total_weight;

  if (cl == 2){

    LAYERS[0].fitness = epsilon;
    LAYERS[1].fitness = 1-epsilon;

    total_weight = FI_function(LAYERS[0].rho.act, coupling)*LAYERS[0].fitness
                  + FI_function(LAYERS[1].rho.act, coupling)*LAYERS[1].fitness;

    if (total_weight <= 0){total_weight = 1;} // This parameter is needed for the start if rho_act_i = 0

    LAYERS[0].weight = FI_function(LAYERS[0].rho.act, coupling)*LAYERS[0].fitness/(float)total_weight;
    LAYERS[1].weight = FI_function(LAYERS[1].rho.act, coupling)*LAYERS[1].fitness/(float)total_weight;

  }else {

    printf("\n");
    printf("   >>> FITNESS MODEL NOT APPLIED <<<\n");
    printf("\n");

    total_weight = 0.0;

    for(i=0;i<cl;i++){
      total_weight += FI_function(LAYERS[i].rho.act, coupling);
    }

    if (total_weight <= 0){total_weight = 1.0;}

    for(i=0;i<cl;i++){
      LAYERS[i].weight = (FI_function(LAYERS[i].rho.act, coupling))/(float)total_weight;
    }
  }
}

void Reactivation(struct network *net){

  int j, l, tmp, index, index_inf, index_act;
  int infected;

  // which susceptible link becomes infected [I,S] -> [I,I]?
  index_inf = net->SIS.pointEA_P[rand()%(net->SIS.count.p2a)];
  infected = net->SIS.EA[1][index_inf];
  net->SIS.NI[net->SIS.Nact] = infected;
  net->SIS.Nact += 1;

  // Which new neighbors become susceptible?
  index = get_index(infected,net->NODES,net->EN.N);
  net->SIS.NS[index] = 1;

  for(j=net->POINT.INI[index];j<=net->POINT.FIN[index];j++){

    if(net->SIS.pointerEA[j] == 0){
      net->SIS.EA[0][net->SIS.Eact] = infected;
      net->SIS.EA[1][net->SIS.Eact] = net->LINKS[j];
      net->SIS.pointerEA[j] = net->SIS.Eact+1;

      index_act = get_index(net->LINKS[j],net->NODES,net->EN.N);

      for(l=net->POINT.INI[index_act];l<=net->POINT.FIN[index_act];l++){
        if(net->LINKS[l] == infected){
          net->SIS.pointerEA[l] = net->SIS.Eact+1;
        }
      }
      net->SIS.Eact += 1;

    }

    else if(net->SIS.pointerEA[j] != 0){
      tmp = net->SIS.pointerEA[j];
      net->SIS.EA[0][tmp-1] = net->SIS.EA[0][net->SIS.Eact-1];
      net->SIS.EA[1][tmp-1] = net->SIS.EA[1][net->SIS.Eact-1];
      net->SIS.pointerEA[j] = 0;
      OverwritePointer(net,infected,net->lenLINKS,tmp,j);
      net->SIS.EA[0][net->SIS.Eact-1] = 0;
      net->SIS.EA[1][net->SIS.Eact-1] = 0;
      net->SIS.Eact -= 1;
    }
  }

}

void Deactivation(struct network *net){

  int j, l, tmp, index, index_rec, index_act, recovered;

  index_rec = rand()%(net->SIS.Nact);
  recovered = net->SIS.NI[index_rec];
  net->SIS.Nact -= 1;
  net->SIS.NI[index_rec] = net->SIS.NI[net->SIS.Nact];
  net->SIS.NI[net->SIS.Nact] = 0;

  index = get_index(recovered,net->NODES,net->EN.N);
  net->SIS.NS[index] = 2;

  for(j=net->POINT.INI[index];j<=net->POINT.FIN[index];j++){
    if(net->SIS.pointerEA[j] == 0){
      // !!!! this is the other way around than in the case of activation,
      // here in EA the infected node [0] is which is in LINKS,
      // while the susceptible becomes the recovered
      net->SIS.EA[0][net->SIS.Eact] = net->LINKS[j];
      net->SIS.EA[1][net->SIS.Eact] = recovered;
      net->SIS.Eact += 1;
      net->SIS.pointerEA[j] = net->SIS.Eact;

      index_act = get_index(net->LINKS[j],net->NODES,net->EN.N);
      for(l=net->POINT.INI[index_act];l<=net->POINT.FIN[index_act];l++){
        if(net->LINKS[l] == recovered){
          net->SIS.pointerEA[l] = net->SIS.Eact;
        }
      }
    }
    else if(net->SIS.pointerEA[j] != 0){
      tmp = net->SIS.pointerEA[j];
      net->SIS.EA[0][tmp-1] = net->SIS.EA[0][net->SIS.Eact-1];
      net->SIS.EA[1][tmp-1] = net->SIS.EA[1][net->SIS.Eact-1];
      net->SIS.pointerEA[j] = 0;
      OverwritePointer(net,recovered,net->lenLINKS,tmp,j);
      net->SIS.EA[0][net->SIS.Eact-1] = 0;
      net->SIS.EA[1][net->SIS.Eact-1] = 0;
      net->SIS.Eact -= 1;
    }
  }
}

void MMactivation(struct network *net){

  int index, index_sus, index_inf, index_act, infected;
  int j, l, tmp;

  index_sus = rand()%(net->SIS.Nsus);
  index_inf = net->SIS.pointNS[index_sus];
  infected = net->NODES[index_inf];
  net->SIS.NS[index_inf] = 1; // by definition of NS only nodes still non active can be chosen
  net->SIS.Nsus -= 1;
  net->SIS.pointNS[index_sus] = net->SIS.pointNS[net->SIS.Nsus];
  net->SIS.NI[net->SIS.Nact] = infected;
  net->SIS.Nact += 1;

  index = get_index(infected,net->NODES,net->EN.N); // Which new neighbors become susceptible?
  for(j=net->POINT.INI[index];j<=net->POINT.FIN[index];j++){

    index_act = get_index(net->LINKS[j],net->NODES,net->EN.N);

    if(net->SIS.pointerEA[j] == 0){
      net->SIS.EA[0][net->SIS.Eact] = infected;
      net->SIS.EA[1][net->SIS.Eact] = net->LINKS[j];
      net->SIS.pointerEA[j] = net->SIS.Eact+1;

      for(l=net->POINT.INI[index_act];l<=net->POINT.FIN[index_act];l++){

        if(net->LINKS[l] == infected){
          net->SIS.pointerEA[l] = net->SIS.Eact+1;
        }
      }
      net->SIS.Eact += 1;
    }
    else if(net->SIS.pointerEA[j] != 0){

      tmp = net->SIS.pointerEA[j];
      net->SIS.EA[0][tmp-1] = net->SIS.EA[0][net->SIS.Eact-1];
      net->SIS.EA[1][tmp-1] = net->SIS.EA[1][net->SIS.Eact-1];
      net->SIS.pointerEA[j] = 0;
      OverwritePointer(net,infected,net->lenLINKS,tmp,j);
      net->SIS.EA[0][net->SIS.Eact-1] = 0;
      net->SIS.EA[1][net->SIS.Eact-1] = 0;
      net->SIS.Eact -= 1;
    }
  }
}

void Activation(struct network *net){

  int index_sus, index, index_act, infected;
  int j, l, tmp;

  index_sus = net->SIS.pointEA_S[rand()%(net->SIS.count.s2a)];
  infected = net->SIS.EA[1][index_sus];
  net->SIS.NI[net->SIS.Nact] = infected;
  net->SIS.Nact += 1;

  index = get_index(infected,net->NODES,net->EN.N);
  net->SIS.NS[index] = 1; // by definition of NS only nodes still non active can be chosen

  net->SIS.Nsus -= 1;
  net->SIS.pointNS[index] = net->SIS.pointNS[net->SIS.Nsus];

  // Which new neighbors become susceptible?


  for(j=net->POINT.INI[index];j<=net->POINT.FIN[index];j++){

    index_act = get_index(net->LINKS[j],net->NODES,net->EN.N);

    if(net->SIS.pointerEA[j] == 0){
      net->SIS.EA[0][net->SIS.Eact] = infected;
      net->SIS.EA[1][net->SIS.Eact] = net->LINKS[j];
      net->SIS.pointerEA[j] = net->SIS.Eact+1;
      for(l=net->POINT.INI[index_act];l<=net->POINT.FIN[index_act];l++){
        if(net->LINKS[l] == infected){
          net->SIS.pointerEA[l] = net->SIS.Eact+1;
        }
      }
      net->SIS.Eact += 1;
    }
    else if(net->SIS.pointerEA[j] != 0){
      tmp = net->SIS.pointerEA[j];
      net->SIS.EA[0][tmp-1] = net->SIS.EA[0][net->SIS.Eact-1];
      net->SIS.EA[1][tmp-1] = net->SIS.EA[1][net->SIS.Eact-1];
      net->SIS.pointerEA[j] = 0;
      OverwritePointer(net,infected,net->lenLINKS,tmp,j);
      net->SIS.EA[0][net->SIS.Eact-1] = 0;
      net->SIS.EA[1][net->SIS.Eact-1] = 0;
      net->SIS.Eact -= 1;
    }
  }
}

float FI_function(float rho, float coupling){

  float fi;

  fi = pow(rho,coupling);

  return fi;
}

bool EpidemicDied(struct network *LAYERS, int cl){

  int i;
  bool dead;

  dead = true;
  i=0;

  while((i<cl)&& dead){
    dead = (LAYERS[i].SIS.Nact == 0);
    i += 1;
  }

  return dead;
}
