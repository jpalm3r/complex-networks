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


int main(){

  char *file, *file1;
  int i;
  char pathFile[4000];
  char *directory = malloc(5000);
  char *folder;
  int* DEGREE_0;
  struct stat dir = {0};

  struct network net1;

  /********************************************************************
  *********************************************************************

  PHASE I: Network file characerization.

    This part of the code is intended to read the source file
    containing the network links in a .txt file with 2 columns,
    each row representing the link and both columns the involved
    nodes. This part returns several vectors containing this
    information.

  *********************************************************************
  *********************************************************************/

  file = malloc(100*sizeof(char));
  file1 = malloc(100*sizeof(char));
  folder = malloc(100*sizeof(char));

  printf("\n");
  printf("\n");
  printf("    ******************************************************\n");
  printf("      N E T W O R K S  6.3:\n");
  printf("\n");
  printf("        - 6.0. Layers, April 2017\n");
  printf("        - 6.1. Structure redefinition, April 2017\n");
  printf("        - 6.2. Structure redefinition, April 2017\n");
  printf("        - 6.3. Separation of edges, April 2017\n");
  printf("        - 6.4. Validation, April 2017\n");
  printf("\n");
  printf("                                      Author: Jaume Palmer\n");
  printf("    ******************************************************\n");
  printf("\n");
  printf("\n");
  printf("      1) Enter the network filename :");
  scanf("%s", file1);
  printf("\n");
  printf("      2) Enter the folder path (./foldername/):");
  scanf("%s", folder);
  printf("\n");

  strcpy(directory,file1);
  if(stat(directory, &dir) == -1)
  {
      mkdir(directory, 0755);
      printf("created directory testdir successfully! \n");
  }

  file = strcat(folder, file1);

  net1.EN = get_E_N(file);

  net1.NODES = malloc(sizeof(int)*net1.EN.N);
  net1.POINT.INI = malloc(net1.EN.N*sizeof(int));
  net1.POINT.FIN = malloc(net1.EN.N*sizeof(int));
  net1.DEGREE = malloc(net1.EN.N*sizeof(int));
  DEGREE_0 = malloc(net1.EN.maxNODE*sizeof(int));

  printf("...getting pointers\n");
  net1.loops = 0;
  get_pointers(&net1, file, DEGREE_0);

  net1.lenLINKS = 2*(net1.EN.E-net1.loops);
  net1.LINKS = (int*)malloc(net1.lenLINKS*sizeof(int));

  InitializeNet(&net1,file, directory);


  /********************************************************************
  *********************************************************************

  PHASE II: SIS model implementation (06/03/17)

  PHASE III: 2 layer multiplex model (24/03/17)

  PHASE IV: Competition (19/04/2017)

  *********************************************************************
  *********************************************************************/

  int j,r,rep;
  int RUNS = 1; int REPS = 1;
  float t, tau, rand1, rand2, t_final;
  float step, delta_step;
  bool keep_running;

  float avg_prv;
  float epsilon = 0.45;
  float coupling = 0.32;
  float SumTotalRate;
  float lambda, MAXlambda; // rate of activation
  float sigma = 1; // rate of recovery
  float mu = 0.1; // rate of spontaneous activation
  float** PREV;
  struct network* LAYERS;
  struct network* LAYERS_0;
  float* G_COEFFS;

  int cl = 2; //Competing networks

  float FRAC[cl];
  float FIT[cl];
  float Dstep[cl];
  float StepNet[cl];

  char** SavingDirs_0 = malloc(cl*sizeof(char*));
  for(i=0;i<cl;i++){SavingDirs_0[i] = malloc(1000*sizeof(char));}

  char** SavingDirs = malloc(cl*sizeof(char*));
  for(i=0;i<cl;i++){SavingDirs[i] = malloc(1000*sizeof(char));}

  FILE** CreatedFiles = malloc(cl*sizeof(FILE*));
  FILE** LayerAverages = malloc(cl*sizeof(FILE*));

  int FileCount[cl];


  for(i=0;i<cl;i++){
    sprintf(SavingDirs_0[i], "./%s/layer_%d",directory,i+1);
    if(stat(SavingDirs_0[i], &dir) == -1){
        mkdir(SavingDirs_0[i], 0755);
    }
  }

  for(i=0;i<cl;i++){
    FileCount[i] = 0;
    sprintf(SavingDirs[i],"./%s/nets",SavingDirs_0[i]);
    if(stat(SavingDirs[i], &dir) == -1){
        mkdir(SavingDirs[i], 0755);
    }
  }

  LAYERS_0 = malloc(cl*sizeof(struct network));
  LAYERS = malloc(cl*sizeof(struct network));

  printf("...starting epidemic spreading\n");

  /***********************************************************************

    Nsus = nodes still only offline (Susceptible)
    net1.SIS.Eact = online nodes susceptible to be activated (Passive)
    net1.SIS.Nact = Active nodes that can become passive (Active)

    Es2a = offline nodes that have neighbours online and can be activated

  ***********************************************************************/

  /***********************************************************************

  In NS:

    '1' = node already in the online (Active or Passive)
    '0' = node susceptible with no neighbours active in the online

  ***********************************************************************/


  printf("...running simulation\n");

  for(i=0;i<cl;i++){
    sprintf(pathFile, "%s/epidemic_output_avg.txt", SavingDirs_0[i]);
    LayerAverages[i] = fopen(pathFile,"w");
  }

  for(i=0;i<cl;i++){

    LAYERS_0[i] = net1;
    LAYERS_0[i].SIS.Nsus = net1.EN.N;

    LAYERS_0[i].SIS.NI = malloc(net1.EN.N*sizeof(int)); // 1D array of infected nodes
    LAYERS_0[i].SIS.pointNS = malloc(sizeof(int)*net1.EN.N);
    LAYERS_0[i].SIS.NS = malloc(sizeof(int)*net1.EN.N);
    LAYERS_0[i].SIS.pointEA_S = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS_0[i].SIS.pointEA_P = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS_0[i].SIS.pointerEA = malloc(net1.lenLINKS*sizeof(int)); //same size as LINKS
    LAYERS_0[i].SIS.EA = malloc(2*sizeof(int*)); // 2D array of active links
    for(j=0;j<2;j++){
      LAYERS_0[i].SIS.EA[j] = malloc(net1.EN.E*sizeof(int));
    }

  }

  //------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------

  for(i=0;i<cl;i++){

    LAYERS[i].SIS.Nsus = net1.EN.N;

    LAYERS[i].NODES = malloc(net1.EN.N*sizeof(int));
    LAYERS[i].LINKS = malloc(net1.lenLINKS*sizeof(int));
    LAYERS[i].DEGREE = malloc(net1.EN.N*sizeof(int));
    LAYERS[i].POINT.INI = malloc(net1.EN.N*sizeof(int));
    LAYERS[i].POINT.FIN = malloc(net1.EN.N*sizeof(int));

    LAYERS[i].SIS.NI = malloc(net1.EN.N*sizeof(int)); // 1D array of infected nodes
    LAYERS[i].SIS.pointNS = malloc(sizeof(int)*net1.EN.N);
    LAYERS[i].SIS.NS = malloc(sizeof(int)*net1.EN.N);
    LAYERS[i].SIS.pointEA_S = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS[i].SIS.pointEA_P = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS[i].SIS.pointerEA = malloc(net1.lenLINKS*sizeof(int)); //same size as LINKS
    LAYERS[i].SIS.EA = malloc(2*sizeof(int*)); // 2D array of active links
    for(j=0;j<2;j++){
      LAYERS[i].SIS.EA[j] = malloc(net1.EN.E*sizeof(int));
    }
  }

  //------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------

  G_COEFFS = malloc((4*cl)*sizeof(float));
  PREV = malloc(cl*sizeof(float*));
  for(i=0;i<cl;i++){
    PREV[i] = malloc(REPS*sizeof(float));
  }


  printf("...>> building initial population                     \n");


  lambda = 0.4; MAXlambda = 1.0;
  // mu = 0.1;
  // sigma = 1;

  // At the beginning there is only one OSN, so initial population is built
  //  with only this



  for(i=0;i<cl;i++){

    printf("\n");
    printf("      %d) Enter the initial fraction of active users in net %d: ",4+i,i+1);
    scanf("%f", &FRAC[i]);
    printf("\n");
  }

  for(i=0;i<cl;i++){FIT[i] = 0.0;}
  for(i=0;i<cl;i++){initial_population(&LAYERS_0[i], FRAC[i], FIT[i]);}


  printf("\r");
  printf("\n");
  printf("\n");
  printf("      %d) Enter the time you want to compute: ",5+(cl-1));
  scanf("%f", &t_final);
  printf("\n");
  if(RUNS>1){
    printf("      -------------------------------------------------------\n");
    printf("      Computing for lambda values between [%f:%f]\n", lambda, MAXlambda);
    printf("      -------------------------------------------------------\n");
    printf("\n");
  }
  else if(RUNS==1){
    REPS = 1;
    printf("      -------------------------------\n");
    printf("      Parameters:\n");
    printf("\n");
    printf("        - lambda = %f\n", lambda);
    printf("        - v = %f\n", lambda/mu);
    printf("        - coupling  = %f\n", coupling);
    printf("\n");
    printf("        - epsilon  = %f\n", epsilon);
    printf("      -------------------------------\n");
    printf("\n");
  }
  printf("...>> spreading infection                             \n");

  clock_t begin = clock();
  srand(time(NULL));

  /***********************************************************************
    BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP
    BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP
    BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP
    BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP
  ***********************************************************************/

  for(i=0;i<cl;i++){
    StepNet[i] = 2; // 2 is the smallest possible net
    Dstep[i] = LAYERS_0[i].EN.N/(float)50;

    FileCount[i] = 0;

    for(j=0;j<REPS;j++){
      PREV[i][j] = 0;
    }
  }

  for(r=1;r<=RUNS;r++){

    for(rep=0;rep<REPS;rep++){

      ResetVectors(LAYERS,LAYERS_0,cl);

      Update_pointNS(LAYERS,cl);
      Update_pointEA_S(LAYERS,cl);
      Compute_Density(LAYERS, cl);
      Update_Weight(LAYERS, cl, epsilon, coupling);
      Update_ActRate(LAYERS, cl, lambda);

      for(i=0;i<cl;i++){
        sprintf(pathFile, "%s/tevol_lambda-%f.txt", SavingDirs_0[i],lambda);
        CreatedFiles[i] = fopen(pathFile,"w");
      }

      step = 0;
      delta_step = t_final/(float)50;

      t = 0;
      keep_running = true;

      // GILLESPIE METHOD
      //    "keep_running" condition is added to be able to
      //    start Gillespie without any active node
      //    and give enough time for the infection to grow

      while(keep_running){

        rand1 = ((double)rand()+1)/ (double)((unsigned)RAND_MAX);
        rand2 = ((double)rand()+1)/ (double)((unsigned)RAND_MAX); //random number between 0 and 1

        //a_react  | a_deact  | a_mmact  | a_act
        Compute_GllspCoeffs(LAYERS, cl, lambda, sigma, mu, G_COEFFS);

        SumTotalRate = 0.0;
        for(i=0;i<cl;i++){
          SumTotalRate += LAYERS[i].SIS.a_i.TOTAL;
        }
        tau = (1/(float)(SumTotalRate))*log(1/rand1);

        for(i=0;i<cl;i++){
          printf("[%d->]  rho_A: %f | rho_P: %f  ",i+1,LAYERS[i].rho.act, LAYERS[i].rho.pas);
        } printf("rep %d/%d. run %d/%d \r",rep+1, REPS, r, RUNS);

        Event_Happens(LAYERS, rand2, G_COEFFS);

        Update_pointNS(LAYERS,cl);
        Update_pointEA_S(LAYERS,cl);
        Compute_Density(LAYERS, cl);
        Update_Weight(LAYERS, cl, epsilon, coupling);
        Update_ActRate(LAYERS, cl, lambda);

        t += tau;

        if(t>=step){
          //only writing the case for one net each run to save time
          for(i=0;i<cl;i++){fprintf(CreatedFiles[i],"%f %f %f\n", t, LAYERS[i].rho.act, LAYERS[i].rho.susc);}

          step += delta_step;
        }


        for(i=0;i<cl;i++){
          if((LAYERS[i].SIS.Nact > StepNet[i])){
            //WriteNet_Act(LAYERS[i],SavingDirs[i],FileCount[i]);
            WriteNet_Total(LAYERS[i],SavingDirs[i],FileCount[i]);
            StepNet[i] += Dstep[i];
            FileCount[i] += 1;
          }
        }

        if(t>1){ keep_running = (t<t_final) && (!EpidemicDied(LAYERS,cl));}
      }

      for(i=0;i<cl;i++){
        fclose(CreatedFiles[i]);

        PREV[i][rep] = LAYERS[i].SIS.Nact/(float)(LAYERS[i].EN.N-LAYERS[i].SIS.Nsus);
      }

    }

    for(i=0;i<cl;i++){
      avg_prv = ComputeMean(PREV[i],REPS);
      fprintf(LayerAverages[i],"%f %f\n",lambda,avg_prv);
    }

    lambda += MAXlambda/(float)RUNS;

  }


  free(G_COEFFS);

  for(i=0;i<cl;i++){
    fclose(LayerAverages[i]);
  }

  if(EpidemicDied(LAYERS,cl)){ printf("Epidemic died.                                           \n");}

  printf("End of SIS run.                                                   \n");
  printf("\n");

  printf("...freeing\n");
  free(DEGREE_0);

  for(i=0;i<cl;i++){
    free(PREV[i]);
  }free(PREV);

  for(i=0;i<cl;i++){

    printf("...freeing 1\r");free(LAYERS[i].SIS.pointerEA);
    printf("...freeing 2\r");free(LAYERS[i].SIS.pointEA_P);
    printf("...freeing 3\r");free(LAYERS[i].SIS.pointEA_S);
    printf("...freeing 4\r");free(LAYERS[i].POINT.INI);
    printf("...freeing 5\r");free(LAYERS[i].POINT.FIN);
    printf("...freeing 6\r");free(LAYERS[i].LINKS);
    printf("...freeing 7\r");free(LAYERS[i].DEGREE);
    printf("...freeing 8\r");free(LAYERS[i].NODES);
    printf("...freeing 9\r");free(LAYERS[i].SIS.NI);
    printf("...freeing 10\r");free(LAYERS[i].SIS.pointNS);
    printf("...freeing 11\r");free(LAYERS[i].SIS.NS);
    printf("...freeing layer %d finished!\n", i+1);
    for(j=0;j<2;j++){
      free(LAYERS[i].SIS.EA[j]);
    }free(LAYERS[i].SIS.EA);
  }
  // Memory leaks !!
  //free(folder);
  //free(file);
  //free(file1);
  //free(directory1);
  //free(directory2);
  //free(directory3);
  //free(directory4);


  clock_t end = clock();
  double time_spent = (double)(end-begin)/CLOCKS_PER_SEC;

  printf("(time spent: %f s)\n", time_spent);
  printf("\n");

  return 0;
}
