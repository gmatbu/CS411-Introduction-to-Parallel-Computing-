#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> 
#include <assert.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

int N, G;
int rank,p;
time_t t;
int root;

// initalize outfiles
FILE * statFile;

// function declaration
void printLocalBlock(int * array, int row, int col);
int * GenerateInitialGoL(int row, int col);
void print_matrix(int **m);
// Need memory access to prev generation array 
// to make changes to it
// addditionally I want to pass by reference
void Simulate(int ** prevGenArr);
void DisplayGoL(int * blockArray);
void ComputeGen(int ** lastGen, int * topRow, int * bottomRow, int ** newGen);

int main(int argc,char *argv[]) {
   struct timeval t1,t2;

   int * array = NULL;
   statFile = fopen("data.csv", "a");
   // time variables
   double totalRuntime;
   root = 0;
   
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_size(MPI_COMM_WORLD,&p);

   srand((unsigned) time(&t));
   printf("my rank=%d\n",rank);
   printf("Rank=%d: number of processes =%d\n",rank,p);

   if (argc >= 2) {
      N = atoi(argv[1]);
      G = atoi(argv[2]);
      printf("Debug: Size of matrix = %d\n", N);
      printf("Debug: Generations = %d\n", G);
   }
    

   // array is editable and has the original copy
   // of allocated memory nby generate initial game of life
   array = GenerateInitialGoL((int)N/p, N);
   printLocalBlock(array, (int)N/p, N);

   Simulate(&array);
   
   free(array);

   MPI_Finalize();
   /*if (rank == root) {
        fprintf(statFile, "\n\n");
   }*/
   fclose(statFile);
   return 0;
}

// initializes only a block of matrix
// row == int N/p
// col == N
// can't generate a dynamically allocated 2D matrix
// as it won't be contiguous in memory 
int * GenerateInitialGoL(int row, int col){
    int i = 0, j = 0;
    // declare N/p rows
    int *m = (int *) malloc(row*col * sizeof(int));
    
    // declare N columns
    /*for(i = 0; i < row; i++){
        m[i] = malloc(col * sizeof(int));
    }*/

    // initialize the block
    for(i = 0; i < row; i++){
        for(j = 0; j < col; j++){
            int seed = rank+1;
            seed = seed*i;
            // even == alive 
            // odd == dead 
            m[i*col + j] = rand_r(&seed) % 10; 
        }
    }

    return m;
}

void printLocalBlock(int * array, int row, int col){

    printf("Rank = %d\n", rank);
    printf("-------------------MATRIX--------------------\n\n");
    // row = int N/p
    // col = N
    for (int i = 0; i < row; i++){
        for (int j = 0; j < col; j++){
                if(j==0){
                    printf("|");
                }
                printf(" %d ", array[i*col + j]);
        }
        printf("|\n");
    }
    printf("\n");
    return;
}

void Simulate(int ** prevGenArr) {
    int * nextGenArr = GenerateInitialGoL((int)N/p, N);
    int * topRow = (int *) malloc(sizeof(int)*N);
    int * bottomRow = (int *) malloc(sizeof(int)*N);
    int prevRank = rank - 1, nextRank = rank + 1;
    MPI_Status statusTop, statusBottom;
    struct timeval gt1, gt2, ct1, ct2, dt1, dt2;
    double generationTime = 0.0, displayTime = 0.0, commTime = 0.0;
    double maxCommTime;
    int displayCounter = 0;

    gettimeofday(&gt1,NULL);
    // simulation
    for(int g=0; g < G; g++, displayCounter++){
        
        gettimeofday(&ct1, NULL);    
        // ensures all processes are excuting
        // the same generation at any given time.
        MPI_Barrier(MPI_COMM_WORLD);
        gettimeofday(&ct2, NULL);
        commTime += (ct2.tv_sec-ct1.tv_sec)*1000000 + (ct2.tv_usec-ct1.tv_usec);
            

        // send and recv
        if (rank == 0) {
            prevRank = p-1;
        }
        if (rank == (p-1)) {
            nextRank = 0;
        }


        gettimeofday(&ct1, NULL);
        // SEND AND RECEIVE ROWS
        // send top row to prevRank and recv bottom row from next rank
        MPI_Sendrecv(*prevGenArr, N, MPI_INT, prevRank, rank, bottomRow, N, MPI_INT, nextRank, MPI_ANY_TAG, MPI_COMM_WORLD, &statusBottom);
        // send bottom row to next rank and recv top row from prev rank
        MPI_Sendrecv(((*prevGenArr)+ (N*(((int)N/p)-1))), N, MPI_INT, prevRank, rank, topRow, N, MPI_INT, nextRank, MPI_ANY_TAG, MPI_COMM_WORLD, &statusBottom);

        gettimeofday(&ct2, NULL);
        commTime += (ct2.tv_sec-ct1.tv_sec)*1000000 + (ct2.tv_usec-ct1.tv_usec);
//        printf("count recvd = %ld\n MPI_SOURCE = %d\n MPI_TAG = %d\n MPI_ERR=%d\n", statusTop._ucount, statusTop.MPI_SOURCE, statusTop.MPI_TAG, statusTop.MPI_ERROR);

      //  printf("print top row for rank %d\n", rank);
      //  printLocalBlock(topRow, 1, N);
      //  printf("print bottom row for rank %d\n", rank);
      //  printLocalBlock(bottomRow, 1, N);

        ComputeGen(prevGenArr, topRow, bottomRow, &nextGenArr);
        if (displayCounter % 3 == 0) {
            gettimeofday(&dt1, NULL);
            // display new generation
            if (rank == 0) {
                printf("---------Displaying Gen %d----------\n\n", g+1);
            }
            DisplayGoL(*prevGenArr);
            gettimeofday(&dt2, NULL);
            displayTime += (dt2.tv_sec-dt1.tv_sec)*1000000 + (dt2.tv_usec-dt1.tv_usec);
        }
    }

    gettimeofday(&gt2,NULL);
    generationTime = (gt2.tv_sec-gt1.tv_sec)*1000000 + (gt2.tv_usec-gt1.tv_usec);


//    MPI_Reduce(&commTime, &maxCommTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0) { 
        printf("generation time before = %lf\n", generationTime);
        printf("comm time before = %lf\n", commTime);
        printf("displaytime before = %lf\n", displayTime);
  //      printf("comm time after = %lf\n", maxCommTime);
    }
    double tRuntime = (generationTime - displayTime);
    double avgRuntime = (generationTime - displayTime)/G;
   if (rank == root) {
        //fprintf(statFile, "G, p, N, Truntime, avg runtime, Tcomm, avg comm, Tcomp, avg comp\n");
        fprintf(statFile, "%d, %d, %d, %lf, %lf, %lf, %lf, %lf, %lf\n", G, p, N, tRuntime, avgRuntime, commTime, commTime/G, tRuntime - commTime, avgRuntime - (commTime/G));
        /*
         fprintf(statFile, "Runtime (Total Runtime - Display Time) = %lf us\n", generationTime - displayTime);
         fprintf(statFile, "Average Runtime = %lf us\n", avgRuntime);  
         fprintf(statFile, "Total Communication Time = %lf us\n", commTime);  
         fprintf(statFile, "Average Communication Time = %lf us\n", commTime/G);  
         fprintf(statFile, "Total Computation Time = %lf us\n", generationTime - displayTime - commTime);  
         fprintf(statFile, "Average Computation Time = %lf us\n", avgRuntime - (commTime/G));  
  */
         }
    // cleanup 
    free(topRow);
    free(bottomRow);
    free(nextGenArr);
    return;
}

void DisplayGoL(int * blockArray){
    // generate nxn matrix to accumulate all the data in
    int * matrix;
    int root = 0;
   
    if(rank == root){
        matrix = (int *)malloc(N*N* sizeof(int));
    }

    MPI_Gather(blockArray, (N* ((int)N/p)), MPI_INT, matrix, (N * ((int)N/p)), MPI_INT, root, MPI_COMM_WORLD);
    
    if (rank == root) {
        printLocalBlock(matrix, N, N);
        free(matrix);
    }
     
    return;
}

void ComputeGen(int ** lastGen, int * topRow, int * bottomRow, int ** newGen) {
    int aliveCount = 0;
    // take care of middle portion of the block
    
    for (int i = 1; i < ((int)N/p)-1; i++){
        for (int j = 1; j < N-1; j++) {
            aliveCount = 8;
            // North
            aliveCount -= (((*lastGen)[(i-1)*N + j]) & 1);
            // South
            aliveCount -= (((*lastGen)[(i+1)*N + j]) & 1);
            // East
            aliveCount -= (((*lastGen)[i*N + j+1]) & 1);
            // West
            aliveCount -= (((*lastGen)[i*N + j-1]) & 1);
            // North East
            aliveCount -= (((*lastGen)[(i-1)*N + j+1]) & 1);
            // North West
            aliveCount -= (((*lastGen)[(i-1)*N + j-1]) & 1);
            // South East
            aliveCount -= (((*lastGen)[(i+1)*N + j+1]) & 1);
            // South West
            aliveCount -= (((*lastGen)[(i+1)*N + j-1]) & 1);
            
            // Fate of cell is decided
            // cell dies of loneliness
            if (aliveCount < 3) {
                (*newGen)[i*N + j] = 1;
            }
            // cell dies of overcrowding
            else if (aliveCount > 5) {
                (*newGen)[i*N + j] = 1;
            }
            // cell is either brought back to life or stays alive
            else {
                (*newGen)[i*N + j] = 0;
            }
        }
    }


    // take care of edges

    // top row and bottom row
    for (int j = 1; j < N-1; j++) { 
        for (int i = 0; i < ((int)N/p); i += ((int)N/p)-1) {
            aliveCount = 8;
            // top row
            if (i -1 < 0){
                // North
                aliveCount -= topRow[j] & 1;
                // South
                aliveCount -= (((*lastGen)[(i+1)*N + j]) & 1);
                // North East
                aliveCount -= topRow[j+1] & 1;
                // North West
                aliveCount -= topRow[j-1] & 1;
                // South East
                aliveCount -= (((*lastGen)[(i+1)*N + j+1]) & 1);
                // South West
                aliveCount -= (((*lastGen)[(i+1)*N + j-1]) & 1);
            }
            else {
                // North
                aliveCount -= (((*lastGen)[(i-1)*N + j]) & 1);       
                // South
                aliveCount -= bottomRow[j] & 1;    
                // North East
                aliveCount -= (((*lastGen)[(i -1)*N + j+1]) & 1);
                // North West
                aliveCount -= (((*lastGen)[(i -1)*N + j-1]) & 1);
                // South East
                aliveCount -= bottomRow[j+1] & 1;
                // South West
                aliveCount -= bottomRow[j-1] & 1;
            }
            // East
            aliveCount -= (((*lastGen)[i*N + j+1]) & 1);
            // West
            aliveCount -= (((*lastGen)[i*N + j-1]) & 1);
        
            // computation
            // Fate of cell is decided
            // cell dies of loneliness
            if (aliveCount < 3) {
                (*newGen)[i*N + j] = 1;
            }
            // cell dies of overcrowding
            else if (aliveCount > 5) {
                (*newGen)[i*N + j] = 1;
            }
            // cell is either brought back to life or stays alive
            else {
                (*newGen)[i*N + j] = 0;
            }
            if (N==p) {
                i += 1;
            }
        } 
    }  



    // left edge and right edge
    for (int i = 1; i < ((int)N/p)-1; i++) {
        for (int j = 0; j < N; j+= N-1) {
           aliveCount = 8;
           // North
           aliveCount -= (((*lastGen)[(i-1)*N + j]) & 1);
           // South
           aliveCount -= (((*lastGen)[(i+1)*N + j]) & 1);
           // East 
           aliveCount -= (((*lastGen)[i*N + ((N + j + 1) % N)]) & 1);
           // West 
           aliveCount -= (((*lastGen)[i*N + ((N + j - 1) % N)]) & 1);
           // North East
           aliveCount -= (((*lastGen)[(i-1)*N + ((N + j + 1) % N)]) & 1);
           // North West
           aliveCount -= (((*lastGen)[(i-1)*N + ((N + j - 1) % N)]) & 1);
           // South East
           aliveCount -= (((*lastGen)[(i+1)*N + ((N + j + 1) % N)]) & 1);
           // South West
           aliveCount -= (((*lastGen)[(i+1)*N + ((N + j - 1) % N)]) & 1);


           // computation
           // Fate of cell is decided
           // cell dies of loneliness
           if (aliveCount < 3) {
               (*newGen)[i*N + j] = 1;
           }
           // cell dies of overcrowding
           else if (aliveCount > 5) {
               (*newGen)[i*N + j] = 1;
           }
           // cell is either brought back to life or stays alive
           else {
               (*newGen)[i*N + j] = 0;
           }
        }
    }

    // take care of corners
    for (int i = 0; i < (int)N/p; i += ((int)N/p)-1) {
        for (int j = 0; j < N; j += N-1) {
            aliveCount = 8;
            // top row
            if (i == 0) {
                // North
                aliveCount -= (topRow[j] & 1);
                // South
                aliveCount -= (((*lastGen)[(i+1)*N + j]) & 1);
                // North East
                aliveCount -= (topRow[(N + j + 1) % N] & 1);
                // North West
                aliveCount -= (topRow[(N + j - 1) % N] & 1);
                // South East
                aliveCount -= (((*lastGen)[(i+1)*N + ((N + j + 1) % N)]) & 1);
                // South West
                aliveCount -= (((*lastGen)[(i+1)*N + ((N + j - 1) % N)]) & 1);
            }
            // bottom row
            else {
                // North
                aliveCount -= (((*lastGen)[(i-1)*N + j]) & 1);
                // South
                aliveCount -= (bottomRow[j] & 1);
                // North East
                aliveCount -= (((*lastGen)[(i-1)*N + ((N + j + 1) % N)]) & 1);
                // North West
                aliveCount -= (((*lastGen)[(i+1)*N + ((N + j - 1) % N)]) & 1);
                // South East
                aliveCount -= (bottomRow[(N + j + 1) % N] & 1);
                // South West
                aliveCount -= (bottomRow[(N + j - 1) % N] & 1);
            }
           // East
           aliveCount -= (((*lastGen)[i*N + ((N + j + 1) % N)]) & 1);
           // West
           aliveCount -= (((*lastGen)[i*N + ((N + j - 1) % N)]) & 1);
        
           // computation
           // Fate of cell is decided
           // cell dies of loneliness
           if (aliveCount < 3) {
               (*newGen)[i*N + j] = 1;
           }
           // cell dies of overcrowding
           else if (aliveCount > 5) {
               (*newGen)[i*N + j] = 1;
           }
           // cell is either brought back to life or stays alive
           else {
               (*newGen)[i*N + j] = 0;
           }
        }
        if (N == p){
             i += 1;
        }
    } 
    // copy newGen into lastGen
    for (int i = 0; i < N * ((int)N/p); i++) {
        (*lastGen)[i] = (*newGen)[i];
    }

    return;
}
