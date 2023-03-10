/*	Cpt S 411, Introduction to Parallel Computing
 *	School of Electrical Engineering and Computer Science
 *	
 *	Example code
 *	Send receive test:
 *   	rank 1 sends to rank 0 (all other ranks sit idle)
 *   	For timing use of C gettimeofday() is recommended.
 * */


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> 
#include <assert.h>
#include <sys/time.h>

int main(int argc,char *argv[])
{

   int rank,p;
   struct timeval t1,t2;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_size(MPI_COMM_WORLD,&p);

   printf("my rank=%d\n",rank);
   printf("Rank=%d: number of processes =%d\n",rank,p);

   // make sure there are at least two processes - send and receive
   assert(p>=2);
  
   int maxSize = 1;
   // number of iterations to record average value
   int loops = 10;

   // get maxSize (provided in MB)
     if (argc >= 2) {
        maxSize = atoi(argv[1]);
        assert(maxSize >= 1);
        printf("Debug: input size is = %d\n", maxSize);
    }

   /*
    * check argument variables
    *
   printf("argc = %d\n", argc);
     printf("argv = ");
   for (int i = 0; i < argc; i++){
           printf("%s ", argv[i]);
   }*/
 

    //initialize buffers
    char *x = (char *) calloc (1024*1025*maxSize, sizeof(char));
    char *y = (char *) calloc (1024*1025*maxSize, sizeof(char));

    // time variable declaration
    double tSend, tRecv;
    // initalize outfiles
    FILE * SendOutfile = fopen("nbSend.txt", "w");
    FILE * RecvOutfile = fopen("nbRecv.txt", "w");
    
    int destRank = 0;
   
    MPI_Request request;
    MPI_Status  status;
    // loop sends data from 1Byte -> maxSize MB (multiply by 1024*1024)
    // send size is a number but unit of data is sent in byte
   for(int sendSize = 1; sendSize <= maxSize*1024*1024; sendSize *= 2){
 
           // process 1 is the sender
            if(rank==1) {
                //x = (char *) calloc (sendSize, sizeof(char));
                gettimeofday(&t1,NULL);
                for (int i = 0; i < loops; i++) {
                    MPI_Send(x,sendSize,MPI_BYTE,destRank,0,MPI_COMM_WORLD);
                }
                gettimeofday(&t2,NULL);
                tSend = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
                fprintf(SendOutfile,"Rank=%d: sent message size %d B; Send time %.2f microsec\n", rank,sendSize,(double)tSend/loops);
                //free(x);
                //x = NULL;
            }
            else if (rank==0) {
                //y = (char *) calloc (sendSize,sizeof(char));
                gettimeofday(&t1,NULL);
                for(int i = 0; i < loops; i++) {
                    MPI_Irecv(y,sendSize,MPI_BYTE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD, &request);
                    MPI_Wait(&request, &status);
                }
                gettimeofday(&t2,NULL);
                tRecv = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);      
                fprintf(RecvOutfile,"Rank=%d: received message size %d B; Recv time %.2f microsec\n",rank,sendSize,(double)tRecv/loops);
                //free(y);
                //y = NULL;
            }

  }

    fclose(SendOutfile);
    fclose(RecvOutfile);
   MPI_Finalize();
}
