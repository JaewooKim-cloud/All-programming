#include <iostream>
#include <array>
#include <cmath>
#include <mpi.h>
#include <ctime>

void genericScan(float *a, const long long n)
{	
  MPI_Status status;
	
  int nprocs,procno,d;

  MPI_Comm_size(MPI_COMM_WORLD,&nprocs); // nprocs : the number of threads
  MPI_Comm_rank(MPI_COMM_WORLD,&procno); // procno : The current id of the thread 
  int p = nprocs; // The number of threads/processes.
  d = int(log2(n));
  float **b, **c;
  b = new float * [d+1];
  c = new float * [d+1];
  
  for ( int i = 0; i < d+1; i++ ) {
	  	b[i] = new float[n];
		c[i] = new float[n];
	}	
  /*
   My code starts here. 
   */
   // up-sweep phase
  for (long long i=(n/p)*procno; i<(n/p)*(procno+1);i++) {
	  b[0][i]=a[i];
  }
		if(procno==0) {		 
  			MPI_Gather ( MPI_IN_PLACE , n/(p) , MPI_FLOAT , b[0] , n/(p) , MPI_FLOAT , 0,MPI_COMM_WORLD );
		}
		else {
			MPI_Gather ( &b[0][procno*n/(p)] , n/(p) , MPI_FLOAT , b[0] , n/(p) , MPI_FLOAT , 0 ,MPI_COMM_WORLD );
		}
		MPI_Bcast(&b[0][0], n, MPI_FLOAT, 0, MPI_COMM_WORLD);
	for (int h=1;h<=d;h++) {
	  
	  if (n<long(pow(2,h))*p) {
		for (long long i=0;i<(n/int(pow(2,h)));i++) {
	  
			b[h][i] = b[h-1][2*i] + b[h-1][2*i+1];
			 
				}
				
			
		   
		  
		} else {
	for (long long i=(n/int(pow(2,h))/p)*procno;i<(n/int(pow(2,h))/p)*(procno+1);i++) {
	  
	  b[h][i] = b[h-1][2*i] + b[h-1][2*i+1]; 
	  
	}    
			if(procno==0) {
	  		MPI_Gather ( MPI_IN_PLACE , n/int(pow(2,h))/p , MPI_FLOAT , b[h] , n/(pow(2,h))/p , MPI_FLOAT , 0 ,MPI_COMM_WORLD );
			}
			else {
			MPI_Gather ( &b[h][(n/int(pow(2,h))/p)*procno] , n/int(pow(2,h))/p , MPI_FLOAT , b[h] , n/(pow(2,h))/p , MPI_FLOAT , 0 ,MPI_COMM_WORLD );
			}
			MPI_Bcast(&b[h][0], n/int(pow(2,h)), MPI_FLOAT, 0, MPI_COMM_WORLD);
			
 		}
  }
 
  // down-sweep phase
  c[d][0]=0;
 for(int h = d-1; h >= 0 ; h--) {
	 
	 if (n<int(pow(2,h))*p) {
	 for (long long i=0;i<(n/int(pow(2,h)));i++) {
		 if (i%2==0) {
			 c[h][i]=c[h+1][i/2];
			 
		 } else {
			 c[h][i]=c[h+1][(i-1)/2] + b[h][i-1];
			 
		 }
	 }
		 
	 } else {
	 for (long long i=(n/int(pow(2,h))/p)*procno;i<(n/int(pow(2,h))/p)*(procno+1);i++) {
		 if (i%2==0) {
			 c[h][i]=c[h+1][i/2];
			 
		 } else {
			 c[h][i]=c[h+1][(i-1)/2] + b[h][i-1];
			 
		 }
	 }
		if(procno==0){
	 	MPI_Gather ( MPI_IN_PLACE , (n/int(pow(2,h))/p) , MPI_FLOAT , c[h] , (n/(pow(2,h))/p) , MPI_FLOAT , 0 ,MPI_COMM_WORLD );
		}
		else {
		MPI_Gather ( & c[h][(n/int(pow(2,h))/p)*procno] , (n/int(pow(2,h))/p) , MPI_FLOAT , c[h] , (n/(pow(2,h))/p) , MPI_FLOAT , 0 ,MPI_COMM_WORLD );
		}
		MPI_Bcast(&c[h][0], n/int(pow(2,h)), MPI_FLOAT, 0, MPI_COMM_WORLD);	 
		
	}
 }
 
 for (long long i=(n/p)*procno; i<(n/p)*(procno+1);i++)  {
	 a[i] = a[i]+c[0][i];
 }
		if(procno==0) {
	 	MPI_Gather ( MPI_IN_PLACE , (n/p) , MPI_FLOAT , a , (n/p) , MPI_FLOAT , 0 ,MPI_COMM_WORLD );
		}
		else{
		MPI_Gather ( & a[(n/p)*procno] , (n/p) , MPI_FLOAT , a , (n/p) , MPI_FLOAT , 0 ,MPI_COMM_WORLD );
		}
	
 // The end of the parallel scan code.
}

int main( int argc, char **argv){	

  MPI_Init(0,0);
  MPI_Comm comm = MPI_COMM_WORLD;
  int nprocs,procno,d;
    std::clock_t start;
    double duration;

  MPI_Comm_size(MPI_COMM_WORLD,&nprocs); // nprocs : the number of threads
  MPI_Comm_rank(MPI_COMM_WORLD,&procno); // procno : The current id of the thread 
  int p = nprocs; // The number of threads/processes.
	
   
	// problem setup
	long long n=1024*1024*1024;
	//n=n*2;
	std::cout << n << "\n";
	if(argc>1) n=atoi(argv[1]);
	float *x = (float *) malloc( sizeof(float)*n) ;

	for(long long i=(n/p)*procno;i<(n/p)*(procno+1);i++) x[i] = i%2;
	
	if(procno==0){	
	
	MPI_Gather ( MPI_IN_PLACE , (n/p) , MPI_FLOAT , x , (n/p) , MPI_FLOAT , 0 ,MPI_COMM_WORLD );	
	}
	else
	{

	MPI_Gather ( &x[(n/p)*procno] , (n/p) , MPI_FLOAT , x , (n/p) , MPI_FLOAT , 0 ,MPI_COMM_WORLD );


	}

	
	// print input
	std::cout << "The initialization has been complete." <<"\n";	 

	// scan
	MPI_Barrier( MPI_COMM_WORLD );
start = std::clock(); // get current time
	genericScan(x,  n);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	MPI_Barrier( MPI_COMM_WORLD );
    std::cout << "Operation took "<< duration << "seconds" << std::endl;

	//print output
	
	if (procno==0){
	for(long long i=0;i<100;i++) std::cout<<x[i]<<" "; 	printf("\n");
	}
	MPI_Finalize(); // The end of parallel computation.
	// clean up
	free(x);

	return 0;
}
	



	


	



