/*
 *Author: Anay Sonawane
 *File name: linear_mpi.c
 *
 *
 *Pipelined Gaussina elimination using MPI
 *
 *
 *
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "mpi.h"


void print_input_matrix(float *A, int N){

int i,j;

for(i=0;i<N;i++){
	for(j=0; j<N;j++)
	{
		printf("%0.3f ",A[i*N+j]);
	}

	printf("\n");
}

}

void print_sol_matrix(float *A, int N){

int i;

for(i=0;i<N;i++){
	printf("%f ",A[i]);
	printf("\n");
}

}


void fill_matrix(float* A, float *x_sol, int N){

int i,j;
float rowsum;

srand(time(NULL));

for(i=0;i<N;i++){
	rowsum = 0.0;
	for(j=0; j<N;j++)	
	{
		if(i != j){
			
			A[i*N+j]= (float)(rand() % 100);
			rowsum += A[i*N+j];
		}
	}
	
	A[i*N+i] = rowsum + 1;
	
	x_sol[i] = (float)(rand() % 100);
}

}


void mult_matrix(float* A, float* x, float* B, int N)
{
	int i,j;
	float sum = 0.0;

	for(i=0;i<N;i++)
	{
		sum = 0.0;
		for(j=0;j<N;j++)
 		{
			sum +=A[i*N+j]*x[j];
		}
		
		B[i] = sum;
	}
}


void back_substitution(float* A, float* x, float* B, int* loc, int N)
{
	int i,j;
	
	for(i=N-1;i>=0;i--)
	{
		x[loc[i]] = B[i]/A[i*N+loc[i]];
		for(j=0;j<N-1;j++)
		{
			B[j] = B[j] - x[loc[i]]*A[j*N+loc[i]];
		}
	}
	
}

void main(int argc, char *argv[])
{

float *A, *x, *x_sol, *B;
float *A_Block, *B_Block, *local_row, *temp_arr, *x_Block;
int* loc;
int msg[2],picked,temp;
float ratio,temp_ratio;

int N,i,j,k,id,p,mag,msg_recvd;
int i_b,root;

struct timeval start_time, stop_time, elapsed_time;
double total_time,perf;
MPI_Status status;

//FILE *fp=fopen("openmp.csv","a+");

MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&p);
MPI_Comm_rank(MPI_COMM_WORLD,&id);

//print_input_matrix(A,N);
//print_sol_matrix(x_sol,N);
//print_sol_matrix(B,N);

////////////////////////Pipelined Gaussian Elimination////////////////////////////////////	

loc = (int*)malloc(N*sizeof(int));
for(i=0;i<N;i++)
{
	loc[i] = i;
}


//Fill the matrix A and B randomly
if(id == 0)
{
	N = atoi(argv[1]);

	A = (float*)malloc(N*N*sizeof(float));
	B = (float*)malloc(N*sizeof(float));
	x = (float*)malloc(N*sizeof(float));
	x_sol = (float*)malloc(N*sizeof(float));
		
	fill_matrix(A,x_sol,N);
	mult_matrix(A,x_sol,B,N);
	
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
else
{
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
	
//Processing blocks for every node
A_Block = (float*)malloc(N*(N/p)*sizeof(float));
B_Block = (float*)malloc((N/p)*sizeof(float));
x_Block = (float*)malloc((N/p)*sizeof(float));
local_row = (float*)malloc((N+1)*sizeof(float));
temp_arr = (float*)malloc((N*N)*sizeof(float));

//Distribute the data to other nodes 
if(id==0)
{
	k = 0;
	for(i=0;i<N;i++)
	{
		if(i%p == id)
		{
			for(j=0;j<N;j++)
			{
				A_Block[k*N+j] = A[i*N+j];
			}
			
				B_Block[k] = B[i];
				k++;
		}
		else
		{
			MPI_Send(&A[i*N],N,MPI_FLOAT,i%p,0,MPI_COMM_WORLD);
			MPI_Send(&B[i],1,MPI_FLOAT,i%p,1,MPI_COMM_WORLD);
		}
	}	
}
else
{
	for(i=0;i<N/p;i++)
	{
		MPI_Recv(&A_Block[i*N],N,MPI_FLOAT,0,0,MPI_COMM_WORLD,&status);
		MPI_Recv(&B_Block[i],1,MPI_FLOAT,0,1,MPI_COMM_WORLD,&status);

	}
}

MPI_Barrier(MPI_COMM_WORLD);

//printf("%d\n",id);
total_time = -MPI_Wtime();

////////////////////////////////////////Pipelining Starts Here///////////////////////////////////////////////////
i_b = 0;
for(i=0;i<N;i++)
{
	if(i%p == id)
	{       
		mag = 0;
		for(j=i;j<N;j++)
		{
			if(fabsf(A_Block[i_b*N+loc[j]]) > mag)
			{
				mag = fabsf(A_Block[i_b*N+loc[j]]);
				picked = j;
			}
		}
	
		temp = loc[i];
		loc[i] = loc[picked];
		loc[picked] = temp;
		
		memcpy(&local_row[0],&A_Block[i_b*N],N*sizeof(float));
		memcpy(&local_row[N],&B_Block[i_b],sizeof(float));
		
		i_b++;
		//Send the information to the next node
		if(id == p-1){
			msg[0] = i;
			msg[1] = picked;
			MPI_Send(&msg[0],2,MPI_INT,0,10,MPI_COMM_WORLD);
			MPI_Send(&local_row[0],N+1,MPI_FLOAT,0,20,MPI_COMM_WORLD);
		}
		else{
			msg[0] = i;
			msg[1] = picked;
			MPI_Send(&msg[0],2,MPI_INT,id+1,10,MPI_COMM_WORLD);
			MPI_Send(&local_row[0],N+1,MPI_FLOAT,id+1,20,MPI_COMM_WORLD);
		}
	
	}
	else
	{
		//Wait to receve information from previous node
		if(id == 0){
			
			//printf("waiting @ %d\n",id);
			MPI_Recv(&msg[0],2,MPI_INT,p-1,10,MPI_COMM_WORLD, &status);
			MPI_Recv(&local_row[0],N+1,MPI_FLOAT,p-1,20,MPI_COMM_WORLD,&status);
			//printf("Data recvd @ %d\n",id);
		}
		else
		{	// printf("waiting @ %d\n",id);
			 MPI_Recv(&msg[0],2,MPI_INT,id-1,10,MPI_COMM_WORLD,&status);
			 MPI_Recv(&local_row[0],N+1,MPI_FLOAT,id-1,20,MPI_COMM_WORLD,&status);
			// printf("Data recvd @ %d\n",id);
		}

		//Immediately send the information to next node
		if(id == p-1){
			if(msg[0]%p != 0)
			{
				//msg[0] = i;
				//msg[1] = picked;
				MPI_Send(&msg[0],2,MPI_INT,0,10,MPI_COMM_WORLD);
				MPI_Send(&local_row[0],N+1,MPI_FLOAT,0,20,MPI_COMM_WORLD);

			}
		}
		else{
			if((msg[0]%p) != (id+1))
			{
				//msg[0] = i;
				//msg[1] = picked;
				MPI_Send(&msg[0],2,MPI_INT,id+1,10,MPI_COMM_WORLD);
				MPI_Send(&local_row[0],N+1,MPI_FLOAT,id+1,20,MPI_COMM_WORLD);

			}
		}


		picked = msg[1];
		
		if(i_b < N/p)
		{
			temp = loc[i];
			loc[i] = loc[picked];
			loc[picked] = temp;
			
		}
	}	

	

	//Eliminate the elements using the pivot column
	for(j=i+1;j<N;j++)
	{
		if(j%p == id){
				
			//if(fabsf(A[i*N+loc[i]])== 0) ratio = 1.0;
			//else
			ratio = A_Block[(j/p)*N + loc[i]]/local_row[loc[i]];

			//printf("%f\t%f\t%f\n",A[j*N+loc[i]],A[i*N+loc[i]],ratio);
			for(k=i;k<N+1;k++)
			{
				if(k==N)
				{
					B_Block[j/p] = B_Block[j/p] - (local_row[k]*ratio);
				}
				else
				{
					A_Block[(j/p)*N + loc[k]] = A_Block[(j/p)*N + loc[k]] - (local_row[loc[k]]*ratio);
				}
			}	
		}
	}
}


/////////////////////////////////////////////Parallel Back Substitution/////////////////////////////////(Not Working)
/*for(i=N-1,i_b=(N/p)-1;i>=0;i--)
{
	
	root = i%p;
	if(root==id)
	{
		x_Block[i_b] = B_Block[i_b]/A_Block[i_b*N+i];
		temp_ratio = x_Block[i_b];
		i_b--;
		for(k=0;k<p;k++)
		{
			if(k != root)MPI_Send(&temp_ratio,1, MPI_FLOAT,k, 9, MPI_COMM_WORLD);
		}
	}
	else
		MPI_Recv(&temp_ratio,1, MPI_FLOAT,root,9, MPI_COMM_WORLD,&status);

	for(j=0;j<N/p;j++)
	{
		B_Block[j] = B_Block[j] - temp_ratio*A_Block[j*N+i];
	}

}*/


//Get the matrix A and B from other nodes
if(id == 0){
	for(i=0,i_b=0;i<N;i++)
	{
		if(i%p==id)
		{
			for(j=0;j<N;j++)
			{
				A[i*N+j] = A_Block[i_b*N+j];
			}
			B[i] = B_Block[i_b];
			//x[i] = x_Block[i_b];
			i_b++;
		}
		else
		{
			MPI_Recv(&A[i*N],N,MPI_FLOAT,i%p,100,MPI_COMM_WORLD,&status);
			MPI_Recv(&B[i],1,MPI_FLOAT,i%p,200,MPI_COMM_WORLD,&status);
			//MPI_Recv(&x[i],1,MPI_FLOAT,i%p,300,MPI_COMM_WORLD,&status);
		}
	}
}
else
{
	for(i=0;i<N/p;i++)
	{
		MPI_Send(&A_Block[i*N],N,MPI_FLOAT,0,100,MPI_COMM_WORLD);
		MPI_Send(&B_Block[i],1,MPI_FLOAT,0,200,MPI_COMM_WORLD);
		//MPI_Send(&x_Block[i],1,MPI_FLOAT,0,300,MPI_COMM_WORLD);
	}
}

MPI_Barrier(MPI_COMM_WORLD);
total_time += MPI_Wtime();
MPI_Finalize();

if(id ==0) {

	///////////////////////////////////////////Back Substitution////////////////////////////////////////////////
	for(i=N-1;i>=0;i--)
	{
		x[i] = B[i]/A[i*N+i];
		//printf("%f\n",x[i]);
		for(j=0;j<N-1;j++)
		{
			B[j] = B[j] - (x[i] * A[j*N+i]);
		}
	}

	perf = ((2*pow(N,3)/3) + (3*pow(N,2)/2) - (13*N/6))/(total_time*pow(10,9));
	printf ("Processes=%d Total_time=%10.6f Peformance=%lf\n",
			p, total_time, perf);

	//print_input_matrix(A,N);
	//print_sol_matrix(x_sol,N);
	//print_sol_matrix(x,N);
}
}

