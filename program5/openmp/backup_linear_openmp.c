#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>


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

void gaussian_elimination(float*A, float*B, int N, int p)
{

int i, j, k, id, *column_loc, *msg, picked=0,temp;
float ratio,mag;

int init=0;

msg = (int*)calloc(p*N,sizeof(int));
column_loc = (int*)malloc(p*sizeof(int));

#pragma omp parallel default (none) private(i,j,k,id,mag,ratio,temp,picked) shared(A,B,N,p,msg,column_loc) num_threads(p)
{
	
	id = omp_get_thread_num();

	for(i=0;i<N;i++)
	{
		if(i%p == id)
		{
			mag = 0;
			for(j=i;j<N;j++)
			{
				if(fabsf(A[i*N+j]) > mag)
				{
					mag = fabsf(A[i*N+j]);
					picked = j;
				}
			}
			
		}
		else
		{	
			//printf("waiting @ %d on %d\n",id,i);
			while(!msg[id*N+i]);
			//printf("entered @ %d on %d\n",id,i);
			msg[id] = 0;
			picked = column_loc[id];
		}

		if(id == p-1){
			if((i%p)!=0)
			{
				msg[i] = 1;
				column_loc[0] = picked;
				//printf("signaled 0 on %d\n",i);
			}
		}
		else{
			if((i%p) != (id+1))
			{
				msg[(id+1)*N+i] = 1;
				column_loc[id+1] = picked;
				//printf("signaled %d on %d\n",id+1,i);
			}
		}
	
		for(j=i+1;j<N;j++)
		{
			if(j%p == id){
				
				//if(fabsf(A[i*N+loc[i]])== 0) ratio = 1.0;
				//else
				ratio = A[j*N + i]/A[i*N + i];

				//printf("%f\t%f\t%f\n",A[j*N+loc[i]],A[i*N+loc[i]],ratio);
				for(k=i;k<N+1;k++)
				{
					if(k==N)
					{
						B[j] = B[j] - (B[i]*ratio);
					}
					else
					{
						A[j*N + k] = A[j*N + k] - (A[i*N + k]*ratio);
					}
				}	
			}
		}
	}
	
}
	
}

void back_substitution(float* A, float* x, float* B, int N)
{
	int i,j;
	
	for(i=N-1;i>=0;i--)
	{
		x[i] = B[i]/A[i*N+i];
		for(j=0;j<N-1;j++)
		{
			B[j] = B[j] - x[i]*A[j*N+i];
		}
	}
	
}

void main(int argc, char *argv[])
{

float *A, *x, *x_sol, *B;
int* loc_global;

int N,i,threads,p;
struct timeval start_time, stop_time, elapsed_time;
double total_time,perf;

//FILE *fp=fopen("openmp.csv","a+");

N = atoi(argv[1]);
p = omp_get_max_threads();

A = (float*)malloc(N*N*sizeof(float));
B = (float*)malloc(N*sizeof(float));
x = (float*)malloc(N*sizeof(float));
x_sol = (float*)malloc(N*sizeof(float));

fill_matrix(A,x_sol,N);
mult_matrix(A,x_sol,B,N);

//print_input_matrix(A,N);
print_sol_matrix(x_sol,N);
//print_sol_matrix(B,N);

gaussian_elimination(A,B,N,p);

back_substitution(A,x,B,N);

//print_input_matrix(A,N);
print_sol_matrix(x,N);

}

