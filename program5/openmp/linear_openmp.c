/*Author : Anay Sonawanr
 *File Name: linear_openmp.c
 *
 *pipelined gaussian elimination
 */
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


//Pipelined Gaussian Elimination Function
void gaussian_elimination(float*A, float*B, int*loc_global, int N, int p)
{

int i, j, k, id, i_b, *loc, *msg, *column_loc, picked=0,temp;
float ratio,mag;

int init=0;

msg = (int*)calloc(N*p,sizeof(int));
column_loc = (int*)malloc(N*p*sizeof(int));

//Start Parallel threads
#pragma omp parallel default (none) private(i,j,k,id,loc,i_b,mag,ratio,temp,picked) shared(A,B,loc_global,N,p,msg,column_loc) num_threads(p)
{
	
	loc = (int*)malloc(N*sizeof(int));

	for(i=0;i<N;i++)
	{
		loc[i] = i;	
	}
	
	id = omp_get_thread_num();

	//Start Guassian Elimination Algorithm
	for(i=0,i_b=0;i<N;i++)
	{
		//Finding the pivot column
		if(i%p == id)
		{
			mag = 0;
			for(j=i;j<N;j++)
			{
				if(fabsf(A[i*N+loc[j]]) > mag)
				{
					mag = fabsf(A[i*N+loc[j]]);
					picked = j;
				}
			}

			i_b++;

			temp = loc[i];
			loc[i] = loc[picked];
			loc[picked] = temp;

			//Send info to the next thread
			if(id == p-1){
				msg[0+i] = 1;
				column_loc[0+i] = picked;
			}
			else{
				msg[(id+1)*N+i] = 1;
				column_loc[(id+1)*N+i] = picked;
			}

			
		}
		else
		{	

			//Recev info from the previous thread
			while(!msg[id*N+i]);
			msg[id*N+i] = 0;
			picked = column_loc[id*N+i];

			//Immediately send the info to next node
			if(id == p-1){
				if((i%p)!=0)
				{
					msg[0+i] = 1;
					column_loc[0+i] = picked;
				}
			}
			else{
				if((i%p) != (id+1))
				{
					msg[(id+1)*N+i] = 1;
					column_loc[(id+1)*N+i] = picked;
				}
			}


			if(i_b < N/p)
			{
				temp = loc[i];
				loc[i] = loc[picked];
				loc[picked] = temp;
			}

		}
		
		//eliminate the elements using pivot columns
		for(j=i+1;j<N;j++)
		{
			if(j%p == id){
				
				//if(fabsf(A[i*N+loc[i]])== 0) ratio = 1.0;
				//else
				ratio = A[j*N + loc[i]]/A[i*N + loc[i]];

				//printf("%f\t%f\t%f\n",A[j*N+loc[i]],A[i*N+loc[i]],ratio);
				for(k=i;k<N+1;k++)
				{
					if(k==N)
					{
						B[j] = B[j] - (B[i]*ratio);
					}
					else
					{
						A[j*N + loc[k]] = A[j*N + loc[k]] - (A[i*N + loc[k]]*ratio);
					}
				}	
			}
		}
	}
	
	//Partial Pivot column
	if(id == 0){	
		for(i=0;i<N;i++)
		{	
			loc_global[i] = loc[i];
		}
	}
}
	
}


//Back substitution to find the solution of matrix
void back_substitution(float* A, float* x, float* B, int* loc, int N)
{
	int i,j;
	
	for(i=N-1;i>=0;i--)
	{
		x[loc[i]] = B[i]/A[i*N+loc[i]];
	
		#pragma omp parallel for private(j)	
		for(j=0;j<N-1;j++)
		{
			B[j] = B[j] - x[loc[i]]*A[j*N+loc[i]];
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
loc_global = (int*)malloc(N*sizeof(int));

fill_matrix(A,x_sol,N);
mult_matrix(A,x_sol,B,N);

printf("Expexted Solution : \n");
print_sol_matrix(x_sol,N);

//To calculaet the time required
gettimeofday(&start_time, NULL);

//Gaussian Elimination
gaussian_elimination(A,B,loc_global,N,p);
//Back Substitution
back_substitution(A,x,B,loc_global,N);

gettimeofday(&stop_time, NULL);

timersub(&stop_time, &start_time, &elapsed_time);

//print solution matrix
printf("Actual Solution :\n");
print_sol_matrix(x,N);

//Performance Calculations
total_time =(double)( elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

perf = ((2*pow(N,3)/3) + (3*pow(N,2)/2) - (13*N/6))/(total_time*pow(10,9));
printf("Threads=%d Total_time=%lf Performance=%lf \n",p,total_time,perf);

}


