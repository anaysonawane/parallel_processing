/*Author : Anay Sonawanr
 *File Name: psrs_openmp.c
 *
 *psrs implementation using OpenMP
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

void print_array(int *A, int N)
{
	int i;
	for(i=0;i<N;i++)
	{
		printf("%d ",A[i]);
	}
	printf("\n");
}

void fill_array(int *A, int N )
{
 int i;

 srand(time(NULL));
 for(i=0;i<N;i++)
 {
	A[i] = rand()%1000;
 }

}

int partition(int *A, int left, int right, int pivot)
{
	int temp;
	int leftP = left-1;
	int rightP = right;

	while(1)
	{
		while(A[++leftP]<pivot){}
		while(rightP>0 && A[--rightP]>pivot){}
		
		if(leftP >= rightP){
			break;
		}
		else
		{
			temp = A[leftP];
			A[leftP] = A[rightP];
			A[rightP] = temp;
		}
	}
	
	temp = A[leftP];
	A[leftP] = A[right];
	A[right] = temp;

	return leftP;
}

void quicksort(int *A, int left, int right)
{
	
	if(right-left <= 0)
	{
		return;
	}
	else
	{
		int pivot = A[right];
		int part = partition(A,left,right,pivot);
		quicksort(A,left,part-1);
		quicksort(A,part+1,right);
	}
	
}

void psrs(int *A, int N, int p)
{
int i,j,k;
int *a;
int left, right,id,offset;
int *pivot_a, *part_size, *part_offset, *a_final;

pivot_a = (int*)malloc(p*p*sizeof(int));
part_offset = (int*)malloc(p*p*sizeof(int));
part_size = (int*)malloc(p*p*sizeof(int));


#pragma omp parallel default(none) private(left,right,id,a_final,i,j,k,offset) shared(A,N,p,a,pivot_a,part_size,part_offset) num_threads(p)
{

	id = omp_get_thread_num();

	a_final = (int*)malloc(N*sizeof(int));

	left = id*N/p;
	right = left+(N/p)-1;	

	quicksort(A,left,right);
	
	for(i=0,j=left;i<p;i++,j=j+p) //get the pivot values from all the process
	{
		pivot_a[p*id+i] = A[j]; 
	}
	
	if(id == p-1){

		quicksort(pivot_a,0,(p*p)-1);
		//print_array(pivot_a,p*p);	
	}	
	
	#pragma omp barrier

	i=left;
	j=p;
	offset = 0;
	for(k=0;k<p;k++)
	{
		if(k < p-1)
		{
			while((A[i] <= pivot_a[j]) && (i<left+N/p))
			{
				i++;
			}	

			part_size[id*p+k] = i-left-offset;
			part_offset[id*p+k] = offset;
		
			offset = i-left;
			j=j+p;
				
		}
		else if(k == p-1)
		{
			part_size[id*p+k] = N/p-offset;
			part_offset[id*p+k] = offset;
		}
	}
		
	#pragma omp barrier

	/*if(id == 0)
	{
		print_array(part_offset,p*p);
		print_array(part_size,p*p);
	}*/

	i=0;
	for(k=0;k<p;k++)
	{
		for(j=part_offset[k*p+id];j<part_size[k*p+id]+part_offset[k*p+id];j++)
		{
			a_final[i] = A[j+k*(N/p)];
			i++;		
		}
		
	}
	
	part_size[id]=i;
	
	quicksort(a_final,0,i-1);
	printf("Sublist from process %d: \n",id);
	print_array(a_final,i);
	printf("==========================================================\n");

	#pragma omp barrier

	if(id == 0)
	{
		offset=0;
		for(i=0;i<p;i++)
		{
			part_offset[i]=offset;
			offset=offset+part_size[i];
		}

		//print_array(part_offset,p);
		//print_array(part_size,p);
	}

	#pragma omp barrier
	
	for(i=part_offset[id],j=0;i<part_offset[id]+part_size[id];j++,i++)
	{
		A[i] = a_final[j];
	}

}
	
}


void main(int argc, char *argv[])
{

int *A;
int N,p;

struct timeval start_time, stop_time, elapsed_time;
double total_time,perf;

FILE *fp=fopen("openmp.csv","a+");

N = atoi(argv[1]);
p = omp_get_max_threads();

A = (int*)malloc(N*sizeof(int));

fill_array(A,N);
printf("Input Array: \n");
print_array(A,N);
printf("========================================================== \n");

//To calculaet the time required
gettimeofday(&start_time, NULL);

psrs(A,N,p);

gettimeofday(&stop_time, NULL);

timersub(&stop_time, &start_time, &elapsed_time);

//print solution matrix
printf("Actual Solution :\n");
print_array(A,N);
printf("========================================================== \n");

//Performance Calculations
total_time =(double)( elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
printf("Threads=%d Total_time=%lf \n",p,total_time);
fprintf(fp,"%d,%d,%lf\n",p,N,total_time);
fclose(fp);

}


