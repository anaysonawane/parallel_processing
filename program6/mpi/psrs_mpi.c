/*Author : Anay Sonawanr
 *File Name: psrs_openmp.c
 *
 *psrs implementation using OpenMP
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"

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

/*void psrs(int *A, int N, int p)
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

	/*i=0;
	for(k=0;k<p;k++)
	{
		for(j=part_offset[k*p+id];j<part_size[k*p+id]+part_offset[k*p+id];j++)
		{
			a_final[i] = A[j+k*(N/p)];
			i++;		
		}
		
	}
	
	quicksort(a_final,0,i-1);
	//print_array(a_final,i);

}
	
}*/


void main(int argc, char *argv[])
{

int *A,*local_A,*sample_A, *final_A, *part_size, *part_size_rcv, *part_offset, *part_offset_rcv,*final_size, *final_offset, *cnt, *disp;
int N,p,id,i,j,k,offset,size;

double elapsed_time, total_time;

MPI_Status status;

FILE *fp;

MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&p);
MPI_Comm_rank(MPI_COMM_WORLD,&id);

if(id == 0)
{
	N = atoi(argv[1]);
	A = (int*)malloc(N*sizeof(int));
	fill_array(A,N);
	printf("Input List:\n");
	print_array(A,N);
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
else
{
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

local_A = (int*)malloc(N/p*sizeof(int));
final_A = (int*)malloc(N*sizeof(int));
sample_A = (int*)malloc(p*sizeof(int));
part_size = (int*)malloc(p*sizeof(int));
part_offset = (int*)malloc(p*sizeof(int));
part_size_rcv = (int*)malloc(p*sizeof(int));
part_offset_rcv = (int*)malloc(p*sizeof(int));
final_size = (int*)malloc(p*sizeof(int));
final_offset = (int*)malloc(p*sizeof(int));

cnt = (int*)malloc(p*sizeof(int));
disp = (int*)malloc(p*sizeof(int));

MPI_Scatter(A,N/p, MPI_INT,local_A,N/p, MPI_INT, 0, MPI_COMM_WORLD);

MPI_Barrier(MPI_COMM_WORLD);

total_time = -MPI_Wtime();
//quick sort at every node
quicksort(local_A,0,N/p-1);
/*printf("Sorted Array from %d :",id);
print_array(local_A,N/p);
*/

if(id == p-1)
{

	int* sample_all = (int*)malloc(p*p*sizeof(int));
	
	k = 0;
	for(i=0;i<p-1;i++)
	{
		MPI_Recv(sample_A,p,MPI_INT,i,0,MPI_COMM_WORLD,&status);
		for(j=0;j<p;j++,k++)
		{
			//printf("%d\n",j);
			sample_all[k] = sample_A[j];
		}
	}

        for(i=k,j=0;i<p+k;i++,j=j+p)
	{
		sample_all[i] = local_A[j]; 
	}

	quicksort(sample_all,0,(p*p)-1);

	/*printf("Sorted samples :");
	print_array(sample_all,p*p);
	*/
	for(i=0,j=p;i<p-1;i++,j=j+p)
	{
		sample_A[i] = sample_all[j];
	}
	
	MPI_Bcast(sample_A, p-1, MPI_INT, p-1, MPI_COMM_WORLD);
	
	/*printf("Pivots: ");
	print_array(sample_A,p-1);
	*/	
	free(sample_all);
}
else
{
	//Smple sorted array
	for(i=0,j=0;i<p;i++,j=j+p)
	{
		sample_A[i] = local_A[j]; 
	}

	MPI_Send(sample_A,p,MPI_INT,p-1,0,MPI_COMM_WORLD);

	MPI_Bcast(sample_A, p-1, MPI_INT, p-1, MPI_COMM_WORLD);	
	
}

	i=0;
	offset = 0;
	for(k=0;k<p;k++)
	{
		if(k < p-1)
		{
			while((local_A[i] <= sample_A[k]) && (i<N/p))
			{
				i++;
			}	

			part_size[k] = i-offset;
			part_offset[k] = offset;
		
			offset = i;
				
		}
		else if(k == p-1)
		{
			part_size[k] = N/p-offset;
			part_offset[k] = offset;
		}
	}


	/*printf("Offset Array from %d :",id);
	print_array(part_offset,p);
	printf("Size Array from %d :",id);
	print_array(part_size,p);
	*/

	for(i=0;i<p;i++)
	{
		cnt[i] = 1;
		disp[i] = i;
	}

	MPI_Alltoallv(part_size,cnt,disp,MPI_INT,
			part_size_rcv,cnt,disp,MPI_INT,MPI_COMM_WORLD);

	free(cnt);
	free(disp);

	offset = 0;
	for(i=0;i<p;i++)
	{
		part_offset_rcv[i] = offset;
		offset = offset+part_size_rcv[i];
	}
	
	/*printf("Rcv Offset Array from %d :",id);
	print_array(part_offset_rcv,p);
	printf("Rcv Size Array from %d :",id);
	print_array(part_size_rcv,p);
	*/

	MPI_Alltoallv(local_A,part_size,part_offset,MPI_INT,
			final_A,part_size_rcv,part_offset_rcv,MPI_INT,MPI_COMM_WORLD);

	free(part_size);
	free(part_offset);
	
	quicksort(final_A,0,part_size_rcv[p-1]+part_offset_rcv[p-1]-1);
	
	//if(id==0)
//	{
	printf("Final Array From %d :",id);
	print_array(final_A, part_size_rcv[p-1]+part_offset_rcv[p-1]);
//	}
	
	if(id == 0)
	{
		i=0;
		final_size[i]=part_size_rcv[p-1]+part_offset_rcv[p-1];
		
		for(i = 1;i<p;i++)
		{
			MPI_Recv(&size,1,MPI_INT,i,0,MPI_COMM_WORLD,&status);
			final_size[i]=size;
		}

		offset=0;
		for(i=0;i<p;i++)
		{
			final_offset[i]=offset;
			offset = offset + final_size[i];
		}
		
		size = final_size[0];
		
	}
	else
	{
		size = part_size_rcv[p-1]+part_offset_rcv[p-1];
		MPI_Send(&size,1,MPI_INT,0,0,MPI_COMM_WORLD);
		
	}

	free(part_size_rcv);
	free(part_offset_rcv);
	
	MPI_Gatherv(final_A,size,MPI_INT,A,final_size,final_offset,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	
	total_time += MPI_Wtime();
	if(id == 0){
		fp = fopen("mpi.csv","a+");
                printf("Final Sorted Array :\n");
                print_array(A,N);
                printf ("Processes=%d Total_time=%10.6f\n", p, total_time);
                fprintf(fp,"%d,%d,%lf\n", p,N,elapsed_time);
                fclose(fp);
	}
	MPI_Finalize();

}


