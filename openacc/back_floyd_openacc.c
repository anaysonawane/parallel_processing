#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define MIN(a,b) ((a)<(b)?(a):(b))

void compute_shortest_paths(int[][6],int);

int main(int argc, char *argv[])
{

int n=6, i, j;
int arr[6][6];
struct timeval start_time, stop_time, elapsed_time;

//printf("Enter the dimensions of matrix: \n");
//scanf("%d",&n);

//arr = (int*)malloc(n*n*sizeof(int));

for (i=0;i<n;i++){

	for(j=0;j<n;j++){
		if(i==j){
			arr[i][j]=0;
		}
		else if((i==j-1) || (j==i -1)){
			arr[i][j] = 1;
		}
		else{
			arr[i][j]=n;
		}
	}
}


//Print the array
/*printf("input:\n");
for(i=0;i<n;i++)
{
	for(j=0;j<n;j++)
	{
		printf("%d ", arr[i*n+j]);
	}
	printf("\n");
}*/


gettimeofday(&start_time, NULL);

compute_shortest_paths(arr,n);

gettimeofday(&stop_time, NULL);

timersub(&stop_time, &start_time, &elapsed_time);
printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

printf("output: \n");
for(i=0;i<n;i++)
{
	for(j=0;j<n;j++)
	{
		printf("%d ", arr[i][j]);
	}
	printf("\n");
}



}


void compute_shortest_paths (int a[][6], int n)
{
   int	i, j, k, dummy;
  // #pragma acc data copy(a)
   #pragma acc kernels
   #pragma acc loop auto
   for (k = 0; k < n; k++)
   {
      #pragma acc loop auto  
      for (i = 0; i < n; i++){
         #pragma acc loop auto
         for (j = 0; j < n; j++)
         	a[i][j]= MIN(a[i][j], a[i][k] + a[k][j]);
      }
   }	
}
