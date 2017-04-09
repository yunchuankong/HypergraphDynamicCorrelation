
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>

/*
void make3(int *tuples, int *x, int *theOne, int *nog, int *not){ 
// Only record gene indices, ignore LA scores
    int i,j,row=0,n,where;
	where=*theOne;
    n=*nog-where-1;
	for (i=0;i<n;i++){ // i+1 > *theOne
		for (j=0;j<=i;j++){
			if (x[ i*(i+1)/2+j ] == 1){
                tuples[row*2] = j+1+where;
                tuples[row*2+1] = i+2+where;
                row++;
			}
		}
	}
}
*/

void make3(int *tuples, int *x, int *theOne, int *nog, int *not){ 
// Only record gene indices, ignore LA scores
    int i,j,row=0,n;
	n=*nog-1;
	for (i=*theOne;i<n;i++){ // i+1 > *theOne
		for (j=*theOne;j<=i;j++){
			if (x[ i*(i+1)/2+j ] == 1){
                tuples[row*2] = j+1;
                tuples[row*2+1] = i+2;
                row++;
			}
		}
	}
}


int cmpfunc (const void * a, const void * b) // for qsort()
{
   return ( *(int*)a - *(int*)b );
}

void count(int *tri, int *n, int *flow, int *N){
	int i, nog, ind[3];
	nog=*N;
	for (i=0;i<*n;i++){
		ind[0]=tri[i*3+0]-1;
		ind[1]=tri[i*3+1]-1;
		ind[2]=tri[i*3+2]-1;
		// qsort(ind,3,sizeof(int),cmpfunc);
		flow[ ind[0]*nog+ind[1] ]+=1;
		flow[ ind[0]*nog+ind[2] ]+=1;
		flow[ ind[1]*nog+ind[0] ]+=1;
		flow[ ind[1]*nog+ind[2] ]+=1;
		flow[ ind[2]*nog+ind[0] ]+=1;
		flow[ ind[2]*nog+ind[1] ]+=1;
	}
}

void shrink(int *tri, int *n, int *adj, int *lab, int *k){
    int i, x, y, z, m;
    int index[3]; 
	m=*k; // # of clusters + 1
    for(i = 0; i < *n; i++){ // enumerate each row of triplets
		for (x=0;x<(m-1);x++){
			if (lab[ (tri[i*3+0]-1)*m+x ]!=0){
				for (y=0;y<(m-1);y++){
					if (lab[ (tri[i*3+1]-1)*m+y ]!=0){
						for (z=0;z<(m-1);z++){
							if (lab[ (tri[i*3+2]-1)*m+z ]!=0){
								index[0]=x;
								index[1]=y;
								index[2]=z;
								qsort(index,3,sizeof(int),cmpfunc);
								adj[ index[0]*(m-1)*(m-1)+index[1]*(m-1)+index[2] ]+=1;	
								if (i%100000==0){
									printf("shirinking %d: %d,%d,%d\n",i,index[0],index[1],index[2]);
								}
							}
						} // for z
					}
				} // for y
			}
		} // for x
    } // for i
}

void fullyConnect(int *adj, int *labSum, int *k){
    int m, x, y, z;
    m=*k; // # of clusters
    for (x=0;x<m;x++){
        for (y=0;y<=x;y++){
            for (z=0;z<=y;z++){
                if (x!=y && y!=z){
                    adj[ z*m*m+y*m+x ]=labSum[x]*labSum[y]*labSum[z];
                }
                else if (x==y && x!=z){
                    adj[ z*m*m+y*m+x ]=labSum[x]*(labSum[x]-1)/2*labSum[z];   
                }
                else if (x==z && x!=y){
                    adj[ z*m*m+y*m+x ]=labSum[x]*(labSum[x]-1)/2*labSum[y];   
                }
                else if (x!=y && y==z){
                    adj[ z*m*m+y*m+x ]=labSum[y]*(labSum[y]-1)/2*labSum[x];   
                }
                else {
                    adj[ z*m*m+y*m+x ]=labSum[x]*(labSum[x]-1)*(labSum[x]-2)/6;
                }
            } // for z
        } // for y
    } // for x    
}

void normalize(double *net, double *full, int *k, double *adj){
    int x, y, z, m;
    m=*k; // # of clusters
    for (x=0;x<m;x++){
        for (y=0;y<=x;y++){
            for (z=0;z<=y;z++){
                if (full[ z*m*m+y*m+x ]>0.1){ // comparision to zero of double type
                    adj[ z*m*m+y*m+x ]=net[ z*m*m+y*m+x ]/full[ z*m*m+y*m+x ];
                }
            }
        }
    }        
}



