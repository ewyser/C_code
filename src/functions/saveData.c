#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "../include/saveData.h"
void saveData(DAT *data, char *filename, int dim1, int dim2){
     FILE *fp=fopen(filename,"w");
     
     if(dim2==1){
         for(int n=0;n<dim1;n++) {
             fprintf(fp,"%f\n",data[n]);
         }
     }
     else if(dim2==3){
         for(int n=0;n<dim1;n++) {
             fprintf(fp,"%f,%f,%f\n",data[n+0*dim1],data[n+1*dim1],data[n+2*dim1]);
         }
     }
     else if(dim2==6){
         for(int n=0;n<dim1;n++) {
             fprintf(fp,"%f,%f,%f,%f,%f,%f\n",data[0+dim2*n],data[1+dim2*n],data[2+dim2*n],data[3+dim2*n],data[4+dim2*n],data[5+dim2*n]);
         }        
     }
     fclose(fp);
}