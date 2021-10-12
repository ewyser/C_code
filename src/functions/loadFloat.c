#include "stdlib.h"
#include "stdio.h"
DAT* loadFloat(int size, char* filename){
        // creater pointer to file & open file
        DAT *ptr=malloc(size*sizeof(DAT));
        FILE  *fid = fopen(filename, "r");
        if(fid==NULL){
            printf("\n! %s not found -> program killed !\n",filename);
            fclose(fid);
            exit(1);
        }
        // read data from file & store value at &par
        for(int i=0;i<size;i++){
            fscanf(fid,"%lf" ,&ptr[i]);
        }
        // clear & free
        fclose(fid);
        return(ptr);    
}
// if fp == 0 then read integer and if fp == 1 then read float