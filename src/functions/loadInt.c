#include "stdlib.h"
#include "stdio.h"
int* loadInt(int size, char* filename){
        // creater pointer to file & open file
        int *ptr=malloc(size*sizeof(int));
        FILE  *fid = fopen(filename, "r");
        if(fid==NULL){
            printf("\n! %s not found -> program killed !\n",filename);
            fclose(fid);
            exit(-1);
        }
        // read data from file & store value at &par
        for(int i=0;i<size;i++){
            fscanf(fid,"%d" ,&ptr[i]);
        }
        // clear & free
        fclose(fid);
        return(ptr);    
}
// if fp == 0 then read integer and if fp == 1 then read float