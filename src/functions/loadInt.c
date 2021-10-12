#include "stdlib.h"
#include "stdio.h"
int* loadInt(int size, char* filename){
        // creater pointer to file & open file
        int *ptr=malloc(size*sizeof(int));
        FILE  *fid = fopen(filename, "r");
        if(fid==NULL){
            printf("\n importing: %s not found ...",filename);
            printf("\n          |-> program killed");
            fclose(fid);
            printf("\n press any key to continue...");  
            getchar();  
            exit(1);
        }
        // read data from file & store value at &par
        printf("\n importing: %s",filename);
        for(int i=0;i<size;i++){
            fscanf(fid,"%d" ,&ptr[i]);
        }
        // clear & free
        fclose(fid);
        return(ptr);    
}