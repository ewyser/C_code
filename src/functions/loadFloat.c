#include "stdlib.h"
#include "stdio.h"
DAT* loadFloat(int size, char* filename){
        // creater pointer to file & open file
        DAT *ptr=malloc(size*sizeof(DAT));
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
            if(sizeof(DAT)==4){
                fscanf(fid,"%f",&ptr[i]);    
            }
            else if(sizeof(DAT)==8){
                fscanf(fid,"%lf",&ptr[i]);
            }
            
        }
        // clear & free
        fclose(fid);
        return(ptr);    
}