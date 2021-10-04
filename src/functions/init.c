#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "../include/type_t.h"
#include "../include/macros.h"

void init(mesh_t* meD,point_t* mpD){
    // mesh & material point geometric properties
        DAT  *par = (DAT*)malloc(12*sizeof(DAT));
        FILE *fid = fopen("param.dat", "rb");                                                      ;\
                    fread(par, sizeof(DAT),12,fid); 

        mpD->nmp   = (int) par[0 ];
        meD->nn    = (int) par[1 ];
        meD->nno[3]= (int) par[2 ];
        meD->h[0]  = (DAT) par[3 ];
        meD->h[1]  = (DAT) par[4 ];
        meD->h[2]  = (DAT) par[5 ];
        meD->min[0]= (DAT) par[6 ];
        meD->min[1]= (DAT) par[7 ];
        meD->min[2]= (DAT) par[8 ];
        meD->nno[0]= (int) par[9 ];
        meD->nno[1]= (int) par[10];
        meD->nno[2]= (int) par[11];
        meD->nno[3]= meD->nno[0]*meD->nno[1]*meD->nno[2];
        meD->nel[0]= meD->nno[0]-1;
        meD->nel[1]= meD->nno[1]-1;
        meD->nel[2]= meD->nno[2]-1;
        meD->nel[3]= meD->nel[0]*meD->nel[1]*meD->nel[2];  
        free(par);
        fclose(fid);
    // mesh
        zero_t(meD,pel            ,2*meD->nno[3]      ,DAT);    
        load_t(meD,e2n ,"e2n.dat" ,meD->nn*meD->nel[3],int);
        zero_t(meD,mn             ,  meD->nno[3]      ,DAT);
        load_t(meD,xn  ,"xn.dat"  ,3*meD->nno[3]      ,DAT);
        zero_t(meD,pn             ,3*meD->nno[3]      ,DAT);
        zero_t(meD,pn             ,3*meD->nno[3]      ,DAT);
        zero_t(meD,fen            ,3*meD->nno[3]      ,DAT);
        zero_t(meD,fin            ,3*meD->nno[3]      ,DAT);
        zero_t(meD,fn             ,3*meD->nno[3]      ,DAT);
        zero_t(meD,an             ,3*meD->nno[3]      ,DAT);
        zero_t(meD,vn             ,3*meD->nno[3]      ,DAT);
        zero_t(meD,un             ,3*meD->nno[3]      ,DAT);
        load_t(meD,bc  ,"bcs.dat" ,3*meD->nno[3]      ,int);
    // material point
        load_t(mpD,mp  ,"mp.dat"  ,mpD->nmp           ,DAT);
        load_t(mpD,vol ,"vol.dat" ,mpD->nmp           ,DAT);
        load_t(mpD,cohp,"cohp.dat",mpD->nmp           ,DAT);
        load_t(mpD,phip,"phip.dat",mpD->nmp           ,DAT);
        zero_t(mpD,epII           ,mpD->nmp           ,DAT);

        load_t(mpD,xp  ,"xp.dat"  ,3*mpD->nmp         ,DAT);
        zero_t(mpD,vp             ,3*mpD->nmp         ,DAT);    
        zero_t(mpD,up             ,3*mpD->nmp         ,DAT);    
        load_t(mpD,lp  ,"lp.dat"  ,3*mpD->nmp         ,DAT);   

        zero_t(mpD,sig            ,6*mpD->nmp         ,DAT);
        zero_t(mpD,eps            ,6*mpD->nmp         ,DAT);
        zero_t(mpD,dev            ,6*mpD->nmp         ,DAT);
        zero_t(mpD,dF             ,9*mpD->nmp,         DAT);
        zero_t(mpD,ome            ,3*mpD->nmp         ,DAT);

        zero_t(mpD,p2e            ,        mpD->nmp   ,int);
        zero_t(mpD,p2n            ,meD->nn*mpD->nmp   ,int);   

        zero_t(mpD,N              ,meD->nn*mpD->nmp   ,DAT);
        zero_t(mpD,dNx            ,meD->nn*mpD->nmp   ,DAT);
        zero_t(mpD,dNy            ,meD->nn*mpD->nmp   ,DAT);
        zero_t(mpD,dNz            ,meD->nn*mpD->nmp   ,DAT);         
    /*
    for(int k=0; k<nmp; k++){
        printf("\n xp = [%f,%f,%f]",mpD->xp[k+0*nmp],mpD->xp[k+1*mpD->nmp],mpD->xp[k+2*mpD->nmp]);
    } 
    */                                                                  
}