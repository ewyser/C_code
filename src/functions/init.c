#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "../include/type_t.h"
#include "../include/macros.h"
#include "../include/init.h"
#include "../include/loadInt.h"
#include "../include/loadFloat.h"

void init(mesh_t *meD,point_t *mpD){
    // mesh & material point geometric properties
        DAT *par   = malloc(18*sizeof(DAT));
        par        = loadFloat(18,"param.txt");
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

        meD->xB[0]= (DAT) par[12];
        meD->xB[1]= (DAT) par[13];
        meD->xB[2]= (DAT) par[14];
        meD->xB[3]= (DAT) par[15];
        meD->xB[4]= (DAT) par[16];
        meD->xB[5]= (DAT) par[13];

        meD->nno[3]= meD->nno[0]*meD->nno[1]*meD->nno[2];
        meD->nel[0]= meD->nno[0]-1;
        meD->nel[1]= meD->nno[1]-1;
        meD->nel[2]= meD->nno[2]-1;
        meD->nel[3]= meD->nel[0]*meD->nel[1]*meD->nel[2];  
        free(par);
    // mesh
        zero_t(meD,pel ,2*meD->nno[3]   ,DAT);    
        zero_t(meD,mn  ,1*meD->nno[3]   ,DAT);
        zero_t(meD,pn  ,3*meD->nno[3]   ,DAT);
        zero_t(meD,pn  ,3*meD->nno[3]   ,DAT);
        zero_t(meD,fen ,3*meD->nno[3]   ,DAT);
        zero_t(meD,fin ,3*meD->nno[3]   ,DAT);
        zero_t(meD,fn  ,3*meD->nno[3]   ,DAT);
        zero_t(meD,an  ,3*meD->nno[3]   ,DAT);
        zero_t(meD,vn  ,3*meD->nno[3]   ,DAT);
        zero_t(meD,un  ,3*meD->nno[3]   ,DAT);
        meD->e2n = malloc((meD->nn*meD->nel[3])*sizeof(int));
        meD->e2n = loadInt((meD->nn*meD->nel[3]),"e2n.txt");
        meD->bc  = malloc((3*meD->nno[3])*sizeof(int));
        meD->bc  = loadInt(3*meD->nno[3],"bcs.txt");
        meD->xn  = malloc((3*meD->nno[3])*sizeof(DAT));
        meD->xn  = loadFloat(3*meD->nno[3],"xn.txt");
    // material point
        zero_t(mpD,epII,mpD->nmp        ,DAT);
        zero_t(mpD,vp  ,3*mpD->nmp      ,DAT);    
        zero_t(mpD,up  ,3*mpD->nmp      ,DAT);     
        zero_t(mpD,sig ,6*mpD->nmp      ,DAT);
        zero_t(mpD,eps ,6*mpD->nmp      ,DAT);
        zero_t(mpD,dev ,6*mpD->nmp      ,DAT);
        zero_t(mpD,dF  ,9*mpD->nmp      ,DAT);
        zero_t(mpD,ome ,3*mpD->nmp      ,DAT);
        zero_t(mpD,p2e ,1*mpD->nmp      ,int);
        zero_t(mpD,p2n ,meD->nn*mpD->nmp,int);   
        zero_t(mpD,N   ,meD->nn*mpD->nmp,DAT);
        zero_t(mpD,dNx ,meD->nn*mpD->nmp,DAT);
        zero_t(mpD,dNy ,meD->nn*mpD->nmp,DAT);
        zero_t(mpD,dNz ,meD->nn*mpD->nmp,DAT);         
        mpD->mp  = malloc(mpD->nmp*sizeof(DAT));
        mpD->mp  = loadFloat(mpD->nmp,"mp.txt");
        mpD->vol = malloc(mpD->nmp*sizeof(DAT));
        mpD->vol = loadFloat(mpD->nmp,"vol.txt");
        mpD->cohp= malloc(mpD->nmp*sizeof(DAT));
        mpD->cohp= loadFloat(mpD->nmp,"cohp.txt");
        mpD->phip= malloc(mpD->nmp*sizeof(DAT));
        mpD->phip= loadFloat(mpD->nmp,"phip.txt");
        mpD->xp  = malloc(3*mpD->nmp*sizeof(DAT));
        mpD->xp  = loadFloat(3*mpD->nmp,"xp.txt");
        mpD->lp  = malloc(3*mpD->nmp*sizeof(DAT));
        mpD->lp  = loadFloat(3*mpD->nmp,"lp.txt");
}