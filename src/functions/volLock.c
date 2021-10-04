#include "../include/type_t.h"
void volLock(mesh_t* meD, point_t* mpD){
    int iD;
    int nmp = mpD->nmp,
        nel = meD->nel[3];
    DAT pr;
    // initialize element pressure
    for(unsigned int e=0;e<nel;e++){
        meD->pel[0+2*e] = (DAT)0.0;
        meD->pel[1+2*e] = (DAT)0.0;
    }
    // accumulate material point pressure on element
    for(unsigned int p=0;p<nmp;p++){
        iD          = mpD->p2e[p]                 ;
        pr          = -(mpD->sig[0+6*p]+mpD->sig[1+6*p]+mpD->sig[2+6*p])/(DAT)3.0;
        mpD->dev[0+6*p ] = mpD->sig[0+6*p]+pr               ;
        mpD->dev[1+6*p ] = mpD->sig[1+6*p]+pr               ;
        mpD->dev[2+6*p ] = mpD->sig[2+6*p]+pr               ;
        mpD->dev[3+6*p ] = mpD->sig[3+6*p]+(DAT)0.0         ;
        mpD->dev[4+6*p ] = mpD->sig[4+6*p]+(DAT)0.0         ;
        mpD->dev[5+6*p ] = mpD->sig[5+6*p]+(DAT)0.0         ;
        meD->pel[0+2*iD]+= pr*mpD->vol[p]                   ;
        meD->pel[1+2*iD]+= mpD->vol[p]                      ;
    }
    // assign element pressure to material point
    for(unsigned int p=0;p<nmp;p++){
        iD          = mpD->p2e[p]                 ;
        pr          = meD->pel[0+2*iD]*1.0/meD->pel[1+2*iD] ;
        mpD->sig[0+6*p]  = mpD->dev[0+6*p]-pr               ;
        mpD->sig[1+6*p]  = mpD->dev[1+6*p]-pr               ;
        mpD->sig[2+6*p]  = mpD->dev[2+6*p]-pr               ;
        mpD->sig[3+6*p]  = mpD->dev[3+6*p]                  ;
        mpD->sig[4+6*p]  = mpD->dev[4+6*p]                  ;
        mpD->sig[5+6*p]  = mpD->dev[5+6*p]                  ;
    }
}