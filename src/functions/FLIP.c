#include "math.h"
#include "../include/type_t.h"
void FLIP(mesh_t* meD, point_t* mpD,DAT dt){
    // update material point velocity
    int id,iDx,iDy,iDz;
    int nmp = mpD->nmp, 
         nn = meD->nn, 
         no = meD->nno[3];
    DAT dvpx,dvpy,dvpz,dxp,dyp,dzp;
    for(int p=0;p<nmp;p++){
        dvpx = dvpy = dvpz = dxp = dyp = dzp = (DAT)0.0;
        for(int n=0;n<nn;n++){
            id    = p+n*nmp        ;
            iDx   = mpD->p2n[id]+0*no   ;
            iDy   = mpD->p2n[id]+1*no   ;
            iDz   = mpD->p2n[id]+2*no   ;
            dvpx += (mpD->N[id]*meD->an[iDx]);
            dvpy += (mpD->N[id]*meD->an[iDy]);
            dvpz += (mpD->N[id]*meD->an[iDz]);
            dxp  += (mpD->N[id]*meD->vn[iDx]);
            dyp  += (mpD->N[id]*meD->vn[iDy]);
            dzp  += (mpD->N[id]*meD->vn[iDz]);
        }
        mpD->vp[p+0*nmp] += (dt*dvpx);
        mpD->vp[p+1*nmp] += (dt*dvpy);
        mpD->vp[p+2*nmp] += (dt*dvpz);
        mpD->xp[p+0*nmp] += (dt*dxp );
        mpD->xp[p+1*nmp] += (dt*dyp );
        mpD->xp[p+2*nmp] += (dt*dzp );
    }   
}