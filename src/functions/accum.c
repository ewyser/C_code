#include "math.h"
#include "../include/type_t.h"
#include "../include/accum.h"
void accum(mesh_t* meD,point_t* mpD,DAT g){
    int id,iD, iDx, iDy, iDz;
    int nmp = mpD->nmp, 
         nn = meD->nn, 
         no = meD->nno[3];
    DAT cache;
    // initialize
    for(int n=0;n<(3*no);n++){
        if(n<no){
            meD->mn[n] = 0.0;
        }
        meD->pn [n] = meD->fen[n] = meD->fin[n] = 0.0;
    }
    // accumulate material point contributions on nodes
    for(int n=0;n<nn; n++){
        for(int p=0;p<nmp;p++){
            // indexing & caching
            id        = p+n*nmp;
            iD        = (mpD->p2n[id]);
            iDx       = iD+0*no;
            iDy       = iD+1*no;
            iDz       = iD+2*no;
            cache     = mpD->N[p+n*nmp]*mpD->mp[p];
            // accumulation
            meD->mn [iD ] += (cache);
            meD->pn [iDx] += (cache*mpD->vp[p+0*nmp]);
            meD->pn [iDy] += (cache*mpD->vp[p+1*nmp]);
            meD->pn [iDz] += (cache*mpD->vp[p+2*nmp]);
            meD->fen[iDz] -= (cache*g          );
            meD->fin[iDx] += mpD->vol[p]*(mpD->dNx[id]*mpD->sig[0+6*p]+mpD->dNy[id]*mpD->sig[3+6*p]+mpD->dNz[id]*mpD->sig[5+6*p]);
            meD->fin[iDy] += mpD->vol[p]*(mpD->dNx[id]*mpD->sig[3+6*p]+mpD->dNy[id]*mpD->sig[1+6*p]+mpD->dNz[id]*mpD->sig[4+6*p]);
            meD->fin[iDz] += mpD->vol[p]*(mpD->dNx[id]*mpD->sig[5+6*p]+mpD->dNy[id]*mpD->sig[4+6*p]+mpD->dNz[id]*mpD->sig[2+6*p]);
        }
    }
}