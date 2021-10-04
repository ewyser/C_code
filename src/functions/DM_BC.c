#include "../include/type_t.h"
void DM_BC(mesh_t* meD, point_t* mpD, DAT dt){
    int id,iDx, iDy, iDz;
    int nmp = mpD->nmp, 
         nn = meD->nn, 
         no = meD->nno[3];
    DAT cache,m,duxp,duyp,duzp;
    // initialize nodal momentum
    for(int n=0;n<(3*no);n++){
        meD->pn[n] = meD->un[n] = (DAT)0.0;
    }
    // accumulate material point momentum
    for(int n=0;n<nn;n++){
        for(int p=0;p<nmp;p++){
            // indexing & caching
            id            = p+n*nmp                 ;
            iDx           = mpD->p2n[id]+0*no       ;
            iDy           = mpD->p2n[id]+1*no       ;
            iDz           = mpD->p2n[id]+2*no       ;
            cache         = mpD->N[id]*mpD->mp[p]   ;
            // accumulation
            meD->pn[iDx] += (cache*mpD->vp[p+0*nmp]);
            meD->pn[iDy] += (cache*mpD->vp[p+1*nmp]);
            meD->pn[iDz] += (cache*mpD->vp[p+2*nmp]);
        }
    }
    // solve for nodal incremental displacement
    for(int n=0;n<no;n++){
        iDx = n+0*no ;
        iDy = n+1*no ;
        iDz = n+2*no ;
        if(meD->mn[n]>(DAT)0.0){
            m            = ((DAT)1.0/(DAT)meD->mn[n])         ;
            meD->un[iDx] = dt*meD->pn[iDx]*m*(DAT)meD->bc[iDx];
            meD->un[iDy] = dt*meD->pn[iDy]*m*(DAT)meD->bc[iDy];
            meD->un[iDz] = dt*meD->pn[iDz]*m*(DAT)meD->bc[iDz];
        }
    }
    // update material point displacement
    for(int p=0;p<nmp;p++){
        duxp = duyp = duzp = (DAT)0.0;
        for(int n=0;n<nn;n++){
            iDx   = mpD->p2n[p+n*nmp]+0*no        ;
            iDy   = mpD->p2n[p+n*nmp]+1*no        ;
            iDz   = mpD->p2n[p+n*nmp]+2*no        ;
            duxp += (mpD->N[p+n*nmp]*meD->un[iDx]);
            duyp += (mpD->N[p+n*nmp]*meD->un[iDy]);
            duzp += (mpD->N[p+n*nmp]*meD->un[iDz]);
        }
        mpD->up[p+0*nmp] += (duxp);
        mpD->up[p+1*nmp] += (duyp);
        mpD->up[p+2*nmp] += (duzp);
    } 
}