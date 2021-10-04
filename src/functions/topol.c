#include "math.h"
#include "../include/type_t.h"
#include "../include/topol.h"
void topol(mesh_t* meD,point_t* mpD){
    DAT dx    = meD->h[0];
    DAT dy    = meD->h[1];
    DAT dz    = meD->h[2];
    DAT xnmin = meD->min[0];
    DAT ynmin = meD->min[1];
    DAT znmin = meD->min[2];
    int nex   = meD->nel[0];
    int nez   = meD->nel[2];
    int nel   = meD->nel[3];
    int nn    = meD->nn;
    int nmp   = mpD->nmp;
    int id;    
    // caching
    DAT one_dx = ((DAT)1.0/(DAT)dx);
    DAT one_dy = ((DAT)1.0/(DAT)dy);
    DAT one_dz = ((DAT)1.0/(DAT)dz);
    // mp-to-element-to-node topology
    for(int i=0;i<nmp;i++){
        mpD->p2e[i] = (int)((floor((mpD->xp[i+2*nmp]-znmin)*one_dz))+(nez)*floor((mpD->xp[i+0*nmp]-xnmin)*one_dx)+nez*nex*(floor((mpD->xp[i+1*nmp]-ynmin)*one_dy)));
        id          = mpD->p2e[i];
        for(int j=0;j<nn;j++){
            mpD->p2n[i+j*nmp] = meD->e2n[id+j*nel];
        }
    }
}