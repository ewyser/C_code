#include "math.h"
#include "../include/type_t.h"
#include "../include/basis.h"
#include "../include/whichType.h"
#include "../include/NdN.h"
void basis(mesh_t *meD, point_t *mpD){
    DAT Nx_dNx[2],Ny_dNy[2],Nz_dNz[2],xmin,xmax;
    int iD,type;
    int nmp = mpD->nmp, 
         nn = meD->nn, 
         no = meD->nno[3];
    for(int n=0;n<nn;n++){
        for(int p=0;p<nmp;p++){
            iD  = p+n*nmp;       
            // x-basis and x-derivative
            type = 0;
            xmin = meD->xB[0];
            xmax = meD->xB[1];
            type = whichType(meD->xn[mpD->p2n[iD]+0*no], xmin, xmax, meD->h[0]);
            NdN(Nx_dNx,(mpD->xp[p+0*nmp]-meD->xn[mpD->p2n[iD]+0*no])*(DAT)1.0/(DAT)meD->h[0],type);
            // y-basis and y-derivative
            type = 0;
            xmin = meD->xB[2];
            xmax = meD->xB[3];
            type = whichType(meD->xn[mpD->p2n[iD]+1*no], xmin, xmax, meD->h[1]);
            NdN(Ny_dNy,(mpD->xp[p+1*nmp]-meD->xn[mpD->p2n[iD]+1*no])*(DAT)1.0/(DAT)meD->h[1],type);
            // z-basis and y-derivative
            type = 0;
            xmin = meD->xB[4];
            xmax = meD->xB[5];
            type = whichType(meD->xn[mpD->p2n[iD]+2*no], xmin, xmax, meD->h[2]);   
            NdN(Nz_dNz,(mpD->xp[p+2*nmp]-meD->xn[mpD->p2n[iD]+2*no])*(DAT)1.0/(DAT)meD->h[2],type);
            // convolution of basis
            mpD->N[iD]   = Nx_dNx[0]*Ny_dNy[0]*Nz_dNz[0];
            mpD->dNx[iD] = Nx_dNx[1]*Ny_dNy[0]*Nz_dNz[0];
            mpD->dNy[iD] = Nx_dNx[0]*Ny_dNy[1]*Nz_dNz[0];
            mpD->dNz[iD] = Nx_dNx[0]*Ny_dNy[0]*Nz_dNz[1];
        }
    }
}