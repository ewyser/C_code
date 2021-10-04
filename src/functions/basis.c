#include "math.h"
#include "../include/type_t.h"
#include "../include/basis.h"
#include "../include/NdN.h"
void basis(mesh_t* meD, point_t* mpD){
    DAT Nx_dNx[2],Ny_dNy[2],Nz_dNz[2];
    int iD;
    int nmp = mpD->nmp, 
         nn = meD->nn, 
         no = meD->nno[3];
    for(int n=0;n<nn;n++){
        for(int p=0;p<nmp;p++){
            iD  = p+n*nmp;       
            // x-basis and x-derivative
            NdN(Nx_dNx,(mpD->xp[p+0*nmp]-meD->xn[mpD->p2n[iD]+0*no]),mpD->lp[p+0*nmp],meD->h[0]);
            // y-basis and y-derivative   
            NdN(Ny_dNy,(mpD->xp[p+1*nmp]-meD->xn[mpD->p2n[iD]+1*no]),mpD->lp[p+1*nmp],meD->h[1]);
            // y-basis and y-derivative   
            NdN(Nz_dNz,(mpD->xp[p+2*nmp]-meD->xn[mpD->p2n[iD]+2*no]),mpD->lp[p+2*nmp],meD->h[2]);
            // convolution of basis
            mpD->N[iD]   = Nx_dNx[0]*Ny_dNy[0]*Nz_dNz[0];
            mpD->dNx[iD] = Nx_dNx[1]*Ny_dNy[0]*Nz_dNz[0];
            mpD->dNy[iD] = Nx_dNx[0]*Ny_dNy[1]*Nz_dNz[0];
            mpD->dNz[iD] = Nx_dNx[0]*Ny_dNy[0]*Nz_dNz[1];
        }
    }
}