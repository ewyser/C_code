#include "math.h"
#include "../include/type_t.h"
#include "../include/strains.h"
void strains(mesh_t* meD, point_t* mpD,DAT dt){
    int id,iDx,iDy,iDz;
    int nmp = mpD->nmp, 
         nn = meD->nn, 
         no = meD->nno[3];
    DAT J;
    // update material point state variables
    for(int p=0;p<nmp;p++){
        // compute incremental deformation gradient
        mpD->dF[0+p*9] = mpD->dF[4+p*9] = mpD->dF[8+p*9] = (DAT)1.0;
        mpD->dF[1+p*9] = mpD->dF[2+p*9] = mpD->dF[3+p*9] = mpD->dF[5+p*9] = mpD->dF[6+p*9] = mpD->dF[7+p*9] = (DAT)0.0;
        for(int n=0;n<nn;n++){
            id              = p+n*nmp          ;
            iDx             = mpD->p2n[id]+0*no     ;
            iDy             = mpD->p2n[id]+1*no     ;
            iDz             = mpD->p2n[id]+2*no     ;
            mpD->dF[0+p*9] += (mpD->dNx[id]*meD->un[iDx]);
            mpD->dF[1+p*9] += (mpD->dNy[id]*meD->un[iDx]);
            mpD->dF[2+p*9] += (mpD->dNz[id]*meD->un[iDx]);
            
            mpD->dF[3+p*9] += (mpD->dNx[id]*meD->un[iDy]);
            mpD->dF[4+p*9] += (mpD->dNy[id]*meD->un[iDy]);
            mpD->dF[5+p*9] += (mpD->dNz[id]*meD->un[iDy]);
            
            mpD->dF[6+p*9] += (mpD->dNx[id]*meD->un[iDz]);
            mpD->dF[7+p*9] += (mpD->dNy[id]*meD->un[iDz]);
            mpD->dF[8+p*9] += (mpD->dNz[id]*meD->un[iDz]);
        }
        // compute incremental strains
        mpD->eps[0+p*6] = mpD->dF[0+p*9]-(DAT)1.0;
        mpD->eps[1+p*6] = mpD->dF[4+p*9]-(DAT)1.0;
        mpD->eps[2+p*6] = mpD->dF[8+p*9]-(DAT)1.0;
        mpD->eps[3+p*6] = mpD->dF[1+p*9]+mpD->dF[3+p*9];
        mpD->eps[4+p*6] = mpD->dF[7+p*9]+mpD->dF[5+p*9];
        mpD->eps[5+p*6] = mpD->dF[6+p*9]+mpD->dF[2+p*9];
        // compute incremental rotation
        mpD->ome[0+p*3] = (DAT)0.5*(mpD->dF[3+p*9]-mpD->dF[1+p*9]);
        mpD->ome[1+p*3] = (DAT)0.5*(mpD->dF[5+p*9]-mpD->dF[7+p*9]);
        mpD->ome[2+p*3] = (DAT)0.5*(mpD->dF[2+p*9]-mpD->dF[6+p*9]);
        // update material point volume and domain lengths
        J          = (DAT)1.0+(mpD->eps[0+p*6]+mpD->eps[1+p*6]+mpD->eps[2+p*6]);
        mpD->vol[p]= J*mpD->vol[p];
        if(sizeof(DAT)==(int)8){
            J           = pow(J,(DAT)0.3333);
        }
        else if(sizeof(DAT)==(int)4){
            J           = powf(J,(DAT)0.3333);
        }
        mpD->lp[p+0*nmp] = J*mpD->lp[p+0*nmp];
        mpD->lp[p+1*nmp] = J*mpD->lp[p+1*nmp];
        mpD->lp[p+2*nmp] = J*mpD->lp[p+2*nmp];
    }
}