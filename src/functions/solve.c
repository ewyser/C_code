#include "math.h"
#include "../include/type_t.h"
void solve(mesh_t* meD, DAT dt){
    int iDx,iDy,iDz;
    int nmp = mpD.nmp, 
         nn = meD->nn, 
         no = meD->nno[3];
    DAT dmp,m,fx,fy,fz,px,py,pz,vx,vy,vz;
    // initialize
    for(int n=0;n<3*no;n++){
        meD->fn[n] = meD->an[n] = meD->vn[n] = 0.0;
    }
    // solve momentum equation on the background mesh
    for(int n=0;n<no;n++){
        iDx = n+0*no;
        iDy = n+1*no;
        iDz = n+2*no;
        if(meD->mn[n]>0.0){
            fx = (meD->fen[iDx]-meD->fin[iDx]);
            fy = (meD->fen[iDy]-meD->fin[iDy]);
            fz = (meD->fen[iDz]-meD->fin[iDz]);
            px = meD->pn[iDx];
            py = meD->pn[iDy];
            pz = meD->pn[iDz];
            vx = px*((DAT)1.0)/((DAT)meD->mn[n]);
            vy = py*((DAT)1.0)/((DAT)meD->mn[n]);
            vz = pz*((DAT)1.0)/((DAT)meD->mn[n]);
            dmp= sqrt(fx*fx+fy*fy+fz*fz);
            fx = fx-D*dmp*((DAT)vx/(DAT)fabs(vx));
            if(fabs(vx)<1E-7){
                fx = (meD->fen[iDx]-meD->fin[iDx]);
            }
            fy = fy-D*dmp*((DAT)vy/(DAT)fabs(vy));
            if(fabs(vy)<1E-7){
                fy = (meD->fen[iDy]-meD->fin[iDy]);
            }
            fz = fz-D*dmp*((DAT)vz/(DAT)fabs(vz));
            if(fabs(vz)<1E-7){
                fz = (meD->fen[iDz]-meD->fin[iDz]);
            }
            m            = ((DAT)1.0)/((DAT)meD->mn[n])  ;
            meD->an[iDx] = (fx      )*m*(DAT)meD->bc[iDx];
            meD->an[iDy] = (fy      )*m*(DAT)meD->bc[iDy];
            meD->an[iDz] = (fz      )*m*(DAT)meD->bc[iDz];
            meD->vn[iDx] = (px+dt*fx)*m*(DAT)meD->bc[iDx];
            meD->vn[iDy] = (py+dt*fy)*m*(DAT)meD->bc[iDy];
            meD->vn[iDz] = (pz+dt*fz)*m*(DAT)meD->bc[iDz];
        }
    }
}