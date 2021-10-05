#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "../include/type_t.h"
#include "../include/CFL.h"
DAT CFL(mesh_t *meD, point_t *mpD, DAT yd, DAT tg, DAT tw){
    int nmp = mpD->nmp;
    DAT cmax,vxpmax,vypmax,vzpmax,dt;
    DAT dx = meD->h[0],
        dy = meD->h[1],
        dz = meD->h[2];
    // determine suitable dt from CFL condition
    vxpmax = fabs(mpD->vp[0+0*nmp]);
    vypmax = fabs(mpD->vp[0+1*nmp]);
    vzpmax = fabs(mpD->vp[0+2*nmp]);
    for(int p=1;p<nmp;p++){
        if(vxpmax<fabs(mpD->vp[p+0*nmp])){
            vxpmax=fabs(mpD->vp[p+0*nmp]);
        }
        if(vypmax<fabs(mpD->vp[p+1*nmp])){
            vypmax=fabs(mpD->vp[p+1*nmp]);
        }
        if(vzpmax<fabs(mpD->vp[p+2*nmp])){
            vzpmax=fabs(mpD->vp[p+2*nmp]);
        }
    }
    cmax =  (DAT)dx/(DAT)(yd + vxpmax);
    if(cmax>((DAT)dy/(DAT)(yd + vypmax))){
        cmax =  (DAT)dy/(DAT)(yd + vypmax);
    }
    else if(cmax>((DAT)dz/(DAT)(yd + vzpmax))){
        cmax =  (DAT)dz/(DAT)(yd + vzpmax);
    }
    dt = 0.5*cmax;
    return(dt);
}