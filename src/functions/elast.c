#include "../include/type_t.h"
void elast(point_t* mpD, DAT* Del, DAT dt){
    int nmp = mpD->nmp;
    DAT sxx0,syy0,szz0,sxy0,syz0,sxz0,dexx,deyy,dezz,dexy,deyz,dexz,doxy,doyz,doxz;
    for(int p=0;p<nmp;p++){
        // store previous stresses
        sxx0   = mpD->sig[0+6*p];
        syy0   = mpD->sig[1+6*p];
        szz0   = mpD->sig[2+6*p];
        sxy0   = mpD->sig[3+6*p];
        syz0   = mpD->sig[4+6*p];
        sxz0   = mpD->sig[5+6*p];
        // strain rate
        dexx   = mpD->eps[0+6*p]/(DAT)dt;
        deyy   = mpD->eps[1+6*p]/(DAT)dt;
        dezz   = mpD->eps[2+6*p]/(DAT)dt;
        dexy   = mpD->eps[3+6*p]/(DAT)dt;
        deyz   = mpD->eps[4+6*p]/(DAT)dt;
        dexz   = mpD->eps[5+6*p]/(DAT)dt;
        // spin
        doxy   = mpD->ome[0+3*p]/dt;
        doyz   = mpD->ome[1+3*p]/dt;
        doxz   = mpD->ome[2+3*p]/dt;
        // update objective stress        
        mpD->sig[0+6*p] +=  (DAT)2.0*dt*(sxy0*doxy+sxz0*doxz);
        mpD->sig[1+6*p] += -(DAT)2.0*dt*(sxy0*doxy-syz0*doyz);
        mpD->sig[2+6*p] += -(DAT)2.0*dt*(sxz0*doxz+syz0*doyz);
        mpD->sig[3+6*p] +=         dt*(doxy*(syy0-sxx0)+syz0*doxz+sxz0*doyz);
        mpD->sig[4+6*p] +=         dt*(doyz*(szz0-syy0)-sxy0*doxz-sxz0*doxy);
        mpD->sig[5+6*p] +=         dt*(doxz*(szz0-sxx0)+syz0*doxy-sxy0*doyz);
        // incremental strain
        mpD->sig[0+6*p] += dt*(Del[0]*dexx+Del[6 ]*deyy+Del[12]*dezz+Del[18]*dexy+Del[24]*deyz+Del[30]*dexz);
        mpD->sig[1+6*p] += dt*(Del[1]*dexx+Del[7 ]*deyy+Del[13]*dezz+Del[19]*dexy+Del[25]*deyz+Del[31]*dexz);
        mpD->sig[2+6*p] += dt*(Del[2]*dexx+Del[8 ]*deyy+Del[14]*dezz+Del[20]*dexy+Del[26]*deyz+Del[32]*dexz);
        mpD->sig[3+6*p] += dt*(Del[3]*dexx+Del[9 ]*deyy+Del[15]*dezz+Del[21]*dexy+Del[27]*deyz+Del[33]*dexz);
        mpD->sig[4+6*p] += dt*(Del[4]*dexx+Del[10]*deyy+Del[16]*dezz+Del[22]*dexy+Del[28]*deyz+Del[34]*dexz);
        mpD->sig[5+6*p] += dt*(Del[5]*dexx+Del[11]*deyy+Del[17]*dezz+Del[23]*dexy+Del[29]*deyz+Del[35]*dexz);   
    }
}