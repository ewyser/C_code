#include "math.h"
#include "../include/type_t.h"
#include "../include/DPPlast.h"
void DPPlast(point_t *mpD, DAT Hp, DAT cohr, DAT Kc, DAT Gc, DAT psi){
    int nmp = mpD->nmp;
    DAT c,Pr,tensile,sxx,syy,szz,sxy,syz,sxz,devxx,devyy,devzz,devxy,devyz,devxz,tau,eta,etaB,xi,sigm,fs,ft,tauP,alpP,h,dlam,tauN,PrN,dep;
    tensile = (DAT)0.0;
    for(int p=0; p<nmp; p++){
        c  = mpD->cohp[p]+(Hp*mpD->epII[p]);
        if(c<cohr){c = cohr;}
        sxx   = mpD->sig[0+6*p] ;
        syy   = mpD->sig[1+6*p] ;
        szz   = mpD->sig[2+6*p] ;
        sxy   = mpD->sig[3+6*p] ;
        syz   = mpD->sig[4+6*p] ;
        sxz   = mpD->sig[5+6*p] ;
        Pr    = (sxx+syy+szz)*((DAT)1.0/(DAT)3.0);
        devxx = sxx - Pr     ; 
        devyy = syy - Pr     ;
        devzz = szz - Pr     ;
        devxy = sxy          ;
        devyz = syz          ;
        devxz = sxz          ;
        tau   = sqrt((DAT)0.5*(devxx*devxx+devyy*devyy+devzz*devzz)+devxy*devxy+devyz*devyz+devxz*devxz);
        eta   = ((DAT)6.0*  sin(mpD->phip[p]))*(DAT)1.0/((DAT)sqrt((DAT)3.0)*((DAT)3.0+sin(mpD->phip[p])));
        etaB  = ((DAT)6.0*  sin(psi         ))*(DAT)1.0/((DAT)sqrt((DAT)3.0)*((DAT)3.0+sin(psi         )));
        xi    = ((DAT)6.0*c*cos(mpD->phip[p]))*(DAT)1.0/((DAT)sqrt((DAT)3.0)*((DAT)3.0+sin(mpD->phip[p])));        
    
        
        //sigm  = min(tensile,xi*(DAT)1.0/((DAT)eta));
        sigm  = xi*(DAT)1.0/((DAT)eta);
        fs    = tau+eta*Pr-xi;
        ft    = Pr-sigm         ;
        tauP  = xi-eta*sigm  ;
        alpP  = sqrt((DAT)1.0+eta*eta)-eta;
        h     = tau-tauP-alpP*(Pr-sigm);
        if((fs>(DAT)0.0 && Pr<sigm)||(h>(DAT)0.0 && Pr>=sigm)){
            dlam             = fs*(DAT)1.0/((DAT)Gc+Kc*eta*etaB)           ;
            PrN              = Pr-Kc*etaB*dlam                             ;
            tauN             = xi-eta*PrN                                  ;
            mpD->sig [0+6*p] = devxx*((DAT)tauN/((DAT)tau))+PrN            ;
            mpD->sig [1+6*p] = devyy*((DAT)tauN/((DAT)tau))+PrN            ;
            mpD->sig [2+6*p] = devzz*((DAT)tauN/((DAT)tau))+PrN            ;
            mpD->sig [3+6*p] = devxy*((DAT)tauN/((DAT)tau))                ;
            mpD->sig [4+6*p] = devyz*((DAT)tauN/((DAT)tau))                ;
            mpD->sig [5+6*p] = devxz*((DAT)tauN/((DAT)tau))                ;
            dep              = (dlam*sqrt((DAT)0.33333333+(DAT)0.22222222*etaB*etaB));
            mpD->epII[p    ]+= dep                                         ;
        }
        if((h<=0.0)&&(Pr>=sigm)){
            dlam             = (Pr-sigm)*((DAT)1.0/((DAT)Kc))        ;
            mpD->sig[0+6*p] += (sigm-Pr)                             ;
            mpD->sig[1+6*p] += (sigm-Pr)                             ;
            mpD->sig[2+6*p] += (sigm-Pr)                             ;
            mpD->sig[3+6*p] += (DAT)0.0                              ;
            mpD->sig[4+6*p] += (DAT)0.0                              ;
            mpD->sig[5+6*p] += (DAT)0.0                              ;
            dep              = (sqrt(2.0)*dlam*((DAT)1.0/((DAT)3.0)));
            mpD->epII[p    ]+= dep                                   ;
        }  
    }
}