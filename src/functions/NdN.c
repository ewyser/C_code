#include "math.h"
#include "../include/NdN.h"
void NdN(DAT* N_dN,DAT xi,DAT lp,DAT dx){
    N_dN[0] = N_dN[1] = (DAT)0.0;
    if( fabs(xi)< (   lp)                     ){N_dN[0]=1.0-(4.0*(xi*xi)+((2.0*lp)*(2.0*lp)))*((DAT)1.0/((DAT)8.0*dx*lp));
                                                N_dN[1]=(-8.0*xi)*((DAT)1.0/((DAT)8.0*dx*lp));
    }
    if((fabs(xi)>=(   lp))&&(fabs(xi)<(dx-lp))){N_dN[0]=1.0-(fabs(xi)*((DAT)1.0/(DAT)dx));
                                                if(xi<0.0){N_dN[1]=( (DAT)1.0/(DAT)dx);}
                                                if(xi>0.0){N_dN[1]=(-(DAT)1.0/(DAT)dx);}
    }
    if((fabs(xi)>=(dx-lp))&&(fabs(xi)<(dx+lp))){N_dN[0]=((dx+lp-fabs(xi))*(dx+lp-fabs(xi)))*((DAT)1.0/((DAT)4.0*dx*lp));
                                                if(xi<0.0){N_dN[1]=(  dx+lp-fabs(xi))*((DAT)1.0/((DAT)2.0*dx*lp)) ;}
                                                if(xi>0.0){N_dN[1]=(-(dx+lp-fabs(xi))*((DAT)1.0/((DAT)2.0*dx*lp)));}
    }
}