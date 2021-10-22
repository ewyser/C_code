#include "math.h"
#include "../include/NdN.h"
void NdN(DAT *N_dN, DAT xi, int type){
    N_dN[0] = N_dN[1] = (DAT)0.0;
    if(type==1){
        if((-2.0<=xi) && (xi<=-1.0)){
            N_dN[0] = 1.0/6.0*(xi*xi*xi)+1.0*(xi*xi)+2.0*xi+4.0/3.0;
            N_dN[1] = 3.0/6.0*(xi*xi   )+2.0*(xi   )+2.0           ;
        }
        else if((-1.0<=xi) && (xi<=0.0)){
            N_dN[0] = -1.0/6.0*(xi*xi*xi)           +1.0*xi+1.0;
            N_dN[1] = -3.0/6.0*(xi*xi   )           +1.0       ;
        }    
        else if( (0.0<=xi) && (xi<=1.0)){
            N_dN[0] =  1.0/6.0*(xi*xi*xi)           -1.0*xi+1.0;
            N_dN[1] =  3.0/6.0*(xi*xi   )           -1.0       ;
        }
        else if( (1.0<=xi) && (xi<=2.0)){
            N_dN[0] = -1.0/6.0*(xi*xi*xi)+1.0*(xi*xi)-2.0*xi+4.0/3.0;
            N_dN[1] = -3.0/6.0*(xi*xi   )+2.0*(xi   )-2.0           ;
        }    
    }












    else if(type==2){
        if((-1.0<=xi) && (xi<=0.0)){
            N_dN[0] = -1.0/3.0*(xi*xi*xi)-1.0*(xi*xi)       +2.0/3.0;
            N_dN[1] = -3.0/3.0*(xi*xi   )-2.0*(xi   )               ;
        }
        else if((0.0<=xi) && (xi<=1.0)){
            N_dN[0] =  1.0/2.0*(xi*xi*xi)-1.0*(xi*xi)       +2.0/3.0;
            N_dN[1] =  3.0/2.0*(xi*xi   )-2.0*(xi   )               ;
        }
        else if((1.0<=xi) && (xi<=2.0)){
            N_dN[0] = -1.0/6.0*(xi*xi*xi)+1.0*(xi*xi)-2.0*xi+4.0/3.0;
            N_dN[1] = -3.0/6.0*(xi*xi   )+2.0*(xi   )-2.0           ;
        }    
    }














    else if(type==3){
        if((-2.0<=xi) && (xi<=-1.0)){
            N_dN[0] = 1.0/6.0*(xi*xi*xi)+1.0*(xi*xi)+2.0*xi+4.0/3.0;
            N_dN[1] = 3.0/6.0*(xi*xi   )+2.0*(xi   )+2.0           ;
        }
        else if((-1.0<=xi) && (xi<=0.0)){
            N_dN[0] = -1.0/2.0*(xi*xi*xi)-1.0*(xi*xi)+2.0/3.0;
            N_dN[1] = -3.0/2.0*(xi*xi   )-2.0*(xi   )        ;
        }
        else if( (0.0<=xi) && (xi<=1.0)){
            N_dN[0] = 1.0/2.0*(xi*xi*xi)-1.0*(xi*xi)+2.0/3.0;
            N_dN[1] = 3.0/2.0*(xi*xi   )-2.0*(xi   )        ;
        }            
        else if( (1.0<=xi) && (xi<=2.0)){
            N_dN[0] = -1.0/6.0*(xi*xi*xi)+1.0*(xi*xi)-2.0*xi+4.0/3.0;
            N_dN[1] = -3.0/6.0*(xi*xi   )+2.0*(xi   )-2.0           ;            
        }            
    }















    else if(type==4){
            if((-2.0<=xi) && (xi<=-1.0)){
            N_dN[0] =  1.0/6.0*(xi*xi*xi)+1.0*(xi*xi)+2.0*xi+4.0/3.0;
            N_dN[1] =  3.0/6.0*(xi*xi   )+2.0*(xi   )+2.0           ;
        }
        else if((-1.0<=xi) && (xi<=0.0)){
            N_dN[0] = -1.0/2.0*(xi*xi*xi)-1.0*(xi*xi)+2.0/3.0;
            N_dN[1] = -3.0/2.0*(xi*xi   )-2.0*xi           ;
        }
        else if( (0.0<=xi) && (xi<=1.0)){
            N_dN[0] =  1.0/3.0*(xi*xi*xi)-1.0*(xi*xi)+2.0/3.0;
            N_dN[1] =  3.0/3.0*(xi*xi   )-2.0*(xi   )        ;
        }        
    }
    else{
            N_dN[0] = 0.0;
            N_dN[1] = 0.0;
    }
}