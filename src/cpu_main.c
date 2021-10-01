#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "time.h"


typedef struct point{
    // POINT DEFINITION
    // scalar-related
        int nmp    ; // number of material point
    // tensor-related    
        DAT* mp    ; // mass
        DAT* vol   ; // volume
        DAT* cohp  ; // initial cohesion
        DAT* phip  ; // friction angle
        DAT* epII  ; // equivalent plastic strain

        DAT* xp    ; // coordinate
        DAT* vp    ; // velocity
        DAT* up    ; // displacement
        DAT* lp    ; // domain lengths

        DAT* sig   ; // cauchy stress tensor
        DAT* eps   ; // infinitesimal strain tensor
        DAT* dev   ; // deviatoric stress tensor
        DAT* dF    ; // increment deformation gradient tensor
        DAT* ome   ; // spin tensor

        int* p2e   ; // point-to-element topology 
        int* p2n   ; // point-to-node topology

        DAT* N     ; // basis function
        DAT* dNx   ; // x-derivative
        DAT* dNy   ; // y-derivative
        DAT* dNz   ; // z-derivative
} point_t;
point_t mpD;
typedef struct mesh{
    // MESH DEFINITION
    // element-related
        DAT* pel   ;
        int* e2n   ;
        DAT  h[3]  ; //[dx,dy,dz]
    // nodes-related
        int  nno[4]; //[nnx ,nny ,nnz ,no]
        int  nel[4]; //[nelx,nely,nelz,no]
        int  nn    ; // -
        DAT  min[3]; //[xnmin,ynmin,znmin]
        DAT* mn    ;
        DAT* xn    ;
        DAT* pn    ;
        DAT* fen   ;
        DAT* fin   ;
        DAT* fn    ;
        DAT* an    ;
        DAT* vn    ;
        DAT* un    ;
        int* bc    ;
} mesh_t;
mesh_t meD;
#include "macros.h"
#include "cpu_functions.h"

int main(){
    init(&meD,&mpD);

    printf("nel = [%d,%d,%d]",meD.nno[0],meD.nno[1],meD.nno[2]);

    // Timer       
    clock_t start, end;
    double CPUinfo[3];
    // Import
    #include "I_O.h"
    // Solver 
    DAT DT   = 0.0;
    int cout_u = 0;
    start = clock();
    while(tw<t){
        if(tw>te && cout_u==0){
            for(int p=0; p<3*mpD.nmp; p++){
                mpD.up[p]=0.0;
            }
            cout_u++;
        }
        // get adaptative dt
        dt = CFL(&meD,&mpD,yd,tg,tw);
        // linear gravitational increase
        g  = getG(tw,tg);
        //
        topol(&meD,&mpD);
        basis(&meD,&mpD);
        accum(&meD,&mpD,g);
        solve(&meD,dt);
        FLIP(&meD,&mpD,dt);
        DM_BC(&meD,&mpD,dt);
        strains(&meD,&mpD,dt);
        elast(&mpD,Del_h,dt);
        if(tw>te){
            DPPlast(&mpD,Hp,cohr,Kc,Gc,psi0);
        }
        volLock(&meD,&mpD);
        // update time & iteration
        tw+=dt;
        it++;
        //printf("\n  workload = %.2f %%",100*tw/t);
        DT+=dt;
        if(((int)FPS>(int)0)&&(DT>=((DAT)1.0/(DAT)FPS))){
            // display workload
            printf("\n  workload = %.2f %%",100*tw/t);
            // reset time interval counter DT
            DT = 0.0;
        }
    } 
    end = clock();
    CPUinfo[0] = (double)(((double)(end-start))/CLOCKS_PER_SEC);
    CPUinfo[1] = it/CPUinfo[0];
    CPUinfo[2] = (IO*it*sizeof(DAT)/(1024*1024*1024*CPUinfo[0]));    
    printf("\n-----------------------------------");
    printf("\nCPU summary: MTPeff = %.2f [GB/s]",CPUinfo[2]);
    printf("\n-----------------------------------");  
    printf("\n  CPU time is %.2f s\n  after %d iterations \n  i.e., %.2f it/s\n",CPUinfo[0],it,CPUinfo[1]);
    // save data
    saveData(mpD.epII,"epII.txt",nmp,1);
    saveData(mpD.lp  ,"lp.txt"  ,nmp,3);
    saveData(mpD.xp  ,"xp.txt"  ,nmp,3);
    saveData(mpD.up  ,"up.txt"  ,nmp,3);
    saveData(mpD.sig ,"sig.txt" ,nmp,6);
}
