#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "time.h"


#include "struct_mesh_mp_mat.h"
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

    varSave(mpD.lp,nmp,3,"domain_length",it);
    varSave(mpD.xp,nmp,3,"xp",it);
}
