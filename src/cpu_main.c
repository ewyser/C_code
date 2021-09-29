#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "time.h"
#include "macros.h"
#include "cpu_functions.h"

int main(){
    // Timer       
    clock_t start, end;
    double CPUinfo[3];
    // Import
    #include "I_O.h"
    // Solver 
    DAT DT   = 0.0;
    int count = 0;
    start = clock();
    while(tw<t){
        if(tw>te && count==0){
            for(int p=0; p<3*nmp; p++){
                up_h[p]=0.0;
            }
            count++;
        }
        // get adaptative dt
        dt = CFL(vp_h,dx,dy,dz,yd,tg,tw,nmp);
        // linear gravitational increase
        g  = getG(tw,tg);
        //
        topol(xp_h,p2e_h,p2n_h,e2n_h,xnmin,ynmin,znmin,dx,dy,dz,nmp,nn,nex,ney,nez,nel);
        basis(xp_h,xn_h,N_h,dNx_h,dNy_h,dNz_h,p2n_h,lp_h,dx,dy,dz,nmp,nn,no);
        accum(mn_h,pn_h,fen_h,fin_h,N_h,dNx_h,dNy_h,dNz_h,mp_h,vp_h,sig_h,vol_h,p2n_h,g,nmp,nn,no);
        solve(fn_h,fen_h,fin_h,mn_h,an_h,pn_h,vn_h,bcs_h,dt,no);
        FLIP(an_h,vn_h,N_h,vp_h,xp_h,p2n_h,dt,nmp,nn,no);
        DM_BC(un_h,pn_h,mn_h,N_h,mp_h,vp_h,up_h,bcs_h,dt,p2n_h,nmp,nn,no);
        strains(un_h,dNx_h,dNy_h,dNz_h,dF_h,eps_h,ome_h,lp_h,vol_h,p2n_h,dt,nmp,nn,no);
        elast(sig_h,eps_h,ome_h,Del_h,nmp,dt);
        if(tw>te){
            DPPlast(sig_h,cohp_h,phip_h,epII_h,Hp,cohr,Kc,Gc,psi0,nmp);
        }
        volLock(pel_h,sig_h,dev_h,vol_h,p2e_h,nmp,nel);
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
    saveData(epII_h,"epII.txt",nmp,1);
    saveData(lp_h  ,"lp.txt"  ,nmp,3);
    saveData(xp_h  ,"xp.txt"  ,nmp,3);
    saveData(up_h  ,"up.txt"  ,nmp,3);
    saveData(sig_h ,"sig.txt" ,nmp,6);

//FILE* fidw=fopen("CPUinfo.dat", "wb");
//fwrite(CPUinfo, 4*sizeof(DAT), 1, fidw);
//fclose(fidw);
}
