#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "time.h"

#include "include/type_t.h"
#include "include/macros.h"
#include "include/functions.h"

int main(){
    printf("\no---------------------------------------------o");
    printf("\n|              ** ep23De v1.0 **              |");
    printf("\no---------------------------------------------o");
    #include "include/io.h"
    point_t *mpD=(point_t *)malloc(sizeof(point_t));
    mesh_t  *meD=(mesh_t  *)malloc(sizeof(point_t));
    init(meD,mpD);
    int nmp = mpD->nmp, 
         nn = meD->nn,
         nel= meD->nel[3], 
         no = meD->nno[3];   
    printf("\n nel = [%d,%d,%d]",meD->nel[0],meD->nel[1],meD->nel[2]);
    printf("\n nno = [%d,%d,%d]",meD->nno[0],meD->nno[1],meD->nno[2]);
    printf("\n nmp = %d",mpD->nmp);
    printf("\no---------------------------------------------o"); 
    // Timer       
    clock_t start, end;
    double CPUinfo[3];
    // Solver 
    int flag_u = 0;
    start = clock();
    while(tw<t){
        if(tw>te && flag_u==0){
            for(int p=0; p<3*mpD->nmp; p++){
                mpD->up[p]=0.0;
            }
            flag_u++;
        }
        // get adaptative dt
        dt = CFL(meD,mpD,yd,tg,tw);
        // linear gravitational increase
        g  = getG(tw,tg);
        //
        topol(meD,mpD);
        basis(meD,mpD);
        accum(meD,mpD,g);
        solve(meD,mpD,dt);
        FLIP(meD,mpD,dt);
        DM_BC(meD,mpD,dt);
        strains(meD,mpD,dt);
        elast(mpD,Del,dt);
        if(tw>te){
            DPPlast(mpD,Hp,cohr,Kc,Gc,psi0);
        }
        volLock(meD,mpD);
        // update time & iteration
        tw+=dt;
        it++;
        if((it%25)==0){
            // display workload
            printf("\n (%.2f %%);\tit = %d [-],\tt = %.2f [s]",100*tw/t,it, tw);
        }
    } 
    
    end = clock();
    CPUinfo[0] = (double)(((double)(end-start))/CLOCKS_PER_SEC);
    CPUinfo[1] = it/CPUinfo[0];
    CPUinfo[2] = (IO*it*sizeof(DAT)/(1024*1024*1024*CPUinfo[0]));
    printf("\no---------------------------------------------o");
    printf("\n Performance: %.2f [GB/s] and %.2f [it/s]",CPUinfo[2],CPUinfo[1]);
    printf("\n CPU time   : %.2f [s]                   ",CPUinfo[0]           );
    printf("\n n_it       : %d [-]                     ",it                   );
    // save data
    saveData(mpD->epII,"epII.txt",nmp,1);
    saveData(mpD->lp  ,"lp.txt"  ,nmp,3);
    saveData(mpD->xp  ,"xp.txt"  ,nmp,3);
    saveData(mpD->up  ,"up.txt"  ,nmp,3);
    saveData(mpD->sig ,"sig.txt" ,nmp,6);

    //varSave(mpD->lp,nmp,3,"domain_length",it);
    //varSave(mpD->xp,nmp,3,"xp",it);
    // free & clear
        free(meD->pel) ;   
        free(meD->e2n) ;   
        free(meD->mn)  ;    
        free(meD->xn)  ;    
        free(meD->pn)  ;    
        free(meD->fen) ;   
        free(meD->fin) ;  
        free(meD->fn)  ;   
        free(meD->an)  ;    
        free(meD->vn)  ;    
        free(meD->un)  ;    
        free(meD->bc)  ;    
    free(meD)          ;
        free(mpD->mp)  ;
        free(mpD->vol) ;    
        free(mpD->cohp);   
        free(mpD->phip);  
        free(mpD->epII);  
        free(mpD->xp)  ;  
        free(mpD->vp)  ;   
        free(mpD->up)  ;  
        free(mpD->lp)  ;   
        free(mpD->sig) ;  
        free(mpD->eps) ;   
        free(mpD->dev) ;  
        free(mpD->dF)  ;    
        free(mpD->ome) ;  
        free(mpD->p2e) ;  
        free(mpD->p2n) ;  
        free(mpD->N)   ;  
        free(mpD->dNx) ;
        free(mpD->dNy) ;
        free(mpD->dNz) ;
    free(mpD)          ;
}
