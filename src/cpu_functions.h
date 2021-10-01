// function declaration
void init(mesh_t* meD,point_t* mpD);
void saveData(DAT* data, char* filename, int dim1, int dim2);
DAT CFL(mesh_t* meD, point_t* mpD, DAT yd, DAT tg, DAT tw);
DAT getG(DAT tw, DAT tg);
void topol(mesh_t* meD,point_t* mpD);
void NdN(DAT* N_dN,DAT xi,DAT lp,DAT dx);
void basis(mesh_t* meD, point_t* mpD);
void accum(mesh_t* meD,point_t* mpD,DAT g);
void solve(mesh_t* meD, DAT dt);
void FLIP(mesh_t* meD, point_t* mpD,DAT dt);
void DM_BC(mesh_t* meD, point_t* mpD, DAT dt);
void strains(mesh_t* meD, point_t* mpD,DAT dt);
void elast(point_t* mpD, DAT* Del, DAT dt);
void DPPlast(point_t* mpD, DAT Hp, DAT cohr, DAT Kc, DAT Gc, DAT psi);
void volLock(mesh_t* meD, point_t* mpD);
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
// function definition
void init(mesh_t* meD,point_t* mpD){
    // mesh & material point geometric properties
        DAT  *par = (DAT*)malloc(12*sizeof(DAT));
        FILE *fid = fopen("param.dat", "rb");                                                      ;\
                    fread(par, sizeof(DAT),12,fid); 

        mpD->nmp   = (int) par[0 ];
        meD->nn    = (int) par[1 ];
        meD->nno[3]= (int) par[2 ];
        meD->h[0]  = (DAT) par[3 ];
        meD->h[1]  = (DAT) par[4 ];
        meD->h[2]  = (DAT) par[5 ];
        meD->min[0]= (DAT) par[6 ];
        meD->min[1]= (DAT) par[7 ];
        meD->min[2]= (DAT) par[8 ];
        meD->nno[0]= (int) par[9 ];
        meD->nno[1]= (int) par[10];
        meD->nno[2]= (int) par[11];
        meD->nno[3]= meD->nno[0]*meD->nno[1]*meD->nno[2];
        meD->nel[0]= meD->nno[0]-1;
        meD->nel[1]= meD->nno[1]-1;
        meD->nel[2]= meD->nno[2]-1;
        meD->nel[3]= meD->nel[0]*meD->nel[1]*meD->nel[2];  
        free(par);
        fclose(fid);
    // mesh
        zero4struct(meD,pel            ,2*meD->nno[3]      ,DAT);    
        load2struct(meD,e2n ,"e2n.dat" ,meD->nn*meD->nel[3],int);
        zero4struct(meD,mn             ,  meD->nno[3]      ,DAT);
        load2struct(meD,xn  ,"xn.dat"  ,3*meD->nno[3]      ,DAT);
        zero4struct(meD,pn             ,3*meD->nno[3]      ,DAT);
        zero4struct(meD,pn             ,3*meD->nno[3]      ,DAT);
        zero4struct(meD,fen            ,3*meD->nno[3]      ,DAT);
        zero4struct(meD,fin            ,3*meD->nno[3]      ,DAT);
        zero4struct(meD,fn             ,3*meD->nno[3]      ,DAT);
        zero4struct(meD,an             ,3*meD->nno[3]      ,DAT);
        zero4struct(meD,vn             ,3*meD->nno[3]      ,DAT);
        zero4struct(meD,un             ,3*meD->nno[3]      ,DAT);
        load2struct(meD,bc  ,"bcs.dat" ,3*meD->nno[3]      ,int);
    // material point
        load2struct(mpD,mp  ,"mp.dat"  ,mpD->nmp           ,DAT);
        load2struct(mpD,vol ,"vol.dat" ,mpD->nmp           ,DAT);
        load2struct(mpD,cohp,"cohp.dat",mpD->nmp           ,DAT);
        load2struct(mpD,phip,"phip.dat",mpD->nmp           ,DAT);
        zero4struct(mpD,epII           ,mpD->nmp           ,DAT);

        load2struct(mpD,xp  ,"xp.dat"  ,3*mpD->nmp         ,DAT);
        zero4struct(mpD,vp             ,3*mpD->nmp         ,DAT);    
        zero4struct(mpD,up             ,3*mpD->nmp         ,DAT);    
        load2struct(mpD,lp  ,"lp.dat"  ,3*mpD->nmp         ,DAT);   

        zero4struct(mpD,sig            ,6*mpD->nmp         ,DAT);
        zero4struct(mpD,eps            ,6*mpD->nmp         ,DAT);
        zero4struct(mpD,dev            ,6*mpD->nmp         ,DAT);
        zero4struct(mpD,dF             ,9*mpD->nmp,         DAT);
        zero4struct(mpD,ome            ,3*mpD->nmp         ,DAT);

        zero4struct(mpD,p2e            ,        mpD->nmp   ,int);
        zero4struct(mpD,p2n            ,meD->nn*mpD->nmp   ,int);   

        zero4struct(mpD,N              ,meD->nn*mpD->nmp   ,DAT);
        zero4struct(mpD,dNx            ,meD->nn*mpD->nmp   ,DAT);
        zero4struct(mpD,dNy            ,meD->nn*mpD->nmp   ,DAT);
        zero4struct(mpD,dNz            ,meD->nn*mpD->nmp   ,DAT);         
    /*
    for(int k=0; k<nmp; k++){
        printf("\n xp = [%f,%f,%f]",mpD->xp[k+0*nmp],mpD->xp[k+1*mpD->nmp],mpD->xp[k+2*mpD->nmp]);
    } 
    */                                                                  
}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
void saveData(DAT* data, char* filename, int dim1, int dim2){
     FILE *fp=fopen(filename,"w");
     
     if(dim2==1){
         for(int n=0;n<dim1;n++) {
             fprintf(fp,"%f\n",data[n]);
         }
     }
     else if(dim2==3){
         for(int n=0;n<dim1;n++) {
             fprintf(fp,"%f,%f,%f\n",data[n+0*dim1],data[n+1*dim1],data[n+2*dim1]);
         }
     }
     else if(dim2==6){
         for(int n=0;n<dim1;n++) {
             fprintf(fp,"%f,%f,%f,%f,%f,%f\n",data[0+dim2*n],data[1+dim2*n],data[2+dim2*n],data[3+dim2*n],data[4+dim2*n],data[5+dim2*n]);
         }        
     }
     fclose(fp);
}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
DAT CFL(mesh_t* meD, point_t* mpD, DAT yd, DAT tg, DAT tw){
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
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
DAT getG(DAT tw, DAT tg){
    DAT g = 0.0;
    if(tw<=tg){
        g = 9.81*tw*((DAT)1.0/(DAT)tg);
    }
    else{
        g = 9.81;
    }
    return g;
}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
void topol(mesh_t* meD,point_t* mpD){
    DAT dx    = meD->h[0];
    DAT dy    = meD->h[1];
    DAT dz    = meD->h[2];
    DAT xnmin = meD->min[0];
    DAT ynmin = meD->min[1];
    DAT znmin = meD->min[2];
    int nex   = meD->nel[0];
    int nez   = meD->nel[2];
    int nel   = meD->nel[3];
    int nn    = meD->nn;
    int nmp   = mpD->nmp;
    int id;    
    // caching
    DAT one_dx = ((DAT)1.0/(DAT)dx);
    DAT one_dy = ((DAT)1.0/(DAT)dy);
    DAT one_dz = ((DAT)1.0/(DAT)dz);
    // mp-to-element-to-node topology
    for(int i=0;i<nmp;i++){
        mpD->p2e[i] = (int)((floor((mpD->xp[i+2*nmp]-znmin)*one_dz))+(nez)*floor((mpD->xp[i+0*nmp]-xnmin)*one_dx)+nez*nex*(floor((mpD->xp[i+1*nmp]-ynmin)*one_dy)));
        id          = mpD->p2e[i];
        for(int j=0;j<nn;j++){
            mpD->p2n[i+j*nmp] = meD->e2n[id+j*nel];
        }
    }
}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
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
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
void basis(mesh_t* meD, point_t* mpD){
    DAT Nx_dNx[2],Ny_dNy[2],Nz_dNz[2];
    int iD;
    int nmp = mpD->nmp, 
         nn = meD->nn, 
         no = meD->nno[3];
    for(int n=0;n<nn;n++){
        for(int p=0;p<nmp;p++){
            iD  = p+n*nmp;       
            // x-basis and x-derivative
            NdN(Nx_dNx,(mpD->xp[p+0*nmp]-meD->xn[mpD->p2n[iD]+0*no]),mpD->lp[p+0*nmp],meD->h[0]);
            // y-basis and y-derivative   
            NdN(Ny_dNy,(mpD->xp[p+1*nmp]-meD->xn[mpD->p2n[iD]+1*no]),mpD->lp[p+1*nmp],meD->h[1]);
            // y-basis and y-derivative   
            NdN(Nz_dNz,(mpD->xp[p+2*nmp]-meD->xn[mpD->p2n[iD]+2*no]),mpD->lp[p+2*nmp],meD->h[2]);
            // convolution of basis
            mpD->N[iD]   = Nx_dNx[0]*Ny_dNy[0]*Nz_dNz[0];
            mpD->dNx[iD] = Nx_dNx[1]*Ny_dNy[0]*Nz_dNz[0];
            mpD->dNy[iD] = Nx_dNx[0]*Ny_dNy[1]*Nz_dNz[0];
            mpD->dNz[iD] = Nx_dNx[0]*Ny_dNy[0]*Nz_dNz[1];
        }
    }
}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
void accum(mesh_t* meD,point_t* mpD,DAT g){
    int id,iD, iDx, iDy, iDz;
    int nmp = mpD->nmp, 
         nn = meD->nn, 
         no = meD->nno[3];
    DAT cache;
    // initialize
    for(int n=0;n<(3*no);n++){
        if(n<no){
            meD->mn[n] = 0.0;
        }
        meD->pn [n] = meD->fen[n] = meD->fin[n] = 0.0;
    }
    // accumulate material point contributions on nodes
    for(int n=0;n<nn; n++){
        for(int p=0;p<nmp;p++){
            // indexing & caching
            id        = p+n*nmp;
            iD        = (mpD->p2n[id]);
            iDx       = iD+0*no;
            iDy       = iD+1*no;
            iDz       = iD+2*no;
            cache     = mpD->N[p+n*nmp]*mpD->mp[p];
            // accumulation
            meD->mn [iD ] += (cache);
            meD->pn [iDx] += (cache*mpD->vp[p+0*nmp]);
            meD->pn [iDy] += (cache*mpD->vp[p+1*nmp]);
            meD->pn [iDz] += (cache*mpD->vp[p+2*nmp]);
            meD->fen[iDz] -= (cache*g          );
            meD->fin[iDx] += mpD->vol[p]*(mpD->dNx[id]*mpD->sig[0+6*p]+mpD->dNy[id]*mpD->sig[3+6*p]+mpD->dNz[id]*mpD->sig[5+6*p]);
            meD->fin[iDy] += mpD->vol[p]*(mpD->dNx[id]*mpD->sig[3+6*p]+mpD->dNy[id]*mpD->sig[1+6*p]+mpD->dNz[id]*mpD->sig[4+6*p]);
            meD->fin[iDz] += mpD->vol[p]*(mpD->dNx[id]*mpD->sig[5+6*p]+mpD->dNy[id]*mpD->sig[4+6*p]+mpD->dNz[id]*mpD->sig[2+6*p]);
        }
    }
}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
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
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
void FLIP(mesh_t* meD, point_t* mpD,DAT dt){
    // update material point velocity
    int id,iDx,iDy,iDz;
    int nmp = mpD->nmp, 
         nn = meD->nn, 
         no = meD->nno[3];
    DAT dvpx,dvpy,dvpz,dxp,dyp,dzp;
    for(int p=0;p<nmp;p++){
        dvpx = dvpy = dvpz = dxp = dyp = dzp = (DAT)0.0;
        for(int n=0;n<nn;n++){
            id    = p+n*nmp        ;
            iDx   = mpD->p2n[id]+0*no   ;
            iDy   = mpD->p2n[id]+1*no   ;
            iDz   = mpD->p2n[id]+2*no   ;
            dvpx += (mpD->N[id]*meD->an[iDx]);
            dvpy += (mpD->N[id]*meD->an[iDy]);
            dvpz += (mpD->N[id]*meD->an[iDz]);
            dxp  += (mpD->N[id]*meD->vn[iDx]);
            dyp  += (mpD->N[id]*meD->vn[iDy]);
            dzp  += (mpD->N[id]*meD->vn[iDz]);
        }
        mpD->vp[p+0*nmp] += (dt*dvpx);
        mpD->vp[p+1*nmp] += (dt*dvpy);
        mpD->vp[p+2*nmp] += (dt*dvpz);
        mpD->xp[p+0*nmp] += (dt*dxp );
        mpD->xp[p+1*nmp] += (dt*dyp );
        mpD->xp[p+2*nmp] += (dt*dzp );
    }   
}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
void DM_BC(mesh_t* meD, point_t* mpD, DAT dt){
    int id,iDx, iDy, iDz;
    int nmp = mpD->nmp, 
         nn = meD->nn, 
         no = meD->nno[3];
    DAT cache,m,duxp,duyp,duzp;
    // initialize nodal momentum
    for(int n=0;n<(3*no);n++){
        meD->pn[n] = meD->un[n] = (DAT)0.0;
    }
    // accumulate material point momentum
    for(int n=0;n<nn;n++){
        for(int p=0;p<nmp;p++){
            // indexing & caching
            id            = p+n*nmp                 ;
            iDx           = mpD->p2n[id]+0*no       ;
            iDy           = mpD->p2n[id]+1*no       ;
            iDz           = mpD->p2n[id]+2*no       ;
            cache         = mpD->N[id]*mpD->mp[p]   ;
            // accumulation
            meD->pn[iDx] += (cache*mpD->vp[p+0*nmp]);
            meD->pn[iDy] += (cache*mpD->vp[p+1*nmp]);
            meD->pn[iDz] += (cache*mpD->vp[p+2*nmp]);
        }
    }
    // solve for nodal incremental displacement
    for(int n=0;n<no;n++){
        iDx = n+0*no ;
        iDy = n+1*no ;
        iDz = n+2*no ;
        if(meD->mn[n]>(DAT)0.0){
            m            = ((DAT)1.0/(DAT)meD->mn[n])         ;
            meD->un[iDx] = dt*meD->pn[iDx]*m*(DAT)meD->bc[iDx];
            meD->un[iDy] = dt*meD->pn[iDy]*m*(DAT)meD->bc[iDy];
            meD->un[iDz] = dt*meD->pn[iDz]*m*(DAT)meD->bc[iDz];
        }
    }
    // update material point displacement
    for(int p=0;p<nmp;p++){
        duxp = duyp = duzp = (DAT)0.0;
        for(int n=0;n<nn;n++){
            iDx   = mpD->p2n[p+n*nmp]+0*no        ;
            iDy   = mpD->p2n[p+n*nmp]+1*no        ;
            iDz   = mpD->p2n[p+n*nmp]+2*no        ;
            duxp += (mpD->N[p+n*nmp]*meD->un[iDx]);
            duyp += (mpD->N[p+n*nmp]*meD->un[iDy]);
            duzp += (mpD->N[p+n*nmp]*meD->un[iDz]);
        }
        mpD->up[p+0*nmp] += (duxp);
        mpD->up[p+1*nmp] += (duyp);
        mpD->up[p+2*nmp] += (duzp);
    } 
}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
void strains(mesh_t* meD, point_t* mpD,DAT dt){
    int id,iDx,iDy,iDz;
    int nmp = mpD->nmp, 
         nn = meD->nn, 
         no = meD->nno[3];
    DAT J;
    // update material point state variables
    for(int p=0;p<nmp;p++){
        // compute incremental deformation gradient
        mpD->dF[0+p*9] = mpD->dF[4+p*9] = mpD->dF[8+p*9] = (DAT)1.0;
        mpD->dF[1+p*9] = mpD->dF[2+p*9] = mpD->dF[3+p*9] = mpD->dF[5+p*9] = mpD->dF[6+p*9] = mpD->dF[7+p*9] = (DAT)0.0;
        for(int n=0;n<nn;n++){
            id              = p+n*nmp          ;
            iDx             = mpD->p2n[id]+0*no     ;
            iDy             = mpD->p2n[id]+1*no     ;
            iDz             = mpD->p2n[id]+2*no     ;
            mpD->dF[0+p*9] += (mpD->dNx[id]*meD->un[iDx]);
            mpD->dF[1+p*9] += (mpD->dNy[id]*meD->un[iDx]);
            mpD->dF[2+p*9] += (mpD->dNz[id]*meD->un[iDx]);
            
            mpD->dF[3+p*9] += (mpD->dNx[id]*meD->un[iDy]);
            mpD->dF[4+p*9] += (mpD->dNy[id]*meD->un[iDy]);
            mpD->dF[5+p*9] += (mpD->dNz[id]*meD->un[iDy]);
            
            mpD->dF[6+p*9] += (mpD->dNx[id]*meD->un[iDz]);
            mpD->dF[7+p*9] += (mpD->dNy[id]*meD->un[iDz]);
            mpD->dF[8+p*9] += (mpD->dNz[id]*meD->un[iDz]);
        }
        // compute incremental strains
        mpD->eps[0+p*6] = mpD->dF[0+p*9]-(DAT)1.0;
        mpD->eps[1+p*6] = mpD->dF[4+p*9]-(DAT)1.0;
        mpD->eps[2+p*6] = mpD->dF[8+p*9]-(DAT)1.0;
        mpD->eps[3+p*6] = mpD->dF[1+p*9]+mpD->dF[3+p*9];
        mpD->eps[4+p*6] = mpD->dF[7+p*9]+mpD->dF[5+p*9];
        mpD->eps[5+p*6] = mpD->dF[6+p*9]+mpD->dF[2+p*9];
        // compute incremental rotation
        mpD->ome[0+p*3] = (DAT)0.5*(mpD->dF[3+p*9]-mpD->dF[1+p*9]);
        mpD->ome[1+p*3] = (DAT)0.5*(mpD->dF[5+p*9]-mpD->dF[7+p*9]);
        mpD->ome[2+p*3] = (DAT)0.5*(mpD->dF[2+p*9]-mpD->dF[6+p*9]);
        // update material point volume and domain lengths
        J          = (DAT)1.0+(mpD->eps[0+p*6]+mpD->eps[1+p*6]+mpD->eps[2+p*6]);
        mpD->vol[p]= J*mpD->vol[p];
        if(sizeof(DAT)==(int)8){
            J           = pow(J,(DAT)0.3333);
        }
        else if(sizeof(DAT)==(int)4){
            J           = powf(J,(DAT)0.3333);
        }
        mpD->lp[p+0*nmp] = J*mpD->lp[p+0*nmp];
        mpD->lp[p+1*nmp] = J*mpD->lp[p+1*nmp];
        mpD->lp[p+2*nmp] = J*mpD->lp[p+2*nmp];
    }
}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
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
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
void DPPlast(point_t* mpD, DAT Hp, DAT cohr, DAT Kc, DAT Gc, DAT psi){
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
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
void volLock(mesh_t* meD, point_t* mpD){
    int iD;
    int nmp = mpD->nmp,
        nel = meD->nel[3];
    DAT pr;
    // initialize element pressure
    for(unsigned int e=0;e<nel;e++){
        meD->pel[0+2*e] = (DAT)0.0;
        meD->pel[1+2*e] = (DAT)0.0;
    }
    // accumulate material point pressure on element
    for(unsigned int p=0;p<nmp;p++){
        iD          = mpD->p2e[p]                 ;
        pr          = -(mpD->sig[0+6*p]+mpD->sig[1+6*p]+mpD->sig[2+6*p])/(DAT)3.0;
        mpD->dev[0+6*p ] = mpD->sig[0+6*p]+pr               ;
        mpD->dev[1+6*p ] = mpD->sig[1+6*p]+pr               ;
        mpD->dev[2+6*p ] = mpD->sig[2+6*p]+pr               ;
        mpD->dev[3+6*p ] = mpD->sig[3+6*p]+(DAT)0.0         ;
        mpD->dev[4+6*p ] = mpD->sig[4+6*p]+(DAT)0.0         ;
        mpD->dev[5+6*p ] = mpD->sig[5+6*p]+(DAT)0.0         ;
        meD->pel[0+2*iD]+= pr*mpD->vol[p]                   ;
        meD->pel[1+2*iD]+= mpD->vol[p]                      ;
    }
    // assign element pressure to material point
    for(unsigned int p=0;p<nmp;p++){
        iD          = mpD->p2e[p]                 ;
        pr          = meD->pel[0+2*iD]*1.0/meD->pel[1+2*iD] ;
        mpD->sig[0+6*p]  = mpD->dev[0+6*p]-pr               ;
        mpD->sig[1+6*p]  = mpD->dev[1+6*p]-pr               ;
        mpD->sig[2+6*p]  = mpD->dev[2+6*p]-pr               ;
        mpD->sig[3+6*p]  = mpD->dev[3+6*p]                  ;
        mpD->sig[4+6*p]  = mpD->dev[4+6*p]                  ;
        mpD->sig[5+6*p]  = mpD->dev[5+6*p]                  ;
    }
}
//-----------------------------------------------------------------------//