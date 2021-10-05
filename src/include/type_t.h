typedef struct point{
    // POINT DEFINITION
    // scalar-related
        int nmp    ; // number of material point
    // tensor-related    
        DAT *mp    ; // mass
        DAT *vol   ; // volume
        DAT *cohp  ; // initial cohesion
        DAT *phip  ; // friction angle
        DAT *epII  ; // equivalent plastic strain

        DAT *xp    ; // coordinate
        DAT *vp    ; // velocity
        DAT *up    ; // displacement
        DAT *lp    ; // domain lengths

        DAT *sig   ; // cauchy stress tensor
        DAT *eps   ; // infinitesimal strain tensor
        DAT *dev   ; // deviatoric stress tensor
        DAT *dF    ; // increment deformation gradient tensor
        DAT *ome   ; // spin tensor

        int *p2e   ; // point-to-element topology 
        int *p2n   ; // point-to-node topology

        DAT *N     ; // basis function
        DAT *dNx   ; // x-derivative
        DAT *dNy   ; // y-derivative
        DAT *dNz   ; // z-derivative
} point_t;

typedef struct mesh{
    // MESH DEFINITION
    // element-related
        DAT *pel   ;
        int *e2n   ;
        DAT  h[3]  ; //[dx,dy,dz]
    // nodes-related
        int  nno[4]; //[nnx ,nny ,nnz ,no]
        int  nel[4]; //[nelx,nely,nelz,no]
        int  nn    ; // -
        DAT  min[3]; //[xnmin,ynmin,znmin]
        DAT *mn    ;
        DAT *xn    ;
        DAT *pn    ;
        DAT *fen   ;
        DAT *fin   ;
        DAT *fn    ;
        DAT *an    ;
        DAT *vn    ;
        DAT *un    ;
        int *bc    ;
} mesh_t;

typedef struct mat{
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
} mat_t;
mat_t phys;