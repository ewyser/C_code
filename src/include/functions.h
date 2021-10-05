// function prototype
void init(mesh_t* meD,point_t* mpD);
void saveData(DAT* data, char* filename, int dim1, int dim2);
DAT CFL(mesh_t* meD, point_t* mpD, DAT yd, DAT tg, DAT tw);
DAT getG(DAT tw, DAT tg);
void topol(mesh_t* meD,point_t* mpD);
void NdN(DAT* N_dN,DAT xi,DAT lp,DAT dx);
void basis(mesh_t* meD, point_t* mpD);
void accum(mesh_t* meD, point_t* mpD,DAT g);
void solve(mesh_t* meD, point_t* mpD,DAT dt);
void FLIP(mesh_t* meD, point_t* mpD,DAT dt);
void DM_BC(mesh_t* meD, point_t* mpD, DAT dt);
void strains(mesh_t* meD, point_t* mpD,DAT dt);
void elast(point_t* mpD, DAT* Del, DAT dt);
void DPPlast(point_t* mpD, DAT Hp, DAT cohr, DAT Kc, DAT Gc, DAT psi);
void volLock(mesh_t* meD, point_t* mpD);