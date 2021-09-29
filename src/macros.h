// header macros.h
#define PI (DAT)3.1415926535897931
#define IO (DAT)(86.0*nmp+44.0*nmp*nn+26.0*no+2.0*nel)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define zeros(A,nx,ny,type)      type *A##_h                                                                          ;\
                                 A##_h = (type*)malloc(nx*ny*sizeof(type))                                            ;
                                 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define load(A,nx,ny,Aname,type) type *A##_h                                                                          ;\
                                 A##_h = (type*)malloc(nx*ny*sizeof(type))                                            ;\
                                 FILE* A##fid=fopen(Aname, "rb")                                                      ;\
                                 fread(A##_h, sizeof(type), nx*ny, A##fid)                                            ;\
                                 fclose(A##fid)                                                                       ;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define save(A,nx,ny,Aname)      FILE* A##fidw=fopen(Aname, "wb")                                                     ;\
                                 fwrite(A##_h, sizeof(DAT), ((nx)*(ny)), A##fidw)                                     ;\
                                 fclose(A##fidw)                                                                      ;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define free_all(A)              free(A##_h);                                                                         ;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                              
#define D2HD(A,nx,ny,name,type,it,sim) char str##A[100],sim##A[100],fname##A[100] = name                              ;\
                                       sprintf(str##A,"_%d",it)                                                       ;\
                                       sprintf(sim##A,"_%d",sim)                                                      ;\
                                       strcat(fname##A,str##A)                                                        ;\
                                       strcat(fname##A,sim##A)                                                        ;\
                                       strcat(fname##A,".txt")                                                        ;\
                                       D2H (A,nx,ny,type)                                                             ;\
                                       save(A,nx,ny,fname##A)                                                         ;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  