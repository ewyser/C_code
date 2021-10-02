// header macros.h
#define PI (DAT)3.1415926535897931
#define IO (DAT)(86.0*nmp+44.0*nmp*nn+26.0*no+2.0*nel)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                              

#define free_all(A)                                         ;\
    free(A##_h)                                             ;
//----------------------------------------------------------// 
#define zeros(A,nx,ny,type)                                 ;\
    type *A##_h                                             ;\
    A##_h = (type*)malloc(nx*ny*sizeof(type))               ;\
//----------------------------------------------------------//
#define load(A,nx,ny,Aname,type)                            ;\
    type *A##_h                                             ;\
    A##_h = (type*)malloc(nx*ny*sizeof(type))               ;\
    FILE* A##fid=fopen(Aname, "rb")                         ;\
    fread(A##_h, sizeof(type), nx*ny, A##fid)               ;\
    fclose(A##fid)                                          ;\
//----------------------------------------------------------//
#define load2struct(strct,A,Aname,dim,type){                ;\
    strct->A = (type*)malloc(dim*sizeof(type))              ;\
    FILE* A##fid=fopen(Aname, "rb")                         ;\
    fread(strct->A, sizeof(type), dim, A##fid)              ;\
    fclose(A##fid)                                          ;\
}                                                           ;
//----------------------------------------------------------//
#define zero4struct(strct,A,dim,type){                      ;\
    strct->A = (type*)malloc(dim*sizeof(type))              ;\
}                                                           ;          
//----------------------------------------------------------//
#define varSave(A,dim1,dim2,name,it){                       ;\
    int it_size = snprintf(NULL,0,"%d",(int)it)             ;\
    int na_size = ((int)sizeof(name)+1)+it_size             ;\
    char filename[na_size+it_size]                          ;\
    sprintf(filename,"%s_%d.txt",name,it)                   ;\
    printf("sizeof(A##) = %d, name = %s",na_size,filename)  ;\
    saveData(A,filename,dim1,dim2)                          ;\
}                                                           ;
//----------------------------------------------------------//


//char greeting[] = "Hello";
//printf("%d",(int)sizeof(greeting));

                                        //                                        char fname[] = name;\
//                                        char str[];\
//                                        sprintf(str,"_%d",(int)it);\
//                                        strcat(fname,str);\