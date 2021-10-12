// PHYSICS
    DAT *phys = malloc(12*sizeof(DAT));
         phys = loadFloat(12,"phys.txt");
    DAT g     = (DAT) phys[0 ];
    DAT rho0  = (DAT) phys[1 ];
    DAT psi0  = (DAT) phys[2 ];
    DAT nu    = (DAT) phys[3 ];
    DAT E     = (DAT) phys[4 ];
    DAT Kc    = (DAT) phys[5 ];
    DAT Gc    = (DAT) phys[6 ];
    DAT cohr  = (DAT) phys[7 ];
    DAT Hp    = (DAT) phys[8 ];
    DAT t     = (DAT) phys[9 ];
    DAT te    = (DAT) phys[10];
    DAT tg    = (DAT) phys[11];
    DAT dt    = (DAT) 0.0;
    DAT yd    = (DAT) sqrt((Kc+1.333*Gc)*(DAT)1.0/(DAT)rho0);
    DAT tw    = (DAT) 0.0;
    int it    = (int) 1;
    DAT Del[36];
        Del[0 ] = Kc+4.0/3.0*Gc; 
        Del[1 ] = Kc-2.0/3.0*Gc;
        Del[2 ] = Kc-2.0/3.0*Gc;
        Del[3 ] = 0.0          ;
        Del[4 ] = 0.0          ;
        Del[5 ] = 0.0          ;

        Del[6 ] = Kc-2.0/3.0*Gc;
        Del[7 ] = Kc+4.0/3.0*Gc; 
        Del[8 ] = Kc-2.0/3.0*Gc;
        Del[9 ] = 0.0          ;
        Del[10] = 0.0          ;
        Del[11] = 0.0          ;

        Del[12] = Kc-2.0/3.0*Gc;
        Del[13] = Kc-2.0/3.0*Gc;
        Del[14] = Kc+4.0/3.0*Gc; 
        Del[15] = 0.0          ;
        Del[16] = 0.0          ;
        Del[17] = 0.0          ;

        Del[18] = 0.0          ;
        Del[19] = 0.0          ;
        Del[20] = 0.0          ; 
        Del[21] = Gc           ;
        Del[22] = 0.0          ;
        Del[23] = 0.0          ;

        Del[24] = 0.0          ;
        Del[25] = 0.0          ;
        Del[26] = 0.0          ; 
        Del[27] = 0.0          ;
        Del[28] = Gc           ;
        Del[29] = 0.0          ;

        Del[30] = 0.0          ;
        Del[31] = 0.0          ;
        Del[32] = 0.0          ; 
        Del[33] = 0.0          ;
        Del[34] = 0.0          ;
        Del[35] = Gc           ;
    free(phys);        