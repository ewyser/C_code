%% MADMAX
% Copyright (C) 2020: Emmanuel Wyser      , emmanuel.wyser[at]unil.ch
%                     Yury Alkhimenkov    , yury.alkhimenkov[at]unil.ch
%                     Michel Jaboyedoff   , michel.jaboyedoff[at]unil.ch
%                     Yury Y. Podladchikov, yury.podladchikov[at]unil.ch
% -------------------------------------------------------------------------%
% version    : v1.0
% date       : february, 2021
% description: explicit mpm (GIMPM) solver based on an updated Lagrangian
% frame for elasto-plastic problem discretized on a 4-noded quadrilateral 
% mesh
% -------------------------------------------------------------------------%
clear,close,clf                                                           ;%
opengl hardware                                                           ;%
delete('*.dat','*.txt','*.exe','*.avi','*.mat','*.lib','*.exp');
% set precision arithmetic
typeD = 'double';
numel = [80];
for sim=1:length(numel)
    disp('------------------------')                                      ;%
    disp(['Run ',num2str(sim),': nel = ',num2str(numel(sim)),''])         ;%
    disp('------------------------')                                      ;%
    %% NON-DIMENSIONAL CONSTANT
    ar      = 0.8                                                         ;% aspect ratio thickness/width
    nu      = 0.3                                                         ;% poisson ratio
    ni      = 2                                                           ;% number of mp in h(1) direction
    %---------------------------------------------------------------------%
    %% INDEPENDANT PHYSICAL CONSTANT
    g       = 9.81                                                        ;% gravitationnal acceleration [m/s^2]
    E       = 1.0e6                                                       ;% young's modulus             [Pa]
    Gc      = E/(2*(1+nu))                                                ;% shear modulus               [Pa]
    Kc      = E/(3*(1-2*nu))                                              ;% bulk modulus                [Pa]
    rho0    = 2700                                                        ;% density                     [kg/m^3]
    n0      = 0.0                                                         ;% initial porosity            [-]
    yd      = sqrt((Kc+4/3*Gc)/rho0)                                      ;% elastic wave speed          [m/s]
    coh0    = 20.0e3                                                      ;% cohesion                    [Pa]
    phi0    = 20.0*pi/180                                                 ;% friction angle              [Rad]
    psi0    =  0.0*pi/180                                                 ;%
    H       = -60e3                                                       ;% softening modulus           [Pa]
    cohr    =  4.0e3                                                      ;% residual cohesion           [Pa]
    phir    = 7.0*pi/180                                                  ;% residual friction angle     [Rad]
    t       = 15.0                                                        ;% simulation time             [s]
    te      = 8.0                                                         ;% elastic loading             [s]
    tg      = te/1.5                                                      ;% gravity increase            [s]
    %---------------------------------------------------------------------%
    %% MESH & MP INITIALIZATION
    lx      = 64.1584                                                     ;% mesh size along x direction
    ly      = lx/4                                                        ;% mesh size along y direction
    lz      = 12.80                                                       ;% mesh size along z direction
    [meD,bc]= meSetup(numel(sim),lx,ly,lz,typeD)                          ;% - see function
    [mpD]   = mpSetup(meD,ni,lz,coh0,cohr,phi0,phir,n0,rho0,typeD)        ;% - see function
    % ISOTROPIC ELASTIC MATRIX
    Del     = [ Kc+4/3*Gc,Kc-2/3*Gc,Kc-2/3*Gc,0.0,0.0,0.0;...
                Kc-2/3*Gc,Kc+4/3*Gc,Kc-2/3*Gc,0.0,0.0,0.0;...
                Kc-2/3*Gc,Kc-2/3*Gc,Kc+4/3*Gc,0.0,0.0,0.0;...
                0.0      ,0.0      ,0.0      ,Gc ,0.0,0.0;...                             
                0.0      ,0.0      ,0.0      ,0.0,Gc ,0.0;...
                0.0      ,0.0      ,0.0      ,0.0,0.0,Gc ;]               ;%     
    Hp      = H*meD.h(1)                                                  ;%
    %---------------------------------------------------------------------%
    %% MPM DM ALGORITHM EXPLICIT SOLVER FOR INFINITESIMAL STRAIN
    run('export');  
end


% writematrix([x p0],'myData.txt','Delimiter',';'); 
% type myData.txt;
