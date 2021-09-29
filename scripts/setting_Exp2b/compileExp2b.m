%% MATLAB SCRIPT TO RUN COMPILATION FOR GPU AND CPU EXECUTABLES (.EXE)
% delete what needs to be deleted
delete('*.dat','*.txt','*.out','*.exe','*.avi','*.mat','*.lib','*.exp');
% set precision arithmetic
typeD = 'double';
% compile cpu code
if(ismac || isunix)
    call = ['clang ../../src/cpu_main.c -o cpu.out -O3 '];
elseif(ispc)
    call = ['clang ../../src/cpu_main.c -o cpu.exe -O3 '];
else
    disp('Outer space OS');
end
if(typeD=='single')
    arg  = ['-DDAT=float  -DSIM=1 -DFPS=5 -DD=0.1'];
    system([call arg]); 
elseif(typeD=='double')
    arg  = ['-DDAT=double -DSIM=1 -DFPS=5 -DD=0.1'];
    system([call arg]);
end