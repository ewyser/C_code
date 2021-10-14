%% EXPORT
%NODES
xn  = [meD.x,meD.y,meD.z];
fid = fopen( 'xn.txt','w'); fprintf(fid,'%f\n',xn(:));fclose(fid);


fid = fopen( 'e2n.txt','w'); fprintf(fid,'%d\n',meD.e2N(:)-1);fclose(fid);

%PHYSICS
fid = fopen( 'cohp.txt','w'); fprintf(fid,'%f\n',mpD.coh(:));fclose(fid);
fid = fopen( 'phip.txt','w'); fprintf(fid,'%f\n',mpD.phi(:));fclose(fid);
%PARAMETERS
xp  = mpD.x;
vol = mpD.V;
lp  = mpD.l;

p   = [mpD.n;meD.nn;meD.no;meD.h(:);min(meD.x);min(meD.y);min(meD.z);meD.nnx;meD.nny;meD.nnz];
fid = fopen( 'param.txt','w'); fprintf(fid,'%f\n',p(:));fclose(fid);
p   = [g;rho0;psi0;nu;E;Kc;Gc;cohr;Hp;t;te;tg];
fid = fopen( 'phys.txt','w'); fprintf(fid,'%f\n',p(:));fclose(fid);
fid = fopen( 'mp.txt','w'); fprintf(fid,'%f\n',mpD.m(:));fclose(fid);
fid = fopen( 'xp.txt','w'); fprintf(fid,'%f\n',xp(:));fclose(fid);
fid = fopen( 'vol.txt','w'); fprintf(fid,'%f\n',vol(:));fclose(fid);
fid = fopen( 'lp.txt','w'); fprintf(fid,'%f\n',lp(:));fclose(fid);

bcx = int32(bc.x);
bcy = int32(bc.y);
bcz = int32(bc.z);
bcx1 = int32(bcx)+0*meD.no;
% bcx2 = int32(bcx)+1*meD.no;
% bcx3 = int32(bcx)+2*meD.no;
% bcy1 = int32(bcy)+0*meD.no;
bcy2 = int32(bcy)+1*meD.no;
% bcy3 = int32(bcy)+2*meD.no;
% bcz1 = int32(bcz)+0*meD.no;
% bcz2 = int32(bcz)+1*meD.no;
bcz3 = int32(bcz)+2*meD.no;
BC  = ones(meD.no*3,1,'int32');
BC([bcx1;bcy2;bcz3]) = 0;
fid = fopen( 'bcs.txt','w'); 
fprintf(fid,'%d\n',BC(:));
fclose(fid);


