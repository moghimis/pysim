function mat2roms(x,y,h,fout,hmin)
%
% mat2roms(x,y,h,fout,hmin)
%

if(exist('hmin')~=1)
  hmin=-inf;
end

rho.x=x';
rho.y=y';  
rho.depth=h';
if(prod(size(hmin))>1)
  rho.mask=hmin';
else
  rho.mask = ones(size(h'));
  rho.mask(h'<=hmin)=0;
  rho.depth(rho.mask==0)=hmin;
end
spherical='F';
projection='mercator';
save('tmp.mat','rho','spherical','projection')
mat2roms_mw('tmp.mat',fout);
!rm tmp.mat
