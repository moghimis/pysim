function simpleWaves(fin,fout)
%
% simpleWaves(fin,fout)
%
% Very simple wave model, to produce an input file for ROMS.  Assumes the
% bathymetry is alongshore-uniform, and uses TG1986 wave breaking.
% Offshore wave conditions are hard-coded.
%
% fin: input grid file
% fout: output waves file
%

%----------------------------------------------------
% hard-coded offshore wave conditions
%----------------------------------------------------

H0=1;  % wave height in meters
T=7.5;    % wave period
theta0=10;  % offshore wave angle in degrees

%----------------------------------------------------
% read input grid (assume alongshore-uniform, constant grid)
%----------------------------------------------------

x=nc_varget(fin,'x_rho'); x=x(1,:); dx=x(2)-x(1);
h=nc_varget(fin,'h'); h=h(1,:);

%----------------------------------------------------
% compute wave transformation using TG1986
%----------------------------------------------------

% solve dispersion relationship
sigma=2*pi/T;
g=9.8;
k=fsolve(@(k)sigma^2-g*k.*tanh(k.*h),sigma./sqrt(g*h),...
         optimset('Display','off'));

% Snell's Law for wave angle (radians)
c=sigma./k;
theta=asin((c./c(end))*sind(theta0));

% TG1986 for wave shoaling/breaking
gamma=0.45;
B=0.8;
rho=1030;
f=1/T;
TGcoef = 3*sqrt(pi)*rho*g*B^3*f./(16*gamma^4*h.^5);
nx=length(x);
H(nx)=H0;
cg=.5*c.*(1 + 2*k.*h./sinh(2*k.*h));
Ef(length(x)) = (1/8)*rho*g*H0^2*cosd(theta0)*cg(end);
for i=(nx-1):-1:1
  Ef(i) = Ef(i+1) - dx*TGcoef(i+1)*H(i+1)^7;
  H(i) = sqrt(Ef(i)/((1/8)*rho*g*cg(i)*cos(theta(i))));
end

% define some other variables needed by ROMS
Dissip_break=TGcoef.*H.^7*dx;
Uwave_rms=1/sqrt(2)*c.*H./(2*h);
theta=90-rad2deg(theta);  % convert to ROMS angle convention

%----------------------------------------------------
% generate output file (wave input for ROMS)
%----------------------------------------------------

% open file
nc = netcdf(fout,'clobber');
nc.type = ncchar('ROMS FORCING file');
nc.history = ncchar(['Created ' datestr(now)]);

% define grid dimensions
h=nc_varget(fin,'h');
[ny,nx]=size(h);
nc('xi_rho' ) = nx;  
nc('eta_rho') = ny;  
nc('wave_time')=1;

% variable definitions...

nc{'wave_time'} = ncdouble('wave_time');
nc{'wave_time'}.long_name = ncchar('wave field time'); 
nc{'wave_time'}.units = ncchar('Julian day'); 
nc{'wave_time'}.field = ncchar('wave_time, scalar, series'); 

nc{'Dissip_break'} = ncdouble('wave_time','eta_rho','xi_rho');
nc{'Dissip_break'}.long_name = ncchar('Dissip_break'); 
nc{'Dissip_break'}.units = ncchar('Watts meter-2'); 
nc{'Dissip_break'}.field = ncchar('Dissip_break, scalar, series'); 

nc{'Dissip_wcap'} = ncdouble('wave_time','eta_rho','xi_rho');
nc{'Dissip_wcap'}.long_name = ncchar('Dissip_wcap'); 
nc{'Dissip_wcap'}.units = ncchar('Watts meter-2'); 
nc{'Dissip_wcap'}.field = ncchar('Dissip_wcap, scalar, series'); 

nc{'Hwave'} = ncdouble('wave_time','eta_rho','xi_rho');
nc{'Hwave'}.long_name = ncchar('Hwave'); 
nc{'Hwave'}.units = ncchar('meter'); 
nc{'Hwave'}.field = ncchar('Hwave, scalar, series'); 

nc{'Pwave_top'} = ncdouble('wave_time','eta_rho','xi_rho');
nc{'Pwave_top'}.long_name = ncchar('Pwave_top'); 
nc{'Pwave_top'}.units = ncchar('second'); 
nc{'Pwave_top'}.field = ncchar('Pwave_top, scalar, series'); 

nc{'Pwave_bot'} = ncdouble('wave_time','eta_rho','xi_rho');
nc{'Pwave_bot'}.long_name = ncchar('Pwave_bot'); 
nc{'Pwave_bot'}.units = ncchar('second'); 
nc{'Pwave_bot'}.field = ncchar('Pwave_bot, scalar, series'); 

nc{'Uwave_rms'} = ncdouble('wave_time','eta_rho','xi_rho');
nc{'Uwave_rms'}.long_name = ncchar('Uwave_rms'); 
nc{'Uwave_rms'}.units = ncchar('meter second-1'); 
nc{'Uwave_rms'}.field = ncchar('Uwave_rms, scalar, series'); 

nc{'Dwave'} = ncdouble('wave_time','eta_rho','xi_rho');
nc{'Dwave'}.long_name = ncchar('Dwave'); 
nc{'Dwave'}.units = ncchar('degrees'); 
nc{'Dwave'}.field = ncchar('Dwave, scalar, series'); 

nc{'Lwave'} = ncdouble('wave_time','eta_rho','xi_rho');
nc{'Lwave'}.long_name = ncchar('Lwave'); 
nc{'Lwave'}.units = ncchar('meter'); 
nc{'Lwave'}.field = ncchar('Lwave, scalar, series'); 

% assign data
nc{'wave_time'}(:)=0;
nc{'Dissip_break'}(:)=repmat(Dissip_break,ny,1);
nc{'Dissip_wcap' }(:)=0;
nc{'Hwave'       }(:)=repmat(H,ny,1);
nc{'Pwave_top'   }(:)=T;
nc{'Pwave_bot'   }(:)=T;
nc{'Uwave_rms'   }(:)=repmat(Uwave_rms,ny,1);
nc{'Dwave'       }(:)=repmat(theta,ny,1);;
nc{'Lwave'       }(:)=repmat(2*pi./k,ny,1);

% close netCDF file
if(~isempty(close(nc)))
  disp(' ## Unable to close the ROMS output file.')
end

