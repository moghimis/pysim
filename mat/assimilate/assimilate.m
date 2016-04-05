function out=assimilate(membersDir,outdir)
%
% out=assim_iter(membersDir,outdir)
%
% Executes assimilation for standardized iterative assimilation directory
% structure
%
% i.e. member <n>'s model output can be found in:
%        membersDir/member<n>/currents/ocean_avg.nc
%
% The posterior bathymetry files are stored in netcdf format in
%        <outdir>/netcdf/bath<n>.nc
%
% If outdir=='', then no output is written
%
warning off
addpath ~/work/coawstRoms/myCOAWST/tools
warning on

matlabpool

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

if(~exist('outdir')) outdir=''; end
if(membersDir(end)~='/') membersDir(end+1)='/'; end
if(length(outdir)>0 & outdir(end) ~='/') outdir(end+1) ='/'; end

%-----------------------------------------------------
% load prior ensemble
%-----------------------------------------------------

% read roms output files for each member
disp('loading prior ensemble')
snap=@(data)squeeze(data(end,:,:));
f=dir([membersDir 'member*']);
N=length(f);
failind=[];
for i=1:N
  if(floor(i/N*10)>floor((i-1)/N*10))
    disp([num2str(floor(i/N*10)*10) '%'])
  end
  memberID(i)=str2num(f(i).name(7:end));
  try
    fn=[membersDir f(i).name '/currents/' subdir 'ocean_his.nc'];
    D(:,:,i)=snap(nc_varget(fn,'Dwave'));
    H(:,:,i)=snap(nc_varget(fn,'Hwave'));
    m(:,:,i)=nc_varget(fn,'mask_rho');
    fn=[membersDir f(i).name '/currents/' subdir 'ocean_avg.nc'];
    u(:,:,i)=uv2rho(snap(nc_varget(fn,'ubar')+nc_varget(fn,'ubar_stokes')),'u');
    v(:,:,i)=uv2rho(snap(nc_varget(fn,'vbar')+nc_varget(fn,'vbar_stokes')),'v');
    us(:,:,i)=uv2rho(snap(nc_varget(fn,'ubar_stokes')),'u');
    vs(:,:,i)=uv2rho(snap(nc_varget(fn,'vbar_stokes')),'v');
    z(:,:,i)=snap(nc_varget(fn,'zeta'));
    h(:,:,i)=nc_varget(fn,'h');
catch
  failind=[failind i];
end

% in case some members failed for whatever reason: discard them
if(length(failind)>0)
  disp(['WARNING WARNING WARNING: the following members failed:'])
  for i=1:length(failind)
    disp(['member' num2str(memberID(failind(i)))])
  end
  keepind=setdiff(1:length(f),failind);
  u=u(:,:,keepind);
  v=v(:,:,keepind);
  us=us(:,:,keepind);
  vs=vs(:,:,keepind);
  D=D(:,:,keepind);
  H=H(:,:,keepind);
  z=z(:,:,keepind);
  m=m(:,:,keepind);
  h=h(:,:,keepind);
  f=f(keepind);
  memberID=memberID(keepind);
end

% load the model grid (is the same for all members)
fn=[membersDir f(1).name '/currents/' subdir 'ocean_his.nc'];
x=nc_varget(fn,'x_rho'); x=x(1,:)';
y=nc_varget(fn,'y_rho'); y=y(:,1);

% reject members which contain islands
bcdata=load([membersDir '../../../bcData.mat']);
zz=z; zz(m==0)=0;
htot=h+bcdata.tide+zz;
htot(htot<0.1)=nan;
[ny,nx,N]=size(htot);
keepind=ones(1,N);
for n=1:N
  for j=1:ny
    if(min(find(~isnan(htot(j,:,n))))<max(find(isnan(htot(j,:,n)))))
      keepind(n)=0;
      break;
    end
  end
end
disp(['throwing out ' num2str(sum(1-keepind)) ' members due to islands'])
disp(['keeping ' num2str(sum(keepind)) ' members'])
keepind=find(keepind==1);
u=u(:,:,keepind);
v=v(:,:,keepind);
us=us(:,:,keepind);
vs=vs(:,:,keepind);
D=D(:,:,keepind);
H=H(:,:,keepind);
z=z(:,:,keepind);
m=m(:,:,keepind);
h=h(:,:,keepind);
f=f(keepind);
memberID=memberID(keepind);
[ny,nx,N]=size(h);

%-----------------------------------------------------
% observations: struct 'meas'
%-----------------------------------------------------

% standard directory structure: observational data file is in parent
% directory and has name obsData.mat
load([membersDir '../../obsData.mat'])

% interpolate ensemble to obs-points
f=fields(meas);
for i=1:length(f)
  obs=getfield(meas,f{i});

  % special case: frequency-wavenumber data comes from dispersion reln,
  % not from swan/shorecirc output files.  Also note, need to estimate
  % wavenumber based on total water depth, not bathymetry; and, need to
  % include the effect of currents.
  if(f{i}=='k')
    if(~isempty(meas.k.data))
      bcdata=load([membersDir '../../../bcData.mat']);
      clear hi
      sigma=2*pi*obs.f;
      g=9.8126;

      % define dispersion relationship factors, following KD86 model in
      % Catalan and Haller (2008)
      if(1)
        DD=@(kk,hh)(8+cosh(4*kk.*hh)-2*tanh(kk.*hh).^2)...
          ./(8*sinh(kk.*hh).^4);
        f1=@(kk,hh)tanh(kk.*hh).^5;
        f2=@(kk,hh)(kk.*hh./sinh(kk.*hh)).^4;
        ee=@(kk,HH)kk.*HH/2;
      else
        warning('using linear dispersion, no currents')
        DD=@(kk,hh)0;
        f1=@(kk,hh)0;
        f2=@(kk,hh)0;
        ee=@(kk,HH)0;
      end
      Uk=@(uu,vv,kk,gam)sqrt(uu.^2+vv.^2).*kk.*cosd(gam);

      % interpolate for variables in dispersion reln
      for n=1:N
        hi(:,n)=interp2(x,y,h(:,:,n)+z(:,:,n),obs.x,obs.y)+bcdata.tide;
        ui(:,n)=interp2(x,y,u(:,:,n)-us(:,:,n),obs.x,obs.y);
        vi(:,n)=interp2(x,y,v(:,:,n)-vs(:,:,n),obs.x,obs.y);
        Hi(:,n)=interp2(x,y,H(:,:,n),obs.x,obs.y);
        alpha=interp2(x,y,90-D(:,:,n),obs.x,obs.y);
        beta=rad2deg(atan2(-vi(:,n),-ui(:,n)));
        gamma(:,n)=beta-alpha;  % angle btwn wave and current
      end
      ui(isnan(ui))=0;
      vi(isnan(vi))=0;
      Hi(isnan(Hi))=0;
      gamma(isnan(gamma))=90;
      kguess=approxDispersion(repmat(sigma,[1 N]),hi);

      % discard locations for which no ensemble wavenumber estimate is
      % possible
      ind=find(~isnan(sigma.*sum(Hi,2).*sum(hi,2).*sum(kguess,2)) & ...
               sum(hi<0.1,2)==0 & sum(Hi==0,2)==0 & ...
               sum(kguess<0,2)==0);
      tmp=fields(obs);
      for j=1:length(tmp)
        data=getfield(obs,tmp{j});
        data=data(ind);
        obs=setfield(obs,tmp{j},data);
      end
      sigma =sigma(ind)   ;
      kguess=kguess(ind,:);
      hi    =hi(ind,:)    ;
      Hi    =Hi(ind,:)    ;
      ui    =ui(ind,:)    ;
      vi    =vi(ind,:)    ;
      gamma =gamma(ind,:) ;

      % solve dispersion relation
      disp('solving dispersion relation')
      model=nan(length(obs.x),N);
      sigma=repmat(sigma,[1 N]);
      parfor j=1:length(model(:))
        solveme=@(k)(sigma(j)-Uk(ui(j),vi(j),k,gamma(j)))^2 ...
                - g*k.*tanh( k.*hi(j)+f2(k,hi(j)).*ee(k,Hi(j)) ) ...
                .*(1+f1(k,hi(j)).*ee(k,Hi(j)).^2.*DD(k,hi(j)));
        model(j)=fzero(solveme,kguess(j),optimset('Display','off'));
      end
      obs.model=model;

    else  % catch for nodata in obs struct
      obs.model=[];
    end

  % if not f-k data, just interpolate from model output
  else
    % assign variable name for a given field
    if(f{i}=='r')  % IR currents are v-velocity
      vname='v';
    elseif(f{i}=='h')  % bathymetry data
      vname='h';
    elseif(f{i}=='a')  % wave angle data
      vname='a';
    else  % all other data are v-velocity
      vname='v';
    end
    eval(['data=' vname ';']);
    if(~isempty(obs.x))
      for n=1:N
        obs.model(:,n)=interp2(x,y,data(:,:,n),obs.x,obs.y);
      end
    else   % catch for nodata in obs-struct
      obs.model=[];
    end
  end

  meas=setfield(meas,f{i},obs);

end  % loop on fields in obs-struct (i.e. loop on obs types)

%-----------------------------------------
% merge together measurement types into one big vector
%-----------------------------------------

% save original structure
measStruct=meas;

% remove null fields
f=fields(meas);
for i=1:length(f)
  if(isempty(getfield(getfield(meas,f{i}),'data')))
    meas=rmfield(meas,f{i});
  end
end

% for non-wavenumber observations, add a fake "frequency" field with a
% unique ID.  This will be used to assign "groups" of variables which
% have sptially correlated observation error
f=fields(meas);
for i=1:length(f)
  if(f{i}~='k')
    this=getfield(meas,f{i});
    if(f{i}=='v')
      this.f=zeros(size(this.data))-1;
    elseif(f{i}=='r')
      this.f=zeros(size(this.data))-2;
    else
      error(['unknown data type ' f{i}])
    end
    meas=setfield(meas,f{i},this);
  end
end

% lump observations together
meas2.x=[];
meas2.y=[];
meas2.data=[];
meas2.s=[];
meas2.model=[];
meas2.f=[];
f=fields(meas);
for i=1:length(f)
  thismeas=getfield(meas,f{i});
  meas2.x=[meas2.x; thismeas.x];
  meas2.y=[meas2.y; thismeas.y];
  meas2.data=[meas2.data; thismeas.data];
  meas2.s=[meas2.s; thismeas.s.*ones(size(thismeas.data))];
  meas2.model=[meas2.model; thismeas.model];
  meas2.f=[meas2.f; thismeas.f];
end
meas=meas2;
clear meas2

% remove any non-real data points
if(~isempty(meas.data))
  validind=find(~isnan(sum(meas.model,2).*meas.data) & ...
                sum(imag(meas.model),2)==0 & ...
                imag(meas.data)==0);
  f=fields(meas);
  for i=1:length(f)
    data=getfield(meas,f{i});
    if(isvector(data))
      data=data(validind); data=data(:);
    else
      data=data(validind,:);
    end
    meas=setfield(meas,f{i},data);
  end
end

% auto-QC step: remove any measurements which disagree with model by more
% than X stdev's
if(length(meas.data)>0)
  qcFact=3;
  ind=find(abs(meas.data-mean(meas.model,2))<=sqrt(meas.s.^2+var(meas.model,[],2))*qcFact);
  ntoss=length(meas.data)-length(ind);
  if(ntoss>0)
    disp(['auto-QC: throwing out ' num2str(ntoss) ' observations'])
    f=fields(measStruct);
    for i=1:length(f)
      obs=getfield(measStruct,f{i});
      tmp=sum(abs(obs.data-mean(obs.model,2))>...
              sqrt(obs.s.^2+var(obs.model,[],2))*qcFact);
      disp(['        (' num2str(tmp) ' out of ' ...
            num2str(length(obs.data)) ' from ' f{i} ')'])
    end
    meas.x=meas.x(ind);
    meas.y=meas.y(ind);
    meas.data=meas.data(ind,:);
    meas.s=meas.s(ind,:);
    meas.model=meas.model(ind,:);
    meas.f=meas.f(ind);
  else
    disp(['no auto-QC tripped, keeping ' num2str(length(meas.data)) ...
          ' observations'])
  end
end

%-----------------------------------------
% assimilate
%-----------------------------------------

% special case: no data to assimilate
if(isempty(meas.data))
  disp('WARNING: no data to assimilate, setting posterior=prior')
  ha=h;
else

% vectorized model
[ny,nx,N]=size(h);
hvec=reshape(h,[nx*ny N]);
uvec=reshape(u,[nx*ny N]);
vvec=reshape(v,[nx*ny N]);

% covariances
disp('computing observation covariance')
warning('using obs error spatial covariance')
% Observation errors: simulate spatially correlated errors with specified
% variance (ensemble 'meas.noise').  Apply regularization to the
% sample covariance matrix to estimate Cdd.  Also, retain the
% sample for later use when perturbing the ensemble
% of observations during the update step (see below)
%
% To simulate the spatially correlated noise, I am using the following
% method:
%
%  a) generate realizations of noise having the desired spatial
%  correlation scales, and *unit* standard deviation.  Stored in 'meas.noise'
%
%  b) scale the unit-stdev spatially correlated noise such that stdev is
%  equal to the desired value at each obs-point.
%
%  c) estimate covariance from the above ensemble of correlated noise.
%  Use shrinkage to regularize the result.
%
% TODO:
%
% - consider including correlation in errors between different frequency
% bands in cBathy.  This version treats each frequency band as an
% independent observation
%
% - This version looks for frequency bands, then assigns independent
% spatially correlated errors for each band.  A more elegant way would be to
% develop observation errors separately for each instrument type, since not
% all instruments have frequency bands
%
warning off
addpath ~/work/randomField/%generate[12]D.m
warning on
% define a hi-res grid on which to compute the covariances
dxc=2; dyc=2;
xc=(min(meas.x)-dxc):dxc:(max(meas.x)+dxc);
yc=(min(meas.y)-dyc):dyc:(max(meas.y)+dyc);
if(mod(length(yc),2)==1) yc=[yc yc(end)+dyc]; end
nxc=length(xc);
nyc=length(yc);
% initialize
meas.noise=nan*meas.model;
% for each frequency...
f=unique(meas.f);
nf=length(f);
for c=1:nf
  disp(['frequency ' num2str(c) ' of ' num2str(nf)])
  ind=find(meas.f==f(c));

  % vbar observations, flagged using "frequency = -1"
  if(f(c)==-1)
    Ldy=15;  % cutoff at 30m vbar window (Ldy is window radius)
    % for each alongshore transect in vbar observation set...
    xx=unique(meas.x(ind));
    for i=1:length(xx)
      % define observation indeces for this vbar transect, indy.
      % Note, eventually indy will touch every vbar observation, because
      % of loop over transects
      indy=find(meas.x(ind)==xx(i));
      indy=ind(indy);
      % compute N realizations of spatially correlated noise with unit
      % stdev, for this transect
      dd=generate1D(nyc,dyc,Ldy,1,N);
      % interpolate the correlated noise to the obs-points, for this set
      % of observations
      for n=1:N
        meas.noise(indy,n)=interp1(yc,dd(:,n),meas.y(indy));
      end
    end  % loop on alongshore transects

  % IR-PIV observations
  elseif(f(c)==-2)
    % from C. Chickadel: "The analysis window size was twice the output grid
    % spacing (50% overlap)".  For SZO, the grid spacing was 8x8 meters.
    % Note Ld[xy] are radius, not diameter, so they should be set equal to
    % grid spacing.
    Ldx=8;
    Ldy=8;
    % generate N realizations of spatially correlated noise
    dd=generate2Dp(nxc,nyc,dxc,dyc,Ldx,Ldy,1,N);
    % interpolate the correlated noise to the obs-points
    for n=1:N
      meas.noise(ind,n)=interp2(xc,yc,dd(:,:,n),meas.x(ind),meas.y(ind));
    end

  % wavenumber band observations
  elseif(f(c)>0)
    Ldx=20;  % cutoff at 20m window size in x
    Ldy=50;  % cutoff at 50m window size in y
    % generate N realizations of spatially correlated noise
    dd=generate2Dp(nxc,nyc,dxc,dyc,Ldx,Ldy,1,N);
    % interpolate the correlated noise to the obs-points
    for n=1:N
      meas.noise(ind,n)=interp2(xc,yc,dd(:,:,n),meas.x(ind),meas.y(ind));
    end

  % other observations, not yet supported
  else
    error('unknown observation type')
  end

end  % loop on frequency bands and/or observation types
% scale the correlated noise to match the desired observation stdev
meas.noise=meas.noise.*repmat(meas.s./std(meas.noise,[],2),[1 N]);
% estimate covariance matrix from ensemble of noise realizations
Cdd=myCov(meas.noise,meas.noise);
%Cdd=omega_vec(Cdd,covDist(meas.x,meas.y,meas.x,meas.y),30);   % localize
C0=diag(meas.s.^2);
delta=0.1;
Cdd=delta*C0+(1-delta)*Cdd;  % regularize (shrinkage method)

% compute model covariances
Chv=myCov(hvec,meas.model);
LCvvL=myCov(meas.model,meas.model);

% localize model covariances
L=75;
disp(['localizing covariances, localization length ' ...
      num2str(L) 'm'])
[xg,yg]=meshgrid(x,y);
Chv=omega_vec(Chv,covDist(xg(:),yg(:),meas.x,meas.y),L);
Cdd=omega_vec(Cdd,covDist(meas.x,meas.y,meas.x,meas.y),L);
LCvvL=omega_vec(LCvvL,covDist(meas.x,meas.y,meas.x,meas.y),L);

% assimilate for posterior ensemble.  Add random noise to
% observations to ensure correct posterior ensemble covariance
disp('assimilating for posterior ensemble')
ChvinvC=Chv*inv(LCvvL+Cdd);
for n=1:N
  %meas.noise=meas.s.*randn(size(meas.s));
  dh(:,:,n)=reshape(ChvinvC*(meas.data+meas.noise(:,n)-meas.model(:,n)),[ny nx]);
  ha(:,:,n)=h(:,:,n)+dh(:,:,n);
end

end  % gate for special no-data case

% note: masked points should never get updated, they should always remain
% dry by definition
ha(isnan(ha))=-999;

% % test code: avoid making extrapolated depth corrections for points outside
% % of the observation domain
% disp('masking out corrections outside of observation window')
% clear dst
% pp=convhull(meas.x,meas.y);
% [xg,yg]=meshgrid(x,y);
% for i=1:nx
%   for j=1:ny
%     dst(j,i)=distancePointPoly([xg(j,i) yg(j,i)],[meas.x(pp) meas.y(pp)]);
%   end
% end
% ind=inpolygon(xg,yg,meas.x(pp),meas.y(pp));
% dst(ind)=-dst(ind);
% decayLen=100;
% w=-.5*sin(pi*dst/decayLen)+.5;
% w(dst<=-decayLen/2)=1;
% w(dst>=decayLen/2)=0;
% dh=ha-h;
% [ny,nx,N]=size(h);
% dh=dh.*repmat(w,[1 1 N]);
% ha=h+dh;

% write the posterior ensemble to disk
if(~isempty(outdir))
  [xg,yg]=meshgrid(x,y);
  disp('writing posterior ensemble to disk')
  ncoutdir=[outdir 'netcdf/'];
  datoutdir=[outdir 'dat/'];
  mkdir(ncoutdir)
  mkdir(datoutdir)
  for n=1:N
    if(floor(n/N*10)>floor((n-1)/N*10))
      disp([num2str(floor(n/N*10)*10) '%'])
    end
    fout=[ncoutdir 'bath' num2str(memberID(n)) '.nc'];
    mat2roms(xg,yg,ha(:,:,n),fout,-99)
    fout=[datoutdir 'bath' num2str(memberID(n)) '.dat'];
    hn=flipud(fliplr(ha(:,:,n)));
    save('-ascii',fout,'hn')
  end
end

% define outputs
%h(m==0)=nan;
%ha(m==0)=nan;
out.x=x;
out.y=y;
out.prior.h=nanmean(h,3);
out.prior.u=nanmean(u,3);
out.prior.v=nanmean(v,3);
out.prior.hstd=nanstd(h,[],3);
out.post.h =nanmean(ha,3);
out.post.hstd=nanstd(ha,[],3);
out.meas=measStruct;

matlabpool close
