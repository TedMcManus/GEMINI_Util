function [] = autocatalogue(ref,outdir,pars,noFields,noParticles)

arguments
    ref (1,1) string
    outdir (1,1) string
    pars (1,1) struct
    noFields (1,1) logical = false
    noParticles (1,1) logical = false
end
outdir=char(outdir);
ref=char(ref);

% Unpackage the structure 
flagdirich=pars.flagdirich;
arc=pars.arc;
field_cadence=pars.field_cadence;
particle_cadence=pars.particle_cadence;
isMoving=pars.isMoving;
E0floor=pars.E0floor;
Qfloor=pars.Qfloor;
ExBg=pars.ExBg;
EyBg=pars.EyBg;
numcores=pars.numcores;
runname=pars.runname;
runname=char(runname);
flattener=pars.flattener;
xy=pars.xy;
invert=pars.invert;
%% Values to be hard-coded:

x2d = 1.84742e-2; % deg km^-1 Converts west-east distances to degrees. Evaluated at r = 7375.36 km (~1000 km alt)
y2d = 7.76854e-3; % deg km^-1 Converts south-north distances to degrees. Evaluated at r = 7375.36 km (~1000 km alt)

%REFERENCE GRID TO USE

if exist('ref','var')
    a = char(ref);
    %remove the trailing slash if it exists (convenience code)
    if a(length(a))=='/'
        ind = length(a)-1;
        a = a(1:ind);
        ref = char(a);
    end
    %set the location of config and grid files
    direcconfig = ref;
    direcgrid = ref;
else
    error('please specify reference sim')
end


%CREATE SOME SPACE FOR OUTPUT FILES
if ~exist('outdir','var')
    error('please specify output directory')
else
    a = char(outdir);
    if a(length(a))=='/' %take off trailing slash
        ind = length(a)-1;
        a = a(1:ind);
        outdir = char(a);
    end
    out = string(outdir);
    if ~exist(out,'dir')
        mkdir(out) %make the output directory
    end
end


%READ IN THE SIMULATION INFORMATION
if (~exist('ymd0','var'))
    %built-in gemini function loads config information into a struct
    confdat=gemini3d.read.config(direcconfig);
    fprintf('Input config.dat file loaded.\n');
end

%Set up time indexing
ymd0=confdat.ymd; %year, month, day of sim
UTsec0=confdat.UTsec0; %initial timestep
tdur=confdat.tdur;

%Load the simulation grid
%use a built-in Gemini function to load the grid into a struct
xg=gemini3d.read.grid(direcgrid);
%lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
fprintf('Grid loaded.\n');

%Create the grid
MLAT=90-squeeze(xg.theta(1,:,:))*180/pi; %magnetic latitude
MLON=squeeze(xg.phi(1,:,:))*180/pi; %magnetic longitude
%lengths of mlat and mlon arrays
llon=xg.lx(2);
llat=xg.lx(3);

%grab two slices of the data and calculate their mean
mlon=MLON(:,1);
mlat=MLAT(1,:);
mlonmean=mean(mlon);
mlatmean=mean(mlat);


%Time indexing
dt = 1;
tmin=0;
tmax=tdur;
run_time=tmin:dt:tmax;
lt=length(run_time);


%Terrible terrible code to make stuff 3D 
%Here, we're turning our 2-D arrays (MLAT and MLON) into 3-d arrays where
%the dimensions are lon x lat x time

%preallocate the arrays
mlat_tmp=ones(llon,llat,lt);
mlon_tmp=ones(llon,llat,lt);
for i = 1:lt
    %store a copy of MLAT at every single timestep (the grid does not vary
    %with time, so all timesteps are identical)
    mlat_tmp(:,:,i)=MLAT;
    mlon_tmp(:,:,i)=MLON;
end
MLAT=mlat_tmp;
MLON=mlon_tmp;

%Now make the most stupid, ridiculuous array you could ever imagine: An
%array with dimensions lat x lon x time (so like 144x216x301), where every
%slice on the third dimension (time) contains a 144x216 array filled with
% copies of the current timestep

IT=ones(llon,llat,lt);
for i = 1:lt
    IT(:,:,i)=IT(:,:,i).*run_time(i); %yep, it's that ridiculous
end


if matches(arc,'gudda') || matches(arc,'')
    % Values to be input by the user:
    [mapQit,mapE0it,mapJ,mapU,mapV] = GUDDA_map(MLON,MLAT,IT,mlonmean,mlatmean,x2d,y2d,pars);
    mapJ=mapJ.*10^-6; %This line of code is the most important one in the whole dang function, it turns current into amps rather than uA
elseif matches(arc,'gmdda')
    [mapQit,mapE0it,mapJ,mapU,mapV] = GMDDA_map(MLON,MLAT,IT,mlonmean,mlatmean,x2d,y2d,pars);
    mapJ=mapJ.*10^-6;
elseif matches(arc,'angle')
    [mapQit,mapE0it,mapJ,mapU,mapV] = ANGLE_map(MLON,MLAT,IT,mlonmean,mlatmean,pars);
    mapJ=mapJ.*10^-6;
elseif matches(arc,'gudda_uneven')
    [mapQit,mapE0it,mapJ,mapU,mapV] = GUDDA_map_new(MLON,MLAT,IT,mlonmean,mlatmean,x2d,y2d,pars);
    mapJ=mapJ.*10^-6.*-1;
elseif matches(arc,'STEVE')
    %For STEVE runs. The precip is defined independent of the fields
    [mapE0it,mapQit]=STEVE_particles_new(pars.STEVE.center,pars.STEVE.E0max,pars.STEVE.E0min,...
        pars.STEVE.Qmax,pars.STEVE.Qmin,ref,outdir,pars.STEVE.mlatsig,true);
    %This is a field definition that is just for STEVE runs. It creates
    %a large gaussian envelope at the target mlat value, and performs a
    %recursive guess-and-check search
    isMoving=false;
    pot_out = flow_spec(ref,pars.STEVE.wtarget,pars.STEVE.displace,pars.STEVE.mlonsig,...
        pars.STEVE.index,pars.STEVE.flagdirich,pars.STEVE.ExBG,...
        pars.STEVE.EyBG,pars.STEVE.lbound,pars.STEVE.ubound,pars.STEVE.steps,xg,pars.STEVE.vtarg);
elseif matches(arc,'Archer')
    %For STEVE runs. The precip is defined independent of the fields
    [mapE0it,mapQit]=STEVE_particles_new(pars.Archer.center,pars.Archer.E0max,pars.Archer.E0min,...
        pars.Archer.Qmax,pars.Archer.Qmin,ref,outdir,pars.Archer.mlatsig,true);
    %This is a field definition that is just for STEVE runs. It creates
    %a large gaussian envelope at the target mlat value, and performs a
    %recursive guess-and-check search
    isMoving=false;
    mapU = flow_spec_archer(ref,pars.Archer.wtarget,pars.Archer.displace,pars.Archer.mlonsig,xg,abs(pars.Archer.vtarg));
    if pars.Archer.vtarg<0
        mapU=-mapU;
    end
    mapU=repmat(mapU,1,1,301);
    mapV=zeros(size(mapU));
    mapJ=zeros(size(mapU));
elseif matches(arc,'maeve')
    [X2,X3,IT] = ndgrid(xg.x2(3:end-2),xg.x3(3:end-2),0:1:confdat.tdur);
    [mapQit,mapE0it,mapJ,mapU,mapV] = MAEVE_map(X2,X3,IT,pars);
end

itplot=10;
pltsv = 0; % Save plot?
fntsz=10;

if ~exist('mapJ','var')&&~exist('mapU','var')&&~exist('mapV','var')
    mapJ=zeros(size(mapE0it));
    mapU=zeros(size(mapE0it));
    mapV=zeros(size(mapE0it));
end

%plot the results
figure(1)
set(gcf,'units','normalized','Position',[0.01 0.1 0.98 0.8])
set(gcf,'PaperPosition',[0,0,6,4],'PaperUnits','inches')
t = subplot(2,3,1);

pcolor(MLON(:,:,itplot),MLAT(:,:,itplot),mapQit(:,:,itplot));
set(gca,'FontSize',fntsz)
shading flat
c = colorbar;
title('Energy Flux, Qit','FontSize',fntsz)
xlabel('mlon [°]','FontSize',fntsz)
ylabel('mlat [°]','FontSize',fntsz)
ylabel(c,'mW m^{-2}','FontSize',fntsz)

t = subplot(2,3,2);
pcolor(MLON(:,:,itplot),MLAT(:,:,itplot),mapE0it(:,:,itplot));
set(gca,'FontSize',fntsz)
shading flat
c = colorbar;
title('Average Energy, E_0','FontSize',fntsz)
xlabel('mlon [°]','FontSize',fntsz)
ylabel('mlat [°]','FontSize',fntsz)
ylabel(c,'eV','FontSize',fntsz)

t = subplot(2,3,3);
pcolor(MLON(:,:,itplot),MLAT(:,:,itplot),mapJ(:,:,itplot)); %Force sign convention
set(gca,'FontSize',fntsz)
shading flat
c = colorbar;
title('Downward FAC, J_{||}','FontSize',fntsz)
xlabel('mlon [°]','FontSize',fntsz)
ylabel('mlat [°]','FontSize',fntsz)
ylabel(c,'A m^{-2}','FontSize',fntsz)
colormap(t,colorcet('D1A'))

t = subplot(2,3,4);
pcolor(MLON(:,:,itplot),MLAT(:,:,itplot),mapU(:,:,itplot));
set(gca,'FontSize',fntsz)
shading flat
c = colorbar;
title('Eastward Flow, U','FontSize',fntsz)
xlabel('mlon [°]','FontSize',fntsz)
ylabel('mlat [°]','FontSize',fntsz)
ylabel(c,'m s^{-1}','FontSize',fntsz)

t = subplot(2,3,5);
pcolor(MLON(:,:,itplot),MLAT(:,:,itplot),mapV(:,:,itplot));
set(gca,'FontSize',fntsz)
shading flat
c = colorbar;
title('Northward Flow, V','FontSize',fntsz)
xlabel('mlon [°]','FontSize',fntsz)
ylabel('mlat [°]','FontSize',fntsz)
ylabel(c,'m s^{-1}','FontSize',fntsz)
if pltsv
    print(arc + "_" + string(itplot) + "s.png",'-dpng')
end


%% Define proper E-field field names and values
if ~noFields
    %set up time indexing
    dt = field_cadence;
    tmin=0;
    tmax=tdur;
    field_time=tmin:dt:tmax;
    lt=length(field_time);

    ymd=ymd0;
    UTsec0=35850;
    if UTsec0==36000
        UTsec=UTsec0+field_time-tdur/2;     %time given in file is the seconds from beginning of hour
    else
        UTsec=UTsec0+field_time;
    end

    %this code was a nightmare to write, and it converts the year, month, day,
    %and UTsec values into Matlab's "datetime" format. I will not commment it
    %because I don't want to understand it again.
    UThrs=UTsec/3600;
    expdate=cat(2,repmat(ymd,[lt,1]),UThrs(:),zeros(lt,1),zeros(lt,1));
    hours=mod(expdate(:,4),1);
    expdate(:,4)=expdate(:,4)-hours;
    expdate(:,5)=hours.*60;
    minute=mod(expdate(:,5),1);
    expdate(:,5)=expdate(:,5)-minute;
    expdate(:,6)=minute.*60;
    E_times=datetime(int32(expdate));

    %background field data, usually just set to zero
    Exit=ones(llon,llat).*ExBg;
    Eyit=ones(llon,llat).*EyBg;

    %Now make it 3-D for consistency
    Ex_t=ones(llon,llat,lt);
    Ey_t=ones(llon,llat,lt);
    for i = 1:lt
        Ex_t(:,:,i)=Exit;
        Ey_t(:,:,i)=Eyit;
    end
    Exit=Ex_t;
    Eyit=Ey_t;

    %% FIELD BOUNDARY CONDITIONS

    %preallocate the arrays to save time
    Vminx1it=zeros(llon,llat,lt); %minimum values of field input (usually just zeros)
    Phi=zeros(llon,llat,lt); %I use this as a temporary variable to hold data. This should prob be changed
    Vmaxx1it=zeros(llon,llat,lt); %The final field input, as it will be written to the input file
    if matches(arc,'STEVE') %then the field data is unchanging, so we just save a copy of it into each timestep
        for i=1:lt
            Vmaxx1it(:,:,i)=pot_out;
        end
    else
        if flagdirich==0 %then we're using current as an input
            for i = 2:lt
                Vmaxx1it(:,:,i)=mapJ(:,:,field_time(i));
            end
        else
            if numcores>1 %use parallel processing. Requires an add-on, but is much faster
                try parpool(numcores);
                catch('Parallel processing not enabled \n');
                end
            end
            if isMoving
                for it=1:lt
                    if it==1 %we always write zeros into the first timestep. IDK why
                        Phivals=zeros(llon,llat);
                        Phi(:,:,it)=Phivals;
                        disp(['it= ',num2str(field_time(it)),'. ',num2str(lt-it),' timesteps remaining']);
                    else
                        tic %for timing

                        %Here's a slight issue - we have a flow structure that we
                        %want, but GEMINI only understands potential. Potential is the integral of flow,
                        % so we need to use a numerical integrator in order to get the potential structure.
                        %This is Jules' script, so contact him if you need to understand it.
                        if numcores>1
                            [Phivals,~,~,~]=flow2phi_new(mapU(:,:,it*10-10),mapV(:,:,it*10-10),xy(1),xy(2),xg,outdir,64,false,true,false,true,false,it);
                        else
                            [Phivals,~,~,~]=flow2phi_new(mapU(:,:,it*10-10),mapV(:,:,it*10-10),xy(1),xy(2),xg,outdir,64,false,true,false,false,false,it);
                        end
                        Phi(:,:,it)=Phivals;
                        out_time = toc;
                        disp(['it= ',num2str(field_time(it)),'. ',num2str(lt-it),' timesteps remaining. \nReconstructed in ',num2str(out_time),'seconds']);
                    end
                    Vmaxx1it=Phi;
                end
            else
                %The run isn't moving. Save time by using the integrator with high fidelity on one timestep
                %and just copying it across the entire run
                disp('WARNING: Solving for stationary boundary conditions')
                tic %for timing
                if numcores>1
                    [Phivals,~,~,~]=flow2phi_new(mapU(:,:,5*10-10),mapV(:,:,5*10-10),xy(1),xy(2),xg,outdir,64,false,true,false,true,false,5);
                else
                    [Phivals,~,~,~]=flow2phi_new(mapU(:,:,5*10-10),mapV(:,:,5*10-10),xy(1),xy(2),xg,outdir,64,false,true,false,false,false,5);
                end
                for iter=1:lt
                    Vmaxx1it(:,:,iter)=Phivals;
                end
            end
            Vmaxx1it(:,:,1)=zeros(llon,llat);
        end
    end
    %idk what these do, but they're all zero
    Vminx2ist=zeros(llat,lt);
    Vmaxx2ist=zeros(llat,lt);
    Vminx3ist=zeros(llon,lt);
    Vmaxx3ist=zeros(llon,lt);
    flagdirich=ones(lt).*flagdirich;
    if strcmp(invert,'Y')
        Vmaxx1it=-Vmaxx1it;
    end

    %% Prepare Fields
    %all the field data is stored in a structure with very particular names, which is then passed to a default
    %gemini script in order to be processed and written out in the correct format
    E=struct();
    E.flagdirich=flagdirich;
    E.Exit=Exit;
    E.Eyit=Eyit;
    E.Vminx1it=Vminx1it;
    E.Vmaxx1it=Vmaxx1it;
    E.Vminx2ist=Vminx2ist;
    E.Vmaxx2ist=Vmaxx2ist;
    E.Vminx3ist=Vminx3ist;
    E.Vmaxx3ist=Vmaxx3ist;
    E.mlon=double(mlon);
    E.mlat=double(mlat);
    E.llon=double(llon);
    E.llat=double(llat);
    E.times=E_times;
end
%% Prepare Precip
if ~noParticles
    %time indexing setup
    dt = particle_cadence;
    tmin=0;
    tmax=tdur;
    particle_time=tmin:dt:tmax;
    lt=length(particle_time);
    ymd=ymd0;
    if UTsec0==36000
        UTsec=UTsec0+particle_time-tdur/2;     %time given in file is the seconds from beginning of hour
    else
        UTsec=UTsec0+particle_time;
    end
    %again, do the datetime conversion
    UThrs=UTsec/3600;
    expdate=cat(2,repmat(ymd,[lt,1]),UThrs(:),zeros(lt,1),zeros(lt,1));
    hours=mod(expdate(:,4),1);
    expdate(:,4)=expdate(:,4)-hours;
    expdate(:,5)=hours.*60;
    minute=mod(expdate(:,5),1);
    expdate(:,5)=expdate(:,5)-minute;
    expdate(:,6)=minute.*60;
    P_times=datetime(int32(expdate));

    %Now set up the struct in the way Gemini likes it
    P=struct;
    P.mlon=double(mlon);
    P.mlat=double(mlat);
    P.llon=double(llon);
    P.llat=double(llat);

    %Define a "floor" for precip. All values lower than the floor are lifted up to the background value
    E0it=zeros(llon,llat,lt);
    Qit=zeros(llon,llat,lt);

    for it=1:lt
        current_time=particle_time(it);
        if current_time==0
            current_time=1; %1-based indexing has its drawbacks
        end
        %use logical indexing & replace all values less than E0floor or Qfloor with E0floor or Qfloor, respectively
        tmpvar=mapE0it(:,:,current_time);
        tmpvar(tmpvar<E0floor)=E0floor;
        E0it(:,:,it)=tmpvar;
        tmpvar=mapQit(:,:,current_time);
        tmpvar(tmpvar<Qfloor)=Qfloor;
        Qit(:,:,it)=tmpvar;
    end

    if flattener
        E0it=ones(size(E0it)).*E0floor;
        Qit=ones(size(Qit)).*Qfloor;
    end
    P.E0it=E0it;
    P.Qit=Qit;
    P.times=P_times;

    %% Write the inputs
    fprintf('\n----------------------Writing Inputs----------------------\n\n');
    fprintf('Particles: ');
    pdir=[outdir,'/',runname,'_particles'];
    mkdir(pdir);
    gemini3d.write.precip(P,pdir,'h5');
    save([pdir,'/particles.mat'],'mlon','mlat','MLAT','MLON','E0it','Qit','P_times','P');
end

if ~noFields
    fprintf('\nFields: ');
    fdir = [outdir,'/',runname,'_fields'];
    mkdir(fdir);
    gemini3d.write.Efield(E,fdir,'h5');
    save([fdir,'/fields.mat'],'mlon','mlat','MLAT','MLON','Exit','Eyit','Vminx*','Vmax*','E_times','mapU','mapV','mapJ','E');
    flagdirich=flagdirich(1);
end
pars.flagdirich=flagdirich(1);
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate')); %close the parallel pool if we're using multithreading
end
fclose('all'); %close all files, just in case
save([outdir,'/setup_meta.mat'],'pars');

%% Input Functions
    % Growth-phase Unipolar Double Discrete Arc (GUDDA)
    function [Qit,E0it,J,U,V] = GUDDA_map(mlon,mlat,it,mlonmean,mlatmean,x2d,y2d,pars)
        par = pars.gudda;
        mlon = mlon - par.driftE*x2d*it;
        mlat = mlat - par.driftN*y2d*it;
        offset = 0*(3*par.sheetwidth/4-par.arcsep/2)*y2d; % move southern arc to equatorward edges of current sheets
        Jpk = 1e3*par.Kpk/par.sheetwidth; % convert line integrated current to peak current density
        mlatctr = mlatmean+(par.spanNS/2).*tanh((mlon-mlonmean)/par.mlonsig); % arc contour definition
        dmlatctrdmlon = (par.spanNS/2).*sech((mlon-mlonmean)/par.mlonsig).^2/par.mlonsig;% used to define tangent/normal to contour.
        dydx = (x2d/y2d).*dmlatctrdmlon;  % Derivative defined in physical space
        %         dydx0 = (x2d/y2d)*(par.spanNS/2)/par.mlonsig; % slope at center
        wfac = sqrt(1+dydx.^2); % width and separation factor to avoid pinching at steep curves
        %         wfac = 1; %%% CHANGE
        %         sfac = (sqrt(1+dydx0^2)-1)/dydx0;
        arcsep = par.arcsep*y2d;
        arcwidth = par.arcwidth*y2d*wfac;
        flowwidth = par.flowwidth*y2d*wfac;
        sheetwidth  = par.sheetwidth*y2d;
        Qit = (par.Qitpk + par.Qitslope*(mlon-mlonmean)).*(...
            gaussian(mlat-mlatctr,-arcsep/2,arcwidth)...
            +gaussian(mlat-mlatctr, arcsep/2,arcwidth)...
            );%.*gaussian(mlon,mlonmean,par.mlonsig);
        E0it = par.E0itpk.*(...
            gaussian(mlat-mlatctr,-arcsep/2,arcwidth)...
            +gaussian(mlat-mlatctr, arcsep/2,arcwidth)...
            );%.*gaussian(mlon,mlonmean,par.mlonsig);
        J = (Jpk + par.Jslope*(mlon-mlonmean)).*(...
            tanh(50*flowwidth.*(mlat-mlatctr-(-3/2)*sheetwidth-offset))...
            -2*tanh(50*flowwidth.*(mlat-mlatctr-(-1/2)*sheetwidth-offset))...
            +tanh(50*flowwidth.*(mlat-mlatctr-( 1/2)*sheetwidth-offset))...
            )/2;%.*gaussian(mlon,mlonmean,par.mlonsig);
        J = J.*(it>2); % Turn on current 2 seconds after precip has settled -- jvi: should this value be in some init file?
        F = (par.Fpk + 1e3*par.Jslope*(mlon-mlonmean)).*(...
            gaussian(mlat-mlatctr,(-3/2)*sheetwidth+offset,flowwidth)... % flow magnitude
            -gaussian(mlat-mlatctr,(-1/2)*sheetwidth+offset,flowwidth)...
            +gaussian(mlat-mlatctr,( 1/2)*sheetwidth+offset,flowwidth)...
            );%.*gaussian(mlon,mlonmean,par.mlonsig);
        U = F./sqrt(1+dydx.^2); % ~eastward flow
        V = F.*dydx./sqrt(1+dydx.^2); % ~northward flow
        function [y] = gaussian(x,pos,fwhm)
            y = exp(-(x-pos).^2./(0.36*(fwhm).^2));
        end
    end
    
    % Growth-phase Multipolar Double Discrete Arc (GMDDA)
    function [Qit,E0it,J,U,V] = GMDDA_map(mlon,mlat,it,mlonmean,mlatmean,x2d,y2d,pars)
        par = pars.gmdda;
        par.sheetwidth = par.arcsep/2; % ensure precip in correct current sheet
        offset = par.sheetwidth*y2d/4; % move arcs to equatorward edges of current sheets
        Jpk = 1e3*par.Kpk/par.sheetwidth; % convert line integrated current to peak current density
        mlatctr = mlatmean+tanh((mlon-mlonmean)/par.mlonsig);
        dmlatctrdmlon = sech((mlon-mlonmean)/par.mlonsig).^2/par.mlonsig;% used to define tangent/normal to contour.
        dydx = (x2d/y2d).*dmlatctrdmlon;  % Derivative defined in physical space
        Qit = par.Qitpk.*(gaussian(mlat-par.drift*y2d*it-mlatctr,-par.arcsep*y2d/2,par.arcwidth*y2d)...
            +gaussian(mlat-par.drift*y2d*it-mlatctr, par.arcsep*y2d/2,par.arcwidth*y2d)...
            );%.*gaussian(mlon,mlonmean,par.mlonsig);
        E0it = par.E0itpk.*(gaussian(mlat-par.drift*y2d*it-mlatctr,-par.arcsep*y2d/2,par.arcwidth*y2d)...
            +gaussian(mlat-par.drift*y2d*it-mlatctr, par.arcsep*y2d/2,par.arcwidth*y2d)...
            );%.*gaussian(mlon,mlonmean,par.mlonsig);
        J = Jpk.*(tanh(50*par.flowwidth*y2d*(mlat-par.drift*y2d*it-mlatctr-(-5/2)*par.sheetwidth*y2d-offset))...
            -2*tanh(50*par.flowwidth*y2d*(mlat-par.drift*y2d*it-mlatctr-(-3/2)*par.sheetwidth*y2d-offset))...
            +2*tanh(50*par.flowwidth*y2d*(mlat-par.drift*y2d*it-mlatctr-(-1/2)*par.sheetwidth*y2d-offset))...
            -2*tanh(50*par.flowwidth*y2d*(mlat-par.drift*y2d*it-mlatctr-( 1/2)*par.sheetwidth*y2d-offset))...
            +tanh(50*par.flowwidth*y2d*(mlat-par.drift*y2d*it-mlatctr-( 3/2)*par.sheetwidth*y2d-offset))...
            )/2;%.*gaussian(mlon,mlonmean,par.mlonsig);
        J = J.*(it>2);
        F = par.Fpk.*(gaussian(mlat-par.drift*y2d*it-mlatctr,(-5/2)*par.sheetwidth*y2d+offset,par.flowwidth*y2d)...
            -gaussian(mlat-par.drift*y2d*it-mlatctr,(-3/2)*par.sheetwidth*y2d+offset,par.flowwidth*y2d)...
            +gaussian(mlat-par.drift*y2d*it-mlatctr,(-1/2)*par.sheetwidth*y2d+offset,par.flowwidth*y2d)...
            -gaussian(mlat-par.drift*y2d*it-mlatctr,( 1/2)*par.sheetwidth*y2d+offset,par.flowwidth*y2d)...
            +gaussian(mlat-par.drift*y2d*it-mlatctr,( 3/2)*par.sheetwidth*y2d+offset,par.flowwidth*y2d)...
            );%.*gaussian(mlon,mlonmean,par.mlonsig);
        U = F./sqrt(1+dydx.^2);
        V = F.*dydx./sqrt(1+dydx.^2);
        function [y] = gaussian(x,pos,fwhm)
            y = exp(-(x-pos).^2./(0.36*(fwhm)^2));
        end
    end
    
    % Matt example: Angle
    function [Qit,E0it,J,U,V] = ANGLE_map(mlon,mlat,it,mlonmean,mlatmean,pars)
        par = pars.angle;
        displace = 10*par.mlatsig;
        mlatctr = mlatmean+displace*tanh((mlon-mlonmean)/(par.mlonsig));
        Qit  = par.Qitpk.*exp(-(mlon-mlonmean).^2/2/par.mlonsig^2).*exp(-(mlat-mlatctr-1.5*par.mlatsig).^2/2/par.mlatsig^2);
        E0it = par.E0itpk.*ones(size(mlon));
        J =   par.Jpk.*exp(-(mlon-mlonmean).^2/2/par.mlonsig^2).*exp(-(mlat-mlatctr-1.5*par.mlatsig).^2/2/par.mlatsig^2);
        J = J-par.Jpk.*exp(-(mlon-mlonmean).^2/2/par.mlonsig^2).*exp(-(mlat-mlatctr+1.5*par.mlatsig).^2/2/par.mlatsig^2);
        J = J.*(it>2);
        U = zeros(size(mlon));
        V = zeros(size(mlon));
    end
    
    % Map input
    function [Qit,E0it,J,U,V] = INTRP_map(mlon,mlat,it,map)
        Qit  = interp3(map.MLAT,map.MLON,map.IT,map.Qit ,mlat,mlon,it);
        E0it = interp3(map.MLAT,map.MLON,map.IT,map.E0it,mlat,mlon,it);
        J    = interp3(map.MLAT,map.MLON,map.IT,map.J   ,mlat,mlon,it);
        U    = interp3(map.MLAT,map.MLON,map.IT,map.U   ,mlat,mlon,it);
        V    = interp3(map.MLAT,map.MLON,map.IT,map.V   ,mlat,mlon,it);
    end
    
    % Magnetospheric Accelerated Energetic inverted-V Electrons (MAEVE)
    function [Qit,E0it,J,U,V] = MAEVE_map(x2,x3,it,pars)
        p = pars.maeve;
        x2 = x2 - p.driftE*it;
        x3 = x3 - p.driftN*it;
        c = (p.ctr_spn/2)*tanh(2*p.ctr_slp*x2/p.ctr_spn);
        dcdx = p.ctr_slp*sech(2*p.ctr_slp*x2/p.ctr_spn).^2;
    %     s = sqrt(1+dcdx.^2);
        s = 1;
        b = bar(x2,p.bar_pos+p.bar_vel*it,p.bar_frc,p.bar_gsl); % loading bar
        d = (2-p.dim_frc*(1-tanh(2*(it-p.dim_del)/p.dim_tim)))/2; % dimming
        J_amp = p.K_amp/p.J_wth;
        Qit = (p.Q_amp_h-p.Q_amp_l-p.Q_floor)*d.*b.*...
              sheet(x3,c+p.Q_wth_l/2+p.Q_off_l+p.Q_off_h*p.Q_wth_l/2,p.Q_wth_h*s,p.Q_gsl_h)...
             +(p.Q_amp_l-p.Q_floor)*...
              sheet(x3,c+p.Q_wth_l/2+p.Q_off_l,p.Q_wth_l*s,p.Q_gsl_l)...
             +p.Q_floor;
        E0it = (p.E_amp_h-p.E_amp_l-p.E_floor)*d.*b.*...
               sheet(x3,c+p.Q_wth_l/2+p.Q_off_l+p.Q_off_h*p.Q_wth_l/2,p.E_wth_h*s,p.E_gsl_h)...
              +(p.E_amp_l-p.E_floor)*...
               sheet(x3,c+p.Q_wth_l/2+p.Q_off_l,p.E_wth_l*s,p.E_gsl_l)... % h and l in same pos as Q
              +p.E_floor;
        J = J_amp*(...
             sheet(x3,c+p.J_wth/2,p.J_wth*s,p.J_gsl)...
            -sheet(x3,c-p.J_wth/2,p.J_wth*s,p.J_gsl)...
            );
        F = p.F_amp*(...
             sheet(x3,c-p.J_wth,p.F_wth*s,p.F_gsl)...
            -sheet(x3,c          ,p.F_wth*s,p.F_gsl)...
            +sheet(x3,c+p.J_wth,p.F_wth*s,p.F_gsl)...
            );
        U = F.*(1/sqrt(1+dcdx.^2));
        V = F.*(dcdx./sqrt(1+dcdx.^2));
    end

    % Basic Functions
    function [v] = sheet(x3,pos,wdth,gsl)
        v = (tanh(2*(x3-pos+wdth/2)./(gsl*wdth))-tanh((x3-pos-wdth/2)./(gsl*wdth)))/2;
    end
    function [v] = bar(x2,pos,frac,gsl)
        v = (2-frac*(1-tanh(2*(x2-pos)/gsl)))/2;
    end
    
%% Other functions
    % user query function
    function [v] = usrq(q,defv)
        vi = input(q);
        if isempty(vi)
            v = defv;
        else
            v = vi;
        end
    end
end