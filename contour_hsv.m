function [] = contour_hsv(files,altref,start,dur,stop,plotmode)

arguments
    files (1,1) string %This is the directory with the data
    altref (1,1) double = 175 %Altitude for center panel
    start (1,1) double = 0 %start time
    dur (1,1) double = 10 %Cadence
    stop  (1,1) double = 300 %timestep at which to halt
    plotmode (1,1) string = 'both' %Plotting mode (defaults to running auto and standard)
end
files=char(files);
direc=files;

%files=dir([direc,'/input_arrays/*_particles/*.mat']);
files=dir([direc,'/*_particles/*.mat']);
load([files.folder,'/',files.name],'MLAT','MLON','Qit');
MLATp=MLAT(:,:,end);
MLONp=MLON(:,:,end);

sized=14;    % label size
titlesized=16;
files=direc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(plotmode,'auto')
    direc = files;
    folder = 'contour_HSV';
    if ~exist([direc,'/',folder],'dir')
        mkdir(direc,folder);
    end
    auto=1;
elseif strcmp(plotmode,'standard')
    direc = files;
    folder = 'contour_HSV';
    if ~exist([direc,'/',folder],'dir')
        mkdir(direc,folder);
    end
    auto=0;
elseif strcmp(plotmode,'BOTH_Auto')
    direc = files;
    folder = 'contour_auto';
    if ~exist([direc,'/',folder],'dir')
        mkdir(direc,folder);
    end
    auto=1;
elseif strcmp(plotmode,'BOTH_Standard')
    direc = files;
    folder = 'contour_standard';
    if ~exist([direc,'/',folder],'dir')
        mkdir(direc,folder);
    end
    auto=0;
end

if strcmp(plotmode,'both')
    fprintf('Running Autocolor \n');
    contour_hsv(files,altref,start,dur,stop,'BOTH_Auto')
    fprintf('Running Standard \n');
    contour_hsv(files,altref,start,dur,stop,'BOTH_Standard')
    return
end

%}
plotmode='auto';
if (~exist('ymd0','var'))
    confdat=gemini3d.read.config(files);
    xg=gemini3d.read.grid(files);
    fprintf('Input config file loaded.\n');
end
%}
ymd0=confdat.ymd;
UTsec0=confdat.UTsec0;
tdur=confdat.tdur;
dtout=confdat.dtout;
flagoutput=confdat.flagoutput;

dt = dtout;
tmin=0;
tmax=tdur;
times=start:dur:stop;
lt=length(times);

ymd=ymd0;
UTsec=UTsec0+times;     %time given in file is the seconds from beginning of hour
UThrs=UTsec/3600;
expdate=cat(2,repmat(ymd,[lt,1]),UThrs(:),zeros(lt,1),zeros(lt,1));
hours=mod(expdate(:,4),1);
expdate(:,4)=expdate(:,4)-hours;
expdate(:,5)=hours.*60;
minute=mod(expdate(:,5),1);
expdate(:,5)=expdate(:,5)-minute;
expdate(:,6)=minute.*60;
dates=datetime(int32(expdate));
names=gemini3d.datelab(dates);
if strcmp(confdat.file_format,"raw")
    confdat.file_format="dat";
end
names=names+"."+confdat.file_format;
%{
if contains(names(1),string(confdat.UTsec0))
    firststr=char(names(1));
    firststr(end-1-length(char(confdat.file_format)))='1';
    names(1)=string(firststr);
end
%}
names=char(names);

for run_num=1
    if strcmp(plotmode,'BOTH_Auto')
        plotmode = 'auto';
    elseif strcmp(plotmode,'BOTH_stndard')
        plotmode = 'standard';
    end
    direc=files;
    clear SIGP_all SIGH_all UTsec_all
end
it=0;
for UTsec=start:dur:stop
    it=it+1;
    %for alt = 100:100:1000
    temp=names(:,:,it);
    savename=temp(1:end-length(char(confdat.file_format))-8);
    if exist([direc,'/',folder,'/',savename,'_hsv.png'],'file')
        disp('File located, moving to next iteration.');
        continue
    end
    data = gemini3d.read.frame([files,'/',names(:,:,it)]);
    disp(['Processing file ',names(:,:,it)]);
    %filename=data.filename;
    %time=data.time;
    %ns=data.ns;
    %vs1=data.vs1;
    %Ts=data.Ts;
    J1=data.J1;
    J2=data.J2;
    J3=data.J3;
    v2=data.v2;
    v3=data.v3;
    Phitop=data.Phitop;
    %lxs=data.lxs;
    ne=data.ne;
    Ti=data.Ti;
    Te=data.Te;
    v1=data.v1;
    lx1=xg.lx(1); 
    lx2=xg.lx(2); 
    lx3=xg.lx(3);
    
    
    
    % Conductivity flag
    doconductivity=0;
    
    % Plotflag
    plotflag=1;
    
    % Altitude range
    altrangehere=[0 500]; %%%% NOTE THAT THIS DOES NOT CHANGE THE HSV PANEL...
    ix2=floor(lx2/2);
    ix3=floor(lx3/2);
    
    % Plot colorbar limits
    v2range=[-10000 10000];   % N-S
    saturated=600;
    %v1range=[-10 10];
    nerange=[9 12];
    %Tirange=[100 1800];
    %Terange=[100 1800];
    JPrange=[-20 20];
    %JHrange=[-9 -4];
    Jparrange=[-20 20];
    Phitoprange=[-5000 5000];
    
    if doconductivity==1
        IntConrangeH=[-110 0];
        IntConrangeP=[0 25];
        sigfac=[0 3];
        sighal=[-4.5e-3 0];
        sigped=[0 8e-4];
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%% Reconstruct conductivities if flag is set
    if doconductivity==1
        dat=gemini3d.read.config(direc);
        [sigP,sigH,sig0,SIGP,SIGH]=gemscr.postprocess.conductivity_reconstruct(xg,ymd,UTsec,activ,ne,Ti,Te,v1);
        SIGP_all(:,:,it)=SIGP;
        SIGH_all(:,:,it)=abs(SIGH); % corrects for Z's sigm convention
        UTsec_all(it)=UTsec;
        clear sigP sigH sig0 SIGP SIGH
        it=it+1;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Calculate JP JH and Jpar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % RESOLVE CURRENTS INTO PEDERSEN AND HALL COMPONENTS
    % lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
    vperp2=squeeze(v2(floor(lx1/2),:,:));
    vperp2=repmat(vperp2,[1,1,lx1]);
    vperp3=squeeze(v3(floor(lx1/2),:,:));
    vperp3=repmat(vperp3,[1,1,lx1]);
    v=cat(4, vperp2, vperp3, zeros(lx2,lx3,lx1));           %create a vector for the ExB drift in the curvilinear basis (i.e. e1,e2,e3 basis vectors), permuted as 231 so that the parallel direction is the third dimension
    magvperp=sqrt(sum(v.^2,4));
    evperp=v./repmat(magvperp,[1,1,1,3]);                %unit vector for ExB drift (curv. basis)
    e1curv=cat(4,zeros(lx2,lx3,lx1),zeros(lx2,lx3,lx1),ones(lx2,lx3,lx1));                         %unit vector (x1-direction) in the curvilinear basis
    
    B=cat(4,zeros(lx2,lx3,lx1),zeros(lx2,lx3,lx1),ones(lx2,lx3,lx1).*permute(xg.Bmag,[2,3,1]));
    E=cross(-1*v,B,4);                                                                             %electric field
    magE=sqrt(sum(E.^2,4));
    eE=E./repmat(magE,[1,1,1,3]);                                                                  %unit vector in the direction of the electric field
    
    J=cat(4,permute(J2,[2,3,1]),permute(J3,[2,3,1]),permute(J1,[2,3,1]));     %current density vector, permuted as 231
    JH=dot(J,-1*evperp,4);                                                    %projection of current in -evperp direction
    JHvec=-1*evperp.*repmat(JH,[1,1,1,3]);
    JP=dot(J,eE,4);                                                           %project of current in eE unit vector direction (direction of electric field)
    JPvec=eE.*repmat(JP,[1,1,1,3]);
    Jfac=cat(4,zeros(lx2,lx3,lx1),zeros(lx2,lx3,lx1),permute(J1,[2,3,1]));    %field aligned current vector
    
    % PERMUTE THE P,H,FAC CURRENT BACK TO WHAT THE MODEL NORMALLY USES
    JH=permute(JHvec,[3,1,2,4]);
    JP=permute(JPvec,[3,1,2,4]);
    Jfac2=permute(Jfac,[3,1,2,4]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% All the things to be plotted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % JPplotting
    JPplot=sqrt(sum(JPvec.^2,4));
    JPplot=ipermute(JPplot,[2,3,1]);
    parmP=log10(JPplot);
    
    % JHplotting
    JHplot=sqrt(sum(JHvec.^2,4));
    JHplot=ipermute(JHplot,[2,3,1]);
    parmH=log10(JHplot);
    
    % Jfac plotting
    parm=(Jfac2(:,:,:,3));
    
    % Velocity plotting
    parmV1=v1;
    parmV2=v2;
    parmV3=v3;
    
    % Electron density plotting
    parmne=log10(ne);
    
    % Temperature plotting
    parmTi=Ti;
    parmTe=Te;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Setting up the plotting grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %SIZE OF SIMULATION -- adjusts for ghost cells
    inds1=3:lx1+2;  % alt
    inds2=3:lx2+2;  % E
    inds3=3:lx3+2;  % N
    Re=6370e3;
    
    %SIZE OF PLOT GRID THAT WE ARE INTERPOLATING ONTO
    meantheta=mean(xg.theta(:));
    y=-1*(xg.theta-meantheta);   %this is a mag colat. coordinate and is only used for defining grid in linspaces below, runs backward from north distance, hence the negative sign
    x=xg.x2(inds2)/Re/sin(meantheta);
    z=xg.alt/1e3;
    lxp=200;
    lyp=200;
    lzp=600;
    minx=min(x(:));
    maxx=max(x(:));
    miny=min(y(:));
    maxy=max(y(:));
    minz=min(z(:));
    maxz=max(z(:));
    xp=linspace(minx,maxx,lxp);     %eastward distance (rads.)
    yp=linspace(miny,maxy,lyp);     %should be interpreted as northward distance (in rads.).  Irrespective of ordering of xg.theta, this will be monotonic increasing!!!
    zp=linspace(minz,maxz,lzp)';     %altitude (kilometers)
    
    MLAT=90-squeeze(xg.theta(1,:,:))*180/pi;
    MLON=squeeze(xg.phi(1,:,:))*180/pi;
    mlat=linspace(min(min(MLAT)),max(max(MLAT)),lyp);      % this way it should auto scale if we change the resolution for publication or just 'cause
    mlon=linspace(min(min(MLON)),max(max(MLON)),lxp);      % i.e. to 200 x 200 or whatever lxp x lyp is...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Slicing up the data for plotting
    %%%%%% I removed the unused slices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %NORTH/EAST SLICE grid space
    zp2=[altref-10,altref,altref+10];
    lzp2=numel(zp2);
    [X2,Y2,Z2]=meshgrid(xp,yp,zp2*1e3);       %lat./lon. meshgrid, need 3D since and altitude slice cuts through all 3 dipole dimensions
    x1plot=Z2(:);      %upward distance
    x2plot=X2(:)*Re*sin(meantheta);     %eastward distance - needs to be fixed to be theta-dependent (e.g. sin theta)
    x3plot=Y2(:)*Re;   %northward distance;
    x3interp=xg.x3(inds3);             %this is northward distance - again backwards from yp
    x3interp=x3interp(:);   %interp doesn't like it unless this is a column vector
    x1plot=single(x1plot);
    
    parmtmp=permute(parmne,[3,2,1]);     %so north dist, east dist., alt.
    parmp2ne=interp3(xg.x2(inds2),x3interp,xg.x1(inds1),parmtmp,x2plot,x3plot,x1plot);
    parmp2ne=reshape(parmp2ne,[lyp,lxp,lzp2]);    %slice expects the first dim. to be "y"
    
    
parmtmp=permute(v2,[3,2,1]);     %so north dist, east dist., alt.
parmp3V2=interp3(xg.x2(inds2),x3interp,xg.x1(inds1),parmtmp,x2plot,x3plot,x1plot);
parmp3V2=reshape(parmp3V2,[lyp,lxp,lzp2]);    %slice expects the first dim. to be "y"

parmtmp=permute(v3,[3,2,1]);     %so north dist, east dist., alt.
parmp3V3=interp3(xg.x2(inds2),x3interp,xg.x1(inds1),parmtmp,x2plot,x3plot,x1plot);
parmp3V3=reshape(parmp3V3,[lyp,lxp,lzp2]);    %slice expects the first dim. to be "y"

        %EDIT - Ted McManus
        if UTsec==0
            tind=1;
        else
            tind=UTsec;
        end
    parmQBG=Qit(:,:,tind); %Pick the right timestep
    sarr = size(Qit);
    if (any(sarr(1:2)~=xg.lx(2:3)))
        parmQBG=zeros(xg.lx(2:3));
    end
    full = permute(repmat(parmQBG, [1 1 xg.lx(1)]), [3 1 2]); %GEMINI loves this garbage
    
    parmtmp=permute(full,[3,2,1]);     %so north dist, east dist., alt.
    parmp2QBG=interp3(xg.x2(inds2),x3interp,xg.x1(inds1),parmtmp,x2plot,x3plot,x1plot);
    parmp2QBG=reshape(parmp2QBG,[lyp,lxp,lzp2]); 
    
    newfull = permute(repmat(Phitop, [1 1 xg.lx(1)]), [3 1 2]); %GEMINI loves this garbage
    parmtmp=permute(newfull,[3,2,1]);     %so north dist, east dist., alt.
    parmp2Phitop=interp3(xg.x2(inds2),x3interp,xg.x1(inds1),parmtmp,x2plot,x3plot,x1plot);
    parmp2Phitop=reshape(parmp2Phitop,[lyp,lxp,lzp2]);    %slice expects the first dim. to be "y"
    
    %{
    parmtmp=permute(v3,[3,2,1]);     %so north dist, east dist., alt.
    parmp2v3=interp3(xg.x2(inds2),x3interp,xg.x1(inds1),parmtmp,x2plot,x3plot,x1plot);
    parmp2v3=reshape(parmp2v3,[lyp,lxp,lzp2]);    %slice expects the first dim. to be "y"
    
    parmtmp=permute(v2,[3,2,1]);     %so north dist, east dist., alt.
    parmp2v2=interp3(xg.x2(inds2),x3interp,xg.x1(inds1),parmtmp,x2plot,x3plot,x1plot);
    parmp2v2=reshape(parmp2v2,[lyp,lxp,lzp2]);    %slice expects the first dim. to be "y"
    %}
    parmtmp=permute(parm,[3,2,1]);     %so north dist, east dist., alt.
    parmp2Jfac=interp3(xg.x2(inds2),x3interp,xg.x1(inds1),parmtmp,x2plot,x3plot,x1plot);
    parmp2Jfac=reshape(parmp2Jfac,[lyp,lxp,lzp2]); 
  
    %NORTH/UP SLICE grid space
    [Y3,Z3]=meshgrid(yp,zp*1e3);
    x1plot=Z3(:);      %upward distance
    x3plot=Y3(:)*Re;   %northward distance;
    
    % interpolating the data onto the slice
    %{
    parmtmp=squeeze(parmP(:,ix2,:));     %so north dist, east dist., alt.
    parmp3P=interp2(xg.x3(inds3),xg.x1(inds1),parmtmp,x3plot,x1plot);
    parmp3P=reshape(parmp3P,[lzp,lyp]);    %slice expects the first dim. to be "y"
    
    parmtmp=squeeze(parmH(:,ix2,:));     %so north dist, east dist., alt.
    parmp3H=interp2(xg.x3(inds3),xg.x1(inds1),parmtmp,x3plot,x1plot);
    parmp3H=reshape(parmp3H,[lzp,lyp]);    %slice expects the first dim. to be "y"
    %}
    parmtmp=squeeze(parm(:,ix2,:));     % FAC    %so north dist, east dist., alt.
    parmp3=interp2(xg.x3(inds3),xg.x1(inds1),parmtmp,x3plot,x1plot);
    parmp3=reshape(parmp3,[lzp,lyp]);    %slice expects the first dim. to be "y"
    
    parmtmp=squeeze(parmV2(:,ix2,:));     %so north dist, east dist., alt.
    parmp3V2_a=interp2(xg.x3(inds3),xg.x1(inds1),parmtmp,x3plot,x1plot);
    parmp3V2_a=reshape(parmp3V2_a,[lzp,lyp]);    %slice expects the first dim. to be "y"
    
    parmtmp=squeeze(parmV3(:,ix2,:));     %so north dist, east dist., alt.
    parmp3V3_a=interp2(xg.x3(inds3),xg.x1(inds1),parmtmp,x3plot,x1plot);
    parmp3V3_a=reshape(parmp3V3_a,[lzp,lyp]);    %slice expects the first dim. to be "y"
    
    parmtmp=squeeze(parmne(:,ix2,:));     %so north dist, east dist., alt.
    parmp3ne=interp2(xg.x3(inds3),xg.x1(inds1),parmtmp,x3plot,x1plot);
    parmp3ne=reshape(parmp3ne,[lzp,lyp]);    %slice expects the first dim. to be "y"
    %{
    parmtmp=squeeze(parmTe(:,ix2,:));     %so north dist, east dist., alt.
    parmp3Te=interp2(xg.x3(inds3),xg.x1(inds1),parmtmp,x3plot,x1plot);
    parmp3Te=reshape(parmp3Te,[lzp,lyp]);    %slice expects the first dim. to be "y"
    %}
    %Following convention

    parmtmp=squeeze(full(:,ix2,:));     %so north dist, east dist., alt.
    parmQBG=interp2(xg.x3(inds3),xg.x1(inds1),parmtmp,x3plot,x1plot);
    parmQBG=reshape(parmQBG,[lzp,lyp]);    %slice expects the first dim. to be "y"
    
    % More fussy grid stuff
    % CONVERT ANGULAR COORDINATES TO MLAT,MLON
    yp2=yp*Re/1e3; %(km)
    [yp2,inds]=sort(yp2);
    xp2=xp*Re*sin(meantheta)/1e3;    %eastward ground distance (km)
    [xp2,inds2]=sort(xp2);
    
    
    % COMPUTE SOME BOUNDS FOR THE PLOTTING
    %{
    minxp=min(xp2(:));
    maxxp=max(xp2(:));
    minyp=min(yp2(:));
    maxyp=max(yp2(:));
    %}
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % For the HSV panel -- have to do this one first or else there is colormap contamination
    koutu = parmp3V2_a(1:53,:);                               % E-W velocity
    koutv = parmp3V3_a(1:53,:);                               % N-S velocity
    vectordirection = ((atan2(koutv,koutu))*57.3+180)/360.; % combine into direction
    vectormagnitude = sqrt(koutu.^2+koutv.^2)/(saturated);  % calculate magnitude
    vectormagnitude(find(vectormagnitude > 1)) = 1;         % adjust to saturate brightness
    map=colorcet('C2','shift',63/256);                      % define new colormap - west is blue
    
    if ~exist('colorWheel_friendly_simple_w.jpg')
        % Makes the jpg if needed, added map passing
        % colorCircle_friendly_simplified(map)   % shade to black - no markings
        colorCircle_friendly_simplified_w(map)   % shade to white - no markings
    end
    
    % Generate the hsv image
    robin=showangularim(vectordirection.*2*pi,colormap(map),'amp',vectormagnitude,'bw',1,'fig',0);
    
    % Define hsv panel axis numbers
    xwords=mlat;
    ywords=zp(1:53);
    hi2=imref2d(size(robin));
    hi2.XWorldLimits = [xwords(1) xwords(end)];
    hi2.YWorldLimits = [ywords(1) ywords(end)];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Plot all the things!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(1);
    set(gcf,'PaperPosition',[0,0,22,8],'PaperUnits','inches')
    %% VERTICAL FLOW HSV
    % Finally plot the hsv panel
    cat2=subplot(3,3,7);
    himage=imshow(robin,hi2); % the plotting part is here
    hold on
    [~,c]=contour(mlat,zp,parmQBG,1,'black');c.LineWidth=.5;
    hold off
    set(gca,'FontSize',sized)
    axis xy
    axis normal               % makes it the same size as the rest of the panels
    colormap(cat2,map)        % this sets the colorbar to be hsv colors
    clb2=colorbar;
    ylabel('Altitude (km)','FontSize',sized)
    xlabel('Mag. Lat.','FontSize',sized)
    set(clb2,'ytick',[0 0.25 .5 0.75 1])  % manually label the directions on the colorbar
    set(clb2,'yticklabel',{'W','S','E','N','W'})
    ylabel(clb2,['Saturated at ',num2str(round(saturated,0)),' m/s'],'FontSize',sized)
    set(gca,'TickDir','in')
    caxis([0,1]);
    
    % Add the little colorwheel
    axes('pos',[0 0 .1 .1])
    % imshow('colorWheel_friendly_simple.jpg')  % No labels
    % % imshow('colorWheel_friendly.jpg')       % With labels
    imshow('colorWheel_friendly_simple_w.jpg')  % No labels - shades to white
    % imshow('colorWheel_friendly_w.jpg')       % With labels - shades to white
    
% For the HSV panel -- have to do this one first or else there is colormap contamination
koutu = parmp3V2(:,:,2);                               % E-W velocity, the "2" index number is taking the slice at altref
koutv = parmp3V3(:,:,2);                               % N-S velocity
vectordirection = ((atan2(koutv,koutu))*57.3+180)/360.; % combine into direction
vectormagnitude = sqrt(koutu.^2+koutv.^2)/(saturated);  % calculate magnitude
vectormagnitude(find(vectormagnitude > 1)) = 1;         % adjust to saturate brightness
map=colorcet('C2','shift',63/256);                      % define new colormap - west is blue

% Generate the hsv image
robin=showangularim(vectordirection.*2*pi,colormap(map),'amp',vectormagnitude,'bw',1,'fig',0);

% Define hsv panel axis numbers - uncomment this section to make pretty plot
xwords=mlon;
ywords=mlat;
hi2=imref2d(size(robin));   % I'm missing toolbox that contains imref2d so I have to comment this section out
hi2.XWorldLimits = [xwords(1) xwords(end)];
hi2.YWorldLimits = [ywords(1) ywords(end)];

%% HORIZONTAL FLOW HSV
cat2=subplot(3,3,[1 4]);
% set(gcf,'PaperPosition',[0,0,6,6],'PaperUnits','inches')
himage=imshow(robin,hi2); % the plotting part is here - uncomment for pretty plot
    hold on
    contour(mlon,mlat,parmp2QBG(:,:,2),1,'red','LineWidth',.5);
    hold off
% himage=imshow(robin); % the plotting part is here
set(himage,'alphadata',~isnan(robin(:,:,2)));
set(gca,'FontSize',sized)
axis xy
axis normal               % makes it the same size as the rest of the panels
colormap(cat2,map)        % this sets the colorbar to be hsv colors
clb2=colorbar;
ylabel('Mag. Lat.','FontSize',sized)
xlabel('Mag. Lon.','FontSize',sized)
set(clb2,'ytick',[0 0.25 .5 0.75 1])  % manually label the directions on the colorbar
set(clb2,'yticklabel',{'W','S','E','N','W'})
ylabel(clb2,['Saturated at ',num2str(round(saturated,0)),' m/s'],'FontSize',sized)
title('Flow (m/s)');
caxis([0,1]);
set(gca,'TickDir','in')
    %{
    cat1=subplot(3,3,4);
    h=imagesc(mlat,zp,parmp3V2);
    hold on
    [~,c]=contour(mlat,zp,parmQBG,3,'black');c.LineWidth=.5;
    hold off
    set(h,'alphadata',~isnan(parmp3V2));
    set(gca,'FontSize',sized)
    axis xy
    clb=colorbar;
    colormap(cat1,'parula')
    ylim(altrangehere)
    xlabel('Mag. Lat.','FontSize',sized)
    ylabel('Altitude (km)','FontSize',sized);
    ylabel(clb,'E-W Ion Velocity (m/s)','FontSize',sized)
    if auto==1
        caxis('auto')
    else
        caxis(v3range)
    end
    %}
    
    cat1=subplot(3,3,[2 5]);
    plot_array=parmp2Phitop(:,:,2);
    h=imagesc(mlon,mlat,plot_array);
    hold on
    contour(mlon,mlat,parmp2QBG(:,:,2),1,'red','LineWidth',.5);
    hold on
    set(h,'alphadata',~isnan(parmp2Phitop(:,:,2)));
    set(gca,'FontSize',sized)
    hold on
    plot([(max(mlon)+min(mlon))/2,(max(mlon)+min(mlon))/2],[min(mlat),max(mlat)],'r--','LineWidth',1); % red line to show slice location
    hold off
    axis xy
    clb=colorbar;
    colormap(cat1,'parula')
    ylabel('Mag. Lat.','FontSize',sized)
    xlabel('Mag. Lon.','FontSize',sized);
    ylabel(clb,['Potential (Volts)'],'FontSize',sized)
    title('Potential (Volts)')
    high=max(plot_array(:));
    low=min(plot_array(:));
    if high>abs(low)
        autoscale=[-high,high];
    else
        autoscale=[-abs(low),abs(low)];
    end
    if auto==1
        caxis(autoscale);
    else
        caxis(Phitoprange);
    end
   

    %% JFAC Horizontal
    cat1=subplot(3,3,[3 6]);
    plotarr=parmp2Jfac(:,:,2).*10^6;
    h=imagesc(mlon,mlat,-plotarr);
    hold on
    contour(mlon,mlat,parmp2QBG(:,:,2),1,'LineColor','black','LineWidth',.5);
    hold off
    set(h,'alphadata',~isnan(parmp2Jfac(:,:,2)))
    set(gca,'FontSize',sized)
    axis xy
    clb=colorbar;
    colormap(cat1,colorcet('D1A'));
    %% 
    ylabel('Mag. Lat.','FontSize',sized)
    xlabel('Mag. Lon.','FontSize',sized);
    ylabel(clb,'J_{fac} (uA/m^2)');
    title(['Field-Aligned Current (uA/m^2)']);
    
    high=max(plotarr(:));
    low=min(plotarr(:));
    if high>abs(low)
        autoscale=[-high,high];
    elseif high<abs(low)
        autoscale=[-abs(low),abs(low)];
    end
    if auto==1
        caxis(autoscale);
    else
        caxis(JPrange);
    end
    %% NE
    cat1=subplot(3,3,8);
    h=imagesc(mlat,zp,parmp3ne);
    hold on
    [~,c]=contour(mlat,zp,parmQBG,1,'black');c.LineWidth=.5;
    hold off
    set(h,'alphadata',~isnan(parmp3ne));
    set(gca,'FontSize',sized)
    axis xy
    clb=colorbar;
    colormap(cat1,'parula')
    ylim(altrangehere)
    xlabel('Mag. Lat.','FontSize',sized)
    ylabel('Altitude (km)','FontSize',sized);
    ylabel(clb,'Electron Density (m^{-3})','FontSize',sized)
    if auto==1
        caxis([min(parmp3ne(:)),max(parmp3ne(:))]);
    else
        caxis(nerange);
    end
    
    %% JFAC
    newmap=colorcet('D1A');       % this one is my favorite (not anymore hehehehe)
    plotarr=-parmp3*10^6;
    cat3=subplot(3,3,9);
    h=imagesc(mlat,zp,plotarr);   % forced the sign convention
    hold on
    [~,c]=contour(mlat,zp,parmQBG,1,'black');c.LineWidth=.5;
    hold off
    set(h,'alphadata',~isnan(parmp3));
    set(gca,'FontSize',sized)
    axis xy
    clb=colorbar;
    colormap(cat3,newmap)
    ylim(altrangehere)
    xlabel('Mag. Lat.','FontSize',sized)
    ylabel('Altitude (km)','FontSize',sized);
    ylabel(clb,'J_{fac} (uA/m^2)','FontSize',sized)
    
    high=max(plotarr(:));
    low=min(plotarr(:));
    if high>abs(low)
        autoscale=[-high,high];
    elseif high<abs(low)
        autoscale=[-abs(low),abs(low)];
    end
    if auto==1
        caxis(autoscale);
    else
        caxis(Jparrange);
    end
    
    %% EDIT: Overall title that shows time of slices and slice altitude
    t = datenum(ymd(1), ymd(2), ymd(3), UThrs, 0, 0);
    ttxt_prime = {[char(dates(it)),' UT at ',num2str(altref),'m alt']};
    sgtitle(ttxt_prime,'fontsize',titlesized,'fontweight','bold');
    
    
    % set(gcf, 'Position', [100 100 1400 500])
    
    if plotflag
        %if ~exist([direc,'/multipanel/'])
        %system(['mkdir ',direc,'/multipanel/']);
        %end
        if true
            temp=names(:,:,it);
            savename=temp(1:end-length(char(confdat.file_format))-8);
            saveas(gcf,[direc,'/',folder,'/',savename,'_hsv.png']);
        end
    end
    
    pause(0.3)
    
    if doconductivity==1
        save([direc,'/',direc(13:end),'_SIG.mat'],'SIGP_all','SIGH_all','UTsec_all','xg','-v7.3')
    end
    %end
end
 
end

function h = show(im, varargin)
        
        Octave = exist('OCTAVE_VERSION', 'builtin') == 5; % Are we running under Octave
        
        s = warning('query','all');                 % Record existing warning state.
        warn_state = onCleanup (@() warning(s));    % Restore warnings at the end
        warning('off');                             % Turn off warnings that might arise if image
        % has to be rescaled to fit on screen
        % Set defaults
        figNo = -1;  % Default indicates create new figure
        Title = '';
        clim = [];
        colourmap = gray(256);
        
        % Overwrite defaults with anything we recognize in the arguments
        %{
for n = 1:length(arg);
    if isscalar(arg{n})
        figNo = arg{n};
        
    elseif ischar(arg{n})
        Title = arg{n};
        
    elseif isnumeric(arg{n}) && all(size(arg{n}) == [1,2])
        clim = arg{n};
        
    elseif isnumeric(arg{n}) && size(arg{n},1) > 1 && size(arg{n},2) == 3
        colourmap = arg{n};
        
    else
        error('Unable to process arguments');
    end
end
        %}
        
        % Check case where im is an image filename rather than image data
        if ~isnumeric(im) && ~islogical(im)
            if isempty(Title)
                Title = im;        % Use file name for title
            end
            im = imread(im);
        elseif isempty(Title)
            Title = inputname(1);  % Use variable name of image data for title
        end
        
        sze = max(size(im));       % Maximum dimension of image
        
        if figNo > 0               % We have a valid figure number
            figure(figNo);         % Reuse or create a figure window with this number
            subplot('position',[0 0 1 1]); % Use the whole window
        elseif figNo == -1
            figNo = figure;        % Create new figure window
            subplot('position',[0 0 1 1]); % Use the whole window
        end
        
        if ndims(im) == 2          % Apply colour map
            if isempty(clim)
                imagesc(im);
            else
                imagesc(im, clim);
            end
            
            colormap(colourmap);
        else
            imshow(im(:,:,1:3));   % Display as RGB (ignore any alpha channel)
        end
        
        if figNo == 0              % Assume we are trying to do a subplot
            figNo = gcf;           % Get the current figure number
            %axis('image'); axis('off');
            title(Title);          % Use a title rather than rename the figure
        else
            %axis('image'); axis('off');
            set(figNo,'name', ['  ' Title])
            
            % If not running Octave and size of image > 500 plot image at 1:1
            % resolution. Otherwise we let imagesc use its default scaling.
            if ~Octave && sze > 500
                truesize(figNo);
            end
        end
        
        if nargout == 1
            h = figNo;
        end
    end
    function n = normalise(im, reqmean, reqvar)
        
        if ~(nargin == 1 | nargin == 3)
            error('No of arguments must be 1 or 3');
        end
        
        if nargin == 1   % Normalise 0 - 1
            if ndims(im) == 3         % Assume colour image
                hsv = rgb2hsv(im);
                v = hsv(:,:,3);
                v = v - min(v(:));    % Just normalise value component
                v = v/max(v(:));
                hsv(:,:,3) = v;
                n = hsv2rgb(hsv);
            else                      % Assume greyscale
                if ~isa(im,'double'), im = double(im); end
                n = im - min(im(:));
                n = n/max(n(:));
            end
            
        else  % Normalise to desired mean and variance
            
            if ndims(im) == 3         % colour image?
                error('cannot normalise colour image to desired mean and variance');
            end
            
            if ~isa(im,'double'), im = double(im); end
            im = im - mean(im(:));
            im = im/std(im(:));      % Zero mean, unit std dev
            
            n = reqmean + im*sqrt(reqvar);
        end
    end

function rgbim = showangularim(varargin)
        
        %     [ang, amp, map, cycle, bw, fig] = parseinputs2(varargin{:});
        p = inputParser;
        
        numericORlogical = @(x) isnumeric(x) || islogical(x);
        
        % The first arguments are the image of angular data and the colour map.
        addRequired(p, 'ang', @isnumeric);
        addRequired(p, 'map', @isnumeric);
        
        % Optional parameter-value pairs and their defaults
        addParameter(p, 'amp', [], @isnumeric);
        addParameter(p, 'cycle', 2*pi, @isnumeric);
        addParameter(p, 'bw', 0, numericORlogical);
        addParameter(p, 'fig', -1, @isnumeric);
        
        parse(p, varargin{:});
        
        ang = p.Results.ang;
        map = p.Results.map;
        amp = p.Results.amp;
        cycle = p.Results.cycle;
        bw = p.Results.bw;
        fig = p.Results.fig;
        if fig < 0,  fig = figure; end
        
        if ~isempty(amp) && ~all(size(amp)==size(ang))
            error('Amplitude data must be same size as angular data');
        end
        % Apply colour map to angular data.  Some care is needed with this.  Unlike
        % normal 'linear' data one cannot apply shifts and/or rescale values to
        % normalise them.  The raw angular data values have to be respected.
        
        ang = mod(ang, cycle);   % Ensure data values are within range 0 - cycle
        rgbim = applycolourmap(ang, map, [0 cycle]);
        
        if ~isempty(amp)  % Display image with rgb values modulated by amplitude
            
            amp = normalise(amp);  % Enforce amplitude  0 - 1
            
            if ~bw  % Modulate rgb values by amplitude fading to black
                for n = 1:3
                    rgbim(:,:,n) = rgbim(:,:,n).*amp;
                end
                
            else  % Modulate rgb values by amplitude fading to white
                for n = 1:3
                    rgbim(:,:,n) = 1 - (1 - rgbim(:,:,n)).*amp;
                end
            end
        end
        
        if fig
            show(rgbim,fig)
        end
        
        % If function was not called with any output arguments clear rgbim so that
        % it is not printed on the screen.
        if ~nargout
            clear('rgbim')
        end
        
        
    end
    function rgbim = applycolourmap(im, map, range)
        
        [ncolours,chan] = size(map);
        assert(chan == 3, 'Colourmap must have 3 columns');
        assert(ndims(im) == 2, 'Image must be single channel');
        if ~isa(im,'double'), im = double(im); end
        
        if ~exist('range', 'var') || isempty(range)
            range = [min(im(:)) max(im(:))];
        end
        assert(range(1) < range(2), 'range(1) must be less than range(2)');
        
        [rows,cols] = size(im);
        
        % Convert image values to integers that can be used to index into colourmap
        im = round( (im-range(1))/(range(2)-range(1)) * (ncolours-1) ) + 1;
        
        mask = isnan(im);
        
        im(mask) = 1;           % Set any Nan entries to 1 and
        im(im < 1) = 1;         % clamp out of range entries.
        im(im > ncolours) = ncolours;
        
        rgbim = zeros(rows,cols,3);
        
        rgbim(:,:,1) = ~mask.*reshape(map(im,1), rows, cols);
        rgbim(:,:,2) = ~mask.*reshape(map(im,2), rows, cols);
        rgbim(:,:,3) = ~mask.*reshape(map(im,3), rows, cols);
    end
    function colorCircle_friendly_simplified(map)
        
        hsvmap = rgb2hsv(map);
        
        %creating HSV vector of hues
        h = hsvmap(:,1);
        s = hsvmap(:,2);
        v = hsvmap(:,3);
        
        %converting the hsv vector of hues to RGB vector at saturation=1
        r=15;
        rgb = zeros(r,3,360);
        
        for j=1:360
            for l = 1:r %l identifies which circle you are in
                ind = floor(j*256/360);
                if ind == 0
                    ind = ind +1;
                end
                %         rgb(l,:,j) = hsv2rgb([h(ind), s(ind), (1-(l/r))*v(ind)]);
                rgb(l,:,j) = hsv2rgb([h(ind), (1-(l/r))*s(ind), v(ind)]);
            end
        end
        
        %plotting in polar scatter plot for a few different radii
        fig = figure;
        pax = polaraxes;
        
        % pax.ThetaZeroLocation='left';%change 0 degrees to be at the top
        % %pax.ThetaDir='clockwise';%change the direction of increasing degrees to be clockwise
        pax.GridColor = 'none';%getting rid of internal gridlines
        pax.RColor = 'none';%getting rid of radial labels (not relevant)
        thetaticks([0 90 180 270])
        thetaticklabels({'','','',''})
        set(gcf,'PaperPosition',[0,0,2,2],'PaperUnits','inches')
        set(gca,'FontSize',14)
        pax.ThetaZeroLocation='left';
        pax.ThetaTick=[0 45 90 135 180 225 270 315];  % manually label the directions on the colorbar
        pax.ThetaTickLabel=({'','','','','','','',''});
        ax1 = gca;
        ax1_pos = ax1.Position;
        ax1.RColor=[1 1 1];
        
        slice = squeeze(rgb(1,:,:));
        hold on
        for k=1:r
            slice = squeeze(rgb(k,:,:));
            polarscatter(linspace(0,2*pi,360), (1-k/r).*ones(1,360),20,slice.','filled');
        end
        thetas = linspace(0,2*pi);
        polarplot(thetas, 1*ones(size(thetas)), 'LineWidth',4, 'Color', 'White')
        ax2 = polaraxes('Position',ax1_pos,'Color','None');
        ax2.ThetaTick=[0 45 90 135 180 225 270 315];  % manually label the directions on the colorbar
        ax2.ThetaTickLabel=({'E','NE','N','NW','W','SW','S','SE'});
        ax2.RTick=[0 .25 .5 .75 1];
        ax2.RTickLabel=({'0%','25%','50%','75%','100%'}); % percentage of saturated value
        ax2.RColor=[1 1 1];
        ax2.RAxisLocation=67.5;
        set(ax2,'FontSize',12,'LineWidth',2)
        hold off
        
        saveas(fig, 'colorWheel_friendly_w.jpg')
        clf
        
        %plotting in polar scatter plot for a few different radii
        pax = polaraxes;
        
        % pax.ThetaZeroLocation='left';%change 0 degrees to be at the top
        % %pax.ThetaDir='clockwise';%change the direction of increasing degrees to be clockwise
        pax.GridColor = 'none';%getting rid of internal gridlines
        pax.RColor = 'none';%getting rid of radial labels (not relevant)
        thetaticks([0 90 180 270])
        thetaticklabels({'','','',''})
        set(gcf,'PaperPosition',[0,0,2,2],'PaperUnits','inches')
        set(gca,'FontSize',14)
        pax.ThetaZeroLocation='left';
        pax.ThetaTick=[0 45 90 135 180 225 270 315];  % manually label the directions on the colorbar
        pax.ThetaTickLabel=({'','','','','','','',''});
        ax1 = gca;
        ax1_pos = ax1.Position;
        ax1.RColor=[1 1 1];
        
        slice = squeeze(rgb(1,:,:));
        hold on
        for k=1:r
            slice = squeeze(rgb(k,:,:));
            polarscatter(linspace(0,2*pi,360), (1-k/r).*ones(1,360),20,slice.','filled');
        end
        thetas = linspace(0,2*pi);
        polarplot(thetas, 1*ones(size(thetas)), 'LineWidth',4, 'Color', 'White')
        ax2 = polaraxes('Position',ax1_pos,'Color','None');
        ax2.ThetaTick=[0 45 90 135 180 225 270 315];  % manually label the directions on the colorbar
        ax2.ThetaTickLabel=({'','','','','','','',''});
        ax2.RTick=[0 .25 .5 .75 1];
        ax2.RTickLabel=({'','','','',''}); % percentage of saturated value
        ax2.RColor=[1 1 1];
        ax2.RAxisLocation=67.5;
        set(ax2,'FontSize',12,'LineWidth',2)
        hold off
        
        saveas(fig, 'colorWheel_friendly_simple_w.jpg')
        clf
    end
    function colorCircle_friendly_simplified_w(map)
        dir=linspace(0,2*pi,360);
        mag=linspace(1,0,100);
        r=length(mag);
        [Dir,Mag]=meshgrid(dir,mag);
        
        robin=showangularim(Dir,colormap(map),'amp',Mag,'bw',1,'fig',0);
        robin=permute(robin,[1,3,2]);
        
        fig = figure;
        pax = polaraxes;
        pax.GridColor = 'none';%getting rid of internal gridlines
        pax.RColor = 'none';%getting rid of radial labels (not relevant)
        thetaticks([0 90 180 270])
        thetaticklabels({'','','',''})
        set(gcf,'PaperPosition',[0,0,2,2],'PaperUnits','inches')
        set(gca,'FontSize',14)
        pax.ThetaZeroLocation='left';
        pax.ThetaTick=[0 45 90 135 180 225 270 315];  % manually label the directions on the colorbar
        pax.ThetaTickLabel=({'','','','','','','',''});
        ax1 = gca;
        ax1_pos = ax1.Position;
        ax1.RColor=[1 1 1];
        
        slice = squeeze(robin(1,:,:));
        hold on
        for k=1:r
            slice = squeeze(robin(k,:,:));
            polarscatter(linspace(0,2*pi,360),(1-k/r).*ones(1,360),20,slice.','filled');
        end
        hold off
        
        print(fig,'colorWheel_friendly_simple_w.jpg','-djpeg')
        % clf
        
        
        %plotting in polar scatter plot for a few different radii
        fig = figure;
        pax = polaraxes;
        
        pax.GridColor = 'none';%getting rid of internal gridlines
        pax.RColor = 'none';%getting rid of radial labels (not relevant)
        thetaticks([0 90 180 270])
        thetaticklabels({'','','',''})
        set(gcf,'PaperPosition',[0,0,2,2],'PaperUnits','inches')
        set(gca,'FontSize',14)
        pax.ThetaZeroLocation='left';
        pax.ThetaTick=[0 45 90 135 180 225 270 315];  % manually label the directions on the colorbar
        pax.ThetaTickLabel=({'','','','','','','',''});
        ax1 = gca;
        ax1_pos = ax1.Position;
        ax1.RColor=[1 1 1];
        
        slice = squeeze(robin(1,:,:));
        hold on
        for k=1:r
            slice = squeeze(robin(k,:,:));
            polarscatter(linspace(0,2*pi,360), (1-k/r).*ones(1,360),20,slice.','filled');
        end
        thetas = linspace(0,2*pi);
        
        polarplot(thetas, 1*ones(size(thetas)), 'LineWidth',4, 'Color', 'White')
        ax2 = polaraxes('Position',ax1_pos,'Color','None');
        ax2.ThetaTick=[0 45 90 135 180 225 270 315];  % manually label the directions on the colorbar
        ax2.ThetaTickLabel=({'E','NE','N','NW','W','SW','S','SE'});
        ax2.RTick=[0 .25 .5 .75 1];
        ax2.RTickLabel=({'0%','25%','50%','75%','100%'}); % percentage of saturated value
        % ax2.RColor=[0.0 0 0];
        ax2.RAxisLocation=67.5;
        set(ax2,'FontSize',12,'LineWidth',2)
        hold off
        
        print 'colorWheel_friendly_w.jpg' '-djpeg'
        % clf
        
    end

