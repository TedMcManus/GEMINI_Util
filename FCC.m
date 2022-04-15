function [] = FCC(direc,start,dur,stop)

tic
%% SIMULATIONS LOCAITONS
flagplot=1;  %save flow, Jfac, conductance plot?
savedata=1;  %save the above data for students?
flagplotC=1; %save conductance only plot?
flagplotT=1; %save continuity terms plot?

loc=['~/../public_html/Gemini3D'];
direc=char(direc);
if strcmp(direc(1:2),'./')
    direc=direc(3:end);
end


thedate='20150201';
%thedate='20150620';

if flagplot==1
    outputdir=[direc,filesep,'FlowCurrentConductancePlots']; %folder for plots
    if ~exist(outputdir)
        system(['mkdir ',outputdir]); %make folder for plots if needed
    end
end
if flagplotC==1
    outputdirC=[direc,filesep,'ConductancePlots']; %folder for plots
    if ~exist(outputdirC)
        system(['mkdir ',outputdirC]); %make folder for plots if needed
    end
end
if flagplotT==1
    outputdirT=[direc,filesep,'ContinuityTerms']; %folder for plots
    if ~exist(outputdirT)
        system(['mkdir ',outputdirT]); %make folder for plots if needed
    end
end

[loc,RunName,~]=fileparts(direc);

xg = gemini3d.read.grid(direc);
cfg = gemini3d.read.config([direc]);
ymd0=cfg.ymd;
thedate=char(strrep(strjoin(string(ymd0)),' ',''));
if (length(char(string(ymd0(2)))))==1
    thedate=[thedate(1:4),'0',thedate(5:end)];
end

%% LOAD THE SIMULATION DATA CLOSEST TO THE REQUESTED TIME
for UTsec=start:dur:stop
    % for UTsec=35920:10:36150
    % for UTsec=36000  %use if only a specific timestep is needed

    time = datetime(ymd0) + seconds(UTsec);

    dat = gemini3d.read.frame([direc,filesep,thedate,'_',num2str(UTsec),'.000000.h5']);

    [sigP,sigH,sig0,SIGP,SIGH,incap,INCAP] = gemscr.postprocess.conductivity_reconstruct(time,dat,cfg, xg); %ram hog and sloooooow
    SIGH=-SIGH;                       %set to K sign
    Jpar=-squeeze(dat.J1(end,:,:));   %set to K sign
    vEW=squeeze(dat.v2(end,:,:));
    vNS=squeeze(dat.v3(end,:,:));
    v2=dat.v2;
    v3=dat.v3;
    try
        E0=h5read([direc,filesep,'STEVE_particles/',thedate,'_',num2str(UTsec),'.000000.h5'],'/E0p');
        Q=h5read([direc,filesep,'STEVE_particles/',thedate,'_',num2str(UTsec),'.000000.h5'],'/Qp');
    catch
        E0=h5read([direc,filesep,'OSSE_particles/',thedate,'_',num2str(UTsec),'.000000.h5'],'/E0p');
        Q=h5read([direc,filesep,'OSSE_particles/',thedate,'_',num2str(UTsec),'.000000.h5'],'/Qp');
    end
    MLAT=90-squeeze(xg.theta(1,:,:))*180/pi;
    MLON=squeeze(xg.phi(1,:,:))*180/pi;

    %% Setting up the plotting grid
    altref=500; % in km
    saturated=1000; % in m/s
    sized=12; % font size

    %SIZE OF SIMULATION -- auto adjusting for ghost cells
    inds1=3:size(v2,1)+2;  % alt
    inds2=3:size(v2,2)+2;  % E
    inds3=3:size(v2,3)+2;  % N
    Re=6370e3;

    %SIZE OF PLOT GRID THAT WE ARE INTERPOLATING ONTO
    meantheta=mean(xg.theta(:));
    y=-1*(xg.theta-meantheta);   %this is a mag colat. coordinate and is only used for defining grid in linspaces below, runs backward from north distance, hence the negative sign
    x=xg.x2(inds2)/Re/sin(meantheta);
    z=xg.alt/1e3;

    % resolution of new plotting grid
    lxp=150;
    lyp=150;
    lzp=180;

    % define new plotting grid
    minx=min(x(:));  maxx=max(x(:));
    miny=min(y(:));  maxy=max(y(:));
    minz=min(z(:));  maxz=max(z(:));
    xp=linspace(minx,maxx,lxp);     %eastward distance (rads.)
    yp=linspace(miny,maxy,lyp);     %should be interpreted as northward distance (in rads.).  Irrespective of ordering of xg.theta, this will be monotonic increasing!!!
    zp=linspace(minz,maxz,lzp)';     %altitude (kilometers)

    % nice axis numbers
    MLAT=90-squeeze(xg.theta(1,:,:))*180/pi;
    MLON=squeeze(xg.phi(1,:,:))*180/pi;
    mlat=linspace(min(min(MLAT)),max(max(MLAT)),lyp);      % this way it should auto scale if we change the resolution for publication or just 'cause
    mlon=linspace(min(min(MLON)),max(max(MLON)),lxp);      % i.e. to 200 x 200 or whatever lxp x lyp is...


    %% Slicing up the data for plotting
    %NORTH/EAST SLICE grid space
    zp2=[altref-10,altref,altref+10];
    lzp2=numel(zp2);
    [X2,Y2,Z2]=meshgrid(xp,yp,zp2*1e3);       %lat./lon. meshgrid, need 3D since and altitude slice cuts through all 3 dipole dimensions
    x1plot=cast(Z2(:),'single');                      %upward distance
    x2plot=X2(:)*Re*sin(meantheta);    %eastward distance - needs to be fixed to be theta-dependent (e.g. sin theta)
    x3plot=Y2(:)*Re;                   %northward distance;
    x3interp=xg.x3(inds3);                %this is northward distance - again backwards from yp
    x3interp=x3interp(:);              %interp doesn't like it unless this is a column vector

    parmtmp=permute(v2,[3,2,1]);     %so north dist, east dist., alt.
    parmp3V2=interp3(xg.x2(inds2),x3interp,xg.x1(inds1),parmtmp,x2plot,x3plot,x1plot);
    parmp3V2=reshape(parmp3V2,[lyp,lxp,lzp2]);    %slice expects the first dim. to be "y"

    parmtmp=permute(v3,[3,2,1]);     %so north dist, east dist., alt.
    parmp3V3=interp3(xg.x2(inds2),x3interp,xg.x1(inds1),parmtmp,x2plot,x3plot,x1plot);
    parmp3V3=reshape(parmp3V3,[lyp,lxp,lzp2]);    %slice expects the first dim. to be "y"


    %% CALCULATE HSV FOR HORIZONTAL SLICE
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
    hi2=imref2d(size(robin));
    hi2.XWorldLimits = [xwords(1) xwords(end)];
    hi2.YWorldLimits = [ywords(1) ywords(end)];


    lbound_UD=1;
    ubound_UD=max(size(MLAT(1,:)));

    lbound_LR=1;
    ubound_LR=max(size(MLAT(:,1)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Plot all the things!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PLOT THE DIFFERENT CONTRIBUTIONS TO FAC - CONVERTED TO K's SIGN CONVENTION
    figure(24);
    set(gcf,'PaperPosition',[0 0 14 10]);

    cat2=subplot(423);
    himage=imshow(robin,hi2);
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
    title('Flow from model','FontSize',sized)

    cat4=subplot(421);
    %pcolor(MLON,MLAT,vEW) % E-W
    vEW=vEW(lbound_LR:ubound_LR,lbound_UD:ubound_UD);
    pcolor(MLON(lbound_LR:ubound_LR,lbound_UD:ubound_UD),MLAT(lbound_LR:ubound_LR,lbound_UD:ubound_UD),vEW);

    set(gca,'FontSize',sized)
    shading flat
    clb=colorbar;
    colormap(cat4,'parula')
    ylabel('Mag. Lat.','FontSize',sized)
    xlabel('Mag. Lon.','FontSize',sized)
    ylabel(clb,'E-W Velocity (m/s)','FontSize',sized)
    caxis([-max(max(abs(vEW))) max(max(abs(vEW)))])
    title('E-W Velocity from Model')

    cat4=subplot(422);
    %pcolor(MLON,MLAT,vNS) % N-S
    vNS=vNS(lbound_LR:ubound_LR,lbound_UD:ubound_UD);
    pcolor(MLON(lbound_LR:ubound_LR,lbound_UD:ubound_UD),MLAT(lbound_LR:ubound_LR,lbound_UD:ubound_UD),vNS);
    set(gca,'FontSize',sized)
    shading flat
    clb=colorbar;
    colormap(cat4,'parula')
    ylabel('Mag. Lat.','FontSize',sized)
    xlabel('Mag. Lon.','FontSize',sized)
    ylabel(clb,'N-S Velocity (m/s)','FontSize',sized)
    caxis([-max(max(abs(vNS))) max(max(abs(vNS)))])
    title('N-S Velocity from Model')

    cat3=subplot(424);
    %pcolor(MLON,MLAT,Jpar./1e-6);
    Jpar=Jpar(lbound_LR:ubound_LR,lbound_UD:ubound_UD);
    pcolor(MLON(lbound_LR:ubound_LR,lbound_UD:ubound_UD),MLAT(lbound_LR:ubound_LR,lbound_UD:ubound_UD),Jpar./1e-6);
    
    set(gca,'FontSize',sized)
    shading flat
    % axis square;
    ylabel('Mag. Lat.','FontSize',sized)
    xlabel('Mag. Lon.','FontSize',sized)
    title('J_{||} from model','FontSize',sized)
    caxis([-max(max(abs(Jpar./1e-6))) max(max(abs(Jpar./1e-6)))]);
    clb=colorbar;
    map2=colorcet('D1A');
    colormap(cat3,map2)
    ylabel(clb,'Current (\muA/m^2)','FontSize',sized)

    cat4=subplot(425);
    %pcolor(MLON,MLAT,SIGP);
    SIGP=SIGP(lbound_LR:ubound_LR,lbound_UD:ubound_UD);
    pcolor(MLON(lbound_LR:ubound_LR,lbound_UD:ubound_UD),MLAT(lbound_LR:ubound_LR,lbound_UD:ubound_UD),SIGP);
    
    set(gca,'FontSize',sized)
    shading flat
    % axis square;
    ylabel('Mag. Lat.','FontSize',sized)
    xlabel('Mag. Lon.','FontSize',sized)
    title('\Sigma_P from model','FontSize',sized)
    % colorcet('D1A')
    colormap(cat4,'parula')
    clb=colorbar;
    ylabel(clb,'Conductance (mhos)','FontSize',sized)
    % cax=[-10e-6 10e-6];  % hardcoded range, should change
    % caxis(cax);

    cat5=subplot(426);
    %pcolor(MLON,MLAT,SIGH);

    SIGH=SIGH(lbound_LR:ubound_LR,lbound_UD:ubound_UD);
    pcolor(MLON(lbound_LR:ubound_LR,lbound_UD:ubound_UD),MLAT(lbound_LR:ubound_LR,lbound_UD:ubound_UD),SIGH);
    set(gca,'FontSize',sized)
    shading flat
    % axis square;
    ylabel('Mag. Lat.','FontSize',sized)
    xlabel('Mag. Lon.','FontSize',sized)
    title('\Sigma_H from model','FontSize',sized)
    % caxis(cax);
    colormap(cat5,'parula')
    clb=colorbar;
    ylabel(clb,'Conductance (mhos)','FontSize',sized)

    cat5=subplot(427);
    %pcolor(MLON,MLAT,Q);

    Q=Q(lbound_LR:ubound_LR,lbound_UD:ubound_UD);
    pcolor(MLON(lbound_LR:ubound_LR,lbound_UD:ubound_UD),MLAT(lbound_LR:ubound_LR,lbound_UD:ubound_UD),Q);
    set(gca,'FontSize',sized)
    shading flat
    % axis square;
    ylabel('Mag. Lat.','FontSize',sized)
    xlabel('Mag. Lon.','FontSize',sized)
    title('Q model input','FontSize',sized)
    % caxis(cax);
    colormap(cat5,'parula')
    clb=colorbar;
    ylabel(clb,'Energy Flux (mW/m^2)','FontSize',sized)

    cat5=subplot(428);
    %pcolor(MLON,MLAT,E0);

    E0=E0(lbound_LR:ubound_LR,lbound_UD:ubound_UD);
    pcolor(MLON(lbound_LR:ubound_LR,lbound_UD:ubound_UD),MLAT(lbound_LR:ubound_LR,lbound_UD:ubound_UD),E0);
    set(gca,'FontSize',sized)
    shading flat
    % axis square;
    ylabel('Mag. Lat.','FontSize',sized)
    xlabel('Mag. Lon.','FontSize',sized)
    title('E0 model input','FontSize',sized)
    % caxis(cax);
    colormap(cat5,'parula')
    clb=colorbar;
    ylabel(clb,'Characteristic Energy (eV)','FontSize',sized)

    if(flagplot)
        print('-dpng',[outputdir,filesep,'FlowCurrentConductances_',num2str(UTsec),'.png'],'-r300');
    end
    if(savedata)
        save([outputdir,filesep,'FlowCurrentConductances_',num2str(UTsec),'.mat'],'Jpar','vEW','vNS','SIGH','SIGP','Q','E0','MLON','MLAT','UTsec','RunName')
    end




    figure(25);
    set(gcf,'PaperPosition',[0 0 8 3]);

    mer=subplot(121);
    pcolor(MLON,MLAT,SIGP);
    shading flat
    % axis square;
    xlabel('Mag. Dist. E (km)');
    ylabel('Mag. N-S (km)');
    title('\Sigma_P from model')
    % colorcet('D1A')
    colormap(mer,'parula')
    clb=colorbar;
    ylabel(clb,'Conductance (mhos)')
    % cax=[-10e-6 10e-6];  % hardcoded range, should change
    % caxis(cax);

    mar=subplot(122);
    pcolor(MLON,MLAT,SIGH);
    shading flat
    % axis square;
    xlabel('Mag. Dist. E (km)');
    ylabel('Mag. N-S (km)');
    title('\Sigma_H from model')
    % caxis(cax);
    colormap(mer,'parula')
    clb=colorbar;
    ylabel(clb,'Conductance (mhos)')


    if(flagplotC)
        print('-dpng',[outputdirC,filesep,'Conductances_',num2str(UTsec),'.png'],'-r300');
    end





    % COMPUTE VARIOUS TERMS OF THE CURRENT CONTINUITY EQUATION - MATT's SIGN CONVENTION!!!
    x2=xg.x2(3:end-2);    %strip off ghost cells
    x3=xg.x3(3:end-2);
    x1=xg.x1(3:end-2);
    lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);

    [v,E]=gemscr.postprocess.Efield(xg, dat.v2, dat.v3);                    %Electric field and drift vectors;
    E2=E(:,:,:,2);
    E3=E(:,:,:,3);

    divE=divergence(x2,x3,squeeze(E2(end,:,:))',squeeze(E3(end,:,:))');     %all of the transposing in what follow is just a way for me to deal with matlab's x vs. y expectations in terms of ordering of indices
    divE=divE';
    JdivE=-1*(SIGP.*divE);

    [gradSIGP2,gradSIGP3]=gradient(SIGP',x2,x3);
    gradSIGP2=gradSIGP2';                            %go back to GEMINI index ordering to avoid confusion
    gradSIGP3=gradSIGP3';
    JgradSIGP=-1*(gradSIGP2.*squeeze(E2(end,:,:))+gradSIGP3.*squeeze(E3(end,:,:)));

    bhat=cat(4,-1*ones(lx1,lx2,lx3),zeros(lx1,lx2,lx3),zeros(lx1,lx2,lx3));    %unit vector along field in curvilinear basis - opposite of the z-direction in my simulation since northern hemisphere
    [gradSIGH2,gradSIGH3]=gradient(abs(SIGH'),x2,x3);                          %abs due to my wacky sign convention of positive x1 is up...
    gradSIGH=cat(3,gradSIGH2',gradSIGH3');           %3rd index is components, go back to index ordering used in GEMINI with tranposes...
    bhatxE=cross(bhat,E,4);                          %cross product to be taken along 4th dim of arrays
    bhatxE=squeeze(bhatxE(end,:,:,2:3));             %should be constant along x1, just use the final cell - also only need to x2 and x3 components since field-line integrated...
    JgradSIGH=-1*dot(gradSIGH,bhatxE,3);

    Jfac=cat(4, dat.J1,zeros(lx1,lx2,lx3),zeros(lx1,lx2,lx3));    %field aligned current vector - from current_decompose.m


    % PLOT THE DIFFERENT CONTRIBUTIONS TO FAC - CONVERTED TO K's SIGN CONVENTION
    figure(26);
    set(gcf,'PaperPosition',[0 0 12 3]);

    jmax=max([max(Jfac(:)),max(JdivE(:)),max(JgradSIGP(:)),max(JgradSIGH(:))]);
    jmin=min([min(Jfac(:)),min(JdivE(:)),min(JgradSIGP(:)),min(JgradSIGH(:))]);
    jmin=abs(jmin);
    total_max=max([jmin,jmax]);
    cax_arr=[-1e-5,1e-5];
    cax_small=[-2e-6,2e-6];

    hat5=subplot(235);
    pcolor(x2/1e3,x3/1e3,-squeeze(Jfac(end,lbound_LR:ubound_LR,lbound_UD:ubound_UD,1))');
    shading flat
    % axis square;
    xlabel('Mag. Dist. E (km)');
    ylabel('Mag. N-S (km)');
    title('J_{||} from model')
    map10=colorcet('D1A');
    colormap(hat5,map10)
    clb=colorbar;
    ylabel(clb,'Current (A/m^{2})')
    %cax=[-40e-6 40e-6];  % hardcoded range, should change
    caxis(cax_arr);

    hat5=subplot(231);
    pcolor(x2/1e3,x3/1e3,-JdivE(lbound_LR:ubound_LR,lbound_UD:ubound_UD)');
    shading flat
    % axis square;
    xlabel('Mag. Dist. E (km)');
    ylabel('Mag. N-S (km)');
    title('\Sigma_P \nabla \cdot E')
    caxis(cax_arr);
    map10=colorcet('D1A');
    colormap(hat5,map10)
    clb=colorbar;
    ylabel(clb,'Current (A/m^{2})')

    hat5=subplot(232);
    pcolor(x2/1e3,x3/1e3,-JgradSIGP(lbound_LR:ubound_LR,lbound_UD:ubound_UD)');
    shading flat
    % axis square;
    xlabel('Mag. Dist. E (km)');
    ylabel('Mag. N-S (km)');
    title('\nabla \Sigma_P \cdot E')
    c=max(max(-JdivE));
    cax2=[-c c];
    caxis(cax2);
    map10=colorcet('D1A');
    colormap(hat5,map10)
    clb=colorbar;
    ylabel(clb,'Current (A/m^{2})')
    caxis(cax_arr);

    hat5=subplot(233);
    pcolor(x2/1e3,x3/1e3,-JgradSIGH(lbound_LR:ubound_LR,lbound_UD:ubound_UD)');
    shading flat
    % axis square;
    xlabel('Mag. Dist. E (km)');
    ylabel('Mag. N-S (km)');
    title('\nabla \Sigma_H \cdot (b \times E)')
    %caxis(cax);
    caxis(cax_small);
    map10=colorcet('D1A');
    colormap(hat5,map10)
    clb=colorbar;
    ylabel(clb,'Current (A/m^{2})')

    hat5=subplot(234);
    pcolor(x2/1e3,x3/1e3,-(JdivE(lbound_LR:ubound_LR,lbound_UD:ubound_UD)+JgradSIGP(lbound_LR:ubound_LR,lbound_UD:ubound_UD)+JgradSIGH(lbound_LR:ubound_LR,lbound_UD:ubound_UD))');
    shading flat
    % axis square;
    xlabel('Mag. Dist. E (km)');
    ylabel('Mag. N-S (km)');
    title('J_{||} from sum of sources')
    caxis(cax_arr);
    map10=colorcet('D1A');
    colormap(hat5,map10)
    clb=colorbar;
    ylabel(clb,'Current (A/m^{2})')
    if(flagplotT)
        print('-dpng',[outputdirT,filesep,'JfacDecomp_',num2str(UTsec),'.png'],'-r300');
    end





end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% I'm including the hsv functions below for convience...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function h = show(im, varargin)
        arg=varargin;
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HSV panel colormap alternatives (including perceptually uniform and
% colorblind-friendly) provided by:
% Peter Kovesi. "Good Colour Maps: How to Design Them."
% arXiv:1509.03700 [cs.GR] 2015.
% "Many widely used colour maps have perceptual flat spots that can hide
% features as large as 10% of your total data range. They may also have
% points of locally high colour contrast leading to the perception of false
% features in your data when there are none. MATLAB's 'hot', 'jet', and
% 'hsv' colour maps suffer from these problems. "
% https://www.peterkovesi.com/matlabfns/index.html#colour
% Needed functions available online: applycolourmap.m, colorcet.m,
% normalise.m, show.m, showangularim.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

