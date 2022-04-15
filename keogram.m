function [] = keogram(direc,outfile,start,dur,stop)
tic
%% OVERVIEW
%This ungodly monster plotter needs to plot
%Density (150 km, 200 km, 300 km, 800km)
%Flow in x2 and x3 directions
%Current
%SigH and SigP
%%
direc = char(direc);
outfile = char(outfile);
%thedate='20150620';

confdat=gemini3d.read.config(direc);
ymd0=confdat.ymd;
thedate=char(strrep(strjoin(string(ymd0)),' ',''));
if (length(char(string(ymd0(2)))))==1
    thedate=[thedate(1:4),'0',thedate(5:end)];
end
UTsec0=confdat.UTsec0;
tdur=confdat.tdur;
dtout=confdat.dtout;
flagoutput=confdat.flagoutput;

xg = gemini3d.read.grid(direc);
cfg = gemini3d.read.config(direc);

dens_95=zeros(xg.lx(3),31);
dens_120=dens_95;
dens_300=dens_95;
dens_800=dens_95;
V_2=dens_95;
V_3=dens_95;
J_par=dens_95;
SigmaH=dens_95;
SigmaP=dens_95;

it=0;
start=start+cfg.UTsec0;
stop=stop+cfg.UTsec0;
broken=false;


MLAT=90-squeeze(xg.theta(1,:,:))*180/pi;
MLON=squeeze(xg.phi(1,:,:))*180/pi;
mlat=MLAT(1,:);
time=[1:31];
time=time*10;

[~,lbound]=min(abs(mlat-52));
[~,ubound]=min(abs(mlat-58));
[TIME,LAT]=meshgrid(time,mlat(lbound:ubound));

%% LOAD THE SIMULATION DATA
for UTsec= start:dur:stop
    if exist([direc,filesep,outfile],'file')
        load([direc,filesep,outfile]);
        disp('Located data file, skipping data creation');
        broken=true;
        break
    end
    % for UTsec=35920:10:36150
    % for UTsec=36000  %use if only a specific timestep is needed
    it=it+1;
    time = datetime(ymd0) + seconds(UTsec);

    dat = gemini3d.read.frame([direc,filesep,thedate,'_',num2str(UTsec),'.000000.h5']);

    [~,~,~,SIGP,SIGH,~,~] = gemscr.postprocess.conductivity_reconstruct(time,dat,cfg, xg); %ram hog and sloooooow
    SIGH=-SIGH;                       %set to K sign
    Jpar=-squeeze(dat.J1(end,:,:));   %set to K sign
    v2=squeeze(dat.v2(end,:,:));
    v3=squeeze(dat.v3(end,:,:));

    %% Load Density
    ne=dat.ne;
    parmne=log10(ne);
    keoslice_95=squeeze(parmne(8,75,:)); %255 lon and 150 km
    dens_95(:,it)=keoslice_95;
    keoslice_120=squeeze(parmne(21,75,:)); %255 lon and 200 km
    dens_120(:,it)=keoslice_120;
    keoslice_300=squeeze(parmne(96,75,:)); %255 lon and 300 km
    dens_300(:,it)=keoslice_300;
    keoslice_800=squeeze(parmne(157,75,:)); %255 lon and 800 km
    dens_800(:,it)=keoslice_800;

    %% Load Current
    Jpar=Jpar(75,:);
    J_par(:,it)=Jpar;

    %% Load Flows
    V_2(:,it)=v2(75,:);
    V_3(:,it)=v3(75,:);

    %% Load Conductances
    SigmaH(:,it)=SIGH(75,:);
    SigmaP(:,it)=SIGP(75,:);

end
%     MLAT=90-squeeze(xg.theta(1,:,:))*180/pi;
%     MLON=squeeze(xg.phi(1,:,:))*180/pi;
%     mlat=MLAT(1,:);
%     time=[1:31];
%     time=time*10;
%     [LAT,T]=meshgrid(time,mlat);
%
%     figure(1);
%     set(gcf, 'Position',  [0, 0, 2000, 1000])
%
%     subplot(4,2,1);
%     pcolor(LAT,T,dens_150);shading flat;
%     ylabel('Mag Lat');
%     xlabel('Elapsed Time (sec)');
%     title('Density at 150 km');
%     colormap('Parula');
%     colorbar;
%
%     subplot(4,2,2);
%     pcolor(LAT,T,dens_200);shading flat;
%     ylabel('Mag Lat');
%     xlabel('Elapsed Time (sec)');
%     title('Density at 200 km');
%     colorbar;
%
%
%     subplot(4,2,3);
%     pcolor(LAT,T,dens_300);shading flat;
%     ylabel('Mag Lat');
%     xlabel('Elapsed Time (sec)');
%     title('Density at 300 km');
%     colorbar;
%
%
%     jax=subplot(4,2,4);
%     pcolor(LAT,T,J_par);shading flat;
%     ylabel('Mag Lat');
%     xlabel('Elapsed Time (sec)');
%     title('Current (\muA)');
%     J=abs(Jpar);
%     jmax=max(J(:));
%     colormap(jax,colorcet('D1A'));
%     caxis([-jmax,jmax]);
%     colorbar;
%
%
%     subplot(4,2,7);
%     pcolor(LAT,T,V_2);shading flat;
%     ylabel('Mag Lat');
%     xlabel('Elapsed Time (sec)');
%     title('Flow in x2 (km/s)');
%     colorbar;
%
%     subplot(4,2,8);
%     pcolor(LAT,T,V_3);shading flat;
%     ylabel('Mag Lat');
%     xlabel('Elapsed Time (sec)');
%     title('Flow in x3 (km/s)');
%     colorbar;
%
%
%     subplot(4,2,5);
%     pcolor(LAT,T,SigmaH);shading flat;
%     ylabel('Mag Lat');
%     xlabel('Elapsed Time (sec)');
%     title('Hall Conductance');
%     colorbar;
%
%     subplot(4,2,6);
%     pcolor(LAT,T,SigmaP);shading flat;
%     ylabel('Mag Lat');
%     xlabel('Elapsed Time (sec)');
%     title('Pedersen Conductance');
%     colorbar;
%
if(~broken)


    figure(1);
    set(gcf, 'Position',  [0, 0, 2000, 1000])

    dmax=max([max(dens_95(:)),max(dens_120(:)),max(dens_300(:)),max(dens_800(:))]);
    dmin=min([min(dens_95(:)),min(dens_120(:)),min(dens_300(:)),min(dens_800(:))]);

    jmax=max([max(Jpar(:)),max(SIGH(:)),max(SIGP(:))]);

    fsized=14;
    subplot(4,2,7);
    pcolor(TIME,LAT,dens_95(lbound:ubound,:));shading flat;
    ylabel('Mag Lat');
    xlabel('Elapsed Time (sec)');
    title('\textbf{Density at 95 km (\(m^{-3}\))}','Interpreter','latex',FontSize=fsized)
    colormap('Parula');
    colorbar;
    caxis([dmin,dmax]);

    subplot(4,2,5);
    pcolor(TIME,LAT,dens_120(lbound:ubound,:));shading flat;
    ylabel('Mag Lat');
    xlabel('Elapsed Time (sec)');
    title('\textbf{Density at 120 km (\(m^{-3}\))}','Interpreter','latex',FontSize=fsized)
    colorbar;
    caxis([dmin,dmax]);


    subplot(4,2,3);
    pcolor(TIME,LAT,dens_300(lbound:ubound,:));shading flat;
    ylabel('Mag Lat');
    xlabel('Elapsed Time (sec)');
    title('\textbf{Density at 300 km (\(m^{-3}\))}','Interpreter','latex',FontSize=fsized)
    colorbar;
    caxis([dmin,dmax]);

    subplot(4,2,1);
    pcolor(TIME,LAT,dens_800(lbound:ubound,:));shading flat;
    ylabel('Mag Lat');
    xlabel('Elapsed Time (sec)');
    title('\textbf{Density at 800 km (\(m^{-3}\))}','Interpreter','latex',FontSize=fsized)
    colorbar;
    caxis([dmin,dmax]);


    ax8=subplot(4,2,8);
    pcolor(TIME,LAT,V_2(lbound:ubound,:));shading flat;
    ylabel('Mag Lat');
    xlabel('Elapsed Time (sec)');
    title('\textbf{Flow in East (\(km/s\))}','Interpreter','latex',FontSize=fsized)
    colorbar;
    colormap(ax8,colorcet('L9'));

    jax=subplot(4,2,4);
    pcolor(TIME,LAT,SigmaH(lbound:ubound,:));shading flat;
    ylabel('Mag Lat');
    xlabel('Elapsed Time (sec)');
    title('\textbf{Hall Conductance (\(mhos\))}','Interpreter','latex',FontSize=fsized)
    colormap(jax,colorcet('D1A'));
    colorbar;
    caxis([-jmax,jmax]);

    jax=subplot(4,2,6);
    pcolor(TIME,LAT,SigmaP(lbound:ubound,:));shading flat;
    ylabel('Mag Lat');
    xlabel('Elapsed Time (sec)');
    %title('Pedersen Conductance (^{\muA}/_{m^2})');
    title('\textbf{Pedersen Conductance (\(mhos\))}','Interpreter','latex',FontSize=fsized)
    colormap(jax,colorcet('D1A'));
    colorbar;
    caxis([-jmax,jmax]);

    jax=subplot(4,2,2);
    pcolor(TIME,LAT,J_par(lbound:ubound,:).*10^6);shading flat;
    ylabel('Mag Lat');
    xlabel('Elapsed Time (sec)');
    title('\textbf{Parallel Current} (\(\mu A/m^{2}\))','Interpreter','latex',FontSize=fsized)
    colormap(jax,colorcet('D1A'));
    caxis([-jmax,jmax]);
    colorbar;



    save([direc,filesep,outfile])
    saveas(gcf,[direc,filesep,'keogram.png']);
end
toc

