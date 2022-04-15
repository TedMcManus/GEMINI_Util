%% Creates a 4-panel movie with density cuts at various altitudes and a current panel
%Assumes the current matlab workspace has a struct called "out" which is the  
%result of running the loadvar script. Out must contain both density and current data
%This version of the script lets the user clip the latitude range of the output plots

%Name the video!
vidfile = VideoWriter('./Dens_Movies/Alt.mp4','MPEG-4');

%Lat ranges
latlow=63;
lathigh=67;

MLAT=90-squeeze(xg.theta(1,:,:))*180/pi;
MLON=squeeze(xg.phi(1,:,:))*180/pi;
llon=xg.lx(2);
llat=xg.lx(3);
mlon=MLON(:,1);
mlat=MLAT(1,:);
mlonmean=mean(mlon);
mlatmean=mean(mlat);
mlatspan=max(MLAT(:))-min(MLAT(:));
mlonspan=max(MLON(:))-min(MLON(:));

Tsized=20; %Title font
Fsized=10; %Font size
T_Sized=22; %Overall title size

UTsec0=36000;
dtout=10;
vidfile.FrameRate = 5;
open(vidfile);

%Density and current data
dens = out.ne;
J1=out.J1;

[~,ind]=min(abs(xg.alt(:,1,1)-95000));
dens_95=log10(squeeze(dens(ind,:,:,:)));

 [~,ind]=min(abs(xg.alt(:,1,1)-120000));
dens_120=log10(squeeze(dens(ind,:,:,:)));

 [~,ind]=min(abs(xg.alt(:,1,1)-300000));
dens_300=log10(squeeze(dens(ind,:,:,:)));

 [~,ind]=min(abs(xg.alt(:,1,1)-800000));
dens_800=log10(squeeze(dens(ind,:,:,:)));

J1=squeeze(J1(end,:,:,:));

%calculate the cumulative colorbar for all the data by taking the max and min of all the data
%this ranges over all the timesteps so it's nice and consistent
maxarr=[max(dens_95(:)),max(dens_120(:)),max(dens_300(:)),max(dens_800(:))];
minarr=[min(dens_95(:)),min(dens_120(:)),min(dens_300(:)),min(dens_800(:))];
max_=max(maxarr(:));
min_=min(minarr(:));
cbar=[log10(min_),log10(max_)];

%current colorbar is separate
JMax=max(abs(J1(:)));
Jbar=[-JMax,JMax]*10^6;

[~,ind]=min(abs(mlat-latlow));
mltlow=ind;

%Clip all the data to the specified ranges
[~,ind]=min(abs(mlat-lathigh));
mlthigh=ind;


MLAT=MLAT(:,mltlow:mlthigh);
MLON=MLON(:,mltlow:mlthigh);
dens_95=dens_95(:,mltlow:mlthigh,:);
dens_120=dens_120(:,mltlow:mlthigh,:);
dens_300=dens_300(:,mltlow:mlthigh,:);
dens_800=dens_800(:,mltlow:mlthigh,:);
J1=J1(:,mltlow:mlthigh,:);




f = figure('units','normalized','outerposition',[0 0 1 1]);

%This part gets a little gnarly. We're making a 5-panel plot where every panel is at a different altitude.
%We're also doing manual positioning to make sure everything fits on the plot. 
for it = 1:length(out.times)

    ax1=subplot(5,1,1); %create a grid of 5 subplots
    curr=squeeze(J1(:,:,it))*10^6; %get the current and scale it to be in microamps
    pcolor(ax1,MLON,MLAT,curr);shading flat; %make a plot
    colormap(ax1,colorcet('D1A'));
    %xlabel(ax1,'Mag lon','Interpreter','latex','FontSize',Fsized);
    %ylabel(ax1,'Mag lat','Interpreter','latex','FontSize',Fsized);
    colorbar(ax1); %make a colorbar 
    ax1.XTickLabel=[]; %remove the tick labels bc they're ugly
    ax1.YTickLabel=[];
    caxis(ax1,Jbar);
    %axis square;
    set(gca, 'PlotBoxAspectRatio', [3.3,1,1]); %manually set aspect ratio for plots
    axis tight; %reduce space between subplots
    yyaxis right; %target the right side of the plots for the label
    %I'm using LaTeX to format all the plots so they look the same. 
    %The weird syntax in the labels is all LaTeX garbage
    ylabel(ax1,'Topside','FontSize',Fsized,'Interpreter','latex',Color='black');
    a = gca;
    a.YTickLabel={};

    %Here, we create a colorbar and position/size it manually
    %I used the property inspector in the plot window to get these numbers
    J_Color=colorbar();
    set(J_Color, 'Position', [0.665126012092531,0.108212809917355,0.01425378548896,0.817431646354749],'units','normalized');
    ylabel(J_Color,'Current\((\mu A/m^2)\)','Interpreter','latex','FontSize',Fsized);

    %The other plots are pretty similar 

    ax1=subplot(5,1,5);
    pcolor(ax1,MLON,MLAT,log10(dens_95(:,:,it)));shading flat;
    xlabel(ax1,'Mag lon','Interpreter','latex','FontSize',Fsized);
    ylabel(ax1,'Mag lat','Interpreter','latex','FontSize',Fsized);
    %t=title(ax1,'95 km','FontSize',Fsized,'Interpreter','latex');
    %colorbar(ax1);
    caxis(cbar);
    %axis square;
    set(gca, 'PlotBoxAspectRatio', [3.3,1,1]);
    axis tight;
    yyaxis right;
    ylabel('95 km','color','black','FontSize',Fsized,'Interpreter','latex');
    a = gca;
    a.YTickLabel={};
    colormap(ax1,"parula");

    ax1=subplot(5,1,4);
    pcolor(ax1,MLON,MLAT,log10(dens_120(:,:,it)));shading flat;
    %xlabel(ax1,'Mag lon','Interpreter','latex','FontSize',Fsized);
    %ylabel(ax1,'Mag lat','Interpreter','latex','FontSize',Fsized);
    %title(ax1,'120 km','FontSize',Fsized,'Interpreter','latex');
    %colorbar(ax1);
    ax1.XTickLabel=[];
    ax1.YTickLabel=[];
    caxis(cbar);
    %axis square;
    set(gca, 'PlotBoxAspectRatio', [3.3,1,1]);
    axis tight;
    yyaxis right;
    ylabel('120 km','color','black','FontSize',Fsized,'Interpreter','latex');
    a = gca;
    a.YTickLabel={};

    ax1=subplot(5,1,3);
    pcolor(ax1,MLON,MLAT,log10(dens_300(:,:,it)));shading flat;
    %xlabel(ax1,'Mag lon','Interpreter','latex','FontSize',Fsized);
    %ylabel(ax1,'Mag lat','Interpreter','latex','FontSize',Fsized);
    %title(ax1,'300 km','FontSize',Fsized,'Interpreter','latex');
    %colorbar(ax1);
    ax1.XTickLabel=[];
    ax1.YTickLabel=[];
    caxis(cbar);
    %axis square;
    set(gca, 'PlotBoxAspectRatio', [3.3,1,1]);
    axis tight;
    yyaxis right;
    ylabel('300 km','color','black','FontSize',Fsized,'Interpreter','latex');
    a = gca;
    a.YTickLabel={};

    ax1=subplot(5,1,2);
    pcolor(ax1,MLON,MLAT,log10(dens_800(:,:,it)));shading flat;
    %xlabel(ax1,'Mag lon','Interpreter','latex','FontSize',Fsized);
    %ylabel(ax1,'Mag lat','Interpreter','latex','FontSize',Fsized);
    %title(ax1,'800 km','FontSize',Fsized,'Interpreter','latex');
    %colorbar(ax1);
    ax1.XTickLabel=[];
    ax1.YTickLabel=[];
    caxis(cbar);
    %axis square;
    set(gca, 'PlotBoxAspectRatio', [3.3,1,1]);
    axis tight;
    yyaxis right;
    ylabel('800 km','color','black','FontSize',Fsized,'Interpreter','latex');
    a = gca;
    a.YTickLabel={};

    %Here, we're creating and positioning the cumulate colorbar 
    h=colorbar;
    set(h, 'Position', [0.360929771293373,0.110278925619835,0.01425378548896,0.814407909074143],'units','normalized');
    ylabel(h,'Density\((m^{-3})\)','Interpreter','latex','FontSize',Fsized);
    UT=UTsec0+(out.times(it)-1)*dtout;
    
    %Cumulative plot title 
    sgtitle([char(string(UT)),' UTSEC'],'FontSize',Tsized,'Interpreter','latex');
    F =getframe(figure(1));
    writeVideo(vidfile,F);
    clf;
end

close(vidfile);
fclose('all');
