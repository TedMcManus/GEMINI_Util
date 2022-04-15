function [] = plot_multirun(inp,time,xg,legendtitle,legendnames,mode)
%% This function is designed to plot the data gathered by load_multirun
% First, run data=load_multirun(folderlist), where folderlist is a string array containing
% names of folders which contain GEMINI output. Then, run this function w/ inp as the output
% data struct from load_multirun, time as the timestep to plot (usually a number from 1 to 31),
% the grid data, a legend title, and names for each run in the list. If mode is 'auto', we
% autoscale the colorbar ranges, otherwise, everything is on the same scale

%% Stuff you might want to change
latlow=53; %low latitude plot limit
lathigh=57; %high latitude plot limit

Fsized=16; %label font size
Fsized_leg=14; %legend font size
Tsized=24; %title font size

%% Stuff you won't want to change

names=fieldnames(inp);

%for some reason we need to index backwards. Idk why but it works
order=linspace(length(legendnames),1,length(legendnames));

hold on;
clf %just in case

%We're plotting with respect to mag lat
MLAT=90-squeeze(xg.theta(1,:,:))*180/pi;
mlat=MLAT(1,:);

%clip the lat ranges (useful for visual clarity)
[~,ind]=min(abs(mlat-latlow));
lbound_left=ind;

[~,ind]=min(abs(mlat-lathigh));
lbound_right=ind;

prange=[lbound_left:lbound_right];
mlatp=mlat(prange);

if ~exist("legendnames",'var')
    legendnames=names;
end

if ~exist("mode",'var')
    mode='auto';
end

%these are dummy variables that get overwritten on the first loop
dmax_max=0;
dmin_min=1000000000;

if ~strcmp(mode,'auto')
    for i = order
        %Load all the relevant data
        current_struct=inp.(char(names(i)));
        dens_95=current_struct.('dens_95');
        dens_120=current_struct.('dens_120');
        dens_300=current_struct.('dens_300');
        dens_800=current_struct.('dens_800');

        %and clip it to the latitude range
        dens_95p=dens_95(prange,time);
        dens_120p=dens_120(prange,time);
        dens_300p=dens_300(prange,time);
        dens_800p=dens_800(prange,time);

        %find the overall max for this sim
        dmax=max([max(dens_95p(:)),max(dens_120p(:)),max(dens_300p(:)),max(dens_800p(:))]);
        dmin=min([min(dens_95p(:)),min(dens_120p(:)),min(dens_300p(:)),min(dens_800p(:))]);

        %and use that classic max algorithm to find the max over all sims
        %we use this to give an overall colorbar for the plots if we aren't autoscaling
        if dmax>dmax_max
            dmax_max=dmax;
        end
        if dmin_min>dmin
            dmin_min=dmin;
        end
    end
    dmax=dmax_max;
    dmin=dmin_min;
end

figure(1)
for i = order
    %grab all the density data
    current_struct=inp.(char(names(i)));
    dens_95=current_struct.('dens_95');
    dens_120=current_struct.('dens_120');
    dens_300=current_struct.('dens_300');
    dens_800=current_struct.('dens_800');
    
    %and clip it
    dens_95p=dens_95(prange,time);
    dens_120p=dens_120(prange,time);
    dens_300p=dens_300(prange,time);
    dens_800p=dens_800(prange,time);

    %make a nice figure, tighten the axes, and create subplots
    set(gcf, 'Position',  [0, 0, 2000, 1000]);
    a=subplot(2,2,1);
    axis tight;

    %Plot the densities in altitude order and name them 
    hold(a,'on');
    plot(a,mlatp,dens_95p,'DisplayName',char(legendnames(i)));
    if ~strcmp(mode,'auto')
        ylim([dmin,dmax]);
    end

    b=subplot(2,2,2);
    axis tight;
    hold(b,'on');
    plot(b,mlatp,dens_120p,'DisplayName',char(legendnames(i)));
    if ~strcmp(mode,'auto')
        ylim([dmin,dmax]);
    end

    c=subplot(2,2,3);
    axis tight;
    hold(c,'on');
    plot(c,mlatp,dens_300p,'DisplayName',char(legendnames(i)));
    if ~strcmp(mode,'auto')
        ylim([dmin,dmax]);
    end

    d=subplot(2,2,4);
    axis tight;
    hold(d,'on');
    plot(d,mlatp,dens_800p,'DisplayName',char(legendnames(i)));
    if ~strcmp(mode,'auto')
        ylim([dmin,dmax]);
    end
end
lgd=legend;

%Make labels and titles using latex interpreter
title(a,'95 km','Interpreter','latex','FontSize',Tsized);
xlabel(a,'Mag Lat (deg)','Interpreter','latex','FontSize',Fsized);
ylabel(a,'Log\(_{10}\) of Density \((m^{-3})\)  ','Interpreter','latex','FontSize',Fsized);

title(b,'120 km','Interpreter','latex','FontSize',Tsized);
xlabel(b,'Mag Lat (deg)','Interpreter','latex','FontSize',Fsized);
ylabel(b,'Log\(_{10}\) of Density \((m^{-3})\)  ','Interpreter','latex','FontSize',Fsized);

title(c,'300 km','Interpreter','latex','FontSize',Tsized);
xlabel(c,'Mag Lat (deg)','Interpreter','latex','FontSize',Fsized);
ylabel(c,'Log\(_{10}\) of Density \((m^{-3})\)  ','Interpreter','latex','FontSize',Fsized);

title(d,'800 km','Interpreter','latex','FontSize',Tsized);
xlabel(d,'Mag Lat (deg)','Interpreter','latex','FontSize',Fsized);
ylabel(d,'Log\(_{10}\) of Density \((m^{-3})\)  ','Interpreter','latex','FontSize',Fsized);

sgtitle('Density Comparison for Winter Non-Counterstreaming Events','Interpreter','latex','FontSize',Tsized)


figure(2)
for i = order

    %load in relevant vars
    current_struct=inp.(char(names(i)));
    Jpar=current_struct.('Jpar');
    SigmaH=current_struct.('SigmaH');
    SigmaP=current_struct.('SigmaP');
    V_2=current_struct.('V_2');

    %and clip them to size
    Jparp=Jpar(prange,time);
    SigmaHp=SigmaH(prange,time);
    SigmaPp=SigmaP(prange,time);
    V2p=V_2(prange,time);

    
    set(gcf, 'Position',  [0, 0, 2000, 1000]);
    lgd2=legend;

    %Now plot flow and conductances, labeling them properly
    a=subplot(2,2,1);
    axis tight;
    hold(a,'on');
    plot(a,mlatp,Jparp*10^6,'DisplayName',char(legendnames(i)));

    b=subplot(2,2,2);
    axis tight;
    hold(b,'on');
    plot(b,mlatp,SigmaHp,'DisplayName',char(legendnames(i)));

    c=subplot(2,2,3);
    axis tight;
    hold(c,'on');
    plot(c,mlatp,SigmaPp,'DisplayName',char(legendnames(i)));

    d=subplot(2,2,4);
    axis tight;
    hold(d,'on');
    plot(d,mlatp,-V2p,'DisplayName',char(legendnames(i)));

end

axis tight
title(a,'Parallel Current','Interpreter','latex','FontSize',Tsized);
xlabel(a,'Mag Lat (deg)','Interpreter','latex','FontSize',Fsized);
ylabel(a,'Current \( \left ( \mu A/m^{-3} \right ) \)','Interpreter','latex','FontSize',Fsized);

title(b,'Hall Conductance','Interpreter','latex','FontSize',Tsized);
xlabel(b,'Mag Lat (deg)','Interpreter','latex','FontSize',Fsized);
ylabel(b,'Conductivity (mhos)','Interpreter','latex','FontSize',Fsized);

title(c,'Pedersen Conductance','Interpreter','latex','FontSize',Tsized);
xlabel(c,'Mag Lat (deg)','Interpreter','latex','FontSize',Fsized);
ylabel(c,'Conductivity (mhos)','Interpreter','latex','FontSize',Fsized);

title(d,'Cross-Track Westwards Flow','Interpreter','latex','FontSize',Tsized);
xlabel(d,'Mag Lat (deg)','Interpreter','latex','FontSize',Fsized);
ylabel(d,'Velocity (m/s)','Interpreter','latex','FontSize',Fsized);


%Legend titles and positioning go here b/c it works...
%You might want to alter this, I clearly fine-tuned it by hand with the property editor
title(lgd2,legendtitle,'Interpreter','latex','FontSize',Fsized_leg)

title(lgd,legendtitle,'Interpreter','latex','FontSize',Fsized_leg)

set(lgd2,'Position',[0.010262323417659,0.431308819957644,0.112136176519731,0.079419193556814]);
set(lgd,'Position',[0.010262323417659,0.431308819957644,0.112136176519731,0.079419193556814]);

%overall title for plot 2
sgtitle('Current/Conductivity Comparison for Winter Non-Counterstreaming Events','Interpreter','latex','FontSize',Tsized)

end


