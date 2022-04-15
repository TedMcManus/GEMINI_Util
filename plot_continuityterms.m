%% Make plots of current from various terms of the current continuity equation
% shows overall parallel current and the Pedersen and Hall sources that make it up

%% READ THIS INFO BEFORE USING:

% First, you're gonna want to run dat=gemini3d.read.frame(frame), where frame is the filename of the frame you
% want to investigate

% Next, load THERUN/FlowCurrentConductancePlots/FlowCurrentConductances_TIMEHERE.mat, where THERUN is the name
% of the run you're looking at, ane TIMEHERE is the timestep that you want to investigate. If this folder does not exist, 
% run process(THERUN) on the simulation's output folder to generate standard plot outputs, which invcludes the aforementioned .mat file

%% Things you might want to change 
Tsized=24; %Title font size
Fsized=16; %Label font size
Lsized=12; %Legend text size

latlow=53; %Low-latitude cut 
lathigh=57; %high-latitude cut
% The plot's X axis ranges from latlow to lathigh

overall_title='Title'; %overall plot title

%% COMPUTE VARIOUS TERMS OF THE CURRENT CONTINUITY EQUATION (code by Meghan Burleigh)
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

%% Now do some plotting. Thanks Meghan! 

gemgrid;

% Create X-axis
[~,ind]=min(abs(mlat-latlow));
low=ind;

[~,ind]=min(abs(mlat-lathigh));
high=ind;

X=mlat(low:high); %This is what we're plotting over
slice=xg.lx(2)/2;

%Current from the gradient of Hall conductivity
Sig_H=squeeze(-JgradSIGH(slice,low:high)*10^6);
%Current from the gradient of Pedersen conductivity
Sig_P=squeeze(-JgradSIGP(slice,low:high)*10^6);
%Current from Pedersen conductivity times divergence of E
Div_E=squeeze(-JdivE(slice,low:high)*10^6);

%Find an overall max and min for a colorbar
arrmax=[max(Sig_H(:)),max(Sig_P(:)),max(Div_E(:))];
maximum=max(arrmax(:));

arrmin=[min(Sig_H(:)),min(Sig_P(:)),min(Div_E(:))];
minimum=min(arrmin(:));

%% Make the plot

% Create a figure and axes
fig = figure();
ax1 = axes();

% Choose colors for the various flavors of current
axColors = [
    .3        .6         .9;        % Div E is BLUEISH
    1        0    0; % Overall current is RED
    20/255       100/255    40/255;    % Hall current  is GREEN
    .7        .1         .9];       % Pedersen is PURPLE

axis(ax1, 'tight') %eliminate the space on the size
hold on

% Plot all the variables, color them, and create legend names
plot(X, Sig_H, '-','LineWidth',2, 'Color', axColors(3,:),'DisplayName','\nabla \Sigma_H \cdot (b \times E)');
plot(X, Sig_P, '-','LineWidth',2,'Color', axColors(1,:),'DisplayName','\nabla \Sigma_P \cdot E');
plot(X, Div_E, '-','LineWidth',2, 'Color', axColors(4,:),'DisplayName','\Sigma_P \nabla \cdot E');
plot(X, Sig_H+Sig_P+Div_E, '-','LineWidth',2, 'Color', axColors(2,:),'DisplayName','J_{||} from sum of sources')

% Font sizes and titles
ax1.FontSize=Fsized;
lgd=legend;
lgd.FontSize=Lsized;
xlabel('Magnetic Latitude (degrees)');
ylabel('Current (\mu A)')
title(overall_title,'FontSize',Tsized);

