MLAT=90-squeeze(xg.theta(1,:,:))*180/pi;
MLON=squeeze(xg.phi(1,:,:))*180/pi;
mlat=MLAT(1,:);
time=[1:31];
time=time*10;
[LAT,T]=meshgrid(time,mlat);

figure(1);
set(gcf, 'Position',  [0, 0, 2000, 1000])

subplot(4,2,1);
pcolor(LAT,T,dens_150);shading flat;
ylabel('Latitude');
xlabel('Elapsed Time (sec)');
title('Density at 150 km');
colormap('Parula');
colorbar;

subplot(4,2,2);
pcolor(LAT,T,dens_200);shading flat;
ylabel('Latitude');
xlabel('Elapsed Time (sec)');
title('Density at 200 km');
colorbar;


subplot(4,2,3);
pcolor(LAT,T,dens_300);shading flat;
ylabel('Latitude');
xlabel('Elapsed Time (sec)');
title('Density at 300 km');
colorbar;


jax=subplot(4,2,4);
pcolor(LAT,T,J_par);shading flat;
ylabel('Latitude');
xlabel('Elapsed Time (sec)');
title('Current (\muA)');
J=abs(Jpar);
jmax=max(J(:));
colormap(jax,colorcet('D1A'));
caxis([-jmax,jmax]);
colorbar;


subplot(4,2,7);
pcolor(LAT,T,V_2);shading flat;
ylabel('Latitude');
xlabel('Elapsed Time (sec)');
title('Flow in x2 (km/s)');
colorbar;

subplot(4,2,8);
pcolor(LAT,T,V_3);shading flat;
ylabel('Latitude');
xlabel('Elapsed Time (sec)');
title('Flow in x3 (km/s)');
colorbar;


subplot(4,2,5);
pcolor(LAT,T,SigmaH);shading flat;
ylabel('Latitude');
xlabel('Elapsed Time (sec)');
title('Hall Conductance');
colorbar;

subplot(4,2,6);
pcolor(LAT,T,SigmaP);shading flat;
ylabel('Latitude');
xlabel('Elapsed Time (sec)');
title('Pedersen Conductance');
colorbar;

%saveas(gcf,[direc,filesep,'keogram.png']);