%Create the grid variables and put them in the current workplace
%this is a nice little shorthand - load the grid, then call this and 
%you have a bunch of useful grid variables
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