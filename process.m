function[]=process(direc,start,dur,stop)
%Ideally just run process('./data_directory') and everything works ;)

arguments
    direc (1,1) string %This is the directory with the data
    start (1,1) double = 0 %start time
    dur (1,1) double = 10 %Cadence
    stop  (1,1) double = 300 %timestep at which to halt
end
    

if ~exist('direc','var')
	fprintf('direc start dur stop')
end

if (~exist('ymd0','var'))
    confdat=gemini3d.read.config(direc);
    fprintf('Input config file loaded.\n');
end

UTsec0=confdat.UTsec0;

disp('Beginning HSV plotting routine');
HSV(direc,175,start,dur,stop,'both');

disp('Beginning HSV countour plotting routines');
contour_hsv(direc,300,start,dur,stop,'both');

disp('Beginning keogram plotting routine');
keogram(direc,'keogram_data.mat',start,dur,stop);

disp('Beginning FCC plotting routine');
FCC(direc,UTsec0,dur,UTsec0+stop);

direc=char(direc);

if strcmp(direc(1:2),'./')
	direc=direc(3:end);
end

outdir=['~/../public_html/Gemini3D/',direc];
mkdir(outdir);
system(['scp -r ./',direc,'/* ',outdir,'/']);

fprintf('%s\n','BAZN0NGA');
