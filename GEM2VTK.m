function [] = GEM2VTK(inp,times,folder,filename,xg)
%This is a script for writing a time series of VTK files
%Can be run in two modes
%% MODE 1: instruct is a structure from loadvar
%This assumes that you have loaded all the arrays you want into a structrure using "instruct=loadvar(options)"
%ex: var2VTK(structure,'./folder','base_filename',xg)

%% MODE 2: instruct is a folder containing data
%This loads all the data from the folder into a structure and then writes the data into a VTK series
%ex: var2VTK('./inputfolder','./outputfolder','base_filename',xg)

arguments
    inp ; % The result of running loadvar, or a directory where data can be found
    times = [1:31] % array of timesteps
    folder string = [inp,filesep,'VTKdata']; %Which folder to write the files to
    filename string = "vtkseries"; %Base filename
    xg = gemini3d.read.grid(inp); % A grid structure from running xg = gemini3d.read.grid(gridlocation)
end
tic;
%files are formatted like filename_1, filename_2,...,filename_N
%The time spacing need not be uniform, but it must be increasing for Paraview to read it

folder=char(folder);
filename=char(filename);

if ischar(class(inp))
    try
        inp=loadvar(inp,times,["density","flow","current","temperature"],xg,[inp,filesep,'VTK_DATA.mat']);
    catch 'Input folder not located';
    end
end

%% These are important user paramters

geogrid=0; %should we use a geographic space? if zero, use lat/lon
smallgrid=1; %rescale into (unphysical) box bounded between 0 and 1?
%If you don't do this, you will not have a nice square grid in Paraview
%However, turning this on is unrealistic compared to the real grid spacing

%% Load in relevant data
[~,errcase]=mkdir(folder);
if strcmp(errcase,'Directory already exists.')
disp(['A VTK Timeseries named ',filename, ' already exists in this location. It will be overwritten'])
    %This will be thrown if we try to write a folder that already exists, which implies we've run this
    %script before and we're going to be overwriting data. The user should know that. They can probably
    % stop us in time :)
end


%MLAT/MLON grid
MLAT=90-squeeze(xg.theta(:,:,:))*180/pi;
MLON=squeeze(xg.phi(:,:,:))*180/pi;
x=MLON;
y=MLAT;
z=xg.alt;
z=z/10000; %To make things visually appealing, alt goes by 10km increments
%otherwise, the data is an insanely tall and thin rectangle

times=inp.times;

if smallgrid %Make the grid small
    x=rescale(x);
    y=rescale(y);
    z=rescale(z);
end

if geogrid %use glat/glon
    x=xg.glon;
    y=xg.glat;
    z=xg.alt;
    z=z/10000; %To make things visually appealing, alt goes by 10km increments
    if smallgrid
        x=rescale(x);
        y=rescale(y);
        z=recale(z);
    end
end

nr_of_elements=numel(x);

tarr=times; %holdover variable from copy-paste heheh
%need not be uniformly spaced, but must be strictly increasing (I think)

vars=fieldnames(inp);

Scalars=struct; %store data
scalarnames={}; %store names
counter=0;

%Did we add scalars or vectors? This will go to one if we end up adding any
%this little bit of annoyance keeps us from accessing an empty struct
Scalar_added=0;
Vec_added=0;

if any(contains(vars,"ne")) %We're adding electron density
    Scalar_added=1; %we added some scalars
    counter=counter+1;
    X = inp.ne;
    name='Density';
    scalarnames(counter)={name};
    Scalars=setfield(Scalars,name,X); %frick you I won't do what you tell me
end

if any(contains(vars,"Te")) %We're adding electron temperature
    Scalar_added=1; %we added some scalars
    counter=counter+1;
    stop = 0;
    X = inp.Te;
    name='Temperature';
    scalarnames(counter)={name};
    Scalars=setfield(Scalars,name,X); %hehe
end

Vectors=struct; %Stores the data
vectornames={}; %Stores the field names
counter=0;
if any(contains(vars,"V1")) %We're adding flows
    Vec_added=1; %Have we added some vectors?
    counter=counter+1;
    [X] = [inp.V2;inp.V3;inp.V1];
    name='Flow';
    vectornames(counter)={name};
    Vectors=setfield(Vectors,name,X); %bad practice? yes. Does it work? yes.
end

if any(contains(vars,"J1")) %We're adding current
    Vec_added=1; %Have we added some vectors?
    counter=counter+1;
    [X] = [inp.J2;inp.J3;inp.J1];
    name='Current';
    vectornames(counter)={name};
    Vectors=setfield(Vectors,name,X); %bad practice? yes. Does it work? yes.
end


counter=0; %again, to keep track of iterations in the for loop
for i=tarr
    % Some code is from vtkwrite.m from the file exchange, some code is back-engineered from the VTK
    % specification document, some is fresh off the noggin. This is the result of an immense amount of
    % suffering. Hopefully nobody ever needs to read this comment, but if you do, appreciate my Herculean
    % efforts to slay the exception-rasing hydra that is Paraview.

    counter=counter+1;
    current_name=[folder,filesep,filename,'_',num2str(i),'.vtk']; %format the filename
    fid=fopen(current_name,'w');
    %Binary file header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'VTK file courtesy of Ted McManus and His Electric Minions\n'); %I think you can change this? Don't really know why it's needed
    %This stuff you definitely cannot change without someone throwing a fit
    fprintf(fid, 'BINARY\n\n');
    fprintf(fid, 'DATASET STRUCTURED_GRID\n');
    fprintf(fid, ['DIMENSIONS ' num2str(size(x,1)) ' ' num2str(size(x,2)) ' ' num2str(size(x,3)) '\n']);
    fprintf(fid, ['POINTS ' num2str(nr_of_elements) ' float\n']);
    fclose(fid);

    %Put in the grid
    fid = fopen(current_name, 'a');
    fwrite(fid, [reshape(x,1,nr_of_elements);  reshape(y,1,nr_of_elements); reshape(z,1,nr_of_elements)],'float','b');
    fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements)]);
    fclose(fid);

    if Vec_added
        for iter=[1:length(vectornames)]
            fid = fopen(current_name, 'a');
            X=getfield(Vectors,char(vectornames(iter))); %If you know a better way, by all means
            xlen = length(x(:,1,1));

            %This part bears explaining. If we have a time array=[1,5,10] and a NxNxNx30 data array
            %we want to take data(:,:,:,1), data(:,:,:,5) and data(:,:,:,10)
            %However, if we have t=[1,5,10] and an NxNxNx3 array, we can
            %assume that the array has been prepared such that data(N,N,N,3) is actually t=10.
            %In case 1, we want to index with the time index, but in case 2
            %we want to index sequentially. This is important because the loadvar script prepares data
            %in the case 2 format.

            if length(tarr)~=length(X(1,1,1,:))
                Xdim = X(1:xlen,:,:,i);
                Ydim = X(xlen+1:2*xlen,:,:,i);
                Zdim = X(xlen*2+1:xlen*3,:,:,i);
            else
                Xdim = X(1:xlen,:,:,counter);
                Ydim = X(xlen+1:2*xlen,:,:,counter);
                Zdim = X(xlen*2+1:xlen*3,:,:,counter);
            end

            %Print out the vectors as binary data
            fprintf(fid, ['\n\nVECTORS ', char(vectornames(iter)), ' float\n']); %binary header
            fwrite(fid, [reshape(Xdim,1,nr_of_elements);  reshape(Ydim,1,nr_of_elements); reshape(Zdim,1,nr_of_elements)],'float','b'); %binary dat
            fclose(fid);
        end
    end
    if Scalar_added
        for iter=[1:length(scalarnames)]
            %open file to append
            fid = fopen(current_name, 'a');
            X=getfield(Scalars,string(char(scalarnames(iter)))); %At least I didn't use eval
            %See the long comment above
            if length(tarr)~=length(X(1,1,1,:))
                X=X(:,:,:,i);
            else
                X=X(:,:,:,counter);
            end
            OutArr = reshape(X,1,nr_of_elements);

            %Annoying binary header
            fprintf(fid, ['\n\nSCALARS ',char(scalarnames(iter)),' float\n']);
            fprintf(fid, 'LOOKUP_TABLE default\n');
            fwrite (fid, OutArr,'float','b'); %binary data
            fclose(fid);
        end
    end
end

fclose('all');

fprintf(['And so it ends: not with a bang, but with an "Elapsed time is ',num2str(toc),' seconds"\n']);