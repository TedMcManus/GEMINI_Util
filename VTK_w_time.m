%Another nice little psuedo-function for VTK.
%This is a script for writing a time series of VTK files
%with the same variables in each. Using the loadvar script
%makes this easier, as it collects a bunch of gemini
%timesteps into a neat struct. However, this can write any
%properly sized arrays from your workspace into VTK time series format

%This assumes that you have loaded all the arrays you want into your
%workspace. Using "out=loadvar(options)" works

geogrid=1; %should we use a geographic space?
smallgrid=0; %rescale into (unphysical) box bounded between 0 and 1?
%If you don't do this, you have to recenter your view in Paraview
%However, turning this on is unrealistic compared to the real grid spacing


filename=input('base filename: ','s');
%files are formatted like filename_1, filename_2,...,filename_N
%The time spacing need not be uniform, but it must be increasing
folder=input('output directory: ','s');
%Folder is designed to be unique.
if ~exist(folder,'dir')
    mkdir(folder)
else
    warn=input(['I am about delete everything from folder ',folder,'. Kill the process to save your files :}']); %necessary? no.
    %you could maybe find a workaround to this. I made each VTK time series
    %have its own folder so fwrite doesn't complain about the file existing
    %when it tries to open it
    delete([folder,'/*']);
end

%MLAT/MLON grid
MLAT=90-squeeze(xg.theta(:,:,:))*180/pi;
MLON=squeeze(xg.phi(:,:,:))*180/pi;
x=MLON;
y=MLAT;
z=xg.alt;
z=z/10000; %To make things visually appealing, alt goes by 10km increments

if smallgrid %Make the grid small
    x=rescale(x);
    y=rescale(y);
    z=rescale(z);
end
tic

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

tarr=input('Enter timesteps as an array (ex: [2,9,21,30])');
%need not be uniformly spaced, but must be strictly increasing (I think)

inp = input('Add vector data [Y/N]? ','s');
if strcmp(inp,'Y')
    stop = 0;
    Vectors=struct; %Stores the data
    vectornames={}; %Stores the field names
    counter = 0; %keep track of how many times we've looped
    Vec_added=1; %Have we added some vectors?
    while stop==0
        counter=counter+1;
        [X] = input('Enter a list of 3 4-D arrays in format [X;Y;Z]');
        name=input('Name your vector: ','s');
        vectornames(counter)={name};
        Vectors=setfield(Vectors,name,X); %bad practice? yes. Does it work? yes.
        inp2=input('Continue [Y/N]?:','s');
        if strcmp(inp2,'Y')
            stop = 0; %continue
        elseif strcmp(inp2,'N')
            stop = 1; %break
        end
    end
else
    Vec_added=0; %we didn't add vectors
end

inp = input('Add scalar data [Y/N]? ','s');
if strcmp(inp,'Y')
    Scalar_added=1; %we added some scalars
    stop = 0;
    Scalars=struct; %store data
    scalarnames={}; %store names
    counter = 0; %keep track of iterations
    while stop==0
        counter=counter+1;
        X = input('Enter a 4-D array: ');
        name=input('Name your scalars: ','s');
        scalarnames(counter)={name};
        Scalars=setfield(Scalars,name,X); %hehe
        inp2=input('Continue [Y/N]?:','s');
        if strcmp(inp2,'Y')
            stop = 0; %continue
        elseif strcmp(inp2,'N')
            stop = 1; %break
        end
    end
else
    Scalar_added=0; %no scalars added
end

counter=0; %again, to keep track of iterations in the for loop
for i=tarr
    counter=counter+1;
    current_name=[folder,filesep,filename,'_',num2str(i),'.vtk']; %format the filename
    fid=fopen(current_name,'w');
    %Binary file header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'VTK from Matlab\n'); %I think you can change this? Don't really know why it's needed
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
            %in the case 2 format

            try
                if length(tarr)~=length(X(1,1,1,:))
                    Xdim = X(1:xlen,:,:,i);
                    Ydim = X(xlen+1:2*xlen,:,:,i);
                    Zdim = X(xlen*2+1:xlen*3,:,:,i);
                else
                    Xdim = X(1:xlen,:,:,counter);
                    Ydim = X(xlen+1:2*xlen,:,:,counter);
                    Zdim = X(xlen*2+1:xlen*3,:,:,counter);
                end
            catch 'If you are seeing this error message, make sure your vectors are added as [x;y;x], NOT [x,y,z]. The semicolon is not optional';
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

toc