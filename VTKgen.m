%This whole thing is like a pseudo-function but it's written as a
%script to allow it to access the base workspace, which should 
%contain an xg and all the stuff you want to load. This is a script for doing a 
%fast&dirty write of a single timestep to VTK. For importing multiple
%timesteps, VTK_w_time provides a more elegant solution

filename=input('Output file name: ','s');
%Name of the output file

if exist('Jfac_all','var')
    JPar = Jfac_all;
    JH = Normalize(JH_all);
    JP = Normalize(JP_all);
end

%% Set up the grid
MLAT=90-squeeze(xg.theta(:,:,:))*180/pi;
MLON=squeeze(xg.phi(:,:,:))*180/pi;
x=MLON;
y=MLAT;
z=xg.alt;
z=z/10000;
%This part makes the grid square from 0 to 1. 
%It could be disabled but it gets kinda weird to work with 
%in paraview and I'm lazy 
x=rescale(x);
y=rescale(y);
z=rescale(z);
tic

nr_of_elements=numel(x);
fid = fopen(filename, 'w'); 

%ASCII file header
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'BINARY\n\n');
fprintf(fid, 'DATASET STRUCTURED_GRID\n');
fprintf(fid, ['DIMENSIONS ' num2str(size(x,1)) ' ' num2str(size(x,2)) ' ' num2str(size(x,3)) '\n']);
fprintf(fid, ['POINTS ' num2str(nr_of_elements) ' float\n']);
fclose(fid);

%append binary x,y,z grid data
fid = fopen(filename, 'a'); 
fwrite(fid, [reshape(x,1,nr_of_elements);  reshape(y,1,nr_of_elements); reshape(z,1,nr_of_elements)],'float','b');
fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements)]);

inp = input('Add vector data? \n [Y/N]','s');
if strcmp(inp,'Y')
    i = 0;
    while i == 0
        [X] = input('Enter a list of 3 arrays in format [X;Y;Z] \n These will define the data\n');
        %MAKE SURE you use a semicolon for this to work
        %Now split the big array into its X, Y,and Z components
        xlen = length(x(:,1,1));
        Xdim = X(1:xlen,:,:);
        Ydim = X(xlen+1:2*xlen,:,:);
        Zdim = X(xlen*2+1:xlen*3,:,:);
        inp3 = input('Name your vector field! \n','s');
        fprintf(fid, ['\n VECTORS ', inp3, ' float\n']); %binary header
        %What follows is how VTK likes its vector data
        fwrite(fid, [reshape(Xdim,1,nr_of_elements);  reshape(Ydim,1,nr_of_elements); reshape(Zdim,1,nr_of_elements)],'float','b'); %binary data
        inp2 = input('Add more data? [Y/N]\n','s');
        if strcmp(inp2,'Y')
            i = 0; %continue
        elseif strcmp(inp2,'N')
            i = 1; %break
        end
    end
end

inp = input('Add scalar data? \n [Y/N]\n Enter ALL to load all timesteps \n','s');
if strcmp(inp,'Y')
    i = 0;
    while i == 0
        X = input('Enter an array whose dimensions are consistent with the grid \n');
        oof = false; %The array is the wrong size
        try OutArr = reshape(X,1,nr_of_elements);
        catch('Matrix dimensions inconsistent with grid');
            disp('');
            oof = true; %Let the user try again
        end
        if ~oof
            inp3 = input('Name your vector field! \n','s');
            %Scalar binary header is two lines, you can specify a lookup
            %table if you want, I don't
            fprintf(fid, ['\nSCALARS ',inp3,' float\n']); 
            fprintf(fid, 'LOOKUP_TABLE default\n'); 
            %big ol' line of numbers
            fwrite (fid, OutArr,'float','b'); %binary data
            inp2 = input('Add more data? [Y/N]\n','s');
            if strcmp(inp2,'Y')
                i = 0; %continue
            elseif strcmp(inp2,'N')
                i = 1; %break
            end
        end
    end
end


fclose(fid);
toc