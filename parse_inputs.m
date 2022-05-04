function [pars] = parse_inputs(infile,outfolder,confloc)
%% Parse an text file and create a GEMINI input folder

arguments
    infile string % Path to the text file
    outfolder string % Path to folder that will hold GEMINI inputs. If it doesn't exist, we create it
    confloc char % Where is the grid and config
end

if ~exist(outfolder,'dir')
    mkdir(outfolder)
end

%Copy over the grid, config, and initial conditions
if exist(confloc,'dir')
    try
        copyfile([confloc,filesep,'simgrid.h5'],outfolder);
        copyfile([confloc,filesep,'simsize.h5'],outfolder);
        copyfile([confloc,filesep,'config.nml'],outfolder);
        copyfile([confloc,filesep,'initial_conditions.h5'],outfolder);
    catch 
        try
            copyfile([confloc,filesep,'inputs/simgrid.h5'],outfolder);
            copyfile([confloc,filesep,'inputs/simsize.h5'],outfolder);
            copyfile([confloc,filesep,'inputs/config.nml'],outfolder);
            copyfile([confloc,filesep,'inputs/initial_conditions.h5'],outfolder);
        catch error('Configuration files not found');
        end
    end
else
    error('Configuration files not found')
end

% We know this exists because we just made it
confdat=gemini3d.read.config(outfolder);

% Name the run according to the run name in the config file
p=char(confdat.prec_dir);
ind=find(p=='_');
runname=p(1:ind-1);


%Create table with the inputs
tmp=readtable(infile);
%Convert to arrays
params=tmp.Parameter;
vals=tmp.value;
%This is the output
pars=struct;

%Indexing - start to the index of 'GLOBALS' is the map-specific stuff, 
%and 'GLOBALS' to 'ENDFILE' is the global stuff
ind=find(ismember(params,'GLOBALS'));
endind=find(ismember(params,'ENDFILE'));

%The first parameter is the name of the arc. Create a substruct that is named 
%with the name of the arc (i.e, create "pars.gudda" for a gudda run)
arc=char(params(1));
pars.(arc)=struct;
pars.arc=arc;


%Make a structure with name-value pairs from the text file
for i = 2:ind-1 
    pars.(arc).(char(params(i)))=vals(i);
end

for i = ind:endind-1
    pars.(char(params(i)))=vals(i);
end

%Take out the useless stuff
pars=rmfield(pars,'GLOBALS');

%You can change this if you want
pars.runname=runname;

%X and Y need to be in an array. Do that. 
pars.xy=[pars.x,pars.y];
pars=rmfield(pars,'x');
pars=rmfield(pars,'y');

%Now save the ouptuts
if ~exist(outfolder,'dir')
    mkdir(outfolder)
end
save([char(outfolder),'/setup_meta.mat'],"pars");

%Now run the automated input generator to make the field and particle inputs
autocatalogue(outfolder,outfolder,pars,false,false);
