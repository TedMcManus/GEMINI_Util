function [pars_meta] = run_series(ref,outputnames)
%% Creates a series of runs based on some reference simulation
% Prompts the user with relevant parameters, which can be changed by passing in
% an array of new values.
arguments
    ref (1,1) string %reference sim
    outputnames (1,:) string  %string array containing all the names of the runs we want to do
end
%EXAMPLE: [p_meta]=run_series('./mybaserun',["baserun_LOBG","baserun_HIBG","baserun_NOMINALBG"])
ref=char(ref);

num_runs=length(outputnames);

for i = 1:num_runs
    %Copy over the grid, config, and initial conditions
    confloc=ref;
    outfolder=char(outputnames(i));
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
            catch 
                error('Configuration files not found in reference directory');
            end
        end
    else
        error('Configuration files not found')
    end
end

%load the file with all the parameters stored in it (called "setup_meta")
try
    load([ref,'/setup_meta.mat'],'pars');
catch
    warning('setup file not found - running legacy version');
    name=input('Prefix for the _fields folder (eg OSSE in OSSE_fields): ','s');
    %If we don't find the setup_meta file, run a simpler "legacy" version of the run_series script
    %The legacy version lets the user change global parameters and scale/offset feautures
    run_series_legacy(ref,outputnames,name);
    return
end

%Print the "global" parameters (the parameters not tied to the specific input generator)
disp('------------------------------------------------------GLOBAL PARAMETERS------------------------------------------------------');
pars %#ok<NOPRT>
%Print the "map-specific" parameters (the parameters tied to the specific input generator)
disp('------------------------------------------------------MAP-SPECIFIC PARAMETERS------------------------------------------------------');
pars.(pars.arc)

%Make a big struct of structs, each struct in this big struct represents a run
pars_meta=struct;
for i = 1:num_runs
    pars_meta.(['run_',char(string(i))])=pars;
end

%Let the user skip particle or field data creation
F_or_P=input(['Do you want to copy the field or particle inputs to the new run unchanged? ' ...
    '\nWrite "fields" or "f" to copy field data and "particles" or "p" to copy particle data' ...
    '\nEnter to change both\n'],'s');

if strcmp(F_or_P,'f')||strcmp(F_or_P,'fields')||strcmp(F_or_P,'F')||strcmp(F_or_P,'Fields')
    noFields=true;
else
    noFields=false;
end

if strcmp(F_or_P,'p')||strcmp(F_or_P,'particles')||strcmp(F_or_P,'P')||strcmp(F_or_P,'Particles')
    noParticles=true;
else
    noParticles=false;
end

%Let the user alter values
%If the user is running multiple runs, the parameters for each run are entered in an array
glob=input('Do you want to alter a global parameter? [Y/N]\n','s');
if(strcmp(glob,'Y')||strcmp(glob,'y'))
    stop=0;
    while stop==0
        parm=input('Which global parameter?\n','S');
        val=input('To which new value (enter an array for multiple runs) \n');
        for i = 1:length(val)
            pars_meta.(['run_',char(string(i))]).(parm)=val(i);
        end
        inp2=input('Continue [Y/N]?:','s');
        if strcmp(inp2,'Y')
            stop = 0; %continue
        elseif strcmp(inp2,'N')
            stop = 1; %break
        end
    end
end

map_spec=input('Do you want to alter a map=specific parameter? [Y/N] \n','s');
if(strcmp(map_spec,'Y')||strcmp(map_spec,'y'))
    stop=0;
    while stop==0
        parm=input('Which map-specific parameter?\n','S');
        val=input('To which new value (enter an array for multiple runs) \n');
        for i = 1:length(val)
            pars_meta.(['run_',char(string(i))]).(pars_meta.(['run_',char(string(i))]).arc)=pars.(pars.arc);
            pars_meta.(['run_',char(string(i))]).(pars_meta.(['run_',char(string(i))]).arc).(parm)=val(i);
        end
        inp2=input('Continue [Y/N]?: ','s');
        if strcmp(inp2,'Y')
            stop = 0; %continue
        elseif strcmp(inp2,'N')
            stop = 1; %break
        end
    end
end

%Now run the automated run generator as many times as we need to
for i = 1:num_runs
    autocatalogue(ref,outputnames(i),pars_meta.(['run_',char(string(i))]),noFields,noParticles);
end
if strcmp(F_or_P,'f')||strcmp(F_or_P,'fields')||strcmp(F_or_P,'F')||strcmp(F_or_P,'Fields')
    warning('Overwriting fields')
    for i = 1:num_runs
        copyfile([ref,'/*_fields'],outputnames(i));
    end
end

if strcmp(F_or_P,'p')||strcmp(F_or_P,'particles')||strcmp(F_or_P,'P')||strcmp(F_or_P,'Particles')
    warning('Overwriting particles')
    for i = 1:num_runs
        copyfile([ref,'/*_fields'],outputnames(i));
    end
end