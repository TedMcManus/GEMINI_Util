function [E_meta,P_meta] = run_series_legacy(ref,outputnames,name)
%% Permutes a set of GEMINI inputs to create a sequence of simulations specified by the user
% This funciton is a liiiitle gnarly b/c of how often it uses dynamic field names... 
% But I will do my best :)

arguments
    ref (1,1) string %reference sim
    outputnames (1,:) string  %string array containing all the names of the runs we want to do
    name (1,1) string %prefix for the folders (eg STEVE in STEVE_fields)
end

ref=char(ref);
name=char(name);
num_runs=length(outputnames);

% Check if the user wants to permute the field inputs
Fields=input('Alter field inputs? [Y/N]\n','s');

if strcmp(Fields,'y')||strcmp(Fields,'Y') %then the user wants to alter fields... here we go

    %Find the field data input folder and load its data
    files= dir([ref,'/*_fields/*.mat']);
    load([files.folder,filesep,files.name]);

    %load the grid and whatnot
    xg=gemini3d.read.grid(ref);
    gemgrid;
    llon=xg.lx(2);
    llat=xg.lx(3);

    %read a config file and set up time indexing
    confdat=gemini3d.read.config(ref);
    ymd0=confdat.ymd;
    UTsec0=confdat.UTsec0;
    tdur=confdat.tdur;
    dtout=confdat.dtout;
    flagoutput=confdat.flagoutput;

    dt = dtout;
    tmin=0;
    tmax=tdur;
    times=tmin:dt:tmax; %this is the array of field data timesteps relative to the first timestep
    lt=length(times);

    % Create a nice little array of datetime objects b/c we might need it
    ymd=confdat.ymd;
    UTsec=UTsec0+times;     %time given in file is the seconds from beginning of hour
    UThrs=UTsec/3600;
    expdate=cat(2,repmat(ymd,[lt,1]),UThrs(:),zeros(lt,1),zeros(lt,1));
    hours=mod(expdate(:,4),1);
    expdate(:,4)=expdate(:,4)-hours;
    expdate(:,5)=hours.*60;
    minute=mod(expdate(:,5),1);
    expdate(:,5)=expdate(:,5)-minute;
    expdate(:,6)=minute.*60;
    dates=datetime(int32(expdate));
    names=gemini3d.datelab(dates);

    % Read out the flagdirich variable from an output file
    % This tells GEMINI whether we're using potential or FAC as an input
    files= dir([ref,'/*_fields/*00.h5']);
    flagdirich=h5read([files(end).folder,filesep,files(end).name],'/flagdirich');
    flagdirich=ones(301).*double(flagdirich);

    %Create the struct that the input writing scripts need to create an input folder
    E=struct();
    E.flagdirich=flagdirich;
    E.Exit=Exit;
    E.Eyit=Eyit;
    E.Vminx1it=Vminx1it;
    E.Vmaxx1it=Vmaxx1it;
    E.Vminx2ist=Vminx2ist;
    E.Vmaxx2ist=Vmaxx2ist;
    E.Vminx3ist=Vminx3ist;
    E.Vmaxx3ist=Vmaxx3ist;
    E.mlon=double(mlon);
    E.mlat=double(mlat);
    E.llon=double(llon);
    E.llat=double(llat);
    E.times=E_times;

    %Here's a fun one: E_meta is a struct of structs 
    %for every simulation we wanna do, we're gonna get a struct inside this struct
    E_meta=struct;
    for i = 1:num_runs
        E_meta.(['run_',char(string(i))])=E; %put an E inside the meta E
    end

    % Let the user change some stuff
    choice = input(['\nAvailable options for field data:\n' ...
        '1: Change background fields\n' ...
        '2: Scale/offset potential or FAC\n' ...
        '3: Both\n'])
    %If the user wants to change the background values
    if choice==1||choice==3
        %show them what the values were
        disp(['ExBG: ',num2str(Exit(1,1)),' EyBG: ',num2str(Eyit(1,1))]);
        %ask them for new values
        Ex=input('\nArray of new ExBG values: ');
        Ey=input('\nArray of new EyBG values: ');
        %and write them into the E_meta struct
        for i = 1:length(Ex)
            E_meta.(['run_',char(string(i))]).('Exit')=ones(size(E.Exit)).*Ex(i);
            E_meta.(['run_',char(string(i))]).('Eyit')=ones(size(E.Exit)).*Ey(i);
        end
    end
    
    %if the user wants to scale/offset the values, let them!
    if choice==2||choice==3
        mult=input('\nMultiply the drivers by the array of values: ');
        add=input('\nAnd add the array: ');
        for i = 1:length(add)
            E_meta.(['run_',char(string(i))]).Vmaxx1it=(E.Vmaxx1it.*mult(i))+add(i); %and do the math while we store it into the big struct
        end
    end
    disp('field inputs being copied')
    for i = 1:num_runs
        %now, for every run we want to do, write out the data
        mkdir([char(outputnames(i)),filesep,char(name),'_fields']);
        gemini3d.write.Efield(E_meta.(['run_',char(string(i))]),[char(outputnames(i)),filesep,char(name),'_fields'],'h5');
        Exit=E_meta.(['run_',char(string(i))]).Exit;
        Eyit=E_meta.(['run_',char(string(i))]).Eyit;
        Vmaxx1it=E_meta.(['run_',char(string(i))]).Vmaxx1it;
        save([char(outputnames(i)),filesep,char(name),'_fields','/fields.mat'],'mlon','mlat','MLAT','MLON','Exit','Eyit','Vminx*','Vmax*','E_times');
    end
else %just copy the inputs over 
    disp('field inputs being copied')
    for i = 1:num_runs
        copyfile([ref,'/*_fields'],outputnames(i));
    end
end

% and now we go on to the particle inputs
Particles=input('Alter particle inputs? [Y/N]\n','s');
if strcmp(Particles,'y')||strcmp(Particles,'Y')
    %if the user wants to alter particles... let's do it!

    %Load in all the reference stuff from the "ref" folder
    files= dir([ref,'/*_particles/*.mat']);
    load([files.folder,filesep,files.name]);
    xg=gemini3d.read.grid(ref);
    gemgrid;
    llon=xg.lx(2);
    llat=xg.lx(3);
    P_meta=struct;
    P=struct;
    P.mlon=double(mlon);
    P.mlat=double(mlat);
    P.llon=double(llon);
    P.llat=double(llat);
    P.E0it=E0it;
    P.Qit=Qit;
    P.times=P_times;

    %Create a meta struct, as above for the fields
    for i = 1:num_runs
        P_meta.(['run_',char(string(i))])=P;
    end
    % Let the user know what they can do 
    choice = input(['\nAvailable options for particle data:\n' ...
        '1: Change background Q and E0 (retaining values above background)\n' ...
        '2: Completely flatten the Q and E0 to background values\n']);

    if choice==1 %the user wants to change the background values
        disp(['E0BG: ',num2str(min(E0it(:))),' QBG: ',num2str(min(Qit(:)))]);
        EBG=input('\nArray of new E0BG values: ');
        QBG=input('\nArray of new QBG values: ');

        % Go over every sim in the series of sims and alter the background precip to the specified values
        for i = 1:length(EBG)
            E0_tmp=E0it;
            E0_tmp(E0_tmp==min(E0_tmp(:)))=EBG(i);
            E0_tmp(E0_tmp<EBG(i))=EBG(i);
            Q_tmp=Qit;
            Q_tmp(Q_tmp==min(Q_tmp(:)))=QBG(i);
            Q_tmp(Q_tmp<QBG(i))=QBG(i);
            P_meta.(['run_',char(string(i))]).('E0it')=E0_tmp;
            P_meta.(['run_',char(string(i))]).('Qit')=Q_tmp;
        end
    end

    if choice==2 %the user wants to flatten the arrays to a background value
        disp(['E0BG: ',num2str(min(E0it(:))),' QBG: ',num2str(min(Qit(:)))]);
        EBG=input('\nArray of new E0BG values: ');
        QBG=input('\nArray of new QBG values: ');
        % erase all the data and flatten it out
        for i = 1:length(EBG)
            P_meta.(['run_',char(string(i))]).('E0it')=ones(size(E0it)).*EBG(i);
            P_meta.(['run_',char(string(i))]).('Qit')=ones(size(Qit)).*QBG(i);
        end
    end
    for i = 1:num_runs
        %now, for every run we want to do, write out the data
        mkdir([char(outputnames(i)),filesep,char(name),'_particles']);
        gemini3d.write.precip(P_meta.(['run_',char(string(i))]),[char(outputnames(i)),filesep,char(name),'_particles'],'h5');
        E0it=P_meta.(['run_',char(string(i))]).E0it;
        Qit=P_meta.(['run_',char(string(i))]).Qit;
        P=P_meta.(['run_',char(string(i))]);
        save([char(outputnames(i)),filesep,char(name),'_particles','/particles.mat'],'mlon','mlat','MLAT','MLON','E0it','Qit','P_times','P');
    end
end
