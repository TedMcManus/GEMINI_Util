function [out] = loadvar(direc,times,vars,xg,outputfile)

arguments
    direc char %location of data files
    times double %times you want to load
    vars string %variables you want to load
    xg (1,1) struct %the grid
    outputfile char = 'none' %Do you want to save the output to a separate .mat file?
end
%Times is an increasing array of doubles; [1,4,7,30] works
%Vars is a string array that contains the names of variables you want. 
%OPTIONS- density,temp/temperature,flow,current


startdir=pwd;
%go to the right place

cd(direc);

lx1=xg.lx(1);
lx2=xg.lx(2);
lx3=xg.lx(3);

is_density=any(contains(vars,"density"));
if is_density
    ne=zeros(lx1,lx2,lx3,length(times));
end

is_temp=any([contains(vars,"temperature"),contains(vars,"temp")]);
if is_temp
    Te=zeros(lx1,lx2,lx3,length(times));
end

is_flow=any([contains(vars,"flow")]);
if is_flow
    V1=zeros(lx1,lx2,lx3,length(times));
    V2=V1;
    V3=V1;
end

is_current=any([contains(vars,"current")]);
if is_current
    J1=zeros(lx1,lx2,lx3,length(times));
    J2=V1;
    J3=V1;
end


files=dir('./*_fields/fields.mat'); 
load([files(1).folder,filesep,files(1).name],'E_times') %Find and load time data
counter=0; %iteration number
for i = times %timestep number
    
    %Get the filename
    frame=gemini3d.datelab(E_times(i));
    frame=[char(frame),'.h5'];
    counter=counter+1;
    
    %Read all the data we want
    if is_density
        tmpvar=h5read(frame,'/nsall');
        ne(:,:,:,counter)=tmpvar(:,:,:,7);
    end
    if is_temp
        tmpvar=h5read(frame,'/Tsall');
        Te(:,:,:,counter)=tmpvar(:,:,:,7);
    end
    if is_flow
        tmpvar=h5read(frame,'/vs1all');
        V1(:,:,:,counter)=mean(tmpvar,4);
        tmpvar=h5read(frame,'/v2avgall');
        V2(:,:,:,counter)=tmpvar(:,:,:);
        tmpvar=h5read(frame,'/v3avgall');
        V3(:,:,:,counter)=tmpvar(:,:,:);
    end
    if is_current
        tmpvar=h5read(frame,'/J1all');
        J1(:,:,:,counter)=tmpvar;
        tmpvar=h5read(frame,'/J2all');
        J2(:,:,:,counter)=tmpvar(:,:,:);
        tmpvar=h5read(frame,'/J3all');
        J3(:,:,:,counter)=tmpvar(:,:,:);
    end
    
end

%You don't see the setfield use
out=struct();
if exist('ne','var')
    out=setfield(out,"ne",ne);
end

if exist('Te','var')
    out=setfield(out,"Te",Te);
end

if exist('V1','var')
    out=setfield(out,"V1",V1);
    out=setfield(out,"V2",V2);
    out=setfield(out,"V3",V3);
end

if exist('J1','var')
    out=setfield(out,"J1",J1);
    out=setfield(out,"J2",J2);
    out=setfield(out,"J3",J3);
end

out.times=times;

cd(startdir);
if ~strcmp(outputfile,'none')
    save(outputfile,'out','-v7.3'); %save the data to a HUGE .mat file
end
end