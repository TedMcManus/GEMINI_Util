function [output_struct] = load_multirun(runlist)
%% Loads a slice from the keogram data for a series of runs (assuming they have keogram plots)
% this is useful for doing comparison plots

tic

%This will store all the data in the end
output_struct=struct();

for i = 1:length(runlist)
    dataloc=char(runlist(i));
    if strcmp(dataloc(1:2),'./')
        dataloc=dataloc(3:end);
    end
    if exist([dataloc,filesep,'keogram_data.mat'],'file') 
        load([dataloc,filesep,'keogram_data.mat']);
        %Use dynamic field indexing to put all the data from keogram_data.mat into a struct
        output_struct.(dataloc).('dens_95')=dens_95;
        output_struct.(dataloc).('dens_120')=dens_120;
        output_struct.(dataloc).('dens_300')=dens_300;
        output_struct.(dataloc).('dens_800')=dens_800;
        output_struct.(dataloc).('V_2')=V_2;
        output_struct.(dataloc).('V_3')=V_3;
        output_struct.(dataloc).('Jpar')=J_par;
        output_struct.(dataloc).('SigmaH')=SigmaH;
        output_struct.(dataloc).('SigmaP')=SigmaP;
    end
end

toc %output is always appreciated - it lets the user know how much time they're wasting


