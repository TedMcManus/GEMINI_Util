%% Example script for making scatter plots

% The folders that contain the data
fldrs=["folder1","folder2","folder3"];

% The overall title
bigtitle="This is an overall title interpreted in $\LaTeX$";

% The title of the legend
lgdtitle="This is the legend title interpreted in $\LaTeX$";

% Names of the legend
LNames = ["Run 1","Run 2","Run 3"];

% Icons for each scatter
icons=[".","x","o"];

% Where is the grid
xg=gemini3d.read.grid('path_to_grid');

% Some garbage to make things work
s=size(icons);
s=s(2);
tmp="";
ind=0;
for i= 0:s-1
    ind=ind+1;
    tmp(ind)=icons(s-i);
end
icons=tmp;

tmp="";
ind=0;
for i= 0:s-1
    ind=ind+1;
    tmp(ind)=LNames(s-i);
end
LNames=tmp;

%% This part is the section that actually runs the scatter plot

% TO SAVE ALL THE DATA TO A FILE NAMED './SCATTER_DATA.MAT'
% plot_scatters(fldrs,xg,bigtitle,lgdtitle,LNames,'save',null(1),'./scatter_data',icons,"")

% TO LOAD ALL THE DATA FROM A FILE NAMED './SCATTER_DATA.MAT'
% plot_scatters(fldrs,xg,bigtitle,lgdtitle,LNames,'load',null(1),'./scatter_data',icons,"")

% TO LOAD ALL THE DATA FROM A STRUCTURE IN YOUR WORKSPACE CALLED "SAVE_STRUCT"
% plot_scatters(fldrs,xg,bigtitle,lgdtitle,LNames,'load',save_struct,'./nothing',icons,"")


%% Anatomy of the function
%plot_scatters(
% fldrs, (the folders)

% xg, (the grid)

% bigtitle, (the overall title)

% lgdtitle, (title of the legend)

% LNames, (names of each dot on the plot)

% savecase, ('load' will attempt to load data from the "dataloc" file or the "data" positional argument.
%            'save' will attempt to save data to a file called 'dataloc.mat')

% data, (a structure with all the relevant data. This will be created by running plot_scatter with the 'load'
%        option and it will be saved to a file called "dataloc")  

% dataloc, (either the place to load data from or the place to save data to)
% icons, (the type of dot on the scatter for each run)

% colors ('auto' will automatically do colors. You can also pass a string array with a color for each run
%         individually)
%)

