clear all
close all
clc

path2root = '..';

addpath(fullfile(path2root,'src','common'))
addpath(fullfile(path2root,'src','lib'))
addpath(fullfile(path2root,'src','model'))
addpath(fullfile(path2root,'plot_dir'))

%% Path to the result mat-files 
% needs to be set for plotSweepResults.m to work correctly.
% /!\ You should not use a remote location here, the files are too large.
% Instead, download the res.mat on your computer and indicate the path to
% the local copy. Alternatively, you can create a new res.mat by running
% paramsweep.m 

path2res = fullfile(path2root,'data','paramsweep');

%% Path to the figure folder for saving figures
% set here the path where you want to save the figures produced by the
% scripts
if ispc()
% path2fig = '\\oscar\macintosh hd\Users\lucas\Documents\EPFL\PDM\fig';
end
path2fig = fullfile(path2root,'fig');

if ~isfolder(path2fig)
    mkdir(path2fig)
end

save(fullfile(path2root,'paths'),'path2fig','path2root','path2res')