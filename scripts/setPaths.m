clear all
close all
clc

path2root = '..';

addpath(fullfile(path2root,'src','common'))
addpath(fullfile(path2root,'src','lib'))
addpath(fullfile(path2root,'src','model'))
addpath(fullfile(path2root,'plot_dir'))

path2res = '\\sti1files.epfl.ch\unfold\unfold-commun\temp_lucas';
% path2res = fullfile(path2root,'data','paramsweep');

% set here the path where you want to save the figures produced by the
% scripts
path2fig = '\\oscar\macintosh hd\Users\lucas\Documents\EPFL\PDM\fig';
% path2fig = fullfile(path2root,'fig');

if ~isfolder(path2fig)
    mkdir(path2fig)
end

save(fullfile(path2root,'paths'),'path2fig','path2root','path2res')