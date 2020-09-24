clear all
close all
clc

addpath(fullfile('..','src','common'))
addpath(fullfile('..','src','lib'))
addpath(fullfile('..','src','model'))
addpath(fullfile('..','plot_dir'))

path2root = '..';
path2res = '\\sti1files.epfl.ch\unfold\unfold-commun\temp_lucas';
% set here the path where you want to save the figures produced by the
% scripts
% path2fig = '\\oscar\macintosh hd\Users\lucas\Documents\EPFL\PDM\fig';
path2fig = fullfile(path2root,'fig');

if ~isfolder(path2fig)
    mkdir(path2fig)
end

save(fullfile(path2root,'paths'),'path2fig','path2root','path2res')