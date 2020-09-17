clear all
close all
clc

addpath(fullfile('..','src','common'))
addpath(fullfile('..','src','lib'))
addpath(fullfile('..','src','model'))
addpath(fullfile('..','plot_dir'))

path2root = '..';

% set here the path where you want to save the figures produced by the
% scripts
% path2fig = '\\oscar\macintosh hd\Users\lucas\Documents\EPFL\PDM\fig';
path2fig = fullfile(path2root,'fig');

save(fullfile(path2root,'paths'),'path2fig')