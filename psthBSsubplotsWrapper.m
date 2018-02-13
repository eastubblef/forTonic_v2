function [] = psthBSsubplotsWrapper(fname, behavFile, units, alignVar)
%% This mfile calls psthBSsubplots

fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/laserChunks/untagged';
fname = strcat(fpath, '/', '151105_untagged_pop.mat');                           
behavFile = strcat(fpath, '/', '151105_untagged_bh.mat');

load(fname);
pop = PopData;
units = 3:length(pop.session.unit);
% varargin = units;

load(behavFile);
behavior = ContData.behavior;

 wheelsLfirstValid = behavior.wheelsLfirstValid;
 wheelsRfirstValid = behavior.wheelsRfirstValid;
 wheelsILfirstValid = behavior.wheelsILfirstValid;
 wheelsIRfirstValid = behavior.wheelsIRfirstValid;
 rewTimesValid = behavior.rewTimesValid;

 alignVar = wheelsLfirstValid;
 alignVar2 = wheelsRfirstValid;

 
 [] = psthBSsubplots2(fname, behavFile, units, alignVar)

