%% Script to load the standard head-fixed behavior rig recording session - Beth's local version

%% Basic data loading
chunkDuration = sprintf('t:%g:%g',0,300); %300 seconds is tagging recording time (and same time clustered)

% Load the continuous data recording channels. Examples:
% filenamestr = '/Users/dudmanj/Desktop/m42/m42-20110427-wide-plustim-008.ns5'
% dataSpks = openNSx('report','read',filenamestr, 't:0:300', 'sec');

% filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/150529ipsi/150529_3004.ns5';
filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/test/BethTest001.ns5';
% filenamestr = '/Users/stubblefielde/Desktop/m42/m42-20110427-wide-plustim-008.ns5'
dataSpks = openNSxBS('report','read',filenamestr, 't:0:300', 'sec', 'p:double');

% Load the continuous behavior monitor channels
% filenamestr = '/Users/dudmanj/Desktop/m42/m42-20110427-wide-plustim-008.ns4'
% dataAinp = openNSx('report','read',filenamestr,'e:137:141','t:0:300', 'sec');
% filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/150529ipsi/150529_3004.ns4';
filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/test/BethTest001.ns4';
% filenamestr = '/Users/stubblefielde/Desktop/m42/m42-20110427-wide-plustim-008.ns4'
dataAinp = openNSxBS('report','read',filenamestr,'e:1:128','t:0:300', 'sec', 'p:double');

% Load behavior event times and digital stamps
% filenamestr = '/Users/dudmanj/Desktop/m42/m42-20110427-wide-plustim-008.nev'
% dataEvents = openNEV(filenamestr,'report','read','nomat');
% filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/150529ipsi/150529_3004.nev';
filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/test/BethTest001.nev';
% filenamestr = '/Users/stubblefielde/Desktop/m42/m42-20110427-wide-plustim-008.nev'
dataEvents = openNEVBS(filenamestr,'report','read','nomat');

openNSxBS('report','read','c:\data\sample.ns5', 'e:15:30', 't:3:10', 'min', 'p:int16');