%% This script is for plotting the average wheel trajectory aligned to event times (i.e. movement onset)

% fpath = '/Volumes/My Passport for Mac/171013/behaveNewUnits';
%Open the relavent figure for the overlay, then run this:
figure;
hold on;
% yFactor = 1e-4;  %to position the wheel data on the y axis of existing figure
yFactor = 1;       %to position the wheel data on the y axis of new figure

fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/171013';
filestr = '171013002.ns4';
% behavFile = strcat(fpath, '/', '171013_bhCue4.mat'); nope
% behavFile = strcat(fpath, '/', '171013_bhCue3.mat'); nope
% behavFile = strcat(fpath, '/', '171013_bhCue5.mat'); nope
% behavFile = strcat(fpath, '/', '171013_bhCue6.mat'); nope
behavFile = strcat(fpath, '/', '171013_bhCueNS4'); %variable, but aligned the best
% behavFile = strcat(fpath, '/', '171013_bhAllNS4.mat'); %variable, but aligned the best

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/171012/newBehaveUnits';
% filestr = '171012001.ns4';
% behavFile = strcat(fpath, '/', '171012_bhCue.mat');

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170118/behaveNewUnits';
% filestr = '170118002.ns4';
% % behavFile = strcat(fpath, '/', '170118_bhCue2.mat');
% behavFile = strcat(fpath, '/', '170118_bhCueMv.mat');


% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170111';
% filestr = '170111002.ns4';
% % behavFile = strcat(fpath, '/', '170111_bhCue.mat');
% behavFile = strcat(fpath, '/', '170111_bhCAllNS4.mat');

% behavFile = strcat(fpath, '/', '170111_bhCue2.mat');
% behavFile = strcat(fpath, '/', '170111_bhCueMv.mat');

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170112/behaveNew';
% filestr = '170112002.ns4';
% behavFile = strcat(fpath, '/', '170112_bhCue.mat'); %works! so why no mvmnt neural activity?

% filestr2 = strcat(fpath, '/', filestr);
% load(behavFile)

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/newBehaveUnits/behaveChunks';
% % cd fpath;
% % fname = strcat(fpath, '/', '151105_all_pop.mat');                          
% % behavFile = strcat(fpath, '/', '151105_all_bh.mat');
% fname = strcat(fpath, '/', '151105_update_pop.mat');                          
% behavFile = strcat(fpath, '/', '151105_update_bh.mat');
% filestr = '151105tag1001.ns4';

filestr2 = strcat(fpath, '/', filestr);
load(behavFile)

[data] = openNSx(filestr2);              %MUST LOAD THE openNSx.m file from 2016 for this to work (currently in forTonic_v2 folder, mac HD)

% LrewTrialStarts  = ContData.behavior.LrewTrialStarts;  %for 151105
% RrewTrialStarts  = ContData.behavior.RrewTrialStarts;  %for 151105
% eventTimesIpsi = LrewTrialStarts; %for 151105
% eventTimesContra = RrewTrialStarts; %for 151105

% 
% fastCueIpsi = export_mvmt_data_for_psths.fastCueL;
% fastCueContra = export_mvmt_data_for_psths.fastCueR;
% slowCueIpsi = export_mvmt_data_for_psths.slowCueL;
% slowCueContra = export_mvmt_data_for_psths.slowCueR;
%
%% From bhAllNS4.mat (uses CSV-aligned to NS4 data for trial starts; uses CSV for reaction times
% All_L_trialStarts = export_mvmt_data_for_psths.All_L_trialStarts;  %for 170111 & 170118, difficult sessions
% All_R_trialStarts = export_mvmt_data_for_psths.All_R_trialStarts;

% All_L_trialStarts_valid = export_mvmt_data_for_psths.All_L_trialStarts_valid;  %for 170111 & 170118, difficult sessions
% All_R_trialStarts_valid = export_mvmt_data_for_psths.All_R_trialStarts_valid;

% eventTimesIpsi = All_L_trialStarts;
% eventTimesContra = All_R_trialStarts;
 
% LrewTrueTrialStarts = export_mvmt_data_for_psths.LrewTrueTrialStarts;
% RrewTrueTrialStarts = export_mvmt_data_for_psths.RrewTrueTrialStarts;

%% from bhCueNS4.mat - uses the NS4 file only for trial start
Lcue = export_mvmt_data_for_psths.Lcue;  %for NS4 aligned only (171012 & 171013) from bhCueNS4.mat
Rcue = export_mvmt_data_for_psths.Rcue;

%For all trial starts
eventTimesIpsi = Lcue;
eventTimesContra = Rcue;


%For correct trial only starts:
% eventTimesIpsi = LrewTrueTrialStarts;
% eventTimesContra = RrewTrueTrialStarts;

% eventTimesIpsi = slowCueIpsi;
% eventTimesContra = slowCueContra;

% eventWin = [1e3 2e3];
eventWin = [2e3 2e3];

%% Continuous data:

% chan.wheel = 131;                                   
% wheelMvs1 = find(data.MetaTags.ChannelID==chan.wheel);
% wheelMvs2 = double(data.Data(wheelMvs1,:)); 
% wheelMvs = decimate(wheelMvs2,10); 

chan.start = 129;                                   
start1 = find(data.MetaTags.ChannelID==chan.start);
start2 = double(data.Data(start1,:)); 
starts = decimate(start2,10); 

window = eventWin;
indicesIpsi = eventTimesIpsi %(12:13);
indicesContra = eventTimesContra %(12:13);
source = starts;

[sinkIpsi] = TNC_ExtTrigWins(source,indicesIpsi,window);  %plot the ipsi wheel movements
[sinkContra] = TNC_ExtTrigWins(source,indicesContra,window);
% wheel_postSgolay = sgolayfilt(source, 9, 101);

% figure; hold on;
% plot(sinkIpsi.avg,'-','LineWidth',2); title('171013')  
sinkIpsi.avg = sinkIpsi.avg * yFactor;
sinkIpsi.err = sinkIpsi.err * yFactor;
%     shadedErrorBar(-window(1):window(2), sinkIpsi.avg, sinkIpsi.err,'bl');  %blue for left; for 151105 dataset = ipsi
%     xlabel('ms', 'FontSize', 24); ylabel('wheel velocity', 'FontSize', 24);
%     xlabel('ms'); ylabel('wheel velocity');
    plot(-window(1):window(2), sinkIpsi.avg, 'bl');
    %     legend('ipsi = blue');
            ax = gca; 
%             ax.FontSize = 24;


% plot(sinkContra.avg, '-','LineWidth',2);
% figure; hold on;
sinkContra.avg = sinkContra.avg * yFactor;
sinkContra.err = sinkContra.err * yFactor;

%     shadedErrorBar(-window(1):window(2), sinkContra.avg, sinkContra.err,'r');  %blue for left; for 151105 dataset = ipsi
    plot(-window(1):window(2), sinkContra.avg, 'r');

% xlabel('ms', 'FontSize', 24); ylabel('wheel velocity', 'FontSize', 24);
%  xlabel('ms'); ylabel('wheel velocity');

% legend('ipsi = blue');
            ax = gca; 
%             ax.FontSize = 24;

% source is the continuous data.
% indices are movement stamps
% window is [1000 2000]