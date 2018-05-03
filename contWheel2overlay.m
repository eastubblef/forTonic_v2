%% This script is for plotting the average wheel trajectory aligned to event times (i.e. movement onset)

% fpath = '/Volumes/My Passport for Mac/171013/behaveNewUnits';
%Open the relavent figure for the overlay, then run this:
hold on;
yFactor = 1e-4;  %to position the wheel data on the y axis

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/171013';
% filestr = '171013002.ns4';
% % behavFile = strcat(fpath, '/', '171013_bh3Rxn.mat');
% behavFile = strcat(fpath, '/', '171013_bh5Rxn.mat');

fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170111';
filestr = '170111002.ns4';
behavFile = strcat(fpath, '/', '170111_bh5Rxn.mat'); %works without offset


filestr2 = strcat(fpath, '/', filestr);
load(behavFile)
% 
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105';
% % cd fpath;
% fname = strcat(fpath, '/', '151105_all_pop.mat');                          
% behavFile = strcat(fpath, '/', '151105_all_bh.mat');
% filestr = '151105tag1001.ns4';
% filestr2 = strcat('/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/newBehaveUnits/behaveChunks', '/', filestr);
% load(behavFile)

[data] = openNSx(filestr2);              %MUST LOAD THE openNSx.m file from 2016 for this to work (currently in forTonic_v2 folder, mac HD)

% wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;  %2017 data
% wheelsRfirstValid = export_mvmt_data_for_psths.wheelRfirstValid;
wheelALL_LfirstValid = export_mvmt_data_for_psths.wheelALL_LfirstValid;
wheelALL_RfirstValid = export_mvmt_data_for_psths.wheelALL_RfirstValid;

eventTimesIpsi = wheelALL_LfirstValid;
eventTimesContra = wheelALL_RfirstValid;

% eventWin = [1e3 2e3];
eventWin = [2e3 2e3];

chan.wheel = 131;                                   
wheelMvs1 = find(data.MetaTags.ChannelID==chan.wheel);
wheelMvs2 = double(data.Data(wheelMvs1,:)); 
wheelMvs = decimate(wheelMvs2,10); 
% figure; hold on; plot(wheelMvs, 'k'); xlim([1e5 1.5e5]);    %longer trace snippet filtered

wheel_golay = sgolayfilt(wheelMvs, 9, 101);
% figure; hold on; plot(wheel_golay, 'k'); xlim([1e5 1.5e5]); %longer trace snippet filtered

% wheel_golay = sgolayfilt(wheelMvs, 9, 101);
sourceIpsi = eventTimesIpsi;
sourceContra = eventTimesContra;
% source = wheel_golay;

window = eventWin;
indicesIpsi = wheelALL_LfirstValid;
indicesContra = wheelALL_RfirstValid;
source = wheelMvs;

[sinkIpsi] = TNC_ExtTrigWins(source,indicesIpsi,window);  %plot the ipsi wheel movements
[sinkContra] = TNC_ExtTrigWins(source,indicesContra,window);
% wheel_postSgolay = sgolayfilt(source, 9, 101);

% figure; hold on;
% plot(sinkIpsi.avg,'-','LineWidth',2); title('171013')  
sinkIpsi.avg = sinkIpsi.avg * yFactor;
sinkIpsi.err = sinkIpsi.err * yFactor;
    shadedErrorBar(-window(1):window(2), sinkIpsi.avg, sinkIpsi.err,'bl');  %blue for left; for 151105 dataset = ipsi
%     xlabel('ms', 'FontSize', 24); ylabel('wheel velocity', 'FontSize', 24);
%     xlabel('ms'); ylabel('wheel velocity');

    %     legend('ipsi = blue');
            ax = gca; 
%             ax.FontSize = 24;


% plot(sinkContra.avg, '-','LineWidth',2);
% figure; hold on;
sinkContra.avg = sinkContra.avg * yFactor;
sinkContra.err = sinkContra.err * yFactor;

    shadedErrorBar(-window(1):window(2), sinkContra.avg, sinkContra.err,'r');  %blue for left; for 151105 dataset = ipsi

% xlabel('ms', 'FontSize', 24); ylabel('wheel velocity', 'FontSize', 24);
%  xlabel('ms'); ylabel('wheel velocity');

% legend('ipsi = blue');
            ax = gca; 
%             ax.FontSize = 24;

% source is the continuous data.
% indices are movement stamps
% window is [1000 2000]