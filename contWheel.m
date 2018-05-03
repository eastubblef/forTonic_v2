%% This script is for plotting the average wheel trajectory aligned to event times (i.e. movement onset)

fpath = '/Volumes/My Passport for Mac/171013/behaveNewUnits';
% cd fpath;
filestr = '171013002.ns4';
behavFile = strcat(fpath, '/', '171013_bh3Rxn.mat');
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

eventWin = [1e3 2e3];

chan.wheel = 131;                                   
wheelMvs1 = find(data.MetaTags.ChannelID==chan.wheel);
wheelMvs2 = double(data.Data(wheelMvs1,:)); 
wheelMvs = decimate(wheelMvs2,10); 
% figure; hold on; plot(wheelMvs, 'k'); xlim([1e5 1.5e5]);    %longer trace snippet filtered

wheel_golay = sgolayfilt(wheelMvs, 9, 101);
figure; hold on; plot(wheel_golay, 'k'); xlim([1e5 1.5e5]); %longer trace snippet filtered

% wheel_golay = sgolayfilt(wheelMvs, 9, 101);
sourceIpsi = eventTimesIpsi;
sourceContra = eventTimesContra;
% source = wheel_golay;

window = eventWin;
% indicesIpsi = wheelALL_LfirstValid;
% indicesContra = wheelALL_RfirstValid;
indices= wheelALL_LfirstValid;
source = wheelMvs;
[sink] = TNC_ExtTrigWins(source,indices,window);  %plot the ipsi wheel movements
sinkIpsi = sink;

clear indices;
indices = wheelALL_RfirstValid;
[sink] = TNC_ExtTrigWins(source,indices,window);
sinkContra = sink;
% wheel_postSgolay = sgolayfilt(source, 9, 101);
figure; hold on;
plot(sinkIpsi.avg, ); title('171013')  %xlim([-100 2000]);
plot(sinkContra.avg); %xlim([-100 2000]);
xlabel('ms', 'FontSize', 24); ylabel('wheel velocity', 'FontSize', 24);
legend('show');
            ax = gca; 
            ax.FontSize = 24;

% source is the continuous data.
% indices are movement stamps
% window is [1000 2000]