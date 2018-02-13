%% This script is to be run to call TNC_MoverBehaviorExtractBS_v6 & TNC_ConvertTSDtoPopDataBS
% EAS 11.23.16

%% Major update 11.22.16 - to be used for behavior/spike alignment calling tb structure,created by behaveLoad6working.m

% USE THIS FOR ERROR TIME-OUT SESSIONS! USE PREVIOUS VERSION IF NO ERROR TIME-OUTS
% Inputs:
% filenamestr: path to the file, for example, '/Users/dudmanj/Desktop/m42/m42-20110427-wide-plustim-008.ns5'
% targetName: will save the data as this, can be whatever you want 
% dataRate: data type, needs to be 'ns4', 'ns3' or 'ns2' depending on acquisition file
% chan: for each channel, there will be a number associated with it. Need to find chan.rew, chan.thr, chan.lick, chan.x, chan.y. 
% determine these #s in: HeadFixedBehaviorLoadingScript- specifically, Load the continuous behavior monitor channels section
%      w/in dataAinp data file, look under dataAinp.MetaTags.ChannelID. Those #s are assoc. with specific blackrock channels.
% must previously run behaveLoad6working.m for extracting the relevant behavior from csv files; this is a _tb.mat file with a "taskbase" structure for further alignment to neural data 

% FILES NEEDED: _tb.mat file (TS of behavioral data), TNC_MoveBehaveExtractContinuous, TNC_ConvertTSDtoPopDataBS (TSs of neural data)
% UPDATED 11.26.15 for Vgattwo and Vgatthree recording/tagging/behavior experiments
% ainp1: trial start & solenoid input
% ainp2: laser - chan.stim = ContData.behavior.threshInds (after running TNC_MoverBehaviorExtractBS.m) = TSs for laser pulses.
% ainp3: motor input

% NOTE: 
% IMPORTANT FOR FIRST TIME USE: When calling openNEV.m - "BasicHeader" var can only be read when file read & write priveledges are turned on for the folder in which the data (.nev) are stored

%% Define vars for behavioral input (as the continuous signal to Blackrock
filenamestr = '~';                    
dataRate = 'ns4';    %10kS/sec
% dataRate = 'ns5';  %30kS/sec
% dataRate = 'ns2';  %1kS/sec: likely not fast enough to keep up w/ solenoid valve discharge < 1ms
% dataRate = 'ns6';  %raw

%% UPDATED 11.26.15 for trial start/solenoid input as new ain1; motor was for ain3. Relevant to Vgattwo and Vgatthree recordings
chan.rew = 129;                         % solenoid input & trial starts:   Ainp = 1
chan.thr = 130;                         % for laser - used for tagging:    Ainp = 2
chan.y = 131;                           % now for lick (post 5.1.16):      Ainp = 3
chan.x = [];                

%% Continuous Data structure is returned with ContData.behavior.threshInds = TS for behavioral events (i.e., laser pulses) 

% targetName = '161015New';
% targetName = 'test3';
targetName = '160505New';

elimTime = 270;                                                             % (s) [for 90 s chunks, and elimination of first # chunks, if needed]

%% for loading the proper nsx file:

initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/160505';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/test';

%% Load the extracted behavior file structure:

% taskbaseFile = 'mVgatfourtag_2016_05_05_16_tb.mat';
taskbaseFile = 'mVgatsix_tb.mat';
% taskbaseFile = 'test3';

taskbase = strcat(initPath,'/', taskbaseFile);
filestr = taskbase;
   
%% Load the proper tsd file: 

fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/160505';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/test';

[source] = strcat(fpath, '*_tsd.mat');
sessNum = 1;

[ContData, filenamestrB] = TNC_MoveBehaveExtractContinuous1(taskbase, filenamestr, targetName, dataRate, chan, elimTime, fpath);  %This will save the extracted recording data of this file to the day of interest (ex: 150529ipsi) as a bh.mat file

%% Get TSs for sorted (.ns5) neuronal spikes

% [PopData] = TNC_ConvertTSDtoPopDataBS(source, sessNum, fpath);
% 
% % Save the PopData structure so there's no need to re-extract sorted units:
%  cd(fpath);
%  save([targetName '_pop.mat'],'PopData');                              
%  disp(['saved as ' targetName '_pop.mat']);
    
%% OLD: Now plot movement of wheel from ContData.behavior structure: apply additional filters for plotting - no longer needed 1.26.16

trialStartsNsol = single(ContData.behavior.trialStrtsNSolData);
trialStartsNsol_sgolay = single(ContData.behavior.trialStrtsNSolData_sgolayF);
% 
% filtMin = 10;
% trialStartsNsol = single(ContData.behavior.trialStrtsNSolData);
% for w = 1:length(trialStartsNsol)/filtMin
%     newW(w) = filtMin*w -9;
%     trialStartsNsolFiltMin = trialStartsNsol(newW);
% end

% % % apply 1 kHz filter (needed for more instantaneous fluctuations) - but keep the first voltage; Note: index numbers will be shifted to the right by 1000 
% % filtLess = 1000;
% % trialStartsNsol = single(ContData.behavior.trialStrtsNSolData);
% % for w = 1:length(trialStartsNsol)/filtLess
% %     newW(w) = filtLess*w - 999;
% %     trialStartsNsolFiltLess = trialStartsNsol(newW);
% % end
% 
figure; hold on;
plot(trialStartsNsol_sgolay,'b','LineWidth',1,'Color',[0 0 0]); title(['160505, sgoLayFilt']);
    ylabel('/10 = mV');  %ylim([-1200 500]);                                  
    xlabel('time (x = x0,000 ms)'); %xlim([-1200 500]);
hold off;
% % 
% figure; hold on;
% plot(trialStartsNsol,'b','LineWidth',1,'Color',[0.5 0.5 0.5]); title(['161005, no filt']);
% hold off;
trialStartsNsol = trialStartsNsol(1:200000);
figure; hold on;
plot(trialStartsNsol,'b','LineWidth',1,'Color',[0.5 0.5 0.5]); title(['160505, no filt']);
% hold off;

