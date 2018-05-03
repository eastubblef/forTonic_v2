%% This script is to be run to call TNC_MoverBehaviorExtractBS_v6 & TNC_ConvertTSDtoPopDataBS
% updated 6.2017
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
chan.thr = 130;                         % for laser or lick:               Ainp = 2; laser and lick 3.20.17
chan.y = 131;                           % for wheel                         % now for lick (post 5.1.16):      Ainp = 3; nope wheel plugged in again (1.11.17)
chan.x = [];                

%% Continuous Data structure is returned with ContData.behavior.threshInds = TS for behavioral events (i.e., laser pulses) 

% targetName = '170112newBehave';
% targetName = '170321behave';
% targetName = '170508Behave';
% targetName = '170501_testBehave';
% targetName = '170627behaveSC3';
% targetName = '170718behaveSC2';
% targetName = '170118behave';
% targetName = '170111behave';
targetName = '171012';
% targetName = '170428';
% targetName = '170111newBehave';
% targetName = '170118_newBehave';

elimTime = 270;                                                             % (s) [for 90 s chunks, and elimination of first # chunks, if needed]

%% for loading the proper nsx file:

% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170112';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170111';
% initPath  = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170321';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/test/testLaser';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170426/SC2';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170427';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/test/testVgat';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/test/170508/SC1';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/test/170508';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/safeKeepingLaserLast';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170527';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/SC2/170602';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/SC3';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/SC2';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170118';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/160505';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/171012';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/forWhisper';
% initPath = '/Volumes/My Passport for Mac/170111/newBehaveUnits';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170112/shanks3&4/4thPass/behaveSegs/behaveNewUnits';
initPath = '/Volumes/My Passport for Mac/171012/newBehaveUnits';
% initPath = '/Volumes/My Passport for Mac/Vgatfive/170118/behaveNewUnits';
%% Load the extracted behavior file structure:

% taskbaseFile = 'mVgatseven_tb.mat';
% taskbaseFile = 'mVgatfive_tb.mat';
% taskbaseFile = '170501test';
% taskbaseFile = '170508test';
% taskbaseFile = '170428';

% taskbaseFile = '170111';
taskbaseFile = '171012';
% taskbaseFile = '151105';
% taskbaseFile = '170602';
% taskbaseFile = '170627';
% taskbaseFile = '170718';
% taskbaseFile = '170118';
% taskbaseFile = '170112';

taskbase = strcat(initPath,'/', taskbaseFile);
filestr = taskbase;
   
%% Load the proper tsd file: 

fpath = initPath; %easiest to have it in same folder
[source] = strcat(fpath, '*_tsd.mat');
sessNum = 1;

[ContData, filenamestrB] = TNC_MoveBehaveExtractMid(taskbase, filenamestr, targetName, dataRate, chan, elimTime, fpath);  %This will save the extracted recording data of this file to the day of interest (ex: 150529ipsi) as a bh.mat file
% [ContData, filenamestrT] = TNC_MoverBehaviorExtractBSwheel_2(taskbase, filenamestr, targetName, dataRate, chan, elimTime, fpath);  %This will save the extracted recording data of this file to the day of interest (ex: 150529ipsi) as a bh.mat file

%% Get TSs for sorted (.ns5) neuronal spikes

[PopData] = TNC_ConvertTSDtoPopDataBS(source, sessNum, fpath);

% Save the PopData structure so there's no need to re-extract sorted units:
%  cd(fpath);
 save([targetName '_pop.mat'],'PopData');                              
 disp(['saved as ' targetName '_pop.mat']);
    
%% OLD: Now plot movement of wheel from ContData.behavior structure: apply additional filters for plotting - no longer needed 1.26.16

% trialStartsNsol = single(ContData.behavior.trialStrtsNSolData(80000:200000));
% trialStartsNsol_sgolay = single(ContData.behavior.trialStrtsNSolData_sgolayF);
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
% figure; hold on;
% plot(trialStartsNsol_sgolay,'b','LineWidth',1,'Color',[0 0 0]); title(['161015, sgoLayFilt']);
%     ylabel('/10 = mV');  %ylim([-1200 500]);                                  
%     xlabel('time (x = x0,000 ms)'); %xlim([-1200 500]);
% hold off;
% 
% figure; hold on;
% plot(trialStartsNsol,'b','LineWidth',1,'Color',[0.5 0.5 0.5]); title(['161015, no filt']);
% hold off;

