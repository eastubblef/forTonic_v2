%% This script is to be run to call TNC_MoverBehaviorExtractBS_v6 & TNC_ConvertTSDtoPopDataBS

%% Major update 11.03.16 - to be used for behavior/spike alignment calling taskbase structure that is created by behaveLoad6.m

% USE THIS FOR ERROR TIME-OUT SESSIONS! USE PREVIOUS VERSION IF NO ERROR TIME-OUTS
% filenamestr: path to the file, for example, '/Users/dudmanj/Desktop/m42/m42-20110427-wide-plustim-008.ns5'
% targetName: will automatically save the data as this, can be whatever you want 
% dataRate: data type, needs to be 'ns4', 'ns3' or 'ns2' depending on acquisition file
% chan: for each channel, there will be a number associated with it. They'll be different depending on your setup. Need to find chan.rew, chan.thr, chan.lick, chan.x, chan.y. 
% should be able to figure numbers out by running HeadFixedBehaviorLoadingScript- specifically the Load the continuous behavior monitor channels section
%      within the dataAinp data file, look under dataAinp.MetaTags.ChannelID. Those numbers are the different numbers associated with specific blackrock channels.
% must previously run behaveLoad5.m for extracting the relevant behavior from csv files; this is a _tb.mat file with a "taskbase" structure for alignment to neural data 

% Summary:
% Files needed: _tb.mat file (TS of behavioral data), TNC_MoverBehaviorExtractBS_v6, openNSxBS, openNEVBS, TNC_ConvertTSDtoPopDataBS (TSs of neural data)
% PRIOR TO 11.1.15: 
% Beth's Blackrock settings (analog inputs):
% ainp1: lick input with sampling rate of 10 kS/s, HP filtered at 100 Hz with a threshold of -14ish (chan ID = 129)
% ainp2: laser input filtered with sampling rate of 10 kS/s, thresholded to 18 mV (chan ID = 130)
% ainp3: motor input ("extract spikes" is unchecked as of 6/23/15 - (chanID = 131)); up is R movement, down is L

% UPDATED 11.26.15 for Vgattwo and Vgatthree recording/tagging/behavior experiments
% ainp1: trial start & solenoid valve input
% ainp2: laser - will become chan.stim = ContData.behavior.threshInds (after running TNC_MoverBehaviorExtractBS.m) = TSs for laser pulses.
% ainp3: motor input

% NOTE: 
% IMPORTANT FOR FIRST TIME USE: When calling openNEV.m - "BasicHeader" var can only be read when file read & write priveledges are turned on for the folder in which the data (.nev) are stored

%% Define vars for behavioral input (as the continuous signal to Blackrock
filenamestr = '~';                    
dataRate = 'ns4';    %10kS/sec
% dataRate = 'ns5';  %30kS/sec
% dataRate = 'ns2';  %1kS/sec
% dataRate = 'ns6';  %raw

%% UPDATED 11.26.15 for trial start/solenoid input as new ain1; motor in for ain3. Relevant to Vgattwo and Vgatthree recordings
chan.thr = 130;                         % for laser - used for tagging; Ainp = 2
chan.y = 131;                           % now for lick (post 5.1.16); was for motor; Ainp = 3
chan.x = [];                
chan.rew = 129;                         % solenoid input and trial starts from Ainp = 1

%% Continuous Data structure is returned with ContData.behavior.threshInds = TS for behavioral events (i.e., laser pulses) 

% targetName = '160505';
% targetName = '151119_untag';
% targetName = '160505_all';
% targetName = '161119test3';
% targetName = '161119behaveNew';
% targetName = '161118test';
% targetName = '51105UltraNew'; funky!
targetName = '161015';

elimTime = 270;                                                             % (s) [for 90 s chunks, and elimination of first # chunks, if needed]

%% for loading the proper nsx file:

initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/161015';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/test';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/161005';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better';

%% Load the extracted behavior file structure:
% taskbaseFile = 'mVgattwotag_2015_11_18_172457_tb.mat';
% taskbaseFile = 'mVgatfourtag_2016_05_05_16_tb.mat';

taskbaseFile = 'mVgatsix_tb.mat';

taskbase = strcat(initPath,'/', taskbaseFile);
filestr = taskbase;
   
%% Load the proper event files & tsd file: 

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/151119good/behaveChunks/tagged';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/160505/behaveChunks';
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/161015';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/test';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/161005';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better';

[source] = strcat(fpath, '*_tsd.mat');
sessNum = 1;

[ContData, filenamestrB] = TNC_MoveBehaveExtractNew_v6(taskbase, filenamestr, targetName, dataRate, chan, elimTime, fpath);  %This will save the extracted recording data of this file to the day of interest (ex: 150529ipsi) as a bh.mat file

%% Get TSs for sorted (.ns5) neuronal spikes

[PopData] = TNC_ConvertTSDtoPopDataBS(source, sessNum, fpath);

% Save the PopData structure so there's no need to re-extract sorted units:
 cd(fpath);
 save([targetName '_pop.mat'],'PopData');                              
 disp(['saved as ' targetName '_pop.mat']);
    
%% OLD: Now plot movement of wheel from ContData.behavior structure: apply additional filters for plotting - no longer needed 1.26.16

% % apply 10 kHz filter - overkill
% % filtMore = 10000;
% % wheelDataBig = single(ContData.behavior.sLeverData(2,:));
% % for w = 1:length(wheelDataBig)/filtMore
% %     newW(w) = filtMore*w;
% %     wheelDataFiltMore = wheelDataBig(newW);
% % end
% 
% % apply 1 kHz filter (needed for more instantaneous fluctuations) - but keep the first voltage; Note: index numbers will be shifted to the right by 1000 
% filtLess = 1000
% wheelDataBig = single(ContData.behavior.sLeverData(2,:));
% for w2 = 1:length(wheelDataBig)/filtLess
%     newW2(w2) = filtLess*w2 - 999;
%     wheelDataFiltLess = wheelDataBig(newW2);
% end
% 
% % apply .1 kHz filter (needed for more instantaneous fluctuations) - but keep the first voltage, index 1, 101, 201... (shifted to the right by 100)
% filtWayLess = 100;
% wheelDataBig = single(ContData.behavior.sLeverData(2,:));
% for w3 = 1:length(wheelDataBig)/filtWayLess
%     newW3(w3) = filtWayLess*w3 - 99;
%     wheelDataFiltWayLess = wheelDataBig(newW3);
% end
% 
% figure; hold on;
% plot(wheelDataFiltLess,'b','LineWidth',1,'Color',[0.5 0.5 0.5]); title(['1 kHz additional; wheel movements, R = up ']);
% hold off;
% 
% figure; hold on;
% plot(wheelDataFiltWayLess,'b','LineWidth',1,'Color',[0.5 0.5 0.5]); title(['.1 kHz additional; wheel movements, R = up ']);
% hold off;
% 
% %Find "center" of raw wheel voltages
% centers = median(wheelDataBig); 
% wheelMin = min(wheelDataBig);    %just to see
% wheelMax = max(wheelDataBig);    %just to see
% 
% %Movements >25ish seem "real" for 151105 data:
% % tHold = 25;
% % for l = 1:length(wheelDataBig(end-1))
% %     moveL = find(wheelDataBig(l+1) < wheelDataBig(l));
% % end
% 
% %% pseudocoding: keep it simple:
% % zeel = [0 2 3 4 3 2 0]
% % 
% % neel = [1 10 20 30 40 50 60]
% % newNeel = [];
% % for z = 1:length(zeel)
% %     if zeel(z) > zeel(z+1)
% %         if z+1 == length(zeel)
% %             break;
% %         else
% %             
% %             newNeel(z) = neel(z);  %zero vals will hold index positions in which case is not true
% %         
% %         end
% %     end
% % end
% 
