%% This script is to be run to call Josh's TNC_MoverBehaviorExtract & TNC_ConvertTSDtoPopData

% filenamestr: path to the file, for example, '/Users/dudmanj/Desktop/m42/m42-20110427-wide-plustim-008.ns5'
% targetName: will automatically save the data as this, can be whatever you want 
% dataRate: data type, needs to be 'ns4', 'ns3' or 'ns2' depending on file
% chan: for each channel, there will be a number associated with it. They'll be different depending on your setup. Need to find chan.rew, chan.thr, chan.lick, chan.x, chan.y. 
% should be able to figure numbers out by running HeadFixedBehaviorLoadingScript- specifically the Load the continuous behavior monitor channels section
% within the dataAinp data file, look under dataAinp.MetaTags.ChannelID. Those numbers are the different numbers associated with specific blackrock channels.
% 
% Summary:
% Files from TONIC_v2 needed: TNC_MoverBehaviorExtract, HeadFixedBehaviorLoadingScript, openNSx, openNEV
% Use HeadFixedBehaviorLoadingScript to get information for the chan variable in TNC_MoverBehaviorExtract
% Create chan variable to be used in TNC_MoverBehaviorExtract- you'll need chan.rew, chan.thr, chan.lick, chan.x, chan.y (those are the variables currently in TNC_MoverBehaviorExtract)
%

% Notes:
% Beth's Blackrock settings (analog inputs):
% ainp1: lick input with sampling rate of 10 kS/s, HP filtered at 100 Hz with a threshold of -14ish (chan ID = 129)
% ainp2: laser input filtered with sampling rate of 10 kS/s, thresholded to 18 mV (chan ID = 130)
% ainp3: motor input ("extract spikes" is unchecked as of 6/23/15 - (chanID = 131)); up is R movement, down is L

% UPDATED 11.26.15 for Vgattwo and Vgatthree recording/tagging/behavior experiments
% ainp1: trial start & solenoid valve input
% ainp2: laser
% ainp3: motor input

% Also, need chan.stim = ContData.behavior.threshInds (after running TNC_MoverBehaviorExtractBS.m)
% ContData.behavor.threshInds will be the TS (ms) for behavioral events (i.e., laser pulse in).

% Calls TNC_ConvertTSDtoPopDataBS for TS of NEURAL data
% PopData will be the structure of TS for actual sorted (raw, .ns5) spikes by unit # per shank per session. Source calls the tsd.mat file.
% NO LONGER calls raster.m for final plotting of spikes' TSs to aligned behavioral event... use another mfile (oneUnitAlignSpike2BS)for that.

% File naming: 151105: 1001 = acquired data; 1001_bh.mat = laser-tagged analyzed data; 1002_bh.mat = ns4 channel analyzed data

% NOTE: 
% IMPORTANT FOR FIRST TIME USE: When calling openNEV.m - "BasicHeader" var can only be read when file read & write priveledges are turned on for the folder in which the data (.nev) are stored

%% New issue: 1.14.16
% 1001 in fname indicates tagging only (i.e. laser in only); 
% 1002 or 2002 fnames will be behaviorally-relavent 
%   -151105_90chunk1_3_1002.ns5_shank1_tsd.mat, saved from TNC_SS_GUI 
%   -151105_shank1_seg3_tsd.mat, saved from ExtendManualSortBS2


%% Begin - must define vars for behavioral input to Blackrock
% Inputs
filenamestr = '~';                     % BS prefers uiget - called in TNC_Mover openNSxBS & openNevBS (ease of going around others' naming strategies)
dataRate = 'ns5';

%% UPDATED 11.26.15 for trial start/solenoid input as new ain1; motor in for ain3.
chan.thr = 130;                         % for laser - used for tagging; Josh says corresponds to Ainp = 2
chan.y = 131;                           % for motor - also not used for tagging; Josh says corresponds to Ainp = 3
chan.x = [];
% chan.lick = 129                       % lick input - piezo not used currently
chan.rew = 129;                         % solenoid input and trial starts from Ainp = 1

%% Continuous Data structure is returned with ContData.behavior.threshInds = TS for behavioral events (i.e., laser pulses) 

% targetName = '151105_all';
% targetName = '161015_all';
% targetName = '151105_newUnits';
% targetName = '151105_newUnits2BehaveRxn';
% targetName = '151119_newUnits2Behave';
% targetName = '151118_newUnits2Behave';
% targetName = '151106_newUnits2BehaveRxn';
% targetName = '151104_newUnits2BehaveRxn';

targetName = '151119_tag3';


if numel(targetName == 11)
    filename = targetName(1:end - 5);
else filename = targetName(1:end - 6);                                      % Otherwise the file will have 12 elements at the end
end

elimTime = 270;                                                             % (s) [for 90 s chunks, and elimination of first # chunks, if needed]

%% Load the proper nsx file:

% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/wfManSortTest/u1_5_9'
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/behaveChunks';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/161015';
% initPath = '/Volumes/My Passport for Mac/Vgatthree_updated/151105/2ndpass/better/behaveChunks';
% initPath ='/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/newBehaveUnits/behaveChunks';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/behave/untag/newBehaveUnits/untagged';
initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/behave/tag';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks/newBehaveUnits/behaveChunks';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151106/behave/newUnits';
% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151104';

taskbaseFile = '*tb.mat';
taskbase = strcat(initPath,'/', taskbaseFile);
filestr = taskbase;
   
%% Load the proper tsd file: (ex: 151105_shank1_seg3_tsd.mat from ExtendManualSortBS2)
% fpath = filestr;
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/behaveChunks'
% fpath = '/Volumes/My Passport for Mac/Vgatthree_updated/151105/2ndpass/better/behaveChunks';
fpath = initPath;
[source] = strcat(fpath, '*_tsd.mat');
sessNum = 1;

[ContData, filenamestrB] = TNC_MoveBehaveExtractNew_v1(taskbase, filenamestr, targetName, dataRate, chan, elimTime, fpath);  %This will save the extracted recording data of this file to the day of interest (ex: 150529ipsi) as a bh.mat file

%% Get TSs for sorted (.ns5) neuronal spikes

[PopData] = TNC_ConvertTSDtoPopDataBS(source, sessNum, fpath);

% Save the PopData structure so there's no need to re-extract sorted units:
 cd(fpath);
 save([targetName '_pop.mat'],'PopData');                              
 disp(['saved as ' targetName '_pop.mat']);
    
%% Now plot movement of wheel from ContData.behavior structure: apply additional filters for plotting

% apply 10 kHz filter - overkill
% filtMore = 10000;
% wheelDataBig = single(ContData.behavior.sLeverData(2,:));
% for w = 1:length(wheelDataBig)/filtMore
%     newW(w) = filtMore*w;
%     wheelDataFiltMore = wheelDataBig(newW);
% end

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
