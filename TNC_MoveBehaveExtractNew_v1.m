
 function [ContData, filenamestrB] = TNC_MoveBehaveExtractNew_v1(taskbase, filenamestr, targetName, dataRate, chan, elimTime, fpath)
%% This is EAS's version of Josh's from Tonic_v2 

% UPDATED LAST: 
% 2.2.18  - incorporate slow rxn time extractions
% 3.21.16 - 3.29.16 use taskbase file to align TS of spikes to trial start & L v R movements (behaveLoad3)
% 3.2.16  - incorporated behavLoad.m for loading behavioral data from Beth's (updated) "taskbase" structure, based on csv file

%  ContData.behavior.threshInds = TS for laser pulse-on (divide by 30, 30kHz = sampling frequency for spikes)
% l. 51-76: chan.rew (129), now ONLY reflects every other TS (NOT proper sol reward inds) - same with trial starts

% RE-ESTABLISH ACTIVE INPUTS FOR 11.26.15 DATA (1.14.16 & 1.28.16): 
% motor in = chan.y (131) captured on ns4 file, sampling rate was 10kHz; divide by 10)

% Also saves a bh.mat file that could be loaded into the next mfile for TS from recordings: TNC_ConvertTSDtoPopDataBS
% Currently, this mfile is called from MoverScript.m

%% FILE NAMES
    
    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Beginning file: ' filenamestr(1,1:length(filenamestr)-3)]);
    disp(' ');disp(' ');  
    dispOn = 0;
    
%% Load event data: BS updated 2.27.16
    newPath = fpath;
    cd(newPath);
    if exist('filenamestr', 'var');
        [filenamestrE, path] = uigetfile('*.nev*','select the .nev file', fpath)                
    end    
    NEV = openNEVBS2(filenamestrE,'read','nosave','nomat', 'report');       %Fucking Kian; see this m file for stupid updates.        
    dataEvents = NEV;
    
    startsRewIndsTmp = find(dataEvents.Data.Spikes.Electrode == chan.rew);                    %inds of trial start and solenoids detected (ch. 129); get rid of first val, since it's a pulse detected when Orienter is started
    startsRewTsTmp = dataEvents.Data.Spikes.TimeStamp(startsRewIndsTmp(2:end))./30;                 % don't ./30 yet. will need to parse out trial start from solenoid here
    startsRewTsRoundTmp = round(dataEvents.Data.Spikes.TimeStamp(startsRewIndsTmp(2:end))./30);     % This is the value to compare - don't ./30 yet.     
    validStrtsRews1 = find(diff(startsRewTsRoundTmp)>50);                                     % eliminate pulses detected twice
    
    validStrtsRews = startsRewTsRoundTmp(validStrtsRews1);                                    % true nev clock from true 1st trial start
    startsRewZeroNev = validStrtsRews-validStrtsRews(1);                                      % sanity check: zero out the clock for comparison to taskbase
    
    % Tangent to investigate the diffs bt. ea. line (clock alignment troubleshooting):
    diffTmp = diff(startsRewZeroNev);
    itiNev = 3000;    
    solBufferNev = 290;                                                                       % seems to have an upper limit of 3289
    solMinusStarts = diffTmp < itiNev+solBufferNev & diffTmp > 3213;                          % give time buffer - usually in the positive direction for solenoid clicks to next trial start
    numSol2starts = numel(find(solMinusStarts == 1));                                         % compare to csv file
        
%% Load taskbase structure for behavioral data: BS updated 3.2.16 
%                      
    if exist('filenamestr', 'var');
        [filenamestrT, path] = uigetfile('*_tb.mat*','select the taskbase file', fpath)       % get the actual # correct trials (and TS) from csv file
    end
    taskbase = strcat(path, filenamestrT);
    load(taskbase);
 
    startTimesCsv = taskbase.trialStartTimes';
    startsCentsCsv = taskbase.events';
    startsRewsCsv = taskbase.startsRews';                                                            % trial starts and solenoid discharges from behavior file
    rewTimesCsv = taskbase.rewTimes';
    % Zero out csv trial starts for comparison/alignment with nev file
     validStrtsRewsCsv1 = find(diff(startsRewsCsv)>50);                                               % account for very low vals: rx time that is virtually impossible due to mouse likely already moving
     validStrtsRewsCsv = startsRewsCsv(validStrtsRewsCsv1);                                           % align this clock to nev
     startsRewsZeroCsv = validStrtsRewsCsv-validStrtsRewsCsv(1);                                      % zero out TS 1; first trial start = 0
    
    % Tangent to investigate the diffs bt. ea. line: (clock alignment troubleshooting)
    diffStartsRewsCsv = diff(startsRewsZeroCsv);                                                     % refine the 3 s rule for the tighter csv file:
    meanDiffCsv = mean(diffStartsRewsCsv);
    solBufferCsv = 80;
    iti = 3000;
    solMinusStartsCSV = diffStartsRewsCsv < iti+solBufferCsv & diffStartsRewsCsv > iti-solBufferCsv; % give time buffer 
    numSol2startsCsv = numel(find(solMinusStartsCSV == 1));                                          % should be the num of csv rewarded trials; compare to nev

%% 3.3.16 - Parse out solenoid/trial start info 
    
    sol_diffs = diffTmp(solMinusStarts == 1);                                   % 2nd trial start is first; these are the nev diffs bt sol and next trial start (for uninterrupted trials);
    meanSol_diffs = mean(sol_diffs); 
    newRewTS = [];
    for i = 2:length(validStrtsRews)
         if solMinusStarts(i) == 1
            newRewTS(i) = validStrtsRews(i);                                    % vector (with inds) for when TS reflect a solenoid click
         end                                                                    % zeros are either trial starts or abberrant pulse into nev
         if i == length(solMinusStarts)
             break
         end
    end
    solsMinusTS1 = find(newRewTS > 10);
    solsMinusTS = newRewTS(solsMinusTS1);                                       % these are trial starts, only after rew trials, starting from second trial
    ContData.behavior.solsMinusTS = solsMinusTS;                                % previously plotted these on rasters

%% Tangent: compare the mean(diff(csv = validStrtsRewsCsv :nev = validStrtsRews)) TS.

    newRewTScsv = [];
    for i = 2:length(validStrtsRewsCsv)
         if solMinusStartsCSV(i) == 1
            newRewTScsv(i) = validStrtsRewsCsv(i);                          % vector (with inds) for when TS reflect a solenoid click
         end                                                                % zeros are either trial starts or abberrant pulse into nev
         if i == length(solMinusStartsCSV)
             break
         end
    end
    solsMinusTS1csv = find(newRewTScsv > 10);
    solsMinusTScsv = newRewTScsv(solsMinusTS1csv);                          % these are trial starts, only after rew trials, starting from second trial
    clockRunNev = solsMinusTS(end) - solsMinusTS(1) ;                       % should be ~ 3 s (3000 ms) diff due to offset in favor of nev file?
    clockRunCsv = solsMinusTScsv(end) - solsMinusTScsv(1);
    clockRunDiff = clockRunNev - clockRunCsv;
    clockDiffs = mean(diff(validStrtsRewsCsv));

%% Align clocks to csv trial starts:

    % startTimesCsv = startTimesCsv - startTimesCsv(1);                         % includes correct and re-started trials
    LRtrial = taskbase.LRtrial;
    leftIndsCenteredTimes = taskbase.leftIndsCenteredTimes;                     % rewarded trials only
    rightIndsCenteredTimes = taskbase.rightIndsCenteredTimes;                   % rewarded trials only
    leftIncorrectTimes = taskbase.leftIncorrectTimes;                           % incorrect threshold was reached
    rightIncorrectTimes = taskbase.rightIncorrectTimes;                         % incorrect threshold was reached
    trialStrts = [];
    for i = 1:length(startsRewsZeroCsv)
        for j = 1:length(startsRewZeroNev)
            if i == 1
                trialStrts(i) = 1;                                               % report a 1 for the first trial alignment (and since the first trial start for both clocks (by design) is defined as aligned) 
            else
                if startsRewZeroNev(j) > startsRewsZeroCsv(i)-60 && startsRewZeroNev(j) < startsRewsZeroCsv(i)+60; % creates a window to search within 40 ms +/-
                    trialStrts(j) = startsRewZeroNev(j);                         % indexed positions are true to the zeroed (to first true trial start) nev clock
                end                                                             
            end                                                                  % then, after this, index to trial starts csv and find L mvments v R mvments per trial
        end
    end
    trialStartsPre = find(trialStrts >= 1);  %may need to use the indices for rxn time extractions 2.2.18
    trialStarts = trialStrts(trialStartsPre);

    lCenteredTimesCsv = taskbase.leftCenteredTimes;
    lCenteredZeroCsv = lCenteredTimesCsv - startTimesCsv(1);                     % zero-out l/r center times csv vector to index back into zeroed nev file
    
    lIncorrectZeroCsv = leftIncorrectTimes - startTimesCsv(1);
    
    % sanity check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    lIndsCenteredTimesZero = [];                                                 % zero-out l center times, keeping indexed position
    for x = 1:length(leftIndsCenteredTimes)
        if leftIndsCenteredTimes(x) > 1
            lIndsCenteredTimesZero(x) = leftIndsCenteredTimes(x) - startTimesCsv(1);
        end
        if leftIndsCenteredTimes(x) == 0
            lIndsCenteredTimesZero(x) = 0;                                       % this vector should be the same as lCenteredZeroCsv except at correct index position for trial number
        end
    end    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % L rewarded trials' start times
    lCenters = [];
    for a = 1:length(trialStarts)
        for b = 1:length(lCenteredZeroCsv)
            if a == numel(trialStarts)
                break
            end
            if b == numel(lCenteredZeroCsv)
                break
            end
            if lCenteredZeroCsv(b) > trialStarts(a) && lCenteredZeroCsv(b) < trialStarts(a+1) && lCenteredZeroCsv(b) < lCenteredZeroCsv(b+1) % just to be sure
               lCenters(a) = trialStarts(a);   %has trial inds
            end
        end
    end
    
    %modified 2.2.18:
    lCenters1 = lCenters > 0; LrewTrialStartsZeroed = lCenters(lCenters1);
    lCenters1Inds = 1:length(lCenters1); LrewTrialStartsZeroedInds = lCenters1Inds(lCenters1); %these are L trial inds
    
    validsCorrection = single(validStrtsRews(1));
    trialStartsValid = trialStarts + validsCorrection;
   
    trialInds = 1:length(trialStartsValid);             %to keep track of trialInds
    LrewTrialsStarts = LrewTrialStartsZeroed + validsCorrection;
    
    % L incorrect trials' start times
    lIncorrects = [];
    for c = 1:length(trialStarts)
        for d = 1:length(lIncorrectZeroCsv)
            if c == numel(trialStarts)
                break
            end
            if d == numel(lIncorrectZeroCsv)
                break
            end
            if lIncorrectZeroCsv(d) > trialStarts(c) && lIncorrectZeroCsv(d) < trialStarts(c+1) && lIncorrectZeroCsv(d) < lIncorrectZeroCsv(d+1)
               lIncorrects(c) = trialStarts(c);
            end
        end
    end
    lIncorrects1 = lIncorrects > 0; LincorrectsTrialStartsZeroed = lIncorrects(lIncorrects1);
    LincorrectTrialsStarts = LincorrectsTrialStartsZeroed + validsCorrection;
   
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % right side trials
    rCenteredTimesCsv = taskbase.rightCenteredTimes;
    rCenteredZeroCsv = rCenteredTimesCsv - startTimesCsv(1);
    rIncorrectZeroCsv = rightIncorrectTimes - startTimesCsv(1);

    % r rewarded trials' start times
    rCenters = [];
    for m = 1:length(trialStarts)
        for n = 1:length(rCenteredZeroCsv)
            if m == numel(trialStarts)
                break
            end
            if n == numel(rCenteredZeroCsv)
                break
            end
            if rCenteredZeroCsv(n) > trialStarts(m) && rCenteredZeroCsv(n) < trialStarts(m+1) && rCenteredZeroCsv(n) < rCenteredZeroCsv(n+1) % just to be sure
               rCenters(m) = trialStarts(m);
            end
        end
    end
    
    %modified 2.2.18:
    rCenters1 = rCenters > 0; RrewTrialStartsZeroed = rCenters(rCenters1);
    rCenters1Inds = 1:length(rCenters1); RrewTrialStartsZeroedInds = rCenters1Inds(rCenters1); %these are L trial inds
    RrewTrialsStarts = RrewTrialStartsZeroed + validsCorrection; 
    
    %r incorrect trials' start times
    rIncorrects = [];
    for o = 1:length(trialStarts)
        for p = 1:length(rIncorrectZeroCsv)
            if o == numel(trialStarts)
                break
            end
            if p == numel(rIncorrectZeroCsv)
                break
            end
            if rIncorrectZeroCsv(p) > trialStarts(o) && rIncorrectZeroCsv(p) < trialStarts(o+1) && rIncorrectZeroCsv(p) < rIncorrectZeroCsv(p+1)
                rIncorrects(o) = trialStarts(o);
            end
        end
    end
    rIncorrects1 = rIncorrects > 0; RIncorrectTrialStartsZeroed = rIncorrects(rIncorrects1);
    RincorrectTrialsStarts = RIncorrectTrialStartsZeroed + validsCorrection;

    ContData.behavior.LrewTrialStarts = LrewTrialsStarts;
    ContData.behavior.trialStartsValid = trialStartsValid;
    ContData.behavior.RrewTrialStarts = RrewTrialsStarts;
    ContData.behavior.LincorrectTrialStarts = LincorrectTrialsStarts;
    ContData.behavior.RincorrectTrialStarts = RincorrectTrialsStarts;
    
    rewardedTrialsNev = numel(ContData.behavior.LrewTrialStarts) + numel(ContData.behavior.RrewTrialStarts);
    disp(['Number of rewarded trials (.nev): ' num2str(rewardedTrialsNev)]);
    disp(' ');disp(' '); 
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Load the taskbase data relevent to wheel movement extractions L v R. Could compare to ns4 file. 
    % Note: these values still retain 0s to hold indices where L or Rtrials were either not rewarded, or were rewarded to other side.
    
    % Correct L and R trials:
    wheelLfirstTime = taskbase.wheelLfirstTime;
    wheelLfirst = taskbase.wheelLfirst;
    wheelRfirst = taskbase.wheelRfirst;
    wheelRfirstTime = taskbase.wheelRfirstTime;

    %Modified 2.2.18: keep track of the trial indices
    % Remove the empty indices and extract just the first vals needed and use LrewTrialStarts/RewTrialStarts for alignment
    wheelLfirstT = wheelLfirstTime > 0; wheelLfirsts = wheelLfirstTime(wheelLfirstT);
    wheelLTinds = 1:length(wheelLfirstT); wheelLfirstTinds = wheelLTinds(wheelLfirstT); %thse are the trial inds for Lmvmnts
    wheelLfirstsZero = wheelLfirsts - startTimesCsv(1);
    
    wheelRfirstT = wheelRfirstTime > 0; wheelRfirsts = wheelRfirstTime(wheelRfirstT);
    wheelRTinds = 1:length(wheelRfirstT); wheelRfirstTinds = wheelRTinds(wheelRfirstT); %these are the trial indices for R mvmnts
    wheelRfirstsZero = wheelRfirsts - startTimesCsv(1); 

    %Added 2.2.18:
    maxTrialInds = max(trialInds);
    maxWheelLfirstT = length(wheelLfirstT);
    if maxTrialInds > maxWheelLfirstT
        LTrialInds = NaN(maxTrialInds,1);
    else LTrialInds = NaN(maxWheelLfirstT,1);
    end
    LTrialIndsMax1 = 1:length(LTrialInds);
    LTrialInds = LTrialIndsMax1(wheelLfirstT == 1);  %use this to test against the rxMat
    
    maxWheelRfirstT = length(wheelRfirstT);
    if maxTrialInds > maxWheelRfirstT
        RTrialInds = NaN(maxTrialInds,1);
    else RTrialInds = NaN(maxWheelRfirstT,1);
    end
    RTrialIndsMax1 = 1:length(RTrialInds);
    RTrialInds = RTrialIndsMax1(wheelRfirstT == 1);  %use this to test against the rxMat
      
    extractWheelCsvLR = numel(wheelLfirsts) + numel(wheelRfirsts);
 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ pick up here
    % Incorrect L and R trials:
    wheelILfirstTime = taskbase.wheelILfirstTime;
    wheelILfirst = taskbase.wheelILfirst;
    wheelIRfirst = taskbase.wheelIRfirst;
    wheelIRfirstTime = taskbase.wheelIRfirstTime;

    rightIncorrectTimes = taskbase.rightIncorrectTimes
    wheelILfirstT = wheelILfirstTime > 0; wheelILfirsts = wheelILfirstTime(wheelILfirstT);
    wheelILfirstsZero = wheelILfirsts - startTimesCsv(1);
    wheelIRfirstT = wheelIRfirstTime > 0; wheelIRfirsts = wheelIRfirstTime(wheelIRfirstT);
    wheelIRfirstsZero = wheelIRfirsts - startTimesCsv(1);

    extractWheelCsvILnIR = numel(wheelILfirsts) + numel(wheelIRfirsts);
    
%% aligned (to nev-based trial starts) wheel movements that were based on csv file
%  May want to consider including this in the overall trial alignment for absolute certainty - rx times (diffL, diffR) are on super-fast side for some trials (25 ms)
    
    % Use this for most sessions (151105 incl.) correct L movements: eliminates first L movements that don't fit to an extracted trial start into Blackrock
    wheelLfirstsAligned = []; wheelLfirstsAlignedInds = [];
    for i = 1:length(LrewTrialStartsZeroed)
        if numel(LrewTrialStartsZeroed) < numel(wheelLfirstsZero) && i == numel(LrewTrialStartsZeroed)
            break
        else if i == numel(wheelLfirstsZero)
                break
            else if wheelLfirstsZero(i) > LrewTrialStartsZeroed(i) && wheelLfirstsZero(i) < LrewTrialStartsZeroed(i+1) && wheelLfirstsZero(i) < wheelLfirstsZero(i+1) % just to be sure
                    wheelLfirstsAlignedInds(i) = wheelLfirstTinds(i);      %added 2.2.18
                    wheelLfirstsAligned(i) = wheelLfirstsZero(i);
                else if wheelLfirstsZero(i) < LrewTrialStartsZeroed(i) && wheelLfirstsZero(i+1) > LrewTrialStartsZeroed(i) && wheelLfirstsZero(i+1) < LrewTrialStartsZeroed(i+1)
                        wheelLfirstsAlignedInds(i) = wheelLfirstTinds(i+1); %added 2.2.18
                        wheelLfirstsAligned(i) = wheelLfirstsZero(i+1);
                    end
                end
            end
        end
    end 
        
    wheelILfirstsAligned = [];
    for j = 1:length(LincorrectsTrialStartsZeroed)
        if numel(LincorrectsTrialStartsZeroed) < numel(wheelILfirstsZero)&& j == numel(LincorrectsTrialStartsZeroed)
            break
        else if j == numel(wheelILfirstsZero)
                break
            else if wheelILfirstsZero(j) > LincorrectsTrialStartsZeroed(j) && wheelLfirstsZero(j) < LincorrectsTrialStartsZeroed(j+1) && wheelILfirstsZero(j) < wheelILfirstsZero(j+1)
                    wheelILfirstsAligned(j) = wheelILfirstsZero(j);
                else if wheelILfirstsZero(j) < LincorrectsTrialStartsZeroed(j) && wheelILfirstsZero(j+1) > LincorrectsTrialStartsZeroed(j) && wheelILfirstsZero(j+1) < LincorrectsTrialStartsZeroed(j+1)
                        wheelILfirstsAligned(j) = wheelILfirstsZero(j+1);
                    end
                end
            end
        end
    end
      
%     % Use for 151106? Alternative for alignment of L movement vectors if above doesn't work       
%     wheelLfirstsAligned = []; wheelLfirstsAlignedInds = [];
%     for i = 1:length(LrewTrialStartsZeroed)
%         for j = 1:length(wheelLfirstsZero)
%             if i == numel(LrewTrialStartsZeroed)
%                 break
%             end
%             if j == numel(wheelLfirstsZero)
%                 break
% %             else if wheelLfirstsZero(j) < LrewTrialStartsZeroed(i) && wheelLfirstsZero(j) > LrewTrialStartsZeroed(i-1)
% %                     wheelLfirstsAligned(j) = wheelLfirstsZero(j-1);
% %                else if wheelLfirstsZero(j+1) < LrewTrialStartsZeroed(i+1) && wheelLfirstsZero(j+1) > LrewTrialStartsZeroed(i+1)
% %                     wheelLfirstsAligned(j) = wheelLfirstsZero(j-1);
%                  else if wheelLfirstsZero(j) > LrewTrialStartsZeroed(i) && wheelLfirstsZero(j) < LrewTrialStartsZeroed(i+1) && wheelLfirstsZero(j) < wheelLfirstsZero(j+1) % just to be sure
%                          wheelLfirstsAligned(j) = wheelLfirstsZero(j);
%                          wheelLfirstsAlignedInds(j) = wheelLfirstTinds(j);      %added 2.2.18
% %                     end
%                 end
%             end
%         end
%     end
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
    % R movements: works better this way for this file's R movements. Think about this.
    wheelRfirstsAligned = []; wheelRfirstsAlignedInds = [];
    for i = 1:length(RrewTrialStartsZeroed)
        for j = 1:length(wheelRfirstsZero)
            if i == numel(RrewTrialStartsZeroed)
                break
            end
            if j == numel(wheelRfirstsZero)
                break
            end
            if wheelRfirstsZero(j) > RrewTrialStartsZeroed(i) && wheelRfirstsZero(j) < RrewTrialStartsZeroed(i+1) && wheelRfirstsZero(j) < wheelRfirstsZero(j+1) % just to be sure
                wheelRfirstsAlignedInds(j) = wheelRfirstTinds(j);
                wheelRfirstsAligned(j) = wheelRfirstsZero(j);
            end
        end
    end
    
    wheelIRfirstsAligned = [];
    for k = 1:length(RIncorrectTrialStartsZeroed)
        for l = 1:length(wheelIRfirstsZero)
            if k == numel(RIncorrectTrialStartsZeroed)
                break
            end
            if l == numel(wheelIRfirstsZero)
                break
            end
            if wheelIRfirstsZero(l) > RIncorrectTrialStartsZeroed(k) && wheelIRfirstsZero(l) < RIncorrectTrialStartsZeroed(k+1) && wheelIRfirstsZero(l) < wheelIRfirstsZero(l+1)
                wheelIRfirstsAligned(l) = wheelIRfirstsZero(l);
            end
        end
    end
        
    % Alternative 2 for alignment of R movement vectors if above doesn't work
%         if wheelRfirstsZero(i) > RrewTrialStartsZeroed(i) && wheelRfirstsZero(i) < RrewTrialStartsZeroed(i+1) && wheelRfirstsZero(i) < wheelRfirstsZero(i+1) % just to be sure
%             wheelRfirstsAligned(i) = wheelRfirstsZero(i);
%             
%         else if wheelRfirstsZero(i) < RrewTrialStartsZeroed(i) && wheelRfirstsZero(i+1) > RrewTrialStartsZeroed(i) && wheelRfirstsZero(i+1) < RrewTrialStartsZeroed(i+1)
%                 wheelRfirstsAligned(i) = wheelRfirstsZero(i+1);
%             end
%         end
%     end
    
%     if isempty(wheelIRfirstsAligned)
%         for i = 1:length(RIncorrectTrialStartsZeroed)
%             if i == numel(RIncorrectTrialStartsZeroed)
%                 break
%             end
% 
%             if wheelIRfirstsZero(i) > RIncorrectTrialStartsZeroed(i) && wheelIRfirstsZero(i) < RIncorrectTrialStartsZeroed(i+1) && wheelIRfirstsZero(i) < wheelIRfirstsZero(i+1) % just to be sure
%                 wheelIRfirstsAligned(i) = wheelIRfirstsZero(i);
% 
%             else if wheelIRfirstsZero(i) < RIncorrectTrialStartsZeroed(i) && wheelIRfirstsZero(i+1) > RIncorrectTrialStartsZeroed(i) && wheelIRfirstsZero(i+1) < RIncorrectTrialStartsZeroed(i+1)
%                     wheelIRfirstsAligned(i) = wheelIRfirstsZero(i+1);
%                 end
%             end
%         end
%     end

    wheelsLfirstValid = wheelLfirstsAligned + validsCorrection;
    wheelsRfirstValid = wheelRfirstsAligned + validsCorrection;
    wheelsILfirstValid = wheelILfirstsAligned + validsCorrection;
    wheelsIRfirstValid = wheelIRfirstsAligned + validsCorrection;
    
%      export_mvmt_data_for_psths.wheelALL_LfirstValid = sort(horzcat(wheelIRfirstsAligned, wheelLfirstsAligned));
%      export_mvmt_data_for_psths.wheelALL_RfirstValid = sort(horzcat(wheelILfirstsAligned, wheelRfirstsAligned));
    wheelALL_LfirstValid = sort(horzcat(wheelsIRfirstValid, wheelsLfirstValid));
    wheelALL_RfirstValid = sort(horzcat(wheelsILfirstValid, wheelsRfirstValid));
    
    %sanity check - ensure that initial movements were after trial starts (i.e. diffL & diffR should be positive values)
    %Also could use these as rxn times for correct L and R-mv trials
    for j = 1:length(LrewTrialsStarts)
        if j == length(wheelsLfirstValid)
            break
        else
            diffL(j) = wheelsLfirstValid(j) - LrewTrialsStarts(j);  %2.2.18: could pull out these inds
        end
     end
     for j = 1:length(RrewTrialsStarts)
        if j == length(wheelsRfirstValid)
            break
        else
            diffR(j) = wheelsRfirstValid(j) - RrewTrialsStarts(j);  %2.2.18: could pull out these inds
        end
     end
     
     
     %New 2.2.18: These are the slow rxn time trials' mvmnt onset timestamps:
     wheelsLslow1 = find(diffL >= 400);
     wheelsLslow = wheelsLfirstValid(wheelsLslow1);
     
     wheelsRslow1 = find(diffR >= 400);
     wheelsRslow = wheelsRfirstValid(wheelsRslow1);
     
%      meanDiffL = mean(diffL); minDiffL = min(diffL); maxDiffL = max(diffL); medianDiffL = median(diffL);
%      meanDiffR = mean(diffR); minDiffR = min(diffR); maxDiffR = max(diffR); medianDiffR = median(diffR);
     
    %  wheelsRfirst = wheelRfirstsZero + validsCorrection; 
%      checkLmv = wheelsLfirstValid - single(trialStartsValid(1));                             % should be nearly the same as the zeroed L movements
%      checkRmv = wheelsRfirstValid - single(trialStartsValid(1));                             % should be nearly the same as the zeroed R movements
    %  ContData.behavior.wheelsLfirst = wheelsLfirst;
    %  ContData.behavior.wheelsRfirst = wheelsRfirst;
         
     ContData.behavior.wheelsLfirstValid = wheelsLfirstValid;
     ContData.behavior.wheelsRfirstValid = wheelsRfirstValid;
     ContData.behavior.wheelsILfirstValid = wheelsILfirstValid;
     ContData.behavior.wheelsIRfirstValid = wheelsIRfirstValid;
     
     ContData.behavior.wheelsLslowRxn = wheelsLslow;
     ContData.behavior.wheelsRslowRxn = wheelsRslow;
    
%% LOAD CONTINUOUS LEVER (WHEEL) DATA
    
 switch dataRate
    
    case 'ns4'    
%          filenamestrE = [filenamestr(1,1:length(filenamestr)-3) dataRate]      %BS change to accept my own files
        [filenamestrE, path] = uigetfile('*.ns4','select the .ns4 file', fpath); %BS works, but unneccessary extra opening of file by hand

        Ns4DATA = openNSxBS('report','read',filenamestrE);
        clear leverData sLeverData tmpLeverData
        
          xChan = find(Ns4DATA.MetaTags.ChannelID==chan.x);                    %my x channel is empty  
          yChan = find(Ns4DATA.MetaTags.ChannelID==chan.y);                        %Motor in
          lChan = find(Ns4DATA.MetaTags.ChannelID==chan.lick);
        
         leverData(1,:) = decimate(Ns4DATA.Data(xChan,:),10);
         leverData(2,:) = decimate(Ns4DATA.Data(yChan,:),10);                     %Correct for 10kS/s rate
         rawLick = decimate(Ns4DATA.Data(lChan,:),10);

    case 'ns3'
        Ns4DATA = openNSx('report','read',filenamestrE);
        clear leverData sLeverData tmpLeverData
        leverData(1,:) = decimate(Ns4DATA.Data(chan.x,:),2);
        leverData(2,:) = decimate(Ns4DATA.Data(chan.y,:),2);
        rawLick = decimate(Ns4DATA.Data(3,:),10);
        
        
%     case 'ns2'
% %         filenamestrE = [filenamestr(1,1:length(filenamestr)-3) dataRate]  %BS
%         filenamestrE = uigetfile('*.ns4','select the .ns4 file', dataRate);
%         Ns2DATA = openNSx('report','read',filenamestrE);
% 
%         clear leverData sLeverData tmpLeverData
%         xChan = find(Ns2DATA.MetaTags.ChannelID==chan.x);                   %BS
%         yChan = find(Ns2DATA.MetaTags.ChannelID==chan.y);
%         lChan = find(Ns2DATA.MetaTags.ChannelID==chan.lick);
%         leverData(1,:)  = Ns2DATA.Data(xChan,:);                          %BS
%         leverData(2,:)  = Ns2DATA.Data(yChan,:);
%         rawLick         = Ns2DATA.Data(lChan,:);
% 
% %         leverData(1,:)  = Ns4DATA.Data(xChan,:);                          %BS
% %         leverData(2,:)  = Ns4DATA.Data(yChan,:);
% %         rawLick         = Ns4DATA.Data(lChan,:);

end
    
%     sLeverData(1,:) = sgolayfilt(leverData(1,:),9,101);  %motor %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
%     sLeverData(2,:) = sgolayfilt(leverData(2,:),9,101);  %motor
% 
%     sWheelData = sLeverData;
%     ContData.behavior.sWheelData = sWheelData;
     
%     ContData.behavior.sLeverData = sLeverData;           %number of data points are = numel(sWheelData)/10 (and rounded up) for correction of 10kHz sampling rate

%     difLick = diff(rawLick);
%     evLick=zeros(1,numel(rawLick));
% 
%     ContData.behavior.evLick    = evLick;
%     ContData.behavior.rawLick   = rawLick;


%% 2.2.18: New for pulling out reaction times from the tb file for 151104-141119 data sets
csvRxnTimes = taskbase.rxCondensed;
trialStartsPre = trialStartsPre';
rxnMat = csvRxnTimes(trialStartsPre,:);

slow = find(rxnMat(:,2) >= 400);

slowRxnMat = rxnMat(slow,:);                %first col = trial number; 2nd col is rxn time > 400 ms; 3rd col is L/R correct/incorrect trials

findLSlows = slowRxnMat(:,3) == -1; 
LcorrSlowRxn = slowRxnMat(findLSlows,:);    %slow rxn times that were correct L

findRSlows = slowRxnMat(:,3) == 1;
RcorrSlowRxn = slowRxnMat(findRSlows,:);

%Determine if the known rxn-time trials correspond to L/R mv indices: 
% sameLTrial = any(LTrialInds == LcorrSlowRxn(:,1));
% LTrialInds = LTrialInds(sameLTrial);
% sameLTrialValid = any(wheelLfirstsAlignedInds == LcorrSlowRxn(:,1));
% LTrialSlows = wheelLfirstsAlignedInds(sameLTrialValid); %these are the slow L correct trials 
% ContData.behavior.wheelsLfirstSlowRxn = wheelsLfirstValid(LTrialSlows); %nope

% sameRTrial = any(RTrialInds == RcorrSlowRxn(:,1));
% RTrialInds = RTrialInds(sameRTrial);
% sameRTrialValid = any(wheelRfirstsAlignedInds == RcorrSlowRxn(:,1));
% RTrialSlows = wheelRfirstsAlignedInds(sameRTrialValid); %these are the slow L correct trials
% ContData.behavior.wheelsRfirstSlowRxn = wheelsRfirstValid(RTrialSlows); %nope


%% BS 2.1.16 - based on ns4 file - commented out 3.21.16 since I can use csv (taskbase file now)
% 
% % Find any leftward movement (-going voltages) or rightward movement (+)
%   % Find "center" of raw wheel voltages
%     centers = median(sWheelData); 
%     wheelMin = min(sWheelData);    %just to see
%     wheelMax = max(sWheelData);    %just to see
% 
%     wheelData = single(sWheelData); 
%     wheelData = sWheelData(2,:);
%     wheelDataInds = [];
%     for w = 1:length(wheelData);
%         wheelDataInds(:,w) = w;
%     end
% 
%     newWheelMat(1,:) = wheelDataInds;
%     newWheelMat(2,:) = wheelData;
% 
%     % L movement extractions:
%     allL_wheel = [];
%     allL_wheelInds = [];
%     for z = 1:length(wheelData)
%         if wheelData(z) > wheelData(z+1)
%             if z+1 == length(wheelData)
%                 break;
%             else
%                 allL_wheel(z) = wheelData(z);
%                 allL_wheelInds(z) = wheelDataInds(z);
%             end
%         end
%     end
% 
%     % Eliminate wheel recordings that occurred before first trial start (.nev)
% %     validWheel_L = find(allL_wheelInds' >= trialStartsNew(1));            % These are recording times(ms) when wheel is moving left; pull out first vals...
%     validWheel_L = find(allL_wheelInds' >= ContData.behavior.trialStartsValid(1));
%     
%     validsL = [];
%     validsL_first = [];
%     for v = 1:length(validWheel_L)
%         if v == numel(validWheel_L)
%             break;
%         else
%             validsL(v+1) = validWheel_L(v+1) > validWheel_L(v) + 5;         % Cut-off is 5ms;
%         end
%     end
% 
%     validsL_first = validWheel_L(validsL == 1);                             % Here are the first TSs for leftward vectors
%     
%     %Write the vectors for all L wheel movements to the ContData structure
% %     ContData.behavior.validWheel_L = validWheel_L;
%     ContData.behavior.newWheelMat = newWheelMat;
% %     ContData.behavior.allL_wheel = single(allL_wheel);                    % Great, but these aren't necessary
% %     ContData.behavior.allL_wheelInds = single(allL_wheelInds);
%     
%     allL_wheelIndsShort = find(allL_wheelInds>0);
% %     ContData.behavior.allL_wheelIndsShort = single(allL_wheelIndsShort); 
% %     validWheel_L = find(allL_wheelIndsShort' >= trialStartsNew(1));       %These are recording times(ms) when wheel is moving left (without the zeros)
%     ContData.behavior.validsL_first = validsL_first;
% 
%     %R movement extractions:
%     allR_wheel = [];
%     allR_wheelInds = [];
%     zrlim = numel(wheelData(1:end-10))
% 
%     for zr = 1:length(wheelData)
%         if wheelData(zr+1) > wheelData(zr)
%             if zr == zrlim;
%                 break;
%             else
%                 
%                 allR_wheel(zr) = wheelData(zr);
%                 allR_wheelInds(zr) = wheelDataInds(zr);
%                 
%             end
%         end
%     end
% 
%     % Eliminate wheel recordings that occurred before first trial start (.nev)
%     validWheel_R = find(allR_wheelInds' >= ContData.behavior.trialStartsValid(1));  %These are recording times(ms) when wheel is moving left; pull out first vals...
%     
%     validsR = [];
%     validsR_first = [];
%     for vr = 1:length(validWheel_R)
%         if vr == numel(validWheel_R)
%             break;
%         else
%             validsR(vr+1) = validWheel_R(vr+1) > validWheel_R(vr) + 5;              %cut-off is 5ms
%         end
%     end
% 
%     validsR_first = validWheel_R(validsR == 1);                                     %Here are the first TSs for rightward vectors
%     ContData.behavior.validsR_first = validsR_first;
    
 %% EXTRACT VELOCITY DATA FROM CONT LEVER DATA. BS: motor in: sLeverV = unfiltered velocities; sLeverVm = filtered
% 
%     numSamples = size(ContData.behavior.sLeverData,2);
%     tmpLeverData(1,:) = sgolayfilt(ContData.behavior.sLeverData(1,:),3,151); 
%     tmpLeverData(2,:) = sgolayfilt(ContData.behavior.sLeverData(2,:),3,151);
% 
%     sLeverV = zeros(1,numSamples);
% 
%     disp(' ');disp(' ');disp('Extracting velocity...');
% 
%     dX = diff(tmpLeverData(1,:));
%     dY = diff(tmpLeverData(2,:));
% 
%     sLeverV = sqrt( dX.^2 + dY.^2 );
% 
%     disp(' ');disp(' Complete. ');disp(' ');
% 
%     ContData.behavior.sLeverV = sLeverV;
%     ContData.behavior.sLeverVm = sgolayfilt(ContData.behavior.sLeverV,3,501); %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
%     
% %     figure(1); plot(sLeverV);
%     
%     clear sLeverV;
% 
% %% FIND MOVEMENTS
%     method = 'vel';
%     numC = 5;
%     clear progSt* reach
% 
%     currX       = ContData.behavior.sLeverData(1,:);
%     currY       = ContData.behavior.sLeverData(2,:);
%     currV       = ContData.behavior.sLeverVm; 
% 
%     pre     = 10;
%     post    = 10;       
%     minSpace = 250;
%     count = 1;
% 
%     % threshold the velocities
%     switch dataRate
%         case 'ns4'                                                            %currV are the velocity values
%             allValidSamps   = find(currV(1:numel(currV)-(minSpace.*2))>0.35); %find point where filtered velocities begin at -500 and >.35 (will be that index and above)
%         case 'ns2'
%             allValidSamps   = find(currV(1:numel(currV)-(minSpace.*2))>1.5);
%     end
%     
%     if numel(allValidSamps)>0
% 
%         switch method
% 
%             case 'vel'
% 
%     %             progStartTMP(count,1)  = currStamps(allValidSamps(1));
%                 progStartTMP(count,2)  = currX(allValidSamps(1));
%                 progStartTMP(count,3)  = currY(allValidSamps(1));           %find currY positions for velocities >.35
%                 progStartTMP(count,4)  = allValidSamps(1);
% 
%                 for j=2:numel(allValidSamps)                                %allValidSamps are the indexed velocities >.35
% 
%                     if allValidSamps(j)>allValidSamps(j-1)+minSpace         %if the indexed velocities are increasing by at least 250
% 
%                         postN = find(currV(allValidSamps(j-1):allValidSamps(j))<1,1,'first');  %postN = find the first velocity values at those positions
%                         if numel(postN)<1
%                             figure(3); clf;                                                         %if postN has less than 1 value, plot                
%                             plot(progStartTMP(count,4):allValidSamps(j),currV(progStartTMP(count,4):allValidSamps(j)),'k');  
%                             hold on; plot(progStartTMP(count,4),currV(progStartTMP(count,4)),'bo',allValidSamps(j-1),currV(allValidSamps(j-1)),'ro');
%                             pause(0.1);
% 
%                             disp('Cannot find stop');                                               %and set postN = post (which was defined as 10)
%                             postN=post;
%                         end
%     %                     progStopTMP(count,1)   = currStamps(allValidSamps(j-1)+post);
%                         progStopTMP(count,2)   = currX(allValidSamps(j-1)+postN);                %But if postN has >1 value, create a matrix (progStopTMP) and find currY at index allValidSamps at the previous positions and add 10 
%                         progStopTMP(count,3)   = currY(allValidSamps(j-1)+postN);
%                         progStopTMP(count,4)   = allValidSamps(j-1)+postN;
%                         count                  = count+1;
% 
%                         preN = find(currV(allValidSamps(j-1) : allValidSamps(j))<1,1,'last');
%                         if numel(preN)<1
%                             disp('Cannot find start');
%     %                     progStartTMP(count,1)  =
%     %                     currStamps(allValidSamps(j)-pre)                                        %Now create the startTMP matrix
%                             progStartTMP(count,2)  = currX(allValidSamps(j)-pre);
%                             progStartTMP(count,3)  = currY(allValidSamps(j)-pre);
%                             progStartTMP(count,4)  = allValidSamps(j)-pre;                    
%                         else
%     %                     progStartTMP(count,1)  = currStamps(allValidSamps(j-1)+preN);
%                             progStartTMP(count,2)  = currX(allValidSamps(j-1)+preN);
%                             progStartTMP(count,3)  = currY(allValidSamps(j-1)+preN);
%                             progStartTMP(count,4)  = allValidSamps(j-1)+preN;
%                         end
%                     end
% 
%                     if j==numel(allValidSamps)
%     %                     post = find(currV(allValidSamps(j):allValidSamps(j)+minSpace)<0.5,1,'first');
%     %                     progStopTMP(count,1)   = currStamps(allValidSamps(j)+post);
%                         progStopTMP(count,2)   = currX(allValidSamps(j)+post);
%                         progStopTMP(count,3)   = currY(allValidSamps(j)+post);
%                         progStopTMP(count,4)   = allValidSamps(j)+post;
%                     end
% 
%                 end
% 
%                 count = 1;
%                 for k = 1:size(progStartTMP,1)
% 
%                     if k==1
%                         reach.init = 1;
%                     end
% 
%                     % reaches must be at least 50 ms long
%                     if progStopTMP(k,4)-progStartTMP(k,4)>=90 & progStartTMP(k,4)>minSpace
% 
%                         trajAngle   = atan2(progStopTMP(k,3)-progStartTMP(k,3),progStopTMP(k,2)-progStartTMP(k,2));
% 
%                         if (pdist2([progStopTMP(k,2),progStopTMP(k,3)],[mean(currX),mean(currY)]) > pdist2([progStartTMP(k,2),progStartTMP(k,3)],[mean(currX),mean(currY)]))
%                             reach.out(count) = 1;
%                         else
%                             reach.out(count) = 0;
%                         end
%                         velTraj = ContData.behavior.sLeverV(progStartTMP(k,4) : progStopTMP(k,4));
%                         xVals = ContData.behavior.sLeverData(1,progStartTMP(k,4) : progStopTMP(k,4));
%                         yVals = ContData.behavior.sLeverData(2,progStartTMP(k,4) : progStopTMP(k,4));
% 
%                         reach.start(count,:)  = progStartTMP(k,:);
%                         reach.stop(count,:)   = progStopTMP(k,:);
%                         reach.angle(count,1)  = trajAngle;
%                         reach.dist(count,1)   = trapz(velTraj);
%                         reach.dist(count,2)   = pdist2(progStartTMP(k,2:3) , progStopTMP(k,2:3));
% 
%                         tmp = findpeaks(velTraj);
%                         reach.numpks(count,1) = numel(tmp.loc);
%                         reach.dur(count,1)    = progStopTMP(k,4) - progStartTMP(k,4);
%                         reach.vel(count,1)   = max(velTraj);
%                         reach.vel(count,2)   = trapz(velTraj) ./ reach.dur(count,1);
%                         reach.vel(count,3)   = var(velTraj);
%                         reach.vel(count,4)   = find(velTraj==max(velTraj),1);
% 
%                         reach.acc(count,1)   = max(diff(velTraj));
%                         reach.acc(count,2)   = mean(diff(velTraj));
%                         reach.acc(count,3)   = max(diff(velTraj(1:90))); % max in first 90 ms of movement
% 
%                         reach.tort(count,1)  = reach.dist(count,1) ./ pdist2([progStopTMP(k,2),progStopTMP(k,3)],[progStartTMP(k,2),progStartTMP(k,3)]);
% 
%                         % find max displacement of the reach
%                         xVals = xVals - xVals(1);
%                         yVals = yVals - yVals(1);
%                         
%                         
%                         
%                         valThrInd = find(ContData.behavior.threshInds>reach.start(count,4) & ContData.behavior.threshInds<reach.stop(count,4));
%                         if numel(valThrInd)>0
%                             reach.rewarded(count)= 1;
%                         else
%                             reach.rewarded(count)= 0;                       
%                         end
% 
%                         count                 = count+1;
% 
%                     end            
%                 end
% 
% 
%         end
% 
%         disp(['Num valid reaches: ' num2str(count-1)]);
% 
%         level = ones(1,size(progStartTMP,1)).*0.35;
%         figure(3); clf;
%         plot(ContData.behavior.sLeverV,'k'); title(['sLeverV']);hold on;
%         plot(progStartTMP(:,4),level,'r^'); title(['reach starts']);
%         plot(progStopTMP(:,4),level,'bo');  title(['reach stops']);
% 
%         [vals,inds] = sort(reach.dur,'ascend');    
% 
%         reach.numReaches = size(reach.start,1);
% 
%         dims = floor(sqrt(reach.numReaches));
%         scaler = 3000;
% 
%         if dispOn
%             
%             figure(5); clf;
%         
%             for p = 1:1:dims.^2
% 
%                 l = inds(p);
% 
%                 offset = l.*ones(1,numel(reach.start(l,4):reach.stop(l,4)));
%                 [x,y] = ind2sub(dims,p);
% 
%                 figure(5);
%                 plot((x.*scaler),(y.*scaler),'k.','MarkerSize',8);
%                 hold on;
%                 plot(currX(reach.start(l,4):reach.stop(l,4))-currX(reach.start(l,4))+(x.*scaler),currY(reach.start(l,4):reach.stop(l,4))-currY(reach.start(l,4))+(y.*scaler) , 'k','LineWidth',1,'Color',[p/reach.numReaches 0.67-(p/reach.numReaches)*0.67 1-(p/reach.numReaches)]);
%                 if i==dims.^2
%                     axis tight; axis off;
%                 end
% 
%                 startInd = round(reach.start(l,4)./2);
%                 stopInd = round(reach.stop(l,4)./2);
% 
%             end
%         end
%         
%         if reach.numReaches>5
%                 figure(7); clf;
%                 subplot(311);
%                 p = polyfit(reach.dist(:,1),reach.vel(:,1),1);
%                 plot(reach.dist,reach.vel(:,1),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
%                 xlabel('Total distance');
%                 ylabel('Peak Velocity');
% 
%                 subplot(312);
%                 p = polyfit(reach.dist(:,1),reach.dur,2);
%                 plot(reach.dist,reach.dur,'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
%                 xlabel('Total distance');
%                 ylabel('Movement Duration (ms)');
% 
%                 subplot(313);
%                 p = polyfit(reach.dist(:,1),reach.acc(:,3),1);
%                 plot(reach.dist,reach.acc(:,3),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
%                 xlabel('Total distance');
%                 ylabel('Initial acceleration');
% 
%                 figure(8); clf;
%                 rose(reach.angle);
% 
%                 drawnow;
%         end
%         
%     else
%         
%         reach.numReaches=0;
%     
%     end
%     
% 
%     
% %% DIMENSION REDUCTION OF REACH PARAMETERS
% 
% if reach.numReaches > 5
%     clear matFor*
% 
%     matForDR = [reach.vel reach.acc reach.tort reach.dur reach.dist];
% 
%     for j=1:size(matForDR,2)
%         matForCov(:,j) = (matForDR(:,j) - mean(matForDR(:,j))) ./ std(matForDR(:,j)) ;
%     end
% 
%     figure(11);
%     imagesc(corr(matForCov),[0 1]);
%     map = colormap(TNC_CreateRBColormap(1024,'cpb'));
% %     colormap(1-map);
% 
%     [eV,eD] = eig(cov(matForCov));
%     figure(12); plot(cumsum(diag(eD))./sum(diag(eD)),'ko-');
% 
%     for m=1:size(matForDR,1)
%         reach.pca(1,m) = dot(matForCov(m,:),eV(:,size(matForCov,2))');
%         reach.pca(2,m) = dot(matForCov(m,:),eV(:,size(matForCov,2)-1)');
%         reach.pca(3,m) = dot(matForCov(m,:),eV(:,size(matForCov,2)-2)');
%     end
% end
% 
% %% CALCULATE META-REACH PARAMETERS
%     if reach.numReaches > 5
%         [winReach] = TNC_ReachVigorWindow(reach.pca(1,:),reach.numReaches,9); % capture fluctuations along the maximally variant dimension
%         ContData.behavior.winReach = winReach;
%     end
%     
% %% WAS THE REACH REWARDED?
% if reach.numReaches > 5
%     for p=1:reach.numReaches
% 
%         tmp = find(ContData.behavior.threshInds>reach.start(p,4) & ContData.behavior.threshInds<reach.stop(p,4));
%         if numel(tmp)>0
%             reach.rewarded(p)=1;
%         else
%             reach.rewarded(p)=0;
%         end
%         
%     end
% end    
% %% WRITE THE REACH STRUCTURE OF THE CONTDATA STRUCTURE
%     ContData.behavior.reach = reach;
    
%% COMPLETED ALL ANALYSIS. SAVE AND START OVER

    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Completed file: ' filenamestr(1,1:length(filenamestr)-3)]);
    
    filenamestrB = strcat(newPath, filenamestrE);
     
     cd(fpath);
     
     save([targetName '_bh.mat'],'ContData');                            
     disp(['saved as ' targetName '_bh.mat']);
                                                                     
    disp('%-------------------------------------------------------------------');
    disp(' ');
    
