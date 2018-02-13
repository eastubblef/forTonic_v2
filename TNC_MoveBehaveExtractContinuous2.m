
 function [ContData, filenamestrB] = TNC_MoveBehaveExtractContinuous2(taskbase, filenamestr, targetName, dataRate, chan, elimTime, fpath)
%% This is Beth's working version of Josh's extraction code from Tonic_v2 

% UPDATED LAST:
% 3.08.17  - pull out fast rx times v slow ones for movement alignment
% 11.28.16 - consider aligning a few trials out for re-starts on csv clock(l. 50ish); usually using 1st TS of CSV works, but find first trial after an obvious restart on the ns4 voltage trace
% 11.24.16 - set T-hold for ns4 continuous signal's trial start pulses (l. 84ish) 
% 11.22.16 - eff the NEV signals' abberant and finicky nature; use ns4 file to pick up trial start and solenoid pulses for clock alignment to behavior
% 11.17.16 - solenoid TSs from NEV file vary widely. See line 98-102!
% 11.3.16 - to fit with better behaveLoad6.m file for taskbase extraction
% 3.21.16 - 3.29.16 use taskbase file to align TS of spikes to trial start & L v R movements (behaveLoad5)
% ll.51-52 are lickInds that I commented out without piezo-data 11.26.15
% ll. 79 & 84, 113-117 commented out for chan.lick 11.26.15

%  Note: removed the following vars since that is empty in my Blackrock system (7.16.15): NEVfiles
%  chan.x
%  ContData.behavior.threshInds = TS for laser pulse-on (divide by 30, 30kHz = sampling frequency for spikes)
% l. 51-76: chan.rew (129), funky: (NOT only proper sol reward inds) - mixed with trial starts

% saves a bh.mat file for continuous behavioral data
% this mfile is called from MoverScript_continuous.m

%% FILE NAMES
    
    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Beginning file: ' filenamestr(1,1:length(filenamestr)-3)]);
    disp(' ');disp(' ');  

    dispOn = 0;
    
    newPath = fpath;
    cd(newPath);
    
%% Load taskbase structure for behavioral data 
                      
    if exist('filenamestr', 'var');
        [filenamestrT, path] = uigetfile('*_tb.mat*','select the taskbase file', fpath)       % get the actual # correct trials (and TS) from csv file
    end
    taskbase = strcat(path, filenamestrT);
    load(taskbase);
    startTimesCsv = taskbase.trialStartTimes';
    startsRewsCsv = taskbase.startsRews';                                                            % trial starts and solenoid discharges from behavior file
    rewTimesCsv = taskbase.rewTimes';
    
    diffStrtsCSV = diff(startTimesCsv);                                                              % use for alignment to ns4 files
    
    % Zero out csv trial starts for comparison/alignment with nev file
     validStrtsRewsCsv1 = find(diff(startsRewsCsv)>50);                                               % account for very low vals: rx time that is virtually impossible due to mouse likely already moving
     validStrtsRewsCsv = startsRewsCsv(validStrtsRewsCsv1);                                           % align this clock to nev
     startsRewsZeroCsv = validStrtsRewsCsv-validStrtsRewsCsv(1);                                      % zero out TS 1; first trial start = 0
     
     rewTimesCsvZero = rewTimesCsv-validStrtsRewsCsv(1);  
%      rewTimesZeroCsv = rewTimesCsv-validStrtsRewsCsv(1);  
     startTimesZeroCsv = startTimesCsv - startTimesCsv(1);
       
     %for particularly hard-to-align sessions (due to restarts), check this too:
%      startTimesZeroCsv = startTimesCsv - startTimesCsv(2);
%      startTimesZeroCsv = startTimesZeroCsv(2:end);                              %must get zero as first element
     

%% LOAD EVENT DATA 
    newPath = fpath;
    cd(newPath);

%      if exist('filenamestr', 'var');
         [filenamestrE, path] = uigetfile('*.nev*','select the .nev file', fpath)               
%      end
    
    dataEvents = openNEVBS_laser(filenamestrE,'read','nosave','nomat', 'report');                    

%% LOAD CONTINUOUS DATA
    
 switch dataRate
    
    case 'ns4'    
        [filenamestrE, path] = uigetfile('*.ns4','select the .ns4 file', fpath); 
        Ns4DATA = openNSxBS('report','read',filenamestrE);

        % How Josh's variables go down        
%         xChan = find(Ns4DATA.MetaTags.ChannelID==chan.x);                 %my x channel is empty  
%         yChan = find(Ns4DATA.MetaTags.ChannelID==chan.y);                 %motor in; then lick; now motor
%         lChan = find(Ns4DATA.MetaTags.ChannelID==chan.lick);
        thrChan = find(Ns4DATA.MetaTags.ChannelID==chan.thr);
        yChan =   find(Ns4DATA.MetaTags.ChannelID==chan.y);                 %lick channel - nope. Motor in
        trialStrtsChan = find(Ns4DATA.MetaTags.ChannelID==chan.rew);        %trial start/sol channel
%         leverData(1,:) = decimate(Ns4DATA.Data(xChan,:),10);
%         leverData(2,:) = decimate(Ns4DATA.Data(yChan,:),10);              %correct for 10kS/s rate
%         rawLick = decimate(Ns4DATA.Data(lChan,:),10);
        
        wheel = decimate(Ns4DATA.Data(yChan,:),10);                         
        trialStrtsNSol = decimate(Ns4DATA.Data(trialStrtsChan,:),10);       %convert these to ms (/10)
        lick = decimate(Ns4DATA.Data(thrChan,:),10);
        
 end     
%         sLeverData(1,:) = sgolayfilt(leverData(1,:),9,101);                   %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
     trialStrtsNSolData_sgolayF = sgolayfilt(trialStrtsNSol,9,101);             %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
     trialStrtsNSolData = trialStrtsNSol;                                       %likely too raw

     ContData.behavior.trialStrtsNSolData = trialStrtsNSolData; 
     figure; hold on;
     plot(trialStrtsNSolData(1:150000));
%      plot(trialStrtsNSolData(60000:100000));
%      plot(trialStrtsNSolData);
     
%      figure; hold on;
%      plot(lick(1:150000));
%      plot(lick);

%      figure; hold on;
%      wheel_sgolayFilt = sgolayfilt(wheel,9,101);
%      plot(wheel_sgolayFilt(1:50000), 'r');  %but this will essentialy apply a filter that was already applied during line noise cancellation

%      figure; hold on;
     plot(wheel(1:150000));
%      plot(wheel);

%      ContData.behavior.trialStrtsNSolData_sgolayF = trialStrtsNSolData_sgolayF;  %#points are = numel(sWheelData)/10 (and rounded up) for correction of 10kHz sampling rate

     startsRewNs4 = trialStrtsNSolData;
     
%% New 5.1.17: in order to (more) reliably estimate trial starts
%  xThold will be the length of time that the trial start "top of the gate" length is relevant - always ~ 102 ms
%  
%      xThold = 100;  %161015 session
%      yTholdMin = 240;  yTholdMax = 320      %161005 session
%      xThold = 100;  xTholdMax = 115         %161005 session
%      yTholdMin = 200;  yTholdMax = 900;     %170112 session
%      xThold = 100;  %xTholdMax = 115        %170112 session   
     yTholdMin = 200;  yTholdMax = 600;       %170111 session
%      yTholdMin = 5500;  yTholdMax = 5600;   %170321 session
%      xThold = 100;  %xTholdMax = 115          %170112 session   
%      yTholdMin = 150;  yTholdMax = 900;     %170501testBoth session
     xThold = 90;                           %maybe use a lower t-hold value since the first few may not get captured
     firstT = find(startsRewNs4 > yTholdMin & startsRewNs4 < yTholdMax);    %find the tall gaps = ea. trial start; need a Thold > any sol signal
     vecT = [];
     for i = 1:length(firstT)
         if i+1 == length(firstT)
             break
         end
         if firstT(i+1) == firstT(i) + 1
             vecT(i) = firstT(i);                                           %when the sliding vector breaks, a zero will appear;
         end
     end                                                                    %Find the span of >0 vals = 113 (ms); the "top of the gate" will last approx. 102 ms
     Tvec = [];
     for i = 1:length(vecT)
         if i+xThold == length(vecT)
             break
         end
         if vecT(i+xThold) == vecT(i) + xThold;
             Tvec(i:i+xThold) = vecT(i):vecT(i) + xThold;
         end
     end                                                                    %first trial should be the 1st element of second vector in which 100ish vals > 0;
     FTinds = Tvec > 0;
     newVecT = [];
     for i = 1:length(FTinds)
         if FTinds(i) == 0 && FTinds(i+1) == 1
            newVecT(i) = FTinds(i+1);
         end
     end  
     Tinds = find(newVecT == 1);                                            %these should be trial start times from 2nd element on.
%      firstTrialInd = Tinds(2);
%      firstTrialStartNS4 = vecT(firstTrialInd + 1);
     allTrialsNS4 = vecT(Tinds+1);                                          
     diffTrialStartsNS4 = diff(allTrialsNS4');

    %% Here comes the variability among datasets; figure out more elegant way to do this by directly clicking on the voltage trace?
     %Or use diffTrials CSV compared to NS4:
     
     if diffStrtsCSV(1) >= diffTrialStartsNS4(3)-3 || diffStrtsCSV(1) <= diffTrialStartsNS4(3)+3;      %give it 3 ms buffer
         refinedTrialStartsNS4 = allTrialsNS4(3:end);
         if diffStrtsCSV(1) >= diffTrialStartsNS4(2)-3 || diffStrtsCSV(1) <= diffTrialStartsNS4(2)+3
             refinedTrialStartsNS4 = allTrialsNS4(2:end);
         else if diffStrtsCSV(1) >= diffTrialStartsNS4(1)-3 || diffStrtsCSV(1) <= diffTrialStartsNS4(1)+3
                 refinedTrialStartsNS4 = allTrialsNS4(1:end);
             end
         end
     end
    
    diffRefinedTrialStartsNS4 = diff(refinedTrialStartsNS4');  
    forCompareNS4lowLim = diffRefinedTrialStartsNS4(1:10)-1;
    forCompareNS4highLim = diffRefinedTrialStartsNS4(1:10)+1;
    forCompareCSV = diffStrtsCSV(1:10); 
    
    compareDiffs = forCompareNS4 == forCompareCSV | forCompareNS4 == forCompareCSV+1 | forCompareNS4 == forCompareCSV-1; 

    if compareDiffs(1) == 1
     trialStartsZeroNS4 = refinedTrialStartsNS4' - refinedTrialStartsNS4(1)';           
    end

%% confirm no extraneous trials/left out trials:     
    clear i
    trialStrts = [];
    for i = 1:length(startTimesZeroCsv)
        for j = 1:length(trialStartsZeroNS4)
            if i == 1
                trialStrts(i) = 1;                                               % report a 1 for the first trial alignment (and since the first trial start for both clocks (by design) is defined as aligned) 
            else
                if trialStartsZeroNS4(j) > startTimesZeroCsv(i)-50 && trialStartsZeroNS4(j) < startTimesZeroCsv(i)+50; % creates a window to search within 40 ms +/-
                    trialStrts(j) = trialStartsZeroNS4(j);                         % indexed positions are true to the zeroed (to first true trial start) nev clock
                end                                                             
            end                                                                  % then, after this, index to trial starts csv and find L mvments v R mvments per trial
        end
    end
    trialStartsPre = find(trialStrts >= 1);
    trialStarts = trialStrts(trialStartsPre);         % both clocks aligned
%     trialStartsWInds(1,:) = 1:length(trialStarts);  % updated for Rx time/moveOn alignment to keep up with trial#
%     trialStartsWInds(2,:) = trialStarts;            % updated for Rx time/moveOn alignment to keep up with trial#
    
    ContData.behavior.trialStarts = trialStarts;
    
    leftIncorrectTimes = taskbase.leftIncorrectTimes;
    leftIndsCenteredTimes = taskbase.leftIndsCenteredTimes;

    lCenteredTimesCsv = taskbase.leftCenteredTimes;
    lCenteredZeroCsv = lCenteredTimesCsv - startTimesCsv(1);                     % zero-out l/r center times csv vector to index back into zeroed nev file
    lIncorrectZeroCsv = leftIncorrectTimes - startTimesCsv(1);
    
    % alignment of solenoid clicks/reward times from csv to blackrock
    rewTimes = [];
    for m = 1:length(trialStarts)   %ns4
        for n = 1:length(rewTimesCsvZero)
            if m == numel(trialStarts) %ns4
                break
            end
            if n == numel(rewTimesCsvZero)
                break
            end
            if rewTimesCsvZero(n) > trialStarts(m) && rewTimesCsvZero(n) < trialStarts(m+1) && rewTimesCsvZero(n) < rewTimesCsvZero(n+1) % just to be sure
                rewTimes(n) = rewTimesCsvZero(n);
            end
        end
    end
%     validsCorrection = single(validStrtsRews(1));
%     validsCorr = single(validStrtsRewsCsv(1));
%     rewTimesValid = rewTimes + validsCorr;
   rewTimesValid = rewTimes;   %ns4
   ContData.behavior.rewTimes = rewTimes;
   
    % sanity check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    lIndsCenteredTimesZero = [];                                                 % zero-out l center times, keeping indexed position
    for x = 1:length(leftIndsCenteredTimes);
        if leftIndsCenteredTimes(x) > 1
            lIndsCenteredTimesZero(x) = leftIndsCenteredTimes(x) - startTimesCsv(1);
        end
        if leftIndsCenteredTimes(x) == 0
            lIndsCenteredTimesZero(x) = 0;                                       % this vector should be the same as lCenteredZeroCsv except at correct index position for trial number
        end
    end
    
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
               lCenters(a) = trialStarts(a);
            end
        end
    end
    lCenters1 = lCenters > 0; LrewTrialStartsZeroed = lCenters(lCenters1);
 
    %~~~~~~~~~~~~~~~Be aware that these may change if I correct for multiple timestamps into CSV~~~~~~~~~~~~~~~~~~~~~~~~  
%     validsCorrection = single(validStrtsRews(1));
%     trialStartsValid = trialStarts + validsCorrection;
    LrewTrialsStarts = LrewTrialStartsZeroed;
    LrewTrialStarts = LrewTrialStartsZeroed;

    % L incorrect trials' start times
    lIncorrects = [];
    for c = 1:length(trialStarts)
        for d = 1:length(lIncorrectZeroCsv)
            if c == numel(trialStarts)
                break
             else if d == numel(lIncorrectZeroCsv)
                    break
                else if lIncorrectZeroCsv(d) > trialStarts(c) && lIncorrectZeroCsv(d) < trialStarts(c+1) && lIncorrectZeroCsv(d) < lIncorrectZeroCsv(d+1)
                        lIncorrects(c) = trialStarts(c);
                    end
                end
            end
        end
    end
    lIncorrects1 = lIncorrects > 0; LincorrectsTrialStartsZeroed = lIncorrects(lIncorrects1);
%     LincorrectTrialsStarts = LincorrectsTrialStartsZeroed + validsCorrection;
    LincorrectTrialsStarts = LincorrectsTrialStartsZeroed;
   
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % right side trials
    rCenteredTimesCsv = taskbase.rightCenteredTimes;
    rightIncorrectTimes = taskbase.rightIncorrectTimes;
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
    rCenters1 = rCenters > 0; RrewTrialStartsZeroed = rCenters(rCenters1);
%     RrewTrialsStarts = RrewTrialStartsZeroed + validsCorrection; 
    RrewTrialsStarts = RrewTrialStartsZeroed; 

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
%     RincorrectTrialsStarts = RIncorrectTrialStartsZeroed + validsCorrection;
    RincorrectTrialsStarts = RIncorrectTrialStartsZeroed;

    ContData.behavior.LrewTrialStarts = LrewTrialsStarts;
%     ContData.behavior.trialStartsValid = trialStartsValid;
    ContData.behavior.trialStartsValid = trialStarts;
    ContData.behavior.RrewTrialStarts = RrewTrialsStarts;
    ContData.behavior.LincorrectTrialStarts = LincorrectTrialsStarts;
    ContData.behavior.RincorrectTrialStarts = RincorrectTrialsStarts;
%     ContData.behavior.rewTimesValid = rewTimesValid;
%     
%     rewardedTrialsNev = numel(ContData.behavior.LrewTrialStarts) + numel(ContData.behavior.RrewTrialStarts);
%     disp(['Number of rewarded trials (.nev): ' num2str(rewardedTrialsNev)]);
%     disp(' ');disp(' '); 
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Load the taskbase data relevent to wheel movement extractions L v R. Could compare to ns4 file. 
    % Note: these values still retain 0s to hold indices where L or Rtrials were either not rewarded, or were rewarded to other side.
    
    % Correct L and R trials:
    wheelLfirstTime = taskbase.wheelLfirstTime;
    wheelLfirst = taskbase.wheelLfirst;
    wheelRfirst = taskbase.wheelRfirst;
    wheelRfirstTime = taskbase.wheelRfirstTime;
    
    wheelLfirstTimeFastest = taskbase.wheelLfirstTimeFastest;
    wheelRfirstTimeFastest = taskbase.wheelRfirstTimeFastest;
    wheelLfirstTimeMids = taskbase.wheelLfirstTimeMids;
    wheelRfirstTimeMids = taskbase.wheelRfirstTimeMids;
    wheelLfirstTimeSlow = taskbase.wheelLfirstTimeSlow;
    wheelRfirstTimeSlow = taskbase.wheelRfirstTimeSlow;

    % Remove the empty indices and extract just the first vals needed and use LrewTrialStarts/RewTrialStarts for alignment
    wheelLfirstT = wheelLfirstTime > 0; wheelLfirsts = wheelLfirstTime(wheelLfirstT);
    wheelLfirstsZero = wheelLfirsts - startTimesCsv(1);
    
    wheelLfirstTFastest = wheelLfirstTimeFastest > 0; wheelLfirstsFastest = wheelLfirstTimeFastest(wheelLfirstTFastest);
    wheelLfirstsFastestZero = wheelLfirstsFastest - startTimesCsv(1);

    wheelLfirstTMids = wheelLfirstTimeMids > 0; wheelLfirstsMids = wheelLfirstTimeMids(wheelLfirstTMids);
    wheelLfirstsMidsZero = wheelLfirstsMids - startTimesCsv(1);

    wheelRfirstT = wheelRfirstTime > 0; wheelRfirsts = wheelRfirstTime(wheelRfirstT);
    wheelRfirstsZero = wheelRfirsts - startTimesCsv(1); 
    
    wheelRfirstTFastest = wheelRfirstTimeFastest > 0; wheelRfirstsFastest = wheelRfirstTimeFastest(wheelRfirstTFastest);
    wheelRfirstsFastestZero = wheelRfirstsFastest - startTimesCsv(1); 

    wheelRfirstTMids = wheelRfirstTimeMids > 0; wheelRfirstsMids = wheelRfirstTimeMids(wheelRfirstTMids);
    wheelRfirstsMidsZero = wheelRfirstsMids - startTimesCsv(1); 

    extractWheelCsvLR = numel(wheelLfirsts) + numel(wheelRfirsts);
%     disp(['Number of L/R wheel movement extractions (csv): ' num2str(extractedWheelCsv);
%     disp(' '); disp(' ');
 
%% Get this working later:
    % Incorrect L and R trials:
    wheelILfirstTime = taskbase.wheelILfirstTime;
    wheelILfirst = taskbase.wheelILfirst;
    wheelIRfirst = taskbase.wheelIRfirst;
    wheelIRfirstTime = taskbase.wheelIRfirstTime;

    rightIncorrectTimes = taskbase.rightIncorrectTimes;
    wheelILfirstT = wheelILfirstTime > 0; wheelILfirsts = wheelILfirstTime(wheelILfirstT);
    wheelILfirstsZero = wheelILfirsts - startTimesCsv(1);
    wheelIRfirstT = wheelIRfirstTime > 0; wheelIRfirsts = wheelIRfirstTime(wheelIRfirstT);
    wheelIRfirstsZero = wheelIRfirsts - startTimesCsv(1);

    extractWheelCsvILnIR = numel(wheelILfirsts) + numel(wheelIRfirsts);
    
%% aligned (to ns4-based trial starts) wheel movements that were based on csv file
%  May want to consider including this in the overall trial alignment for absolute certainty - rx times (diffL, diffR) are on super-fast side for some trials (25 ms)
    
    % correct L movements: eliminates first L movements that don't fit to an extracted trial start into Blackrock
    wheelLfirstsAligned = [];
    for i = 1:length(LrewTrialStartsZeroed)
        if numel(LrewTrialStartsZeroed) < numel(wheelLfirstsZero) && i == numel(LrewTrialStartsZeroed)
            break
        else if i == numel(wheelLfirstsZero)
                break
            else if wheelLfirstsZero(i) > LrewTrialStartsZeroed(i) && wheelLfirstsZero(i) < LrewTrialStartsZeroed(i+1) && wheelLfirstsZero(i) < wheelLfirstsZero(i+1) % just to be sure
                    wheelLfirstsAligned(i) = wheelLfirstsZero(i);
                else if wheelLfirstsZero(i) < LrewTrialStartsZeroed(i) && wheelLfirstsZero(i+1) > LrewTrialStartsZeroed(i) && wheelLfirstsZero(i+1) < LrewTrialStartsZeroed(i+1)
                        wheelLfirstsAligned(i) = wheelLfirstsZero(i+1);
                    end
                end
            end
        end
    end
    wheelLfirstsFastestAligned = [];
    for i = 1:length(LrewTrialStartsZeroed)
        if numel(LrewTrialStartsZeroed) < numel(wheelLfirstsFastestZero) && i == numel(LrewTrialStartsZeroed)
            break
        else if i == numel(wheelLfirstsFastestZero)
                break
            else if wheelLfirstsFastestZero(i) > LrewTrialStartsZeroed(i) && wheelLfirstsFastestZero(i) < LrewTrialStartsZeroed(i+1) && wheelLfirstsFastestZero(i) < wheelLfirstsFastestZero(i+1) % just to be sure
                    wheelLfirstsFastestAligned(i) = wheelLfirstsFastestZero(i);
                else if wheelLfirstsFastestZero(i) < LrewTrialStartsZeroed(i) && wheelLfirstsFastestZero(i+1) > LrewTrialStartsZeroed(i) && wheelLfirstsFastestZero(i+1) < LrewTrialStartsZeroed(i+1)
                        wheelLfirstsFastestAligned(i) = wheelLfirstsFastestZero(i+1);
                    end
                end
            end
        end
    end
    wheelLfirstsMidsAligned = [];
    for i = 1:length(LrewTrialStartsZeroed)
        if numel(LrewTrialStartsZeroed) < numel(wheelLfirstsMidsZero) && i == numel(LrewTrialStartsZeroed)
            break
        else if i == numel(wheelLfirstsMidsZero)
                break
            else if wheelLfirstsMidsZero(i) > LrewTrialStartsZeroed(i) && wheelLfirstsMidsZero(i) < LrewTrialStartsZeroed(i+1) && wheelLfirstsMidsZero(i) < wheelLfirstsMidsZero(i+1) % just to be sure
                    wheelLfirstsMidsAligned(i) = wheelLfirstsMidsZero(i);
                else if wheelLfirstsMidsZero(i) < LrewTrialStartsZeroed(i) && wheelLfirstsMidsZero(i+1) > LrewTrialStartsZeroed(i) && wheelLfirstsMidsZero(i+1) < LrewTrialStartsZeroed(i+1)
                        wheelLfirstsMidsAligned(i) = wheelLfirstsMidsZero(i+1);
                    end
                end
            end
        end
    end

    
    wheelILfirstsAligned1 = [];
    for j = 1:length(LincorrectsTrialStartsZeroed)
        if numel(LincorrectsTrialStartsZeroed) < numel(wheelILfirstsZero)&& j == numel(LincorrectsTrialStartsZeroed)
            break
        end
        if j == numel(wheelILfirstsZero)
            break
        end
        if wheelILfirstsZero(j) > LincorrectsTrialStartsZeroed(j) && wheelLfirstsZero(j) < LincorrectsTrialStartsZeroed(j+1) && wheelILfirstsZero(j) < wheelILfirstsZero(j+1)
            wheelILfirstsAligned1(j) = wheelILfirstsZero(j);
        else if wheelILfirstsZero(j) < LincorrectsTrialStartsZeroed(j) && wheelILfirstsZero(j+1) > LincorrectsTrialStartsZeroed(j) && wheelILfirstsZero(j+1) < LincorrectsTrialStartsZeroed(j+1)
                wheelILfirstsAligned1(j) = wheelILfirstsZero(j+1);
            end
        end
    end 
    wheelILfirstsAligned2 = wheelILfirstsAligned1 > 1;
    wheelILfirstsAligned = wheelILfirstsAligned1(wheelILfirstsAligned2);
    
% Alternative for alignment of L movement vectors if above doesn't work 
  if numel(wheelLfirstsAligned) < numel(LrewTrialStartsZeroed - 50)         % if the previous loop left out 50 of the L rew trials, do this loop
    wheelLfirstsAligned = [];
    for i = 1:length(LrewTrialStartsZeroed)
        for j = 1:length(wheelLfirstsZero)
            if i == numel(LrewTrialStartsZeroed)
                break
            end
            if j == numel(wheelLfirstsZero)
                 else if wheelLfirstsZero(j) > LrewTrialStartsZeroed(i) && wheelLfirstsZero(j) < LrewTrialStartsZeroed(i+1) && wheelLfirstsZero(j) < wheelLfirstsZero(j+1) % just to be sure
                        wheelLfirstsAligned(j) = wheelLfirstsZero(j);
                end
            end
         end
      end
  end        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
    % R movements: works better this way for some file's R movements. Think about this.
    wheelRfirstsAligned = [];
    for i = 1:length(RrewTrialStartsZeroed)
        for j = 1:length(wheelRfirstsZero)
            if i == numel(RrewTrialStartsZeroed)
                break
            end
            if j == numel(wheelRfirstsZero)
                break
            end
             if wheelRfirstsZero(j) > RrewTrialStartsZeroed(i) && wheelRfirstsZero(j) < RrewTrialStartsZeroed(i+1) && wheelRfirstsZero(j) < wheelRfirstsZero(j+1) % just to be sure
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
%         
% % Alternative 2 for alignment of R movement vectors if above doesn't work
% %     if numel(wheelRfirstsAligned) < numel(RrewTrialStartsZeroed/2)         % if the previous loop left out 50 of the L rew trials, do this loop
% %         wheelRfirstsAligned = [];
% %         for i = 1:length(RrewTrialStartsZeroed)
% %             if wheelRfirstsZero(i) > RrewTrialStartsZeroed(i) && wheelRfirstsZero(i) < RrewTrialStartsZeroed(i+1) && wheelRfirstsZero(i) < wheelRfirstsZero(i+1) % just to be sure
% %                 wheelRfirstsAligned(i) = wheelRfirstsZero(i);
% %             else if wheelRfirstsZero(i) < RrewTrialStartsZeroed(i) && wheelRfirstsZero(i+1) > RrewTrialStartsZeroed(i) && wheelRfirstsZero(i+1) < RrewTrialStartsZeroed(i+1)
% %                     wheelRfirstsAligned(i) = wheelRfirstsZero(i+1);
% %                 end
% %             end
% %         end
% %     end   
%     if numel(wheelIRfirstsAligned) < numel(RIncorrectTrialStartsZeroed/2) 
%         for i = 1:length(RIncorrectTrialStartsZeroed)
%             if i == numel(RIncorrectTrialStartsZeroed)
%                 break
%             end
%             if wheelIRfirstsZero(i) > RIncorrectTrialStartsZeroed(i) && wheelIRfirstsZero(i) < RIncorrectTrialStartsZeroed(i+1) && wheelIRfirstsZero(i) < wheelIRfirstsZero(i+1) % just to be sure
%                 wheelIRfirstsAligned(i) = wheelIRfirstsZero(i);
% 
%             else if wheelIRfirstsZero(i) < RIncorrectTrialStartsZeroed(i) && wheelIRfirstsZero(i+1) > RIncorrectTrialStartsZeroed(i) && wheelIRfirstsZero(i+1) < RIncorrectTrialStartsZeroed(i+1)
%                     wheelIRfirstsAligned(i) = wheelIRfirstsZero(i+1);
%                 end
%             end
%         end
%     end 
%     wheelsLfirstValid = wheelLfirstsAligned + validsCorrection;
%     wheelsRfirstValid = wheelRfirstsAligned + validsCorrection;
    
    wheelsLfirstValid = wheelLfirstsAligned;
    wheelsRfirstValid = wheelRfirstsAligned;

%     wheelsILfirstValid = wheelILfirstsAligned + validsCorrection;
%     wheelsIRfirstValid = wheelIRfirstsAligned + validsCorrection;
    
    wheelsILfirstValid = wheelILfirstsAligned;
    wheelsIRfirstValid = wheelIRfirstsAligned;

    %  wheelsRfirst = wheelRfirstsZero + validsCorrection; 
    %  checkLmv = wheelsLfirst - single(trialStartsValid(1));                             % should be nearly the same as the zeroed L movements
    %  checkRmv = wheelsRfirst - single(trialStartsValid(1));                             % should be nearly the same as the zeroed R movements
    %  ContData.behavior.wheelsLfirst = wheelsLfirst;
    %  ContData.behavior.wheelsRfirst = wheelsRfirst;
     ContData.behavior.wheelsLfirstValid = wheelsLfirstValid;
     ContData.behavior.wheelsRfirstValid = wheelsRfirstValid;
     ContData.behavior.wheelsILfirstValid = wheelsILfirstValid;
     ContData.behavior.wheelsIRfirstValid = wheelsIRfirstValid;
     ContData.behavior.wheelsAll_L = sort(horzcat(wheelsLfirstValid, wheelsIRfirstValid));
     ContData.behavior.wheelsAll_R = sort(horzcat(wheelsRfirstValid, wheelsILfirstValid));
    
     
    %sanity check - ensure that initial movements were after trial starts (i.e. diffL & diffR = rx times; should be positive values)
    for j = 1:length(LrewTrialsStarts)
        if j == length(wheelsLfirstValid)
            break
        else
            diffL(j) = wheelsLfirstValid(j) - LrewTrialsStarts(j);
        end
     end
     for j = 1:length(RrewTrialsStarts)
        if j == length(wheelsRfirstValid)
            break
        else
            diffR(j) = wheelsRfirstValid(j) - RrewTrialsStarts(j);
        end
     end 

     %% Rx times - option for removal and focusing on short/long trials: is SNr involved more/less for either?
     
     % Get a sense of "rx times" for the two directions: will be longer since these movements were defined as the first one that led to subsequent centering
%      meanDiffL = mean(diffL); minDiffL = min(diffL); maxDiffL = max(diffL); medianDiffL = median(diffL);
%      meanDiffR = mean(diffR); minDiffR = min(diffR); maxDiffR = max(diffR); medianDiffR = median(diffR);
  
     trialMat = taskbase.trialMat;        % Readout: % sanity matrix: 1.trial ind, 2.trialStarts, 3.Lfirst rewarded movement, 4.Lcentered, 5.Rfirst rewarded movement, 6.Rcentered, 7.ILmoveOn, 8.barOffL, 9.IRmoveOn, 10.barOffR, 11.trial end times 12.rx times
     rxMat = taskbase.rxMat;              % Readout: 1.trial #, 2.rxCL, 3.rxCR, 4.rxIL, 5.rxIR
     rxCondensed = taskbase.rxCondensed;  % Readout: 1.trial#, 2.rxTime, 3.CL = -1; CR = 1; IL = 0; IR = 0; Should have all info needed - above 2 are likely overkill



     
    %% This code is now more detailed & should be used in MoverScript_continuous1 calling TNC_MoveBehaeExtractContinuous1.m:
    
%      sLeverData(1,:) = sgolayfilt(leverData(1,:),9,101);                      %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
%      trialStrtsNSolData_sgolayF = sgolayfilt(trialStrtsNSol,9,101);             %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
%      trialStrtsNSolData = trialStrtsNSol;                                       %likely too raw
% 
%      ContData.behavior.trialStrtsNSolData = trialStrtsNSolData;                  %number of data points are = numel(sWheelData)/10 (and rounded up) for correction of 10kHz sampling rate
%      ContData.behavior.trialStrtsNSolData_sgolayF = trialStrtsNSolData_sgolayF;  %number of data points are = numel(sWheelData)/10 (and rounded up) for correction of 10kHz sampling rate
% 
%      startsRewNs4 = trialStrtsNSolData_sgolayF;
%      firstT = find(startsRewNs4 > 550);
     
%      vecT = [];
%      for i = 1:length(firstT)
%          if firstT(i+1) == firstT(i) + 1
%              vecT(i) = firstT(i);      %when the sliding vector breaks, a zero will appear; I will want the next value for the "start" of 1st trial start
%              if i+1 == length(firstT)
%                  break
%              end
%          end
%      end                                %I want the first value of the 2nd increasing-by-one vector
%       
%      FTinds = find(vecT== 0);
%      firstTrialInd = firstT(FTinds + 1);                                    %Boom! Note that this trial index is its timestamp
%      firstTrialIndBest = firstTrialInd - 15;                                %Will be a value further back, nearer to pulse onset 
%      
%      firstS = find(startsRewNs4 > 280 & startsRewNs4 < 400);                %first solenoid click will be > 300 && < 400 on y axis:
%      vecS = [];
%      for i = 1:length(firstS)
% %          for t = 1:length(firstTrialInd);
% %              if firstS(i) > firstTrialInd                                 %for test3, 26160 is the ind I want (just after 3rd zero); using this line too soongives too many zeros
%                  if firstS(i+1) == firstS(i) + 1
%                      vecS(i) = firstS(i);                                   %when the sliding vector breaks, a zero will appear; I will want the next value for the "start" of 1st trial start
%                      if i+1 == length(firstS)
%                          break
%                      end
%                  end
%      end
%      
%      FSinds = find(vecS == 0);
%      firstSolInd1 = FSinds(4);                                              %should be 1 ind before the solenoid
%      firstSolInd = vecS(firstSolInd1+1);                                    %Boom! 26160
%      firstSolIndBest = firstSolInd - 9;                                     %Will be a value further back, nearer to pulse onset (-9 should be ~ 100mv)
%      
%      %Would be nice to perform a double-check: use ITI & work back
%      itibuff = 100;
%      findFirstSitiHi = firstSolInd + (3000+itibuff);
%      findFirstSitiLo = firstSolInd + (3000-itibuff);

%older:
%     difLick = diff(rawLick);
%     evLick=zeros(1,numel(rawLick));
% 
%     ContData.behavior.evLick    = evLick;
%     ContData.behavior.rawLick   = rawLick;

%% BS 2.1.16 - for motor in, based on ns4 file - commented out 3.21.16 since I can use csv (taskbase file now)
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

     
     
%      %Would be nice to perform a double-check: use ITI & work back
%      itibuff = 100;
%      findFirstSitiHi = firstSolInd + (3000+itibuff);
%      findFirstSitiLo = firstSolInd + (3000-itibuff);
    

% % %% WRITE THE REACH STRUCTURE OF THE CONTDATA STRUCTURE
% %     ContData.behavior.reach = reach;
%     
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
    
