
 function [ContData, filenamestrB] = TNC_MoveBehaveExtractContinuous1(taskbase, filenamestr, targetName, dataRate, chan, elimTime, fpath)
%% This is Beth's version of Josh's extraction code from Tonic_v2 

% UPDATED LAST:
% 11.29.16 - called by moverScript_continuous1 for plotting the voltage trace
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
    
%% Load event data - use the openNEVlaser version since it's cleaner
    newPath = fpath;
    cd(newPath);
%     if exist('filenamestr', 'var');
%         [filenamestrE, path] = uigetfile('*.nev*','select the .nev file', fpath)                
%     end
%     
%     dataEvents = openNEVBS_laser(filenamestrE,'read','nosave','nomat', 'report');
%     
%     startsRewIndsTmp = find(dataEvents.Data.Spikes.Electrode == chan.rew);                    % inds of trial start and solenoids detected (ch. 129); get rid of first val, since it's a pulse detected when Orienter is started
%     startsRewTsTmp = dataEvents.Data.Spikes.Timestamps(startsRewIndsTmp(2:end))./30;                 % don't ./30 yet? will need to parse out trial start from solenoid here
%     startsRewTsRoundTmp = round(dataEvents.Data.Spikes.Timestamps(startsRewIndsTmp(2:end))./30);     % This is the value to compare     
%     diffFilteredNEV = diff(startsRewTsRoundTmp);
%     diffFilteredNEV = single(diffFilteredNEV);
%     validStrtsRews1 = find(diff(startsRewTsRoundTmp)>50);                                     % eliminate pulses detected twice
% %     validStrtsRews1 = find(diff(startsRewTsRoundTmp)>60);                                     % trial start pulse width ~ 56 ms
%         validStrtsRews = startsRewTsRoundTmp(validStrtsRews1);                                    % true nev clock from true 1st trial start
%         startsRewZeroNev = single(validStrtsRews-validStrtsRews(1));                                      % zero out the clock for comparison to taskbase
%         diffTmpNEV = single(diff(startsRewZeroNev));
% 
% %     searchStrtsRews1 = find(diff(startsRewTsRoundTmp)>1); %to compare
% %       searchStrtsRews = startsRewTsRoundTmp(searchStrtsRews1);                                   
% %       searchStartsRewZeroNev = searchStrtsRews-searchStrtsRews(1);                                      
% %       searchDiffTmpNEV = diff(searchStartsRewZeroNev);
% %     itiNev = 3000;                                                                            % to assist in troubleshooting
% %     solBufferNev = 290;                                                                       % seems to have max ~3213 & min ~3290
% %     solMinusStarts = diffTmpNEV < itiNev+solBufferNev & diffTmpNEV > 3213;                    % give time buffer - usually in the positive direction for solenoid clicks to next trial start
% %     numSol2starts = numel(find(solMinusStarts == 1));                                         % compare to csv file
%     
% %     figure; plot(diffTmpNEV(1:900));
% %     figure; plot(searchDiffTmpNEV(1:900));
% %     figure; plot(diffFilteredNEV(1:700)); 
% 
%% Load taskbase structure for behavioral data 
                      
    if exist('filenamestr', 'var');
        [filenamestrT, path] = uigetfile('*_tb.mat*','select the taskbase file', fpath)       % get the actual # correct trials (and TS) from csv file
    end

    taskbase = strcat(path, filenamestrT);
    load(taskbase);
 
    startTimesCsv = taskbase.trialStartTimes';
    startsRewsCsv = taskbase.startsRews';                                                            % trial starts and solenoid discharges from behavior file
    rewTimesCsv = taskbase.rewTimes';
    
    % Zero out csv trial starts for comparison/alignment with nev file
     validStrtsRewsCsv1 = find(diff(startsRewsCsv)>50);                                               % account for very low vals: rx time that is virtually impossible due to mouse likely already moving
     validStrtsRewsCsv = startsRewsCsv(validStrtsRewsCsv1);                                           % align this clock to nev
     startsRewsZeroCsv = validStrtsRewsCsv-validStrtsRewsCsv(1);                                      % zero out TS 1; first trial start = 0
     
     rewTimesZeroCsv = rewTimesCsv-validStrtsRewsCsv(1);    
%      
%     % Tangent to investigate the diffs bt. ea. line: (clock alignment troubleshooting)
%     diffStartsRewsCsv = diff(startsRewsZeroCsv);                                                     % refine the 3 s rule for the tighter csv file:
% %     meanDiffCsv = mean(diffStartsRewsCsv);
% %     solBufferCsv = 80;
% %     iti = 3000;
% %     solMinusStartsCSV = diffStartsRewsCsv < iti+solBufferCsv & diffStartsRewsCsv > iti-solBufferCsv; % give time buffer 
% %     numSol2startsCsv = numel(find(solMinusStartsCSV == 1));                                          % should be the num of csv rewarded trials; compare to nev
% 
% %% Parse out solenoid/trial start info: hack leftover from troubleshooting clock alignment
% % These solenoid times can shift > 300 ms up/down, depending on the NEV file
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %For the test dataset:
% % solBufferNev = 170; %NEV: add for high end = 3190
% % iti = 3020;         %NEV: low end
% % itiCsvVal = 3020;   %Csv: low end
% % solBufferCsv = 250; %Csv: add for high end = 3270
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %For 161015 dataset & 161119 test3:
% % solBufferNev = 50; %NEV: add for high end = 3190 161015 dataset
% % % solBufferNev = 200; %NEV: add for high end = 3190 161119 test3
% % iti = 3020;         %NEV: low end
% % itiCsvVal = 3020;   %Csv: low end
% % solBufferCsv = 75; %Csv: add for high end = 3270
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %for 161005 dataset & old 151105:
% solBufferNev = 270; %NEV: add for high end = 3190
% iti = 3020;         %NEV: low end
% itiCsvVal = 3020;   %Csv: low end
% solBufferCsv = 40; %Csv: add for high end = 3270
%  
% %% Align clocks to csv trial starts:
% 
%     % startTimesCsv = startTimesCsv - startTimesCsv(1);                         % includes correct and re-started trials
%     LRtrial = taskbase.LRtrial;
%     rewTimesCsvZero = rewTimesCsv - taskbase.trialStartTimes(1);
% %     LtrialStartsZero = taskbase.LtrialStarts - taskbase.trialStartTimes(1);
% %     RtrialStartsZero = taskbase.RtrialStarts - taskbase.trialStartTimes(1);
%     leftIndsCenteredTimes = taskbase.leftIndsCenteredTimes;                     % rewarded trials only
%     rightIndsCenteredTimes = taskbase.rightIndsCenteredTimes;                   % rewarded trials only
%     leftIncorrectTimes = taskbase.leftIncorrectTimes;                           % incorrect threshold was reached
%     rightIncorrectTimes = taskbase.rightIncorrectTimes;                         % incorrect threshold was reached
%    
% %% Added this section to correct for extraneous pulses into Blackrock during error time-outs:
% 
% % need to id incorrect trials in order to subtract off pulses:
% % a = find(diffTmpNEV > 1015 & diffTmpNEV < 1070);  % diff max & min of the repeating error-trial pulses
% % b = diffTmpNEV(a);
% % diffTmpNEV = [240 100 1029 1028 1026 1028 3244];
% % startsRewZero = [240 340 1360 3000 4010 5020 8010 10000];
% 
% newStartsRewZeroNev = [];
% if numel(diffTmpNEV) < numel(startsRewZeroNev)
%     diffTmpNEV(end+1) = NaN;
% end
% 
% i = 1;
% for d = 1:length(startsRewZeroNev)
%     d = i;
%     if d >= numel(startsRewZeroNev)
%         break
%     else
%         if i == 1
%             newStartsRewZeroNev(1) = startsRewZeroNev(1);
%             i = i+1;
%         else if diffTmpNEV(i)>1010 && diffTmpNEV(i)<1090 && diffTmpNEV(i+1)>1010 && diffTmpNEV(i+1)<1090 && diffTmpNEV(i+2)>1010 && diffTmpNEV(i+2)<1090 && diffTmpNEV(i+3)>1010 && diffTmpNEV(i+3)<1090
%                 newStartsRewZeroNev(i) = startsRewZeroNev(i);  % the 1 sec TS always occur in 4s or 3s
%                 newStartsRewZeroNev(i+1:i+4) = NaN;
%                 i = i+5;
%             else if diffTmpNEV(i)>1010 && diffTmpNEV(i)<1090 && diffTmpNEV(i+1)>1010 && diffTmpNEV(i+1)<1090 && diffTmpNEV(i+2)>1010 && diffTmpNEV(i+2)<1090
%                     newStartsRewZeroNev(i) = startsRewZeroNev(i);
%                     newStartsRewZeroNev(i+1:i+3) = NaN;
%                     i = i+4;
%                 else newStartsRewZeroNev(i) = startsRewZeroNev(i);
%                     i = i+1;
%                 end
%             end
%         end
%     end
% end
% numExtras = numel(find(isnan(newStartsRewZeroNev)));
% newStartsRewZeroNEV1 = isfinite(newStartsRewZeroNev);
% newStartsRewZeroNEV = newStartsRewZeroNev(newStartsRewZeroNEV1);
% 
% diffNewStartsRewZeroNEV = diff(newStartsRewZeroNEV);
% 
% %Correct for strange, short (100ms) repeating TSs in NEV file that occurred bt. first trial start (csv(1) = 0) and first solenoid click (csv(2))
% clear d;
% clear i;
% firstSolcsv = startsRewsZeroCsv(2);
% firstSolcsvBuffmin = firstSolcsv - solBufferCsv;
% firstSolcsvBuffmax = firstSolcsv + solBufferCsv;
% firstSolnev = find(newStartsRewZeroNEV >= firstSolcsvBuffmin & newStartsRewZeroNEV <= firstSolcsvBuffmax);  %usually the 9th element since I do the first trial
% firstSolNEV = newStartsRewZeroNEV(firstSolnev);  % this NEV value is the first solenoid click.
% 
% if numel(firstSolNEV > 1)
%     firstSolNEV = firstSolNEV(1);
% else if numel(firstSolNEV == 1)
%     firstSolNEV = firstSolNEV;
%     end
% end
% 
% newStartsRewZeroNEVsol = [];                     % get rid of initial extraneous pulses between 1st trial TS (0) and first solenoid click
% i = 1;
% for d = 1:length(newStartsRewZeroNEV)
%     d = i;
%     if d >= numel(newStartsRewZeroNEV)
%         break
%     else
%         if i == 1
%             newStartsRewZeroNEVsol(1) = newStartsRewZeroNEV(i);
%             i = i+1;
%         else if i < firstSolnev && i > 1
%                 newStartsRewZeroNEVsol(i) = NaN;
%                 i = i+1;
%             else newStartsRewZeroNEVsol(i) = newStartsRewZeroNEV(i);
%                 i = i+1;
%             end
%         end
%     end
% end
% newStartsRewZeroNEVnew = isfinite(newStartsRewZeroNEVsol);
% newStartsRewZeroNEV2 = newStartsRewZeroNEVsol(newStartsRewZeroNEVnew);
% diffTmpNEV2 = diff(newStartsRewZeroNEV2);
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %hack for the 161015 dataset:
% % solBufferNev = 38; %add for high end
% % iti = 3252;        %low end
% % itiCsvVal = 3020;
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %hack for the 161005 dataset:
% solBufferNev = 90; %add for high end
% iti = 3200;        %low end
% itiCsvVal = 3020;
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %hack for the 161119test dataset:
% % solBufferNev = 96; %add for high end
% % iti = 3100;        %low end
% % itiCsvVal = 3020;
% 
% itiNEV = find(diffTmpNEV2 >= iti & diffTmpNEV2 <= iti+solBufferNev);                    % add an index of 1 to these for newStrtRewZeroNEV2
% itiCSV = find(diffStartsRewsCsv >= itiCsvVal & diffStartsRewsCsv <= iti+solBufferNev);  % add an index of 1 to these for diffStrtsRewsCsv
% firstITInev = newStartsRewZeroNEV2(itiNEV(1)+1);                                        % Yes! works for 2 datasets (161015 and 161118test)
% 
% clear d; clear i;
% newStartsRewZeroNEViti = [];                                                            % get rid of initial extraneous pulses between 1st solenoid click and next trial start
% i = 1;
% for d = 1:length(newStartsRewZeroNEV2)
%     d = i;
%     if d >= numel(newStartsRewZeroNEV2)
%         break
%     else
%         if i == 1 || i == 2
%             newStartsRewZeroNEViti(i) = newStartsRewZeroNEV2(i);
%             i = i+1;
%         else if i >= 3 && i < itiNEV(1)+1
%                 newStartsRewZeroNEViti(i) = NaN;
%                 i = i+1;
%             else if i == itiNEV(1)+1  
%                 newStartsRewZeroNEViti(i) = newStartsRewZeroNEV2(i);
%                 i = i+1;
%                 else if i > itiNEV(1)+1  
%                   newStartsRewZeroNEViti(i) = newStartsRewZeroNEV2(i);
%                   i = i+1;
%                     end
%                 end
%             end
%         end
%     end
% end
% newStartsRewZeroNEVnew = isfinite(newStartsRewZeroNEViti);
% newStartsRewZeroNEV3 = newStartsRewZeroNEViti(newStartsRewZeroNEVnew);
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % end the awesome code that should pull out abberant pulses into NEV file
% %% Resume
% clear i
%     trialStrts = [];
%     for i = 1:length(startsRewsZeroCsv)
%         for j = 1:length(newStartsRewZeroNEV2)
%             if i == 1
%                 trialStrts(i) = 1;                                               % report a 1 for the first trial alignment (and since the first trial start for both clocks (by design) is defined as aligned) 
%             else
%                 if newStartsRewZeroNEV2(j) > startsRewsZeroCsv(i)-40 && newStartsRewZeroNEV2(j) < startsRewsZeroCsv(i)+40; % creates a window to search within 40 ms +/-
%                     trialStrts(j) = newStartsRewZeroNEV2(j);                         % indexed positions are true to the zeroed (to first true trial start) nev clock
%                 end                                                             
%             end                                                                  % then, after this, index to trial starts csv and find L mvments v R mvments per trial
%         end
%     end
%     trialStartsPre = find(trialStrts >= 1);
%     trialStarts = trialStrts(trialStartsPre);
% 
%     lCenteredTimesCsv = taskbase.leftCenteredTimes;
%     lCenteredZeroCsv = lCenteredTimesCsv - startTimesCsv(1);                     % zero-out l/r center times csv vector to index back into zeroed nev file
%     lIncorrectZeroCsv = leftIncorrectTimes - startTimesCsv(1);
%     
%     % alignment of solenoid clicks/reward times from csv to blackrock
%     rewTimes = [];
%     for m = 1:length(trialStarts)
%         for n = 1:length(rewTimesCsvZero)
%             if m == numel(trialStarts)
%                 break
%             end
%             if n == numel(rewTimesCsvZero)
%                 break
%             end
%             if rewTimesCsvZero(n) > trialStarts(m) && rewTimesCsvZero(n) < trialStarts(m+1) && rewTimesCsvZero(n) < rewTimesCsvZero(n+1) % just to be sure
%                 rewTimes(n) = rewTimesCsvZero(n);
%             end
%         end
%     end
%     validsCorrection = single(validStrtsRews(1));
%     rewTimesValid = rewTimes + validsCorrection;
% 
%     % sanity check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     lIndsCenteredTimesZero = [];                                                 % zero-out l center times, keeping indexed position
%     for x = 1:length(leftIndsCenteredTimes);
%         if leftIndsCenteredTimes(x) > 1
%             lIndsCenteredTimesZero(x) = leftIndsCenteredTimes(x) - startTimesCsv(1);
%         end
%         if leftIndsCenteredTimes(x) == 0
%             lIndsCenteredTimesZero(x) = 0;                                       % this vector should be the same as lCenteredZeroCsv except at correct index position for trial number
%         end
%     end    
%     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     % L rewarded trials' start times
%     lCenters = [];
%     for a = 1:length(trialStarts)
%         for b = 1:length(lCenteredZeroCsv)
%             if a == numel(trialStarts)
%                 break
%             end
%             if b == numel(lCenteredZeroCsv)
%                 break
%             end
%             if lCenteredZeroCsv(b) > trialStarts(a) && lCenteredZeroCsv(b) < trialStarts(a+1) && lCenteredZeroCsv(b) < lCenteredZeroCsv(b+1) % just to be sure
%                lCenters(a) = trialStarts(a);
%             end
%         end
%     end
%     lCenters1 = lCenters > 0; LrewTrialStartsZeroed = lCenters(lCenters1);
% %     validsCorrection = single(validStrtsRews(1));
%     trialStartsValid = trialStarts + validsCorrection;
%     LrewTrialsStarts = LrewTrialStartsZeroed + validsCorrection;
%     
%     % L incorrect trials' start times
%     lIncorrects = [];
%     for c = 1:length(trialStarts)
%         for d = 1:length(lIncorrectZeroCsv)
%             if c == numel(trialStarts)
%                 break
%              else if d == numel(lIncorrectZeroCsv)
%                     break
%                 else if lIncorrectZeroCsv(d) > trialStarts(c) && lIncorrectZeroCsv(d) < trialStarts(c+1) && lIncorrectZeroCsv(d) < lIncorrectZeroCsv(d+1)
%                         lIncorrects(c) = trialStarts(c);
%                     end
%                 end
%             end
%         end
%     end
%     lIncorrects1 = lIncorrects > 0; LincorrectsTrialStartsZeroed = lIncorrects(lIncorrects1);
%     LincorrectTrialsStarts = LincorrectsTrialStartsZeroed + validsCorrection;
%    
%     %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     % right side trials
%     rCenteredTimesCsv = taskbase.rightCenteredTimes;
%     rCenteredZeroCsv = rCenteredTimesCsv - startTimesCsv(1);
%     rIncorrectZeroCsv = rightIncorrectTimes - startTimesCsv(1);
% 
%     % r rewarded trials' start times
%     rCenters = [];
%     for m = 1:length(trialStarts)
%         for n = 1:length(rCenteredZeroCsv)
%             if m == numel(trialStarts)
%                 break
%             end
%             if n == numel(rCenteredZeroCsv)
%                 break
%             end
%             if rCenteredZeroCsv(n) > trialStarts(m) && rCenteredZeroCsv(n) < trialStarts(m+1) && rCenteredZeroCsv(n) < rCenteredZeroCsv(n+1) % just to be sure
%                rCenters(m) = trialStarts(m);
%             end
%         end
%     end
%     rCenters1 = rCenters > 0; RrewTrialStartsZeroed = rCenters(rCenters1);
%     RrewTrialsStarts = RrewTrialStartsZeroed + validsCorrection; 
%     
%     %r incorrect trials' start times
%     rIncorrects = [];
%     for o = 1:length(trialStarts)
%         for p = 1:length(rIncorrectZeroCsv)
%             if o == numel(trialStarts)
%                 break
%             end
%             if p == numel(rIncorrectZeroCsv)
%                 break
%             end
%             if rIncorrectZeroCsv(p) > trialStarts(o) && rIncorrectZeroCsv(p) < trialStarts(o+1) && rIncorrectZeroCsv(p) < rIncorrectZeroCsv(p+1)
%                 rIncorrects(o) = trialStarts(o);
%             end
%         end
%     end
%     rIncorrects1 = rIncorrects > 0; RIncorrectTrialStartsZeroed = rIncorrects(rIncorrects1);
%     RincorrectTrialsStarts = RIncorrectTrialStartsZeroed + validsCorrection;
% 
%     ContData.behavior.LrewTrialStarts = LrewTrialsStarts;
%     ContData.behavior.trialStartsValid = trialStartsValid;
%     ContData.behavior.RrewTrialStarts = RrewTrialsStarts;
%     ContData.behavior.LincorrectTrialStarts = LincorrectTrialsStarts;
%     ContData.behavior.RincorrectTrialStarts = RincorrectTrialsStarts;
%     ContData.behavior.rewTimesValid = rewTimesValid;
%     
%     rewardedTrialsNev = numel(ContData.behavior.LrewTrialStarts) + numel(ContData.behavior.RrewTrialStarts);
%     disp(['Number of rewarded trials (.nev): ' num2str(rewardedTrialsNev)]);
%     disp(' ');disp(' '); 
%     
%     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     % Load the taskbase data relevent to wheel movement extractions L v R. Could compare to ns4 file. 
%     % Note: these values still retain 0s to hold indices where L or Rtrials were either not rewarded, or were rewarded to other side.
%     
%     % Correct L and R trials:
%     wheelLfirstTime = taskbase.wheelLfirstTime;
%     wheelLfirst = taskbase.wheelLfirst;
%     wheelRfirst = taskbase.wheelRfirst;
%     wheelRfirstTime = taskbase.wheelRfirstTime;
% 
%     % Remove the empty indices and extract just the first vals needed and use LrewTrialStarts/RewTrialStarts for alignment
%     wheelLfirstT = wheelLfirstTime > 0; wheelLfirsts = wheelLfirstTime(wheelLfirstT);
%     wheelLfirstsZero = wheelLfirsts - startTimesCsv(1);
%     wheelRfirstT = wheelRfirstTime > 0; wheelRfirsts = wheelRfirstTime(wheelRfirstT);
%     wheelRfirstsZero = wheelRfirsts - startTimesCsv(1); 
%     
%     extractWheelCsvLR = numel(wheelLfirsts) + numel(wheelRfirsts);
% %     disp(['Number of L/R wheel movement extractions (csv): ' num2str(extractedWheelCsv);
% %     disp(' '); disp(' ');
%  
%     % Incorrect L and R trials:
%     wheelILfirstTime = taskbase.wheelILfirstTime;
%     wheelILfirst = taskbase.wheelILfirst;
%     wheelIRfirst = taskbase.wheelIRfirst;
%     wheelIRfirstTime = taskbase.wheelIRfirstTime;
% 
%     rightIncorrectTimes = taskbase.rightIncorrectTimes
%     wheelILfirstT = wheelILfirstTime > 0; wheelILfirsts = wheelILfirstTime(wheelILfirstT);
%     wheelILfirstsZero = wheelILfirsts - startTimesCsv(1);
%     wheelIRfirstT = wheelIRfirstTime > 0; wheelIRfirsts = wheelIRfirstTime(wheelIRfirstT);
%     wheelIRfirstsZero = wheelIRfirsts - startTimesCsv(1);
% 
%     extractWheelCsvILnIR = numel(wheelILfirsts) + numel(wheelIRfirsts);
%     
% %% aligned (to nev-based trial starts) wheel movements that were based on csv file
% %  May want to consider including this in the overall trial alignment for absolute certainty - rx times (diffL, diffR) are on super-fast side for some trials (25 ms)
%     
%     % correct L movements: eliminates first L movements that don't fit to an extracted trial start into Blackrock
%     wheelLfirstsAligned = [];
%     for i = 1:length(LrewTrialStartsZeroed)
%         if numel(LrewTrialStartsZeroed) < numel(wheelLfirstsZero) && i == numel(LrewTrialStartsZeroed)
%             break
%         else if i == numel(wheelLfirstsZero)
%                 break
%             else if wheelLfirstsZero(i) > LrewTrialStartsZeroed(i) && wheelLfirstsZero(i) < LrewTrialStartsZeroed(i+1) && wheelLfirstsZero(i) < wheelLfirstsZero(i+1) % just to be sure
%                     wheelLfirstsAligned(i) = wheelLfirstsZero(i);
%                 else if wheelLfirstsZero(i) < LrewTrialStartsZeroed(i) && wheelLfirstsZero(i+1) > LrewTrialStartsZeroed(i) && wheelLfirstsZero(i+1) < LrewTrialStartsZeroed(i+1)
%                         wheelLfirstsAligned(i) = wheelLfirstsZero(i+1);
%                     end
%                 end
%             end
%         end
%     end
%     
%     wheelILfirstsAligned = [];
%     for j = 1:length(LincorrectsTrialStartsZeroed)
%         if numel(LincorrectsTrialStartsZeroed) < numel(wheelILfirstsZero)&& j == numel(LincorrectsTrialStartsZeroed)
%             break
%         end
%         if j == numel(wheelILfirstsZero)
%             break
%         end
%         if wheelILfirstsZero(j) > LincorrectsTrialStartsZeroed(j) && wheelLfirstsZero(j) < LincorrectsTrialStartsZeroed(j+1) && wheelILfirstsZero(j) < wheelILfirstsZero(j+1)
%             wheelILfirstsAligned(j) = wheelILfirstsZero(j);
%         else if wheelILfirstsZero(j) < LincorrectsTrialStartsZeroed(j) && wheelILfirstsZero(j+1) > LincorrectsTrialStartsZeroed(j) && wheelILfirstsZero(j+1) < LincorrectsTrialStartsZeroed(j+1)
%                 wheelILfirstsAligned(j) = wheelILfirstsZero(j+1);
%             end
%         end
%     end
%       
% % Alternative for alignment of L movement vectors if above doesn't work 
%   if numel(wheelLfirstsAligned) < numel(LrewTrialStartsZeroed - 50)         % if the previous loop left out 50 of the L rew trials, do this loop
%     wheelLfirstsAligned = [];
%     for i = 1:length(LrewTrialStartsZeroed)
%         for j = 1:length(wheelLfirstsZero)
%             if i == numel(LrewTrialStartsZeroed)
%                 break
%             end
%             if j == numel(wheelLfirstsZero)
%                  else if wheelLfirstsZero(j) > LrewTrialStartsZeroed(i) && wheelLfirstsZero(j) < LrewTrialStartsZeroed(i+1) && wheelLfirstsZero(j) < wheelLfirstsZero(j+1) % just to be sure
%                         wheelLfirstsAligned(j) = wheelLfirstsZero(j);
%                 end
%             end
%          end
%       end
%   end        
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
%     % R movements: works better this way for some file's R movements. Think about this.
%     wheelRfirstsAligned = [];
%     for i = 1:length(RrewTrialStartsZeroed)
%         for j = 1:length(wheelRfirstsZero)
%             if i == numel(RrewTrialStartsZeroed)
%                 break
%             end
%             if j == numel(wheelRfirstsZero)
%                 break
%             end
%              if wheelRfirstsZero(j) > RrewTrialStartsZeroed(i) && wheelRfirstsZero(j) < RrewTrialStartsZeroed(i+1) && wheelRfirstsZero(j) < wheelRfirstsZero(j+1) % just to be sure
%                     wheelRfirstsAligned(j) = wheelRfirstsZero(j);
%              end
%         end
%     end
%     
%     wheelIRfirstsAligned = [];
%     for k = 1:length(RIncorrectTrialStartsZeroed)
%         for l = 1:length(wheelIRfirstsZero)
%             if k == numel(RIncorrectTrialStartsZeroed)
%                 break
%             end
%             if l == numel(wheelIRfirstsZero)
%                 break
%             end
%             if wheelIRfirstsZero(l) > RIncorrectTrialStartsZeroed(k) && wheelIRfirstsZero(l) < RIncorrectTrialStartsZeroed(k+1) && wheelIRfirstsZero(l) < wheelIRfirstsZero(l+1)
%                 wheelIRfirstsAligned(l) = wheelIRfirstsZero(l);
%             end
%         end
%     end
%         
% % Alternative 2 for alignment of R movement vectors if above doesn't work
% %     if numel(wheelRfirstsAligned) < numel(RrewTrialStartsZeroed/2)         % if the previous loop left out 50 of the L rew trials, do this loop
% %         wheelRfirstsAligned = [];
% %         for i = 1:length(RrewTrialStartsZeroed)
% %             if wheelRfirstsZero(i) > RrewTrialStartsZeroed(i) && wheelRfirstsZero(i) < RrewTrialStartsZeroed(i+1) && wheelRfirstsZero(i) < wheelRfirstsZero(i+1) % just to be sure
% %                 wheelRfirstsAligned(i) = wheelRfirstsZero(i);
% % 
% %             else if wheelRfirstsZero(i) < RrewTrialStartsZeroed(i) && wheelRfirstsZero(i+1) > RrewTrialStartsZeroed(i) && wheelRfirstsZero(i+1) < RrewTrialStartsZeroed(i+1)
% %                     wheelRfirstsAligned(i) = wheelRfirstsZero(i+1);
% %                 end
% %             end
% %         end
% %     end
%     
%     if numel(wheelIRfirstsAligned) < numel(RIncorrectTrialStartsZeroed/2) 
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
% 
%     wheelsLfirstValid = wheelLfirstsAligned + validsCorrection;
%     wheelsRfirstValid = wheelRfirstsAligned + validsCorrection;
%     wheelsILfirstValid = wheelILfirstsAligned + validsCorrection;
%     wheelsIRfirstValid = wheelIRfirstsAligned + validsCorrection;
%     
%     %sanity check - ensure that initial movements were after trial starts (i.e. diffL & diffR should be positive values)
%     for j = 1:length(LrewTrialsStarts)
%         if j == length(wheelsLfirstValid)
%             break
%         else
%             diffL(j) = wheelsLfirstValid(j) - LrewTrialsStarts(j);
%         end
%      end
%      for j = 1:length(RrewTrialsStarts)
%         if j == length(wheelsRfirstValid)
%             break
%         else
%             diffR(j) = wheelsRfirstValid(j) - RrewTrialsStarts(j);
%         end
%      end 
%      % Get a sense of "rx times" for the two directions: will be longer since these movements were defined as the first one that led to subsequent centering
% %      meanDiffL = mean(diffL); minDiffL = min(diffL); maxDiffL = max(diffL); medianDiffL = median(diffL);
% %      meanDiffR = mean(diffR); minDiffR = min(diffR); maxDiffR = max(diffR); medianDiffR = median(diffR);
%      
%     %  wheelsRfirst = wheelRfirstsZero + validsCorrection; 
%     %  checkLmv = wheelsLfirst - single(trialStartsValid(1));                             % should be nearly the same as the zeroed L movements
%     %  checkRmv = wheelsRfirst - single(trialStartsValid(1));                             % should be nearly the same as the zeroed R movements
%     %  ContData.behavior.wheelsLfirst = wheelsLfirst;
%     %  ContData.behavior.wheelsRfirst = wheelsRfirst;
%          
%      ContData.behavior.wheelsLfirstValid = wheelsLfirstValid;
%      ContData.behavior.wheelsRfirstValid = wheelsRfirstValid;
%      ContData.behavior.wheelsILfirstValid = wheelsILfirstValid;
%      ContData.behavior.wheelsIRfirstValid = wheelsIRfirstValid;
    
%% LOAD CONTINUOUS DATA
    
 switch dataRate
    
    case 'ns4'    
        [filenamestrE, path] = uigetfile('*.ns4','select the .ns4 file', fpath); 
        Ns4DATA = openNSxBS('report','read',filenamestrE);

        % How Josh's variables go down        
%         xChan = find(Ns4DATA.MetaTags.ChannelID==chan.x);                 %my x channel is empty  
%         yChan = find(Ns4DATA.MetaTags.ChannelID==chan.y);                 %motor in; now is lick
%         lChan = find(Ns4DATA.MetaTags.ChannelID==chan.lick);
        yChan = find(Ns4DATA.MetaTags.ChannelID==chan.y);                   %lick channel
        trialStrtsChan = find(Ns4DATA.MetaTags.ChannelID==chan.rew);        %trial start/sol channel

%         leverData(1,:) = decimate(Ns4DATA.Data(xChan,:),10);
%         leverData(2,:) = decimate(Ns4DATA.Data(yChan,:),10);              %correct for 10kS/s rate
%         rawLick = decimate(Ns4DATA.Data(lChan,:),10);
%         lick = decimate(Ns4DATA.Data(yChan,:),10);                        %get this enabled, Beth!
        trialStrtsNSol = decimate(Ns4DATA.Data(trialStrtsChan,:),10);       %convert these to ms (/10)

%     case 'ns2'
%         filenamestrE = uigetfile('*.ns2','select the .ns2 file', dataRate);
%         Ns2DATA = openNSx('report','read',filenamestrE);
%         clear leverData sLeverData tmpLeverData
% 
%         yChan = find(Ns2DATA.MetaTags.ChannelID==chan.y);
%         lick = decimate(Ns4DATA.Data(yChan,:),10);
        
 end
%     
%      sLeverData(1,:) = sgolayfilt(leverData(1,:),9,101);                      %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
     trialStrtsNSolData_sgolayF = sgolayfilt(trialStrtsNSol,9,101);             %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
     trialStrtsNSolData = trialStrtsNSol;                                       %likely too raw

     ContData.behavior.trialStrtsNSolData = trialStrtsNSolData;                  %number of data points are = numel(sWheelData)/10 (and rounded up) for correction of 10kHz sampling rate
     ContData.behavior.trialStrtsNSolData_sgolayF = trialStrtsNSolData_sgolayF;  %number of data points are = numel(sWheelData)/10 (and rounded up) for correction of 10kHz sampling rate

%      startsRewNs4 = trialStrtsNSolData_sgolayF;
%      firstT = find(startsRewNs4 > 550);
%      
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
%     
% 








%older:
% %     difLick = diff(rawLick);
% %     evLick=zeros(1,numel(rawLick));
% % 
% %     ContData.behavior.evLick    = evLick;
% %     ContData.behavior.rawLick   = rawLick;
% 
% %% BS 2.1.16 - based on ns4 file - commented out 3.21.16 since I can use csv (taskbase file now)
% % 
% % % Find any leftward movement (-going voltages) or rightward movement (+)
% %   % Find "center" of raw wheel voltages
% %     centers = median(sWheelData); 
% %     wheelMin = min(sWheelData);    %just to see
% %     wheelMax = max(sWheelData);    %just to see
% % 
% %     wheelData = single(sWheelData); 
% %     wheelData = sWheelData(2,:);
% %     wheelDataInds = [];
% %     for w = 1:length(wheelData);
% %         wheelDataInds(:,w) = w;
% %     end
% % 
% %     newWheelMat(1,:) = wheelDataInds;
% %     newWheelMat(2,:) = wheelData;
% % 
% %     % L movement extractions:
% %     allL_wheel = [];
% %     allL_wheelInds = [];
% %     for z = 1:length(wheelData)
% %         if wheelData(z) > wheelData(z+1)
% %             if z+1 == length(wheelData)
% %                 break;
% %             else
% %                 allL_wheel(z) = wheelData(z);
% %                 allL_wheelInds(z) = wheelDataInds(z);
% %             end
% %         end
% %     end
% % 
% %     % Eliminate wheel recordings that occurred before first trial start (.nev)
% % %     validWheel_L = find(allL_wheelInds' >= trialStartsNew(1));            % These are recording times(ms) when wheel is moving left; pull out first vals...
% %     validWheel_L = find(allL_wheelInds' >= ContData.behavior.trialStartsValid(1));
% %     
% %     validsL = [];
% %     validsL_first = [];
% %     for v = 1:length(validWheel_L)
% %         if v == numel(validWheel_L)
% %             break;
% %         else
% %             validsL(v+1) = validWheel_L(v+1) > validWheel_L(v) + 5;         % Cut-off is 5ms;
% %         end
% %     end
% % 
% %     validsL_first = validWheel_L(validsL == 1);                             % Here are the first TSs for leftward vectors
% %     
% %     %Write the vectors for all L wheel movements to the ContData structure
% % %     ContData.behavior.validWheel_L = validWheel_L;
% %     ContData.behavior.newWheelMat = newWheelMat;
% % %     ContData.behavior.allL_wheel = single(allL_wheel);                    % Great, but these aren't necessary
% % %     ContData.behavior.allL_wheelInds = single(allL_wheelInds);
% %     
% %     allL_wheelIndsShort = find(allL_wheelInds>0);
% % %     ContData.behavior.allL_wheelIndsShort = single(allL_wheelIndsShort); 
% % %     validWheel_L = find(allL_wheelIndsShort' >= trialStartsNew(1));       %These are recording times(ms) when wheel is moving left (without the zeros)
% %     ContData.behavior.validsL_first = validsL_first;
% % 
% %     %R movement extractions:
% %     allR_wheel = [];
% %     allR_wheelInds = [];
% %     zrlim = numel(wheelData(1:end-10))
% % 
% %     for zr = 1:length(wheelData)
% %         if wheelData(zr+1) > wheelData(zr)
% %             if zr == zrlim;
% %                 break;
% %             else
% %                 
% %                 allR_wheel(zr) = wheelData(zr);
% %                 allR_wheelInds(zr) = wheelDataInds(zr);
% %                 
% %             end
% %         end
% %     end
% % 
% %     % Eliminate wheel recordings that occurred before first trial start (.nev)
% %     validWheel_R = find(allR_wheelInds' >= ContData.behavior.trialStartsValid(1));  %These are recording times(ms) when wheel is moving left; pull out first vals...
% %     
% %     validsR = [];
% %     validsR_first = [];
% %     for vr = 1:length(validWheel_R)
% %         if vr == numel(validWheel_R)
% %             break;
% %         else
% %             validsR(vr+1) = validWheel_R(vr+1) > validWheel_R(vr) + 5;              %cut-off is 5ms
% %         end
% %     end
% % 
% %     validsR_first = validWheel_R(validsR == 1);                                     %Here are the first TSs for rightward vectors
% %     ContData.behavior.validsR_first = validsR_first;
%     
%  %% EXTRACT VELOCITY DATA FROM CONT LEVER DATA. BS: motor in: sLeverV = unfiltered velocities; sLeverVm = filtered
% % 
% %     numSamples = size(ContData.behavior.sLeverData,2);
% %     tmpLeverData(1,:) = sgolayfilt(ContData.behavior.sLeverData(1,:),3,151); 
% %     tmpLeverData(2,:) = sgolayfilt(ContData.behavior.sLeverData(2,:),3,151);
% % 
% %     sLeverV = zeros(1,numSamples);
% % 
% %     disp(' ');disp(' ');disp('Extracting velocity...');
% % 
% %     dX = diff(tmpLeverData(1,:));
% %     dY = diff(tmpLeverData(2,:));
% % 
% %     sLeverV = sqrt( dX.^2 + dY.^2 );
% % 
% %     disp(' ');disp(' Complete. ');disp(' ');
% % 
% %     ContData.behavior.sLeverV = sLeverV;
% %     ContData.behavior.sLeverVm = sgolayfilt(ContData.behavior.sLeverV,3,501); %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
% %     
% % %     figure(1); plot(sLeverV);
% %     
% %     clear sLeverV;
% % 
% % %% FIND MOVEMENTS
% %     method = 'vel';
% %     numC = 5;
% %     clear progSt* reach
% % 
% %     currX       = ContData.behavior.sLeverData(1,:);
% %     currY       = ContData.behavior.sLeverData(2,:);
% %     currV       = ContData.behavior.sLeverVm; 
% % 
% %     pre     = 10;
% %     post    = 10;       
% %     minSpace = 250;
% %     count = 1;
% % 
% %     % threshold the velocities
% %     switch dataRate
% %         case 'ns4'                                                            %currV are the velocity values
% %             allValidSamps   = find(currV(1:numel(currV)-(minSpace.*2))>0.35); %find point where filtered velocities begin at -500 and >.35 (will be that index and above)
% %         case 'ns2'
% %             allValidSamps   = find(currV(1:numel(currV)-(minSpace.*2))>1.5);
% %     end
% %     
% %     if numel(allValidSamps)>0
% % 
% %         switch method
% % 
% %             case 'vel'
% % 
% %     %             progStartTMP(count,1)  = currStamps(allValidSamps(1));
% %                 progStartTMP(count,2)  = currX(allValidSamps(1));
% %                 progStartTMP(count,3)  = currY(allValidSamps(1));           %find currY positions for velocities >.35
% %                 progStartTMP(count,4)  = allValidSamps(1);
% % 
% %                 for j=2:numel(allValidSamps)                                %allValidSamps are the indexed velocities >.35
% % 
% %                     if allValidSamps(j)>allValidSamps(j-1)+minSpace         %if the indexed velocities are increasing by at least 250
% % 
% %                         postN = find(currV(allValidSamps(j-1):allValidSamps(j))<1,1,'first');  %postN = find the first velocity values at those positions
% %                         if numel(postN)<1
% %                             figure(3); clf;                                                         %if postN has less than 1 value, plot                
% %                             plot(progStartTMP(count,4):allValidSamps(j),currV(progStartTMP(count,4):allValidSamps(j)),'k');  
% %                             hold on; plot(progStartTMP(count,4),currV(progStartTMP(count,4)),'bo',allValidSamps(j-1),currV(allValidSamps(j-1)),'ro');
% %                             pause(0.1);
% % 
% %                             disp('Cannot find stop');                                               %and set postN = post (which was defined as 10)
% %                             postN=post;
% %                         end
% %     %                     progStopTMP(count,1)   = currStamps(allValidSamps(j-1)+post);
% %                         progStopTMP(count,2)   = currX(allValidSamps(j-1)+postN);                %But if postN has >1 value, create a matrix (progStopTMP) and find currY at index allValidSamps at the previous positions and add 10 
% %                         progStopTMP(count,3)   = currY(allValidSamps(j-1)+postN);
% %                         progStopTMP(count,4)   = allValidSamps(j-1)+postN;
% %                         count                  = count+1;
% % 
% %                         preN = find(currV(allValidSamps(j-1) : allValidSamps(j))<1,1,'last');
% %                         if numel(preN)<1
% %                             disp('Cannot find start');
% %     %                     progStartTMP(count,1)  =
% %     %                     currStamps(allValidSamps(j)-pre)                                        %Now create the startTMP matrix
% %                             progStartTMP(count,2)  = currX(allValidSamps(j)-pre);
% %                             progStartTMP(count,3)  = currY(allValidSamps(j)-pre);
% %                             progStartTMP(count,4)  = allValidSamps(j)-pre;                    
% %                         else
% %     %                     progStartTMP(count,1)  = currStamps(allValidSamps(j-1)+preN);
% %                             progStartTMP(count,2)  = currX(allValidSamps(j-1)+preN);
% %                             progStartTMP(count,3)  = currY(allValidSamps(j-1)+preN);
% %                             progStartTMP(count,4)  = allValidSamps(j-1)+preN;
% %                         end
% %                     end
% % 
% %                     if j==numel(allValidSamps)
% %     %                     post = find(currV(allValidSamps(j):allValidSamps(j)+minSpace)<0.5,1,'first');
% %     %                     progStopTMP(count,1)   = currStamps(allValidSamps(j)+post);
% %                         progStopTMP(count,2)   = currX(allValidSamps(j)+post);
% %                         progStopTMP(count,3)   = currY(allValidSamps(j)+post);
% %                         progStopTMP(count,4)   = allValidSamps(j)+post;
% %                     end
% % 
% %                 end
% % 
% %                 count = 1;
% %                 for k = 1:size(progStartTMP,1)
% % 
% %                     if k==1
% %                         reach.init = 1;
% %                     end
% % 
% %                     % reaches must be at least 50 ms long
% %                     if progStopTMP(k,4)-progStartTMP(k,4)>=90 & progStartTMP(k,4)>minSpace
% % 
% %                         trajAngle   = atan2(progStopTMP(k,3)-progStartTMP(k,3),progStopTMP(k,2)-progStartTMP(k,2));
% % 
% %                         if (pdist2([progStopTMP(k,2),progStopTMP(k,3)],[mean(currX),mean(currY)]) > pdist2([progStartTMP(k,2),progStartTMP(k,3)],[mean(currX),mean(currY)]))
% %                             reach.out(count) = 1;
% %                         else
% %                             reach.out(count) = 0;
% %                         end
% %                         velTraj = ContData.behavior.sLeverV(progStartTMP(k,4) : progStopTMP(k,4));
% %                         xVals = ContData.behavior.sLeverData(1,progStartTMP(k,4) : progStopTMP(k,4));
% %                         yVals = ContData.behavior.sLeverData(2,progStartTMP(k,4) : progStopTMP(k,4));
% % 
% %                         reach.start(count,:)  = progStartTMP(k,:);
% %                         reach.stop(count,:)   = progStopTMP(k,:);
% %                         reach.angle(count,1)  = trajAngle;
% %                         reach.dist(count,1)   = trapz(velTraj);
% %                         reach.dist(count,2)   = pdist2(progStartTMP(k,2:3) , progStopTMP(k,2:3));
% % 
% %                         tmp = findpeaks(velTraj);
% %                         reach.numpks(count,1) = numel(tmp.loc);
% %                         reach.dur(count,1)    = progStopTMP(k,4) - progStartTMP(k,4);
% %                         reach.vel(count,1)   = max(velTraj);
% %                         reach.vel(count,2)   = trapz(velTraj) ./ reach.dur(count,1);
% %                         reach.vel(count,3)   = var(velTraj);
% %                         reach.vel(count,4)   = find(velTraj==max(velTraj),1);
% % 
% %                         reach.acc(count,1)   = max(diff(velTraj));
% %                         reach.acc(count,2)   = mean(diff(velTraj));
% %                         reach.acc(count,3)   = max(diff(velTraj(1:90))); % max in first 90 ms of movement
% % 
% %                         reach.tort(count,1)  = reach.dist(count,1) ./ pdist2([progStopTMP(k,2),progStopTMP(k,3)],[progStartTMP(k,2),progStartTMP(k,3)]);
% % 
% %                         % find max displacement of the reach
% %                         xVals = xVals - xVals(1);
% %                         yVals = yVals - yVals(1);
% %                         
% %                         
% %                         
% %                         valThrInd = find(ContData.behavior.threshInds>reach.start(count,4) & ContData.behavior.threshInds<reach.stop(count,4));
% %                         if numel(valThrInd)>0
% %                             reach.rewarded(count)= 1;
% %                         else
% %                             reach.rewarded(count)= 0;                       
% %                         end
% % 
% %                         count                 = count+1;
% % 
% %                     end            
% %                 end
% % 
% % 
% %         end
% % 
% %         disp(['Num valid reaches: ' num2str(count-1)]);
% % 
% %         level = ones(1,size(progStartTMP,1)).*0.35;
% %         figure(3); clf;
% %         plot(ContData.behavior.sLeverV,'k'); title(['sLeverV']);hold on;
% %         plot(progStartTMP(:,4),level,'r^'); title(['reach starts']);
% %         plot(progStopTMP(:,4),level,'bo');  title(['reach stops']);
% % 
% %         [vals,inds] = sort(reach.dur,'ascend');    
% % 
% %         reach.numReaches = size(reach.start,1);
% % 
% %         dims = floor(sqrt(reach.numReaches));
% %         scaler = 3000;
% % 
% %         if dispOn
% %             
% %             figure(5); clf;
% %         
% %             for p = 1:1:dims.^2
% % 
% %                 l = inds(p);
% % 
% %                 offset = l.*ones(1,numel(reach.start(l,4):reach.stop(l,4)));
% %                 [x,y] = ind2sub(dims,p);
% % 
% %                 figure(5);
% %                 plot((x.*scaler),(y.*scaler),'k.','MarkerSize',8);
% %                 hold on;
% %                 plot(currX(reach.start(l,4):reach.stop(l,4))-currX(reach.start(l,4))+(x.*scaler),currY(reach.start(l,4):reach.stop(l,4))-currY(reach.start(l,4))+(y.*scaler) , 'k','LineWidth',1,'Color',[p/reach.numReaches 0.67-(p/reach.numReaches)*0.67 1-(p/reach.numReaches)]);
% %                 if i==dims.^2
% %                     axis tight; axis off;
% %                 end
% % 
% %                 startInd = round(reach.start(l,4)./2);
% %                 stopInd = round(reach.stop(l,4)./2);
% % 
% %             end
% %         end
% %         
% %         if reach.numReaches>5
% %                 figure(7); clf;
% %                 subplot(311);
% %                 p = polyfit(reach.dist(:,1),reach.vel(:,1),1);
% %                 plot(reach.dist,reach.vel(:,1),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
% %                 xlabel('Total distance');
% %                 ylabel('Peak Velocity');
% % 
% %                 subplot(312);
% %                 p = polyfit(reach.dist(:,1),reach.dur,2);
% %                 plot(reach.dist,reach.dur,'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
% %                 xlabel('Total distance');
% %                 ylabel('Movement Duration (ms)');
% % 
% %                 subplot(313);
% %                 p = polyfit(reach.dist(:,1),reach.acc(:,3),1);
% %                 plot(reach.dist,reach.acc(:,3),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
% %                 xlabel('Total distance');
% %                 ylabel('Initial acceleration');
% % 
% %                 figure(8); clf;
% %                 rose(reach.angle);
% % 
% %                 drawnow;
% %         end
% %         
% %     else
% %         
% %         reach.numReaches=0;
% %     
% %     end
% %     
% % 
% %     
% % %% DIMENSION REDUCTION OF REACH PARAMETERS
% % 
% % if reach.numReaches > 5
% %     clear matFor*
% % 
% %     matForDR = [reach.vel reach.acc reach.tort reach.dur reach.dist];
% % 
% %     for j=1:size(matForDR,2)
% %         matForCov(:,j) = (matForDR(:,j) - mean(matForDR(:,j))) ./ std(matForDR(:,j)) ;
% %     end
% % 
% %     figure(11);
% %     imagesc(corr(matForCov),[0 1]);
% %     map = colormap(TNC_CreateRBColormap(1024,'cpb'));
% % %     colormap(1-map);
% % 
% %     [eV,eD] = eig(cov(matForCov));
% %     figure(12); plot(cumsum(diag(eD))./sum(diag(eD)),'ko-');
% % 
% %     for m=1:size(matForDR,1)
% %         reach.pca(1,m) = dot(matForCov(m,:),eV(:,size(matForCov,2))');
% %         reach.pca(2,m) = dot(matForCov(m,:),eV(:,size(matForCov,2)-1)');
% %         reach.pca(3,m) = dot(matForCov(m,:),eV(:,size(matForCov,2)-2)');
% %     end
% % end
% % 
% % %% CALCULATE META-REACH PARAMETERS
% %     if reach.numReaches > 5
% %         [winReach] = TNC_ReachVigorWindow(reach.pca(1,:),reach.numReaches,9); % capture fluctuations along the maximally variant dimension
% %         ContData.behavior.winReach = winReach;
% %     end
% %     
% % %% WAS THE REACH REWARDED?
% % if reach.numReaches > 5
% %     for p=1:reach.numReaches
% % 
% %         tmp = find(ContData.behavior.threshInds>reach.start(p,4) & ContData.behavior.threshInds<reach.stop(p,4));
% %         if numel(tmp)>0
% %             reach.rewarded(p)=1;
% %         else
% %             reach.rewarded(p)=0;
% %         end
% %         
% %     end
% % end    
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
    
