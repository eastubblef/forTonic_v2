
 function [ContData, filenamestrT] = TNC_MoverBehaviorExtractBSwheel_2(taskbase, filenamestr, targetName, dataRate, chan, elimTime, fpath)
%% This is EAS's version of Josh's from Tonic_v2 

% UPDATED
% 2.8.17 to try and extract wheel data from blackrock rather than .csv file - use MoverScript_continuous to call
% 3.2.16 to incorporate behavLoad.m for loading behavioral data from Beth's (updated) "taskbase" structure

%  Note: removed the following vars since that is empty in my Blackrock system (7.16.15): NEVfiles
%  chan.x, chan.rew, 
%  ContData.behavior.threshInds = TS for laser pulse-on (divide by 30, 30kHz = sampling frequency for spikes)

% RE-ESTABLISH ACTIVE INPUTS FOR 11.26.15 DATA (1.14.16 & 1.28.16): 
% motor in = chan.y (131) captured on ns4 file, sampling rate was 10kHz; divide by 10)

% Also saves a bh.mat file that could be loaded into the next mfile for TS from recordings: TNC_ConvertTSDtoPopDataBS
% chan.rew = 129;                         % solenoid input & trial starts:              Ainp = 1
% chan.thr = 130;                         % for laser - used for tagging:               Ainp = 2
% chan.y = 131;                           % wheel (1.2.17); was for lick(post 5.1.16):  Ainp = 3
% chan.x = [];                

%% FILE NAMES
    
    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Beginning file: ' filenamestr(1,1:length(filenamestr)-3)]);
    disp(' ');disp(' ');  

    dispOn = 0;
    
%% Load event data: 
%     newPath = fpath;
%     cd(newPath);
%     if exist('filenamestr', 'var');
%         [filenamestrE, path] = uigetfile('*.nev*','select the .nev file', fpath)                
%     end    
%     dataEvents = openNEVBS_laser(filenamestrE,'read','nosave','nomat', 'report');              
%     rewIndsTmp = find(dataEvents.Data.Spikes.Electrode == chan.rew);                            %TS of trial start and solenoids detected (ch. 129)
%     rewIndsTmpZero = rewIndsTmp-rewIndsTmp(1);
%     deltaRewIndsTmp = diff(rewIndsTmp);
    
%% Load taskbase structure for behavioral data: updated 3.2.16 - use w/ .nev data for TS of trial starts alignment
                     
    if exist('filenamestr', 'var');
        [filenamestrT, path] = uigetfile('*_tb.mat*','select the taskbase file', fpath)           %get the actual # correct trials (and TS) from csv file
    end

    taskbase = strcat(path, filenamestrT);
    load(taskbase);

    startTimesCsv = taskbase.trialStartTimes';
    startsRewsCsv = taskbase.startsRews';                                  %trial starts and solenoid discharges from behavior file
    rewTimesCsv = taskbase.rewTimes';
    
    % Zero out csv trial starts for comparison/alignment with ns4 file
     validStrtsRewsCsv1 = find(diff(startsRewsCsv)>50);                    %account for very low vals: rx time that is virtually impossible due to mouse likely already moving
     validStrtsRewsCsv = startsRewsCsv(validStrtsRewsCsv1);                %align this clock to nev
     startsRewsZeroCsv = validStrtsRewsCsv-validStrtsRewsCsv(1);           %zero out TS 1; first trial start = 0
     
     rewTimesCsvZero = rewTimesCsv-validStrtsRewsCsv(1);  
%      rewTimesZeroCsv = rewTimesCsv-validStrtsRewsCsv(1);  
     startTimesZeroCsv = startTimesCsv - startTimesCsv(1);
    
%      startTimesZeroCsv = startTimesCsv - startTimesCsv(2);               %for particularly hard-to-align sessions (due to restarts), check this too:
%      startTimesZeroCsv = startTimesZeroCsv(2:end);                       %must get zero as first element

%% LOAD CONTINUOUS WHEEL/TRIAL DATA
    
switch dataRate
    
    case 'ns4'    
        [filenamestrE, path] = uigetfile('*.ns4','select the .ns4 file', fpath); 
 
        Ns4DATA = openNSxBS('report','read',filenamestrE);
%          xChan = find(Ns4DATA.MetaTags.ChannelID==chan.x);               %my x channel is empty  
        yChan = find(Ns4DATA.MetaTags.ChannelID==chan.y);                  %Motor in
        trialStrtsChan = find(Ns4DATA.MetaTags.ChannelID==chan.rew);       %trial start/sol channel

%          leverData(1,:) = decimate(Ns4DATA.Data(xChan,:),10);
        leverData(1,:) = decimate(Ns4DATA.Data(yChan,:),10);               %correct for 10kS/s rate
%         rawLick = decimate(Ns4DATA.Data(lChan,:),10);
        trialStrtsNSol = decimate(Ns4DATA.Data(trialStrtsChan,:),10);      %same; convert these to ms (/10)

%     case 'ns2'
%         filenamestrE = uigetfile('*.ns2','select the .ns4 file', dataRate);
%         Ns2DATA = openNSx('report','read',filenamestrE);
%         clear leverData sLeverData tmpLeverData
%         xChan = find(Ns2DATA.MetaTags.ChannelID==chan.x);                   
%         yChan = find(Ns2DATA.MetaTags.ChannelID==chan.y);
%         lChan = find(Ns2DATA.MetaTags.ChannelID==chan.lick);
%         leverData(1,:)  = Ns2DATA.Data(xChan,:);                       
%         leverData(2,:)  = Ns2DATA.Data(yChan,:);
%         rawLick         = Ns2DATA.Data(lChan,:);       
end

     trialStrtsNSolData_sgolayF = sgolayfilt(trialStrtsNSol,9,101);        %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
     trialStrtsNSolData = trialStrtsNSol;                                  %unfiltered
     sLeverData = sgolayfilt(leverData(1,:),9,101);                        %motor %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
     wheelData = sLeverData;
     
%     wheelFilt = single(wheelData);                                       %plot, just to see 
%     filtMin = 100;     
%     for w = 1:length(wheelFilt)/filtMin
%         newW(w) = filtMin*w -99;
%         wheelFiltMin = wheelFilt(newW);
%     end
%     figure; hold on;
%     plot(wheelFiltMin(1:1000));

     ContData.behavior.trialStrtsNSolDataFilt = trialStrtsNSolData_sgolayF; 
     ContData.behavior.trialStrtsNSolData = trialStrtsNSolData; 
     ContData.behavior.wheelData = wheelData;                              %number of data points are = numel(sWheelData)/10 (and rounded up) for correction of 10kHz sampling rate
     ContData.behavior.wheelDataRaw = leverData;
%     
%      figure; hold on; %plot(trialStrtsNSolData(1:150000));
%        plot(trialStrtsNSolData_sgolayF(80000:130000), 'r');
% %      plot(wheelData(1:200000));                                            %overlay it with trial strts/sol signal
%        plot(wheelData(80000:130000));     
% 
%  %% Find trial starts:
%      
%      startsRewNs4 = trialStrtsNSolData;
%      yTholdMin = 200;  yTholdMax = 900;   %170112 session
%      xThold = 100;  %xTholdMax = 115      %170112 session   
% 
%      firstT = find(startsRewNs4 > yTholdMin & startsRewNs4 < yTholdMax);   %find the tall gaps = ea. trial start; need a Thold > any sol signal
%      vecT = [];
%      for i = 1:length(firstT)
%          if i+1 == length(firstT)
%              break
%          end
%          if firstT(i+1) == firstT(i) + 1
%              vecT(i) = firstT(i);                                          %when the sliding vector breaks, a zero will appear;
%          end
%      end                                                                   %Find the span of >0 vals = 113 (ms)
%      Tvec = [];
%      for i = 1:length(vecT)
%          if i+xThold == length(vecT)
%              break
%          end
%          if vecT(i+xThold) == vecT(i) + xThold;
%              Tvec(i:i+xThold) = vecT(i):vecT(i) + xThold;
%          end
%      end                                                                   %first trial should be the 1st element of second vector in which 100ish vals > 0;
%      FTinds = Tvec > 0;
%      newVecT = [];
%      for i = 1:length(FTinds)
%          if FTinds(i) == 0 && FTinds(i+1) == 1
%           newVecT(i) = FTinds(i+1);
%          end
%      end  
%      Tinds = find(newVecT == 1);                                           %these should be trial start times from 2nd element on.
% %      firstTrialInd = Tinds(2);
% %      firstTrialStartNS4 = vecT(firstTrialInd + 1);
%      allTrialsNS4 = vecT(Tinds+1);
%      
%      %Here comes the variability among datasets; figure out more elegant way to do this
% %      allTrialStartsNS4 = allTrialsNS4(2:end);                            %Boom! Note that this trial index is its timestamp
%      allTrialStartsNS4 = allTrialsNS4(1:end);                              %Boom! Note that this trial index is its timestamp
%      diffTrialStartsNS4 = diff(allTrialStartsNS4);
% %      trialStartsZeroNS4 = allTrialStartsNS4 - allTrialStartsNS4(1);      %for 161015 (and most sessions)
% %      trialStartsZeroNS4 = allTrialStartsNS4 - allTrialStartsNS4(4);      %for 161005 and 160505
% %      trialStartsZeroNS4 = trialStartsZeroNS4(4:end);
%      trialStartsZeroNS4 = allTrialStartsNS4 - allTrialStartsNS4(2);        %for 161005 and 160505
%      trialStartsZeroNS4 = trialStartsZeroNS4(2:end);
% 
%      trueTrialStartsNS4 = allTrialsNS4(2:end);                             %Boom! (Use l. 158 & 164 to determine which val to start with THIS is what I want to align wheel to. 
% 
% %% confirm no extraneous trials/left out trials; alignment of trial starts (.ns4) to .csv clock:     
% %     clear i
% %     trialStrts = [];
% %     for i = 1:length(startTimesZeroCsv)
% %         for j = 1:length(trialStartsZeroNS4)
% %             if i == 1
% %                 trialStrts(i) = 1;                                               % report a 1 for the first trial alignment (and since the first trial start for both clocks (by design) is defined as aligned) 
% %             else
% %                 if trialStartsZeroNS4(j) > startTimesZeroCsv(i)-80 && trialStartsZeroNS4(j) < startTimesZeroCsv(i)+80; % creates a window to search within 40 ms +/-
% %                     trialStrts(j) = trialStartsZeroNS4(j);                         % indexed positions are true to the zeroed (to first true trial start) nev clock
% %                 end                                                             
% %             end                                                                  % then, after this, index to trial starts csv and find L mvments v R mvments per trial
% %         end
% %     end
% %     trialStartsPre = find(trialStrts >= 1);
% %     trialStarts = trialStrts(trialStartsPre);
% %     
% %     ContData.behavior.trialStarts = trialStarts;
% %     
% %     leftIncorrectTimes = taskbase.leftIncorrectTimes;
% %     leftIndsCenteredTimes = taskbase.leftIndsCenteredTimes;
% % 
% %     lCenteredTimesCsv = taskbase.leftCenteredTimes;
% %     lCenteredZeroCsv = lCenteredTimesCsv - startTimesCsv(1);                     % zero-out l/r center times csv vector to index back into zeroed nev file
% %     lIncorrectZeroCsv = leftIncorrectTimes - startTimesCsv(1);
% %     
% %     % alignment of solenoid clicks/reward times from csv to blackrock
% %     rewTimes = [];
% %     for m = 1:length(trialStarts)      %ns4
% %         for n = 1:length(rewTimesCsvZero)
% %             if m == numel(trialStarts) %ns4
% %                 break
% %             end
% %             if n == numel(rewTimesCsvZero)
% %                 break
% %             end
% %             if rewTimesCsvZero(n) > trialStarts(m) && rewTimesCsvZero(n) < trialStarts(m+1) && rewTimesCsvZero(n) < rewTimesCsvZero(n+1) 
% %                 rewTimes(n) = rewTimesCsvZero(n);
% %             end
% %         end
% %     end
% %    rewTimesValid = rewTimes;            %ns4
% %    ContData.behavior.rewTimes = rewTimes;                                   %up to l. 189 from TNC_MoveBehaveExtractContinuous.m
%         
% %% Find any leftward movement vectors (-going voltages) 
%     
%     wheelMin = min(wheelData);    %just to see
%     wheelMax = max(wheelData);    %just to see
%     wheelDataInds = 1:numel(wheelData);
% 
%     % L movement extractions:
%     allL_wheel = [];
%     allL_wheelInds = [];
%     for z = 1:length(wheelData)
%             if z+1 == length(wheelData)
%                 break;
%             else
%             if wheelData(z) > wheelData(z+1)
%                 allL_wheel(z) = wheelData(z);
%                 allL_wheelInds(z) = wheelDataInds(z);
%             end
%         end
%     end
%     
%     %Eliminate wheel movements that occurred before 1st trial start (.ns4)
%     validWheel_L = find(allL_wheelInds >= trueTrialStartsNS4(1));          %these are recording times(ms) when wheel is moving left; pull out first vals...
%     
%     validsL = [];
%     validsL_first = [];
%     for v = 1:length(validWheel_L)
%         if v == numel(validWheel_L)
%             break;
%         else
%             validsL(v+1) = validWheel_L(v+1) > validWheel_L(v) + 20;       %cut-off is 70ms
%         end
%     end
% 
%     validsL_first = validWheel_L(validsL == 1);                            %here are the first TSs for leftward vectors
%     ContData.behavior.wheel_L_first = validsL_first;
%         
% %     allL_wheelIndsShort = find(allL_wheelInds>0);
% %     ContData.behavior.allL_wheelIndsShort = single(allL_wheelIndsShort); 
% %     validWheel_L = find(allL_wheelIndsShort' >= trialStartsNew(1));      %these are recording times(ms) when wheel is moving left (without the zeros)
% 
% %% Find any rightward movement vectors (+going voltages) 
% 
%     %R movement extractions:
%     allR_wheel = [];
%     allR_wheelInds = [];
%     zrlim = numel(wheelData(1:end-10))
% 
%     for zr = 1:length(wheelData)
%         if zr == zrlim;
%                 break;
%         else
%             if wheelData(zr+1) > wheelData(zr)
%                 allR_wheel(zr) = wheelData(zr);
%                 allR_wheelInds(zr) = wheelDataInds(zr); 
%             end
%         end
%     end
% 
%     %Eliminate wheel recordings that occurred before first trial start (.ns4)
%     validWheel_R = find(allR_wheelInds >= trueTrialStartsNS4(1));          %these are recording times(ms) when wheel is moving left; pull out first vals...
%     
%     validsR = [];
%     validsR_first = [];
%     for vr = 1:length(validWheel_R)
%         if vr == numel(validWheel_R)
%             break;
%         else
%             validsR(vr+1) = validWheel_R(vr+1) > validWheel_R(vr) + 20;     %length of vector cut-off is 70ms
% %             validsL(v+1) = validWheel_L(v+1) > validWheel_L(v) + 5; 
%         end
%     end
% 
%     validsR_first = validWheel_R(validsR == 1);                            %here are the first TSs for rightward vectors
%     ContData.behavior.wheel_R_first = validsR_first;
% 
%     
% %% EXTRACT VELOCITY DATA FROM CONT LEVER DATA. motor in: sLeverV = unfiltered velocities; sLeverVm = filtered
% % BUT my signal is already a velocity signal.
% 
% % Get advice from Josh
% 
% %     numSamples = size(ContData.behavior.sLeverData,2);
% %     tmpLeverData(1,:) = sgolayfilt(ContData.behavior.sLeverData(1,:),3,151); 
% %     tmpLeverData(2,:) = sgolayfilt(ContData.behavior.sLeverData(2,:),3,151);
%     numSamples = size(leverData);
%     tmpWheelData = sgolayfilt(leverData(1,:),3,151);
%     
% %     sLeverV = zeros(1,numSamples);
% %     dX = diff(tmpLeverData(1,:));
% %     dY = diff(tmpLeverData(2,:));
% %     sLeverV = sqrt( dX.^2 + dY.^2 );
% %     ContData.behavior.sLeverV = sLeverV;
% %     ContData.behavior.sLeverVm = sgolayfilt(ContData.behavior.sLeverV,3,501); %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
% 
%     wheelV = NaN(numSamples);
%     disp(' ');disp(' ');disp('Extracting velocity...');
%     dX = diff(tmpWheelData(1,:));
%     wheelV = sqrt(dX.^2);
%     
%     disp(' ');disp(' Complete. ');disp(' ');
%     ContData.behavior.wheelV = wheelV;
%     ContData.behavior.wheelVm = sgolayfilt(ContData.behavior.wheelV,3,501);
%     
% %     figure(1); plot(sLeverV);
%     figure; hold on; plot(wheelV);
%     clear wheelV;
% 
% %% FIND MOVEMENTS
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
%% WRITE THE REACH STRUCTURE OF THE CONTDATA STRUCTURE
%     ContData.behavior.reach = reach;
    
%% COMPLETED ALL ANALYSIS. SAVE AND START OVER

    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Completed file: ' filenamestr(1,1:length(filenamestr)-3)]);
    %BS added:
    filenamestrT = strcat(path, filenamestrE);
%     % save the data from this file
     
%      cd(fpath);
     
     save([targetName 'wheel_bh.mat'],'ContData');                              %BS
     disp(['saved as ' targetName 'wheel_bh.mat']);
                                                                
%     save(FR_file{n_i}, 'FR', 'FR_bine');
      save([targetName 'wheel_bh.mat'],'ContData');      
      disp(['saved as ' targetName 'wheel_bh.mat']);  
     
    disp('%-------------------------------------------------------------------');
    disp(' ');
    
