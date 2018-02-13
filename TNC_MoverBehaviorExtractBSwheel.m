function [ContData, filenamestrB] = TNC_MoverBehaviorExtractBSwheel(csvFile, filenamestr, targetName, dataRate, chan, elimTime, fpath)
%% This is BS's version of Josh's from Tonic_v2 
% ll.51-52 are lickInds that I commented out without piezo-data 11.26.15
% ll. 79 & 84, 113-117 commented out for chan.lick 11.26.15
% ll. 53-71 for parsing out solenoid clicks from trial starts for ch.rew 1.28.16

%  updated last 7.16.15; 11.26.15; 1.14.16; 1.28.16; 2.11.16 (OpenNSxNewBS)
%  called 

%  Note: removed the following vars since that is empty in my Blackrock system (7.16.15): NEVfiles
%  chan.x, chan.rew, 
%  ContData.behavior.threshInds = TS for laser pulse-on (divide by 30, 30kHz = sampling frequency for spikes)
% l. 51-76: chan.rew (129), now ONLY reflects every other TS (NOT proper sol reward inds) - same with trial starts; likely polluted with incorrect trial (bar goes off screen)


%  RE-ESTABLISH ACTIVE INPUTS FOR 11.26.15 DATA (1.14.16 & 1.28.16): 
%  motor in = chan.y (131) captured on ns4 file, sampling rate was 10kHz; divide by 10)
% l. 63 for me is not going to give me the number of trials since I do not laser stim every trial

% Also saves a bh.mat file that could be loaded into the next mfile for TS from recordings: TNC_ConvertTSDtoPopDataBS
% Currently, this mfile is called from MoverScript.m

%% FILE NAMES
    
    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Beginning file: ' filenamestr(1,1:length(filenamestr)-3)]);
    disp(' ');disp(' ');  

    dispOn = 0;
    
%% LOAD EVENT DATA - BS made changes 6/24/15
    newPath = fpath;
    cd(newPath);
    if exist('filenamestr', 'var');
        [filenamestrE, path] = uigetfile('*.nev*','select the .nev file', fpath)                %BS
    end
    
    dataEvents = openNEVBS_laser(filenamestrE,'read','nosave','nomat', 'report');               %BS

    rewIndsTmp = find(dataEvents.Data.Spikes.Electrode==chan.rew);                              %trial start, incorrect trial, and solenoid in detected (ch. 129)
%     plot(dataEvents.Data.Spikes.Electrode(1:10000))                                             %BS added to see the voltage trace (all .nev sigs in) w/ time
    
%         % proper filtering of .nev file & eliminate any pulses that were detected twice
%         rewTsTmp = dataEvents.Data.Spikes.Timestamps(rewIndsTmp)./30;                           %BS 1.14.16: will need to parse out trial start from solenoid here
%         rewTsRoundTmp = round(dataEvents.Data.Spikes.Timestamps(rewIndsTmp)./30);               %BS 1.29.16: account for 30kHz filtering for TSs 
% %         validRewInds = find(diff(rewTsTmp)>200);
%         validRewInds = find(diff(rewTsTmp)>50);                                                 %For BS task, 200ms will cause violation of trial start to some solenoid clicks
%         validRewInds = rewTsTmp; 
%         
%      % write valid data to structure
%       ContData.behavior.rewardInds = validRewInds;                   %These should be trial starts, incorrect bar off, and solenoids
%       
% %       %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%       %BS added 1.28.16 - Parse out trialStarts                       
%         starts = single(validRewInds);                                %vals 1,3,5... are trial starts 
%         trialStarts = starts(1:end);                                  %the last trial start will not have an end (end on next-to-last) (i.e. trials will end on odd index#) 
%         newT = [1:length(trialStarts)/2];                  
%         for t = 2:length(newT)  
%             newT(t) = (2*t)-1;
%             trialStartsNew = trialStarts(newT);                       %pull out vals 1,3,5...
%         end                                                           %solenoid var is now the solenoid TS
%        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        %BS added 1.28.16 - Parse out solenoid info
%         solenoid_trialStarts = single(validRewInds);                  %vals 1,3,5... are trial starts; vals 2,4,6 are solenoids
%         solenoid_trialStartsNew = solenoid_trialStarts(1:end-1);      %adjust for the last trial's sol click never happening below 
%         newS = [1:length(solenoid_trialStarts)/2];                  
%         for s = 1:length(newS)                                        
%             newS(s) = (2*s);
%             solenoid = solenoid_trialStartsNew(newS);                 %pull out vals 2,4,6...
%         end                                                           %solenoid var is now the solenoid TS
% %       %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     
%     ContData.behavior.startInds = trialStartsNew;
%     ContData.behavior.rewardInds = solenoid;
%     
%     thrIndsTmp = find(dataEvents.Data.Spikes.Electrode==chan.thr);                               %threshold channel: this var = indices (of TS) for when channel 130 (laser) pulses were detected; 
%     %Basically, this is a vector of TS when laser pulses were detected.
%         % eliminate any pulses that were detected twice
%         thrTsTmp = round(dataEvents.Data.Spikes.Timestamps(thrIndsTmp)./30);                     %divide this by 30 (spike acquisition filtering in kHz) so that vals will be accurate
%         validThrInds = find(diff(thrTsTmp)>10);
%     % write valid data to structure
%     ContData.behavior.threshInds = thrTsTmp(validThrInds);                                       %actual TSs of laser channel in
% 
%     % Sanity check
%     minRecording = ((ContData.behavior.threshInds(end) - ContData.behavior.threshInds(1))/1000)/60;
%     
%     disp(['Number of trial starts (.nev): ' num2str(numel(trialStartsNew))]);
%     disp(' ');disp(' '); 
%        
%     %lickInds = find(dataEvents.Data.Spikes.Electrode==chan.lick);                               %licking - disconnected piezo due to noise 11.26.15
%     %ContData.behavior.lickInds = round(dataEvents.Data.Spikes.Timestamps(lickInds)./30);

%% Load csv file: BS updated 1.29.16 - using as sanity check against .nev data
                     
    if exist('filenamestr', 'var');
        [filenamestrE, path] = uigetfile('*.csv*','select the .csv file', fpath)                %BS get the actual # correct trials from csv file
    end

     data.csvTrials = dlmread(csvFile,',',2,0);
     data.attemptedTrials = data.csvTrials(:,1);
     data.numRewardedTrials = data.csvTrials(end,1);
     data.trialStarts = data.csvTrials(:,2);
     data.trialEnds = data.csvTrials(:,4);
     
     %find the rewarded trials:
     attemptedTrials = data.attemptedTrials;
     trialStarts = data.trialStarts;
     numRew = unique(attemptedTrials);
     
    for r = 1:length(attemptedTrials)
        if r == numel(attemptedTrials)
            break                                                              %note: last trial info will be omitted
        else if attemptedTrials(r) < attemptedTrials(r+1)
                rewardedTrialsStarts(r) = trialStarts(r);  
            else if attemptedTrials(r+1) == attemptedTrials(r)
                    rewardedTrialsStarts(r) = NaN;                             %put a NaN in the 1st slot where trial is repeated
                end
            end
        end
    end


    rewardedTrialStartsOnes = find(isfinite(rewardedTrialsStarts));
    rewardedTrialStarts = rewardedTrialsStarts(rewardedTrialStartsOnes);        %Here are the csv TS for all rewarded trial starts

    numTrials = data.numRewardedTrials;
    data.rewardedTrialStarts = rewardedTrialStarts;                             %BS put actual number of rewarded trials into data struct
    data.numIncorrectTrials = numel(attemptedTrials) - numel(numRew);

     %     % work backwards from last trial
%     numTrials = numel(ContData.behavior.threshInds)                           %from Josh's code
%     
%     for p=numTrials:-1:1
%         ContData.behavior.block(p) = 7-floor(p./15);
%     end
    
    disp(['Number of trial starts (.csv): ' num2str(numel(r))]);
    disp(' ');disp(' '); 
    
    
    %%
    %         % proper filtering of .nev file & eliminate any pulses that were detected twice
%         rewTsTmp = dataEvents.Data.Spikes.Timestamps(rewIndsTmp)./30;                           %BS 1.14.16: will need to parse out trial start from solenoid here
%         rewTsRoundTmp = round(dataEvents.Data.Spikes.Timestamps(rewIndsTmp)./30);               %BS 1.29.16: account for 30kHz filtering for TSs 
% %         validRewInds = find(diff(rewTsTmp)>200);
%         validRewInds = find(diff(rewTsTmp)>50);                                                 %For BS task, 200ms will cause violation of trial start to some solenoid clicks
%         validRewInds = rewTsTmp; 
%         
%      % write valid data to structure
%       ContData.behavior.rewardInds = validRewInds;                   %These should be trial starts, incorrect bar off, and solenoids
%       
% %       %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%       %BS added 1.28.16 - Parse out trialStarts                       
%         starts = single(validRewInds);                                %vals 1,3,5... are trial starts 
%         trialStarts = starts(1:end);                                  %the last trial start will not have an end (end on next-to-last) (i.e. trials will end on odd index#) 
%         newT = [1:length(trialStarts)/2];                  
%         for t = 2:length(newT)  
%             newT(t) = (2*t)-1;
%             trialStartsNew = trialStarts(newT);                       %pull out vals 1,3,5...
%         end                                                           %solenoid var is now the solenoid TS
%        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        %BS added 1.28.16 - Parse out solenoid info
%         solenoid_trialStarts = single(validRewInds);                  %vals 1,3,5... are trial starts; vals 2,4,6 are solenoids
%         solenoid_trialStartsNew = solenoid_trialStarts(1:end-1);      %adjust for the last trial's sol click never happening below 
%         newS = [1:length(solenoid_trialStarts)/2];                  
%         for s = 1:length(newS)                                        
%             newS(s) = (2*s);
%             solenoid = solenoid_trialStartsNew(newS);                 %pull out vals 2,4,6...
%         end                                                           %solenoid var is now the solenoid TS
% %       %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     
%     ContData.behavior.startInds = trialStartsNew;
%     ContData.behavior.rewardInds = solenoid;
%     
%     thrIndsTmp = find(dataEvents.Data.Spikes.Electrode==chan.thr);                               %threshold channel: this var = indices (of TS) for when channel 130 (laser) pulses were detected; 
%     %Basically, this is a vector of TS when laser pulses were detected.
%         % eliminate any pulses that were detected twice
%         thrTsTmp = round(dataEvents.Data.Spikes.Timestamps(thrIndsTmp)./30);                     %divide this by 30 (spike acquisition filtering in kHz) so that vals will be accurate
%         validThrInds = find(diff(thrTsTmp)>10);
%     % write valid data to structure
%     ContData.behavior.threshInds = thrTsTmp(validThrInds);                                       %actual TSs of laser channel in
% 
%     % Sanity check
%     minRecording = ((ContData.behavior.threshInds(end) - ContData.behavior.threshInds(1))/1000)/60;
%     
%     disp(['Number of trial starts (.nev): ' num2str(numel(trialStartsNew))]);
%     disp(' ');disp(' '); 
%        
%     %lickInds = find(dataEvents.Data.Spikes.Electrode==chan.lick);                               %licking - disconnected piezo due to noise 11.26.15
%     %ContData.behavior.lickInds = round(dataEvents.Data.Spikes.Timestamps(lickInds)./30);


%% LOAD CONTINUOUS LEVER (WHEEL) DATA
    
switch dataRate
    
    case 'ns4'    
%          filenamestrE = [filenamestr(1,1:length(filenamestr)-3) dataRate]      %BS change to accept my own files
        [filenamestrE, path] = uigetfile('*.ns4','select the .ns4 file', fpath); %BS works, but unneccessary extra opening of file by hand

        
        Ns4DATA = openNSxBS('report','read',filenamestrE);
%         Ns4DATA = openNSxNewBS('report','read',filenamestrE);               %updated BS since Josh's update 2.10.16

        clear leverData sLeverData tmpLeverData
        
%           xChan = find(Ns4DATA.MetaTags.ChannelID==chan.x);                 %my x channel is empty  
        yChan = find(Ns4DATA.MetaTags.ChannelID==chan.y);                   %Motor in
%         lChan = find(Ns4DATA.MetaTags.ChannelID==chan.lick);
        
%          leverData(1,:) = decimate(Ns4DATA.Data(xChan,:),10);
        leverData(2,:) = decimate(Ns4DATA.Data(yChan,:),10);                %Correct for 10kS/s rate
%         rawLick = decimate(Ns4DATA.Data(lChan,:),10);

    case 'ns3'
        Ns4DATA = openNSx('report','read',filenamestrE);
        clear leverData sLeverData tmpLeverData
        leverData(1,:) = decimate(Ns4DATA.Data(chan.x,:),2);
        leverData(2,:) = decimate(Ns4DATA.Data(chan.y,:),2);
        rawLick = decimate(Ns4DATA.Data(3,:),10);
        
        
    case 'ns2'
%         filenamestrE = [filenamestr(1,1:length(filenamestr)-3) dataRate]  %BS
        filenamestrE = uigetfile('*.ns4','select the .ns4 file', dataRate);
        Ns2DATA = openNSx('report','read',filenamestrE);
        clear leverData sLeverData tmpLeverData
        xChan = find(Ns2DATA.MetaTags.ChannelID==chan.x);                   %BS
        yChan = find(Ns2DATA.MetaTags.ChannelID==chan.y);
        lChan = find(Ns2DATA.MetaTags.ChannelID==chan.lick);
%         leverData(1,:)  = Ns2DATA.Data(xChan,:);                          %BS
        leverData(2,:)  = Ns2DATA.Data(yChan,:);
        rawLick         = Ns2DATA.Data(lChan,:);
        
end
    
    sLeverData(1,:) = sgolayfilt(leverData(1,:),9,101);  %motor %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
    sLeverData(2,:) = sgolayfilt(leverData(2,:),9,101);  %motor

    sWheelData = sLeverData;
    ContData.behavior.sWheelData = sWheelData;
     
    ContData.behavior.sLeverData = sLeverData;          %number of data points are = numel(sWheelData)/10 (and rounded up) for correction of 10kHz sampling rate

%     difLick = diff(rawLick);
%     evLick=zeros(1,numel(rawLick));
% 
%     ContData.behavior.evLick    = evLick;
%     ContData.behavior.rawLick   = rawLick;

    

%% BS 2.1.16

% Find any leftward movement (-going voltages) or rightward movement (+)
  % Find "center" of raw wheel voltages
    centers = median(sWheelData); 
    wheelMin = min(sWheelData);    %just to see
    wheelMax = max(sWheelData);    %just to see

    wheelData = single(sWheelData); 
    wheelData = sWheelData(2,:);
    wheelDataInds = [];
    for w = 1:length(wheelData);
        wheelDataInds(:,w) = w;
    end

    newWheelMat(1,:) = wheelDataInds;
    newWheelMat(2,:) = wheelData;

    % L movement extractions:
    allL_wheel = [];
    allL_wheelInds = [];
    for z = 1:length(wheelData)
        if wheelData(z) > wheelData(z+1)
            if z+1 == length(wheelData)
                break;
            else
                allL_wheel(z) = wheelData(z);
                allL_wheelInds(z) = wheelDataInds(z);
%                 validWheel = find(allL_wheelInds
            end
        end
    end

    %Need to eliminate wheel recordings that occurred before first trial start (.nev)
    validWheel_L = find(allL_wheelInds' >= trialStartsNew(1));              %These are recording times(ms) when wheel is moving left; pull out first vals...
    
    validsL = [];
    validsL_first = [];
    for v = 1:length(validWheel_L)
        if v == numel(validWheel_L)
            break;
        else
            validsL(v+1) = validWheel_L(v+1) > validWheel_L(v) + 5;         %Cut-off is 5ms;
        end
    end

    validsL_first = validWheel_L(validsL == 1);                             %Here are the first TSs for leftward vectors
    
    %Write the vectors for all L wheel movements to the ContData structure
%     ContData.behavior.validWheel_L = validWheel_L;
    ContData.behavior.newWheelMat = newWheelMat;
%     ContData.behavior.allL_wheel = single(allL_wheel);                    %Great, but these aren't necessary
%     ContData.behavior.allL_wheelInds = single(allL_wheelInds);
    
    allL_wheelIndsShort = find(allL_wheelInds>0);
%     ContData.behavior.allL_wheelIndsShort = single(allL_wheelIndsShort); 
%     validWheel_L = find(allL_wheelIndsShort' >= trialStartsNew(1));       %These are recording times(ms) when wheel is moving left (without the zeros)
    ContData.behavior.validsL_first = validsL_first;

    %R movement extractions:
    allR_wheel = [];
    allR_wheelInds = [];
    zrlim = numel(wheelData(1:end-10))

    for zr = 1:length(wheelData)
        if wheelData(zr+1) > wheelData(zr)
%             if zr == numel(wheelData(1:end-10)) takes looop way too long; issue is that if vals keep adding in zero, "exceeds matrix dims"
            if zr == zrlim;
                break;
            else
                
                allR_wheel(zr) = wheelData(zr);
                allR_wheelInds(zr) = wheelDataInds(zr);
                
            end
        end
    end

    %Need to eliminate wheel recordings that occurred before first trial start (.nev)
    validWheel_R = find(allR_wheelInds' >= trialStartsNew(1));              %These are recording times(ms) when wheel is moving left; pull out first vals...
    
    validsR = [];
    validsR_first = [];
    for vr = 1:length(validWheel_R)
        if vr == numel(validWheel_R)
            break;
        else
            validsR(vr+1) = validWheel_R(vr+1) > validWheel_R(vr) + 5;       %cut-off is 5ms
%             validsL(v+1) = validWheel_L(v+1) > validWheel_L(v) + 5; 
        end
    end

    validsR_first = validWheel_R(validsR == 1);                              %Here are the first TSs for rightward vectors
    ContData.behavior.validsR_first = validsR_first;
    
%% EXTRACT VELOCITY DATA FROM CONT LEVER DATA. BS: motor in: sLeverV = unfiltered velocities; sLeverVm = filtered

    numSamples = size(ContData.behavior.sLeverData,2);
    tmpLeverData(1,:) = sgolayfilt(ContData.behavior.sLeverData(1,:),3,151); 
    tmpLeverData(2,:) = sgolayfilt(ContData.behavior.sLeverData(2,:),3,151);

    sLeverV = zeros(1,numSamples);

    disp(' ');disp(' ');disp('Extracting velocity...');

    dX = diff(tmpLeverData(1,:));
    dY = diff(tmpLeverData(2,:));

    sLeverV = sqrt( dX.^2 + dY.^2 );

    disp(' ');disp(' Complete. ');disp(' ');

    ContData.behavior.sLeverV = sLeverV;
    ContData.behavior.sLeverVm = sgolayfilt(ContData.behavior.sLeverV,3,501); %y = sgolayfilt(x,k,f) applies Savitzky-Golay FIR smoothing filter to x. k is the polynomial order and must be < f, frame size.
    
%     figure(1); plot(sLeverV);
    
    clear sLeverV;

%% FIND MOVEMENTS
    method = 'vel';
    numC = 5;
    clear progSt* reach

    currX       = ContData.behavior.sLeverData(1,:);
    currY       = ContData.behavior.sLeverData(2,:);
    currV       = ContData.behavior.sLeverVm; 

    pre     = 10;
    post    = 10;       
    minSpace = 250;
    count = 1;

    % threshold the velocities
    switch dataRate
        case 'ns4'                                                            %currV are the velocity values
            allValidSamps   = find(currV(1:numel(currV)-(minSpace.*2))>0.35); %find point where filtered velocities begin at -500 and >.35 (will be that index and above)
        case 'ns2'
            allValidSamps   = find(currV(1:numel(currV)-(minSpace.*2))>1.5);
    end
    
    if numel(allValidSamps)>0

        switch method

            case 'vel'

    %             progStartTMP(count,1)  = currStamps(allValidSamps(1));
                progStartTMP(count,2)  = currX(allValidSamps(1));
                progStartTMP(count,3)  = currY(allValidSamps(1));           %find currY positions for velocities >.35
                progStartTMP(count,4)  = allValidSamps(1);

                for j=2:numel(allValidSamps)                                %allValidSamps are the indexed velocities >.35

                    if allValidSamps(j)>allValidSamps(j-1)+minSpace         %if the indexed velocities are increasing by at least 250

                        postN = find(currV(allValidSamps(j-1):allValidSamps(j))<1,1,'first');  %postN = find the first velocity values at those positions
                        if numel(postN)<1
                            figure(3); clf;                                                         %if postN has less than 1 value, plot                
                            plot(progStartTMP(count,4):allValidSamps(j),currV(progStartTMP(count,4):allValidSamps(j)),'k');  
                            hold on; plot(progStartTMP(count,4),currV(progStartTMP(count,4)),'bo',allValidSamps(j-1),currV(allValidSamps(j-1)),'ro');
                            pause(0.1);

                            disp('Cannot find stop');                                               %and set postN = post (which was defined as 10)
                            postN=post;
                        end
    %                     progStopTMP(count,1)   = currStamps(allValidSamps(j-1)+post);
                        progStopTMP(count,2)   = currX(allValidSamps(j-1)+postN);                %But if postN has >1 value, create a matrix (progStopTMP) and find currY at index allValidSamps at the previous positions and add 10 
                        progStopTMP(count,3)   = currY(allValidSamps(j-1)+postN);
                        progStopTMP(count,4)   = allValidSamps(j-1)+postN;
                        count                  = count+1;

                        preN = find(currV(allValidSamps(j-1) : allValidSamps(j))<1,1,'last');
                        if numel(preN)<1
                            disp('Cannot find start');
    %                     progStartTMP(count,1)  =
    %                     currStamps(allValidSamps(j)-pre)                                        %Now create the startTMP matrix
                            progStartTMP(count,2)  = currX(allValidSamps(j)-pre);
                            progStartTMP(count,3)  = currY(allValidSamps(j)-pre);
                            progStartTMP(count,4)  = allValidSamps(j)-pre;                    
                        else
    %                     progStartTMP(count,1)  = currStamps(allValidSamps(j-1)+preN);
                            progStartTMP(count,2)  = currX(allValidSamps(j-1)+preN);
                            progStartTMP(count,3)  = currY(allValidSamps(j-1)+preN);
                            progStartTMP(count,4)  = allValidSamps(j-1)+preN;
                        end
                    end

                    if j==numel(allValidSamps)
    %                     post = find(currV(allValidSamps(j):allValidSamps(j)+minSpace)<0.5,1,'first');
    %                     progStopTMP(count,1)   = currStamps(allValidSamps(j)+post);
                        progStopTMP(count,2)   = currX(allValidSamps(j)+post);
                        progStopTMP(count,3)   = currY(allValidSamps(j)+post);
                        progStopTMP(count,4)   = allValidSamps(j)+post;
                    end

                end

                count = 1;
                for k = 1:size(progStartTMP,1)

                    if k==1
                        reach.init = 1;
                    end

                    % reaches must be at least 50 ms long
                    if progStopTMP(k,4)-progStartTMP(k,4)>=90 & progStartTMP(k,4)>minSpace

                        trajAngle   = atan2(progStopTMP(k,3)-progStartTMP(k,3),progStopTMP(k,2)-progStartTMP(k,2));

                        if (pdist2([progStopTMP(k,2),progStopTMP(k,3)],[mean(currX),mean(currY)]) > pdist2([progStartTMP(k,2),progStartTMP(k,3)],[mean(currX),mean(currY)]))
                            reach.out(count) = 1;
                        else
                            reach.out(count) = 0;
                        end
                        velTraj = ContData.behavior.sLeverV(progStartTMP(k,4) : progStopTMP(k,4));
                        xVals = ContData.behavior.sLeverData(1,progStartTMP(k,4) : progStopTMP(k,4));
                        yVals = ContData.behavior.sLeverData(2,progStartTMP(k,4) : progStopTMP(k,4));

                        reach.start(count,:)  = progStartTMP(k,:);
                        reach.stop(count,:)   = progStopTMP(k,:);
                        reach.angle(count,1)  = trajAngle;
                        reach.dist(count,1)   = trapz(velTraj);
                        reach.dist(count,2)   = pdist2(progStartTMP(k,2:3) , progStopTMP(k,2:3));

                        tmp = findpeaks(velTraj);
                        reach.numpks(count,1) = numel(tmp.loc);
                        reach.dur(count,1)    = progStopTMP(k,4) - progStartTMP(k,4);
                        reach.vel(count,1)   = max(velTraj);
                        reach.vel(count,2)   = trapz(velTraj) ./ reach.dur(count,1);
                        reach.vel(count,3)   = var(velTraj);
                        reach.vel(count,4)   = find(velTraj==max(velTraj),1);

                        reach.acc(count,1)   = max(diff(velTraj));
                        reach.acc(count,2)   = mean(diff(velTraj));
                        reach.acc(count,3)   = max(diff(velTraj(1:90))); % max in first 90 ms of movement

                        reach.tort(count,1)  = reach.dist(count,1) ./ pdist2([progStopTMP(k,2),progStopTMP(k,3)],[progStartTMP(k,2),progStartTMP(k,3)]);

                        % find max displacement of the reach
                        xVals = xVals - xVals(1);
                        yVals = yVals - yVals(1);
                        
                        
                        
                        valThrInd = find(ContData.behavior.threshInds>reach.start(count,4) & ContData.behavior.threshInds<reach.stop(count,4));
                        if numel(valThrInd)>0
                            reach.rewarded(count)= 1;
                        else
                            reach.rewarded(count)= 0;                       
                        end

                        count                 = count+1;

                    end            
                end


        end

        disp(['Num valid reaches: ' num2str(count-1)]);

        level = ones(1,size(progStartTMP,1)).*0.35;
        figure(3); clf;
        plot(ContData.behavior.sLeverV,'k'); title(['sLeverV']);hold on;
        plot(progStartTMP(:,4),level,'r^'); title(['reach starts']);
        plot(progStopTMP(:,4),level,'bo');  title(['reach stops']);

        [vals,inds] = sort(reach.dur,'ascend');    

        reach.numReaches = size(reach.start,1);

        dims = floor(sqrt(reach.numReaches));
        scaler = 3000;

        if dispOn
            
            figure(5); clf;
        
            for p = 1:1:dims.^2

                l = inds(p);

                offset = l.*ones(1,numel(reach.start(l,4):reach.stop(l,4)));
                [x,y] = ind2sub(dims,p);

                figure(5);
                plot((x.*scaler),(y.*scaler),'k.','MarkerSize',8);
                hold on;
                plot(currX(reach.start(l,4):reach.stop(l,4))-currX(reach.start(l,4))+(x.*scaler),currY(reach.start(l,4):reach.stop(l,4))-currY(reach.start(l,4))+(y.*scaler) , 'k','LineWidth',1,'Color',[p/reach.numReaches 0.67-(p/reach.numReaches)*0.67 1-(p/reach.numReaches)]);
                if i==dims.^2
                    axis tight; axis off;
                end

                startInd = round(reach.start(l,4)./2);
                stopInd = round(reach.stop(l,4)./2);

            end
        end
        
        if reach.numReaches>5
                figure(7); clf;
                subplot(311);
                p = polyfit(reach.dist(:,1),reach.vel(:,1),1);
                plot(reach.dist,reach.vel(:,1),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
                xlabel('Total distance');
                ylabel('Peak Velocity');

                subplot(312);
                p = polyfit(reach.dist(:,1),reach.dur,2);
                plot(reach.dist,reach.dur,'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
                xlabel('Total distance');
                ylabel('Movement Duration (ms)');

                subplot(313);
                p = polyfit(reach.dist(:,1),reach.acc(:,3),1);
                plot(reach.dist,reach.acc(:,3),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
                xlabel('Total distance');
                ylabel('Initial acceleration');

                figure(8); clf;
                rose(reach.angle);

                drawnow;
        end
        
    else
        
        reach.numReaches=0;
    
    end
    

    
%% DIMENSION REDUCTION OF REACH PARAMETERS

if reach.numReaches > 5
    clear matFor*

    matForDR = [reach.vel reach.acc reach.tort reach.dur reach.dist];

    for j=1:size(matForDR,2)
        matForCov(:,j) = (matForDR(:,j) - mean(matForDR(:,j))) ./ std(matForDR(:,j)) ;
    end

    figure(11);
    imagesc(corr(matForCov),[0 1]);
    map = colormap(TNC_CreateRBColormap(1024,'cpb'));
%     colormap(1-map);

    [eV,eD] = eig(cov(matForCov));
    figure(12); plot(cumsum(diag(eD))./sum(diag(eD)),'ko-');

    for m=1:size(matForDR,1)
        reach.pca(1,m) = dot(matForCov(m,:),eV(:,size(matForCov,2))');
        reach.pca(2,m) = dot(matForCov(m,:),eV(:,size(matForCov,2)-1)');
        reach.pca(3,m) = dot(matForCov(m,:),eV(:,size(matForCov,2)-2)');
    end
end

%% CALCULATE META-REACH PARAMETERS
    if reach.numReaches > 5
        [winReach] = TNC_ReachVigorWindow(reach.pca(1,:),reach.numReaches,9); % capture fluctuations along the maximally variant dimension
        ContData.behavior.winReach = winReach;
    end
    
%% WAS THE REACH REWARDED?
if reach.numReaches > 5
    for p=1:reach.numReaches

        tmp = find(ContData.behavior.threshInds>reach.start(p,4) & ContData.behavior.threshInds<reach.stop(p,4));
        if numel(tmp)>0
            reach.rewarded(p)=1;
        else
            reach.rewarded(p)=0;
        end
        
    end
end    
%% WRITE THE REACH STRUCTURE OF THE CONTDATA STRUCTURE
    ContData.behavior.reach = reach;
    
%% COMPLETED ALL ANALYSIS. SAVE AND START OVER

    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Completed file: ' filenamestr(1,1:length(filenamestr)-3)]);
    %BS added:
    filenamestrB = strcat(newPath, filenamestrE);
%     % save the data from this file
     
     cd(fpath);
     
     save([targetName '_bh.mat'],'ContData');                              %BS
     disp(['saved as ' targetName '_bh.mat']);
                                                                
%     save(FR_file{n_i}, 'FR', 'FR_bine');
%       save([targetName '_bh.mat'],'ContData');      
      disp(['saved as ' targetName '_bh.mat']);  
     
    disp('%-------------------------------------------------------------------');
    disp(' ');
    
