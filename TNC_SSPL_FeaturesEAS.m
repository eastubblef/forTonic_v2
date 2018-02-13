% _________________________________________________________________________
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / Janelia
% 
% BUG REPORTING: josh@dudmanlab.org
% FURTHER INFORMATION: www.dudmanlab.org
% Copyright (C) 2015 by Howard Hughes Medical Institute.
% _________________________________________________________________________
% _________________________________________________________________________
% 
% DETAILS: extract spike waveforms and features suitable for clustering with the TNC_SS_GUI interface
% 
% INPUTS:
% fileName: name of the ns5 file used from processing spike data
% arrayType: name of the array - depends upon function TNC_RemapElecPos
%  *** NOTE: if the array type string is not recognized all sites are processed independently.
% chunk: how many seconds of data should be processed per chunk
% dispOn: display online progress in figures [[[SLOW]]]
% targetName: base name of target files containg the featStruct (*_ft) and sessionStruct (*_ss) outputs
% maxSegs: if 0 process all chunks, otherwise process maxSegs number of chunks
% 
% EXAMPLE:
% [sessionStruct,featStruct] = TNC_SSPL_Features('15_05_27m28003.ns5','NN_h64_dhs',90,0,'g15_05_27m28003?,0);
% 

function [sessionStruct,featStruct] = TNC_SSPL_FeaturesEAS(fileName,arrayType,chunk,dispOn,targetName)

sampPrior = 25;
pathName = './'; % can be changed or otherwise specified if necessary (required for the update to loader v6.2.0 from BlackRock)

%% LOAD A CONFIGURATION FILE IF IT EXISTS

tmp = dir('config_analy*');
if numel(tmp)>0
    S=load(tmp.name);
    % ===================================
    % Choose filtering method
    % ===================================
    par.detect_fmin = S.config_an.detect_fmin;
    par.detect_fmax = S.config_an.detect_fmax;
    par

    % ===================================
    % Choose thresholding method
    % ===================================
    thresh = S.config_an.thresh;

    art_flag = S.config_an.art_flag;
    
else

    %=============================
    % USE ELLIPTICAL FILTER w/ quiroga settings
    par.detect_fmin = 300
    par.detect_fmax = 3000;

    thresh = 4
    art_flag = 0;
    
end

    par.sr = 30000;

%% PRE-PROCESSING
    % GET INFO ABOUT THE CURRENT RECORDING
    Ns5DATA = openNSx('report',[pathName fileName]);
    disp(' ');
    disp(' ');
%     disp(['Data contains ' num2str(numSegs) ' x ' num2str(chunk) ' seconds long segments.']);

    % determine how to calculate segments
    switch numel(chunk)
        case 1
            numSegs = ceil(Ns5DATA.MetaTags.Duration ./ chunk); %%2.6.18 Changed to actually call the duration properly from the field
            timeOffset = 0;
            for k=1:numSegs
                seg_txt{k} = ['t:' num2str((k-1)*chunk) ':' num2str(k*chunk)]
                chunkOff(k) = (k-1)*chunk;
            end
            
        case 3
            numSegs = ceil( (chunk(2)-chunk(1)) ./ chunk(3) );
            timeOffset = chunk(1);
            for k=1:numSegs
                seg_txt{k} = ['t:' num2str(((k-1)*chunk(3))+chunk(1)) ':' num2str((k*chunk(3))+chunk(1))]
                chunkOff(k) = ((k-1)*chunk(3))+chunk(1);
            end
            
        otherwise
            numSegs = size(chunk,1);
            timeOffset = 0;
            disp('***NOTE***');
            disp(' ');
            disp('Assuming 0 timeOffset because user will know how to correct offsets for subsequent data processing');
            disp(' ');
            disp('***NOTE***');
            for k=1:numSegs
                seg_txt{k} = ['t:' num2str(chunk(k,1)) ':' num2str(chunk(k,2))]
                chunkOff(k) = chunk(k,1);
            end
    end
    
    sessionStruct.chunk     = chunkOff;
    sessionStruct.arrayType = arrayType;
    sessionStruct.fileName  = fileName;

    
    
    disp([num2str(numSegs) ' segments to be processed.']);
    disp(' ');
    disp(' ');

    totChar = numel(fileName);
    
    if dispOn==0
        loadBehavior=0;
    else
        d = dir([fileName(1:totChar-4) '_bh.mat']);
        if numel(d)>0
            loadBehavior = 1;
        else 
            disp('No associated behavioral data was found.');
            loadBehavior = 0;
        end
    end
    
%% EXTRACT MULTIDIMENSIONAL EVENTS

    wfs = [];
    wts = [];
    chan = [];
    inds = [];
    vals = [];

    % LOAD BEHAVIORAL DATA TO DISPLAY WITH THE SPIKING DATA
    totChar = numel(fileName);
    ppBehavFileName = [fileName(1:totChar-4) '.mat']
    
    if loadBehavior
        eval(['load ' fileName(1:totChar-4) '.mat']);
        behaviorTrue = 1;
    else
        behaviorTrue = 0;
    end
    
    for k=1:numSegs

    % LOAD ALL CHANNEL DATA 
    % *** NOTE: READ/WRITE TIME SLOWS THIS DOWN, BUT REDUCES MEMORY DEMANDS - CONSIDER CHANGING...
        timeStr = ['t:' num2str((k-1)*chunk) ':' num2str(k*chunk)];
        timeStr = seg_txt{k};
        disp(' ');
        disp(' ');
        disp(['Loading data over the range of ' timeStr]);
        disp(' ');
        disp(' ');
        Ns5DATA = openNSx('report','read',[pathName fileName],timeStr,'sec');
        
    % SMOOTH/FILTER DATA 
        numChan     = numel(find(Ns5DATA.MetaTags.ChannelID<129));

    % SOME CHECKING OF DATA FOR PROBES WHERE MULTIPLE SITES ARE USED FOR SORTING
        if ~strcmp(arrayType,'SingleSites')
            allRecChan = double(Ns5DATA.MetaTags.ChannelID);
            if min(allRecChan)>64 % in case the recording was done in Bank B
                allRecChan = allRecChan-64;
            end
        else
            if k==1 % pre-allocate memory for the structure
                sessionStruct.seg(numSegs).shank(numChan).spikes    = [];
                featStruct.seg(numSegs).shank(numChan).col          = 0;
            end
        end

        [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(1,arrayType);
        numRows = size(electrodeMatrix,1);
        
        if numRows>1

            %________________________________________________________
            %_______CODE FOR MULTISITE ARRAYS________________________
            %________________________________________________________
                    
                %=============================
                % GET THE ELECTRODE MAPPINGS
                [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(1,arrayType);
                
                for j=1:size(electrodeMatrix,2) %shanks                                       
                    
                    wfs = [];
                    wts = [];
                    chan = [];
                    inds = [];
                    vals = [];
                    
                    disp('_________________________________________________________________________');
                    disp('_________________________________________________________________________');
                    disp(['...Seg ' num2str(k) ' of ' num2str(numSegs) ' | shank ' num2str(j)]);
                    disp(' ');

                    spkData = zeros(size(electrodeMatrix,1),size(Ns5DATA.Data,2));
                    
                %=============================
                %=============================
                    disp(['1) Filtering the raw data...']);
                    
                    for i=1:size(electrodeMatrix,1) % Filter is so fast now that it is worse to use parallel filtering unless chunk is long
                        if numel( find(allRecChan==electrodeMatrix(i,j)) ) > 0  
                            spkData(i,:) = TNC_FilterData2(double(Ns5DATA.Data(find(allRecChan==electrodeMatrix(i,j)),:)),par);
                        else
                            disp(['***ALERT*** Channel ' num2str(electrodeMatrix(i,j)) ' was not recorded in this segment...']);
                            spkData(i,:) = zeros(size(Ns5DATA.Data(1,:)));
                        end                        
                    end

                    % COMMON MODE REJECTION
                    avgSeg      = median(spkData,1);
                    for i=1:size(electrodeMatrix,1)
                        spkData(i,:) = (spkData(i,:)-avgSeg) ; % ./ 4; % division by four is to get into the units of uV
                    end
                                                            
                %=============================
                %=============================
                    disp(['2) Detecting events...']);
                    if art_flag==1
                        [problem_samps] = TNC_FindArtifacts(sum(spkData,1),30); %blanks 1 ms of artifact
                        spkData(:,problem_samps)=randn(size(spkData,1),numel(problem_samps));
                    else
                        problem_samps = [];
                    end
                    
                    parfor i=1:size(electrodeMatrix,1) % electrodes
                    
                        %=============================
                        % DETECT EVENTS ON THIS CHANNEL
                        if numel(find(allRecChan==electrodeMatrix(i,j)))==1
                            [parEv(i).events] = TNC_SSPL_EventDetect(spkData(i,:),30,thresh,1);
                        else
                            disp(['***ALERT*** Channel ' num2str(electrodeMatrix(i,j)) ' has 0 events in this segment...']);
                            parEv(i).events.inds = [];
                        end
                    end
                    
                    for i=1:size(electrodeMatrix,1) % electrodes
                        
                        %=============================
                        % ADD TO TOTAL LIST OF EVENTS
                        inds = [inds parEv(i).events.inds];
                        vals = [vals spkData(i,parEv(i).events.inds)];
                        chan = [chan ones(1,numel(parEv(i).events.inds)).*i];
                        
                    end
                
                    %=============================
                    % FIND UNIQUE EVENTS ACROSS THE SHANK (ALL EVENT MUST BE A FIXED TIME APART)
                    [allInds,origIs,uIs]    = unique(inds);
                    allVals                 = vals(origIs);
                    allChan                 = chan(origIs);
                    
                    % considering 12 samples(400 us) to be outside of peak
                    % shift time seen from what is truly the same unit
                    valids = find(diff(allInds)>12);
                    valids = [ 1 , valids + 1]; % deal with the N-1 size of diff, first spike is by definition well separated

                    % for closer indices need to make a decision about
                    % whether its 2 different units that happen to be
                    % nearly simultaneous or just the same unit with
                    % components shifted in time relative to each other
                    invalids = find(diff(allInds)<=12);
                    invalids = invalids + 1; % deal with the N-1 size of diff
                    
                    %=============================
                    % RECONCILE THE INVALID SAMPLES
                    % first look for sequential invalids with no intervening valid indices (these should be treated as a group).
                    % for isolated invalid indices just need to compare those to the previously valid index, but it probably won't matter because the valid index will, by the definition of invalid, capture the invalid waveform.
                    m=2; startK=1; addValids=[];
                    while m < numel(invalids)
                        
                        if invalids(m)-invalids(m-1)==1
                            % in a sequence of short gaps with no
                            % intervening valid gaps
                            
                        else
                            endK = m-1;
                            
                            % evaluate the range invalids(startK)-1 (i.e. last valid) to invalids(endK) to find the
                            values      = allVals(invalids(startK)-1:invalids(endK));
                            chans       = allChan(invalids(startK)-1:invalids(endK));
                            indices     = allInds(invalids(startK)-1:invalids(endK));
                            lSpkLoc     = find(values==min(values)); % look for the largest spike out of the set
                            if lSpkLoc == 1 %largest spike was at the valid point
                                if chans(1)~=chans(2)
                                  addValids = [addValids invalids(startK)];
                                end
                            else
                                addValids = [addValids invalids(lSpkLoc+startK-2)];
                            end
                            
                            startK = m;
                        end
                        
                        m=m+1;

                    end
                
                    if 1
                        eventInds   = unique(allInds(valids)); % use unique to sort the indices
                        disp(['     >>> Total events: ' num2str(size(eventInds)) '...(' num2str(numel(addValids)) ')']);
                    else
                        eventInds   = unique(allInds([valids addValids])); % use unique to sort the indices
                        disp(['     >>> Total events: ' num2str(size(eventInds)) '...(' num2str(numel(addValids)) ')']);
                    end
                    
                %=============================
                %=============================
                    disp('3) Extract & quantify spike waveforms...');
                    [spikes]    = TNC_SSPL_EventExtractME(spkData,eventInds,[sampPrior,45]);
                    
                    if k==1
%                         [features]  = TNC_SSPL_EventQuantME(spikes.wfs, [], round(spikes.inds./30), sampPrior, size(electrodeMatrix,1));
                        [features] = TNC_SSPL_EventQuantME(spikes.wfs, [], round(spikes.inds./30) , sampPrior, size(electrodeMatrix,1));                        
                    else
%                         [features]  = TNC_SSPL_EventQuantME(spikes.wfs, features.pca , round(spikes.inds./30), sampPrior, size(electrodeMatrix,1));                        
                        [features] = TNC_SSPL_EventQuantME(spikes.wfs, featStruct.seg(1).shank(j).proj , round(spikes.inds./30) , sampPrior, size(electrodeMatrix,1));
                    end
                    
                %=============================
                %=============================
                    disp('4) Updating the feature structure with new data...');
                        
                        featStruct.seg(k).shank(j).inds     = spikes.inds + (chunkOff(k).*30000);
                        featStruct.seg(k).shank(j).ts       = round(featStruct.seg(k).shank(j).inds./30);
                        featStruct.seg(k).shank(j).id       = zeros(numel(spikes.inds),1);
                        featStruct.seg(k).shank(j).params   = features.params;
                        featStruct.seg(k).shank(j).proj     = features.proj;
                        
                        if k==1 && j==1
                            featStruct.paramNames       = features.paramNames;
                            featStruct.rSp              = rSp;
                            featStruct.cSp              = cSp;
                            featStruct.electrodeMatrix  = electrodeMatrix;
                            sessionStruct.resolution    = spikes.resolution;
                            sessionStruct.winL          = spikes.winL;
                            sessionStruct.winR          = spikes.winR;
                        end

                        sessionStruct.seg(k).shank(j).wfs   = spikes.wfs;
                        sessionStruct.seg(k).shank(j).inds  = spikes.inds;
                    disp('_________________________________________________________________________');
                    disp(' ');
                    disp(' ');                    
                    disp(' ');
                    
                end

            %________________________________________________________
            %_______CODE FOR ISOLATED SITE ARRAYS____________________
            %________________________________________________________

        else
                
                disp(' ');
                disp(' ');
                disp('_________________________________________________________________________');
                disp('...BEGIN FILTERING DATA');
                        
                %=============================
                % FILTERING DATA
                disp('1) Filtering data on all channels...');                
                    parfor eInd=1:numChan % shanks | sites
                        [bandPassedData] = TNC_FilterData2(double(Ns5DATA.Data(eInd,:)),par);
                        spkData_all(eInd,:) = bandPassedData;                        
                    end
                    avgSeg = median(spkData_all,1);

                %=============================
                % COMMON MODE REJECTION
                disp('1b) Common mode rejection...');                
                    for eInd=1:numChan
                        spkData(eInd,:) = spkData_all(eInd,:) - avgSeg;
                        spks_struct(eInd).spikes = [];
                        spks_struct(eInd).features = [];
                    end
                
                %=============================
                % EXTRACT WAVEFORMS FROM ALL TRODES FOR EACH UNIQUE EVENT
                disp('2+3) Detecting, extracting, quantifying events...');                
                    if art_flag==1
                        [problem_samps] = TNC_FindArtifacts(sum(spkData,1),10); % blanks out 10 samples (0.33 ms) of artifact data 
                        spkData(:,problem_samps)=randn(size(spkData,1),numel(problem_samps));
                    else
                        problem_samps = [];
                    end
                    
                    parfor eInd=1:numChan % shanks | sites
                        [events] = TNC_SSPL_EventDetect(spkData(eInd,:),30,4,1); % basically always has to be 4 to ensure not getting overloaded
%                         [events] = TNC_SSPL_EventDetect(spkData(eInd,:),30,thresh,1);
                        if numel(events.ts)>10
                            [spikes]    = TNC_SSPL_EventExtract(spkData(eInd,:),events.inds,[15,45],dispOn);
                            [features]  = TNC_SSPL_EventQuantSE(spikes,[]);
                            spks_struct(eInd).spikes    = spikes
                            spks_struct(eInd).features  = features
                        end
                    end
                
                %=============================
                % UPDATE STRUCTURE FOR SPIKE SORTING
                disp('4) Updating the feature structure with new data...');
                    for eInd=1:numChan                        
                        featStruct.seg(k).shank(eInd).inds     = spks_struct(eInd).spikes.inds + (chunkOff(k).*30000);
                        featStruct.seg(k).shank(eInd).id       = zeros(numel(spks_struct(eInd).spikes.inds),1);
                        featStruct.seg(k).shank(eInd).ts       = round(featStruct.seg(k).shank(eInd).inds./30)';
                        featStruct.seg(k).shank(eInd).params   = spks_struct(eInd).features.params;

                        if k==1 && eInd==1
                            [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(1,arrayType);
                            featStruct.paramNames = spks_struct(eInd).features.paramNames;
                            featStruct.rSp = rSp;
                            featStruct.cSp = cSp;
                            featStruct.electrodeMatrix = electrodeMatrix;

                            sessionStruct.resolution = spks_struct(eInd).spikes.resolution;
                            sessionStruct.winL       = spks_struct(eInd).spikes.winL;
                            sessionStruct.winR       = spks_struct(eInd).spikes.winR;
                        end

                        sessionStruct.seg(k).shank(eInd).wfs        = spks_struct(eInd).spikes.wfs;
                        sessionStruct.seg(k).shank(eInd).inds       = spks_struct(eInd).spikes.inds;                    
                    end
                    
                disp('_________________________________________________________________________');
                disp(' ');
                disp(' ');                    
                disp(' ');

        end
        
    end

%% SAVE FEATURE STRUCTURE

% clear paramVals featStruct wfStruct newValues
offset                  = 0;
featStruct.chunk        = chunkOff;
featStruct.arrayType    = arrayType;

disp('_________________________________________________________________________');
disp(['save ~/' targetName '_ft featStruct']); 
save([targetName '_ft.mat'],'featStruct');

%% SAVE SESSION STRUCTURE
 
disp(['save ~/' targetName '_ss sessionStruct']);
save([targetName '_ss.mat'],'sessionStruct');

%% CLOSE DOWN THE PARALLEL POOL

% delete(gcp);
disp('_________________________________________________________________________');

%% DEPRECATED CODE
% Ns5DATA.MetaTags.DataDurationSec = 1000;    
% 
% % determine how to calculate segments
%     switch numel(chunk)
%         case 1
%             numSegs = ceil(Ns5DATA.MetaTags.DataDurationSec ./ chunk);
%             timeOffset = 0;
%             for k=1:numSegs
%                 seg_txt{k} = ['t:' num2str((k-1)*chunk) ':' num2str(k*chunk)]
%             end
%             
%         case 3
%             numSegs = ceil( (chunk(2)-chunk(1)) ./ chunk(3) );
%             timeOffset = chunk(1);
%             for k=1:numSegs
%                 seg_txt{k} = ['t:' num2str(((k-1)*chunk(3))+chunk(1)) ':' num2str((k*chunk(3))+chunk(1))]
%             end
%             
%         otherwise
%             numSegs = size(chunk,1);
%             timeOffset = 0;
%             disp('***NOTE***');
%             disp(' ');
%             disp('Assuming 0 timeOffset because user will know how to correct offsets for subsequent data processing');
%             disp(' ');
%             disp('***NOTE***');
%             for k=1:numSegs
%                 seg_txt{k} = ['t:' num2str(chunk(k,1)) ':' num2str(chunk(k,2))]
%             end
%     end
%     disp([num2str(numSegs) ' segments to be processed.']);
%     disp(' ');
%     disp(' ');

%                             disp(['Filtering data on channel ' num2str(electrodeMatrix(i,j)) '|' num2str(find(allRecChan==electrodeMatrix(i,j))) ' ...']);

                        %=============================
                        % USE A WAVELET FILTERING METHOD FOR SPEED
%                             wname = 'db4'; 
%                             maxlevel = 6; % 7.5kHz low pass
%                             [c,l] = wavedec(rawData, maxlevel, wname);
%                             c = wthcoef('a', c, l);
%                             spkData(i,:) = sgolayfilt(waverec(c, l, wname),3,9);


% BUILD THE WF STRUCTURE 
% 
% for j=1:size(sessionStruct.seg(1).shank,2)
%         
%     for i=1:sessionStruct.numSegs
%         
%         wfdata.seg(i).wfs(numel(sessionStruct.seg(i).shank(j).spikes.wfs)).values = zeros(8,61);
% 
%         for k=1:numel(sessionStruct.seg(i).shank(j).spikes.wfs)
%             wfdata.seg(i).wfs(k).values = int16(sessionStruct.seg(i).shank(j).spikes.wfs(k).values);
%         end
%     end
%     
%     % For getting waveform data 
%     wfdata.chunk    = chunk;
%     wfdata.shank    = j;
%     wfdata.filename = fileName;
%     wfdata.datatype = 'int16';
%     
%     disp(['save ' fileName(1:totChar-4) '_shank' num2str(j) '_wf wfdata']);
%     eval(['save ' fileName(1:totChar-4) '_shank' num2str(j) '_wf wfdata']);
% 
%     clear wfdata    
%     
% end


% 
%                         % DISPLAY PHYS AND BEHAVIOR
%                         if dispOn==1
%                             disp('Displaying data for the current channel...');
%                             if i==1
%                                 figure(1); clf;
%                             end
% 
%                             figure(1);
% 
%                             subplot(10,1,1:8);
%                             if rem(j,2)==0
%     %                             plot(spkData(i,:)+((((j-1).*8)+(i-1)).*1000),'Color',[0 0.67 1]); hold on; axis tight;
%                                 plot([events.inds ; events.inds],[spkData(i,events.inds).*0 ; spkData(i,events.inds)]+(((j-1).*8)+(i-1)).*1000,'Color',[0 0.67 1]); hold on; axis off; axis tight;
%                             else
%     %                             plot(spkData(i,:)+((((j-1).*8)+(i-1)).*1000),'Color',[1 0 0]); hold on; axis tight;
%                                 plot([events.inds ; events.inds],[spkData(i,events.inds).*0 ; spkData(i,events.inds)]+(((j-1).*8)+(i-1)).*1000,'Color',[1 0    0]); hold on; axis off; axis tight;
%                             end
% 
%                             subplot(10,1,9);
%                             plot(ContData.behavior.sLeverV((chunk.*1000.*(k-1))+1:chunk.*1000.*k),'k');
%                             axis off; axis tight;
%                             subplot(10,1,10);
%                             plot(ContData.behavior.rawLick((chunk.*1000.*(k-1))+1:chunk.*1000.*k),'k');
%                             axis off; axis tight;
%                                 
%                             drawnow;
%                         end
                    
%                         %=============================
%                         % DISPLAY PHYS AND BEHAVIOR
%                         if dispOn==1
%                             disp('Displaying data for the current channel...');
%                             if i==1
%                                 figure(1); clf;
%                             end
% 
%                             figure(1);
% 
%                             subplot(10,1,1:size(electrodeMatrix,1));
%                             if rem(j,2)==0
%                                 plot([events.inds ; events.inds],[spkData(i,events.inds).*0 ; spkData(i,events.inds)]+(((j-1).*size(electrodeMatrix,1))+(i-1)).*1000,'Color',[0 0.67 1]); hold on; axis off; axis tight;
%                             else
%                                 plot([events.inds ; events.inds],[spkData(i,events.inds).*0 ; spkData(i,events.inds)]+(((j-1).*size(electrodeMatrix,1))+(i-1)).*1000,'Color',[1 0    0]); hold on; axis off; axis tight;
%                             end
% 
%                             if behaviorTrue
%                                 subplot(10,1,9);
%                                 plot(ContData.behavior.sLeverV((chunk.*1000.*(k-1))+1:chunk.*1000.*k),'k');
%                                 axis off; axis tight;
%                                 subplot(10,1,10);
%                                 plot(ContData.behavior.rawLick((chunk.*1000.*(k-1))+1:chunk.*1000.*k),'k');
%                                 axis off; axis tight;
%                             end
%                             
%                             drawnow;
%                         end


