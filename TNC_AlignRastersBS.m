function [response] = TNC_AlignRastersBS(delta,spkStamps,alignStamps,window,rasterFlag,boxcar)
% FUNCTION DETAILS: Using a set of timestamps provided in ALIGNSTAMPS creates a matrix of DATA vectors that span -WINDOW(1,1) to +WINDOW(1,2) around each timestamp.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% 

alignStamps = round(alignStamps);

if size(alignStamps,1) > size(alignStamps,2)
    alignStamps = alignStamps';
end

count = 1;
numStamps = size(alignStamps,2);

% avgRate     = mean(delta(alignStamps(1):alignStamps(numStamps)));
if alignStamps(numStamps) > numel(delta)
    stdRate     = std(delta(alignStamps(1):numel(delta)));     
    meanRate     = mean(delta(alignStamps(1):numel(delta)));     
else
    stdRate     = std(delta(alignStamps(1):alignStamps(numStamps))); 
    meanRate     = mean(delta(alignStamps(1):numel(delta)));     
end

response.image.aligned = zeros(numStamps, window(1,1)+window(1,2)+1);

for i = 1:numStamps

    currStamp = round(alignStamps(1,i));

    % get the list of timestamps within the current window
    if rasterFlag
        validSpks = find(spkStamps>currStamp-window(1,1) & spkStamps<currStamp+window(1,2));
        response.raster.trial(i).ts = spkStamps(validSpks) - currStamp;        
    end
    
    % case 'image' (use all indices by padding)
    if currStamp + window(1,2) < size(delta,2)
        if currStamp-window(1,1) > 0
            response.image.aligned(count,:)   = delta(1,currStamp-window(1,1):currStamp+window(1,2));
            response.image.indicesUsed(count) = i;
        else
%             response.image.aligned(count,:)   = [zeros(1,abs(currStamp-window(1,1)-1)),delta(1,1:currStamp+window(1,2))];
%             response.image.indicesUsed(count) = i;
            count = count-1;
        end
    else
        if currStamp-window(1,1) > 0
            tmp = delta(1,currStamp-window(1,1):length(delta));
            response.image.aligned(count,1:length(tmp))   = tmp;
            response.image.indicesUsed(count) = i;
        end
    end

    count = count+1;
        
end

% write out the response variables...
count = count-1;

% size(response.image.aligned,2)
    
if boxcar == 1
    centers = (0:6:size(response.image.aligned,2)-6) + 3;
    boxcarTMP = sum(response.image.aligned,1)./numStamps;
    for i = 1:length(centers)
        response.image.boxcar(i) = sum(boxcarTMP(centers(i)-2:centers(i)+3));
    end   
end

if stableTrials==-1
    response.image.alignedtime    = -window(1,1):1:window(1,2);
    response.image.psthAVG        = mean(response.image.aligned,1);% - meanRate;
    response.image.psthSEM        = std( response.image.aligned,0,1) ./ sqrt(size(response.image.aligned,1)-1);
else
    response.image.alignedtime    = -window(1,1):1:window(1,2);
    response.image.psthAVG        = mean(response.image.aligned(stableTrials:count,:),1);% - meanRate;
    response.image.psthSEM        = std( response.image.aligned(stableTrials:count,:),0,1) ./ sqrt(size(response.image.aligned(stableTrials:count,:),1)-1);
end

% use the entire left window as an estimate of the mean
avgRate = mean(response.image.psthAVG(1:round(window(1,1)./2)));
if stdRate>0.0001
    response.image.psthZ          = (response.image.psthAVG-avgRate) ./ stdRate;  
    response.image.psthZe         = response.image.psthSEM ./ stdRate;  
    response.noZ                  = 0;
else
    disp(['STD=' num2str( mean(response.image.psthAVG(1,1:window(1,1)),2) ) ' | cannot z-score']);
    response.image.psthZ          = response.image.psthAVG;
    response.image.psthZe         = response.image.psthAVG;
    response.noZ                  = 1;
end

response.meanRate = meanRate;
response.stdRate = stdRate;
