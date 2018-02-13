
function [stimTimeA, stimTimeB] = callsROCpref(pop, fname, varargin, varAlign1, varAlign2)

%% This function extracts the following for inputs into ROC_preference  
% a - number_of_trials x 1 vector of firing rates under Condition A (varAlign1)
% b - number_of_trials x 1 vector of firing rates under Condition B (varAlign2)

fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks';
fname = strcat(fpath, '/', '151118_all_pop.mat');
behavFile = strcat(fpath, '/', '151118_all_bh.mat');

neuronType = 'tag&Untag';

load(fname);
pop = PopData;

units = 1:length(pop.session.unit);
% units = pop.session.unit(1);

varargin = units;

load(behavFile)
behavior = ContData.behavior;

LrewTrialStarts = behavior.LrewTrialStarts;
RrewTrialStarts = behavior.RrewTrialStarts;

LincorrectTrialStarts = behavior.LincorrectTrialStarts;
RincorrectTrialStarts = behavior.RincorrectTrialStarts;
 
wheelsLfirstValid = behavior.wheelsLfirstValid;
wheelsRfirstValid = behavior.wheelsRfirstValid;
wheelsILfirstValid = behavior.wheelsILfirstValid;
wheelsIRfirstValid = behavior.wheelsIRfirstValid;
rewTimesValid = behavior.rewTimesValid;

varAlign1 = LrewTrialStarts;
varAlign2 = RrewTrialStarts;

%% Unpack the timestamps for spikes that correspond to varAlign1:

stimTime = [];  
stim = varAlign1;

stimTime1 = [stim; 1:length(stim)];                                         %Figure out what Yttri did: create a matrix in which the first row is my laser on times, 2nd row is 1-47 (to number the indices of stim?)
if sum(stim > 2000) == 0                                                        %if the sum of (any time is > 2000, it will be a one; if not, it will be zero)... for my purposes, the else will have to appy since this isn't true
   stimTime1 = stimTime1(:, stimTime1(1,:) < 3.5);                              %which can't happen, given above.
else
   stimTime1 = stimTime1(:, stimTime1(1,:) > 2000);                            %else, this will have to happen. Which means nothing changes.
end

for i = 2:length(stimTime1)                                                     %for i=the 2nd element-end of my sol on times
%     if (stimTime1(1, i) - stimTime1(1, i-1) > 200);                           %there should be (100?) ms bt. stim times for wheel movement(?); use the first row; not sure why Yttri used the iterative row
    if (stimTime1(1, i) - stimTime1(1, i-1) > 200);
        stimTime = [stimTime, stimTime1(1, i)];                                 %then the new matrix = [2-last index of laser times]
    end
end
stimTimeA = stimTime;

%% Unpack the units field of the popData structure

tss = nan(length(units), 50000);  %likely high, equivalent to 1 hour @ 25hz    Stands for time stamp spike; %Run : find(diff(tss<0))
for thisunitID = 1:length(units);
    thisunit = units(thisunitID);
    thists = pop.session.unit(thisunit).ts;
    tss(thisunitID, 1:length(thists)) = thists;
end
[ms,ns] = find(isfinite(tss));
tss1 = tss(:, 1:max(ns));

histWidthBack = 2000    
histWidthFront = 1000;

binWidth = 1;
xbins = [-histWidthFront:binWidth:histWidthBack];
dimm=ceil(sqrt(length(units))); %
for k = 1:length(units);                                                                    %set up a loop to hit ea. of 32 units
    tssNew1=tss1(k,:);                                                                      %ea. row through the iterative process
    forHist1=nan(length(stimTimeA),2*histWidthBack);                                         %forHist = nans that have num rows = length(stimTime) and 10000 cols
    stimNum=1;
    for stims = stimTimeA                                                                    %stims is new Var for stimTime
        trigSpikes=tssNew1(ismember(tssNew1,(stims-histWidthFront):(stims+histWidthBack)))-stims; %trigSpikes = ea. units' TS(members of that row's TS, col = laser pulse TS-200ms:TS+200)-laserTS
        if length(trigSpikes)>2                                                             %typically around 10 or so?
            forHist1(stimNum,1:length(trigSpikes))=trigSpikes;                              %forHist(1,col=1:length(trigSpikes)) = trigSpikes
        end
        stimNum=stimNum+1;                                                                  %ea. time thru loop, add another row to forHist = trigSpikes
    end
    [m,n]=find(isfinite(forHist1));                                                         %matrix [m,n] = indices in forHist that are finite
    forHist1=forHist1(:,1:max(n));                                                          %forHist = build the matrix
    
    [~, theRank]=sort(sum(isfinite(forHist1')));                                            %theRank = explanatory
    forHist1=forHist1(fliplr(theRank),:);
    
%     [data1] = rasterBS(forHist1);       %122 x 4 = num of trials x
    %         newData = find(isfinite(data)); %Call to raster - return 'data' for plotting
    %         [x,y] = size(data);
    %         xlim([-histWidthFront, histWidthBack]);
    %         ylabel('Trials');
    %         xlabel('ms');
    %         hold on;
end

%% Unpack the timestamps for spikes that correspond to varAlign2:

stimTime = [];
stim = varAlign2;

stimTime1 = [stim; 1:length(stim)];                                         %Figure out what Yttri did: create a matrix in which the first row is my laser on times, 2nd row is 1-47 (to number the indices of stim?)
if sum(stim > 2000) == 0                                                        %if the sum of (any time is > 2000, it will be a one; if not, it will be zero)... for my purposes, the else will have to appy since this isn't true
    stimTime1 = stimTime1(:, stimTime1(1,:) < 3.5);                              %which can't happen, given above.
else
    stimTime1 = stimTime1(:, stimTime1(1,:) > 2000);                            %else, this will have to happen. Which means nothing changes.
end

for i = 2:length(stimTime1)                                                     %for i=the 2nd element-end of my sol on times
    %     if (stimTime1(1, i) - stimTime1(1, i-1) > 200);                             %there should be (100?) ms bt. stim times for wheel movement(?); use the first row; not sure why Yttri used the iterative row
    if (stimTime1(1, i) - stimTime1(1, i-1) > 200);
        stimTimeB = [stimTime, stimTime1(1, i)];                                 %then the new matrix = [2-last index of laser times]
    end
end

%% Unpack the units field of the popData structure

tss = nan(length(units), 50000);  %likely high, equivalent to 1 hour @ 25hz    Stands for time stamp spike; %Run : find(diff(tss<0))
for thisunitID = 1:length(units);
    thisunit = units(thisunitID);
    thists = pop.session.unit(thisunit).ts;
    tss(thisunitID, 1:length(thists)) = thists;
end
[os,ps] = find(isfinite(tss));
tss2 = tss(:, 1:max(ps));

histWidthBack = 2000
histWidthFront = 1000;

binWidth = 1;
xbins = [-histWidthFront:binWidth:histWidthBack];
dimm=ceil(sqrt(length(units))); %
for k = 1:length(units);                                                                    %set up a loop to hit ea. of 32 units
    tssNew2=tss2(k,:);                                                                          %ea. row through the iterative process
    forHist2=nan(length(stimTimeB),2*histWidthBack);                                          %forHist = nans that have num rows = length(stimTime) and 10000 cols
    stimNum=1;
    for stims = stimTimeB                                                                    %stims is new Var for stimTime
        trigSpikes=tssNew2(ismember(tssNew2,(stims-histWidthFront):(stims+histWidthBack)))-stims; %trigSpikes = ea. units' TS(members of that row's TS, col = laser pulse TS-200ms:TS+200)-laserTS
        if length(trigSpikes)>2                                                             %typically around 10 or so?
            forHist2(stimNum,1:length(trigSpikes))=trigSpikes;                               %forHist(1,col=1:length(trigSpikes)) = trigSpikes
        end
        stimNum=stimNum+1;                                                                  %ea. time thru loop, add another row to forHist = trigSpikes
    end
    [m,n]=find(isfinite(forHist2));                                                          %matrix [m,n] = indices in forHist that are finite
    forHist2=forHist2(:,1:max(n));                                                            %forHist = build the matrix
    
    [~, theRank]=sort(sum(isfinite(forHist2')));                                             %theRank = explanatory
    forHist2=forHist2(fliplr(theRank),:);
    
%     [data2] = rasterBS(forHist2);       %122 x 4 = num of trials x
    %         newData = find(isfinite(data)); %Call to raster - return 'data' for plotting
    %         [x,y] = size(data);
    %         xlim([-histWidthFront, histWidthBack]);
    %         ylabel('Trials');
    %         xlabel('ms');
    %         hold on;
    
end

%% Call the ROC_preference m file:
% INPUTS:  
% a - number_of_trials x 1 vector of firing rates under Condition A
% b - number_of_trials x 1 vector of firing rates under Condition B
% num_repeats - number of times to permute a and b and recalculate preference (used for calculating the p_val associates with pref)


% a = forHist1;  must be a #trials x 1 vector for condition a
% b = forHist2;  must be a #trials x 1 vector for condition b
a = stimTimeA;
b = stimTimeB;
num_repeats = 500;

[pref, p_val] = ROC_preference(a, b, num_repeats)


