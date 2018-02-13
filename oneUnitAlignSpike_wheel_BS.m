

function [forHist, histVals, tss] = oneUnitAlignSpike_wheel_BS(pop, fname, varargin, varAlign) %this does actually plot one neuron's raster at a time!
%% BS modified from oneUnitAlignSpike2BS to align to solenoid trigger/trial starts/
% Further modified for L/R raw movements 2.24.16

%% BS modified from Yttri's mfile to use after TNC_ConvertTSDtoPopDataBS 7/15/15
% One unit at a time is plotted on a raster and psth - multiple per figure 8/2/15
% This is the final step for aligning TSs of neural data to particular behavioral event (i.e. laser pulse on)
% Use this after  MoverScript.m (Already obtained the PopData and ContData structures

% INPUTS: pop = TS from individual neurons' raw data (cluster 1,2,3)... from "PopData" structure (ex: 150529_6008_pop.mat)
%         fname = file that contains PopData neural data structure
%         behavFile = behavioral .ns4 file extracted from blackrock as "ContData" structure 
%         threshInds = ContData.behavior.threshInds contains TS for blackrock channel to which to align neural data (example, laser pulse-in channel) 

% OUTPUTS: is a stim x ts aligned matrix for rasterizing, y vals from aligned histogram, and timestamps
%          tss is a matrix of TSs extracted from the pop field of PopData (tss = time stamp spike)


%%
% Inputs: 
% laserOn = 20; %(For ea. 1Hz stim, laser on duration was _ms)
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/sorted4/151119ipsi/1001'
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/laserChunks/untagged';

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/151119good/behaveChunks/tagged';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151104/1001/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151106/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/160505';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better/laserChunks/untagged';
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/161005';
% behavFile = strcat(fpath, '/', '160505_all_bh.mat');
% fname = strcat(fpath, '/', '160505_all_pop.mat');
% fname = strcat(fpath, '/', '151119_all_pop.mat');
% behavFile = strcat(fpath, '/', '151119_all_bh.mat');
% fname = strcat(fpath, '/', '151118_all_pop.mat');
% behavFile = strcat(fpath, '/', '151118_all_bh.mat');
% behavFile = strcat(fpath, '/', '151105_untagged_bh.mat');
% fname = strcat(fpath, '/', '151105_untagged_pop.mat');
behavFile = strcat(fpath, '/', '161005_bh.mat');
fname = strcat(fpath, '/', '161005_pop.mat');

% neuronType = 'untagged';
% neuronType = 'tag';
neuronType = 'tagged';

load(fname);
pop = PopData;
units = 1:length(pop.session.unit); %use this as the default
% units = pop.session.unit(4:5);

varargin = units;

load(behavFile)
behavior = ContData.behavior;

% validsR_first = behavior.validsR_first';
% validsL_first = behavior.validsL_first';
% rewardTS = behavior.rewTrialStarts;
LrewTrialStarts = behavior.LrewTrialStarts;
RrewTrialStarts = behavior.RrewTrialStarts;

LincorrectTrialStarts = behavior.LincorrectTrialStarts;
RincorrectTrialStarts = behavior.RincorrectTrialStarts;

%  wheelsLfirst = behavior.wheelsLfirst;
%  wheelsRfirst = behavior.wheelsRfirst;

 wheelsLfirstValid = behavior.wheelsLfirstValid;
 wheelsRfirstValid = behavior.wheelsRfirstValid;
 wheelsILfirstValid = behavior.wheelsILfirstValid;
 wheelsIRfirstValid = behavior.wheelsIRfirstValid;
 rewTimesValid = behavior.rewTimesValid;

% varAlign = validsR_first;
% varAlign = validsL_first;

varAlign = wheelsLfirstValid;
% varAlign = wheelsRfirstValid;
% varAlign = wheelsILfirstValid;
% varAlign = wheelsIRfirstValid;
% varAlign = rewTimesValid;
% varAlign = LrewTrialStarts;
% varAlign = RrewTrialStarts;
%  varAlign = LincorrectTrialStarts;
% varAlign = RincorrectTrialStarts;

%% For alignment to solenoid valve ons or trial starts
%Look for all L or R movements from ContData.behavior structure:

stimTime = [];  
stim = varAlign;

stimTime1 = [stim; 1:length(stim)];                                         %Figure out what Yttri did: create a matrix in which the first row is my laser on times, 2nd row is 1-47 (to number the indices of stim?)
if sum(stim > 2000) == 0                                                        %if the sum of (any time is > 2000, it will be a one; if not, it will be zero)... for my purposes, the else will have to appy since this isn't true
   stimTime1 = stimTime1(:, stimTime1(1,:) < 3.5);                              %which can't happen, given above.
else
   stimTime1 = stimTime1(:, stimTime1(1,:) > 2000);                            %else, this will have to happen. Which means nothing changes.
end


for i = 2:length(stimTime1)                                                     %for i=the 2nd element-end of my sol on times
%     if (stimTime1(1, i) - stimTime1(1, i-1) > 200);                             %there should be (100?) ms bt. stim times for wheel movement(?); use the first row; not sure why Yttri used the iterative row
    if (stimTime1(1, i) - stimTime1(1, i-1) > 200);
        stimTime = [stimTime, stimTime1(1, i)];                                 %then the new matrix = [2-last index of laser times]
    end
end

%% Unpack the units field of the popData structure

tss = nan(length(units), 50000);  %likely high, equivalent to 1 hour @ 25hz    Stands for time stamp spike; %Run : find(diff(tss<0))
for thisunitID = 1:length(units);
    thisunit = units(thisunitID);
    thists = pop.session.unit(thisunit).ts;
    tss(thisunitID, 1:length(thists)) = thists;
end
[ms,ns] = find(isfinite(tss));
tss = tss(:, 1:max(ns));

%% Plot one unit on its own raster: stimTime = vector of laser/validsL_first TS; tss = matrix in which ea. row is a neuron's TS.

% histWidthBack = 700
% histWidthFront = 700;

histWidthBack = 2000    
histWidthFront = 1000;

binWidth = 1;
xbins = [-histWidthFront:binWidth:histWidthBack];

% Plot each unit on its own raster:
eachUnit = 0;
% if eachUnit==1                                                                                %He's made this loop nonconsequential
if eachUnit == 0
figure; hold on;
    dimm=ceil(sqrt(length(units))); %
    for k = 1:length(units);                                                                    %set up a loop to hit ea. of 32 units
        tss1=tss(k,:);                                                                          %ea. row through the iterative process
        forHist=nan(length(stimTime),2*histWidthBack);                                          %forHist = nans that have num rows = length(stimTime) and 10000 cols 
        stimNum=1;
        for stims = stimTime                                                                    %stims is new Var for stimTime
            trigSpikes=tss1(ismember(tss1,(stims-histWidthFront):(stims+histWidthBack)))-stims; %trigSpikes = ea. units' TS(members of that row's TS, col = laser pulse TS-200ms:TS+200)-laserTS
            if length(trigSpikes)>2                                                             %typically around 10 or so?
                forHist(stimNum,1:length(trigSpikes))=trigSpikes;                               %forHist(1,col=1:length(trigSpikes)) = trigSpikes
            end
            stimNum=stimNum+1;                                                                  %ea. time thru loop, add another row to forHist = trigSpikes
        end
        [m,n]=find(isfinite(forHist));                                                          %matrix [m,n] = indices in forHist that are finite
        forHist=forHist(:,1:max(n));                                                            %forHist = build the matrix
        
        [~, theRank]=sort(sum(isfinite(forHist')));                                             %theRank = explanatory
        forHist=forHist(fliplr(theRank),:);
        
        subplot(dimm,dimm,k);
        
        [data] = rasterBSwheel(forHist,k);       %122 x 4 = num of trials x      %Updated 10.21.16 for k input
%         newData = find(isfinite(data)); %Call to raster - return 'data' for plotting
%         [x,y] = size(data);
        xlim([-histWidthFront, histWidthBack]);
%         ylim([0 numel(threshInds]);        %since ThreshInds are the TSs at zero - nope; too much.
        
        ylabel('Trials');
        xlabel('ms');
        hold on;
%         l = plot([0,0],[0, length(stimTime(1,:))],'bl', 'LineWidth', 0.5);
%         l = plot([0,0],[0, stimTime(:,1)],'bl', 'LineWidth', 0.5);

%      fnameNew = '151105tag_pop.mat';
%         newfname1 = fnameNew(end-17:end);
%         newfname1 = fnameNew(end-21:end-7);

%         newfnamepre = newfname1(end-17:end-11);
%         newfnamepost = newfname1(end-3:end);
%         newfnamepre = newfname1;

% % %         newfnamemid = newfname1(end-11:end-8);
%         for unitNum = 1:length(k)
%             unitNum = num2str(k);
%             unitname = strcat(neuronType,' u#', unitNum);
% %             newfname = strcat(newfnamepre, '-', unitname);
% %             title(newfname);
%             title(unitname);
%         end
    end
end

%%
% % nbins=(histWidthBack+histWidthFront)/binWidth; %40; %go for 5 ms bins - yttri's
% forHist=nan(length(stimTime),2*histWidthBack);                                              
% stimNum=1;
% for stims = stimTime                                    
%     trigSpikes=tss(ismember(tss,(stims-histWidthFront):(stims+histWidthBack)))-stims;
%     if length(trigSpikes)>2 %typically around 10 or so?
%       forHist(stimNum,1:length(trigSpikes))=trigSpikes;
%     end
%     stimNum=stimNum+1;
% end
% [m,n]=find(isfinite(forHist));
% forHist=forHist(:,1:max(n));                            %Ea row of forHist is a single unit
% 
% [~, theRank]=sort(sum(isfinite(forHist')));             %Sort along cols such that nans are excluded and vals are summed
% forHist=forHist(fliplr(theRank),:);

