

function [forHist, histVals, tss] = allUnitsAlignSpikeBS(pop, fname, varargin, threshInds)% No, this does actually plot one neuron's raster at a time!

%% BS modified from Yttri's mfile to use after TNC_ConvertTSDtoPopDataBS 7/15/15
% This is the final step for aligning TSs of neural data to particular behavioral event (i.e. laser pulse on)
% Use this after  MoverScript.m (Already obtained the PopData and ContData structures

% INPUTS: pop = TS from individual neurons' raw data (cluster 1,2,3)... from "PopData" structure (ex: 150529_6008_pop.mat)
%         fname = file that contains PopData neural data structure
%         behavFile = behavioral .ns4 file extracted from blackrock as "ContData" structure 
%         threshInds = ContData.behavior.threshInds contains TS for blackrock channel to which to align neural data (example, laser pulse-in channel) 

% OUTPUTS: is a stim x ts aligned matrix for rasterizing, y vals from aligned histogram, and timestamps
%          tss is a matrix of TSs extracted from the pop field of PopData

%%
% Inputs: 
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Gad2_tag/150529ipsi/5007';
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Gad2_tag/150606contra/3004'; 

cd(fpath);
fname = strcat(fpath, '/', '150606_3004_pop.mat');
% fname = strcat(fpath, '/', '150529_5007_pop.mat');
load(fname);
pop = PopData;
units = 1:length(pop.session.unit);
varargin = units;

behavFile = strcat(fpath, '/', '150606_3004_bh.mat');
% behavFile = strcat(fpath, '/', '150529_5007_bh.mat');
load(behavFile)
behavior = ContData.behavior;
threshInds = behavior.threshInds;

units=1:length(pop.session.unit);

if nargin ==3
   units=varargin;
end

%% BS - get alignment times from laser pulses                               
% Yttri calls source file here: ('fileX'), filenamestr (.ns2) and if you want, which units
% 
% get alignment times
% TNC_ConvertTSDtoPopData(datanamestr,sessNum)                              %could use .ns5 file
% NsDATA = openNSx('report','read', filenamestr);
% stim=NsDATA.Data(4,:);                                                    %This is Yttri's laser (he's trying to put out times > 500ms)
                                                                            %BS gets laser input from ContData.behavior.threshInds
stimTime = [];                                                                
stim = double(threshInds)';
stimTime1 = [stim; 1:length(stim)];                                         %Figure out what Yttri did: create a matrix in which the first row is my laser on times, 2nd row is 1-47 (to number the indices of stim?)
if sum(stim > 2000) == 0                                                        %if the sum of (any time is > 2000, it will be a one; if not, it will be zero)... for my purposes, the else will have to appy since this isn't true
   stimTime1 = stimTime1(:, stimTime1(1,:) < 3.5);                              %which can't happen, given above.
else
   stimTime1 = stimTime1(:, stimTime1(1,:) > 2000);                            %else, this will have to happen. Which means nothing changes.
end


%BS attempt for this to make sense:
for i = 2:length(stimTime1)                                                     %for i=the 2nd element-end of my laser on times
    if (stimTime1(1, i) - stimTime1(1, i-1) > 900);                             %there should be 1 sec bt. stim times; use the first row; not sure why Yttri used the iterative row
        stimTime = [stimTime, stimTime1(1, i)];                                 %then the new matrix = [2-last index of laser times]
    end
end

%% Unpack the units field of the popData structure

tss = nan(length(units), 50000);  %probably too much, equivalent to 1 hour @ 25hz    Stands for time stamp spike; %Run : find(diff(tss<0))
for thisunitID = 1:length(units);
    thisunit = units(thisunitID);
    thists = pop.session.unit(thisunit).ts;
    tss(thisunitID, 1:length(thists)) = thists;
end
[ms,ns] = find(isfinite(tss));
tss = tss(:, 1:max(ns));

% Unpack one shank at a time from popData 7/23/15
% shanks = pop.session.unit.sh; %nope

%% This plots all units on the same raster: stimTime = vector of laser TS; tss = matrix in which ea. row is a neuron's TS.

% histWidthBack=1000;
histWidthBack = 500
histWidthFront = 200;
binWidth = 1;
xbins = [-histWidthFront:binWidth:histWidthBack];

%This plots all units on the same raster:
eachUnit=0
if eachUnit==1
    figure;
    dimm=ceil(sqrt(length(units)));
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
        rasterBS(forHist);
        xlim([-histWidthFront, histWidthBack])
        ylim([1 length(stimTime(:,1))])
    end
end
% 
% nbins=(histWidthBack+histWidthFront)/binWidth; %40; %go for 5 ms bins - yttri's
forHist=nan(length(stimTime),2*histWidthBack);                                              
stimNum=1;
for stims = stimTime                                    
    trigSpikes=tss(ismember(tss,(stims-histWidthFront):(stims+histWidthBack)))-stims;
    if length(trigSpikes)>2 %typically around 10 or so?
      forHist(stimNum,1:length(trigSpikes))=trigSpikes;
    end
    stimNum=stimNum+1;
end
[m,n]=find(isfinite(forHist));
forHist=forHist(:,1:max(n));

[~, theRank]=sort(sum(isfinite(forHist')));
forHist=forHist(fliplr(theRank),:);

%% Do the plotting
figure()
subplot(211); 

rasterBS(forHist);
newfname1 = fname(end-16:end);
newfnamepre = newfname1(end-16:end-13);
newfnamepost = newfname1(end-6:end);
newfnamemid = newfname1(end-11:end-8);
newfname = strcat(newfnamepre, '-', newfnamemid, newfnamepost);
title(newfname);
ylabel('Trials from multiple units');
xlabel('ms');

subplot(212);

hz = 1000/binWidth;
summer = (max(ms)*max(m));
histVals = hist(forHist(forHist<(histWidthBack+1)),xbins) * hz/summer; 
plot(xbins,histVals)
    % bar(histVals,[-histWidthFront:binWidth:(histWidthBack-binWidth)]) %need to get X values correct if you bar
hold on
    % plot([0,0],[0,round2(max(histVals),1)+5],'r', 'LineWidth', 3);            %Not sure what round2 is about
plot([0,0],[0,round(max(histVals),1)+5],'r', 'LineWidth', 3);

    % ylim=[0,round2(max(histVals),1)+5];
xlabel('ms');
ylabel('Aligned Pop Firing Rate')
hold off

