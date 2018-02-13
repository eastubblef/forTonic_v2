function [forHist, histVals, tss] = alignSpikeBS(pop, filenamestr, varargin, threshInds)

%% BS modified from Yttri's mfile to use after TNC_ConvertTSDtoPopDataBS 7/15/15
% This is the final step for aligning TSs of neural data to particular behavioral event (i.e. laser pulse on)
% Called from MoverScript.m

% INPUTS: pop = TS from individual neurons' raw data (cluster 1,2,3)... formerly from PopData structure
%         filenamestr = behavioral .ns4 file into blackrock (Yttri). BS will make this her laser pulse

% OUTPUTS: is a stim x ts aligned matrix for rasterizing, y vals from aligned histogram, and timestamps
%          tss is a matrix of TSs extracted from the pop field of PopData, formerly

%%
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
% stim=NsDATA.Data(4,:);                                                    %This is Yttri's laser (he's trying to put out times>500ms)
                                                                            %BS gets laser input from ContData.behavior.threshInds
stimTime=[];                                                                
stim = double(threshInds)';
stimTime1=[stim;1:length(stim)];                                            %Figure this out: create a matrix in which the first row is my laser on times, 2nd row is 1-47 (to number the indices of stim?)
if sum(stim>2000)==0                                                        %if the sum of (any time is > 2000, it will be a one; if not, it will be zero)... for my purposes, the else will have to appy since this isn't true
   stimTime1=stimTime1(:,stimTime1(1,:)< 3.5);                              %Which can't happen, given above.
else
    stimTime1=stimTime1(:,stimTime1(1,:)>2000);                             %Else, This will have to happen. Which means nothing changes.
end

% for i=2:length(stimTime1)                                                   %What is this loop doing? for i=the 2nd element-end of my laser on times
%     if (stimTime1(2,i)-stimTime1(2,i-1)>500)                                %if the second row of those times - the i before it (which will = 1 for me) each time through the loop (i.e. will never be > 500)
%        stimTime=[stimTime,stimTime1(2,i)];                                  %Then the new matrix = [47] ? - But it can never happen; therefore matrix is empty
%     end
% end

%BS attempt for this to make sense:
for i=2:length(stimTime1)                                                   %What is this loop doing? for i=the 2nd element-end of my laser on times
    if (stimTime1(1,i)-stimTime1(1,i-1)>900);                               %There should be 1 sec bt. stim times; use the first row; not sure why Yttri used the iterative row
        stimTime=[stimTime,stimTime1(1,i)];                                 %Then the new matrix = [2-last index of laser times]
    end
end

%% unpack popData structure
tss=nan(length(units), 50000);  %probably too much, equivalent to 1 hour @ 25hz
for thisunitID = 1:length(units);
    thisunit=units(thisunitID);
    thists=pop.session.unit(thisunit).ts;
    tss(thisunitID,1:length(thists))=thists;
end
[ms,ns]=find(isfinite(tss));
tss=tss(:,1:max(ns));

histWidthBack=1000;
histWidthFront=200;
binWidth=1;
xbins=[-histWidthFront:binWidth:histWidthBack];

eachUnit=0
if eachUnit==1
    figure;
    dimm=ceil(sqrt(length(units)));
    for k = 1:length(units);
        tss1=tss(k,:);
        forHist=nan(length(stimTime),2*histWidthBack);
        stimNum=1;
        for stims = stimTime
            trigSpikes=tss1(ismember(tss1,(stims-histWidthFront):(stims+histWidthBack)))-stims;
            if length(trigSpikes)>2 %typically around 10 or so?
                forHist(stimNum,1:length(trigSpikes))=trigSpikes;
            end
            stimNum=stimNum+1;
        end
        [m,n]=find(isfinite(forHist));
        forHist=forHist(:,1:max(n));
        
        [~, theRank]=sort(sum(isfinite(forHist')));
        forHist=forHist(fliplr(theRank),:);
        
        subplot(dimm,dimm,k);
        rasterBS(forHist);
        xlim([-histWidthFront, histWidthBack])
        ylim([1 length(stimTime(:,1))])
    end
end

   %nbins=(histWidthBack+histWidthFront)/binWidth; %40; %go for 5 ms bins
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
title(filenamestr);
ylabel('Trials')

subplot(212);

hz=1000/binWidth;
summer=(max(ms)*max(m));
histVals=hist(forHist(forHist<(histWidthBack+1)),xbins)*hz/summer; 
plot(xbins,histVals)
%bar(histVals,[-histWidthFront:binWidth:(histWidthBack-binWidth)]) %need to get X values correct if you bar
hold on
% plot([0,0],[0,round2(max(histVals),1)+5],'r', 'LineWidth', 3);            %Not sure what round2 is about
plot([0,0],[0,round(max(histVals),1)+5],'r', 'LineWidth', 3);

%ylim=[0,round2(max(histVals),1)+5];
xlabel=('ms');
ylabel('Aligned Pop Firing Rate')
hold off
