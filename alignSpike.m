function [forHist, histVals, tss] = alignSpike(pop,filenamestr,varagin)
%output is a stim x ts aligned matrix for rasterizing, y vals from
%aligned histogram, and timestamps
    units=1:length(pop.session.unit);

if nargin ==3
   units=varagin;
end
%need to identify source ('fileX'), filenamestr (.ns2) and if you want, which units

%get alignment times
%TNC_ConvertTSDtoPopData(datanamestr,sessNum) % could use .ns5 file
NsDATA = openNSx('report','read', filenamestr);
stimTime=[];
stim=NsDATA.Data(4,:);
stimTime1=[stim;1:length(stim)];
if sum(stim>2000)==0
   stimTime1=stimTime1(:,stimTime1(1,:)< 3.5);
else
stimTime1=stimTime1(:,stimTime1(1,:)>2000);
end
for i=2:length(stimTime1) 
    if (stimTime1(2,i)-stimTime1(2,i-1)>500)
        stimTime=[stimTime,stimTime1(2,i)];
    end
end

%unpack popData structure
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
raster(forHist);
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

figure()
subplot(211); 
raster(forHist);
title(filenamestr);
ylabel('Trials')
subplot(212);
hz=1000/binWidth;
summer=(max(ms)*max(m));
histVals=hist(forHist(forHist<(histWidthBack+1)),xbins)*hz/summer; 
plot(xbins,histVals)
%bar(histVals,[-histWidthFront:binWidth:(histWidthBack-binWidth)]) %need to get X values correct if you bar
hold on
plot([0,0],[0,round2(max(histVals),1)+5],'r', 'LineWidth', 3);
%ylim=[0,round2(max(histVals),1)+5];
xlabel=('ms');
ylabel('Aligned Pop Firing Rate')
hold off
