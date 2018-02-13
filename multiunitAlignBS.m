function [forHist, histVals, tss] = multiunitAlignBS( filenamestr , ftFile)

%% Yttri's mfile 7/15/15

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
NsDATA = openNSx('report','read', filenamestr);
stimTime=[];
stim=NsDATA.Data(4,:);
stimTime1=[stim;1:length(stim)];
if sum(stim>2000)==0
   stimTime1=stimTime1(:,stimTime1(1,:)< 3.5);
else
stimTime1=stimTime1(:,stimTime1(1,:)>4600);
end
for i=2:length(stimTime1) 
    if (stimTime1(2,i)-stimTime1(2,i-1)>500)
        stimTime=[stimTime,stimTime1(2,i)];
    end
end

ggg=load(ftFile');
allTs=[];
for ii=1:length(ggg.featStruct.seg)
    for jj=1:length(ggg.featStruct.seg(ii).shank)
        allTs=[allTs;ggg.featStruct.seg(ii).shank(jj).ts];
    end
end
tss=allTs;
[ms,ns]=find(isfinite(tss));
tss=tss(:,1:max(ns));

  histWidthBack=800;
   histWidthFront=200;
   binWidth=1;
   xbins=[-histWidthFront:binWidth:histWidthBack];
   %nbins=(histWidthBack+histWidthFront)/binWidth; %40; %go for 5 ms bins
   forHist=nan(length(stimTime),50000);
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
hz=10/binWidth;
summer=mean(nanmean(forHist));%length(forHist(:,1));%(max(ms)*max(m));
histVals=hist(forHist(forHist<(histWidthBack+1)),xbins); 
plot(xbins,histVals/nanmean(histVals(1:100))*100)
%bar(histVals,[-histWidthFront:binWidth:(histWidthBack-binWidth)]) %need to get X values correct if you bar
hold on
plot([0,0],[0,100],'g', 'LineWidth', 3);
%ylim=[0,round2(max(histVals),1)+5];
xlabel=('ms');
ylabel('Aligned MUA % baseline FR')
hold off


end

