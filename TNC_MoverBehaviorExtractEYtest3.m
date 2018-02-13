%function [reach, ContData] = TNC_MoverBehaviorExtractEY(varargin)%,targetName)
func =0;
%func =1;

plotOn=0;

%% HouseKeeping
if func
filenamestr=varargin{1};
if nargin ==2
    g=varargin{2}; plotOn=1;
    if varargin{2}==0
        g=999; plotOn=0;
    end 
end
else plotOn=1
end
% if you use flexistick, use moveavg(X,30) to smooth out
    
date=str2num(filenamestr([1:2,4:5,7:8]));
newJoy=1;
if date<130920
   newJoy=0;
end

    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Beginning file: ' filenamestr(1,1:length(filenamestr)-3)]);
    disp(' ');disp(' ');  

    dispOn = 0;
    clear('reach'); clear('ContData'); clear('currX');
    
%% LOAD EVENT DATA

    Ns4DATA = openNSx('read', filenamestr);%took out 'report';  [ filenamestr, '.ns2']);
 stimTime=[]; rewTime=[];

% Reward
    reward=Ns4DATA.Data(5,:);
rewTime1=[reward;1:length(reward)];
rewTimee=rewTime1(2,rewTime1(1,:)>4500);
rewTime=rewTimee(1);
for i=2:length(rewTimee) 
    if (rewTimee(i)-rewTimee(i-1)>450)
        rewTime=[rewTime,rewTimee(i)];
    end
end
ContData.behavior.threshInds=rewTime';
ContData.behavior.rewardInds = ContData.behavior.threshInds;
    
if length(ContData.behavior.rewardInds)> 20
    numTrials = numel(ContData.behavior.threshInds);
        
        licks = abs(Ns4DATA.Data(1,:)-mean(Ns4DATA.Data(1,:)));
lickThresh=sort(licks);
lickTime1=find(licks>lickThresh(end-4)/5);
lickTimes=lickTime1(1);
for k = 2:length(lickTime1)
   if lickTime1(k)-lickTimes(end)>80 % = to 15 Hz.  way too fast, but reasonable-ish
       lickTimes=[lickTimes, lickTime1(k)];
   end
end
ContData.behavior.lickTimes=lickTimes;

     
     
     valRew=[];  %skipped for now
     justRewardedTrials=0; %just the reaches immediately before reward, or all except ITI reaches
     if justRewardedTrials
    for thisReward2 = 1:length(ContData.behavior.rewardInds)
        thisReward = ContData.behavior.rewardInds(thisReward2);
        if thisReward-3500>0
        valRew=[valRew,thisReward-3500:(thisReward)];
        else
             valRew=[valRew,1:thisReward];
        end
    end   
     else notvalRew=[];
         ttt=1:length(Ns4DATA.Data);
         for thisReward2 = 1:length(ContData.behavior.rewardInds)
        thisReward = ContData.behavior.rewardInds(thisReward2);
        if thisReward+1501>length(Ns4DATA.Data)
           notvalRew=[notvalRew,thisReward:length(Ns4DATA.Data)];
        else
             notvalRew=[notvalRew,thisReward:thisReward+1500];
        end   
         end
         valRew=ttt(ismember(ttt, notvalRew)==0);
     end
        
    
    if numTrials>20 & length(valRew)>1000 
     
stim=Ns4DATA.Data(4,:);
stimTime1=[stim;1:length(stim)];
stimTime1=stimTime1(2,stimTime1(1,:)>median(stim)+std(stim)*4);%30000);
if length(stimTime1)>5
stimTime=stimTime1(1);
for i=2:length(stimTime1) 
    if (stimTime1(i)-stimTime1(i-1)>450)
        stimTime=[stimTime,stimTime1(i)];
    end
end
stimCount=length(stimTime);
else stimTime=[];stimCount=0;
end
if length(stimTime)>150
  stimTime=[];
end
ContData.behavior.stim=stimTime;
 lickAnal; %generates ContData.behavior.licks (count, hz, stimmed)


        clear leverData sLeverData tmpLeverData
        
        if any(diff(ContData.behavior.rewardInds)>1.5e5)
            realEnd=ContData.behavior.rewardInds(find(diff(ContData.behavior.rewardInds)>1.5e5,1))+1000;
            Ns4DATA.Data=Ns4DATA.Data(:,1:realEnd);
            ContData.behavior.rewardInds = ContData.behavior.rewardInds(ContData.behavior.rewardInds<realEnd);
        disp('cut off long break > 2.5 minutes')
        end
        
    
        leverData(1,:) = Ns4DATA.Data(3,:)-median(Ns4DATA.Data(3,:));
        leverData(2,:) = Ns4DATA.Data(2,:)-median(Ns4DATA.Data(2,:));
        if var(leverData(2,:))/ var(leverData(1,:))>2
            leverData=[leverData(2,:); leverData(1,:)];
        end   
        if newJoy==0
        %   leverData=[-leverData(2,:); leverData(1,:)];
              %if abs(mean(leverData(1,:)))<5
           %       leverData(1,:)=leverData(1,:)*33.3;
          %        leverData(2,:)=leverData(2,:)*10;%33.3;
       %end
        end
        
          if mean(leverData(1,:) )<0  %flips so that most movements are positive
            leverData(1,:) =-leverData(1,:);
        end
        if mean(leverData(2,:) )<0
            leverData(2,:) =-leverData(2,:);
        end
    
    if sum(leverData(1,:)>std(leverData(1,:))*2)/sum(leverData(2,:)>std(leverData(2,:))*2) >5 
        leverData(2,:)=zeros(1,length(leverData(2,:)));
    end

       % reset voltage 
       startReset=0;
        if  ismember(str2num(filenamestr(10:11)),[10:13])
          startReset=0
        end
          if startReset==1
pp=ContData.behavior.rewardInds(1:end-1); %find reward times
%if newJoy==0
 %   pp=[pp;ContData.behavior.stim(2:end-1)'-3100];
  %  pp=sort(pp)';
%end
reachDelay=3000;
resetOffset1=leverData(1,pp+reachDelay); %what is the abs position at reset time
resetOffset2=leverData(2,pp+reachDelay);
for pp1=1:length(pp)-1;
    leverData(1,pp(pp1):pp(pp1+1)-1)=leverData(1,pp(pp1):pp(pp1+1)-1)-resetOffset1(pp1);
    leverData(2,pp(pp1):pp(pp1+1)-1)=leverData(2,pp(pp1):pp(pp1+1)-1)-resetOffset2(pp1);
end
% leverData(1,pp(pp1+1):end)=leverData(1,pp(pp1+1):end)-resetOffset1(pp1+1); %not really necessary, but good to be thorough
 %leverData(2,pp(pp1+1):end)=leverData(2,pp(pp1+1):end)-resetOffset2(pp1+1); %not really necessary, but good to be thorough
          end
   
     odd1=leverData(1,:)>(mean(leverData(1,:))+4*std(leverData(1,:)));     
     leverData(1,odd1) = mean(leverData(1,:))+4*std(leverData(1,:));
      odd2=leverData(2,:)>(mean(leverData(2,:))+4*std(leverData(2,:)));     
     leverData(2,odd2) = mean(leverData(2,:))+4*std(leverData(2,:));
     
    sLeverData(1,:) = sgolayfilt(leverData(1,:),9,101);
    sLeverData(2,:) = sgolayfilt(leverData(2,:),9,101);

    ContData.behavior.sLeverData = sLeverData;
%%

%% EXTRACT VELOCITY DATA FROM CONT LEVER DATA

    numSamples = size(ContData.behavior.sLeverData,2);
    tmpLeverData(1,:) = ContData.behavior.sLeverData(1,:);%sgolayfilt(ContData.behavior.sLeverData(1,:),3,151);
    tmpLeverData(2,:) = ContData.behavior.sLeverData(2,:);%sgolayfilt(ContData.behavior.sLeverData(2,:),3,151);
%these extra smoothings were NOT HELPFUL

    if newJoy==0
       %tmpLeverData(1,:)  = sgolayfilt(tmpLeverData(1,:) ,3,51);
       tmpLeverData(1,:)  = moveavg(tmpLeverData(1,:) ,151);
    end
    
    sLeverV = zeros(1,numSamples);

   % disp(' ');disp(' ');disp('Extracting velocity...');

%    dX = diff(tmpLeverData(1,:));
 %   dY = diff(tmpLeverData(2,:));

  %  sLeverV = sqrt( dX.^2 + dY.^2 );
   
    %ok, here's the problem, you've been doing euclidean distance of the
    %VELOCITY! 
  posX=tmpLeverData(1,:);
  posY=tmpLeverData(2,:);
  posXY=abs(posX)+abs(posY);
  dX=zeros(1,length(posX));
  dY=zeros(1,length(posY));
for i = find(abs(posX)>200) %X has a more consistent 0 value
    if i >4
      dX(i)=posX(i)+posX(i-1)-(posX(i-2)+posX(i-3));
      dY(i)=posY(i)+posY(i-1)-(posY(i-2)+posY(i-3));
    end
end
sLeverV=dX+dY;%(abs(dX)+abs(dY)).*(dX>0); %dX+dY  %maybe try .*d>0 too
sLeverV=dX;
%sLeverV= sqrt(dX.^2+dY.^2);
    
% this takes out return velocities!!!
    %sLeverVY=dY; %%%%%
    disp(' ');disp(' Complete. ');disp(' ');

    ContData.behavior.sLeverV = sLeverV;
    ContData.behavior.sLeverVm = sgolayfilt(ContData.behavior.sLeverV,3,501);
     %figure(101); plot(sLeverV);
   
    clear sLeverV;
    clear sLeverY2; %%%%%
%% FIND MOVEMENTS
    method = 'vel';
    numC = 5;
    clear progSt* reach

    currX       = ContData.behavior.sLeverData(1,:);
    currY       = ContData.behavior.sLeverData(2,:);
    currV       = ContData.behavior.sLeverV; %rather than Vm 
   % if any(currX<(median(currX)-std(currX)))
    %    disp('sadfasdfsdfsadfasdfasdfasdfasdfasdf')
    %notGood = currX<(median(currX)-std(currX));
    %currX(notGood)=median(currX);
    %currY(notGood)=median(currY);
    %currV(notGood)=median(currV);
    %end
    
    currX2= currX+currY;%currX;currX2(ismember(1:length(currX), valRew)==0)=median(currX);
    currX2=sqrt(currX.^2 + currY.^2);
    
    % currently, only reaches from 1700 to 250 ms before reward are being
    %considered.  see reach selection criteria (valRew) above
    
    pre     = 100;
    post    = 10;       
    minSpace = 250;
    count = 1;
    scc=sort(currX2);
    scc1=[];
   % for k = 2:length(ContData.behavior.rewardInds)  bad idea
    %  scc1=[scc1,currX(ContData.behavior.rewardInds(k)-3000:ContData.behavior.rewardInds(k))];
    %end
    %scc=sort(scc1);
     scc(scc<median(scc))=median(scc);
    if newJoy==10
cutoff = scc(round(length(scc)*.95));
    else
%cutoff=scc(max(find(diff(scc)>2,5)))
%  cutoff=1400;
cutoff= scc(max(find(sgolayfilt(diff(scc),1,1001)>.02,5)))
  if cutoff>5.5e3
      cutoff=5.5e3
  end
    end

  % figure(124); hold on; plot(scc(1:max(find(diff(moveavg(diff(scc),10000))>1e-4,10))),'k')
 %  plot(max(find(diff(moveavg(diff(scc),10000))>1e-4,10)):length(scc),scc(max(find(diff(moveavg(diff(scc),10000))>1e-4,10)):end),'r')
   
 
    % threshold the velocities
    %allValidSamps   = find(currV(1:numel(currV)-(minSpace.*2))>10); %0.35
    %allValidSamps   = find(currX(1:numel(currX)-(minSpace.*2))>5000);
    if length(cutoff)==0 
       allValidSamps=0
    else allValidSamps   = find(currX2(1:numel(currX2)-(minSpace.*2))>cutoff);
    end 
  if length(allValidSamps)>1000  
allValidSamps = allValidSamps([true,diff(allValidSamps)==1]);

x=allValidSamps;
aaa=[];

 dd=mat2cell(x,1,diff([0,find(diff(x) >3),length(x)]));
 for ii = 1:length(dd)
     if length(dd{ii})>50
         aaa=[aaa,dd{ii}];
     elseif ii>1
         if min(dd{ii})-max(dd{ii-1})<50
             aaa=[aaa,dd{ii}];
         end
     end
 end
 allValidSamps=[1,aaa];
 
 %%
            
  %  allValidSamps   = allValidSamps(ismember(allValidSamps, valRew));
    
    if numel(allValidSamps)>0

        switch method

            case 'vel'

    %             progStartTMP(count,1)  = currStamps(allValidSamps(1));
                progStartTMP(count,2)  = currX(allValidSamps(1));
                progStartTMP(count,3)  = currY(allValidSamps(1));
                progStartTMP(count,4)  = allValidSamps(1);

                for j=2:numel(allValidSamps)

                    if allValidSamps(j)>allValidSamps(j-1)+minSpace
                        
                        postN = find(currV(allValidSamps(j-1):allValidSamps(j))<1.5,1,'first'); %was <1
                        if numel(postN)<1
                            if plotOn
                            figure(213); clf; 
                            plot(progStartTMP(count,4):allValidSamps(j),currV(progStartTMP(count,4):allValidSamps(j)),'k');
                            hold on; plot(progStartTMP(count,4),currV(progStartTMP(count,4)),'bo',allValidSamps(j-1),currV(allValidSamps(j-1)),'ro');
                            pause(0.1);
                            end

                            disp(['Cannot find stop', num2str(j)]);                        
                            postN=post;
                        end
    %                     progStopTMP(count,1)   = currStamps(allValidSamps(j-1)+post);
                        progStopTMP(count,2)   = currX(allValidSamps(j-1)+postN);
                        progStopTMP(count,3)   = currY(allValidSamps(j-1)+postN);
                        progStopTMP(count,4)   = allValidSamps(j-1)+postN;

                        count                  = count+1;

                        preN = find(currV(allValidSamps(j-1) : allValidSamps(j))<1.5,1,'last'); %was <1
                        if numel(preN)<1
                            disp(['Cannot find start',num2str(j)]);
    %                     progStartTMP(count,1)  = currStamps(allValidSamps(j)-pre);
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

               
        durs  = progStopTMP(:,4)-progStartTMP(:,4);
         dwd=durs(durs<nanmean(durs)+1.5*nanstd(durs));
        lowDur=81;
        if  ismember(str2num(filenamestr(10:11)),[10:13])
            lowDur=61;
        end
        goodDurs = durs<(nanmean(durs)+1.5*nanstd(durs)) & durs > 81;%lowDur;  % high was 1200, 100
        % goodDurs = durs<(nanmean(durs)+1.5*nanstd(durs)) & durs >(nanmean(dwd)-.85*nanstd(dwd));
         progStartTMP = progStartTMP (goodDurs,:);
         progStopTMP = progStopTMP (goodDurs,:);
         
%% NUTS N BOLTS
                count = 1;
                for k = 1:size(progStartTMP,1)
                    if k==1
                        reach.init = 1;
                    end

                    % reaches must be at least 50 ms long
                    if progStopTMP(k,4)-progStartTMP(k,4)>=50 & progStartTMP(k,4)>minSpace
%RewStop=int64(ContData.behavior.rewardInds)-int64(progStopTMP(k,4));  %is the reach during the reward period
 %                       if any(RewStop<-200 & RewStop>-1500)==0 % it was <0 and  >-1500
 RewStart=int64(ContData.behavior.rewardInds)-int64(progStartTMP(k,4));  %is the reach during the reward period
                        if any(RewStart<200 & RewStart>-1500)==0 
                        trajAngle   = atan2(progStopTMP(k,3)-progStartTMP(k,3),progStopTMP(k,2)-progStartTMP(k,2));

                        if (pdist2([progStopTMP(k,2),progStopTMP(k,3)],[mean(currX),mean(currY)]) > pdist2([progStartTMP(k,2),progStartTMP(k,3)],[mean(currX),mean(currY)]))
                            reach.out(count) = 1;
                        else
                            reach.out(count) = 0;
                        end
                        velTraj = ContData.behavior.sLeverV(progStartTMP(k,4) : progStopTMP(k,4));
                        xVals = ContData.behavior.sLeverData(1,progStartTMP(k,4) : progStopTMP(k,4))-ContData.behavior.sLeverData(1,progStartTMP(k,4));
                        yVals = ContData.behavior.sLeverData(2,progStartTMP(k,4) : progStopTMP(k,4))-ContData.behavior.sLeverData(2,progStartTMP(k,4));
                        xyVals = sqrt(xVals.^2 +yVals.^2);
                       % figure(k+50); plot(xVals,yVals)
                        %reach traj insert
                      %  currTraj = ([xVals;yVals])';
                    %  trajAnal;
                        
                         reach.start(count,:)  = progStartTMP(k,:);
                        reach.stop(count,:)   = progStopTMP(k,:);
                        reach.angle(count,1)  = trajAngle;
                        reach.dist(count,1)   = trapz(xyVals); %velTraj
                        reach.dist(count,2)   = pdist2(progStartTMP(k,2:3) , progStopTMP(k,2:3));
                        reach.disp(count,1)   = max(xyVals);
                        reach.disp(count,2)   = find(max(xyVals)==xyVals,1)+progStartTMP(k,4);
                        %  figure(222); hold on; plot(xVals+count*7000);

                        tmp = findpeaks(velTraj);
                        reach.numpks(count,1) = numel(tmp.loc);
                        reach.dur(count,1)    = progStopTMP(k,4) - progStartTMP(k,4);
                        reach.vel(count,1)   = max(velTraj);
                        reach.vel(count,2)   = trapz(velTraj) ./ reach.dur(count,1);
                        reach.vel(count,3)   = var(velTraj);
                        reach.vel(count,4)   = find(velTraj==max(velTraj),1);
                        reach.vel(count,5)   = find(velTraj==max(velTraj),1)+reach.start(count,4);
                        
                        Acc=diff(velTraj);
                        for ac= 2:length(Acc)-1
                          if abs(Acc(ac)-Acc(ac-1))>(std(abs(Acc))*4) %big movement
                              if abs(Acc(ac)-Acc(ac+1))>(std(abs(Acc))*4) %subsequent big movement
                                  if (Acc(ac-1)-Acc(ac))*(Acc(ac)-Acc(ac+1))<0 %subsequent movements in opposite directions
                                    Acc(ac-1:ac+1)=0;
                                  end; end;end;
                        end
                               
                        reach.acc(count,1)   = max(Acc);
                        reach.acc(count,2)   = nanmean(Acc);
                        reach.acc(count,3)   = max(Acc(1:50)); % max in first 50 ms of movement
if reach.acc(count,1)<.5
                 %       reach.acc(count,:)=nan;
end
                        reach.acc(count,4)   = find(Acc==max(Acc),1);
                        
                        reach.tort(count,1)  = reach.dist(count,1) ./ pdist2([progStopTMP(k,2),progStopTMP(k,3)],[progStartTMP(k,2),progStartTMP(k,3)]);

                 
                        
                        %valThrInd = find(ContData.behavior.threshInds>reach.start(count,4) & ContData.behavior.threshInds<reach.stop(count,4));
                        %if numel(valThrInd)>0
                        %    reach.rewarded(count)= 1;
                       % else
                       %     reach.rewarded(count)= 0;                       
                      %  end

                        count                 = count+1;
                        else progStartTMP(count,4);
                        %  badStart=[badStart; progStartTMP(k,4), RewStart(RewStart<200 & RewStart>-1500)]
                        end; %end 
                        end;
                end
%%

        end
        
        if count>10
%remove outliers
goodVel = [reach.vel((reach.start(2:end,4)-reach.stop(1:end-1,4))>200,1); reach.vel(end,1)]; %long enough ITI
%goodVel=goodVel(goodVel(:,1)>15);  %way too small of velocities
%if sum(reach.disp(:,1)<1.5e4) <10
 %  goodVel=
%end
goodVel=remOut(goodVel(:,1));
goodDisp=reach.disp(:,1)>(median(reach.disp(:,1))/4);
goodReach=ismember(reach.vel(:,1),goodVel) & ismember(reach.vel(:,1),goodVel) ;

reach.vel=reach.vel(goodReach,:);
reach.start=reach.start(goodReach,:);
reach.stop=reach.stop(goodReach,:);
reach.dur=reach.dur(goodReach,:);
reach.dist=reach.dist(goodReach,:);
reach.acc=reach.acc(goodReach,:);
reach.tort=reach.tort(goodReach,:);
reach.disp=reach.disp(goodReach,:);

velTrace=[];
thisColor=1;
for jj = reach.vel(:,5)'
    traceDur=[(jj-1000):(jj+2000)];
    if jj-1000>0
    if jj+2000<length(currX2)
 % plot(currX2(traceDur),'Color',color(thisColor,:))
  velTrace(thisColor,:)=currX2(traceDur);
  %stimVel(thisColor,:)=currV(traceDur);
  thisColor=thisColor+1;
    end;end
end


ymean = []; xmean = [];
sampleRate = [];
startstop=[reach.start(:,4),progStopTMP(ismember(progStartTMP(:,4), reach.start(:,4)),4)];
%trajAnal2

disp(['Num valid reaches: ' num2str(length(reach.start))]);

        level = ones(1,size(reach.start,1)).*0.35;
        if plotOn
        figure(g); clf;
      %  plot(ContData.behavior.sLeverV,'g'); hold on;
     %   plot((leverData(1,:)-median(leverData(1,:)))/300, 'k'); hold on
  plot(currX2/300,'k'); hold on
  plot(currV/100)
  plot(reach.disp(:,2),reach.disp(:,1)/300,'m*')
        plot(reach.start(:,4),level,'r^','LineWidth',2);
        plot(reach.stop(:,4),level,'bo', 'LineWidth',2);
        plot([0,length(leverData)],[cutoff, cutoff]/300, 'm')
        plot(ContData.behavior.rewardInds,4.5,'g*')
        if length(stimTime)>0
        plot(stimTime,5,'b*')
        end
      plot(reach.vel(:,5),reach.vel(:,1)/100,'c*'); %was 100
title(filenamestr)
 end
        [vals,inds] = sort(reach.dur,'ascend');    

        reach.numReaches = size(reach.start,1);

        dims = floor(sqrt(reach.numReaches));
        scaler = 3000;
       
        
        if dispOn
            
           % figure(125); clf;
        
            for p = 1:1:dims.^2

                l = inds(p);

                offset = l.*ones(1,numel(reach.start(l,4):reach.stop(l,4)));
                [x,y] = ind2sub(dims,p);

               % figure(115);
                %plot((x.*scaler),(y.*scaler),'k.','MarkerSize',8);
                %hold on;
                %plot(currX(reach.start(l,4):reach.stop(l,4))-currX(reach.start(l,4))+(x.*scaler),currY(reach.start(l,4):reach.stop(l,4))-currY(reach.start(l,4))+(y.*scaler) , 'k','LineWidth',1,'Color',[p/reach.numReaches 0.67-(p/reach.numReaches)*0.67 1-(p/reach.numReaches)]);
               % if i==dims.^2
                %    axis tight; axis off;
                %end

                startInd = round(reach.start(l,4)./2);
                stopInd = round(reach.stop(l,4)./2);

            end
        end
        
        if reach.numReaches>555
                figure(117); clf;
                subplot(311);
                p = polyfit(reach.dist(:,1),reach.vel(:,1),1);
                plot(reach.dist(:,1),reach.vel(:,1),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
                xlabel('Total distance');
                ylabel('Peak Velocity');

                subplot(312);
                p = polyfit(reach.dist(:,1),reach.dur,2);
                plot(reach.dist(:,1),reach.dur,'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
                xlabel('Total distance');
                ylabel('Movement Duration (ms)');

                subplot(313);
                p = polyfit(reach.dist(:,1),reach.acc(:,1),1);
                plot(reach.dist(:,1),reach.acc(:,1),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
                xlabel('Total distance');
                ylabel('Peak acceleration');

                %figure(8); clf;
                %rose(reach.angle);

                drawnow;
        end
    else
        reach.numReaches=0; %if count
        end
    else
        
        reach.numReaches=0;
    
    end
    
reach.cutoff=cutoff;
    
%% DIMENSION REDUCTION OF REACH PARAMETERS
if 1==2
if reach.numReaches > 5
    clear matFor*

    matForDR = [reach.vel reach.acc reach.tort reach.dur reach.dist];

    for j=1:size(matForDR,2)
        matForCov(:,j) = (matForDR(:,j) - mean(matForDR(:,j))) ./ std(matForDR(:,j)) ;
    end

   % figure(121);
   % imagesc(corr(matForCov),[0 1]);
    %map = colormap(pink);
    %colormap(1-map);

    [eV,eD] = eig(cov(matForCov));
  %  figure(112); plot(cumsum(diag(eD))./sum(diag(eD)),'ko-');

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
end %if 1==2 line 366
%% WAS THE REACH REWARDED?

if reach.numReaches >0
    for p=1:reach.numReaches  % not particularly good
     tmp=reach.vel(p,5)-int32(ContData.behavior.rewardInds);
      %  tmp = find(ContData.behavior.threshInds>reach.start(p,4) & ContData.behavior.threshInds<reach.stop(p,4));
        if any(tmp<-500 & tmp>-1300)
            reach.rewarded(p)=1;
        else
            reach.rewarded(p)=0;
        end
        
    end 
%% STIM COMPOSITION
stimReach=[];NOstimReach=[];
GoodStim=0;  stimVel=[];
tt=[];
   stimTrace=nan(1,4001);
if length(ContData.behavior.stim)>5
for kk= reach.vel(:,5)'  
    if sum(abs(kk-ContData.behavior.stim)<190)>0
     %   min(abs(kk-ContData.behavior.stim))
    kk1=kk;
    tt=[tt,min(abs(kk-ContData.behavior.stim))];
    if length(kk1)>1
        kk1= kk1(max(kk1));
    end
    stimReach=[stimReach, kk1];
    else NOstimReach=[NOstimReach,kk];
    end
end

color=jet(length(ContData.behavior.stim));
%figure(300+g); hold on; 
thisColor=1;
    stimTrace=zeros(length(ContData.behavior.stim),4001);
for jj = ContData.behavior.stim
    traceDur=[(jj-2000):(jj+2000)];
    if jj-2000>0
    if jj+2000<length(currX2)
 % plot(currX2(traceDur),'Color',color(thisColor,:))
  stimTrace(thisColor,:)=currX2(traceDur);
  stimVel(thisColor,:)=currV(traceDur);
  thisColor=thisColor+1;
    else putIn=currX2(jj-2000:end);
         stimTrace(thisColor,1:length(putIn))=putIn;
  thisColor=thisColor+1;
    end
   else putIn=currX2(1:jj+2000);
         stimTrace(thisColor,1:length(putIn))=putIn;
  thisColor=thisColor+1;
    end
end
stimTrace=stimTrace(:,1:4001);
stimVel=stimVel(:,1:2100);
%semline(stimTrace,'k');
%[maxVal maxInd] = max(stimTrace(:,400:end)');
%figure;hist(maxInd)

stimmed = ismember(reach.vel(:,5), stimReach);
if sum(stimmed)>5
    disp('stimstimstimstimstimstimstimstimstim')
    %% BLAH
ns=[0; stimmed(1:end-1)];
ns=ns==1;
reach.stim.start=reach.start(ismember(reach.vel(:,5), stimReach),:);
reach.stim.stop=reach.stop(ismember(reach.vel(:,5), stimReach),:);
reach.stim.vel=reach.vel(ismember(reach.vel(:,5), stimReach),:);
reach.stim.dur=reach.dur(ismember(reach.vel(:,5), stimReach),:);
reach.stim.dist=reach.dist(ismember(reach.vel(:,5), stimReach),:);
reach.stim.acc=reach.acc(ismember(reach.vel(:,5), stimReach),:);
reach.stim.tort=reach.tort(ismember(reach.vel(:,5), stimReach),:);
reach.stim.disp=reach.disp(ismember(reach.vel(:,5), stimReach),:);


reach.nextstim.start=reach.start(ns,:);
reach.nextstim.stop=reach.stop(ns,:);
reach.nextstim.vel=reach.vel(ns,:);
reach.nextstim.dur=reach.dur(ns,:);
reach.nextstim.dist=reach.dist(ns,:);
reach.nextstim.acc=reach.acc(ns,:);
reach.nextstim.disp=reach.disp(ns,:);
reach.nextstim.tort=reach.tort(ns,:);


reach.NOstim.start=reach.start(ismember(reach.vel(:,5), NOstimReach),:);
reach.NOstim.stop=reach.stop(ismember(reach.vel(:,5), NOstimReach),:);
reach.NOstim.vel=reach.vel(ismember(reach.vel(:,5), NOstimReach),:);
reach.NOstim.dur=reach.dur(ismember(reach.vel(:,5), NOstimReach),:);
reach.NOstim.dist=reach.dist(ismember(reach.vel(:,5), NOstimReach),:);
reach.NOstim.acc=reach.acc(ismember(reach.vel(:,5), NOstimReach),:);
reach.NOstim.disp=reach.disp(ismember(reach.vel(:,5), NOstimReach),:);
reach.NOstim.tort=reach.tort(ismember(reach.vel(:,5), NOstimReach),:);
GoodStim=1
end
end    
    if GoodStim==0
     reach.stim.stop=nan(1,5);
      reach.stim.start=nan(1,5);
    reach.stim.vel=nan(1,5);
reach.stim.dur=nan(1,1);
reach.stim.dist=nan(1,5);
reach.stim.acc=nan(1,5);
reach.stim.tort=nan(1,5);
reach.stim.disp=nan(1,5);
%reach.stim.trace=nan(1,5);

     reach.nextstim.stop=nan(1,5);
      reach.nextstim.start=nan(1,5);
    reach.nextstim.vel=nan(1,5);
reach.nextstim.dur=nan(1,1);
reach.nextstim.dist=nan(1,5);
reach.nextstim.acc=nan(1,5);
reach.nextstim.disp=nan(1,5);
reach.nextstim.tort=nan(1,5);

reach.NOstim.stop=nan(1,5);
reach.NOstim.start=nan(1,5);
 reach.NOstim.vel=nan(1,5);
reach.NOstim.dur=nan(1,1);
reach.NOstim.dist=nan(1,5);
reach.NOstim.acc=nan(1,5);
reach.NOstim.disp=nan(1,5);
reach.NOstim.tort=nan(1,5);
    end

%% WRITE THE REACH STRUCTURE OF THE CONTDATA STRUCTURE
%figure(g+100); hold off;
%plot(reach.vel(:,5), reach.vel(:,1)/9,'ko'); hold on
%plot(reach.stim.vel(:,5), reach.stim.vel(:,1),'bo');
%title(filenamestr)
%lsline;

thisColor=1;
    stimTrace=zeros(length(reach.start(:,4)),4001);
for jj = reach.start(:,4)'%ContData.behavior.stim
    traceDur=[(jj-2000):(jj+2000)];
    if jj-2000>0
    if jj+2000<length(currX2)
 % plot(currX2(traceDur),'Color',color(thisColor,:))
  stimTrace(thisColor,:)=currX2(traceDur);
  thisColor=thisColor+1;
    else putIn=currX2(jj-2000:end);
         stimTrace(thisColor,1:length(putIn))=putIn;
  thisColor=thisColor+1;
    end
   else putIn=currX2(1:jj+2000);
         stimTrace(thisColor,1:length(putIn))=putIn;
  thisColor=thisColor+1;
    end
end
stimTrace=stimTrace(:,1:4001);


reach.st=tt';
reach.stim.trace=stimTrace;
reach.velTrace=velTrace;
reach.stimVel=stimVel;
    ContData.behavior.reach = reach;
    ContData.data=ContData.behavior;
 
%% COMPLETED ALL ANALYSIS. SAVE AND START OVER

    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Completed file: ' filenamestr(1,1:length(filenamestr)-3)]);
    
    % save the data from this file
   % disp(['saved as ' targetName '_bh.mat']);
    %save([targetName '_bh.mat'],'ContData');
    
    disp('%-------------------------------------------------------------------');
    disp(' ');
    else  disp('%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$');
    disp(['NO VALID REACHES : ',filenamestr]);
    reach.numReaches=0; ContData=[];
  end
  else  disp('%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$');
    disp(['NO VALID REACHES : ',filenamestr]);
    reach.numReaches=0; ContData=[];
  end
     else  disp('%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$');
    disp(['NO VALID REACHES : ',filenamestr]);
    reach.numReaches=0; ContData=[];
  end 
   notnow=0;
   if notnow==1
    ll=length(reach.start(:,4));
    trajInd=zeros(ll,2500);
    trajAlign=zeros(ll,2500+1000);
    color=jet(ll);
    for i = 1:ll
        figure(234); hold on;
        thisR=currX(reach.start(i,4):reach.stop(i,4))-currX(reach.start(i,4));
        trajInd(i, 1:length(thisR))=thisR;
        plot(thisR, 'Color', color(i,:))
       figure(100+g); hold on;
         align=find(thisR>2000,1);
         align=find(thisR> max(thisR)*.75,1);
           if align>0
       plot(1001-align:length(thisR)-align+1000,thisR, 'Color', color(i,:)) 
%        trajAlign(i,(1001-align:length(thisR)-align+1000))=thisR;
           end
    end
   %plot(
    %plot(1-1000:length(trajAlign)-1000,mean(trajAlign),'k','LineWidth',3)
    %% comments
    %could align reach to start, slop > threshold, or peak
    %velocity(winner)...
    %surprised reach start is not at the start of reaches / is somewhat
    %noisy
    
   end
    
else
    reach.numReaches=0; ContData=[];
end
if 1==2
ld=13:1.6e5;
count=1;
figure;
for kk = ld
ldd=count/length(ld);
plot(currX(kk),currY(kk),'.', 'Color', [1-ldd 0 ldd]);
count=count+1; hold on; end
hold off
end


