%filenamestr= files(2).name;%'15_05_12m28NPostRecSham.ns2'
fv = openNSx('read', filenamestr);%took out 'report';  [ filenamestr, '.ns2']);

%disp(['Align to ', num2str(find(fv(4,:)>1500)/1000) , 'sec on movie'])
disp('Need a light trigger here to align!')
%movie starts on "file name record", not initial file save
% don't interact with movie menu bar until you're done

fva=(sqrt(fv.Data(2,:)/4000) + sqrt(fv.Data(3,:)/4000)).^2;
%fva=fv.Data(3,:)/4000;
fva1=sgolayfilt(fva-mean(fva(1:100)),3,251);
fvv=[0 0 0];
for i = 4:length(fva1)
    fvv=[fvv, fva1(i)+fva1(i-1)-(fva1(i-2)+fva1(i-3))];
end
fva2=sgolayfilt(fv.Data(2,:)-mean(fv.Data(2,1:100)),3,251)/4000;
fva3=sgolayfilt(fv.Data(3,:)-mean(fv.Data(3,1:100)),3,251)/4000;

licks=fv.Data(1,:);
lickThresh=std(licks)*4;
lickTime1=find(licks>lickThresh);
lickTimes1=lickTime1(1);
for k = 2:length(lickTime1)
   if lickTime1(k)-lickTimes1(end)>80 % = to 15 Hz.  way too fast, but reasonable-ish
       lickTimes1=[lickTimes1, lickTime1(k)];
   end
end
lightOn=abs(fv.Data(4,:))>50;
lightOn1=1;
lightOn1=find(lightOn==1,1);
lightOnAll=find(lightOn==1)-lightOn1;

fva1=fva1(lightOn1:end);
fvv1=sgolayfilt(fvv(lightOn1:end),3,151);
lickTimes1=lickTimes1(lickTimes1>lightOn1);

xwid=3000;
lims=stimLim;%3e4:5.2e4;%10.4e4:12.5e4;%2.5e4;
stimmer=[];
for gg= 1:length(stimIndex)
 stimmer=[stimmer, stimIndex(gg):stimIndex(gg)+450];
end

%xwidF=500; xwidB=2500; lims=lims-xwidF+1;%for the bar
step=20; 
%lims=135001:167000; %135 is approx zero 

figure; hold on
fva=fva1(lims);
fvv=fvv1(lims);
lickTimes=lickTimes1(ismember(lickTimes1,lims))-lims(1);
%fva=(fva-min(fva))/(max(fva)-min(fva));
%fvv=(fvv-min(fvv))/(max(fvv)-min(fvv));
fva=fva-mode(round2(fva,-2));
fva=fva/max(fva);
fvv=fvv/max(fvv);
%plot(zeros(1,length(fva)), 'k--');%'Color', [.2 .2 .2]); hold on
%plot(ones(1,length(fva))*2, 'k--')%,'Color', [0 .6 .2]); hold on
plot(fva,'k', 'LineWidth',2); hold on
plot(fvv+2,'g', 'LineWidth',2); 
%plot(lickTimes,1.01,'bo')
lightOnAll=lightOnAll(ismember(lightOnAll,lims))-lims(1);
plot(stimmer,3.1,'b.')
%plot(lightOnAll,3.1,'b.')
%%

x= 1:length(fva);
%ymax=max(fva)*1.25;
ymax=max(fva)+max(fvv)+1*1.25;
ymin=min(fva)-abs(min(fva)*.40);
ylim([ymin ymax])
nframe=1;

 axis off
%for n=xwidF:step:(x(end)-xwidB) %for the bar
    for n=1:step:(x(end)-xwid)  %for no bar
    %plot(fva(n:(n+xwid-1)),'k');
   % ylim([ymin ymax])
   % xlim([1 xwid])
   %  axis off
   
    xlim([x(n) x(n+xwid-1)])  %for no bar
   %for the bar
   %plot(n , ymin,'g.', 'LineWidth',12)
   %xlim([x(n-xwidF+1) x(n+xwidB)])
          MM(nframe)=getframe;
          nframe=nframe+1;         
end

%%
fps=1000/step;
movie2avi(MM, ['low1newER_',filenamestr(1:end-4), '.avi'],'FPS',fps)
%movie2avi(MM, [filenamestr(1:end-4), '2.avi'],'FPS',fps)
%close; close;
disp('must open avi and resave to work')

