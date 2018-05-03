%Updated in Sep. 2017 for alignments to older tagging/recording data:
%% Must set the timescale offset at line 60ish per session!

% Updated 2.20.18 to show individual units' psth w/in their own figure

%% Load the file of interest:
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170112/shanks3&4/4thPass/behaveSegs/behaveNewUnits';
% csv_data = dlmread('mVgatfive_2017_01_12_165817pXY.csv',',',2,0);
% filestr = '/170112002.ns4';

fpath = '/Volumes/My Passport for Mac/170111/newBehaveUnits';
csv_data = dlmread('mVgatfive_2017_01_11_174943pXY.csv',',',2,0);
filestr = '/170111002.ns4';

% fpath = '/Volumes/My Passport for Mac/Vgatfive/170118/behaveNewUnits';
% csv_data = dlmread('mVgatfive_2017_01_18_163328pXY.csv',',',2,0);
% filestr = '/170118002.ns4';

filestr2 = strcat(fpath, '/', filestr);
[data] = openNSx(filestr2);              %MUST LOAD THE openNSx.m file from 2016 for this to work

% S = load('170112newBehave_pop');
S = load('170111newBehave_pop');
% S = load('170118_newBehave_pop.mat');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PopData = S.PopData;
clear S;

%plot to calculate x_off:
chan.start = 129;                                    %trial start channel (Ain 1)
trialStrtsChan = find(data.MetaTags.ChannelID==chan.start);
trialStrtsChan1 = double(data.Data(trialStrtsChan,:)); 
trialStrts = decimate(trialStrtsChan1,10);           %convert these to ms (/10)
figure; hold on;
plot(trialStrts,'k');

%%
Sigma = 24;
t = 1:1:Sigma.*30;
Gaussian = ( 1./( sqrt(2.*pi.*Sigma.*Sigma) ) ) .* exp( -(t-Sigma.*15).^2 ./ (2.*Sigma).^2 );
integral = trapz(Gaussian);                         %trapz(y) returns the approximate integral of y
Gaussian = Gaussian./integral;
[mapName] = TNC_CreateRBColormap(8,'mbr');          %calls color map function

%%
% figure(3); hold off;
figure; hold off;
% x_off = 107700;   % previously thought for 170112, but next line better
% x_off = 107745;   %this is the offset of the ns4 file when the noise reflects csv start 170112
% x_off = 110510;   %this is the offset of the ns4 file when the noise reflects csv start 170118
x_off = 78935;    %this is the offset of the ns4 file when the noise reflects csv start 170111
% x_off = 128212;   %this is the offset of the ns4 file when the noise reflects csv start 160505

x_vals = (csv_data(:,1)-csv_data(1,1))+x_off; %subtract the x-offset (time) from the csv data's timestamps
d_x_vals = [0 diff(x_vals')];                 %the diff of timestamps from csv file
smth_wheel_v = conv( abs(csv_data(:,2)) , [zeros(1,4) ones(1,5) 0] , 'same' )'; %The convolution of two vecs u (velocity) & v (row of ones padded by 0s), will be the overlap under the points as v slides across u

plot(x_vals , csv_data(:,2) ,'LineWidth',2);  %blue; plots col2 (raw velocity) w/ the time offset calculated from the ns4 file
hold on; 
plot(x_vals , smth_wheel_v ,'LineWidth',2);   %orange; plots the smoothed wheel velocity data (all positive-going); 
plot(x_vals, d_x_vals);                       %yellow; plots the diff of timestamps
 
% plot(x_vals ,[0 diff(-csv_data(:,3)')]./2,'LineWidth',2); 
% plot(x_vals ,csv_data(:,3),'LineWidth',2);
% plot(sgolayfilt((decimate(double(data.Data(3,:)),10)./25),3,101));
% plot(decimate(double(data.Data(1,:)),10)./1e2);

%% find intra movement periods
threshold = 5;
intra_mvmt = find(smth_wheel_v>threshold & d_x_vals<30);  %this is a "meaningful" movement: find the overlapping instances for smoothed velocity > 5 and the diff(timestamps) < 30 ms

% figure(4); hold off;
figure; hold off;

plot(x_vals , csv_data(:,2) ,'LineWidth',1);              %plots csv data in blue; inter-mvment-interval in red (fig. 4)
hold on;
plot(x_vals(intra_mvmt) , csv_data(intra_mvmt,2) ,'.');   %plot the intra-mvment periods on top of the raw velocity (fig 4)

%% separate out into valid movements
mvmt_win = [10,25];
starts = find([0 diff(intra_mvmt)]>10);                                           %specifically finds the gaps between movement vectors & makes the next value = start

starts_valid = find(intra_mvmt(starts)>mvmt_win(1) & intra_mvmt(starts)<(numel(csv_data(:,2))-mvmt_win(2)));

starts_clean = starts(starts_valid);
plot(x_vals(intra_mvmt(starts_clean)) , csv_data(intra_mvmt(starts_clean),2) ,'o');%plots in yellow (still fig 4)

%wheel_moves creates #rows of ea. mvmnt bout's vel: 
%when trial starts & bar = -pos, correct rightward mvment has -vel; 
%when +position, correct leftward mv has + vel.
[wheel_moves] = TNC_ExtTrigWins(csv_data(:,2)',intra_mvmt(starts_clean),mvmt_win); %sink.wins = csv(:,2), indexed by the starts_clean with an expanded window; ea. row shows the indexed vals from 10 backwards/25 ahead of the meaningful index
%170112 session = L hemi; Thus, +vel = L mvment & ipsi

for i=1:size(wheel_moves.wins,1)
    magnitude(i) = trapz(wheel_moves.wins(i,mvmt_win(1):20));  %trapz(y) returns the approximate integral of y
end                                                            %gives a per col. readout of how big ea. movement was
magnitude = magnitude';                 %EAS: find the mvments that had more than 3ish vel readouts (of >10-30, ea) should be ~50 or greater
trueMvsPre = abs(magnitude) > 50;
trueMvsInds = find(trueMvsPre == 1);    %EAS: these are the indices of magnitude that are likely meaningful movements

trueMvsInds_L = find(magnitude < -50);  %EAS: indices of magnitude for large/fast left mvments
trueMvsInds_R = find(magnitude > 50);   %EAS: inds of mag for large/fast right mvments

[vals, inds] = sort(magnitude);         %EAS: this does nothing for me
% figure(5);
figure;
imagesc(wheel_moves.wins(inds,:),[-50 50]); colormap(mapName);           %x is the movement vectors sorted by velocity magnitude 
xlabel('vector window'); ylabel('movement vectors, sorted by vel mag');  %should be separating ipsi v contra (one is above, other below).

%% Relevant for export:
export_mvmt_data_for_psths.times        = x_vals(intra_mvmt(starts_clean))'; %timestamps of the movement vector starts
export_mvmt_data_for_psths.magnitude    = magnitude;
export_mvmt_data_for_psths.mag_sort_i   = inds;

export_mvmt_data_for_psths.Lefts        = trueMvsInds_L                      %magnitude inds for L mvments
export_mvmt_data_for_psths.Rights       = trueMvsInds_R                      %magnitude inds for R mvments
export_mvmt_data_for_psths.LeftsTimes   = export_mvmt_data_for_psths.times(trueMvsInds_L); %timestamps of L mvments
export_mvmt_data_for_psths.RightsTimes  = export_mvmt_data_for_psths.times(trueMvsInds_R); %timestamps of R mvments

%% examine some psths associated with movements
numUnits = numel( PopData.session(1).unit );
% figure(10); clf;
figure; 
for j=1:numUnits  
    
    tmp = PopData.session(1).unit(j).ts;
    delta = zeros(1,ceil(max(tmp)));
    delta(round(tmp)) = 1;

    tmpSmooth = conv(delta,Gaussian,'same');
%     figure(10);
%     figure;  %figure 6: not sure what this is plotting on the left side of the subplot
%     subplot(1,3,[1 2]);
%     plot(tmpSmooth + j*0.1); hold on;
    
%     subplot(1,3,3);
    figure; hold on;
    times = export_mvmt_data_for_psths.times;    %These are the aligned mvmnt times
    Ltimes = export_mvmt_data_for_psths.LeftsTimes;
    Rtimes = export_mvmt_data_for_psths.RightsTimes;
    
    valid_valid = find(times>9.2e4 & times<2.07e6);
    [vals, inds] = sort(magnitude(valid_valid));

    psthWin = [1.0e3,2e3];
    
    %     [sink_tmp]   = TNC_ExtTrigWins(tmpSmooth,times(valid_valid),[750,1000]);
    [sink_tmp]   = TNC_ExtTrigWinsEAS(tmpSmooth,times(valid_valid),[750,1000]);
    if ~isfinite(sink_tmp.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
        j = j+1;
    else
    unitNum = num2str(j);
    unitname = strcat(' u#', unitNum);
   
    %Plot the psth and label the unit plotted: %[-win +win], so mvment onset = ~700
    plot((sink_tmp.avg - mean(sink_tmp.avg(1:500))) + (ones(1,numel(sink_tmp.avg)).*j.*0.01), 'k'); hold on;  
    drawnow;  %plots the movement-aligned related neural activity 
    title(unitname); 
    hold off;
    
    %% Plot the L-movements:
    figure; hold on;
    [sink_tmpL]   = TNC_ExtTrigWinsEAS(tmpSmooth,Ltimes,[1000,2000]); %[-win +win], so mvment onset = ~700
    if ~isfinite(sink_tmpL.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
        j = j+1;
    else
    unitNum = num2str(j);
    unitname = strcat(' u#', unitNum);
   
    %Plot the L-mv psth and label the unit plotted
%     plot((sink_tmpL.avg - mean(sink_tmpL.avg(1:500))) + (ones(1,numel(sink_tmpL.avg)).*j.*0.01), 'bl'); hold on;  
    alignVar1 = Ltimes;
    PopData.currParams.stableTrials = -1;                                         % subtracts the first "trial" when calling AlignRaster

    [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session.unit(j).ts, PopData.currParams.stableTrials, alignVar1, psthWin, 1, 1); 
    shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'bl');  %blue for left; for 151105 dataset = ipsi
    title(unitname); 
    end
    
    %% Plot the R-movements:
    [sink_tmpR]   = TNC_ExtTrigWinsEAS(tmpSmooth,Rtimes,[1000,2000]); %[-win +win], so mvment onset = ~700
    if ~isfinite(sink_tmpR.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
        j = j+1;
    else
    unitNum = num2str(j);
    unitname = strcat(' u#', unitNum);
   
    alignVar2 = Rtimes;
    [respCSS2] = TNC_AlignRasters(tmpSmooth,PopData.session.unit(j).ts, PopData.currParams.stableTrials, alignVar2, psthWin, 1, 1); 
    shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthZ, respCSS2.image.psthZe,'r');  %blue for left; for 151105 dataset = ipsi

    %Plot the L-mv psth and label the unit plotted
%     plot((sink_tmpR.avg - mean(sink_tmpL.avg(1:500))) + (ones(1,numel(sink_tmpR.avg)).*j.*0.01), 'r'); hold on;  
    title(unitname); 
    hold off;
    end
%%
%     if j==5
% %         figure(11);
%         figure;
%         imagesc(sink_tmp.wins(inds,:));
%         
% %         figure(12); hold off;
%         figure; hold off;
%         plot(tmpSmooth.*250);
%         hold on;
%         plot(x_vals , csv_data(:,2) ,'LineWidth',2);        
%         
%     end
    
    
    end
    
end
return

%%