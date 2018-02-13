%Updated in Sep. 2017 for alignments to older tagging/recording data:

csv_data = dlmread('mVgatfive_2017_01_12_165817pXY.csv',',',2,0);
filestr = '/170112002.ns4';
% [data] = openNSx('/170112002.ns4');
% [data] = openNSx(filestr);
[data] = openNSxBS(filestr);

S = load('170112behave_pop.mat');

% csv_data = dlmread('mVgatfive_2017_01_18_163328pXY.csv',',',2,0);
% filestr = '/170118002.ns4';
% % [data] = openNSx('/170118002.ns4');
% [data] = openNSx(filestr);
% S = load('170118behave_pop.mat');

% csv_data = dlmread('mVgatfive_2017_01_11_174943pXY.csv',',',2,0);
% filestr = '/170111002.ns4';
% % [data] = openNSx('/170111002.ns4');
% [data] = openNSx(filestr);
% S = load('170111_pop.mat');

% csv_data = dlmread('mVgatfourtag_2016_05_05_161450pXY.csv',',',2,0);  %R hemisphere recorded
% filestr = '/160505001.ns4';
% [data] = openNSx(filestr);
% % S = load('160505New2_pop.mat');
% S = load('160505New_pop.mat');

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
figure(3); hold off;
% x_off = 107700;   %this is the offset of the ns4 file when the noise reflects csv start 170112
% x_off = 110530;   %this is the offset of the ns4 file when the noise reflects csv start 170118
% x_off = 78935;    %this is the offset of the ns4 file when the noise reflects csv start 170111
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

figure(4); hold off;

plot(x_vals , csv_data(:,2) ,'LineWidth',1);
hold on;
plot(x_vals(intra_mvmt) , csv_data(intra_mvmt,2) ,'.');   %plot the intra-mvment periods on top of the raw velocity

%% separate out into valid movements
mvmt_win = [10,25];
starts = find([0 diff(intra_mvmt)]>10);                   %specifically finds the gaps between movement vectors & makes the next value = start

starts_valid = find(intra_mvmt(starts)>mvmt_win(1) & intra_mvmt(starts)<(numel(csv_data(:,2))-mvmt_win(2)));

starts_clean = starts(starts_valid);
plot(x_vals(intra_mvmt(starts_clean)) , csv_data(intra_mvmt(starts_clean),2) ,'o');

[wheel_moves] = TNC_ExtTrigWins(csv_data(:,2)',intra_mvmt(starts_clean),mvmt_win); %sink.wins = csv(:,2), indexed by the starts_clean with an expanded window; ea. row shows the indexed vals from 10 backwards/25 ahead of the meaningful index

for i=1:size(wheel_moves.wins,1)
    magnitude(i) = trapz(wheel_moves.wins(i,mvmt_win(1):20));  %trapz(y) returns the approximate integral of y
end                                                            %gives a per col. readout of how big ea. movement was

[vals, inds] = sort(magnitude);
figure(5); imagesc(wheel_moves.wins(inds,:),[-50 50]); colormap(mapName);    %x is the movement vectors sorted by velocity magnitude 
    xlabel('vector window'); ylabel('movement vector, sorted by vel mag');

export_mvmt_data_for_psths.times        = x_vals(intra_mvmt(starts_clean))'; %timestamps of the movement vector starts
export_mvmt_data_for_psths.magnitude    = magnitude;
export_mvmt_data_for_psths.mag_sort_i   = inds;

%% examine some psths associated with movements
numUnits = numel( PopData.session(1).unit );
figure(10); clf;
for j=1:numUnits
    
    tmp = PopData.session(1).unit(j).ts;
    delta = zeros(1,ceil(max(tmp)));
    delta(round(tmp)) = 1;

    tmpSmooth = conv(delta,Gaussian,'same');
    
    figure(10);
    subplot(1,3,[1 2]);
    plot(tmpSmooth + j*0.1); hold on;
    
    subplot(1,3,3);
    times = export_mvmt_data_for_psths.times;
    valid_valid = find(times>9.2e4 & times<2.07e6);
    [vals, inds] = sort(magnitude(valid_valid));
    [sink_tmp]   = TNC_ExtTrigWins(tmpSmooth,times(valid_valid),[750,1000]);
    
    plot((sink_tmp.avg - mean(sink_tmp.avg(1:500))) + (ones(1,numel(sink_tmp.avg)).*j.*0.01)); hold on;
    drawnow;
    
    if j==5
       
        figure(11);
        imagesc(sink_tmp.wins(inds,:));
        
        figure(12); hold off;
        plot(tmpSmooth.*250);
        hold on;
        plot(x_vals , csv_data(:,2) ,'LineWidth',2);        
        
    end
title(filestr);
end


%%