%Updated in Sep. 2017 for all movement alignments to older tagging/recording data:
%% Must set the timescale offset at line 60ish per session!

%This script serves as the moveBehaveExtract file to generate a "bh"-like structure for generating psth's aligned to mvment onset.

% Updated 4.27.18: 
    %don't really need the x_offset since this is calculated directly from ns4 file's first trial start
    %cue aligned in figures is from CSV file; rxn time calculated using all csv info put into ns4 file times
    %Works for 171013 file
% Updated 4.12.18 - uses ns4 signal for perfect cue alignment 171012-171013
% Updated 4.5.18 for cue-alignment
% Updated 3.30.18 for longer pre-movement baseline (2 sec) - saves as bh5 file
% Updated 2.20.18 to show individual units' psth w/in their own figure
% Updated 2.21.18 to separate out ipsi/contra movements 
% Updated 2.22.18 to align movements within the trial-start structure - works!
% Updated 3.8.18 to include correct/incorrect movements (ipsi/contra), trial starts (for all these), and rxn times for all these)

% Run moverScript_mid and moveBehaveExtract_mid to generate the Pop Data file for 1701_11, 12, 18 data
%170111 & 170112 sessions = L hemi recording; Thus, +vel = L mvment & ipsi
%170118 session           = R hemi recording; Thus, +vel = L mvment & contra

% Run moverScriptNew & moveBehaveExtract_new for 171012, 13 data
% 171012, 13              = L hemisphere;     Thus, +vel = L mvment & ipsi

%% Load the file of interest:
% threshold = 500;  % all 2017 data - may differ, though, with different pre-amp settings
% threshold = 2.0e3;  % all 2017 data - may differ, though, with different pre-amp settings

% targetName = '171105';
% % fpath = '/Volumes/My Passport for Mac/171013/behaveNewUnits';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/newBehaveUnits/behaveChunks';
% csv_data = dlmread('mVgatthreetag_2015_11_05_141004.csv',',',2,0);
% filestr = '151105tag1001.ns4';
% smooth_off = 0;

targetName = '171012';
% fpath = '/Volumes/My Passport for Mac/171012/newBehaveUnits';
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/171012/newBehaveUnits';
csv_data = dlmread('mSC2Tag_2017_10_12_164510pXY.csv',',',2,0);
filestr = '/171012001.ns4';
smooth_off = 0;
% smooth_off = 800;
% excludeNeurons = [1 5 7 9 14 15 28 33 35]; %should be 9 units

% targetName = '171013';
% % fpath = '/Volumes/My Passport for Mac/171013/behaveNewUnits';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/171013';
% csv_data = dlmread('mSC2tag_2017_10_13_154511pXY.csv',',',2,0);
% filestr = '171013002.ns4';
% smooth_off = 0;
% % smooth_off = 480; %value determined from xvelocity wheel offset in contWheel2.m
% % smooth_off = 590;   %for cue alignment %This and the offset were chosen for bhCue2.mat
% % smooth_off = 750;   %for cue alignment

% targetName = '170111';
% % fpath = '/Volumes/My Passport for Mac/170111/newBehaveUnits';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170111';
% csv_data = dlmread('mVgatfive_2017_01_11_174943pXY.csv',',',2,0);
% filestr = '/170111002.ns4';
% smooth_off = 0;

% targetName = '170112';
% % fpath = '/Volumes/My Passport for Mac/170112/behaveSegs/behaveNewUnits';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170112/behaveNew';
% csv_data = dlmread('mVgatfive_2017_01_12_165817pXY.csv',',',2,0);
% filestr = '/170112002.ns4';
% smooth_off = 0;

% targetName = '170118';
% % fpath = '/Volumes/My Passport for Mac/Vgatfive/170118/behaveNewUnits';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170118/behaveNewUnits';
% csv_data = dlmread('mVgatfive_2017_01_18_163328pXY.csv',',',2,0);
% filestr = '/170118002.ns4';
% smooth_off = 0;

filestr2 = strcat(fpath, '/', filestr);
[data] = openNSx(filestr2);              %MUST LOAD THE openNSx.m file from 2016 for this to work (currently in forTonic_v2 folder, mac HD)

% S = load('170111newBehave_pop');
% S = load('170112newBehave_pop');
% S = load('170118_newBehave_pop.mat');
S = load('171012_pop.mat');
% S = load('171013_pop.mat');
% S = load('151105_update_pop.mat');

rxnTimesFast = 100;
rxnTimesSlowMin = 400;
rxnTimesSlowMax = 600;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PopData = S.PopData;
clear S;

%plot to calculate x_off for trial start from ns4 file: for 170111, 12, 18:, 171012, 
chan.start = 129;                                    %trial start channel (Ain 1)
trialStrtsChan = find(data.MetaTags.ChannelID==chan.start);
trialStrtsChan1 = double(data.Data(trialStrtsChan,:)); 
trialStrts = decimate(trialStrtsChan1,10);           %convert these to ms (/10)

figure; hold on;
plot(trialStrts,'k');

%get wheel data from ns4 file: for 170111, 12, 18:, 171012, 
chan.wheel = 131;                                    %trial start channel (Ain 1)
wheelChan = find(data.MetaTags.ChannelID==chan.wheel);
wheelChan1 = double(data.Data(wheelChan,:)); 
wheelMv1 = decimate(wheelChan1,10);                   %convert these to ms (/10)
plot(wheelMv1,'bl');                                  %plot wheel NS4 velocity signal in blue

% wheel_offset = -80;                                   %find the y-axis offset of the wheel baseline
% wheelMv2 = wheelMv1 - wheel_offset;
% wheel_sgolayF = sgolayfilt(wheelMv2,9,101);           %this filter is good for oscillating baselines
% wheelMvgolay = wheel_sgolayF - wheel_offset;
% wheelMv_true = wheelMvgolay;

%% This is to generate the plot of magnitudes for the movement vectors & set up psths
Sigma = 24;
t = 1:1:Sigma.*30;
Gaussian = ( 1./( sqrt(2.*pi.*Sigma.*Sigma) ) ) .* exp( -(t-Sigma.*15).^2 ./ (2.*Sigma).^2 );
integral = trapz(Gaussian);                         %trapz(y) returns the approximate integral of y
Gaussian = Gaussian./integral;
[mapName] = TNC_CreateRBColormap(8,'mbr');          %calls color map function

%% Set the offset based on figure 1

% x_off = 764510;       %151105
% LtimesMax = 2.2e6;
% RtimesMax = 2.2e6;

% x_off = 78935;        %170111
% LtimesMax = 3.0e6;
% RtimesMax = 3.0e6;

% x_off = 107700;       %previously thought for 170112, but next line better? No
% % x_off = 107000;     %this is the offset of the ns4 file when the noise reflects csv start 170112; 171012
% LtimesMax = 2.1e6;
% RtimesMax = 2.1e6;
 
x_off = 110800;       %171012
LtimesMax = 15e5;
RtimesMax = 15e5;
trialStrtMin = 1e5;

% x_off = 110800;       %171013 bhCue2 works best - movements are right shifted, variable trial starts
% x_off = 102800;       %bhCue5 begin of the csv session start signal - nope; movements are still r-shifted 
% x_off = 103820;       %bhCue6 - end of the csv session start signal
% x_off = 103840;       %past end of the csv session start signal bhCue4
% x_off = 103838;       %past end of the csv session start signal bhCue3
% 

% x_off = 110800;        %171013
% trialStrtMin = 1.1e5;
% LtimesMax = 15760e5;    
% RtimesMax = 15760e5;

% x_off = 110530;        %170118
% LtimesMax = 2.5e6;
% RtimesMax = 2.5e6;
% trialStrtMin = 1.1e5;

%% Align trial start cue using the ns4 file

% figure; hold on;
% plot(trialStrts(trialStrtMin:end),'k');   %nope - I need to retain the indices since they are the timestamp
trialStrt_inds = 1:length(trialStrts);      %retain the indices
trialStrt_true = find(trialStrts>500);
trialStart_almost = find(trialStrt_true > trialStrtMin);
trialStart_almost2 = trialStrt_true(trialStart_almost);
trialStart_diff = diff(trialStart_almost2);
findDiff = find(trialStart_diff > 3000);    %no trial(i+1) can occur < 3000 ms from trial(i)
trialStarts_valid = trialStart_almost2(findDiff+1);

export_mvmt_data_for_psths.trialStarts_validNS4 = trialStarts_valid;

%% Load the taskbase structure for CSV "correct" trial alignments:  2.22.18

t_holdVel = 10;
x_vals = (csv_data(:,1)-csv_data(1,1))+x_off + smooth_off; %add the NS4 x-offset (time) to the csv data's timestamps
d_x_vals = [0 diff(x_vals')];                 %the diff of timestamps from csv file
smth_wheel_v = conv(abs(csv_data(:,2)) , [zeros(1,4) ones(1,5) 0] , 'same' )'; %The convolution of two vecs u (velocity) & v (row of ones padded by 0s), will be the overlap under the points as v slides across u

% plot(x_vals , csv_data(:,2) ,'LineWidth',2);  %blue; plots col2 (raw velocity) w/ the time offset calculated from the ns4 file
hold on; 
plot(x_vals , smth_wheel_v ,'LineWidth',2);     %orange; plots the smoothed wheel velocity data (all positive-going); 

%% find intra movement periods
threshold = 5;  %all 2017 data

intra_mvmt = find(smth_wheel_v>threshold & d_x_vals<30);  %this is a "meaningful" movement: find the overlapping instances for smoothed velocity > 5 and the diff(timestamps) < 30 ms

figure; hold on;
plot(x_vals , csv_data(:,2) ,'LineWidth',1);              %plots csv data in blue; inter-mvment-interval in red (fig. 4)
hold on;
velMvmntsZeroed = x_vals(intra_mvmt);                     %pulls out the ts of CSV file for changing velocities w/ offset factored in
velMvmntsRaw = csv_data(intra_mvmt,2);                    %those raw velocities 
plot(x_vals(intra_mvmt) , csv_data(intra_mvmt,2) ,'.');   %plot the intra-mvment periods on top of the raw velocity (fig 4)
hold off;

%% separate out into valid movements
mvmt_win = [10,25];
starts = find([0 diff(intra_mvmt)]>10);                    %starts = inds of mvment vecs - specifically finds the large gaps between movement vectors & makes the next value = start
% starts = find([0 diff(intra_mvmt)]>5);                    %starts = inds of mvment vecs - specifically finds the large gaps between movement vectors & makes the next value = start
starts_valid = find(intra_mvmt(starts)>mvmt_win(1) & intra_mvmt(starts)<(numel(csv_data(:,2))-mvmt_win(2)));  %these are the csv ts inds of first mvment onset
starts_clean = starts(starts_valid);                       %indices of the starts of meaningful mvments of csv file
figure; hold on; plot(x_vals(intra_mvmt(starts_clean)) , csv_data(intra_mvmt(starts_clean),2) ,'o');%plots in yellow (still fig 4)
%                  plot(x_vals(intra_mvmtNew(starts_cleanNew)) , csv_data(intra_mvmtNew(starts_cleanNew),2) ,'o');%plots in yellow (still fig 4)
%                  plot(x_vals(intra_mvmtNew(starts_cleanNew)) , 0 ,'r.');%plots in yellow (still fig 4)

%wheel_moves creates #rows of ea. mvmnt bout's vel: 
%when trial starts & bar = -pos, correct rightward mvment has -vel;when +position, correct leftward mv has + vel.
[wheel_moves] = TNC_ExtTrigWins(csv_data(:,2)',intra_mvmt(starts_clean),mvmt_win); %sink.wins = csv(:,2), indexed by the starts_clean with an expanded window; ea. row shows the indexed vals from 10 backwards/25 ahead of the meaningful index
% [wheel_movesNew] = TNC_ExtTrigWins(csv_data(:,2)',intra_mvmtNew(starts_cleanNew),mvmt_win); %sink.wins = csv(:,2), indexed by the starts_clean with an expanded window; ea. row shows the indexed vals from 10 backwards/25 ahead of the meaningful index

%170112 session = L hemi; Thus, +vel = L mvment & ipsi

for i=1:size(wheel_moves.wins,1)
    magnitude(i) = trapz(wheel_moves.wins(i,mvmt_win(1):20));  %trapz(y) returns the approximate integral of y
end                                                            %gives a per col. readout of how big ea. movement was
magnitude = magnitude';                 %EAS: find the mvments that had more than 3ish vel readouts (of >10-30, ea) should be ~50 or greater
% t_holdVel = 30;

trueMvsPre = abs(magnitude) > t_holdVel;
trueMvsInds = find(trueMvsPre == 1);    %EAS: these are the indices of magnitude that are likely meaningful movements

trueMvsInds_R = find(magnitude < -t_holdVel);  %EAS: indices of magnitude for large/fast right mvments
trueMvsInds_L = find(magnitude > t_holdVel);   %EAS: inds of mag for large/fast left mvments

trueMv_R_mag = magnitude(trueMvsInds_R);       %Magnitude of all R inds
trueMv_L_mag  = magnitude(trueMvsInds_L);

[vals_vel_R, inds_vel_R] = sort(trueMv_R_mag); %Sorted magnitude of all R mvs 
[vals_vel_L, inds_vel_L] = sort(trueMv_L_mag);

[vals_vel, inds_vel] = sort(magnitude);  %Sorts all inds of the original magnitude (L & R)      
figure;
imagesc(wheel_moves.wins(inds_vel,:),[-50 50]); colormap(mapName);           %x is the movement vectors sorted by velocity magnitude 
xlabel('vector window'); ylabel('movement vectors, sorted by vel mag');  %should be separating ipsi v contra (one is above, other below).

export_mvmt_data_for_psths.times        = x_vals(intra_mvmt(starts_clean))'; %timestamps of the movement vector starts
export_mvmt_data_for_psths.magnitude    = magnitude;
export_mvmt_data_for_psths.mag_sort_i   = inds_vel;

export_mvmt_data_for_psths.Lefts        = trueMvsInds_L;                      %magnitude inds for L mvments
export_mvmt_data_for_psths.Rights       = trueMvsInds_R;                      %magnitude inds for R mvments
export_mvmt_data_for_psths.LeftsTimes   = export_mvmt_data_for_psths.times(trueMvsInds_L); %timestamps of L mvments
export_mvmt_data_for_psths.RightsTimes  = export_mvmt_data_for_psths.times(trueMvsInds_R); %timestamps of R mvments


if exist('csv_data', 'var');
    [filenamestrT, path] = uigetfile('*_tb.mat*','select the taskbase file', fpath)       % get the actual # correct trials (and TS) from csv file
end
taskbase = strcat(path, filenamestrT);
load(taskbase);
startTimesCsv = taskbase.trialStartTimes';
diffStartTimesCsv = diff(startTimesCsv);
diffStartTimesCsv2 = diffStartTimesCsv(2:end);
DifftrialStarts_valid = diff(trialStarts_valid');  %from ns4 file, but the first trial will be missing; 
% [Ids, DiffVals]  = findClosestId2Val(diffStartTimesCsv2(1), DifftrialStarts_valid);
x = DifftrialStarts_valid(1:3) == diffStartTimesCsv2(1:3)+1;
y = DifftrialStarts_valid(1:3) == diffStartTimesCsv2(1:3)-1;
x_vec = find(x == 1);
y_vec = find(y == 1);

if x_vec > 0
    first_alignTS = trialStarts_valid(x_vec);
else if y_vec > 0
        first_alignTS = trialStarts_valid(y_vec);
    
    end
end
first_align = find(trialStarts_valid==first_alignTS);
trueTrialStarts_valid = trialStarts_valid(first_align:end);  %these are the true trial starts ns4
 
startTimesCSVZero1 = startTimesCsv - startTimesCsv(1) + first_alignTS;
trialStarts = startTimesCSVZero1;                             %these are the ns4-aligned CSV starts
  
%% Load the taskbase fields:

%These are correct AND incorrect trials, first mvment onset that lead to centering:
wheelLfirstTimes = taskbase.wheelLfirstTime';
% wheelLfirstZero1 = wheelLfirstTimes-csv_data(1,1) + x_off + smooth_off;
wheelLfirstZero1 = wheelLfirstTimes - startTimesCsv(1) + first_alignTS;
wheelLfirstsZero = unique(wheelLfirstZero1);                %correct for added x_off to zeros

wheelILfirstTimes = taskbase.wheelILfirstTime';
wheelILfirstZero1 = wheelILfirstTimes - startTimesCsv(1) + first_alignTS;
wheelILfirstsZero = unique(wheelILfirstZero1);               %correct for added x_off to zeros

wheelRfirstTimes = taskbase.wheelRfirstTime';
wheelRfirstZero1 = wheelRfirstTimes - startTimesCsv(1) + first_alignTS;
wheelRfirstsZero = unique(wheelRfirstZero1);                 %correct for added x_off to zeros

wheelIRfirstTimes = taskbase.wheelIRfirstTime';
wheelIRfirstZero1 = wheelIRfirstTimes - startTimesCsv(1) + first_alignTS;
wheelIRfirstsZero = unique(wheelIRfirstZero1);               %correct for added x_off to zeros

leftCenteredTimes = taskbase.leftCenteredTimes';
leftCenteredTimesZero1 = leftCenteredTimes - startTimesCsv(1) + first_alignTS;
leftCenteredTimesZero = unique(leftCenteredTimesZero1);

rightCenteredTimes = taskbase.rightCenteredTimes';
rightCenteredTimesZero1 = rightCenteredTimes - startTimesCsv(1) + first_alignTS;
rightCenteredTimesZero = unique(rightCenteredTimesZero1);

leftIncorrectTimes = taskbase.leftIncorrectTimes';
leftIncorrectTimesZero1 = leftIncorrectTimes - startTimesCsv(1) + first_alignTS;
leftIncorrectTimesZero = unique(leftIncorrectTimesZero1);

rightIncorrectTimes = taskbase.rightIncorrectTimes';
rightIncorrectTimesZero1 = rightIncorrectTimes - startTimesCsv(1) + first_alignTS;
rightIncorrectTimesZero = unique(rightIncorrectTimesZero1);

%% Now fit these movements to trial start structure:  added 2.22.18
% L rewarded trials' start times
wheelLtrialStarts = [];
for a = 1:length(trialStarts)
    for b = 1:length(wheelLfirstsZero)
        if a == numel(trialStarts)
            break
        end
        if b == numel(wheelLfirstsZero)
            break
        end
        if wheelLfirstsZero(b) > trialStarts(a) && wheelLfirstsZero(b) < trialStarts(a+1) && wheelLfirstsZero(b) < wheelLfirstsZero(b+1) % just to be sure
           wheelLtrialStarts(a) = trialStarts(a);   %has trial inds
        end
    end
end

wheelILtrialStarts = [];
for a = 1:length(trialStarts)
    for b = 1:length(wheelILfirstsZero)
        if a == numel(trialStarts)
            break
        end
        if b == numel(wheelILfirstsZero)
            break
        end
        if wheelILfirstsZero(b) > trialStarts(a) && wheelILfirstsZero(b) < trialStarts(a+1) && wheelILfirstsZero(b) < wheelILfirstsZero(b+1) % just to be sure
           wheelILtrialStarts(a) = trialStarts(a);   %has trial inds
        end
    end
end

%R rewarded trials' start times
wheelRtrialStarts = [];
for a = 1:length(trialStarts)
    for b = 1:length(wheelRfirstsZero)
        if a == numel(trialStarts)
            break
        end
        if b == numel(wheelRfirstsZero)
            break
        end
        if wheelRfirstsZero(b) > trialStarts(a) && wheelRfirstsZero(b) < trialStarts(a+1) && wheelRfirstsZero(b) < wheelRfirstsZero(b+1) % just to be sure
           wheelRtrialStarts(a) = trialStarts(a);   %has trial inds
        end
    end
end

wheelIRtrialStarts = [];
for a = 1:length(trialStarts)
    for b = 1:length(wheelIRfirstsZero)
        if a == numel(trialStarts)
            break
        end
        if b == numel(wheelIRfirstsZero)
            break
        end
        if wheelIRfirstsZero(b) > trialStarts(a) && wheelIRfirstsZero(b) < trialStarts(a+1) && wheelIRfirstsZero(b) < wheelIRfirstsZero(b+1) % just to be sure
           wheelIRtrialStarts(a) = trialStarts(a);   %has trial inds
        end
    end
end

export_mvmt_data_for_psths.LrewTrialStarts = wheelLtrialStarts;
export_mvmt_data_for_psths.ILrewTrialStarts = wheelILtrialStarts;
export_mvmt_data_for_psths.trialStartsValid = trialStarts;
export_mvmt_data_for_psths.RrewTrialStarts = wheelRtrialStarts;
export_mvmt_data_for_psths.IRrewTrialStarts = wheelIRtrialStarts;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Use this for most sessions' correct L movements: eliminates first L movements that don't fit to an extracted trial start into Blackrock
wheelLfirstsAligned = []; wheelLfirstsAlignedInds = [];LrewTrialStarts2 = [];
LrewTrialStartsZeroed1 = find(wheelLtrialStarts>0);
LrewTrialStartsZeroed = wheelLtrialStarts(LrewTrialStartsZeroed1);
for i = 1:length(LrewTrialStartsZeroed)
    if numel(LrewTrialStartsZeroed) < numel(wheelLfirstsZero) && i == numel(LrewTrialStartsZeroed)
        break
    else if i == numel(wheelLfirstsZero)
            break
        else if wheelLfirstsZero(i) > LrewTrialStartsZeroed(i) && wheelLfirstsZero(i) < LrewTrialStartsZeroed(i+1) % just to be sure
                wheelLfirstsAligned(i) = wheelLfirstsZero(i);
                LrewTrialStarts2(i) = LrewTrialStartsZeroed(i);
            else if wheelLfirstsZero(i) < LrewTrialStartsZeroed(i) && wheelLfirstsZero(i+1) > LrewTrialStartsZeroed(i) && wheelLfirstsZero(i+1) < LrewTrialStartsZeroed(i+1)
%             else if wheelLfirstsZero(i+1) > LrewTrialStartsZeroed(i) && wheelLfirstsZero(i+1) < LrewTrialStartsZeroed(i+1)                    
                    wheelLfirstsAligned(i) = wheelLfirstsZero(i+1);
                    LrewTrialStarts2(i) = LrewTrialStartsZeroed(i);
                end
            end
        end
    end
end

%% For Incorrect R movements:
% Use this for most sessions' correct L movements: eliminates first L movements that don't fit to an extracted trial start into Blackrock
wheelIRfirstsAligned = []; wheelIRfirstsAlignedInds = []; IRTrialStarts2 = [];
IRrewTrialStartsZeroed1 = find(wheelIRtrialStarts>0);
IRrewTrialStartsZeroed = wheelIRtrialStarts(IRrewTrialStartsZeroed1);
for i = 1:length(IRrewTrialStartsZeroed)
    if numel(IRrewTrialStartsZeroed) < numel(wheelIRfirstsZero) && i == numel(IRrewTrialStartsZeroed)
        break
    else if i == numel(wheelIRfirstsZero)
            break
        else if wheelIRfirstsZero(i) > IRrewTrialStartsZeroed(i) && wheelIRfirstsZero(i) < IRrewTrialStartsZeroed(i+1)  % just to be sure
                wheelIRfirstsAligned(i) = wheelIRfirstsZero(i);
                IRTrialStarts2(i) = IRrewTrialStartsZeroed(i);
            else if wheelIRfirstsZero(i) < IRrewTrialStartsZeroed(i) && wheelIRfirstsZero(i+1) > IRrewTrialStartsZeroed(i) && wheelIRfirstsZero(i+1) < IRrewTrialStartsZeroed(i+1)
                    wheelIRfirstsAligned(i) = wheelIRfirstsZero(i+1);
                    IRTrialStarts2(i) = IRrewTrialStartsZeroed(i);
                end
            end
        end
    end
end

%% Use this for most sessions' correct R movements: eliminates first L movements that don't fit to an extracted trial start into Blackrock
wheelRfirstsAligned = []; wheelRfirstsAlignedInds = []; RrewTrialStarts2 = [];
RrewTrialStartsZeroed1 = find(wheelRtrialStarts>0);
RrewTrialStartsZeroed = wheelRtrialStarts(RrewTrialStartsZeroed1);
for i = 1:length(RrewTrialStartsZeroed)
    if numel(RrewTrialStartsZeroed) < numel(wheelRfirstsZero) && i == numel(RrewTrialStartsZeroed)
        break
    else if i == numel(wheelRfirstsZero)
            break
        else if wheelRfirstsZero(i) > RrewTrialStartsZeroed(i) && wheelRfirstsZero(i) < RrewTrialStartsZeroed(i+1)  % just to be sure
                wheelRfirstsAligned(i) = wheelRfirstsZero(i);
                RrewTrialStarts2(i) = RrewTrialStartsZeroed(i);
            else if wheelRfirstsZero(i) < RrewTrialStartsZeroed(i) && wheelRfirstsZero(i+1) > RrewTrialStartsZeroed(i) && wheelRfirstsZero(i+1) < RrewTrialStartsZeroed(i+1)
                    wheelRfirstsAligned(i) = wheelRfirstsZero(i+1);
                    RrewTrialStarts2(i) = RrewTrialStartsZeroed(i);
                end
            end
        end
    end
end

%% For incorrect left trials:
wheelILfirstsAligned = []; wheelILfirstsAlignedInds = []; ILTrialStarts2 = [];
ILrewTrialStartsZeroed1 = find(wheelILtrialStarts>0);
ILrewTrialStartsZeroed = wheelILtrialStarts(ILrewTrialStartsZeroed1);
for i = 1:length(ILrewTrialStartsZeroed)
    if numel(ILrewTrialStartsZeroed) < numel(wheelILfirstsZero) && i == numel(ILrewTrialStartsZeroed)
        break
    else if i == numel(wheelILfirstsZero)
            break
        else if wheelILfirstsZero(i) > ILrewTrialStartsZeroed(i) && wheelILfirstsZero(i) < ILrewTrialStartsZeroed(i+1) % just to be sure
                wheelILfirstsAligned(i) = wheelILfirstsZero(i);
                ILTrialStarts2(i) = ILrewTrialStartsZeroed(i);
            else if wheelILfirstsZero(i) < ILrewTrialStartsZeroed(i) && wheelILfirstsZero(i+1) > ILrewTrialStartsZeroed(i) && wheelILfirstsZero(i+1) < ILrewTrialStartsZeroed(i+1)
                    wheelILfirstsAligned(i) = wheelILfirstsZero(i+1);
                    ILTrialStarts2(i) = ILrewTrialStartsZeroed(i);

                end
            end
        end
    end
end

 export_mvmt_data_for_psths.wheelLfirstValid = wheelLfirstsAligned;
 export_mvmt_data_for_psths.wheelRfirstValid = wheelRfirstsAligned;
 export_mvmt_data_for_psths.wheelILfirstValid = wheelILfirstsAligned;
 export_mvmt_data_for_psths.wheelIRfirstValid = wheelIRfirstsAligned;

 wheelALL_LfirstValid = sort(horzcat(wheelIRfirstsAligned, wheelLfirstsAligned));
 wheelALL_RfirstValid = sort(horzcat(wheelILfirstsAligned, wheelRfirstsAligned));
 export_mvmt_data_for_psths.wheelALL_LfirstValid = sort(horzcat(wheelIRfirstsAligned, wheelLfirstsAligned));
 export_mvmt_data_for_psths.wheelALL_RfirstValid = sort(horzcat(wheelILfirstsAligned, wheelRfirstsAligned));

% All_L_trialStarts = sort(horzcat(LrewTrialStartsZeroed, IRrewTrialStartsZeroed));
% All_R_trialStarts = sort(horzcat(RrewTrialStartsZeroed, ILrewTrialStartsZeroed));

% Use these trial starts:
All_L_trialStarts = sort(horzcat(LrewTrialStarts2, IRTrialStarts2));
All_R_trialStarts = sort(horzcat(RrewTrialStarts2, ILTrialStarts2));
export_mvmt_data_for_psths.All_L_trialStarts = All_L_trialStarts;
export_mvmt_data_for_psths.All_R_trialStarts = All_R_trialStarts;

%% Get wheel velocities for these TS: use the centering time or use the vel already calculated by magnitude code above?
% Better to just extract the velocities from the taskbase file - don't use yet
pooL = ismember(export_mvmt_data_for_psths.LeftsTimes, wheelALL_LfirstValid - 46) | ismember(export_mvmt_data_for_psths.LeftsTimes, wheelALL_LfirstValid - 45) | ismember(export_mvmt_data_for_psths.LeftsTimes, wheelALL_LfirstValid - 44) | ismember(export_mvmt_data_for_psths.LeftsTimes, wheelALL_LfirstValid - 43) | ismember(export_mvmt_data_for_psths.LeftsTimes, wheelALL_LfirstValid - 42)| ismember(export_mvmt_data_for_psths.LeftsTimes, wheelALL_LfirstValid - 41); 
pooR = ismember(export_mvmt_data_for_psths.RightsTimes, wheelALL_RfirstValid - 46) | ismember(export_mvmt_data_for_psths.RightsTimes, wheelALL_RfirstValid - 45) | ismember(export_mvmt_data_for_psths.RightsTimes, wheelALL_RfirstValid - 44) | ismember(export_mvmt_data_for_psths.RightsTimes, wheelALL_RfirstValid - 43) | ismember(export_mvmt_data_for_psths.RightsTimes, wheelALL_RfirstValid - 42)| ismember(export_mvmt_data_for_psths.RightsTimes, wheelALL_RfirstValid - 41); 

findSimL = find(pooL == 1);  %indices of the same values bt. raw mvmnts & trial-aligned mvments
findSimR = find(pooR == 1);
SimL = export_mvmt_data_for_psths.LeftsTimes(pooL); %mvmt times for the above indices that are the same
SimR = export_mvmt_data_for_psths.RightsTimes(pooR);

mag_L = trueMv_L_mag(pooL);
mag_R = trueMv_R_mag(pooR);
[sort_vel_L, indsSort_vel_L] = sort(mag_L);    %Sorted magnitude of all matched L mvs 
[sort_vel_R, indsSort_vel_R] = sort(mag_R);    %Sorted magnitude of all matched R mvs 

%% examine some psths associated with movements & set up matrix for heat plotting all units

numUnits = numel( PopData.session(1).unit );
figure; 
for j=1:numUnits  
    
    tmp = PopData.session(1).unit(j).ts;
    delta = zeros(1,ceil(max(tmp)));
    delta(round(tmp)) = 1;

    tmpSmooth = conv(delta,Gaussian,'same');

    times = export_mvmt_data_for_psths.times;       %These are all mvmnt times
    Ltimes = export_mvmt_data_for_psths.LeftsTimes; %These are all L mvmnt times
    Rtimes = export_mvmt_data_for_psths.RightsTimes;%These are all R mvmnt times
    
    wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;  %These are the trial-aligned L mvment times
    wheelsRfirstValid = export_mvmt_data_for_psths.wheelRfirstValid;  %These are the trial-aligned R mvment times
    
    wheelsAll_LfirstValid =  export_mvmt_data_for_psths.wheelALL_LfirstValid;  %ADDED 3.6.18
    wheelsAll_RfirstValid = export_mvmt_data_for_psths.wheelALL_RfirstValid;   %ADDED 3.6.18
    %% Sort out reaction times: use wheelsAll_RfirstValid and wheelsAll_LfirstValid 3.8.18
    
    % All aligned trials using the CSV times:
    rxnTimesL = wheelsAll_LfirstValid - All_L_trialStarts; %in order of occurence
    rxnTimesR = wheelsAll_RfirstValid - All_R_trialStarts;
    
    [sort_rxn_L, indsSort_rxn_L] = sort(rxnTimesL);    %Sorted reaction times of all L mvs 
    [sort_rxn_R, indsSort_rxn_R] = sort(rxnTimesR);    %Sorted reaction times of all R mvs 

    zoot_L = All_L_trialStarts(indsSort_rxn_L);
    zoot_R = All_R_trialStarts(indsSort_rxn_R);
   
    % Separates out corrects/incorrects, uses the CSV times:
    rxnTimesL_1 = wheelsLfirstValid - LrewTrialStarts2;
    rxnTimesR_1 = wheelsRfirstValid - RrewTrialStarts2;
    rxnTimesIL_1 = wheelILfirstsAligned - ILTrialStarts2;
    if isempty(wheelIRfirstsAligned)
        rxnTimesIR_1 = NaN;
    else
        rxnTimesIR_1 = wheelIRfirstsAligned - IRTrialStarts2;
    end
    
    rxnTimesL = sort(horzcat(rxnTimesL_1, rxnTimesIR_1)); %includes all incorrect/correct trials
    rxnTimesR = sort(horzcat(rxnTimesR_1, rxnTimesIL_1)); %includes all incorrect/correct trials
    
    rxnTimesLslows = find(rxnTimesL > rxnTimesSlowMin & rxnTimesL < rxnTimesSlowMax);
    rxnTimesL400 = rxnTimesL(rxnTimesLslows);
    rxnTimesL400mv = wheelsAll_LfirstValid(rxnTimesLslows);  %includes all incorrect/correct trials
    
    rxnTimesLfast = find(rxnTimesL < rxnTimesFast);
    rxnTimesL200 = rxnTimesL(rxnTimesLfast);
    rxnTimesL200mv = wheelsAll_LfirstValid(rxnTimesLfast);   %includes all incorrect/correct trials

    rxnTimesRslows = find(rxnTimesR > rxnTimesSlowMin & rxnTimesR < rxnTimesSlowMax);
    rxnTimesR400 = rxnTimesR(rxnTimesRslows);
    rxnTimesR400mv = wheelsAll_RfirstValid(rxnTimesRslows);

    rxnTimesRfast = find(rxnTimesR < rxnTimesFast);
    rxnTimesR200 = rxnTimesR(rxnTimesRfast);
    rxnTimesR200mv = wheelsAll_RfirstValid(rxnTimesRfast);    %includes all incorrect/correct trials
    
psthWin = [2.0e3,2e3];
  
%% Plot the L-movement trial start cue, sorted by reaction time
    figure; hold on;
    subplot(2,2,1); hold on;
    
    Ltimes_valid = zoot_L;
    
    [sink_tmpL]   = TNC_ExtTrigWinsEAS(tmpSmooth,Ltimes_valid,[2000,2000]); %[-win +win], so mvment onset = ~700
    if ~isfinite(sink_tmpL.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
        j = j+1;
    else
    unitNum = num2str(j);
    unitname = strcat(' u#', unitNum, 'cue-aligned');
   
    alignVar1 = Ltimes_valid;
    PopData.currParams.stableTrials = -1;                                         % subtracts the first "trial" when calling AlignRaster

    [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session.unit(j).ts, PopData.currParams.stableTrials, alignVar1, psthWin, 1, 1); 
    shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'bl');  %blue for left; for 151105 dataset = ipsi
%     title(unitname); 
    
    Lstarts_trials = respCSS.image.aligned;
    end
        
    %% Plot the R-movements' trial start cue:
    
    Rtimes_valid = zoot_R;
    
    [sink_tmpR]   = TNC_ExtTrigWinsEAS(tmpSmooth,Rtimes_valid,[2000,2000]); %[-win +win], so mvment onset = ~700
    if ~isfinite(sink_tmpR.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
        j = j+1;
    else
    unitNum = num2str(j);
    unitname = strcat(' u#', unitNum, 'cue-aligned', 'all');
    title(unitname); 

    alignVar2 = Rtimes_valid;
    [respCSS2] = TNC_AlignRasters(tmpSmooth,PopData.session.unit(j).ts, PopData.currParams.stableTrials, alignVar2, psthWin, 1, 1); 
    shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthZ, respCSS2.image.psthZe,'r');  %blue for left; for 151105 dataset = ipsi
    hold off;
        
    Lstarts_trials = respCSS.image.aligned;
    Rstarts_trials = respCSS2.image.aligned;

    figure; hold on;
    ipsiTrials = Lstarts_trials;
    title(strcat(' u#', unitNum, 'cue-ipsi'));
    imagesc(flip(ipsiTrials));
%     
    [cmap] = TNC_CreateRBColormap(1024,'gp');
    colormap(cmap);
%     c_ipsi = colorbar;
%     cmin = -0.25; cmax = 0.25; caxis([cmin cmax]);
    colorbar;

    figure; hold on;
%     contraTrials = Rstarts_trials(inds_RMv,:);
    contraTrials = Rstarts_trials;
    title(strcat(' u#', unitNum, 'cue-contra'));
    imagesc(flip(contraTrials));
    
    [cmap] = TNC_CreateRBColormap(1024,'gp');
    colormap(cmap);
%     c_contra = colorbar;
%     cmin = -0.25; cmax = 0.25; caxis([cmin cmax]);
    colorbar;

    end    

% %% Plot the trial-aligned data:
%     % Plot all L-aligned move onset:
%     subplot(2,2,2); hold on;
% %     wheelsLfirstValid_valid = find(wheelsLfirstValid>9.2e4 & wheelsLfirstValid<LtimesMax);
%     LeftsTimes = export_mvmt_data_for_psths.LeftsTimes;
%     Lmoves_valid = find(LeftsTimes > 9.2e4 & LeftsTimes<LtimesMax);
%     
%     [sink_tmpLaligned]   = TNC_ExtTrigWinsEAS(tmpSmooth,LeftsTimes(Lmoves_valid),[2000,2000]); %[-win +win], so mvment onset = ~700
%     if ~isfinite(sink_tmpLaligned.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
%         j = j+1;
%     else
%     unitNum = num2str(j);
%     unitname = strcat(' u#', unitNum, 'mv-aligned, all');
%     alignVar1 = LeftsTimes(Lmoves_valid);
% 
%     PopData.currParams.stableTrials = -1;        
% %     if j == 20  %for the 170112 dataset
% %         break
% %         else
%         [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session.unit(j).ts, PopData.currParams.stableTrials, alignVar1, psthWin, 1, 1); 
%         shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'bl');  %blue for left; for 151105 dataset = ipsi
%         title(unitname); 
% %         end     %for the 170112 dataset
%     end
%     %% Plot the R-movements:
% %     wheelsRfirstValid_valid = find(wheelsRfirstValid>9.2e4 & wheelsRfirstValid<RtimesMax);
% %     RrewTrialStarts = export_mvmt_data_for_psths.RrewTrialStarts;
% %     RrewTrialStarts_valid = find(RrewTrialStarts > 9.2e4 & RrewTrialStarts<RtimesMax);
% 
%     RightsTimes = export_mvmt_data_for_psths.RightsTimes;
%     Rmoves_valid = find(RightsTimes > 9.2e4 & RightsTimes<LtimesMax);
%     
%     [sink_tmpRaligned]   = TNC_ExtTrigWinsEAS(tmpSmooth, RightsTimes(Rmoves_valid),[2000,2000]); %[-win +win], so mvment onset = ~700
%     if ~isfinite(sink_tmpRaligned.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
%         j = j+1;
%     else
%     unitNum = num2str(j);
%     unitname = strcat(' u#', unitNum, 'mv-aligned, all');
%    
% %     alignVar2 = wheelsRfirstValid(wheelsRfirstValid_valid);
%     alignVar2 =  RightsTimes(Rmoves_valid);
%     
%     [respCSS2] = TNC_AlignRasters(tmpSmooth,PopData.session.unit(j).ts, PopData.currParams.stableTrials, alignVar2, psthWin, 1, 1); 
%     shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthZ, respCSS2.image.psthZe,'r');  %blue for left; for 151105 dataset = ipsi
%     title(unitname); 
%     hold off;
%     end
%    
%    export_mvmt_data_for_psths.LmoveAll = LeftsTimes(Lmoves_valid);
%    export_mvmt_data_for_psths.RmoveAll = RightsTimes(Rmoves_valid);
% 
% %    export_mvmt_data_for_psths.LrewTrueTrialStarts = LrewTrialStarts(LrewTrialStarts_valid);
% %    export_mvmt_data_for_psths.RrewTrueTrialStarts = RrewTrialStarts(RrewTrialStarts_valid);
% %  
% %    
% %     %% Rxn time plots:
% %     % Plot the L-movements'cue on (for fast trials):
% %     subplot(2,2,3); hold on;
% %     rxnTimesL200mv_valid = find(rxnTimesL200mv > 9.2e4 & rxnTimesL200mv< LtimesMax);
% %     
% %     [sink_tmpLaligned]   = TNC_ExtTrigWinsEAS(tmpSmooth, All_L_trialStarts(rxnTimesL200mv_valid),[2000,2000]); %[-win +win], so mvment onset = ~700
% %     if ~isfinite(sink_tmpLaligned.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
% %         j = j+1;
% %     else
% %     unitNum = num2str(j);
% %     unitname = strcat(' u#', unitNum, 'cue-aligned rxn < 100 ms');
% %    
% %     %Plot the L-mv psth and label the unit plotted
% % %     plot((sink_tmpL.avg - mean(sink_tmpL.avg(1:500))) + (ones(1,numel(sink_tmpL.avg)).*j.*0.01), 'bl'); hold on;  
% % %     alignVar1 = rxnTimesL200mv(rxnTimesL200mv_valid)
% %     fastCueL = All_L_trialStarts(rxnTimesL200mv_valid);
% %     alignVar1 = fastCueL;
% % 
% %     PopData.currParams.stableTrials = -1;                                         % subtracts the first "trial" when calling AlignRaster
% % 
% %     [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session.unit(j).ts, PopData.currParams.stableTrials, alignVar1, psthWin, 1, 1); 
% %     shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'bl');  %blue for left; for 151105 dataset = ipsi
% %     title(unitname); 
% %     end
% %     
% %    % Plot the R-movements (fast):
% %     rxnTimesR200mv_valid = find(rxnTimesR200mv>9.2e4 & rxnTimesR200mv<RtimesMax);
% % 
% %     [sink_tmpRaligned]   = TNC_ExtTrigWinsEAS(tmpSmooth,All_R_trialStarts(rxnTimesR200mv_valid),[2000,2000]); %[-win +win], so mvment onset = ~700
% %     if ~isfinite(sink_tmpRaligned.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
% %         j = j+1;
% %     else
% %     unitNum = num2str(j);
% %     unitname = strcat(' u#', unitNum, 'cue-aligned rxn < 100 ms');
% %    
% % %     alignVar2 = rxnTimesR200mv(rxnTimesR200mv_valid);
% %     fastCueR = All_R_trialStarts(rxnTimesR200mv_valid);
% %     alignVar2 = fastCueR;
% % 
% % 
% %     [respCSS2] = TNC_AlignRasters(tmpSmooth,PopData.session.unit(j).ts, PopData.currParams.stableTrials, alignVar2, psthWin, 1, 1); 
% %     shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthZ, respCSS2.image.psthZe,'r');  %blue for left; for 151105 dataset = ipsi
% % 
% %     %Plot the L-mv psth and label the unit plotted
% % %     plot((sink_tmpR.avg - mean(sink_tmpL.avg(1:500))) + (ones(1,numel(sink_tmpR.avg)).*j.*0.01), 'r'); hold on;  
% %     title(unitname); 
% %     hold off;
% %     end
% %     
% %     %% Plot slow rxn times:
% %         % Plot the L-movements (slow):
% %     subplot(2,2,4); hold on;
% %     rxnTimesL400mv_valid = find(rxnTimesL400mv > 9.2e4 & rxnTimesL400mv< LtimesMax);
% % 
% %     [sink_tmpLaligned]   = TNC_ExtTrigWinsEAS(tmpSmooth,All_L_trialStarts(rxnTimesL400mv_valid),[2000,2000]); %[-win +win], so mvment onset = ~700
% %     if ~isfinite(sink_tmpLaligned.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
% %         j = j+1;
% %     else
% %     unitNum = num2str(j);
% %     unitname = strcat(' u#', unitNum, 'cue-aligned, 400-600 ms');
% %     
% %     slowCueL = All_L_trialStarts(rxnTimesL400mv_valid);
% %     alignVar1 = slowCueL;
% % 
% %     PopData.currParams.stableTrials = -1;                                         % subtracts the first "trial" when calling AlignRaster
% % 
% %     [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session.unit(j).ts, PopData.currParams.stableTrials, alignVar1, psthWin, 1, 1); 
% %     shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'bl');  %blue for left; for 151105 dataset = ipsi
% %     title(unitname); 
% %     end
% %     
% %    % Plot the R-movements (slow):
% %     rxnTimesR400mv_valid = find(rxnTimesR400mv>9.2e4 & rxnTimesR400mv<RtimesMax);
% % 
% %     [sink_tmpRaligned]   = TNC_ExtTrigWinsEAS(tmpSmooth,All_R_trialStarts(rxnTimesR400mv_valid),[2000,2000]); %[-win +win], so mvment onset = ~700
% %     if ~isfinite(sink_tmpRaligned.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
% %         j = j+1;
% %     else
% %     unitNum = num2str(j);
% %     unitname = strcat(' u#', unitNum, 'cue-aligned, 400-600 ms');
% %    
% % %     alignVar2 = rxnTimesR400mv(rxnTimesR400mv_valid);
% %     slowCueR = All_R_trialStarts(rxnTimesR400mv_valid);
% %     alignVar2 = slowCueR;
% % 
% %     [respCSS2] = TNC_AlignRasters(tmpSmooth,PopData.session.unit(j).ts, PopData.currParams.stableTrials, alignVar2, psthWin, 1, 1); 
% %     shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthZ, respCSS2.image.psthZe,'r');  %blue for left; for 151105 dataset = ipsi
% %     
% %     %Save to the structure:
% %    export_mvmt_data_for_psths.rxnTimesL200mv = rxnTimesL200mv;
% %    export_mvmt_data_for_psths.rxnTimesR200mv = rxnTimesR200mv;
% %    export_mvmt_data_for_psths.rxnTimesL400mv = rxnTimesL400mv;
% %    export_mvmt_data_for_psths.rxnTimesR400mv = rxnTimesR400mv;
% %     
% %    export_mvmt_data_for_psths.fastCueL = fastCueL;
% %    export_mvmt_data_for_psths.fastCueR = fastCueR;
% %    export_mvmt_data_for_psths.slowCueL = slowCueL;
% %    export_mvmt_data_for_psths.slowCueR = slowCueR;
% export_mvmt_data_for_psths.Lcue = L_starts(Ltimes_valid); 
% export_mvmt_data_for_psths.Rcue = R_starts(Rtimes_valid); 
% 
% 
% 
% %     title(unitname); 
% %     hold off;
% %     end
% %    %% Save to the structure:
% %    export_mvmt_data_for_psths.rxnTimesL100mv = rxnTimesL200mv;  %These are for corrects and incorrects
% %    export_mvmt_data_for_psths.rxnTimesR100mv = rxnTimesR200mv;
% %    export_mvmt_data_for_psths.rxnTimesL4to600mv = rxnTimesL400mv;
% %    export_mvmt_data_for_psths.rxnTimesR4to600mv = rxnTimesR400mv;
% % 
% %     end
% %     end
% %     end
end
% export_mvmt_data_for_psths.psthMatIpsi = psthMatIpsi;
% export_mvmt_data_for_psths.psthMatContra = psthMatContra;
% 
% aligned1 = 'Ipsi, cue aligned';
% figure; hold on;
% h_ipsi = imagesc(psthMatIpsi);
% sessionName1 = strcat(filestr, aligned1, ',all trials');
% title(sessionName1);
% colorbar;
% % [cmap] = TNC_CreateRBColormap(1024,'gp');
% % colormap(cmap);
% % c_blcks = colorbar;
% % cmin = 0; cmax = 100; caxis([cmin cmax]);
% 
% aligned2 = 'Contra, cue-aligned';
% figure; hold on;
% % subplot(1,2,2); hold on;
% h_contra = imagesc(psthMatContra);
% sessionName2 = strcat(filestr, aligned2, ',all trials');
% title(sessionName2);
% colorbar;
% % [cmap] = TNC_CreateRBColormap(1024,'gp');
% % colormap(cmap);
% % c_blcks = colorbar;
% % cmin = 0; cmax = 100; caxis([cmin cmax]);
% 
% % figure; hold on;
% % h_hardCode = imagesc([64 64 87 95 89 86]); %randomized % corrects
% 
% % [cmap] = TNC_CreateRBColormap(1024,'gp');
% % colormap(cmap);
% % c_rand = colorbar;
% % cmin = 0; cmax = 100; caxis([cmin cmax]);


    %% Save to structure    
    cd(fpath);
     
    save([targetName '_bhAllNS4.mat'],'export_mvmt_data_for_psths');                            
    disp(['saved as ' targetName '_bhAllNS4.mat']);
                                                                     
    disp('%-------------------------------------------------------------------');
    disp(' ');

return
