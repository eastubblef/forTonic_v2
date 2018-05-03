
%% This script is for aligning on one variable & sorting on the other
% call the bh file of interest
% 4.13.18

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/171013';
% filestr = '171013002.ns4';
% % behavFile = strcat(fpath, '/', '171013_bh3Rxn.mat');
% % behavFile = strcat(fpath, '/', '171013_bh5Rxn.mat');       %movement align only; 2 sec bline
% % behavFile = strcat(fpath, '/', '171013_bhCue.mat');
% % behavFile = strcat(fpath, '/', '171013_bhCue0.mat');       %movement align only; 2 sec bline w/ out convolution
% behavFile = strcat(fpath, '/', '171013_bhCueNS4');           %Works! Best for future data

fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170111';
filestr = '170111002.ns4';
% behavFile = strcat(fpath, '/', '170111_bh5Rxn.mat'); %works without offset
behavFile = strcat(fpath, '/', '170111_bh6RxnSmooth.mat'); %works without offset
S = load('170111newBehave_pop');
LtimesMax = 3.0e6;
RtimesMax = 3.0e6;


% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170118/behaveNewUnits';
% filestr = '170118002.ns4';
% % behavFile = strcat(fpath, '/', '170118_bh5Rxn.mat'); %use this to est. wheel offset
% behavFile = strcat(fpath, '/', '170118_bh6RxnSmooth.mat'); %works without offset
% Note that ipsi = R trials; contra = L

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170112/behaveNew';
% filestr = '170112002.ns4';
% behavFile = strcat(fpath, '/', '170112_bhCue.mat');  %works without offset, but neural activity fucked

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/171012/newBehaveUnits';
% filestr = '171012001.ns4';
% % behavFile = strcat(fpath, '/', '171013_bh3Rxn.mat');
% % behavFile = strcat(fpath, '/', '171012_bhCue.mat');
% behavFile = strcat(fpath, '/', '171012_bh6RxnSmooth.mat');

% % 
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/newBehaveUnits/behaveChunks';
% % cd fpath;
% % fname = strcat(fpath, '/', '151105_all_pop.mat');                          
% % behavFile = strcat(fpath, '/', '151105_all_bh.mat');
% fname = strcat(fpath, '/', '151105_update_pop.mat');                          
% behavFile = strcat(fpath, '/', '151105_update_bh.mat');
% filestr = '151105tag1001.ns4';

%%  MUST LOAD THE openNSx.m file from 2016 for this to work (currently in forTonic_v2 folder, mac HD)

filestr2 = strcat(fpath, '/', filestr);
load(behavFile)
[data] = openNSx(filestr2);              

PopData = S.PopData;
%2017 aligned data - choose the first event to align to
% wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;     %all correct mvments
% wheelsRfirstValid = export_mvmt_data_for_psths.wheelRfirstValid;

All_L_trialStarts = export_mvmt_data_for_psths.All_L_trialStarts;    %for 170111 & 170118
All_R_trialStarts = export_mvmt_data_for_psths.All_R_trialStarts;
eventTimesIpsi = All_L_trialStarts;
eventTimesContra = All_R_trialStarts;

% LrewTrueTrialStarts = export_mvmt_data_for_psths.LrewTrueTrialStarts; %only trial starts that lead to correct trials
% RrewTrueTrialStarts = export_mvmt_data_for_psths.RrewTrueTrialStarts;

% Lcue = export_mvmt_data_for_psths.Lcue;   %For all trial starts: 171012, 171013
% Rcue = export_mvmt_data_for_psths.Rcue;   
% eventTimesIpsi = Lcue;
% eventTimesContra = Rcue;

%% choose event to sort on
wheelALL_LfirstValid = export_mvmt_data_for_psths.wheelALL_LfirstValid;
wheelALL_RfirstValid = export_mvmt_data_for_psths.wheelALL_RfirstValid;

sortTimesIpsi = wheelALL_LfirstValid;
sortTimesContra = wheelALL_RfirstValid;


Sigma = 24;
t = 1:1:Sigma.*30;
Gaussian = ( 1./( sqrt(2.*pi.*Sigma.*Sigma) ) ) .* exp( -(t-Sigma.*15).^2 ./ (2.*Sigma).^2 );
integral = trapz(Gaussian);                         %trapz(y) returns the approximate integral of y
Gaussian = Gaussian./integral;
[mapName] = TNC_CreateRBColormap(8,'mbr');          %calls color map function

% Walk thru ea. unit & build matrix of aligned psths
numUnits = numel( PopData.session(1).unit );
% figure;     

psthWin = [2.0e3,2e3];
psthMat = NaN(length(-psthWin(1):psthWin(2)), length(1:numUnits));

for j=1:numUnits  
    tmp = PopData.session(1).unit(j).ts;
    delta = zeros(1,ceil(max(tmp)));
    delta(round(tmp)) = 1;

    tmpSmooth = conv(delta,Gaussian,'same');
    
%% L-movement trial start cue:
    figure; hold on;
    subplot(2,2,1); hold on;

    Ltimes_valid = find(All_L_trialStarts > 9.2e4 & All_L_trialStarts < LtimesMax);

    [vals, inds] = sort(magnitude(Ltimes_valid));
    [sink_tmpL]   = TNC_ExtTrigWinsEAS(tmpSmooth,Ltimes_valid,[2000,2000]); %[-win +win], so mvment onset = ~700
%     if ~isfinite(sink_tmpL.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
%         j = j+1;
%     else
%     unitNum = num2str(j);
%     unitname = strcat(' u#', unitNum, 'cue-aligned, all');
   
    alignVar(j) = All_L_trialStarts;
    PopData.currParams.stableTrials = -1;          %subtracts the first "trial" when calling AlignRaster

    [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session.unit(j).ts, PopData.currParams.stableTrials, alignVar(j), psthWin, 1, 1); 
%     shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'bl');  %blue for left; for 151105 dataset = ipsi
%     title(unitname); 
    psthMat(j,:) = respCSS.image.psthZ(j)
    
end
    
%     %% Plot the R-movements' trial start cue:
%     Rtimes_valid = find(R_starts >9.2e4 & R_starts < RtimesMax);
%     [vals, inds] = sort(magnitude(Rtimes_valid));
% 
%     [sink_tmpR]   = TNC_ExtTrigWinsEAS(tmpSmooth,R_starts(Rtimes_valid),[2000,2000]); %[-win +win], so mvment onset = ~700
%     if ~isfinite(sink_tmpR.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
%         j = j+1;
%     else
%     unitNum = num2str(j);
%     unitname = strcat(' u#', unitNum, 'cue-aligned,  all');
%    
%     alignVar2 = R_starts(Rtimes_valid);
%     [respCSS2] = TNC_AlignRasters(tmpSmooth,PopData.session.unit(j).ts, PopData.currParams.stableTrials, alignVar2, psthWin, 1, 1); 
%     shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthZ, respCSS2.image.psthZe,'r');  %blue for left; for 151105 dataset = ipsi
%     title(unitname); 
%     hold off;
%     end    

%% Plot the trial-aligned data:
    % Plot all L-movements' trial starts:
    subplot(2,2,2); hold on;
%     wheelsLfirstValid_valid = find(wheelsLfirstValid>9.2e4 & wheelsLfirstValid<LtimesMax);
    LeftsTimes = export_mvmt_data_for_psths.LeftsTimes;
    Lmoves_valid = find(LeftsTimes > 9.2e4 & LeftsTimes<LtimesMax);
    
    [sink_tmpLaligned]   = TNC_ExtTrigWinsEAS(tmpSmooth,LeftsTimes(Lmoves_valid),[2000,2000]); %[-win +win], so mvment onset = ~700
    if ~isfinite(sink_tmpLaligned.wins(1))               %EAS added in case a unit's numel indices ~= numEvents
        j = j+1;
    else
    unitNum = num2str(j);
    unitname = strcat(' u#', unitNum, 'mv-aligned, all');
    alignVar1 = LeftsTimes(Lmoves_valid);

    PopData.currParams.stableTrials = -1;        
%     if j == 20  %for the 170112 dataset
%         break
%         else
        [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session.unit(j).ts, PopData.currParams.stableTrials, alignVar1, psthWin, 1, 1); 
        shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'bl');  %blue for left; for 151105 dataset = ipsi
        title(unitname); 
%         end     %for the 170112 dataset
    end
    %% Plot the R-movements:

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



imagesc(eventTimesIpsi);

