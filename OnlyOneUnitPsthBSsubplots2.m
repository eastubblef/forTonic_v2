
function [] = OnlyOneUnitPsthBSsubplots2(fname, behavFile, units, alignVar)
%% This m file calls Josh's code in order to plot psths (as subplots like oneUnitAlignSpike code for rasters plots)
% this is a standalone; modified for Gidon's R01;
% plots the raw firing rate data (i.e., no z-scoring)

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170112/4thPass/behaveSegs';
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks';

%%
% fname = strcat(fpath, '/', '151105_all_pop.mat');                          
% behavFile = strcat(fpath, '/', '151105_all_bh.mat');
% fname = strcat(fpath, '/', '151105_untagged_pop.mat');                          
% behavFile = strcat(fpath, '/', '151105_untagged_bh.mat');

% fname = strcat(fpath, '/', '151105_1001laser_pop.mat');                          
% behavFile = strcat(fpath, '/', '151105_1001laser_bh.mat');

% behavFile = strcat(fpath, '/', '170112behave_bh.mat');
% popFile = strcat(fpath, '/', '170112behave_pop.mat');
behavFile = strcat(fpath, '/', '151118_all_bh.mat');
popFile = strcat(fpath, '/', '151118_all_pop.mat');

fname = popFile;
load(fname);
pop = PopData;
% units = 1:length(pop.session.unit);
units = pop.session.unit(1,:).ts;  
% units = ans;  

varargin = units;

load(behavFile)
behavior = ContData.behavior;

% LrewTrialStarts = behavior.LrewTrialStarts;
% RrewTrialStarts = behavior.RrewTrialStarts;
 wheelsLfirstValid = behavior.wheelsLfirstValid;
 wheelsRfirstValid = behavior.wheelsRfirstValid;

%  wheelsILfirstValid = behavior.wheelsILfirstValid;
%  wheelsIRfirstValid = behavior.wheelsIRfirstValid;
% rewTimes = behavior.rewTimes;

alignVar1 = wheelsRfirstValid;    %151119: R = ipsi (dotted) for Gidon
alignVar2 = wheelsLfirstValid;    %151119: L = contra (solid) for Gidon


% alignVar1 = wheelsLfirstValid;    %151105: L = ipsi (dotted) for Gidon
% alignVar2 = wheelsRfirstValid;    %151105: R = contra (solid) for Gidon
% alignVar1 = wheelsILfirstValid;
% alignVar2 = wheelsIRfirstValid;
% alignVar1 = rewTimesValid;
% alignVar1 = LrewTrialStarts;
% alignVar2 = RrewTrialStarts;
% alignVar2 = rewTimes; 
% alignVar1 = wheelsILfirstValid;
% alignVar2 = wheelsIRfirstValid;
oneUnit = PopData.session.unit(3).ts;

PopData.currParams.stableTrials = -1;                                         % subtracts the first "trial" when calling AlignRaster

%% 
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 15.4;

[currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);

% for i = 1:length(PopData.session)
    figure; hold on; 
    dimm=ceil(sqrt(length(oneUnit))); %

%     for j = 1:length(units)
        numStamps = length(oneUnit);
        delta = zeros(1,ceil(oneUnit(numStamps)));
        delta(1,(oneUnit)) = 1;
        tmpSmooth = conv(delta,currParams.filter.kernel,'same');

%         psthWin = [.5e3,.5e3]; % best for 151105

        psthWin = [1.0e3,2e3]; % best for 151105
%         psthWin = [1.0e3,1.5e3]; % best for 151105 (Gidon wants Window[1,1] - but doesn't return to baseline)
        
%         subplot(dimm,dimm,j); hold on;        

 %% Call to AlignRasters:
 
        % FIRST plot for ipsi (w/ respect to hemisphere recording) trials
        [respCSS] = TNC_notZscoredAlignRasters(tmpSmooth,oneUnit, PopData.currParams.stableTrials, alignVar1, psthWin, 0, 0); 
%         unitname = num2str(j);
%         dir = plotname
%         plotName = strcat('u#', unitname, '- ', dir);
%             h = plot([-(psthWin(1)) psthWin(2)],[0 0]); 
%             set(h, 'Color',[0.2 0.2 0.2]);

        ylabel('Firing rate', 'FontSize', 22);
        xlabel('ms', 'FontSize', 22);
        ax = gca; 
        ax.FontSize = 22;

            AVG = respCSS.image.psthAVG;
            SEM = respCSS.image.psthSEM;
%             ipsiH = shadedErrorBar(-psthWin(1):psthWin(2), AVG, SEM,'k--');  %blue for left; red for right
            ipsiH = shadedErrorBar(-psthWin(1):psthWin(2), AVG, SEM,'bl--');  %blue for left; red for right

%           set(ipsiH.mainLine, 'Color', [.5, .5, .5]);     

        % SECOND plot (for contra trials)
        [respCSS2] = TNC_notZscoredAlignRasters(tmpSmooth,oneUnit, PopData.currParams.stableTrials, alignVar2, psthWin, 0, 0);
%         shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthAVG, respCSS2.image.psthSEM,'k');  %blue for left; red for right
        shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthAVG, respCSS2.image.psthSEM,'r');  %blue for left; red for right

      ylim([0 0.04]); %best for 151105 raw data (not z-scored)
        
 
     numbins = length(AVG);
     binsize = (psthWin(1) + psthWin(2))/numbins; % this is in seconds; .999 = 1ms
        
    end
% end

%Last line is the key part. The AlignRasters code is looking for:
%[response] = TNC_AlignRasters(continuous spike density function,timestamps for spikes,omit how many trials at beginning?,time stamps to align to,how much time before and adter, want raster plot structured data?, want a boxcar histogram of response?)

% plot the output something like this:

%         plot([-psthWin(1) psthWin(2)],[0 0],'-','Color',[0.2 0.2 0.2]);
%         shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'k');
