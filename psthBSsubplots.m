%% Josh's psth code (modified) 2.26.16: 
function [] = psthBSsubplots(fname, behavFile, alignVar)
%% This m file calls Josh's code in order to plot psths (as subplots like oneUnitAlignSpike code for rasters plots

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/laserChunks/untagged';
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/behaveChunks';
fname = strcat(fpath, '/', '151105_all_pop.mat');
behavFile = strcat(fpath, '/', '151105_all_bh.mat');

% fname = strcat(fpath, '/', '151105_untagged_pop.mat');                           
% behavFile = strcat(fpath, '/', '151105_untagged_bh.mat');

load(fname);
pop = PopData;
units = 1:length(pop.session.unit);
varargin = units;
 
load(behavFile)
behavior = ContData.behavior;

% LrewTrialStarts = behavior.LrewTrialStarts;
% RrewTrialStarts = behavior.RrewTrialStarts;
%  wheelsLfirst = behavior.wheelsLfirst;
%  wheelsRfirst = behavior.wheelsRfirst;

 wheelsLfirstValid = behavior.wheelsLfirstValid;
 wheelsRfirstValid = behavior.wheelsRfirstValid;
 wheelsILfirstValid = behavior.wheelsILfirstValid;
 wheelsIRfirstValid = behavior.wheelsIRfirstValid;
 rewTimesValid = behavior.rewTimesValid;

% alignVar = wheelsLfirst;
%  alignVar = wheelsLfirstValid;
% alignVar = wheelsRfirstValid;
% alignVar = wheelsILfirstValid;
% alignVar = wheelsIRfirstValid;
alignVar = rewTimesValid;

PopData.currParams.stableTrials = -1;                                         % subtracts the first "trial" when calling AlignRaster

%% 
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 15.4;
[currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);

for i = 1:length(PopData.session)
    figure; hold on; 
    plotname = 'solenoid clicks'
    dimm=ceil(sqrt(length(units))); %

    for j = 1:length(units)
        numStamps = length(PopData.session(i).unit(j).ts);
        delta = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
        delta(1,ceil(PopData.session(i).unit(j).ts)) = 1;
        tmpSmooth = conv(delta,currParams.filter.kernel,'same');

%         psthWin = [1.50e3,4.5e3];
        psthWin = [1.0e3,2e3];
        
        subplot(dimm,dimm,j); hold on;

%         [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,ContData.behavior.reach.start(rewReaches,4),psthWin,0,1);
        [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts, PopData.currParams.stableTrials, alignVar, psthWin, 1, 1);
%         [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts, PopData.currParams.stableTrials, validsR_first, psthWin, 1, 1);

        unitname = num2str(j);
        dir = plotname
        plotName = strcat('u#', unitname, '- ', dir);
     
        h(j) = plot([-(psthWin(1)) psthWin(2)],[0 0]); 
        set(h(j), 'Color',[0.2 0.2 0.2]);
        ylabel('Firing Rate'); xlabel('ms'); 
                
        shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'k');  %blue for left; red for right

        title(plotName);
        hold off;
    end
end

%Last line is the key part. The AlignRasters code is looking for:
%[response] = TNC_AlignRasters(continuous spike density function,timestamps for spikes,omit how many trials at beginning?,time stamps to align to,how much time before and adter, want raster plot structured data?, want a boxcar histogram of response?)


% plot the output something like this:

%         plot([-psthWin(1) psthWin(2)],[0 0],'-','Color',[0.2 0.2 0.2]);
%         shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'k');
