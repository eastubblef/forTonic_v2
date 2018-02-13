%% Josh's psth code (modified) 2.26.16: 
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/behaveChunks';
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks';

fname = strcat(fpath, '/', '151118_all_pop.mat');
behavFile = strcat(fpath, '/', '151118_all_bh.mat');
load(fname);
pop = PopData;
units = 1:length(pop.session.unit);
varargin = units;

load(behavFile)
behavior = ContData.behavior;

% validsL_first = behavior.validsL_first'; %for L movements
% validsR_first = behavior.validsR_first'; %for R movements

 wheelsLfirstValid = behavior.wheelsLfirstValid;
 wheelsRfirstValid = behavior.wheelsRfirstValid;

 wheelsILfirstValid = behavior.wheelsILfirstValid;
 wheelsIRfirstValid = behavior.wheelsIRfirstValid;


% PopData.currParams.stableTrials = validsL_first;  %will need to update this once I've aligned to actual trials (csv) - relates to which trials to include
PopData.currParams.stableTrials = -1;               %hmm, not sure if this really subtracts the first "trial" when calling AlignRaster

alignVar = wheelsRfirstValid;

%% 
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 15.4;
[currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);

for i = 1:length(PopData.session)
    for j = 1:length(units)
        numStamps = length(PopData.session(i).unit(j).ts);
        delta = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
        delta(1,ceil(PopData.session(i).unit(j).ts)) = 1;
        tmpSmooth = conv(delta,currParams.filter.kernel,'same');

        psthWin = [1.3e3,4.5e3];

%         [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,ContData.behavior.reach.start(rewReaches,4),psthWin,0,1);
        [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts, PopData.currParams.stableTrials, alignVar, psthWin, 0, 1);
%         [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts, PopData.currParams.stableTrials, validsR_first, psthWin, 1, 1);

        figure(j); hold on;
        unitname = num2str(j);
        fnamepre = fname(end-18:end-13);
        dir = 'L mvment onset'
        plotName = strcat(fnamepre, '-', 'u#', unitname, '- ', dir);
        
        plot([-psthWin(1) psthWin(2)],[0 0],'-','Color',[0.2 0.2 0.2]); ylabel('Normalized Firing Rate'); xlabel('ms'); ylim([-.35 0.35]);
        shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'r');  %blue for left; red for right
        title(plotName);
        hold off;
    end
end

%Last line is the key part. The AlignRasters code is looking for:
%[response] = TNC_AlignRasters(continuous spike density function,timestamps for spikes,omit how many trials at beginning?,time stamps to align to,how much time before and adter, want raster plot structured data?, want a boxcar histogram of response?)


% plot the output something like this:

%         plot([-psthWin(1) psthWin(2)],[0 0],'-','Color',[0.2 0.2 0.2]);
%         shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'k');
