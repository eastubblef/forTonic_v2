
function [alignVarMat] = OnlyOnePsthBSsubplots4concatIpsiL(PopData, behavFile, units, alignVar1, alignVar2)
%% This m file calls Josh's code in order to plot psths (as subplots like oneUnitAlignSpike code for rasters plots)
% this is a standalone; modified for Gidon's R01;
% plots the raw firing rate data (i.e., no z-scoring)

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151106/behave';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better/laserChunks/untagged';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/behave/tag';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170112/shanks3&4/1stPass/behaveSegs1';
%%

% load(fname);
% pop = PopData;
% units = 1:length(pop.session.unit);
% load(behavFile)
% behavior = ContData.behavior;

% FOR 151118-151119:
% alignVar1 = wheelsRfirstValid;    %151119: R = ipsi (dotted) for Gidon
% alignVar2 = wheelsLfirstValid;    %151119: L = contra (solid) for Gidon

% FOR 151104-151106:
% alignVar1 = validsL_first;    %151105: L = ipsi (dotted) for Gidon
% alignVar2 = validsR_first;    %151105: R = contra (solid) for Gidon
% alignVar1 = wheelsRfirstValid;    %151105: L = ipsi (dotted) for Gidon
% alignVar2 = wheelsLfirstValid;    %151105: R = contra (solid) for Gidon

oneUnit = PopData.session.unit(units).ts;  %for 151118 fig
% oneUnit = PopData.session.unit(1).ts;  %for 151105 fig
% oneUnit = Pop.session.unit(4).ts;  %for 151118 fig

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
        
%         subplot(dimm,dimm,j); hold on;        

 %% Call to AlignRasters:
 
        % FIRST plot for ipsi (w/ respect to hemisphere recording) trials
%         [respCSS] = TNC_notZscoredAlignRasters(tmpSmooth,oneUnit, PopData.currParams.stableTrials, alignVar1, psthWin, 0, 0); 
       [respCSS] = TNC_AlignRasters(tmpSmooth,oneUnit, PopData.currParams.stableTrials, alignVar1, psthWin, 1, 1); 
        
%         unitname = num2str(j);
%         dir = plotname
%         plotName = strcat('u#', unitname, '- ', dir);
%         ipsiH = plot([-(psthWin(1)) psthWin(2)],[0 0]); 

        AVG = respCSS.image.psthAVG;
%         AVG = respCSS.image.psthAVG * 100;
        SEM = respCSS.image.psthSEM;
%         ipsiH = shadedErrorBar(-psthWin(1):psthWin(2), AVG, SEM,'k--');  %blue for left; red for right
        ipsiH =  shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'bl-');
        set(ipsiH.mainLine, 'lineWidth', 2);     

% SECOND plot (for contra trials)
%        [respCSS2] = TNC_notZscoredAlignRasters(tmpSmooth,oneUnit, PopData.currParams.stableTrials, alignVar2, psthWin, 0, 0);
       [respCSS2] = TNC_AlignRasters(tmpSmooth, oneUnit, PopData.currParams.stableTrials, alignVar2, psthWin, 1, 1);      
%        AVG2 = respCSS2.image.psthAVG * 100;
%        SEM2 = respCSS2.image.psthSEM;
       contraH =  shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthZ, respCSS2.image.psthZe,'r');
%        contraH = shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthAVG, respCSS2.image.psthSEM,'k');  %blue for left; red for right
%        shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthZ, respCSS2.image.psthZe,'k');  %blue for left; red for right
%        set(contraH.mainLine, 'Color', [0, 0, 0], 'lineWidth', 1.5);     
       set(contraH.mainLine, 'lineWidth', 1.5);     

 
     numbins = length(AVG);
     binsize = (psthWin(1) + psthWin(2))/numbins; % this is in seconds; .999 = 1ms
        
%         ylabel('Norm. firing rate', 'FontSize', 28);
%         ylim([-0.2 1.2]);
        xlabel('ms', 'FontSize', 28);
        ax = gca; 
        ax.FontSize = 28;

        alignVarMat = NaN(1000, 1000);
        alignVar1 = alignVar1';
        alignVar2 = alignVar2';
            alignVarMat(:,1) = alignVar1;
            alignVarMat(:,2) = alignVar2;

end

%Last line is the key part. The AlignRasters code is looking for:
%[response] = TNC_AlignRasters(continuous spike density function,timestamps for spikes,omit how many trials at beginning?,time stamps to align to,how much time before and adter, want raster plot structured data?, want a boxcar histogram of response?)

% plot the output something like this:

%         plot([-psthWin(1) psthWin(2)],[0 0],'-','Color',[0.2 0.2 0.2]);
%         shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'k');
