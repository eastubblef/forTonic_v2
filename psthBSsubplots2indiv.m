
function [] = psthBSsubplots2indiv(fname, behavFile, units, alignVar)
%% This m file calls Josh's code in order to plot psths (as subplots like oneUnitAlignSpike code for rasters plots)
% calls the following mfiles:
% TNC_createGaussian.m - explanatory
% TNC_alignRasters.m -  for creating inputted time window (line 81), setting the inputted zero-alignment of pop data (units) to behavioral data (behavFile)
% corresponding to ipsi v contra "alignVars 1 & 2" (lines 51 & 52) 
% shadedErrorBar.m -  for plotting 

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better/laserChunks/untagged';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151104';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151104/1001/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151106/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170112/1stPass/behaveSegs1'
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/151119good/behaveChunks/tagged';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/160505';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better/laserChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/161005/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/161015';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/160505';
fpath = '/Volumes/My Passport for Mac/170111/newBehaveUnits';
%%
% fname = strcat(fpath, '/', '170112behave_pop.mat');                          
% behavFile = strcat(fpath, '/', '170112behave_bh.mat');

% fname = strcat(fpath, '/', '160505New2_pop.mat');                          
% behavFile = strcat(fpath, '/', '160505New2_bh.mat');
% fname = strcat(fpath, '/', '160505New2_pop.mat');                          
% behavFile = strcat(fpath, '/', '160505New2_bh.mat');

fname = strcat(fpath, '/', '170111newBehave_pop.mat');                          
behavFile = strcat(fpath, '/', '170111_bh3.mat');

load(fname);
pop = PopData;

% allow for plotting only 1 of them
%PopData.session.unit(4).ts  
% units1 = pop.session.unit(6).ts; %151105
units1 = pop.session.unit(12).ts; %170111


% units2 = pop.session.unit(5).ts; 
% units3 = pop.session.unit(17).ts;
% units = 1:length(units1);
% units = 1:length(pop.session.unit);    % plots all of them

units = units1;   %update here

varargin = units; 

load(behavFile)
% behavior = ContData.behavior;
behavior = export_mvmt_data_for_psths

% LrewTrialStarts = behavior.LrewTrialStarts;
% RrewTrialStarts = behavior.RrewTrialStarts;
% LincorrectTrialStarts = behavior.LincorrectTrialStarts;
% RincorrectTrialStarts = behavior.RincorrectTrialStarts;
%  wheelsLfirstValid = behavior.wheelsLfirstValid;
%  wheelsRfirstValid = behavior.wheelsRfirstValid;
 
wheelsLfirstValid = behavior.wheelALL_LfirstValid
wheelsRfirstValid = behavior.wheelALL_RfirstValid

%  wheelsILfirstValid = behavior.wheelsILfirstValid;
%  wheelsIRfirstValid = behavior.wheelsIRfirstValid;
% rewTimes = behavior.rewTimes;

% wheelsAll_R = behavior.wheelsAll_R;
% wheelsAll_L = behavior.wheelsAll_L;

% alignVar1 = wheelsLfirstValid; %ipsi for the 151105 dataset (plots blue)
% alignVar2 = wheelsRfirstValid; %contra for the 151105 dataset (plots red)

alignVar1 = wheelsLfirstValid;  %bl = ipsi
alignVar2 = wheelsRfirstValid;  %r = contra
% alignVar1 = rewTimes;
% alignVar1 = LrewTrialStarts;
% alignVar2 = RrewTrialStarts;
% alignVar1 = wheelsILfirstValid;
% alignVar2 = wheelsIRfirstValid;
% alignVar1 = LincorrectTrialStarts;
% alignVar2 = RincorrectTrialStarts;
% alignVar1 = wheelsAll_L;
% alignVar2 = wheelsAll_R;

PopData.currParams.stableTrials = -1;                                         % subtracts the first "trial" when calling AlignRaster

%% 
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 15.4;
[currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);

% for i = 1:length(PopData.session)
    figure; hold on; 
    plotname = 'rew';
    dimm=ceil(sqrt(length(units))); %

%     for j = 1:length(units)
%         numStamps = length(PopData.session.unit(j).ts);
        numStamps = length(units);
        delta = zeros(1,ceil(units(numStamps)));
        delta(1,ceil(units)) = 1;
        tmpSmooth = conv(delta,currParams.filter.kernel,'same');

%         psthWin = [1.50e3,4.5e3];
        psthWin = [1.0e3,2e3]; % best for 151105 data
%         psthWin = [.5e3, .5e3];
        
%         subplot(dimm,dimm,j); hold on; 
                         
        % second plot
        [respCSS2] = TNC_AlignRasters(tmpSmooth,units, PopData.currParams.stableTrials, alignVar1, psthWin, 1, 0);
%             h(j) = plot([-(psthWin(1)) psthWin(2)],[0 0]); 
%             set(h(j), 'Color',[0.2 0.2 0.2]);
%             ylabel('Firing Rate'); xlabel('ms'); 
        shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthZ, respCSS2.image.psthZe,'bl');  %red for right; for 151105 data = contra
                
        % first plot
        [respCSS] = TNC_AlignRasters(tmpSmooth,units, PopData.currParams.stableTrials, alignVar2, psthWin, 1, 0); 
        unitname = units;
        dir = plotname
%         plotName = strcat('u#', unitname, '- ', dir);
            h = plot([-(psthWin(1)) psthWin(2)], [0 0]); 
             set(h, 'Color',[0.2 0.2 0.2]);
            ylabel('Firing Rate'); xlabel('ms'); 
        shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'r');  %blue for left; for 151105 dataset = ipsi
% 
%         % second plot
%         [respCSS2] = TNC_AlignRasters(tmpSmooth,units, PopData.currParams.stableTrials, alignVar2, psthWin, 1, 0);
% %             h(j) = plot([-(psthWin(1)) psthWin(2)],[0 0]); 
% %             set(h(j), 'Color',[0.2 0.2 0.2]);
% %             ylabel('Firing Rate'); xlabel('ms'); 
%         shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthZ, respCSS2.image.psthZe,'r');  %red for right; for 151105 data = contra
% 
        title('160505u6, rew');
%         hold off;
%     end
% end

%Last line is the key part. The AlignRasters code is looking for:
%[response] = TNC_AlignRasters(continuous spike density function,timestamps for spikes,omit how many trials at beginning?,time stamps to align to,how much time before and adter, want raster plot structured data?, want a boxcar histogram of response?)


% plot the output something like this:

%         plot([-psthWin(1) psthWin(2)],[0 0],'-','Color',[0.2 0.2 0.2]);
%         shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'k');
