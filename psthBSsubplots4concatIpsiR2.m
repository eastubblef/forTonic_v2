
function [psthMat, psthMatStr] = psthBSsubplots4concatIpsiR2(PopData, behavFileStr, units, alignVar1, alignVar2, alignVar1str, alignVar2str, taggedNeurons)
%% This m file calls Josh's code in order to plot psths (as subplots like oneUnitAlignSpike code for rasters plots)
% calls the following mfiles:
% TNC_createGaussian.m - explanatory
% TNC_alignRasters.m -  for creating inputted time window (line 81), setting the inputted zero-alignment of pop data (units) to behavioral data (behavFile)
% corresponding to ipsi v contra "alignVars 1 & 2" (lines 51 & 52) 
% shadedErrorBar.m -  for plotting 

% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/behave/tag';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151104/1001/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151106/behaveChunks';
% fpath = '/Users/stubblefielde/Desktop/behaviorFebMarApr16/Vgatfour/forAnalysis'
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/160505';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better/safeKeepingLaserLast';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170112';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170321';

%%
% fname = strcat(fpath, '/', '151119_tag_pop.mat');
% behavFile = strcat(fpath, '/', '151119_tag_bh.mat');
% fname = strcat(fpath, '/', '151118_all_pop.mat');
% behavFile = strcat(fpath, '/', '151118_all_bh.mat');
% fname = strcat(fpath, '/', '151106_all_pop.mat');                          
% behavFile = strcat(fpath, '/', '151106_all_bh.mat');
% fname = strcat(fpath, '/', '161015behave_pop.mat');                          
% behavFile = strcat(fpath, '/', '161015behave_bh.mat');

%updated behaveLoad6; re-ran MoverScript -> TNC_MoveBehaveExtractNew
% fname = strcat(fpath, '/', '151105_update_pop.mat');                       
% behavFile = strcat(fpath, '/', '151105_update_bh.mat');
% fname = strcat(fpath, '/', '151119_tag_pop.mat');                       
% behavFile = strcat(fpath, '/', '151119_tag_bh.mat');
% behavFile = strcat(fpath, '/', '151119_tagUpdate_bh.mat');
% fname = strcat(fpath, '/', '151105chunk26_1001_3_pop.mat');                       
% behavFile = strcat(fpath, '/', '151105chunk26_1001_3_bh.mat');

% behaveFile = strcat(fpath, '/', '170118behave_bh.mat');
% popFile = strcat(fpath, '/', '170118behave_pop.mat');

% behavFile = strcat(fpath, '/', '170321behave_bh.mat');
% popFile = strcat(fpath, '/', '170321behave_pop.mat');

% behavFile = strcat(fpath, '/', '170112NevBoss_bh.mat');
% popFile = strcat(fpath, '/', '170112NevBoss_pop.mat');
% behavFile = strcat(fpath, '/', '170112behave2_bh.mat');
% popFile = strcat(fpath, '/', '170112behave2_pop.mat');
% behavFile = strcat(fpath, '/', '151119_tag_bh.mat');
% popFile = strcat(fpath, '/', '151119_tag_pop.mat');
% 
% 
% load(popFile);
% pop = PopData;
% 
% % allow for plotting only 1 of them
% units = pop.session.unit;  %for Tonic data
% 
% % units1 = unit(5);
% units = 1:length(units1);
units = 1:length(PopData.session.unit);    % plots all of them
varargin = units;

% load(behavFile)
% behavior = ContData.behavior;

% trialStarts = behavior.trialStarts;
% LrewTrialStarts = behavior.LrewTrialStarts;
% RrewTrialStarts = behavior.RrewTrialStarts;
% % LincorrectTrialStarts = behavior.LincorrectTrialStarts;
% % RincorrectTrialStarts = behavior.RincorrectTrialStarts;
% wheelsLfirstValid = behavior.wheelsLfirstValid;
% wheelsRfirstValid = behavior.wheelsRfirstValid;
% 
% wheelsILfirstValid = behavior.wheelsILfirstValid;
% wheelsIRfirstValid = behavior.wheelsIRfirstValid;
% rewTimes = behavior.rewTimes;
% % 
% wheelsAll_R = behavior.wheelsAll_R;
% wheelsAll_L = behavior.wheelsAll_L;

% %If using wheel data for alignment:
% wheel_L_first = behavior.wheel_L_first;
% wheel_R_first = behavior.wheel_R_first;
% 

% alignVar1 = behavior.wheelsRfirstValid;    %151119: R = ipsi (dotted) for Gidon
% alignVar2 = behavior.wheelsLfirstValid;    %151119: L = contra (solid) for Gidon

% alignVar1 = trialStarts;
% alignVar1 = wheel_L_first;
% alignVar2 = wheel_R_first;

PopData.currParams.stableTrials = -1;                                         % subtracts the first "trial" when calling AlignRaster

%% 
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 15.4;
[currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);

% alignVarMat = NaN(1000, 1000);  %added this section 11.2017 for building a matrix for dim reduction
psthMat = [];
psthMatStr = {}; %fill this with strings for indexing later on

for i = 1:length(PopData.session)
    figure; hold on; 
%     plotname = '151119moveOn';
    plotname = char(behavFileStr);

    dimm=ceil(sqrt(length(units))); %

    for j = 1:length(units)
        numStamps = length(PopData.session(i).unit(j).ts);
        delta = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
        delta(1,ceil(PopData.session(i).unit(j).ts)) = 1;
%         delta(1) = 1;
%         x = 1,ceil(PopData.session(i).unit(j).ts);
%         delta(x) = 1;
        tmpSmooth = conv(delta,currParams.filter.kernel,'same');
        psthWin = [1.0e3,2e3]; % best for 151105 data        
        subplot(dimm,dimm,j); hold on; 
                        
%       first plot
        [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts, PopData.currParams.stableTrials, alignVar1, psthWin, 1, 1); 
        
        unitname = num2str(j);
        shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'bl');  %blue for left; for 151105 dataset = ipsi
%         
        %second plot
        [respCSS2] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts, PopData.currParams.stableTrials, alignVar2, psthWin, 1, 1);
            h(j) = plot([-(psthWin(1)) psthWin(2)],[0 0]); 
            set(h(j), 'Color',[0.2 0.2 0.2]);
            ylabel('Firing Rate'); xlabel('ms'); 
        shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthZ, respCSS2.image.psthZe,'r');  %red for right; for 151105 data = contra
%  
        title(plotname);
        
%         psthAvg1 = respCSS.image.psthAVG';
%         psthAvg2 = respCSS2.image.psthAVG';
        psthAvg1 = respCSS.image.psthZ';
        psthAvg2 = respCSS2.image.psthZ';
    
%         for z = 1:length(psthAvg1)
%             if j == 1
%                 psthMat(z,j) = psthAvg1(z);
%                  psthMatStr{1,j} = char(behavFileStr(1:6));
%                  psthMatStr{2,j} = char(strcat('u', num2str(unitname)));
%                  psthMatStr{3,j} = char(alignVar1str);
% 
%             else
%                 psthAvg1 = respCSS.image.psthZ';
%                 psthMat(z,j+(j-1)) = psthAvg1(z);
%                  psthMatStr{1,j+(j-1)} = char(behavFileStr(1:6));
%                  psthMatStr{2,j+(j-1)} = char(strcat('u', num2str(unitname)));
%                  psthMatStr{3,j+(j-1)} = char(alignVar1str);
%             end
%         end
%         for zz = 1:length(psthAvg2)
%             if j == 1
%                  psthMat(zz,j*2) = psthAvg2(zz); 
%                  psthMatStr{1,j*2} = char(behavFileStr(1:6));
%                  psthMatStr{2,j*2} = char(strcat('u', num2str(unitname)));
%                  psthMatStr{3,j*2} = char(alignVar2str);
% 
%             else 
%                 psthAvg2 = respCSS2.image.psthZ';
%                 psthMat(zz,j*2) = psthAvg2(zz);
%                 psthMatStr{1,j*2} = char(behavFileStr(1:6));
%                 psthMatStr{2,j*2} = char(strcat('u', num2str(unitname)));
%                 psthMatStr{3,j*2} = char(alignVar2str);
% 
%             end
%         end
% %         clear psthAvg1; clear psthAvg2;
%     end
% end
        for z = 1:length(psthAvg1)
            if j == 1
                psthMat(z,j) = psthAvg1(z);
                psthMatStr{1,j} = char(behavFileStr(1:6));
                psthMatStr{2,j} = char(strcat('u', num2str(unitname)));
                psthMatStr{3,j} = char(alignVar1str);
                psthMatStr{4,j} = char('untagged');
                psthMatStr{5,j} = char('rand');
                
            else
                psthAvg1 = respCSS.image.psthZ';
                psthMat(z,j+(j-1)) = psthAvg1(z);
                psthMatStr{1,j+(j-1)} = char(behavFileStr(1:6));
                psthMatStr{2,j+(j-1)} = char(strcat('u', num2str(unitname)));
                psthMatStr{3,j+(j-1)} = char(alignVar1str);
                
                if j == taggedNeurons(1) || j == taggedNeurons(2) || j == taggedNeurons(3)
                    psthMatStr{4,j+(j-1)} = char('tag');
                else
                    psthMatStr{4,j+(j-1)} = char('untagged');
                end
                psthMatStr{5,j} = char('rand');
            end
        end
        for zz = 1:length(psthAvg2)
            if j == 1
                psthMat(zz,j*2) = psthAvg2(zz);
                %                  psthMatStr{1,j*2} = char(behavFileStr(1:6));
                psthMatStr{1,j*2} = char(behavFileStr(1:6));
                psthMatStr{2,j*2} = char(strcat('u', num2str(unitname)));
                psthMatStr{3,j*2} = char(alignVar2str);
                psthMatStr{4,j*2} = char('untagged');
                %                  psthMatStr{5,j*2} = char('tag');
                psthMatStr{5,j*2} = char('rand');
            else
                %                 psthAvg2 = respCSS2.image.psthAVG';
                psthAvg2 = respCSS2.image.psthZ';
                psthMat(zz,j*2) = psthAvg2(zz);           %This should make col 1 = alignVar1; col2 = alignVar2 for ea. unit to build psthMat
                psthMatStr{1,j*2} = char(behavFileStr(1:6));
                psthMatStr{2,j*2} = char(strcat('u', num2str(unitname)));
                psthMatStr{3,j*2} = char(alignVar2str);
                
                if j == taggedNeurons(1) || j == taggedNeurons(2) || j == taggedNeurons(3)
                    psthMatStr{4,j*2} = char('tag');
                else
                    psthMatStr{4,j*2} = char('untagged');
                end
                %                 psthMatStr{4,j*2} = char('untagged');
                %                 psthMatStr{5,j*2} = char('tag');
                psthMatStr{5,j*2} = char('rand');
            end
        end
        %     clear psthAvg1; clear psthAvg2; clear unitname; clear alignVar1str; clear alignVar2str
        
    end
end

%Last line is the key part. The AlignRasters code is looking for:
%[response] = TNC_AlignRasters(continuous spike density function,timestamps for spikes,omit how many trials at beginning?,time stamps to align to,how much time before and adter, want raster plot structured data?, want a boxcar histogram of response?)


% plot the output something like this:

%         plot([-psthWin(1) psthWin(2)],[0 0],'-','Color',[0.2 0.2 0.2]);
%         shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'k');
