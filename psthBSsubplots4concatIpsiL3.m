
function [psthMat, psthMatStr] = psthBSsubplots4concatIpsiL3(PopData, behavFileStr, units, alignVar1, alignVar2, alignVar1str, alignVar2str, taggedNeurons, session)
%% This m file calls Josh's code in order to plot psths (as subplots like oneUnitAlignSpike code for rasters plots)
% calls the following mfiles:
% TNC_createGaussian.m - explanatory
% TNC_alignRasters.m -  for creating inputted time window (line 81), setting the inputted zero-alignment of pop data (units) to behavioral data (behavFile)
% corresponding to ipsi v contra "alignVars 1 & 2" (lines 51 & 52) 
% shadedErrorBar.m -  for plotting 

%% Updated 3.20.18 to generate strings & psths 

% Run this in conjunction with concatPopDataNew_wIncorr_rxnString.m
%% use the L2 version to load the following, but might have string errors?

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

% load(popFile);
% pop = PopData;
%
% % allow for plotting only 1 of them
% units = pop.session.unit;  %for Tonic data
%
% % units1 = unit(5);
% units = 1:length(units1);
if numel(units) == 1
    units = PopData.session.unit(1);
else
    units = 1:length(PopData.session.unit);    % plots all of them
    varargin = units;
end
% load(behavFile)
% behavior = ContData.behavior;

% trialStarts = behavior.trialStarts;

PopData.currParams.stableTrials = -1;                                         % subtracts the first "trial" when calling AlignRaster

%%
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 15.4;
[currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);

% alignVarMat = NaN(1000, 1000);  %added this section 11.2017 for building a matrix for dim reduction
psthMatStr = {}; %fill this with strings
psthMat = [];

% for k = 1:length(taggedNeurons)

for i = 1:length(PopData.session)
    figure; hold on;
    plotname = char(behavFileStr);
    
    dimm=ceil(sqrt(length(units))); %
    
    for j = 1:length(units)
        for k = 1:length(taggedNeurons)
        
        numStamps = length(PopData.session(i).unit(j).ts);
        delta = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
        delta(1,ceil(PopData.session(i).unit(j).ts)) = 1;
        tmpSmooth = conv(delta,currParams.filter.kernel,'same');
%         psthWin = [1.0e3,2e3]; % best for 151105 data
        
        psthWin = [2.0e3,2e3]; % best for 151105 data
        subplot(dimm,dimm,j); hold on;
        
        % first plot
        [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts, PopData.currParams.stableTrials, alignVar1, psthWin, 1, 1);
        unitname = num2str(j);
        shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'bl');  %blue for left; for 151105 dataset = ipsi
        
        % second plot
        [respCSS2] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts, PopData.currParams.stableTrials, alignVar2, psthWin, 1, 1);
        h(j) = plot([-(psthWin(1)) psthWin(2)],[0 0]);
        set(h(j), 'Color',[0.2 0.2 0.2]);
        ylabel('Firing Rate'); xlabel('ms');
        shadedErrorBar(-psthWin(1):psthWin(2), respCSS2.image.psthZ, respCSS2.image.psthZe,'r');  %red for right; for 151105 data = contra
        
        title(behavFileStr(1:6));
        
        psthAvg1 = respCSS.image.psthZ';
        psthAvg2 = respCSS2.image.psthZ';
        
        for z = 1:length(psthAvg1)
            
            if j == 1
                psthMat(z,j) = psthAvg1(z);
                psthMatStr{1,j} = char(behavFileStr(1:6));
                psthMatStr{2,j} = char(strcat('u', num2str(unitname)));
                psthMatStr{3,j} = char(alignVar1str);
                psthMatStr{4,j} = char('untagged');  %As long as no session's first u is tagged
                psthMatStr{5,j} = char('rand');
                
            else
                psthAvg1 = respCSS.image.psthZ';
                psthMat(z,j+(j-1)) = psthAvg1(z);
                psthMatStr{1,j+(j-1)} = char(behavFileStr(1:6));
                psthMatStr{2,j+(j-1)} = char(strcat('u', num2str(unitname)));
                psthMatStr{3,j+(j-1)} = char(alignVar1str);
                %                     psthMatStr{4,j+(j-1)} = char('untagged');
                psthMatStr{5,j+(j-1)} = session;
                
                %                 for k = 1:length(taggedNeurons)
                if taggedNeurons(k) == j
                    psthMatStr{4,j+(j-1)} = char('tag');
                    if k < length(taggedNeurons)
                        k = k+1;
                        break
                        % %                     else psthMatStr{4,j+(j-1)} = char('untagged');
                    end
                end
            end
        end
        
        for zz = 1:length(psthAvg2)
            if j == 1
                psthMat(zz,j*2) = psthAvg2(zz);
                psthMatStr{1,j*2} = char(behavFileStr(1:6));
                psthMatStr{2,j*2} = char(strcat('u', num2str(unitname)));
                psthMatStr{3,j*2} = char(alignVar2str);
%                 psthMatStr{4,j*2} = char('untagged');
                psthMatStr{5,j*2} = session;
            else
                psthAvg2 = respCSS2.image.psthZ';
                psthMat(zz,j*2) = psthAvg2(zz);           %This should make col 1 = alignVar1; col2 = alignVar2 for ea. unit to build psthMat
                psthMatStr{1,j*2} = char(behavFileStr(1:6));
                psthMatStr{2,j*2} = char(strcat('u', num2str(unitname)));
                psthMatStr{3,j*2} = char(alignVar2str);
                %                     psthMatStr{4,j*2} = char('untagged');
                psthMatStr{5,j*2} = session;
                
                if taggedNeurons(k) == j
                    psthMatStr{4,j*2} = char('tag');
                    
                    if k < length(taggedNeurons)
                        k = k+1;
                        break
                    end
                    %                     else psthMatStr{4,j*2} = char('untagged');
                end
              end
           end
        end
    end
end
%% These lines work only for the 171013 session              
%         for z = 1:length(psthAvg1)
%             if j == 1
%                 psthMat(z,j) = psthAvg1(z);
%                 psthMatStr{1,j} = char(behavFileStr(1:6));
%                 psthMatStr{2,j} = char(strcat('u', num2str(unitname)));
%                 psthMatStr{3,j} = char(alignVar1str);
%                 psthMatStr{4,j} = char('untagged');
%                 psthMatStr{5,j} = char('rand');
%                 
%             else
%                 psthAvg1 = respCSS.image.psthZ';
%                 psthMat(z,j+(j-1)) = psthAvg1(z);
%                 psthMatStr{1,j+(j-1)} = char(behavFileStr(1:6));
%                 psthMatStr{2,j+(j-1)} = char(strcat('u', num2str(unitname)));
%                 psthMatStr{3,j+(j-1)} = char(alignVar1str);
%                 
%                 if numel(taggedNeurons) == 2 && j == taggedNeurons(1) || j == taggedNeurons(2) 
%                     psthMatStr{4,j+(j-1)} = char('tag');
%                     else
%                         psthMatStr{4,j+(j-1)} = char('untagged');
%                 end
%                 
%                 psthMatStr{5,j+(j-1)} = char('rand');
%             end
%         end
%         
%         for zz = 1:length(psthAvg2)
%             if j == 1
%                 psthMat(zz,j*2) = psthAvg2(zz);
%                 psthMatStr{1,j*2} = char(behavFileStr(1:6));
%                 psthMatStr{2,j*2} = char(strcat('u', num2str(unitname)));
%                 psthMatStr{3,j*2} = char(alignVar2str);
%                 psthMatStr{4,j*2} = char('untagged');
%                 psthMatStr{5,j*2} = char('rand');
%             else
%                 psthAvg2 = respCSS2.image.psthZ';
%                 psthMat(zz,j*2) = psthAvg2(zz);           %This should make col 1 = alignVar1; col2 = alignVar2 for ea. unit to build psthMat
%                 psthMatStr{1,j*2} = char(behavFileStr(1:6));
%                 psthMatStr{2,j*2} = char(strcat('u', num2str(unitname)));
%                 psthMatStr{3,j*2} = char(alignVar2str);
%                 
%                 if numel(taggedNeurons) == 2 && j == taggedNeurons(1) || j == taggedNeurons(2) 
%                     psthMatStr{4,j*2} = char('tag');
%                     else
%                         psthMatStr{4,j*2} = char('untagged');
%                 end
%                 psthMatStr{5,j*2} = char('rand');
%             end
%         end
%         
%     end
% end

%Last line is the key part. The AlignRasters code is looking for:
%[response] = TNC_AlignRasters(continuous spike density function,timestamps for spikes,omit how many trials at beginning?,time stamps to align to,how much time before and adter, want raster plot structured data?, want a boxcar histogram of response?)


% plot the output something like this:

%         plot([-psthWin(1) psthWin(2)],[0 0],'-','Color',[0.2 0.2 0.2]);
%         shadedErrorBar(-psthWin(1):psthWin(2), respCSS.image.psthZ, respCSS.image.psthZe,'k');
