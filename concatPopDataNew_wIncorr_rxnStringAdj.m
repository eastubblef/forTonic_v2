%This script concatenates the popData.session.units for sessions in which
%the tagged & untagged units are in separate files
%created 3.6.18 to concatenate all new-task-recording files for use to do PCA or calculate diffs in ipsi/contra psths
%% Updated 3.13.18- to include concatenation of the strings for generation of the labels file (to overcome index misinterpretations)
%  Updated for only the units that have movement-related activity(pop2.mat), created by excludePopNeurons.m

% INPUT THE SESSION TYPE (blocks v rand); INPUT TAGGED NEURON UNIT NUMBERS PER SESSION
% concatenate the tagged/untagged data for 170111, 12, 18, 171012, 171013 sessions:
% alignVar 1 = ipsi; alignVar 2 = contra
% filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15'; %save the ouput matrix here
% targetName = '2015_psthMatCorrects4PCA_2';  %the first vrsion was prior to SFN 2017
% 
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/behave/tag';
% fname = strcat(fpath, '/', '151119_tag_pop.mat');                          
% behavFile = strcat(fpath, '/', '151119_tag_bh.mat');

% filenamestr = '/Volumes/My Passport for Mac/concatMvPCANew/lessUnitsNew'; %save the ouput matrix here

filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/concatMvPCANew/lessUnitsNew';
targetName = 'psthMatAdjust_mat_rxn5';  % this is the saved concatenated matrix

% fpath = '/Volumes/My Passport for Mac/171013/behaveNewUnits';
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/171013';
% fname1 = strcat(fpath, '/', '171013_pop.mat');     %should have 17 units; ipsi/contra = 36                     
% behavFile1 = strcat(fpath, '/', '171013_bh3Rxn.mat');
% load(fname1);load(behavFile1);
fname = strcat(fpath, '/', '171013_pop.mat');     %should have 17 units; ipsi/contra = 36                     
behavFile = strcat(fpath, '/', '171013_bh5Rxn.mat');
load(fname);
load(behavFile);

Pop = PopData;

% Load the variables for concatenation
wheelsAll_LfirstValid = export_mvmt_data_for_psths.wheelALL_LfirstValid;  %Added 3.6.18 for correct & incorrects
wheelsAll_RfirstValid = export_mvmt_data_for_psths.wheelALL_RfirstValid;

rxnTimesL100mv = export_mvmt_data_for_psths.rxnTimesL100mv;
rxnTimesR100mv = export_mvmt_data_for_psths.rxnTimesR100mv;
rxnTimesL4to600mv = export_mvmt_data_for_psths.rxnTimesL4to600mv;
rxnTimesR4to600mv = export_mvmt_data_for_psths.rxnTimesR4to600mv;

alignVar1 = rxnTimesL100mv;    %171013: L = ipsi
alignVar2 = rxnTimesR100mv;    %171013: R = contra
alignVar3 = rxnTimesL4to600mv;
alignVar4 = rxnTimesR4to600mv;

alignVar1str = num2str('ipsi fast');
alignVar2str = num2str('contra fast');
alignVar3str = num2str('ipsi slow');
alignVar4str = num2str('contra slow');

units = PopData.session.unit;         
behavFileStr = char(behavFile(end-16:end));
taggedNeurons = [2 16]; session = 'rand';

[psthMat_fast, psthMat_fastStr] = psthBSsubplots4concatIpsiL3(PopData, behavFileStr, units, alignVar1, alignVar2, alignVar1str, alignVar2str, taggedNeurons, session);
[psthMat_slow, psthMat_slowStr] = psthBSsubplots4concatIpsiL3(PopData, behavFileStr, units, alignVar3, alignVar4, alignVar3str, alignVar4str, taggedNeurons, session);

%should concatenate thru all of this session's units & construct these cols: 
%u1 ipsi fast, u1 contra fast, u2 ipsi fast, u2 contra fast, u3... then again for slows
psthMat171013 = horzcat(psthMat_fast, psthMat_slow);  
psthMat171013Str = horzcat({psthMat_fastStr,psthMat_slowStr});
clear alignVar1; clear alignVar2; clear psthMat_fast; clear psthMat_slow;
clear alignVar3; clear alignVar4; clear psthMat_fastStr; clear psthMat_slowStr;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%171012
% fpath = '/Volumes/My Passport for Mac/171012/newBehaveUnits';
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/171012/newBehaveUnits';

% fname = strcat(fpath, '/', '171012_pop.mat');     %should have 36 units total                          
fname = strcat(fpath, '/', '171012_pop2.mat');      %should have all units except for those w/out mv-related activity

% behavFile = strcat(fpath, '/', '171012_bh3Rxn.mat');
behavFile = strcat(fpath, '/', '171012_bh5Rxn.mat');

load(fname);
% Pop = PopData;
PopData = newPopData;
load(behavFile)
wheelsAll_LfirstValid = export_mvmt_data_for_psths.wheelALL_LfirstValid;  %Added 3.6.18
wheelsAll_RfirstValid = export_mvmt_data_for_psths.wheelALL_RfirstValid;

rxnTimesL100mv = export_mvmt_data_for_psths.rxnTimesL100mv;
rxnTimesR100mv = export_mvmt_data_for_psths.rxnTimesR100mv;
rxnTimesL4to600mv = export_mvmt_data_for_psths.rxnTimesL4to600mv;
rxnTimesR4to600mv = export_mvmt_data_for_psths.rxnTimesR4to600mv;

alignVar1 = rxnTimesL100mv;    %171012: L = ipsi
alignVar2 = rxnTimesR100mv;    %171012: R = contra
alignVar3 = rxnTimesL4to600mv;
alignVar4 = rxnTimesR4to600mv;
alignVar1str = num2str('ipsi fast');
alignVar2str = num2str('contra fast');
alignVar3str = num2str('ipsi slow');
alignVar4str = num2str('contra slow');

units = PopData.session.unit;
taggedNeurons = 0; session = 'rand';

%should have 72 cols = fast; 72 cols = slow
behavFileStr = char(behavFile(end-16:end));
[psthMat_fast, psthMat_fastStr] = psthBSsubplots4concatIpsiL3(PopData, behavFileStr, units, alignVar1, alignVar2, alignVar1str, alignVar2str, taggedNeurons, session);
[psthMat_slow, psthMat_slowStr] = psthBSsubplots4concatIpsiL3(PopData, behavFileStr, units, alignVar3, alignVar4, alignVar3str, alignVar4str, taggedNeurons, session);

psthMat171012 = horzcat(psthMat_fast, psthMat_slow);
psthMat171012Str = horzcat({psthMat_fastStr,psthMat_slowStr});

clear alignVar1; clear alignVar2; clear psthMat_fast; clear psthMat_slow;
clear alignVar3; clear alignVar4; clear psthMat_fastStr; clear psthMat_slowStr;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% 170118
% fpath = '/Volumes/My Passport for Mac/Vgatfive/170118/behaveNewUnits';
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170118/behaveNewUnits';

% fname = strcat(fpath, '/', '170118_newBehave_pop.mat');          %17 units                
fname = strcat(fpath, '/', '170118_pop2.mat');          %17 units                

behavFile = strcat(fpath, '/', '170118_bh5Rxn.mat');

load(fname);
load(behavFile)
PopData = newPopData;


% wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;   %for correct trials only
% wheelsRfirstValid = export_mvmt_data_for_psths.wheelRfirstValid;
wheelsAll_LfirstValid = export_mvmt_data_for_psths.wheelALL_LfirstValid;  %Added 3.6.18
wheelsAll_RfirstValid = export_mvmt_data_for_psths.wheelALL_RfirstValid;
rxnTimesL100mv = export_mvmt_data_for_psths.rxnTimesL100mv;
rxnTimesR100mv = export_mvmt_data_for_psths.rxnTimesR100mv;
rxnTimesL4to600mv = export_mvmt_data_for_psths.rxnTimesL4to600mv;
rxnTimesR4to600mv = export_mvmt_data_for_psths.rxnTimesR4to600mv;

alignVar1 = rxnTimesR100mv;    %170118: R = ipsi
alignVar2 = rxnTimesL100mv;    %170118: L = contra
alignVar3 = rxnTimesR4to600mv;
alignVar4 = rxnTimesL4to600mv;
alignVar1str = num2str('ipsi fast');
alignVar2str = num2str('contra fast');
alignVar3str = num2str('ipsi slow');
alignVar4str = num2str('contra slow');

units = PopData.session.unit;
taggedNeurons = [2 5 9]; session = 'rand';

%should have 34 cols = fast; 34 cols = slow
behavFileStr = char(behavFile(end-16:end));
[psthMat_fast, psthMat_fastStr] = psthBSsubplots4concatIpsiR3(PopData, behavFileStr, units, alignVar1, alignVar2, alignVar1str, alignVar2str, taggedNeurons, session);
[psthMat_slow, psthMat_slowStr] = psthBSsubplots4concatIpsiR3(PopData, behavFileStr, units, alignVar3, alignVar4, alignVar3str, alignVar4str, taggedNeurons, session);

psthMat170118 = horzcat(psthMat_fast, psthMat_slow);
psthMat170118Str = horzcat({psthMat_fastStr,psthMat_slowStr});

clear alignVar1; clear alignVar2; clear psthMat_fast; clear psthMat_slow;
clear alignVar3; clear alignVar4; clear psthMat_fastStr; clear psthMat_slowStr;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% fpath = '/Volumes/My Passport for Mac/170111/newBehaveUnits';
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170111';

fname = strcat(fpath, '/', '170111newBehave_pop.mat');         %28 units                 
behavFile = strcat(fpath, '/', '170111_bh5Rxn.mat');

load(fname);
load(behavFile)
% wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;        %For correct trials only
% wheelsRfirstValid = export_mvmt_data_for_psths.wheelRfirstValid;
wheelsAll_LfirstValid = export_mvmt_data_for_psths.wheelALL_LfirstValid;  %Added 3.6.18
wheelsAll_RfirstValid = export_mvmt_data_for_psths.wheelALL_RfirstValid;

rxnTimesL100mv = export_mvmt_data_for_psths.rxnTimesL100mv;
rxnTimesR100mv = export_mvmt_data_for_psths.rxnTimesR100mv;
rxnTimesL4to600mv = export_mvmt_data_for_psths.rxnTimesL4to600mv;
rxnTimesR4to600mv = export_mvmt_data_for_psths.rxnTimesR4to600mv;

alignVar1 = rxnTimesL100mv;    %170111: L = ipsi
alignVar2 = rxnTimesR100mv;    %170111: R = contra
alignVar3 = rxnTimesL4to600mv;
alignVar4 = rxnTimesR4to600mv;
alignVar1str = num2str('ipsi fast');
alignVar2str = num2str('contra fast');
alignVar3str = num2str('ipsi slow');
alignVar4str = num2str('contra slow');

units = PopData.session.unit;         %for Tonic data
taggedNeurons = [6 9 28]; session = 'rand';


%should have 56 cols = fast; 56 cols = slow
behavFileStr = char(behavFile(end-16:end));
[psthMat_fast, psthMat_fastStr] = psthBSsubplots4concatIpsiL3(PopData, behavFileStr, units, alignVar1, alignVar2, alignVar1str, alignVar2str, taggedNeurons, session);
[psthMat_slow, psthMat_slowStr] = psthBSsubplots4concatIpsiL3(PopData, behavFileStr, units, alignVar3, alignVar4, alignVar3str, alignVar4str, taggedNeurons, session);

psthMat170111 = horzcat(psthMat_fast, psthMat_slow);
psthMat170111Str = horzcat({psthMat_fastStr,psthMat_slowStr});

clear alignVar1; clear alignVar2; clear psthMat_fast; clear psthMat_slow;
clear alignVar3; clear alignVar4; clear psthMat_fastStr; clear psthMat_slowStr;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %170112
% fpath = '/Volumes/My Passport for Mac/170112/behaveSegs/behaveNewUnits';
% fname = strcat(fpath, '/', '170112newBehave_pop.mat');             %should have 19 units             
% behavFile = strcat(fpath, '/', '170112_bh.mat');
% 
% load(fname);
% % Pop = PopData;
% load(behavFile)
% 
% wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;
% wheelsRfirstValid = export_mvmt_data_for_psths.wheelRfirstValid;
% 
% alignVar1 = wheelsLfirstValid;    %170112: L = ipsi 
% alignVar2 = wheelsRfirstValid;    %170112: R = contra 
% units = PopData.session.unit;         %for Tonic data
% 
% [psthMat] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar1, alignVar2);
% psthMat170112 = psthMat;
% clear alignVar1; clear alignVar2; clear psthMat;

%% concatenate in this order: 151119tag,151119untag, 151118, 151105, 151104, 151106

psthMatAdjust_rxn5 = horzcat(psthMat171013, psthMat171012, psthMat170118, psthMat170111);
psthMatAdjust_rxn5Str = horzcat({psthMat171013Str}, {psthMat171012Str}, {psthMat170118Str}, {psthMat170111Str});

% should have 117 total neurons; 98 without 170112 
% psthMatCorrects4PCAnew_rxns = horzcat(psthMat171013, psthMat171012, psthMat170118, psthMat170111, psthMat170112);
%     filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15;'
     cd(filenamestr);
     
%      save([targetName],'psthMatCorr_Incorr4mat_rxn_3');
%      save(['strIndices'],'psthMatCorr_Incorr4mat_rxn_3Str');                            

    save([targetName], 'psthMatAdjust_rxn5');
    save(['strIndices5'], 'psthMatAdjust_rxn5Str');
     
    disp(['saved as ' targetName]);
                                                                     
    disp('%-------------------------------------------------------------------');
    disp(' ');


