%This script concatenates the popData.session.units for sessions in which
%the tagged & untagged units are in separate files
%created 3.6.18 to concatenate all new-task-recording files for use to do PCA or calculate diffs in ipsi/contra psths

%concatenate the tagged/untagged data for 170111, 12, 18, 171012, 171013 sessions:
% alignVar 1 = ipsi; alignVar 2 = contra
% filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15'; %save the ouput matrix here
% targetName = '2015_psthMatCorrects4PCA_2';  %the first vrsion was prior to SFN 2017
% 
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/behave/tag';
% fname = strcat(fpath, '/', '151119_tag_pop.mat');                          
% behavFile = strcat(fpath, '/', '151119_tag_bh.mat');

filenamestr = '/Volumes/My Passport for Mac/concatMvPCANew/lessUnitsNew'; %save the ouput matrix here
% filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15'; %save the ouput matrix here for later concat w/ 2015 data
targetName = '2018_psthMatCorr_Incorr_matrix_rxn';  

fpath = '/Volumes/My Passport for Mac/171013/behaveNewUnits';
fname = strcat(fpath, '/', '171013_pop.mat');     %should have 17 units; ipsi/contra = 36                     
behavFile = strcat(fpath, '/', '171013_bh3Rxn.mat');

load(fname);
Pop = PopData;
load(behavFile)

% wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;
% wheelsRfirstValid = export_mvmt_data_for_psths.wheelRfirstValid;
wheelsAll_LfirstValid = export_mvmt_data_for_psths.wheelALL_LfirstValid;  %Added 3.6.18
wheelsAll_RfirstValid = export_mvmt_data_for_psths.wheelALL_RfirstValid;

rxnTimesL100mv = export_mvmt_data_for_psths.rxnTimesL100mv;
rxnTimesR100mv = export_mvmt_data_for_psths.rxnTimesR100mv;
rxnTimesL4to600mv = export_mvmt_data_for_psths.rxnTimesL4to600mv;
rxnTimesR4to600mv = export_mvmt_data_for_psths.rxnTimesR4to600mv;

alignVar1 = rxnTimesL100mv;    %171013: L = ipsi
alignVar2 = rxnTimesR100mv;    %171013: R = contra
alignVar3 = rxnTimesL4to600mv;
alignVar4 = rxnTimesR4to600mv;

units = PopData.session.unit;         

[psthMat_fast] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar1, alignVar2);
[psthMat_slow] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar3, alignVar4);
%should concatenate thru all of this session's units & construct these cols: 
%u1 ipsi fast, u1 contra fast, u2 ipsi fast, u2 contra fast, u3... then again for slows
psthMat171013 = horzcat(psthMat_fast, psthMat_slow);  

clear alignVar1; clear alignVar2; clear psthMat_fast; clear psthMat_slow;
clear alignVar3; clear alignVar4; 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%171012
fpath = '/Volumes/My Passport for Mac/171012/newBehaveUnits';
fname = strcat(fpath, '/', '171012_pop.mat');  %should have 36 units total                          
% behavFile = strcat(fpath, '/', '171012_bh3Rxn.mat');

behavFile = strcat(fpath, '/', '171012_bh4Rxn.mat');

load(fname);
Pop = PopData;
load(behavFile)
% wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;
% wheelsRfirstValid = export_mvmt_data_for_psths.wheelRfirstValid;
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

units = PopData.session.unit;        

%should have 72 cols = fast; 72 cols = slow
[psthMat_fast] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar1, alignVar2);
[psthMat_slow] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar3, alignVar4);
psthMat171012 = horzcat(psthMat_fast, psthMat_slow);

clear alignVar1; clear alignVar2; clear psthMat_fast; clear psthMat_slow;
clear alignVar3; clear alignVar4; 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% 170118
fpath = '/Volumes/My Passport for Mac/Vgatfive/170118/behaveNewUnits';
fname = strcat(fpath, '/', '170118_newBehave_pop.mat');          %17 units                
behavFile = strcat(fpath, '/', '170118_bh3Rxn.mat');

load(fname);
load(behavFile)

% wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;
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

units = PopData.session.unit;         
%should have 34 cols = fast; 34 cols = slow
[psthMat_fast] = psthBSsubplots4concatIpsiR2(PopData, behavFile, units, alignVar1, alignVar2);
[psthMat_slow] = psthBSsubplots4concatIpsiR2(PopData, behavFile, units, alignVar3, alignVar4);
psthMat170118 = horzcat(psthMat_fast, psthMat_slow);
clear alignVar1; clear alignVar2; clear psthMat_fast; clear pstMat_slow;
clear alignVar3; clear alignVar4;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fpath = '/Volumes/My Passport for Mac/170111/newBehaveUnits';
fname = strcat(fpath, '/', '170111newBehave_pop.mat');         %28 units                 
behavFile = strcat(fpath, '/', '170111_bh3Rxn.mat');

load(fname);
load(behavFile)
% wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;
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

units = PopData.session.unit;         %for Tonic data

%should have 56 cols = fast; 56 cols = slow
[psthMat_fast] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar1, alignVar2);
[psthMat_slow] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar3, alignVar4);
psthMat170111 = horzcat(psthMat_fast, psthMat_slow);
clear alignVar1; clear alignVar2; clear psthMat_fast; clear psthMat_slow;
clear alignVar3; clear alignVar4;

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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %151106
% % fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151106/behave';
% % fname = strcat(fpath, '/', '151106_all_pop.mat');                          
% % behavFile = strcat(fpath, '/', '151106_all_bh.mat');
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151106/behave/newUnits';
% fname = strcat(fpath, '/', '151106_newUnits2Behave_pop.mat');         %should have 20 u                 
% behavFile = strcat(fpath, '/', '151106_newUnits2Behave_bh.mat');
% 
% load(fname);
% % Pop = PopData;
% 
% load(behavFile)
% behavior = ContData.behavior;
% 
% wheelsLfirstValid = behavior.wheelsLfirstValid;
% wheelsRfirstValid = behavior.wheelsRfirstValid;
% 
% alignVar1 = wheelsLfirstValid;    %151105: L = ipsi 
% alignVar2 = wheelsRfirstValid;    %151105: R = contra 
% units = PopData.session.unit;         %for Tonic data
% 
% [psthMat] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar1, alignVar2);
% psthMat06 = psthMat;

%% concatenate in this order: 151119tag,151119untag, 151118, 151105, 151104, 151106

psthMatCorr_Incorr4mat_rxn_2 = horzcat(psthMat171013, psthMat171012, psthMat170118, psthMat170111);
% should have 117 total neurons; 98 without 170112 
% psthMatCorrects4PCAnew_rxns = horzcat(psthMat171013, psthMat171012, psthMat170118, psthMat170111, psthMat170112);

%     filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15;'
     cd(filenamestr);
     
     save([targetName],'psthMatCorr_Incorr4mat_rxn_2');                            
     disp(['saved as ' targetName]);
                                                                     
    disp('%-------------------------------------------------------------------');
    disp(' ');


