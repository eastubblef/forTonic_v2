%This script concatenates the popData.session.units for sessions in which
%the tagged & untagged units are in separate files

%concatenate the tagged/untagged data for 170111, 12, 18, 171012, 171013 sessions:
% alignVar 1 = ipsi; alignVar 2 = contra
% filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15'; %save the ouput matrix here
% targetName = '2015_psthMatCorrects4PCA_2';  %the first vrsion was prior to SFN 2017
% 
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/behave/tag';
% fname = strcat(fpath, '/', '151119_tag_pop.mat');                          
% behavFile = strcat(fpath, '/', '151119_tag_bh.mat');

filenamestr = '/Volumes/My Passport for Mac/concatMvPCANew'; %save the ouput matrix here
targetName = '2017_psthMatCorrects4PCA';  

fpath = '/Volumes/My Passport for Mac/171013/behaveNewUnits';
fname = strcat(fpath, '/', '171013_pop.mat');     %should have 17 units; ipsi/contra = 36                     
behavFile = strcat(fpath, '/', '171013_bh.mat');

load(fname);
Pop = PopData;

load(behavFile)

wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;
wheelsRfirstValid = export_mvmt_data_for_psths.wheelRfirstValid;

alignVar1 = wheelsLfirstValid;    %171013: L = ipsi
alignVar2 = wheelsRfirstValid;    %171013: R = contra
units = PopData.session.unit;         

[psthMat] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar1, alignVar2);
psthMat171013 = psthMat;

clear alignVar1; clear alignVar2; clear psthMat;


%171012
fpath = '/Volumes/My Passport for Mac/171012/newBehaveUnits';
fname = strcat(fpath, '/', '171012_pop.mat');  %should have 36 units total                          
behavFile = strcat(fpath, '/', '171012_bh.mat');

load(fname);
Pop = PopData;

load(behavFile)
wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;
wheelsRfirstValid = export_mvmt_data_for_psths.wheelRfirstValid;

alignVar1 = wheelsLfirstValid;    %171012: L = ipsi 
alignVar2 = wheelsRfirstValid;    %171012: R = contra 
units = PopData.session.unit;        

[psthMat] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar1, alignVar2);
psthMat171012 = psthMat;

clear alignVar1; clear alignVar2; clear psthMat;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %concatenate the 151118 tagged data (units 3,4,5):
% % alignVar 1 = ipsi
% clear alignVar1; clear alignVar2; clear units;
% 
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks';
% fname = strcat(fpath, '/', '151118_all_pop.mat');                          
% behavFile = strcat(fpath, '/', '151118_all_bh.mat');
% load(fname);
% Pop = PopData;
% load(behavFile)
% behavior = ContData.behavior;
% 
% alignVar1 = behavior.wheelsRfirstValid;    %151119: R = ipsi 
% alignVar2 = behavior.wheelsLfirstValid;    %151119: L = contra 
% 
% units = 3;
% [alignVarMat] = OnlyOnePsthBSsubplots4concatIpsiR(Pop, behavFile, units, alignVar1, alignVar2);
% alignVarMat18_3 = alignVarMat;

% 170118
fpath = '/Volumes/My Passport for Mac/Vgatfive/170118/behaveNewUnits';
fname = strcat(fpath, '/', '170118_newBehave_pop.mat');          %17 units                
behavFile = strcat(fpath, '/', '170118_bh.mat');

load(fname);
Pop = PopData;

load(behavFile)

wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;
wheelsRfirstValid = export_mvmt_data_for_psths.wheelRfirstValid;

alignVar1 = wheelsRfirstValid;    %170118: R = ipsi 
alignVar2 = wheelsLfirstValid;    %170118: L = contra 
units = PopData.session.unit;         

[psthMat] = psthBSsubplots4concatIpsiR2(PopData, behavFile, units, alignVar1, alignVar2);
psthMat170118 = psthMat;
clear alignVar1; clear alignVar2; clear psthMat;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fpath = '/Volumes/My Passport for Mac/170111/newBehaveUnits';
fname = strcat(fpath, '/', '170111newBehave_pop.mat');         %28 units                 
behavFile = strcat(fpath, '/', '170111_bh.mat');

load(fname);
% Pop = PopData;

load(behavFile)
wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;
wheelsRfirstValid = export_mvmt_data_for_psths.wheelRfirstValid;

alignVar1 = wheelsLfirstValid;    %170111: L = ipsi 
alignVar2 = wheelsRfirstValid;    %170111: R = contra 
units = PopData.session.unit;         %for Tonic data

[psthMat] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar1, alignVar2);
psthMat170111 = psthMat;
clear alignVar1; clear alignVar2; clear psthMat;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%170112
fpath = '/Volumes/My Passport for Mac/170112/behaveSegs/behaveNewUnits';
fname = strcat(fpath, '/', '170112newBehave_pop.mat');             %should have 19 units             
behavFile = strcat(fpath, '/', '170112_bh.mat');

load(fname);
% Pop = PopData;
load(behavFile)

wheelsLfirstValid = export_mvmt_data_for_psths.wheelLfirstValid;
wheelsRfirstValid = export_mvmt_data_for_psths.wheelRfirstValid;

alignVar1 = wheelsLfirstValid;    %170112: L = ipsi 
alignVar2 = wheelsRfirstValid;    %170112: R = contra 
units = PopData.session.unit;         %for Tonic data

[psthMat] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar1, alignVar2);
psthMat170112 = psthMat;
clear alignVar1; clear alignVar2; clear psthMat;

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

psthMatCorrects4PCAnew = horzcat(psthMat171013, psthMat171012, psthMat170118, psthMat170111, psthMat170112);
%  should have 117 total neurons  
% psthMatCorrects4PCAnew_rxns = horzcat(psthMat171013, psthMat171012, psthMat170118, psthMat170111, psthMat170112);

%     filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15;'
     cd(filenamestr);
     
     save([targetName],'psthMatCorrects4PCAnew');                            
     disp(['saved as ' targetName]);
                                                                     
    disp('%-------------------------------------------------------------------');
    disp(' ');


