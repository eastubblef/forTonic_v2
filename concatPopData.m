%This script concatenates the popData.session.units for sessions in which
%the tagged & untagged units are in separate files

%concatenate the tagged data for 151119
% alignVar 1 = ipsi
filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15'; %save the ouput matrix here
targetName = '2015_psthMatCorrects4PCA_2';  %the first vrsion was prior to SFN 2017

fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/behave/tag';
fname = strcat(fpath, '/', '151119_tag_pop.mat');                          
behavFile = strcat(fpath, '/', '151119_tag_bh.mat');

load(fname);
Pop = PopData;

load(behavFile)
behavior = ContData.behavior;

wheelsLfirstValid = behavior.wheelsLfirstValid;
wheelsRfirstValid = behavior.wheelsRfirstValid;

alignVar1 = wheelsRfirstValid;    %151119: R = ipsi 
alignVar2 = wheelsLfirstValid;    %151119: L = contra 
units = PopData.session.unit;         %for Tonic data

[psthMat] = psthBSsubplots4concatIpsiR2(PopData, behavFile, units, alignVar1, alignVar2);
psthMat19tag = psthMat;

clear alignVar1; clear alignVar2; clear psthMat;

%concatenate the UNtagged data for 151119
% alignVar 1 = ipsi
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/behave/untag';
% fname = strcat(fpath, '/', '151119_untag_pop.mat');                          
% behavFile = strcat(fpath, '/', '151119_untag_bh.mat');

%new untagged units
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/behave/untag/newBehaveUnits/untagged';
fname = strcat(fpath, '/', '151119_newUnits2Behave_pop.mat');  %should have 13 units total                          
behavFile = strcat(fpath, '/', '151119_newUnits2Behave_bh.mat');

load(fname);
Pop = PopData;

load(behavFile)
behavior = ContData.behavior;

wheelsLfirstValid = behavior.wheelsLfirstValid;
wheelsRfirstValid = behavior.wheelsRfirstValid;

alignVar1 = wheelsRfirstValid;    %151119: R = ipsi 
alignVar2 = wheelsLfirstValid;    %151119: L = contra 
units = PopData.session.unit;         %for Tonic data

[psthMat] = psthBSsubplots4concatIpsiR2(PopData, behavFile, units, alignVar1, alignVar2);
psthMat19untag = psthMat;

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

%concatenate the tagged data for 151118
% alignVar 1 = ipsi
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks';
% fname = strcat(fpath, '/', '151118_all_pop.mat');                          
% behavFile = strcat(fpath, '/', '151118_all_bh.mat');

fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks/newBehaveUnits/behaveChunks';
fname = strcat(fpath, '/', '151118_newUnits2Behave_pop.mat');          %25 units                
behavFile = strcat(fpath, '/', '151118_newUnits2Behave_bh.mat');

load(fname);
Pop = PopData;

load(behavFile)
behavior = ContData.behavior;

wheelsLfirstValid = behavior.wheelsLfirstValid;
wheelsRfirstValid = behavior.wheelsRfirstValid;

alignVar1 = wheelsRfirstValid;    %151118: R = ipsi 
alignVar2 = wheelsLfirstValid;    %151118: L = contra 
units = PopData.session.unit;         %for Tonic data

[psthMat] = psthBSsubplots4concatIpsiR2(PopData, behavFile, units, alignVar1, alignVar2);
psthMat18 = psthMat;
clear alignVar1; clear alignVar2; clear psthMat;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%151105
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105';
% fname = strcat(fpath, '/', '151105_all_pop.mat');                          
% behavFile = strcat(fpath, '/', '151105_all_bh.mat');
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/newBehaveUnits/behaveChunks';
fname = strcat(fpath, '/', '151105_newUnits2Behave_pop.mat');                          
behavFile = strcat(fpath, '/', '151105_newUnits2Behave_bh.mat');

load(fname);
% Pop = PopData;

load(behavFile)
behavior = ContData.behavior;

wheelsLfirstValid = behavior.wheelsLfirstValid;
wheelsRfirstValid = behavior.wheelsRfirstValid;

alignVar1 = wheelsLfirstValid;    %151105: L = ipsi 
alignVar2 = wheelsRfirstValid;    %151105: R = contra 
units = PopData.session.unit;         %for Tonic data

[psthMat] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar1, alignVar2);
psthMat05 = psthMat;
clear alignVar1; clear alignVar2; clear psthMat;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%151104
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151104';
fname = strcat(fpath, '/', '151104_all_pop.mat');             %should have 22 units             
behavFile = strcat(fpath, '/', '151104_all_bh.mat');

load(fname);
% Pop = PopData;

load(behavFile)
behavior = ContData.behavior;

wheelsLfirstValid = behavior.wheelsLfirstValid;
wheelsRfirstValid = behavior.wheelsRfirstValid;

alignVar1 = wheelsLfirstValid;    %151105: L = ipsi 
alignVar2 = wheelsRfirstValid;    %151105: R = contra 
units = PopData.session.unit;         %for Tonic data

[psthMat] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar1, alignVar2);
psthMat04 = psthMat;
clear alignVar1; clear alignVar2; clear psthMat;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%151106
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151106/behave';
% fname = strcat(fpath, '/', '151106_all_pop.mat');                          
% behavFile = strcat(fpath, '/', '151106_all_bh.mat');
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151106/behave/newUnits';
fname = strcat(fpath, '/', '151106_newUnits2Behave_pop.mat');         %should have 20 u                 
behavFile = strcat(fpath, '/', '151106_newUnits2Behave_bh.mat');

load(fname);
% Pop = PopData;

load(behavFile)
behavior = ContData.behavior;

wheelsLfirstValid = behavior.wheelsLfirstValid;
wheelsRfirstValid = behavior.wheelsRfirstValid;

alignVar1 = wheelsLfirstValid;    %151105: L = ipsi 
alignVar2 = wheelsRfirstValid;    %151105: R = contra 
units = PopData.session.unit;         %for Tonic data

[psthMat] = psthBSsubplots4concatIpsiL2(PopData, behavFile, units, alignVar1, alignVar2);
psthMat06 = psthMat;

%% concatenate in this order: 151119tag,151119untag, 151118, 151105, 151104, 151106

psthMatCorrects4PCA_2 = horzcat(psthMat19tag, psthMat19untag, psthMat18, psthMat05, psthMat04, psthMat06);
    
%     filenamestr = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15;'
     cd(filenamestr);
     
     save([targetName],'psthMatCorrects4PCA_2');                            
     disp(['saved as ' targetName]);
                                                                     
    disp('%-------------------------------------------------------------------');
    disp(' ');


