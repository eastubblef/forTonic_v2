%% Classify PSTHs based upon recording properties
% Updated 2.27.18 for new 2017 datasets
%  Updated 1.31.17 to show PC1 ipsi v. PC1 contra (rather than PC1 v PC3 for both
%  specifically load psths aligned to mvment onset
%  Updated 1.24.17, for 151119, 18, 04, 05, 06 recordings (new untagged neurons added)

clear;
cd '/Volumes/My Passport for Mac/concatMvPCANew/lessUnitsNew';

%New for 49 neurons - Load the older data:

%New for 147 neurons, 1.24.18:
S = load('2018_psthMatCorr_Incorr_matrix_rxn3');
% [num,txt,raw1] = xlsread('2017_psthMatIndsLabel_4sess_rxn.xlsx');
% [num,txt,raw1] = xlsread('2018_psthIndsLabel_rxnNew2.xlsx');
[num,txt,raw1] = xlsread('stringIndex3.xlsx');

%Concatenate the matrices
txtMat = [];
for j = 1:length(txt)
    for i = 1:5                                   %depends on the number of variables
        if ischar(txt{i,j})
           str(i,j) = string(txt{i,j}); 
        else break
        end
    end
end
txtMat = str'; %just put the 3rd col. 1st, copy & paste into the tpose wrksheet.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear num; clear raw1; clear txt;
[num,txt,raw] = xlsread('2018_psthIndsLabel_rxnAll.xlsx', 'tpose');

