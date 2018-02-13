
function [stimTimeA, stimTimeB] = callsROCpref(pop, fname, varargin, varAlign1, varAlign2)

%% This function extracts the following for inputs into ROC_preference  
% a - number_of_trials x 1 vector of firing rates under Condition A (varAlign1)
% b - number_of_trials x 1 vector of firing rates under Condition B (varAlign2)

% calls Josh's code for psth for firing rate to input into the ROC calculation
% find the vector for respCSS1(j).image.alignedtime = the zero (aligned) epoch and extract the respCSS1(j).image.aligned vector (num trials x 1 vector) of firing rates 

%% Inputs

fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks';
fname = strcat(fpath, '/', '151118_all_pop.mat');
behavFile = strcat(fpath, '/', '151118_all_bh.mat');

neuronType = 'tag&Untag';

load(fname);
pop = PopData;
units = 1:length(pop.session.unit);
varargin = units;

load(behavFile)
behavior = ContData.behavior;

%% define the variables for varAlign1 & 2

LrewTrialStarts = behavior.LrewTrialStarts;
RrewTrialStarts = behavior.RrewTrialStarts;
LincorrectTrialStarts = behavior.LincorrectTrialStarts;
RincorrectTrialStarts = behavior.RincorrectTrialStarts;
 
wheelsLfirstValid = behavior.wheelsLfirstValid;
wheelsRfirstValid = behavior.wheelsRfirstValid;
wheelsILfirstValid = behavior.wheelsILfirstValid;
wheelsIRfirstValid = behavior.wheelsIRfirstValid;
rewTimesValid = behavior.rewTimesValid;

varAlign1 = LrewTrialStarts;
varAlign2 = RrewTrialStarts;

PopData.currParams.stableTrials = -1;                                         % subtracts the first "trial" when calling AlignRaster

%% 
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 15.4;
[currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);

for i = 1:length(PopData.session)
%     figure; hold on; 
%     plotname = 'corr L/R trialStart aligned'
    dimm=ceil(sqrt(length(units))); %

    for j = 1:length(units)
        numStamps = length(PopData.session(i).unit(j).ts);
        delta = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
        delta(1,ceil(PopData.session(i).unit(j).ts)) = 1;
        tmpSmooth = conv(delta,currParams.filter.kernel,'same');

        psthWin = [1.0e3,2e3]; 
                
        % condition a
        [respCSS1(j)] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts, PopData.currParams.stableTrials, varAlign1, psthWin, 1, 1); 
%         unitname = num2str(j);


        % condition b
        [respCSS2(j)] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts, PopData.currParams.stableTrials, varAlign2, psthWin, 1, 1);

    end
end

%% Call the ROC_preference m file:
% INPUTS:  
% a - number_of_trials x 1 vector of firing rates under Condition A
% b - number_of_trials x 1 vector of firing rates under Condition B
% num_repeats - number of times to permute a and b and recalculate preference (used for calculating the p_val associates with pref)

for j = 1:length(units)

%     a = respCSS1(j).raster.trial.ts;
%     b = respCSS1(j).raster.trial.ts;
    a = respCSS1(j).image.aligned 
    num_repeats = 500;

    [pref(j), p_val(j)] = ROC_preference(a, b, num_repeats)

end
