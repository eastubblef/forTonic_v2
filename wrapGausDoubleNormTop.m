%% For Paper 2014 plotting the first 4 plots of Fig 5 = histograms; Based on m file called 'wrapGausDoubleNorm.m'
clear;
clear global;

cd '/Users/stubblefielde/Desktop/Felsen Lab/database/updated_14'
singleDoubleDir = '/Users/stubblefielde/Desktop/Felsen Lab/database/updated_14';
[num, txt] = xlsread('Beth_GausFits5.xlsx'); %For single v double better logical arrays

% Pulled values should all be from Beth_GausFits.xlsx
% Based on Beth_databaseLess.xls and Cntrl2.xls

% pvals = 'C:\Felsen Lab docs\Matlab\Data\spikes\FR\200ms\2n4Hz_1pulse\lateFile_pvals2n4' %deleted 1/3/15

%% These 2 files contain the matrices needed:
%Inputs:

cd '/Users/stubblefielde/Desktop/FRs/FRfiles/2n4Hz_1pulse/avg/';
lateFile = '/Users/stubblefielde/Desktop/FRs/FRfiles/2n4Hz_1pulse/avg/lateFile_data2n4';
% lateFile = 'C:\Felsen Lab docs\Matlab\Data\spikes\FR\200 ms\2n4Hz_1pulse\avg\lateFile_data2n4';
load(lateFile);
lateFileCntrl = '/Users/stubblefielde/Desktop/FRs/FRfiles/2n4Hz_1pulse/avg/lateFileCntrl2n4'
% lateFileCntrl = 'C:\Felsen Lab docs\Matlab\Data\spikes\FR\200 ms\2n4Hz_1pulse\avg\lateFileCntrl2n4';
load(lateFileCntrl);
x_vals = [0:0.005:0.1975];                                                  %based on using histc with bin edges of [0:0.005:0.2] (5 msec bins from 0 to 200 msec)

%% Input the proper matrix of 4 latency params
matrixChoice = matLateResAll2n4;
late_hist = matrixChoice;                                                   %This one has 25 neurons; most current 7/23/14
late_histCntrl = matLateCntrl2n4;                                           %This one has 65 neurons; most current 7/23/14

%% 7/23/14 normalize the firing rates for better MSE vals for fits: FRnorm = FR/sum(FR); 
%But don't use norms for binning. To get MSEs use
%wrapperGausFitsDoubVSingle2

%% These inputs will fit latencies' 4 params: single Gaussian!
no_plot_flag = 1                                                            %No plotting

numRes = 25
noRes = 65
numNeurons = 1:numRes;                                                      %specified by Beth_database3.xls
numCntrl = 1:noRes;                                                         %specified by Beth_database3.xls

% Don't include the normalized data for binning:
%   fit_params = 4 element vector consisting of the fitted
%   amplitude1, center1, width, and baseline (in that order)
numParams = 4;                                                              %specified by single Gaussian fits (num of parameters to calc)
longNeuronsFit = nan(length(numNeurons), numParams);
longNeuronsInfo = nan(length(numNeurons), numParams);
longFitSSEsingle = nan(length(numNeurons), 1);
 
cd '/Users/stubblefielde/Desktop/Felsen Lab/mfiles/updated';

for i = 1:length(matrixChoice)                                              %note that late_hist is for binning later (not normalized) average FR
    [fit_params, fit_info] = ConstrainedGaussianFit(x_vals, late_hist(:,i), no_plot_flag);
    longNeuronsFit(i,:) = fit_params;                                       %note that lateHistNorm is NORMALIZED average FR
    longFitSSEsingle(i) = fit_info.gof.sse; 
        if i == length(numNeurons);
           break
    
        end
      
 end
%  longFitMSE = longFitSSEsingle ./ 40; to generate vals
%  in ppt

 %% These inputs will fit latencies' 6 params: double Gaussians!
%  fit_params = 6 element vector consisting of the fitted
%  amplitude1, amplitude2, center1, center2, width, and baseline (in that order)
numParams = 6;                                                              %specified by double Gaussian fits (num of params to fit)                               
longDoubleFit = nan(length(numNeurons), numParams);
longDoubleInfo = nan(length(numNeurons), numParams);
longFitSSEdouble = nan(length(numNeurons), 1);

 for i = 1:length(matrixChoice)   
    [fit_params, fit_info] = ConstrainedDoubleGaussianFit(x_vals, late_hist(:,i), no_plot_flag); %late_hist substitute for binning!
    longDoubleFit(i,:) = fit_params;
    longFitSSEdouble(i) = fit_info.gof.sse;
        if i == length(numNeurons);
           break
    
        end
        
 end
longFitdoubleMSE = longFitSSEdouble ./ 40;  %to generate vals in ppt 
longFitdoubleMSE = longFitdoubleMSE';

%% For control neurons, called & saved from database_2; single fits
numParams = 4;
longCntrlFit = nan(length(numCntrl), numParams);
longCntrlInfo = nan(length(numCntrl), numParams);

% Don't include the normalized data to be fit here:
 for c = 1:length(matLateCntrl2n4)   
    [fit_params, fit_info] = ConstrainedGaussianFit(x_vals, late_histCntrl(:,c), no_plot_flag);
    longCntrlFit(c,:) = fit_params;   
    CntrlFitSSEdouble(i,:) = fit_info.gof.sse;
        if c == length(numCntrl);
           break
    
        end
        
 end

%% Determine which neurons were best fit with single v. double from xlsread
neuronNum = num(1:end, 1);
singleFits = num(1:end, 7);
doubleFits = num(1:end, 8);

%Calculated 14 singleInds and 10 doubleInds
%Now index singleInds and doubleInds into the proper neuron...
longIndFit = NaN(length(numNeurons), 1);
for newi = 1:length(singleFits);
    if singleFits(newi) == 1
        longIndFit(newi) = longNeuronsFit(newi);
    else
        longIndFit(newi) = longDoubleFit(newi);
    end
    
end

%% Plot individual histogram examples (A & B)
PresentationFigSetUp;
hold on;
late_hist = matLateResAll2n4(:,1);
no_plot_flag = 1;
[fit_params, fit_info, fit_object] = ConstrainedGaussianFit(x_vals, late_hist, no_plot_flag);

subplot(2,2,1);
hold on;
    p = plot(x_vals, late_hist, 'ko');
    set(p, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 5);

    p2 = plot(fit_object);
    set(p2, 'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', [.5, .5, .5]);
    legend('off');
    set(p, 'LineWidth', 2);
    set(findobj('Type','line'),'Color', 'k');
    xlabel('Time from first pulse on (s)'); ylabel('Firing rate'); 

hold on;
late_hist = matLateResAll2n4(:, 6);   

[fit_params, fit_info, fit_object] = ConstrainedDoubleGaussianFit(x_vals, late_hist, no_plot_flag);
    
    subplot(2,2,2);
    hold on;
    p3 = plot(x_vals, late_hist, 'ko');
    set(p3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 5);
    p4 = plot(fit_object);
%     set(p2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 2);
    legend('off');
    set(p4, 'LineWidth', 2);
    set(findobj('Type','line'),'Color', 'k');
    xlabel('Time from first pulse on (s)'); ylabel('Firing rate'); 
    set(gca, 'Xlim', [0 0.2]);
     

%% New C. Baseline comparisons - try to do all three: non light-responsive, light-responsive (single), light-responsive (double)

% cntrlBase = longCntrlFit(:,4);
oneBase = longNeuronsFit(:,4);
twoBase = longDoubleFit(:,6);

%These are the vecs to be filled for baselines
singleBase = NaN(length(numNeurons), 1);
doubleBase = NaN(length(numNeurons), 1);

for ibase = 1:length(singleFits)'; 
    if singleFits(ibase) == 1
        singleBase(ibase) = oneBase(ibase)
    end
end

for ibase2 = 1:length(doubleFits);
    if doubleFits(ibase2) == 1
        doubleBase(ibase2) = twoBase(ibase2)
    end
end

%Plotting:
%For control neurons:
edgeCntrlBase = [0:8:160]
midsCntrlBase = conv2(edgeCntrlBase, [1 1], 'valid')/2;
%     subplot(3,2,5);
%     hold on;
    latenciesCntrl_b1 = histc(longCntrlFit(:,4), edgeCntrlBase);     
    latenciesCntrl_b = latenciesCntrl_b1(1:(end-1));
    fractionCntrl_b = latenciesCntrl_b./noRes;

%single fit    
edgeBase = [0:6:120] % 20 bins
midsBase = conv2(edgeBase, [1 1], 'valid')/2;
    latencies_b1 = histc(singleBase, edgeBase);     
    latencies_b = latencies_b1(1:(end-1));
    fraction_b = latencies_b./14;                    %divide by singleFit number
   
%double fit    
    latencies_b2 = histc(doubleBase, edgeBase);
    latencies_b2 = latencies_b2(1:(end-1));
    fraction_b2 = latencies_b2./10                  %divide by doubleFit number
    
    subplot(2, 2, 3);
    hold on;
 
    bMat = horzcat(fractionCntrl_b, fraction_b, fraction_b2);
    h6 = bar(midsBase, bMat, 1.5); 
    
    [pC,tblC,statsC] = anova1(bMat);
    
    xlabel('Baseline (average firing rate)'); ylabel('Fraction of neurons'); 
     set(gca, 'Ylim', [0 0.50], 'Ytick', [0 0.25 0.5], 'YtickLabel', [0 0.25 0.5]);
     set(gca, 'Xtick', [min(edgeBase):40:max(edgeBase)], 'Xticklabel', [min(edgeBase):40:max(edgeBase)]);
     set(gca, 'Xlim', [min(edgeBase) max(edgeBase)]);

    hold on;

%% New D = 3 amplitudes: same color for amp1 and 2 (diff shades); single fit diff color
% Find the single-fit amps
oneAmp = longNeuronsFit(:,1);
singleAmpsFit1 = NaN(length(numNeurons), 1);
for a = 1:length(singleFits)'; 
    if singleFits(a) == 1
        singleAmpsFit1(a) = oneAmp(a)
    end
end

firstAmp = longDoubleFit(:,1);
secondAmp = longDoubleFit(:,2);
doubleAmpsFit1 = NaN(length(numNeurons), 1);
doubleAmpsFit2 = NaN(length(numNeurons), 1);
for newi2 = 1:length(doubleFits)'; 
    if doubleFits(newi2) == 1
        doubleAmpsFit1(newi2) = firstAmp(newi2)
    end
    if doubleFits(newi2) == 1
        doubleAmpsFit2(newi2) = secondAmp(newi2)
    end
end

subplot (2,2,4);
% subplot(4,2,4);
hold on;
%     lateCntrls_amp1 = histc(longCntrlFit(:,1), edgeCntrlAmp);               %calculate nums ber bin using edges
%     lateCntrls_amp = lateCntrls_amp1(1:(end-1));                            %make it 20 elements long
%     fractionCntrls_amp = lateCntrls_amp ./ noRes;                           %sum of y axis components should = 1

edgeAmp = [0:12.5:250];    %20 bins                                        
midsAmp = conv2(edgeAmp, [1 1], 'valid') / 2;
    latencies_amp1 = histc(singleAmpsFit1, edgeAmp);                        %calculate nums per bin using edges
    latencies_amp = latencies_amp1(1:(end-1));                              %correct for last bin
    fraction_amp = latencies_amp ./ 14; 
 
    doubAmpsFit1 = histc(doubleAmpsFit1, edgeAmp);
    doubAmpsFit1a = doubAmpsFit1(1:(end-1));
    fraction_doubAmps1 = doubAmpsFit1a ./10;
    
    doubAmpsFit2 =  histc(doubleAmpsFit2, edgeAmp);
    doubAmpsFit2a = doubAmpsFit2(1:(end-1));
    fraction_doubAmps2 = doubAmpsFit2a ./10;
        
    %plot them side-by-side
        fractionAmpMat = horzcat(fraction_amp, fraction_doubAmps1, fraction_doubAmps2);             
        h = bar(midsAmp, fractionAmpMat, 1.5);
%         h = bar(midsCntrlAmp, fractionAmpMat, 1.2);
%         set(h, 'FaceColor', [.5, .5, .5]);
        xlabel('Amplitude (average firing rate)'); ylabel('Fraction of neurons');
        set(gca, 'XLim', [min(edgeAmp) max(edgeAmp)] , 'Ylim', [0 0.5], 'YTick', [0 0.25 0.5]);
