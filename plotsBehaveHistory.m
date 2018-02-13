%% plotsBehaveHistory
% This calls the blocks and rand mfiles for plotting on one graph
%% Initialize for the blocks file:

% fileDirectory = '/Users/stubblefielde/Desktop/behaviorFebMarApr16/Vgatsix/forAnalysis';
fileDirectory = '/Users/stubblefielde/Desktop/behaviorFebMarApr16/plots';

% posStartsTrueMaxIR = 400;
posStartsTrueMaxIR = 510;
posStartsTrueMax = 200;
posStartsTrueMin = 190;
posStartsTrueMinNeg = -190;
posStartsTrueMaxNeg = -200; 
% posStartsTrueMinNegIL = -400;
posStartsTrueMinNegIL = -510;

%blocked trial constraints:
% trialConstrainMin = 60;
% trialConstrainMax = 170; %best performance within 100ish trials
trialConstrainMin = 77;
trialConstrainMax = 97; %best performance within 100ish trials for L block

% [Prob_Li, Prob_LiMinusOne, Prob_LiMinusTwo, Prob_LiMinusThree, Prob_LiMinusFour, Prob_LiMinusFive, Prob_LiMinusSix, Prob_LiMinusSeven] = behaveLoad5history_blocks(fileDirectory)
[Prob_Li, Prob_LiMinusOne1, Prob_LiMinusTwo1, Prob_LiMinusThree1, Prob_LiMinusFour1, Prob_LiMinusFive1, Prob_LiMinusSix1, Prob_LiMinusSeven1] = behaveLoad5history_blocks(fileDirectory)

%% Blocks: Plot Prob(left trial i) as a function of trial history:
figure; hold on;
%     y = [Prob_Li, Prob_LiMinusOne, Prob_LiMinusTwo, Prob_LiMinusThree, Prob_LiMinusFour, Prob_LiMinusFive, Prob_LiMinusSix, Prob_LiMinusSeven];
Prob_Li_jitter = Prob_Li + 0.01;
y = [Prob_Li_jitter, Prob_LiMinusOne1, Prob_LiMinusTwo1, Prob_LiMinusThree1, Prob_LiMinusFour1, Prob_LiMinusFive1, Prob_LiMinusSix1, Prob_LiMinusSeven1];
h = plot(y, 'k.', 'MarkerSize', 25);
%     ylabel('P choice (Left)', 'FontSize', 20);  %Divided ea. by Prob_Li, so yes, P(L choice) on trial i

    ylabel('P choice (Left_i)', 'FontSize', 22);  
    xlabel('Trial history (i-num trials back) for 1 session', 'FontSize', 22); 
    ylim([0 1.1]);
%     xlabels = single([0 1 2 3 4 5 6 7]);
%     xlabels1 = sort(xlabels, 'descend'); 
%     xlabels1 = single([i i-1 i-2 i-3 i-4 i-5 i-6 i-7]); 
%      set(gca, 'YTick', 0:20:80, 'XTickLabels', xlabels1);
     ax = gca; 
     ax.FontSize = 22;


clear all;

%% Initialize for the randomized file:

fileDirectory = '/Users/stubblefielde/Desktop/behaviorFebMarApr16/Vgatfive/forAnalysis';


posStartsTrueMaxIR = 400;
% posStartsTrueMaxIR = 510;
posStartsTrueMax = 200;
posStartsTrueMin = 190;
posStartsTrueMinNeg = -190;
posStartsTrueMaxNeg = -200; 
posStartsTrueMinNegIL = -400;
% posStartsTrueMinNegIL = -510;

%rand trial constraints:
% trialConstrainMin = 30;
% trialConstrainMax = 190; %best performance within 100ish trials


% [Prob_Li, Prob_LiMinusOne, Prob_LiMinusTwo, Prob_LiMinusThree, Prob_LiMinusFour, Prob_LiMinusFive, Prob_LiMinusSix, Prob_LiMinusSeven] = behaveLoad5history_rand(fileDirectory)
[Prob_Li, Prob_LiMinusOne1, Prob_LiMinusTwo1, Prob_LiMinusThree1, Prob_LiMinusFour1, Prob_LiMinusFive1, Prob_LiMinusSix1, Prob_LiMinusSeven1] = behaveLoad5history_rand(fileDirectory)

Prob_Li_jitter = Prob_Li - 0.01;

%% Rand: Plot Prob(left trial i) as a function of trial history: 

gcf; 
hold on; 
%     y = [Prob_Li, Prob_LiMinusOne, Prob_LiMinusTwo, Prob_LiMinusThree, Prob_LiMinusFour, Prob_LiMinusFive, Prob_LiMinusSix, Prob_LiMinusSeven];
y = [Prob_Li_jitter, Prob_LiMinusOne1, Prob_LiMinusTwo1, Prob_LiMinusThree1, Prob_LiMinusFour1, Prob_LiMinusFive1, Prob_LiMinusSix1, Prob_LiMinusSeven1];        

%     y = [Prob_Li, Prob_Li, Prob_Li, Prob_Li, Prob_Li, Prob_Li, Prob_Li, Prob_Li]; %should get at (Prob choice Li) for Y axis
h = plot(y, 'r.', 'MarkerSize', 25, 'Color', [0.5 0.5 0.5]);
%    ylabel('P choice (Left)', 'FontSize', 20);  %Divided ea. by Prob_Li, so yes, P(L choice) on trial i
    ylabel('P choice (Left_i)', 'FontSize', 22);  
    xlabel('Trial history (i-num trials back)', 'FontSize', 22); 
    ylim([0 1.05]);
%     xlabels = single([0 1 2 3 4 5 6 7]);
%     xlabels1 = sort(xlabels, 'descend'); 
%     xlabels1 = single([i i-1 i-2 i-3 i-4 i-5 i-6 i-7]); 
%     set(gca, 'YTick', 0:.20:1.2, 'XTickLabels', xlabels1);
    set(gca, 'YTick', 0:.20:1.2, 'XTick', 0:1:8);
    ax = gca; 
    set(ax);
    ax.XTickLabelMode = 'manual'
    ax.XTickLabel = {' ', 'i','i-1','i-2','i-3', 'i-4','i-5','i-6','i-7'}
    ax.FontSize = 22;


