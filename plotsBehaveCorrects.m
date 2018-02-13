%% plotsBehaveHistory
% This calls the blocks and rand mfiles for plotting on one graph
%% Initialize for the blocks file:

%7.05 blocks
%7.12 blocks
%7.28 blocks
%7.29 blocks
% Used this and behaveLoad5corrects to arrive at % corrects for the following mice/sessions.

%Blocks:
%Vgatsix:
percCorrBlocks1 = 74;
percCorrBlocks3 = 72;
percCorrBlocks4 = 62;
percCorrBlocks5 = 68; %8/09
% percCorrectsBlocks6 = 
% percCorrectsBlocks7 = 

%Vgatfour:
percCorrBlocks8 = 65;
percCorrBlocks9 = 70;
percCorrBlocks10 = 56;

%Vgattwo:
percCorrBlocks11 = 68;
percCorrBlocks12 = 70;
percCorrBlocks13 = 64;

Vgatsix = horzcat(percCorrBlocks1, percCorrBlocks3, percCorrBlocks4, percCorrBlocks5);
avgBlockNew = mean(Vgatsix);
stdBlockNew = std(Vgatsix);
semBlockNew = stdBlockNew/sqrt(numel(Vgatsix));
E1 = semBlockNew;

VgatBlockOlds = horzcat(percCorrBlocks8, percCorrBlocks9, percCorrBlocks10, percCorrBlocks11, percCorrBlocks12, percCorrBlocks13);
avgBlockOlds = mean(VgatBlockOlds);
stdBlockOlds = std(VgatBlockOlds);
semBlockOlds = stdBlockOlds/sqrt(numel(VgatBlockOlds));
E2 = semBlockOlds;

%% Randomized sessions for percent corrects:
%Vgatthree (10.26-11.05)
%randomized
percCorrRands1 = 49;
percCorrRands2 = 47;
percCorrRands3 = 50;
percCorrRands4 = 46; 
percCorrRands5 = 51; 

%Vgatfive (10.26-11.05)
percCorrRands6 = 63;
percCorrRands7 = 61;
percCorrRands8 = 69;
percCorrRands9 = 62; 
% percCorrRands5 = 51; 


VgatRandOlds = horzcat(percCorrRands1, percCorrRands2, percCorrRands3, percCorrRands4, percCorrRands5);
avgRandOlds = mean(VgatRandOlds);
stdRandOlds = std(VgatRandOlds);
semRandOlds = stdRandOlds/sqrt(numel(VgatRandOlds));
E3 = semRandOlds;

VgatRandNews = horzcat(percCorrRands6, percCorrRands7, percCorrRands8, percCorrRands9);
avgRandNews = mean(VgatRandNews);
stdRandNews = std(VgatRandNews);
semRandNews = stdRandNews/sqrt(numel(VgatRandNews));
E4 = semRandNews;

%% Rand: Plot Prob(left trial i) as a function of trial history: 

figure; 
hold on; 
y1 = avgBlockOlds; y2 = avgBlockNew; y3 = avgRandOlds; y4 = avgRandNews;
x1 = 1; x2 = 2; x3 = 4; x4 = 5;
h1 = errorbar(x1, y1, E1, 's', 'MarkerSize', 7, 'Color', [0 0 0], 'LineWidth', 1);
h2 = errorbar(x2, y2, E2, 'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'LineWidth', 1);
h3 = errorbar(x3, y3, E3, 'o', 'MarkerSize', 7, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
h4 = errorbar(x4, y4, E4, 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 1);
    ylabel('Percent correct', 'FontSize', 20);  
%     xlabel('Blocks', 'FontSize', 22); 
    xlim([0 6]); ylim([0 100]);
    set(gca, 'YTick', 0:10:100, 'XTick', 0:1:6);
    ax = gca; 
    set(ax);
    ax.XTickLabelMode = 'manual'
    ax.XTickLabel = {' ', 'Old Blocks','New', ' ', 'Old rand.', 'New'}
    ax.FontSize = 20;


