%% Classify PSTHs based upon recording properties - this m file corrects for flipping color map for the untagged neurons
%  in the PC1 v 3 scatter

%  Updated 2.27.18 for new 2017 datasets
%  Updated 1.31.17 to show PC1 ipsi v. PC1 contra (rather than PC1 v PC3 for both
%  specifically load psths aligned to mvment onset
%  Updated 1.24.17, for 151119, 18, 04, 05, 06 recordings (new untagged neurons added)
%  Updated 3.6.18 for finding diffs of ipsi/contras for all sessions
clear;

% cd '/Volumes/My Passport for Mac/concatMvPCANew';
cd /Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15
% cd /Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/older mvPCAs

% S = load('2015_psthMatCorrects4PCA.mat'); %for SFN
S = load('2015_psthMatCorrIncorr')';  %includes incorrect trials 3.6.18

% % [num,txt,raw] = xlsread('2015_psthMatIndsLabel2.xlsx','tpose');  %doesn't read xls
% [num,txt,raw] = xlsread('2015_psthMatIndsLabel.xlsx', 'tpose');  %reads older v. for SFN
[num,txt,raw] = xlsread('2018_psth_MatIndsLabelwOldComb.xlsx', 'tpose');  %reads comb. SFN data + new task phys files

%New for 98 neurons, 3.6.18:
T = load('2018_psthMatCorr_Incorr_matrix.mat');  %generated by the concat file that aligns psths to event
% [num,txt,raw1] = xlsread('2018_psth_MatIndsLabel.xlsx');

% Run these lines to generate the tpose version of xls labelInds = txtMat
% % S = load('2015_psthMatCorrects4PCA_2.mat');  %generated by the concat file that aligns psths to event
% [num,txt,raw1] = xlsread('2015_psthMatIndsLabel_newUntags3.xlsx');
% txtMat = [];
% for j = 1:length(txt)
%     for i = 1:3                                   %depends on the number of variables
%         if ischar(txt{i,j})
%            str(i,j) = string(txt{i,j}); 
% %             str(txtMat(i,j)) = string(txt{i,j}); %nope
% %             str(txtMat(i,j)) = txt{i,j};         %nope
% %             txtMat(str(i,j)) = string(txt{i,j}); %just gives me the same thing
% %             str(txtMat(i,j)) = string;           %nope
%         else break
%         end
%     end
% end
% 
% txtMat = str'; %just put the 3rd col. 1st, copy & paste into the tpose wrksheet.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear num; clear raw1; clear txt;
% [num,txt,raw] = xlsread('2018_psth_MatIndsLabel.xlsx', 'tpose');
% clear num; clear raw1; clear txt;
% [num,txt,raw] = xlsread('2018_psth_MatIndsLabel.xlsx', 'tpose');
% clear num; clear raw1; clear txt;
% [num,txt,raw] = xlsread('2018_psth_MatIndsLabel.xlsx', 'tpose');

psthMatCorr_Incorr4diff = S.psthMatCorr_Incorr4diff;                       %these are the 2015 data
psthMatCorr_Incorr4mat = T.psthMatCorr_Incorr4mat;                         %these are the 2017 data

ST.psthMat = horzcat(T.psthMatCorr_Incorr4mat, S.psthMatCorr_Incorr4diff); %order of the spreadsheet for concat of 2017 & 2015 data
save('ST.psthMat.mat');
% save('S.psthMatCorr_incorr4diff');
% mat = S.psthMatCorr_Incorr4mat; %1st row is ipsi, 2nd is contra for 1 unit; 3rd is ipsi, 4th contra for unit2...
% mat = S.psthMatCorr_Incorr4diff; %same format - older GOOD units from 2015

mat = ST.psthMat;
% mat = T.psthMatCorr_Incorr4mat;
% mat = S.psthMatCorr_Incorr4diff;

matrixA(:,1:size(mat,2)./2) = mat(:,1:2:size(mat,2)-1); %these are the ipsi-mvmnt aligned (at 1000 s) 
matrixB(:,1:size(mat,2)./2) = mat(:,2:2:size(mat,2));   %these are the contra-mvmnt aligned (at 1000 s)

matrixIpsi2 = matrixA;
matrixContra2 = matrixB;

save('matrixIpsi2.mat');
save('matrixContra2.mat');

figure(1); 
subplot(131); plot(matrixA(:,1)'); title('u1 mvAlign; ipsi');  %only for unit 1
subplot(132); plot(matrixB(:,1)'); title('u1 mvAlign; contra');
subplot(133); plot((matrixA(:,1)'- matrixB(:,1)').^2); title('u1; ipsi-contra');

% figure(); subplot(121); plot(mat(:,1)'); subplot(122); plot(mat(:,2)')
figure(2); subplot(131); imagesc(matrixA'); title('mvAlign ipsi mat');
subplot(132); imagesc(matrixB'); title('mvAlign contra mat');
subplot(133); imagesc( (matrixA'- matrixB') , [-0.1 0.1]); title('ipsi-contra mat'); colorbar;

%Do SNr neurons show directional (ipsi/contra) tuning? - shows diffs in mvOnset and at reward
tmp = abs( matrixA' - matrixB' ); 
figure(4); clf; subplot(121); plot(mean(mat',1)); xlim([0 3200]);title('avg of all units, mvAlign');
subplot(122); plot(sum(tmp,1)); xlim([0 3200]); %ylim([.010 .0170]); 
title('sum of abs(ipsi-contra)');

figure; hold on; plot(tmp)

%choose example units that show greater ipsi v contra activity: 170111 units 16 & 22

% [num,txt,raw] = xlsread('2015_psthMatIndsLabel_newUntags3.xlsx', 'tpose');
% raw2 = ~nan(raw);


% This is still for generating PCA 3.8.18
for j=1:size(raw,1)
    
    if strfind(raw{j,1},'ipsi')     %if ipsi, write 1 in first col
        taskLogic(j,1)=1;
    else
        taskLogic(j,1)=0;
    end

    if strfind(raw{j,2},'untagged') %if untagged, write 1 in 2nd col
        taskLogic(j,2)=1;           %will be 144 of them; 62 untagged
    else if strfind(raw{j,2}, 'tag')
        taskLogic(j,2)=0;
        end
    end

    if strfind(raw{j,3},'blocks')   %if blocks task, write 1 in 3rd col
        taskLogic(j,3)=1;
    else
        taskLogic(j,3)=0;        
    end
    
    %assign color to the possible 3-var combos: (but this isn't called again)
    %ipsi  tagged    blocks = 5, contra tagged   blocks = 1;
    %ipsi  untagged  blocks = 7, contra untagged blocks = 3;
    %ipsi  tagged    rand   = 4, contra tagged   rand   = 0;
    %ipsi  untagged  rand   = 2, contra untagged rand =   6;
    taskColor(j,1) = base2dec([num2str(taskLogic(j,1)) num2str(taskLogic(j,2)) num2str(taskLogic(j,3))],2);
    
end
labels{1} = 'ipsi=cyan';
labels{2} = 'untagged=cyan';
labels{3} = 'blocks=cyan';
labels2{1} = 'PC1';
labels2{2} = 'PC2' ;
labels2{3} = 'PC3';
labels3{1} = 'ipsi - contra (per neuron)';
labels3{2} = 'untagged=cyan';
labels3{3} = 'blocks=cyan';

%This function returns the low-dimensional representation of the data in the matrix; new mat = mappedA.
%Rows are observations (here: 49 neurons), cols are 1stPCA, 2nd, 3rd:

%mappedA = compute_mapping(psth matrix, type, no_dims, parameters)
%Info on the mapping is returned in the struct "mapping"
% [mappedA,mapping] = compute_mapping(S.psthMatCorrects4PCA','PCA',3);
% [mappedA,mapping] = compute_mapping(S.psthMatCorrects4PCA_2','PCA',3);
[mappedA,mapping] = compute_mapping(S.psthMatCorrects4PCAnew','PCA',3);

mapName = TNC_CreateRBColormap(2,'bo');    %from TONIC:  [mapName] = TNC_CreateRBColormap(numValues,type)
mapName2 = [.5, .5, .5; 0, 0.5, 1];

for jj=1:3
%     figure(1); 
%     figure(jj); subplot(2,3,jj);
    if jj == 3
        figure(jj);
        h = scatter(mappedA(:,1),mappedA(:,3),100,taskLogic(:,jj),'filled', 'MarkerEdgeColor', [0 0 0]); %colormap(3 col matrix = RGB)   
        c = colormap(mapName); 
        xlabel('PC1', 'FontSize', 24); ylabel('PC3', 'FontSize', 24);
        title(labels(jj));
            ax = gca; 
            ax.FontSize = 24;
        hold off;   
        %for all of the 1's, cyan should be plotted, but why doees it flip the color map for the middle-top plot, the colors are flipped?
    end
    
    if jj == 1
        figure(jj);
        mappedA1 = mappedA(:,1); taskLogicIpsi = taskLogic(:,1); taskLogicContra = taskLogic(:,1) == 0;
        mappedIpsi = mappedA1(taskLogicIpsi==1); mappedContra = mappedA1(taskLogicContra == 1);

        %scatter(x,y, size, CData = color);
%         h1 = scatter(mappedA1(taskLogicIpsi == 1), mappedA1(taskLogicIpsi == 0), 100, taskLogic(:,jj), 'filled'); %doesn't work bc taskLogic is the entire 206x1 vector & others are only 1/2
        h1 = scatter(mappedIpsi, mappedContra, 100, 'filled', 'MarkerEdgeColor', [0 0 0]);  %plot PC1 ipsi v PC1 contra
        c = colormap(mapName);
%         h1 = plot(mappedIpsi, mappedContra, '.');  %plot PC1 ipsi v PC1 contra
%         set(h1, 'markerEdgeColor', [0 0 0], 'MarkerSize', '100')

%        group = mappedA1;
%        h1 = gscatter(mappedA1(taskLogicIpsi==1), mappedA1(taskLogicIpsi==0), group, 'br', 'oo', 18);
%         group = mappedIpsi;
%         h1 = gscatter(mappedIpsi, mappedContra, group, 'br', 'oo', 18, 'doleg', 'off');  %plot PC1 ipsi v PC1 contra
        xlabel('PC1 ipsi', 'FontSize', 24); ylabel('PC1 contra', 'FontSize', 24);
        title(labels(jj));
            ax = gca; 
            ax.FontSize = 24;
        hold off;  
    end
    %Make a scatter plot of PCA1 v. PCA3
    if jj == 2
        figure(jj);
        h2 = scatter(mappedA(:,1),mappedA(:,3),100,taskLogic(:,jj),'filled', 'MarkerEdgeColor', [0 0 0]); %colormap(3 col matrix = RGB)
        c2 = colormap(mapName2); 
        xlabel('PC1', 'FontSize', 24); ylabel('PC3', 'FontSize', 24);
        title(labels(jj));
            ax = gca; 
            ax.FontSize = 24;
        
       figure; plot(mapping.M(:,jj)); box off;  %this is the PC of interest for psth's aligned to mvmt onset
        title(labels2(jj));
        xlim([0 3100]); 
        xlabel('mvmt onset = 1000ms');
        hold off;
    end
%     xlabel('PC1'); ylabel('PC3');
%     title(labels(jj));
%     figure(1); subplot(2,3,jj+3);
%     figure(jj); subplot(2,3,jj+3);
end
for jj = 1:3
    figure(jj+3)
    plot(mapping.M(:,jj)); box off; %this is the PC of interest for psth's aligned to mvmt onset
    title(labels2(jj));
    xlim([0 3100]); 
    xlabel('mvmt onset = 1000ms'); 
    
end


%% Plot the marginal histograms of the projections onto PC of interest:
for jj=1:3                              %for ea. column of variables
    posInds = find(taskLogic(:,jj)==1); %1 = ipsi (1st col), or untagged (2nd), or blocks task (3rd)
    negInds = find(taskLogic(:,jj)==0); %0 = contra, tagged, rand task
    
    if jj==1
        figure(7); subplot(3,3,jj); 
%       figure(5); subplot(3,3,jj); 
            %[counts,centers] = hist(x, nbins); %Plot the marginal histograms for the data shown in Figure 1, top panel
            [tmpP,x] = hist(mappedA(posInds,1)-mappedA(negInds,1),-1.5:0.05:0.5); %PC1
            hold off; plot(x,tmpP,'-','LineWidth',2,'color',mapName(1,:));
            hold on; plot([0 0],[0 6],'k-'); box off;   %draw vertical line
            xlim([-1.5 0.5]); ylim([0 20]);
            xlabel('Loadings onto PC1', 'FontSize', 16); ylabel('# Neurons', 'FontSize', 16);
            title(labels3{jj});
         figure(7); subplot(3,3,jj+3); 
%        figure(6); subplot(3,3,jj+3); 

            [tmpP,x] = hist(mappedA(posInds,2)-mappedA(negInds,2),-0.2:0.01:0.2); %PC2
            hold off; plot(x,tmpP,'-','LineWidth',2,'color',mapName(1,:));        
            hold on; plot([0 0],[0 6],'k-'); box off;
            xlabel('Loadings onto PC2', 'FontSize', 16); ylabel('# Neurons', 'Fontsize', 16);

        figure(7); subplot(3,3,jj+6); 
%                 figure(7); subplot(3,3,jj+6); 

            [tmpP,x] = hist(mappedA(posInds,3)-mappedA(negInds,3),-0.2:0.01:0.2);  %PC3
            hold off; plot(x,tmpP,'-','LineWidth',2,'color',mapName(1,:));        
            hold on; plot([0 0],[0 6],'k-'); box off;
            xlabel('Loadings onto PC3','Fontsize', 16); ylabel('# Neurons', 'Fontsize', 16);
            ylim([0 20]);

    else
        [tmpP,x] = hist(mappedA(posInds,1),-1.5:0.05:0.5);
        [tmpN,xx] = hist(mappedA(negInds,1),-1.5:0.05:0.5);
        figure(7); subplot(3,3,jj); 
        hold off; plot(x,tmpP,'-','LineWidth',2,'color',mapName(1,:)); box off;%posInds of mapped A will be blue
        hold on;  plot(xx,tmpN,'-','LineWidth',2,'color',mapName(2,:)); box off;%negInds of mapped A will be gray
        axis tight;
        title(labels3{jj});
        ylim([0 20]);
        
        [tmpP,x] = hist(mappedA(posInds,2),-0.25:0.01:0.25);
        [tmpN,xx] = hist(mappedA(negInds,2),-0.25:0.01:0.25);
        figure(7); subplot(3,3,jj+3); 
        hold off; plot(x,tmpP,'-','LineWidth',2,'color',mapName(1,:)); box off;
        hold on;  plot(xx,tmpN,'-','LineWidth',2,'color',mapName(2,:)); box off;
        axis tight;
        ylim([0 20]);
        
        [tmpP,x] = hist(mappedA(posInds,3),-0.2:0.01:0.2);
        [tmpN,xx] = hist(mappedA(negInds,3),-0.2:0.01:0.2);
        figure(7); subplot(3,3,jj+6); 
        hold off; plot(x,tmpP,'-','LineWidth',2,'color',mapName(1,:)); box off;
        hold on;  plot(xx,tmpN,'-','LineWidth',2,'color',mapName(2,:)); box off;
        axis tight;
        ylim([0 20]);
    end
    
end