%% Classify PSTHs based upon recording properties
%  specifically load psths aligned to mvment onset
%  Updated 1.24.17, for 151119, 18, 04, 05, 06 recordings (new untagged neurons added)

clear;
cd /Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15
% cd /Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/older mvPCAs
% S = load('2015_psthMatCorrects4PCA.mat');
% % [num,txt,raw] = xlsread('2015_psthMatIndsLabel2.xlsx','tpose');  %doesn't read xls
% [num,txt,raw] = xlsread('2015_psthMatIndsLabel.xlsx', 'tpose');  %reads older v. for SFN

%New for 103 neurons, 1.24.18:
S = load('2015_psthMatCorrects4PCA_2.mat');
[num,txt,raw1] = xlsread('2015_psthMatIndsLabel_newUntags2.xlsx');
txtMat = [];
for j = 1:length(txt)
    for i = 1:3;
        if ischar(txt{i,j})
           str(i,j) = string(txt{i,j}); 
%             str(txtMat(i,j)) = string(txt{i,j}); %nope
%             str(txtMat(i,j)) = txt{i,j};         %nope
%             txtMat(str(i,j)) = string(txt{i,j}); %just gives me the same thing
%             str(txtMat(i,j)) = string;           %nope
        else break
        end
    end
end
txtMat = str'; %just put the 3rd col. 1st, copy & paste into the tpose wrksheet.
%xlswrite(filename,A,sheet) writes to the specified worksheet.
% xlswrite('2015_psthMatIndsLabel_newUntags2.xlsx',txtMat,'tpose') %writes to the specified worksheet.
[num,txt,raw] = xlsread('2015_psthMatIndsLabel_newUntags3.xlsx', 'tpose');
% raw2 = ~nan(raw);
for j=1:size(raw,1)
    
    if strfind(raw{j,1},'ipsi')     %if ipsi, write 1 in first col
        taskLogic(j,1)=1;
    else
        taskLogic(j,1)=0;
    end

    if strfind(raw{j,2},'untagged') %if untagged, write 1 in 2nd col
        taskLogic(j,2)=1;           %will be 144 of them; 62 untagged
    else    
        taskLogic(j,2)=0;
    end

    if strfind(raw{j,3},'blocks')   %if blocks task, write 1 in 3rd col
        taskLogic(j,3)=1;
    else
        taskLogic(j,3)=0;        
    end
    
    %assign color to the possible 3-var combos: 
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
labels3{1} = 'IPSI - CONTRA';
labels3{2} = 'untagged=cyan';
labels3{3} = 'blocks=cyan';

%This function returns the low-dimensional representation of the data in the matrix; new mat = mappedA.
%Rows are observations (here: 49 neurons), cols are 1stPCA, 2nd, 3rd:

%mappedA = compute_mapping(psth matrix, type, no_dims, parameters)
%Info on the mapping is returned in the struct "mapping"
% [mappedA,mapping] = compute_mapping(S.psthMatCorrects4PCA','PCA',3);
[mappedA,mapping] = compute_mapping(S.psthMatCorrects4PCA_2','PCA',3);

mapName = TNC_CreateRBColormap(2,'bo');    %from TONIC:  [mapName] = TNC_CreateRBColormap(numValues,type)
% mapName2 = [0 1 1; .5 .5 .5];
for jj=1:3
    figure(1); subplot(2,3,jj);
%      figure(jj); subplot(2,3,jj);

    %Make a scatter plot of PCA1 v. PCA3
    h = scatter(mappedA(:,1),mappedA(:,3),50,taskLogic(:,jj),'filled', 'MarkerEdgeColor', [0 0 0]); %colormap(3 col matrix = RGB)

%      scatter(mappedA(:,1),mappedA(:,3),50,taskLogic(jj,:),'filled'); 
    c = colormap(mapName);                  %for all of the 1's, cyan should be plotted, but why doees it seem like for the middle-top plot, the colors are flipped?
%      colormap(mapName2);                  %THEY ARE FLIPPED

    xlabel('PC1'); ylabel('PC3');
    title(labels(jj));
    figure(1); subplot(2,3,jj+3);
%     figure(jj); subplot(2,3,jj+3);

    plot(mapping.M(:,jj));  %this is the PC of interest for psth's aligned to mvmt onset
%     title(labels2(jj));
    title(labels2(jj));
    xlim([0 3100]); 
    xlabel('mvmt onset = 1000ms'); 
    
end

for jj=1:3                              %for ea. column of variables
    posInds = find(taskLogic(:,jj)==1); %1 = ipsi (1st col), or untagged (2nd), or blocks task (3rd)
    negInds = find(taskLogic(:,jj)==0); %0 = contra, tagged, rand task
    
    if jj==1
        figure(2); subplot(3,3,jj); 
%       figure(5); subplot(3,3,jj); 
            %[counts,centers] = hist(x, nbins); %Plot the marginal
            %histograms for the data shown in Figure 1, top panel
            [tmpP,x] = hist(mappedA(posInds,1)-mappedA(negInds,1),-0.4:0.02:0.4);
            hold off; plot(x,tmpP,'-','LineWidth',2,'color',mapName(1,:));
            hold on; plot([0 0],[0 6],'k-');
            xlabel('Loadings onto PC1'); ylabel('# Neurons');
            title(labels3{jj});
         figure(2); subplot(3,3,jj+3); 
%        figure(6); subplot(3,3,jj+3); 

            [tmpP,x] = hist(mappedA(posInds,2)-mappedA(negInds,2),-0.2:0.01:0.2);
            hold off; plot(x,tmpP,'-','LineWidth',2,'color',mapName(1,:));        
            hold on; plot([0 0],[0 6],'k-');
            xlabel('Loadings onto PC2'); ylabel('# Neurons');

        figure(2); subplot(3,3,jj+6); 
%                 figure(7); subplot(3,3,jj+6); 

            [tmpP,x] = hist(mappedA(posInds,3)-mappedA(negInds,3),-0.2:0.01:0.2);
            hold off; plot(x,tmpP,'-','LineWidth',2,'color',mapName(1,:));        
            hold on; plot([0 0],[0 6],'k-');
            xlabel('Loadings onto PC3'); ylabel('# Neurons');

    else
        [tmpP,x] = hist(mappedA(posInds,1),-0.5:0.05:0.5);
        [tmpN,xx] = hist(mappedA(negInds,1),-0.5:0.05:0.5);
        figure(2); subplot(3,3,jj); 
        hold off; plot(x,tmpP,'-','LineWidth',2,'color',mapName(1,:)); %posInds of mapped A will be blue
        hold on;  plot(xx,tmpN,'-','LineWidth',2,'color',mapName(2,:));%negInds of mapped A will be gray
        axis tight;
        title(labels3{jj});
        
        [tmpP,x] = hist(mappedA(posInds,2),-0.25:0.01:0.25);
        [tmpN,xx] = hist(mappedA(negInds,2),-0.25:0.01:0.25);
        figure(2); subplot(3,3,jj+3); 
        hold off; plot(x,tmpP,'-','LineWidth',2,'color',mapName(1,:));
        hold on;  plot(xx,tmpN,'-','LineWidth',2,'color',mapName(2,:));
        axis tight;
        
        [tmpP,x] = hist(mappedA(posInds,3),-0.2:0.01:0.2);
        [tmpN,xx] = hist(mappedA(negInds,3),-0.2:0.01:0.2);
        figure(2); subplot(3,3,jj+6); 
        hold off; plot(x,tmpP,'-','LineWidth',2,'color',mapName(1,:));
        hold on;  plot(xx,tmpN,'-','LineWidth',2,'color',mapName(2,:));
        axis tight;
    end
    
end