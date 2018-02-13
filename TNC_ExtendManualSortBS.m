% _________________________________________________________________________
% BS updated 1/4/16 - l. 138: index exceeds matrix dimensions (due to diff seg file vs. template)
% BS updated 1/7/16 - l. 190: user must indicate bounds of the pairwise corr w/ template: tried clicking on 0.8, 2000
% BS updated 1/6/15 - l. 241: underfined func or var "target_values"
% BS updated 1/7/16 - after defining bounds of corr., a _shank_tsd.mat file is generated: targets are in black; template is in blue
% BS updated 1/10/15 - ll. 98, 249 generate as many figs as units attmepting to template-match rather than drawing over fig.1
% ________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / Janelia
% 
% BUG REPORTING: josh@dudmanlab.org
% FURTHER INFORMATION: www.dudmanlab.org
% Copyright (C) 2015 by Howard Hughes Medical Institute.
% _________________________________________________________________________
% INPUTS:
% segArray: array of segments that sorting will be applied over
% seed_segs: array of segments that should be used to generate the template
% shank_num: what shank to perform sorting on
% file_Path: basename of the files used ('_ft_tns.mat' '_ss.mat' '_ft.mat' are required)
% updateTnsFlag: 1=overwrite tns file with new sorted indices
% 
% EXAMPLE:
% [shank] = TNC_ExtendManualSort([1:2],[1],5,'DA_F06_20140815-i16',1);
% 
% OUTPUTS:
% Can overwrite the spike ids file (*_tns) used for sorting in GUI

%% Modified BS 12.20.15

function [shank] = TNC_ExtendManualSortBS(segArray,seed_segs,shank_num,file_Path,updateTnsFlag)

% % valid seg range:
% segArray = [1:10]; %first-pass testing with the initial 10chunks that I've extracted (1/5/16)
% seed_segs = [1];
% shank_num = 1;
updateTnsFlag = 0;

segArray = [2];     %first-pass testing with the initial 10chunks that I've extracted (1/5/16)
seed_segs = [1];
shank_num = 1;

% file_PathApply = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree'                                 %will have "chunk" in the name
% file_PathSeed =  '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/sorted/1maxSeg/151104ipsi/1001'  %will have "tag" in the name

file_Path = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/wfManSortTest/' %Put them all in same folder
% segArray = [1:40]; %for Vgatthree 151104 recording; apply over all of 151104chunk2002 segs
% seed_segs = [1];   %for Vgatthree 151104 recording; use 151104tag1001, 1st segment (thus, diff file names - figure out how to incorporate)
% shank_num = 1;     %might as well start on 1

% paths to sessStruct and ft_tns files:
% tns_Path              = [file_Path '_ft_tns.mat']; 
% sessionStruct_Path    = [file_Path '_ss.mat'];
% featData_Path         = [file_Path '_ft.mat'];

% tns_Path              = [file_Path '151105tag1001.ns5_ft_tns.mat'];          %already sorted - use these units' waveforms for template 
% sessionStruct_Path    = [file_Path '151105_90chunk1_10_1002.ns5_ss.mat'];    %not sorted - for matching
% featData_Path         = [file_Path '151105_90chunk1_10_1002.ns5_ft.mat'];    %not sorted - for matching

tns_Path              = [file_Path '151105_90chunk1_3_1002.ns5_ft_tns.mat'];   %already sorted - use these units' waveforms for template 
sessionStruct_Path    = [file_Path '151105_90chunk1_3_1002.ns5_ss.mat'];       %not sorted - for matching
featData_Path         = [file_Path '151105_90chunk1_3_1002.ns5_ft.mat'];       %not sorted - for matching


S = load(sessionStruct_Path);    %sessionStruct has waveforms and inds per seg
sessionStruct=S.sessionStruct;
clear S

S = load(tns_Path);              %this idList consists of cluster #s
idList=S.idList;
clear S

S = load(featData_Path);         %featStruct consists of ids, params (but what do these #s mean?) of ea. shank and seg
featStruct=S.featStruct;
clear S

disp(' ');
disp(['1) Extracting data from shank ' num2str(shank_num) ' using templates defined from segments [' num2str(seed_segs) '] <<< ' file_Path]);

%% Define the template for current unit

% find the number of units in the seed_segs
unit_list_tmp=[];
for i=1:numel(seed_segs)
    unit_list_tmp = [unit_list_tmp unique(idList.seg(seed_segs(1)).shank(shank_num).id)];
end

unit_list = unique(unit_list_tmp(find(unit_list_tmp>0)));                                   %retains u#s from clusters in TNC_SS_GUI

disp(['2) Found ' num2str(numel(unit_list)) ' templates <<< [' num2str(unit_list') ']']);

%% Calculate the inner products for all detected events within the segment range

disp(['3) Extending manually sorted units across all segments...']);

for j=1:numel(unit_list)

%     figure(1); clf;
    figure(j);
    currUnit = unit_list(j);
    
    fprintf('\r');
    disp(['::: Extracting the waveform template for ' num2str(currUnit)]);
    % calculate the template    
    total_events = 0;
    
    for p=1:numel(seed_segs)
        currSeg = seed_segs(p);
        
        ids_for_curr_unit = find(idList.seg(currSeg).shank(shank_num).id==currUnit);   %ids for curr unit are specifically fitting unit 1 to next seg (first pass)
        
        for m=1:numel(ids_for_curr_unit)
            if m==1 & p==1
                template_values = single(sessionStruct.seg(currSeg).shank(shank_num).wfs(ids_for_curr_unit(m)).values);                
            else
                template_values = template_values + single(sessionStruct.seg(currSeg).shank(shank_num).wfs(ids_for_curr_unit(m)).values);
            end
        end
        
        total_events = total_events+numel(ids_for_curr_unit);                   

    end
                                                                                
    % find the 4 channels with the largest amplitudes
    template_amps = min(template_values,[],2);                                  %Why with only 2 clusters are there more than 2 waveforms? (see l. 128)
    [values,indices] = sort( template_amps , 'ascend' );                        %[x,y] = ascending order of template amps for plotting
    chan_per_shank = numel(indices);                                            %each index here is a channel; ea. row of template_values is a channel
    
    % NEED TO DEAL WITH VARIABLE SHANK SIZES HERE
    if chan_per_shank >= 3
        template_vector = [template_values(indices(1),:) template_values(indices(2),:) template_values(indices(3),:)] ./ total_events;
        template_amps = [min(template_values(indices(1),:)) min(template_values(indices(2),:)) min(template_values(indices(3),:))] ./ total_events;        
    else
        template_vector = template_values(indices(1),:) ./ total_events;            %template vector is the 1st avg. waveform (template_values for 1st ind. element)
        template_amps = [min(template_values(indices(1),:))] ./ total_events;        
    end

    % display a quick sanity check
    figure(1); subplot(1,2,1); plot(template_vector,'b','LineWidth',3,'Color',[0.5 0.8 1]); title(['Extracting... ' num2str(currUnit)]); drawnow;
    
    
    fprintf('::: Calculating a heuristic for template match for all segments ');
    for k=1:numel(segArray)
        
        currSeg = segArray(k);
        fprintf('|');
        
        if k==1
            cat_heur = [];
            cat_heurD= [];
        end
        
        heuristics.seg(k).c = zeros(1,numel(sessionStruct.seg(currSeg).shank(shank_num).wfs));    %index exceeds matrix dimensions here for seedSegs of diff file from segArray
        heuristics.seg(k).d = zeros(1,numel(sessionStruct.seg(currSeg).shank(shank_num).wfs));
        
        for m=1:numel(sessionStruct.seg(currSeg).shank(shank_num).wfs)

            tmp = single(sessionStruct.seg(currSeg).shank(shank_num).wfs(m).values);
            
            if chan_per_shank >= 3        
                test_vector = [ tmp(indices(1),:) tmp(indices(2),:) tmp(indices(3),:) ];
                test_amps = [min(tmp(indices(1),:)) min(tmp(indices(2),:)) min(tmp(indices(3),:)) ];
                R = corrcoef(template_vector',test_vector');
                heuristics.seg(k).c(1,m) = R(2,1);
                heuristics.seg(k).d(1,m) = sum( abs(template_amps-test_amps) ) ;                
            else
                test_vector = tmp(indices(1),:) ;
                test_amps = min(tmp(indices(1),:)) ;
                R = corrcoef(template_vector',test_vector');
                heuristics.seg(k).c(1,m) = R(2,1);
                heuristics.seg(k).d(1,m) = sum( abs(template_amps-test_amps) ) ;                
            end

            
        end
        
        cat_heur = [cat_heur heuristics.seg(k).c];
        cat_heurD= [cat_heurD heuristics.seg(k).d];
        
    end
    fprintf('\r');

    [count,values] = hist(cat_heur,2000);
    [countD,valuesD] = hist(cat_heurD,2000);
    
    % display a histogram of projection (correlation) values for thresholding
%     figure(1); subplot(1,2,2); hold off; hist(cat_heur,2000); drawnow;
    figure(1); subplot(1,2,2); hold off; plot(cat_heur,cat_heurD,'k.','MarkerSize',1,'Color',[0.25 0.25 0.25]); axis([0 1 0 max(cat_heurD)]);  
    xlabel('Pairwise correlation with template'); ylabel('Summed error (uV)');
    drawnow;
    
    % ask user to specify the bounds  
    [userThresh,userThreshY] = ginput(1);
    figure(1); subplot(1,2,2); hold on; plot([userThresh userThresh],[0 max(cat_heurD)],'r',[0 1],[userThreshY userThreshY],'r','LineWidth',2); drawnow;

    disp(['::: A threshold of ' num2str(userThresh) ' was chosen by the user']);

%             figure(2); subplot(121); plot(cat_heur,cat_heurD,'k.'); axis([0.5 1 0 200]); subplot(122); hist(cat_heur.*cat_heurD,2000); drawnow;
    
    % walk back through all the segs to get the waveforms and timestamps
    tsTmp = [];
    indTmp = [];
    figure(1);  subplot(1,2,1); hold off; plot(template_vector,'b','LineWidth',3,'Color',[0.5 0.8 1]); hold on; 
        
    for k=1:numel(segArray)    
        
        currSeg = segArray(k);
        
        % take all suprathreshold events
                %         matchedInds = find(heuristics.seg(k).c>=userThresh(1));
                %         tmp_metric = heuristics.seg(k).c.*heuristics.seg(k).d;
                %         matchedInds = find(heuristics.seg(k).c>=userThresh(1) & heuristics.seg(k).d>0);
        matchedInds = find(heuristics.seg(k).c >= userThresh & heuristics.seg(k).d <= userThreshY);
        
        if updateTnsFlag==1
            idList.seg(currSeg).shank(shank_num).id(matchedInds) = currUnit;
        end
        
        % prep data for the tsd format ***NEED TO FIX THIS TO WORK WITH A SINGLE SHANK***
        if size(featStruct.seg(currSeg).shank(shank_num).ts,1) < size(featStruct.seg(currSeg).shank(shank_num).ts,2)
            tsTmp   = [tsTmp featStruct.seg(currSeg).shank(shank_num).ts(matchedInds)+ ((currSeg-1).*featStruct.chunk.*1000)];
        else
            tsTmp   = [tsTmp featStruct.seg(currSeg).shank(shank_num).ts(matchedInds)'+ ((currSeg-1).*featStruct.chunk.*1000)];
        end

        if size(featStruct.seg(currSeg).shank(shank_num).inds,1) < size(featStruct.seg(currSeg).shank(shank_num).inds,2)
            indTmp  = [indTmp featStruct.seg(currSeg).shank(shank_num).inds(matchedInds)];
        else
            indTmp  = [indTmp featStruct.seg(currSeg).shank(shank_num).inds(matchedInds)'];
        end
        
        % update the tns file?
        
        % store the mean extracted waveforms and the input template
        for m=1:numel(matchedInds)
            if m==1
                target_values = single(sessionStruct.seg(currSeg).shank(shank_num).wfs(matchedInds(m)).values);                
            else
                target_values = target_values + single(sessionStruct.seg(currSeg).shank(shank_num).wfs(matchedInds(m)).values);
            end
        end
        
        % display a quick sanity check
        if chan_per_shank >= 3
            target_vector(k,:) = [target_values(indices(1),:) target_values(indices(2),:) target_values(indices(3),:)]  ./ numel(matchedInds) ;
        else
            target_vector(k,:) = [target_values(indices(1),:)]  ./ numel(matchedInds) ;
        end

        figure(2);  subplot(1,2,1); plot(target_vector(k,:),'k'); drawnow;
        

    end
    pause();
    disp(['::: Storing data in .tsd format to the shank datastructure']);

    % store data in a tsd format
    shank.unit(j).ts                    = tsTmp;
    shank.unit(j).inds                  = indTmp;
    shank.unit(j).template_uid          = currUnit;
    if chan_per_shank >= 3   
        shank.unit(j).template_sites        = indices(1:3);
    else
        shank.unit(j).template_sites        = indices(1);
    end
    shank.unit(j).sourceData            = tns_Path;
    shank.unit(j).waveform.template     = template_vector;
    shank.unit(j).waveform.targets      = target_vector;
    shank.unit(j).segs.segArray         = segArray;
    shank.unit(j).segs.seed_segs        = seed_segs;

end

%% WRITE DATA TO DISK
disp('4) Finishing extension of sorted data across all segments in array');
disp('________________________________________________________________________________________________________________________');
disp(' ');
disp(['     save ' file_Path '_shank' num2str(shank_num) 'X_tsd shank']); 
eval(['     save ' file_Path '_shank' num2str(shank_num) 'X_tsd shank']); 

if updateTnsFlag==1
    disp(' ');
    disp(['     save ' tns_Path ' idList']); 
    save(tns_Path,'idList');
end

disp('________________________________________________________________________________________________________________________');


%% DEPRECATED CODE

%             test_vector = [tmp(indices(1),:) tmp(indices(2),:) tmp(indices(3),:) tmp(indices(4),:)];            
%             heuristics.seg(k).d(1,m) = dot(template_vector,test_vector);
%             heuristics.seg(k).d(1,m) = sqrt( mean( (template_vector-test_vector).^2 ) );
%             heuristics.seg(k).d(1,m) = dot(template_amps , test_amps);

%         target_values = target_values';
%         target_vector(k,:) = target_values(1:numel(sessionStruct.seg(currSeg).shank(shank_num).wfs(1).values)) ./ numel(matchedInds);     
%         target_vector(k,:) = [target_values(indices(1),:) target_values(indices(2),:) target_values(indices(3),:) target_values(indices(4),:)]  ./ numel(matchedInds) ;

