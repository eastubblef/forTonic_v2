
% This function can pull out average waveform from clustered unit on the 3 channels with the largest waveforms
% Updated 6.24.16 - works for as many units as are available on the seed segment - for 1 u at a time, plotted more professionally see version 2.0
% BS wrote, in part by fragmenting Josh's TNC_ExtendManualSort code 1.20.16

%% INPUTS:
% segArray:  not needed; for template-matching, use TNC_ExtendManualSortBS2.m
% seed_segs: array of segments that should be used to generate the waveforms
% shank_num: what shank to perform sorting on
% file_Path: basename of the files used ('_ft_tns.mat' '_ss.mat' '_ft.mat' are required)
% updateTnsFlag: 1=overwrite tns file with new sorted indices
% 
% EXAMPLE:
% [shank] = TNC_ExtendManualSort([1:2],[1],5,'DA_F06_20140815-i16',1);
% 
% OUTPUTS:
% Can overwrite the spike ids file (*_tns) used for sorting in GUI

function [shank] = avgWaveforms(seed_segs,shank_num,file_Path,updateTnsFlag)

% % valid seg range for 32 ch. probes (4 shanks, 8 channels each):
% segArray = [1:8];
% seed_segs = [5];
% shank_num = 1;

% paths to sessStruct and ft_tns files:
% tns_Path              = [file_Path '_ft_tns.mat'];
% sessionStruct_Path    = [file_Path '_ss.mat'];
% featData_Path         = [file_Path '_ft.mat'];

updateTnsFlag = 0;    %do not overwrite
% segArray = [3];     
seed_segs = [26];
shank_num = 1;

% file_PathApply = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree'                                 %will have "chunk" in the name
% file_PathSeed =  '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/sorted/1maxSeg/151104ipsi/1001'  %will have "tag" in the name
% file_Path = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/'                                     %Put them all in same folder
file_Path = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better/safeKeepingLaserLast/';

% tns_Path              = [file_Path '160505chunks_1001.ns5_ft_tns.mat'];    %already sorted - use these units' waveforms for template 
% sessionStruct_Path    = [file_Path '160505chunks_1001.ns5_ss.mat'];        %not sorted; for matching - no longer true (1.20.16)
% featData_Path         = [file_Path '160505chunks_1001.ns5_ft.mat'];        %not sorted; for matching - no longer true (1.20.16)

tns_Path              = [file_Path '151105chunks_1001_3.ns5_ft_tns.mat'];  %already sorted - use these units' waveforms for template 
sessionStruct_Path    = [file_Path '151105chunks_1001_3.ns5_ss.mat'];      %not sorted - for matching
featData_Path         = [file_Path '151105chunks_1001_3.ns5_ft.mat'];      %not sorted - for matching

S = load(sessionStruct_Path);
sessionStruct=S.sessionStruct;
clear S

S = load(tns_Path);
idList=S.idList;
clear S

S = load(featData_Path);
featStruct=S.featStruct;
clear S

disp(' ');
disp(['1) Extracting data from shank ' num2str(shank_num) ' using waveforms defined from segments [' num2str(seed_segs) '] <<< ' file_Path]);

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
        
        ids_for_curr_unit = find(idList.seg(currSeg).shank(shank_num).id==currUnit);                                                        %ids for curr unit are found in the tns file (seed file) at these unit's indexed locations
        
        for m=1:numel(ids_for_curr_unit)
            if m==1 & p==1
                template_values = single(sessionStruct.seg(currSeg).shank(shank_num).wfs(ids_for_curr_unit(m)).values);                     %unit 1's template waveform (template_values) is is found in the ss file's wfs indexed to the nth location from the idList          
            else
                template_values = template_values + single(sessionStruct.seg(currSeg).shank(shank_num).wfs(ids_for_curr_unit(m)).values);   %but once iterative >1(ea. time u. 1 appears in the idList), the template_values (waveforms) are added on, in succession
            end
        end
        
        total_events = total_events+numel(ids_for_curr_unit);                  %does numel(id) correpond to #dots in cluster? (see l. 134 - yes)

    end
                                                                                
    % find the 3 channels with the largest amplitudes                           %on the shank of interest for segment to be tested
    template_amps = min(template_values,[],2);                                  %gives a vector from the col. of ea. of 8ch's mins
    [values,indices] = sort( template_amps , 'ascend' );                        %[x,y] = largest to smallest order of template amps for plotting
    chan_per_shank = numel(indices);                                            %each index here is a channel; ea. row of template_values is a channel
    
    % NEED TO DEAL WITH VARIABLE SHANK SIZES HERE
    if chan_per_shank >= 3
        template_vector = [template_values(indices(1),:) template_values(indices(2),:) template_values(indices(3),:)] ./ total_events;
        template_amps = [min(template_values(indices(1),:)) min(template_values(indices(2),:)) min(template_values(indices(3),:))] ./ total_events;        
    else
        template_vector = template_values(indices(1),:) ./ total_events;        %template vector is the largest avg. waveform on 1/8ch. 
        template_amps = [min(template_values(indices(1),:))] ./ total_events;   %this shank's 3 ch. with largest avg. waveforms of this unit     
    end

    % display 3 waveforms for 3 channels
    figure(j); subplot(1,2,1); plot(template_vector,'b','LineWidth',3,'Color',[0.5 0.8 1]); title(['unit ' num2str(currUnit) ' seg ' num2str(seed_segs) ' shank ' num2str(shank_num)]);

    % up to l.142
    
end
