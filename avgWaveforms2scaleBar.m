function [hX, hY, udataX, udataY] = avgWaveforms2scaleBar(scalebarOn, newUnitX, newUnitY, unitNum,laserOn,seed_segs,shank_num,file_Path,updateTnsFlag)

% This function can pull out the average waveform from a clustered unit on the channel w/ the largest waveform
% Updated 6.30.16 - now plots a scalebar (default = 10 uV, 1 ms in length)
% Updated 6.29.16 - funky axes units from the ss file (y/10 = uV; x*10 = uSec) - use new input to make y = uV; x = ms
% Updated 6.24.16 - works for as many units as are available on the seed segment 
% EAS wrote, in part by fragmenting Josh's TNC_ExtendManualSort code 1.20.16

%% INPUTS:
% scalebar          : 1 = plot a scalebar; default is 1 ms long for x-axis; 10uV long for y-axis
% newUnitX, newUnitY: the desired x and y axes units for the waveform plot
% segArray          : not needed; for template-matching, use TNC_ExtendManualSortBS2.m
% seed_segs         : array of segments that should be used to generate the waveforms
% shank_num         : what shank to perform sorting on
% file_Path         : basename of the files used ('_ft_tns.mat' '_ss.mat' '_ft.mat' are required)
% updateTnsFlag     : 1 = overwrite tns file with new sorted indices; NOT recommended
% 
% OUTPUTS:
% Can overwrite the spike ids file (*_tns) used for sorting in GUI

% % valid seg range for 32 ch. probes (4 shanks, 8 channels each):
% segArray = [1:8];
% seed_segs = [5];
% shank_num = 1;
% tns_Path              = [file_Path '_ft_tns.mat'];
% sessionStruct_Path    = [file_Path '_ss.mat'];
% featData_Path         = [file_Path '_ft.mat'];

scalebarOn = 1;          %draw a scalebar on the plot; default x-length = 1 ms
newUnitX = 100;          %desired x axis units - 100 will plot x in ms
newUnitY = 10;           %desired y axis units - 10 will plot y in uV
updateTnsFlag = 0;       %do NOT overwrite the tns file
% segArray = [3];     

% seed_segs = [26];
seed_segs = [23];

% seed_segs = [17]
shank_num = 1;
unitNum = 1;             %cluster number corresponding to individual unit
laserOn = 0;             %if not a laser segment, use 0 (for plotting)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% file_Path = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better/safeKeepingLaserLast/';
file_Path = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better/behaveChunks/';
% file_Path = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151104/1001/laser/';

% Behavior segment
tns_Path              = [file_Path '151105chunks_all.ns5_ft_tns.mat'];    %already sorted - use these units' waveforms for template 
sessionStruct_Path    = [file_Path '151105chunks_all.ns5_ss.mat'];        %not sorted; for matching - no longer true (1.20.16)
featData_Path         = [file_Path '151105chunks_all.ns5_ft.mat'];        %not sorted; for matching - no longer true (1.20.16)

% Laser segment
% tns_Path              = [file_Path '151104chunks_1001.ns5_ft_tns.mat'];  %already sorted - use these units' waveforms for template 
% sessionStruct_Path    = [file_Path '151104chunks_1001.ns5_ss.mat'];      %not sorted - for matching
% featData_Path         = [file_Path '151104chunks_1001.ns5_ft.mat'];      %not sorted - for matching

% tns_Path              = [file_Path '151105chunks_1001_3.ns5_ft_tns.mat'];  %already sorted - use these units' waveforms for template 
% sessionStruct_Path    = [file_Path '151105chunks_1001_3.ns5_ss.mat'];      %not sorted - for matching
% featData_Path         = [file_Path '151105chunks_1001_3.ns5_ft.mat'];      %not sorted - for matching

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

% start by finding the number of units in the seed_segs
unit_list_tmp=[];
for i=1:numel(seed_segs)
    unit_list_tmp = [unit_list_tmp unique(idList.seg(seed_segs(1)).shank(shank_num).id)];
end
unit_list = unique(unit_list_tmp(find(unit_list_tmp>0)));                                   %retains u#s from clusters in TNC_SS_GUI

disp(['2) Found ' num2str(numel(unit_list)) ' templates <<< [' num2str(unit_list') ']']);

%% Calculate the inner products for all detected events within the segment range

disp(['3) Extending manually sorted units across all segments...']);

% for j=1:numel(unit_list)
%     figure(1); clf;
%     figure(j);
figure; hold on;
%     currUnit = unit_list(j);
    currUnit = unitNum
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
                                                                                
    % find the channel with the largest amplitude                            %on the shank of interest for segment to be tested
    template_amps = min(template_values,[],2);                                  %gives a vector from the col. of ea. of 8ch's mins
    [values,indices] = sort( template_amps , 'ascend' );                        %[x,y] = largest to smallest order of template amps for plotting
    chan_per_shank = numel(indices);                                            %each index here is a channel; ea. row of template_values is a channel
    
    %New for one unit plotting of the largest waveform; the first index has the largest amplitude one
    template_vector = template_values(indices(1),:) ./ total_events;        %template vector is the largest avg. waveform on 1/8ch. 
    template_amps = [min(template_values(indices(1),:))] ./ total_events;   %this shank's 3 ch. with largest avg. waveforms of this unit     
    
    ylabel('/10 = uV'); ylim([-1200 500]);                                  %for DataAspectRatio = [1 1]; needed to feed into the scalebar function
    xlabel('*10 = uSec'); xlim([-1200 500]);

    if laserOn == 1
        plot(template_vector,'b','LineWidth',3,'Color',[0.5 0.8 1]); title(['unit ' num2str(currUnit) ' LaserSeg ' num2str(seed_segs) ' shank ' num2str(shank_num)]);
    else
        plot(template_vector,'k','LineWidth',3,'Color',[0 0 0]); title(['unit ' num2str(currUnit) ' BehaveSeg ' num2str(seed_segs) ' shank ' num2str(shank_num)]);
    end

 %since axes are already off by: y/10 = uV; x*10 = uSec coming out of the ss file, make it right:
    if newUnitX == 100
       newXTicks = {-120 -100 -80 -60 -40 -20 0 20 40 60};      %aspectRatio [1 1] compared to yAxis (in uV), but this scales to /10 = ms
       newXTicks = {-12 -10 -8 -6 -4 -2 0 2 4 6};               %this scales to ms
    end
    if newUnitY == 10
       newYTicks = {-120 -100 -80 -60 -40 -20 0 20 40 60};      %scales to microV
    end    
    
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', newXTicks, 'YTickLabelMode', 'manual', 'YTickLabel', newYTicks);
    ylabel('uV'); 
    xlabel('ms'); 
 
%% define scale bar parameters:

    if scalebarOn == 1
       n = newUnitX;    %plot x axis length as 1 ms
       directions = {'northwest','northeast','southeast','southwest'};
       [hX, udataX] = scalebarX('Location', 'northeast','Unit', 'ms', 'scalelength', n);  
    end
    if scalebarOn == 1
       n = newUnitY *10; %plot y axis length as 10 uV
       directions = {'northwest','northeast','southeast','southwest'};
       [hY, udataY] = scalebarY('Location', 'northeast', 'Unit', 'uV', 'scalelength', n);  

    end

    if scalebarOn == 0
        return
    end    
    
end
