%insert starting at line 94   

for p=1:numel(seed_segs)
        currSeg = seed_segs(p);
        
        ids_for_curr_unit = find(idList.seg(currSeg).shank(shank_num).id==currUnit);                                                        
        %ids for curr unit are found in the tns file (seed file) at these unit's indexed locations
        
        for m=1:numel(ids_for_curr_unit)
            if m==1 & p==1
                template_valuesAll = single(sessionStruct.seg(currSeg).shank(shank_num).wfs(ids_for_curr_unit(m)).values);                     
                %unit 1's template waveform (template_values) is is found in the ss file's wfs indexed to the nth location from the idList          
            else
                %template_values = template_values + single(sessionStruct.seg(currSeg).shank(shank_num).wfs(ids_for_curr_unit(m)).values);   
               template_valuesAll = [template_valuesAll;single(sessionStruct.seg(currSeg).shank(shank_num).wfs(ids_for_curr_unit(m)).values)];   
               
                %but once iterative >1(ea. time u. 1 appears in the idList), the template_values (waveforms) are added on, in succession
            end
        end
        
        total_events = total_events+numel(ids_for_curr_unit);                  %does numel(id) correpond to #dots in cluster? (see l. 134 - yes)

    end
    
    siteNum=size(single(sessionStruct.seg(currSeg).shank(shank_num).wfs(ids_for_curr_unit(m)).values),1);

    template_values=[];
    for sites = 1:siteNum  
        whichSites=[sites:siteNum:size(template_valuesAll,1)];
        template_values =  [template_values; nanmean(template_valuesAll(whichSites,:))] ; 
    end
    
      %% find the 1 channel with the largest amplitude                           %on the shank of interest for segment to be tested
   % ...