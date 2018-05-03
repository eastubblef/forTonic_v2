%% This script gets rid of units that don't have movement-related activity
% Calls the existing Pop.mat file and generates a pop2.mat file instead (unwanted units have zeros for psth z-scored values

% fpath = '/Volumes/My Passport for Mac/171012/behaveNewUnits';
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/171012/newBehaveUnits';
fname = strcat(fpath, '/', '171012_pop.mat');     %should have 17 units; ipsi/contra = 36                     
targetName = '171012_pop2.mat';
load(fname);

excludeNeurons = [1 5 7 9 14 15 28 33 35]; %should be 9 units

Pop = PopData;
unit = PopData.session.unit;         

for i = 1:length(unit)
    for j = 1:length(excludeNeurons)
        if i == excludeNeurons(j) 
            newUnit(i).ts = 1;
            i = i+1
        else 
            newUnit(i).ts = unit(i).ts;
        end
    end
end

newPopData.session.unit = newUnit;

save([targetName], 'newPopData');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Next unit
% fpath = '/Volumes/My Passport for Mac/171012/behaveNewUnits';
fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170118/behaveNewUnits';
fname = strcat(fpath, '/', '170118_pop.mat');     %should have 17 units; ipsi/contra = 36                     
targetName = '170118_pop2.mat';
load(fname);

excludeNeurons = [3 4 6]; %should be 9 units

Pop = PopData;
unit = PopData.session.unit;         

for i = 1:length(unit)
    for j = 1:length(excludeNeurons)
        if i == excludeNeurons(j) 
            newUnit(i).ts = 1;
            i = i+1
        else 
            newUnit(i).ts = unit(i).ts;
        end
    end
end

newPopData.session.unit = newUnit;

save([targetName], 'newPopData');

