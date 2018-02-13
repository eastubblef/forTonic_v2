%% Purpose of this mfile is to
% - rename all sorted .ns5 files from *_tsd to tsd.mat
% - move these files into their own folders based on depth (ex: 150603_contra002.ext goes into folder called "002" within 150603contra folder)

%% UPDATED 11.26.15:
% - Files are moved into folder: mouse_name"Sorted"

%%
% path = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Gad2_tag/';
% path = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgat/';
path = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/sorted';

cd(path);
% ext = '150602contra'                                                      %For contra files
% ext = '150228ipsi'
ext = '151105ipsi'

%newDir = '7007'
newDir = '1001'
mkdir(ext, newDir);

%pathString = strcat(path, ext);
pathString = strcat(ext);

cd(pathString);
% nameCat = strcat('*_contra', newDir)                                      %For contra files with "contra" in the name
% nameCat = strcat(ext, '_', newDir);                                       %For ipsi files with "ipsi" in the name
nameCat = strcat(ext(1:end-4), '_', newDir);                                %For ipsi files without "ipsi" in the name

namePlus = strcat(nameCat, '*');

dir namePlus

% d = dir(fullfile(pathString, '*_contra4005*')); %contra files 
d = dir(fullfile(pathString, namePlus));          %works for both files! 

if numel(d) > 0
    for i = 1:numel(d)
        
        disp(['loading data from >>> ' d(i).name]);
        %         s = load(d(i).name); no need to load anything
        %         movefile('filename', newDir);
        movefile(d(i).name, newDir);
    end
else disp(['file not found']);
end


%% Change _tsd files to tsd.mat

cd(newDir)
newPathString = strcat(pathString, '/', newDir);
tsdFiles = dir(fullfile(newPathString, '*_tsd'));

newPathString2 = strcat(newPathString, '/');
% Rename in a LOOP - actually will have to "rename" since "save" renders the file unusable
if numel(tsdFiles.name) > 0
    for n = 1:numel(tsdFiles)
        disp(['converting tsd data to tsd.mat for file: ' tsdFiles(n).name]);

        oldname = [newPathString2 tsdFiles(n).name];
        newname = strcat(newPathString2, tsdFiles(n).name, '.mat');
%         save(newname);
%         delete(oldname); %could be dangerous
        movefile(oldname, newname);
    end 
else disp(['failed to convert tsd data']);
end
disp(['files are ready for further processing']);

