                     
%% LOAD excel file 
newPath = '/Users/stubblefielde/Desktop';
cd(newPath);
    
% [filenamestr, path] = uigetfile('*.xls*','select the .xls file', newPath)              
filenamestr = '32ChMapFlipped4matlab';
[num, txt] = xlsread(filenamestr);

c9 = find(num == C9);

