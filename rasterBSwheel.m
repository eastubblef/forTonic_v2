function [data] = rasterBSwheel(varargin)

%% Yttri's mfile 7/15/15 - So far no mods (BS 7/20/15)
% onUnitAlignSpike_laser & oneUnitAlignSpike_wheel call this function for rasterizing

unitnum = varargin{2};
% unitnum = num2str(k);
unitname = strcat(' u#', unitnum);

check=ishold;
% disp('just checking - input must be row vectors!);
disp('processing')
disp(unitname);
disp(unitnum);


data = varargin{1}; %'data' input always enters this function as a cell array per ea. unit 
Color='k';
% if nargin == 3 
% Color = varargin{2};
% newfig = varargin{3};
% elseif nargin ==2
%     if isstr(varargin{2})
%         Color = varargin{2};
%     else 
%         newfig = varargin{2};
%     end
% end
   
if isstruct(data)
    disp('data!!!')
    tss=nan(length(data.session.unit), 100000);  %probably too much, equivalent to 1 hour @ 25hz
    for thisunit = 1:length(data.session.unit)
        thists=data.session.unit(thisunit).ts;
        tss(thisunit,1:length(thists))=thists;
    end
    [m,n]=find(isfinite(tss));
    data=tss(:,1:max(n));
end
    

if exist('newfig')
    figure(newfig)
end

sized=size(data);
data=flipud(data);
tic
for rows=1:sized(1)
    plot([data(rows,:);data(rows,:)],[ones(1,sized(2))*.9;zeros(1,sized(2))]+rows-1, Color);
    hold on
end

set(gca, 'YTick', []);
hold off
if check
    hold on
end
toc
outs = data;
