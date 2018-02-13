function [PopData] = TNC_ConvertNEVtoPopDataBS(source, sessNum, fpath)

%% This mfile is for gathering sorted NEURONAL spikes for eventual alignment with behavoral events (from TNC_MoverBehaviorExtract, for example).
% INPUTS
% sessNum is the sorted spike unit (cluster 1, 2, 3, etc.) from .ns5 file = raw! (NEV = threshold crossings)
% source is the dir where tsd.mat file is located; tsd files are saved from running: 
%  -TNC_SS_GUI for spike sorting (will have "chunk" in the fname)
%  -TNC_ExtendManualSort (BS2 ver.) for wfs template matching (will have "shank" in the fname)

% OUTPUTS 
% [PopData] is a structure in which the col. "TS" is for each neuron (sorted unit), organized by shank #
% (See Yttri's code alignSpikeBS.m for next step - ea. neuron is diff. row; col is TS of ea. spike for that neuron; take for every unit)

% Currently, this mfile is called from MoverScript.m

% BS Updated last 7.13.15 (for laser/spike - alignment)
% 1.14.15 (for other behavioral events) - say, any movement onset

%%
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/laserChunks/untagged/'
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks/';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/170112/2ndPass/behaveSegs/';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/test/';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/161004/';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16/161012/';
% fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170112/';
fpath = '/Volumes/My Passport for Mac/Vgatthree_updated/151105/2ndpass/better/behaveChunks/';

cd(fpath); 
source = fpath;
sessNum = 1;

ucnt = 0;
figure(1); clf;

% [source] = strcat(fpath, '150529_3004.ns5_shank1_tsd.mat');               %nope, need the dir, not the file
% [fname, path] = uigetfile('*tsd.mat*','select the tsd file', fpath)       %BS added

d = dir([source '*_NEV.mat']); %CAREFUL! THIS WILL LOAD ALL tsd.mat FILES IN THE DIR

if numel(d)>0

    for i=1:numel(d)

        disp(['Loading data from >>> ' d(i).name]);
        s = load(d(i).name);

        numUnits = numel(s.shank.unit);

        for k=1:numUnits
            
            ucnt = ucnt+1;
            
            nonZeros = find(diff(s.shank.unit(k).ts)>1);            
            
            PopData.session(sessNum).unit(ucnt).inds    = s.shank.unit(k).inds(nonZeros);            
            PopData.session(sessNum).unit(ucnt).ts      = round(PopData.session(sessNum).unit(ucnt).inds./30);
%             PopData.session(sessNum).unit(ucnt).ts      = s.shank.unit(k).ts(nonZeros);
            PopData.session(sessNum).unit(ucnt).sh      = str2num(d(i).name(strfind(d(i).name,'shank')+5));
            
            [isi] = TNC_QuantISI(PopData.session(sessNum).unit(ucnt).ts);
            PopData.session(sessNum).unit(ucnt).isi     = isi;

        end

    end

    disp(['Session number ' num2str(sessNum) ' was written to the PopData structure and contains ' num2str(ucnt) ' single units.']);    

    dims = ceil(sqrt(ucnt));
        
    for m=1:ucnt            
        figure(1);
%         subplot(dims,dims,m);
        subplot(1,ucnt,m); hold off;
        PopData.isiDist(m,:) = PopData.session(sessNum).unit(m).isi.hist.logCount;
        semilogx(PopData.session(sessNum).unit(m).isi.hist.logTimes,PopData.session(sessNum).unit(m).isi.hist.logCount,'LineWidth',1,'Color',[(m./ucnt) 0.5 1-(m./ucnt)]); hold on;
        semilogx([1 1],[0 max(PopData.session(sessNum).unit(m).isi.hist.logCount)],'k--');
        title(['sh' num2str(PopData.session(sessNum).unit(m).sh) 'u' num2str(m) '   |   ' num2str(round(1000./PopData.session(sessNum).unit(m).isi.stats.mean)) ' Hz ']);
%         grid on; axis([1e-2 1e4 0 max(PopData.session(sessNum).unit(m).isi.hist.logCount)]); axis off;
    end        
    
%     [ids] = kmeans(PopData.isiDist,3);
    [ids] = (PopData.isiDist);

    [ys,is] = sort(ids);
    PopData.isiClass.y = ys;
    PopData.isiClass.i = is;
    
%     [cMap] = TNC_CreateRBColormap(12,'mbr');
%     figure(2); 
%     subplot(211); imagesc(PopData.isiDist(is,:)); colormap(cMap);
%     subplot(212); imagesc(corr(PopData.isiDist(is,:)'),[0 1]); colormap(cMap);
    
else

    disp('Could not find any sorted files with that name.');

end

