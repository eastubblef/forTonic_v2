
%% Smoothing factor
Sigma = 50;
t = 1:1:Sigma.*30;
Gaussian = ( 1./( sqrt(2.*pi.*Sigma.*Sigma) ) ) .* exp( -(t-Sigma.*15).^2 ./ (2.*Sigma).^2 );
integral = trapz(Gaussian);
Gaussian = Gaussian./integral;
[mapName] = TNC_CreateRBColormap(8,'mbr');

%% Simple joystick analysis + finding reaches
pos.x = sgolayfilt(decimate(Xpos,25),3,501);
pos.y = sgolayfilt(decimate(Ypos,25),3,501);
pos.v = [0 sqrt(diff(pos.x).^2 + diff(pos.y).^2)];
pos.v_win = sgolayfilt( pos.v , 3 , 2001 );
pos.x = pos.x-median(pos.x);
pos.y = pos.y - median(pos.y);
supra_thresh    = find(pos.v_win>0.67);
jumps           = find([1 diff(supra_thresh)]>1);
valids          = [];
max_vel         = [];
direction       = [];
[sink]      = TNC_ExtTrigWins(pos.v_win,supra_thresh(jumps),[1500,3000]);
% [sink_6]    = TNC_ExtTrigWins(tmpSmooth,supra_thresh(jumps),[1500,3000]);
[sink_x]    = TNC_ExtTrigWins(pos.x,supra_thresh(jumps),[1500,3000]);
[sink_y]    = TNC_ExtTrigWins(pos.y,supra_thresh(jumps),[1500,3000]);


% skip the first 20 reaches
offset = 20;

%Find valid reaches
for i=offset:numel(jumps)
    if (sqrt(sink_x.wins(i,2000).^2 + sink_y.wins(i,2000).^2) > sqrt(sink_x.wins(i,1200).^2 + sink_y.wins(i,1200).^2)) && max(sink.wins(i,1:1200))<2
        valids = [valids i];
        this_vel = max(sink.wins(i,1500:2500));
        max_vel = [max_vel this_vel];        
        if sink_y.wins(i,1800) > sink_y.wins(i,1200)
            direction = [direction 1];
            plot(sink_x.wins(i,1000:2000)-sink_x.wins(i,1200),sink_y.wins(i,1000:2000)-sink_y.wins(i,1200),'Color',[mapName(1,:)]); hold on;        
        else
            direction = [direction -1];
            plot(sink_x.wins(i,1000:2000)-sink_x.wins(i,1200),sink_y.wins(i,1000:2000)-sink_y.wins(i,1200),'Color',[mapName(11,:)]); hold on;        
        end
    else
        disp('skipped');
    end
end

%% Get PSTHs for recorded cells

    figure(10); clf;
    
for cell_id = 1:numel(S_clu.cviSpk_clu);

    tmp = double(viTime_spk(S_clu.cviSpk_clu{cell_id}))*1000/25e3;
    delta = zeros(1,ceil(max(tmp)));
    delta(round(tmp)) = 1;

    tmpSmooth = conv(delta,Gaussian,'same');
    
    subplot(1,3,[1 2]);
    plot(tmpSmooth + (ones(1,numel(tmpSmooth)).*cell_id.*0.1)); hold on;
    drawnow;
    
    subplot(1,3,3);
    
    times = supra_thresh(jumps(valids));
    valid_valid = find(times>7e5 & times<ceil(max(tmp)));
    [sink_tmp]    = TNC_ExtTrigWins(tmpSmooth,times(valid_valid),[1000,2500]);
    plot((sink_tmp.avg - mean(sink_tmp.avg(1:500))) + (ones(1,numel(sink_tmp.avg)).*cell_id.*0.01)); hold on;

end

%% Examine neural correlates



for cell_id = 14 %2:14

tmp = double(viTime_spk(S_clu.cviSpk_clu{cell_id}))*1000/25e3;
delta = zeros(1,ceil(max(tmp)));
delta(round(tmp)) = 1;

tmpSmooth = conv(delta,Gaussian,'same');

% [sink]      = TNC_ExtTrigWins(pos.v_win,supra_thresh(jumps),[1500,3000]);
[sink_6]    = TNC_ExtTrigWins(tmpSmooth,supra_thresh(jumps),[1500,3000]);
% [sink_x]    = TNC_ExtTrigWins(pos.x,supra_thresh(jumps),[1500,3000]);
% [sink_y]    = TNC_ExtTrigWins(pos.y,supra_thresh(jumps),[1500,3000]);

figure(1); clf;
valids = [];
max_vel = [];
max_rate = [];
direction = [];
traj_mat = [];

for i=20:285
    if (sqrt(sink_x.wins(i,2000).^2 + sink_y.wins(i,2000).^2) > sqrt(sink_x.wins(i,1200).^2 + sink_y.wins(i,1200).^2)) && max(sink.wins(i,1:1200))<1
        valids = [valids i];
        this_vel = max(sink.wins(i,1500:2500));
        max_vel = [max_vel this_vel];
        this_rate = max(sink_6.wins(i,1750:2400).*1000);
        max_rate = [max_rate this_rate];
        
        if sink_y.wins(i,1800) > sink_y.wins(i,1200)
            direction = [direction 1];
            plot(sink_x.wins(i,1000:2000)-sink_x.wins(i,1200),sink_y.wins(i,1000:2000)-sink_y.wins(i,1200),'Color',[mapName(1,:)]); hold on;        
        else
            direction = [direction -1];
            plot(sink_x.wins(i,1000:2000)-sink_x.wins(i,1200),sink_y.wins(i,1000:2000)-sink_y.wins(i,1200),'Color',[mapName(11,:)]); hold on;        
        end
        
        traj_mat = [traj_mat ; sink_y.wins(i,1000:2000)-sink_y.wins(i,1200)];
        
    else
        disp('skipped');
    end
end

[Vs,Is] = sort(max_vel .* direction);

[mapName] = TNC_CreateRBColormap(8,'mbr');

figure(2); clf;
subplot(121); 
imagesc(sink.wins(valids(Is),:),[0 10]);
colormap(mapName);
subplot(122); 
imagesc(sink_6.wins(valids(Is),:).*1000,[0 40]);
colormap(mapName);

figure(3); clf;
plot((max_vel).*direction,max_rate,'ko');
xlabel('log(velocity)');
ylabel('max(Hz)');

[rho,pval] = corr([max_vel',max_rate']);
[rhoU,pval] = corr([max_vel(find(direction>0))',max_rate(find(direction>0))']);
[rhoD,pval] = corr([max_vel(find(direction<0))',max_rate(find(direction<0))']);

rhoU(2,1)-rhoD(2,1)

% plot PSTHs for each direction
figure(100); 
subplot(7,2,cell_id);
hold off;
shadedErrorBar(-1500:3000,mean(sink_6.wins(valids(vv==1),:),1),std(sink_6.wins(valids(vv==1),:),[],1)./sqrt(numel(valids(vv==1))),{'Color',[mapName(1,:)]});
hold on;
shadedErrorBar(-1500:3000,mean(sink_6.wins(valids(vv==2),:),1),std(sink_6.wins(valids(vv==2),:),[],1)./sqrt(numel(valids(vv==2))),{'Color',[mapName(11,:)]});
plot([0 0], [-0.001 0.02] , 'k--'); axis tight; box off;
drawnow;
end

%%
figure(5);

[binnedData] = TNC_BinAndMean(max_vel,max_rate,5);
plot([ binnedData.bins.center],[ binnedData.bins.avg],'ko-');

