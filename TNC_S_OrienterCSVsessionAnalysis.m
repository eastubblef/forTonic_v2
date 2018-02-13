
csv_data = dlmread('mVgatfive_2017_01_12_165817pXY.csv',',',2,0);
[data] = openNSx('/170112002.ns4');
S = load('170112behave_pop.mat');
PopData = S.PopData;
clear S;

%%
Sigma = 24;
t = 1:1:Sigma.*30;
Gaussian = ( 1./( sqrt(2.*pi.*Sigma.*Sigma) ) ) .* exp( -(t-Sigma.*15).^2 ./ (2.*Sigma).^2 );
integral = trapz(Gaussian);
Gaussian = Gaussian./integral;
[mapName] = TNC_CreateRBColormap(8,'mbr');

%%

figure(3); hold off;

x_off = 107700;
x_vals = (csv_data(:,1)-csv_data(1,1))+x_off;
d_x_vals = [0 diff(x_vals')];
smth_wheel_v = conv( abs(csv_data(:,2)) , [zeros(1,4) ones(1,5) 0] , 'same' )';



plot(x_vals , csv_data(:,2) ,'LineWidth',2);

hold on; 

plot(x_vals , smth_wheel_v ,'LineWidth',2);

plot(x_vals, d_x_vals);

% 
% plot(x_vals ,[0 diff(-csv_data(:,3)')]./2,'LineWidth',2); 
% plot(x_vals ,csv_data(:,3),'LineWidth',2);

% plot(sgolayfilt((decimate(double(data.Data(3,:)),10)./25),3,101));

% plot(decimate(double(data.Data(1,:)),10)./1e2);


%% find intra movement periods

threshold = 5;

intra_mvmt = find(smth_wheel_v>threshold & d_x_vals<30);

figure(4); hold off;

plot(x_vals , csv_data(:,2) ,'LineWidth',1);
hold on;
plot(x_vals(intra_mvmt) , csv_data(intra_mvmt,2) ,'.');

%% separate out into valid movements

mvmt_win = [10,25];
starts = find([0 diff(intra_mvmt)]>10);

starts_valid = find(intra_mvmt(starts)>mvmt_win(1) & intra_mvmt(starts)<(numel(csv_data(:,2))-mvmt_win(2)));

starts_clean = starts(starts_valid);
plot(x_vals(intra_mvmt(starts_clean)) , csv_data(intra_mvmt(starts_clean),2) ,'o');

[wheel_moves] = TNC_ExtTrigWins(csv_data(:,2)',intra_mvmt(starts_clean),mvmt_win);

for i=1:size(wheel_moves.wins,1)
    magnitude(i) = trapz(wheel_moves.wins(i,mvmt_win(1):20));
end

[vals, inds] = sort(magnitude);
figure(5); imagesc(wheel_moves.wins(inds,:),[-50 50]); colormap(mapName); %x is the movement vectors sorted by magnitude 
    xlabel('position'); ylabel('movement vector');

export_mvmt_data_for_psths.times        = x_vals(intra_mvmt(starts_clean))';
export_mvmt_data_for_psths.magnitude    = magnitude;
export_mvmt_data_for_psths.mag_sort_i   = inds;


%% examine some psths associated with movements

numUnits = numel( PopData.session(1).unit );
figure(10); clf;
for j=1:numUnits
    
    tmp = PopData.session(1).unit(j).ts;
    delta = zeros(1,ceil(max(tmp)));
    delta(round(tmp)) = 1;

    tmpSmooth = conv(delta,Gaussian,'same');
    
    figure(10);
    subplot(1,3,[1 2]);
    plot(tmpSmooth + j*0.1); hold on;
    
    subplot(1,3,3);
    times = export_mvmt_data_for_psths.times;
    valid_valid = find(times>9.2e4 & times<2.07e6);
    [vals, inds] = sort(magnitude(valid_valid));
    [sink_tmp]    = TNC_ExtTrigWins(tmpSmooth,times(valid_valid),[750,1000]);
    
    plot((sink_tmp.avg - mean(sink_tmp.avg(1:500))) + (ones(1,numel(sink_tmp.avg)).*j.*0.01)); hold on;
    drawnow;
    
    if j==5
       
        figure(11);
        imagesc(sink_tmp.wins(inds,:));
        
        figure(12); hold off;
        plot(tmpSmooth.*250);
        hold on;
        plot(x_vals , csv_data(:,2) ,'LineWidth',2);        
        
    end

end


%%