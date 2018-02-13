function OptoTagJune2015(spiketimes,tetrodefile,pulsedata, upper_limit)
%% Written by Andrew Wolf - given to BS 6/23/15
%Spiketimes, tetrodefile, and pulsedata should be paths for the appropriate files (see code below); upper_limit is to set upper limit for response(typical 0.01 seconds)
%Want waveforms for spikes w/in 10ms of light pulse and those outside of 10ms so they can be compared

%Load spike times as TS and convert to seconds - may need to use t_to_mat function first
load(spiketimes);
TS = TS./1e4;

%Settings to extract data from tetrode file (from Gidon) using Nlx2MatSpike
FieldSelection = [1 1 1 1 1];
ExtractHeader = 0;
ExtractMode = 1;
ModeArray = [];
Filename = tetrodefile;
[TimeStamps, ScNumbers, CellNumbers, Params, DataPoints] = Nlx2MatSpike(Filename, FieldSelection, ExtractHeader, ExtractMode, ModeArray);

%Convert DataPoints units to microvolts
hold = AD_Volt_Convert(tetrodefile);
DataPoints = hold.MicroVolts;

%Convert imported TimeStamps to seconds
TimeStamps = TimeStamps';
TimeStamps = TimeStamps./1e6;

%Load pulse data taskbase and get trial pulse times, convert to absolute pulse times, combine, and sort. Set upper limit (upper_limit) in function inputs.
load(pulsedata)
vector = cell2mat(taskbase.pulse_on);
pulse1 = taskbase.start_nlx(1:end-1) + vector(:,1);
pulse2 = taskbase.start_nlx(1:end-1) + vector(:,2);
pulse3 = taskbase.start_nlx(1:end-1) + vector(:,3);
pulse4 = taskbase.start_nlx(1:end-1) + vector(:,4);
pulse5 = taskbase.start_nlx(1:end-1) + vector(:,5);
pulse6 = taskbase.start_nlx(1:end-1) + vector(:,6);
pulse7 = taskbase.start_nlx(1:end-1) + vector(:,7);
pulse8 = taskbase.start_nlx(1:end-1) + vector(:,8);
pulses = cat(1,pulse1,pulse2,pulse3,pulse4,pulse5,pulse6,pulse7,pulse8);
sorted = sort(pulses);
start = 0;
sorted_upper = sorted + upper_limit;
sorted_upper = cat(1,start,sorted_upper);

%Spikelist variable preallocation - variable setup to make code run more efficiently (but not necessary with way code is currently written)
spikelist = zeros(20000,1);
spikelist2 = zeros(20000,1);

%Round spikelist variable to two decimal places
spikelist = round(spikelist./.01).*.01;
spikelist2 = round(spikelist2./.01).*.01;

%Need waveforms for spikes w/in 10ms (spikelist) of light pulse and those outside of 10ms (spikelist2)
for pulseindex = 2:length(sorted);
indexspike1 = find(and(TS(:,1) >= sorted(pulseindex - 1), TS(:,1) <= sorted_upper(pulseindex)));
spikelist = cat(1,spikelist,TS(indexspike1));
end

for pulseindex2 = 2:length(sorted);
indexspike2 = find(and(TS(:,1) < sorted(pulseindex2), TS(:,1) > sorted_upper(pulseindex2 - 1)));
spikelist2 = cat(1,spikelist2,TS(indexspike2));
end

%Remove zeros from spikelist and spikelist2
spikelist(any(spikelist==0,2),:)=[];
spikelist2(any(spikelist2==0,2),:)=[];

%Alter structure of DataPoints and reduce to only 1 lead 32x1xY
A = DataPoints(:,1,:);
a = squeeze(A);
a_trans = a';
B = DataPoints(:,2,:);
b = squeeze (B);
b_trans = b';
C = DataPoints(:,3,:);
c = squeeze (C);
c_trans = c';
D = DataPoints(:,4,:);
d = squeeze (D);
d_trans = d';

%Combine TimeStamps and DataPoints for each lead
a_run = horzcat(TimeStamps, a_trans);
b_run = horzcat(TimeStamps, b_trans);
c_run = horzcat(TimeStamps, c_trans);
d_run = horzcat(TimeStamps, d_trans);

%Variable pre-allocation for spike data (following section)
a_val = zeros(10000,32);
b_val = zeros(10000,32);
c_val = zeros(10000,32);
d_val = zeros(10000,32);
a_val2 = zeros(10000,32);
b_val2 = zeros(10000,32);
c_val2 = zeros(10000,32);
d_val2 = zeros(10000,32);

n = 1;

%This is written to pull out the 32 data points for each spike on each lead
for index1 = 1:length(spikelist)
    if  any(a_run(:,1) == spikelist(index1));
        indexa = find(a_run(:,1) == spikelist(index1));
        a_val_temp = a_trans(indexa,:);
        a_val(n,:) = a_val_temp;
        b_val_temp = b_trans(indexa,:);
        b_val(n,:) = b_val_temp;
        c_val_temp = c_trans(indexa,:);
        c_val(n,:) = c_val_temp;
        d_val_temp = d_trans(indexa,:);
        d_val(n,:) = d_val_temp;
        n = n + 1;
    end
end

for index1 = 1:length(spikelist2)
    if  any(a_run(:,1) == spikelist2(index1));
        indexa = find(a_run(:,1) == spikelist2(index1));
        a_val_temp = a_trans(indexa,:);
        a_val2(n,:) = a_val_temp;
        b_val_temp = b_trans(indexa,:);
        b_val2(n,:) = b_val_temp;
        c_val_temp = c_trans(indexa,:);
        c_val2(n,:) = c_val_temp;
        d_val_temp = d_trans(indexa,:);
        d_val2(n,:) = d_val_temp;
        n = n + 1;
    end
end

%Remove zeros from x_val variables
a_val(all(a_val==0,2),:)=[];
b_val(all(b_val==0,2),:)=[];
c_val(all(c_val==0,2),:)=[];
d_val(all(d_val==0,2),:)=[];

a_val2(all(a_val2==0,2),:)=[];
b_val2(all(b_val2==0,2),:)=[];
c_val2(all(c_val2==0,2),:)=[];
d_val2(all(d_val2==0,2),:)=[];

%For plotting purposes for 32 data points
x_val = 1:32;
x_val = x_val';

%Calculating mean data point values for all spikes on each lead
a_valmean = mean(a_val);
a_val2mean = mean(a_val2);
b_valmean = mean(b_val);
b_val2mean = mean(b_val2);
c_valmean = mean(c_val);
c_val2mean = mean(c_val2);
d_valmean = mean(d_val);
d_val2mean = mean(d_val2);

%Cross correlations between light (val) and non-light (val2) spikes on each lead
corrsa = xcorr(a_valmean, a_val2mean, 0, 'coeff');
corrsb = xcorr(b_valmean, b_val2mean, 0, 'coeff');
corrsc = xcorr(c_valmean, c_val2mean, 0, 'coeff');
corrsd = xcorr(d_valmean, d_val2mean, 0, 'coeff');

%energy calculations first - square mean values
a_valmeansq = mean(a_val).^2;
a_val2meansq = mean(a_val2).^2;
b_valmeansq = mean(b_val).^2;
b_val2meansq = mean(b_val2).^2;
c_valmeansq = mean(c_val).^2;
c_val2meansq = mean(c_val2).^2;
d_valmeansq = mean(d_val).^2;
d_val2meansq = mean(d_val2).^2;

%energy calculations second - use matlab trapz
inta1 = trapz(x_val, a_valmeansq');
inta2 = trapz(x_val, a_val2meansq');
intb1 = trapz(x_val, b_valmeansq');
intb2 = trapz(x_val, b_val2meansq');
intc1 = trapz(x_val, c_valmeansq');
intc2 = trapz(x_val, c_val2meansq');
intd1 = trapz(x_val, d_valmeansq');
intd2 = trapz(x_val, d_val2meansq');


x = 1:32;

figure

%Optical ID Session - Waveform Plots (Light Stimulated)
subplot(2,6,1)
plot(x, a_valmean, 'k-')
xlabel('Time')
ylabel('Microvolts')
title('Light-Stim: Lead 1')

subplot(2,6,2)
plot(x, b_valmean, 'k-')
xlabel('Time')
ylabel('Microvolts')
title('Light-Stim: Lead 2')

subplot(2,6,3)
plot(x, c_valmean, 'k-')
xlabel('Time')
ylabel('Microvolts')
title('Light-Stim: Lead 3')

subplot(2,6,4)
plot(x, d_valmean, 'k-')
xlabel('Time')
ylabel('Microvolts')
title('Light-Stim: Lead 4')

%Optical ID Session - Energy Data (Light Stimulated)
ax = subplot(2,6,5);
text(0.1, 1, 'Lead 1');
text(0.1, .75, 'Lead 2');
text(0.1, .50, 'Lead 3');
text(0.1, .25, 'Lead 4');
text(0.5, 1.1, 'Light-Stim: Energies');
text(0.5, 1, num2str(inta1));
text(0.5, 0.75, num2str(intb1));
text(0.5,0.5, num2str(intc1));
text(0.5,0.25,num2str(intd1));
set(ax,'visible','off');

%Optical ID Session - Correlation (to non-light stimulated) Data
ax = subplot(2,6,6);
text(0.1, 1, 'Lead 1');
text(0.1, .75, 'Lead 2');
text(0.1, .50, 'Lead 3');
text(0.1, .25, 'Lead 4');
text(0.5, 1.1, 'Correlations');
text(0.5, 1, num2str(corrsa));
text(0.5, 0.75, num2str(corrsb));
text(0.5,0.5, num2str(corrsc));
text(0.5,0.25,num2str(corrsd));
set(ax,'visible','off');

%Optical ID Session - Waveform Plots (non-light stimulated)
subplot(2,6,7)
plot(x, a_val2mean)
xlabel('Time')
ylabel('Microvolts')
title('No-Light: Lead 1')

subplot(2,6,8)
plot(x, b_val2mean)
xlabel('Time')
ylabel('Microvolts')
title('No-Light: Lead 2')

subplot(2,6,9)
plot(x, c_val2mean)
xlabel('Time')
ylabel('Microvolts')
title('No-Light: Lead 3')

subplot(2,6,10)
plot(x, d_val2mean)
xlabel('Time')
ylabel('Microvolts')
title('No-Light: Lead 4')

%Optical ID Session - Energy Data (non-light stimulated)
ax = subplot(2,6,11);
text(0.1, 1, 'Lead 1');
text(0.1, .75, 'Lead 2');
text(0.1, .50, 'Lead 3');
text(0.1, .25, 'Lead 4');
text(0.5, 1.1, 'No-Light: Energies');
text(0.5, 1, num2str(inta2));
text(0.5, 0.75, num2str(intb2));
text(0.5,0.5, num2str(intc2));
text(0.5,0.25,num2str(intd2));
set(ax,'visible','off');

end
