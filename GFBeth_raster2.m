function [ST, trial_start_times, event_times, trial_numb] = Beth_raster2(spk_times_fn, tb_fn, T);

load(spk_times_fn);

load(tb_fn);

%%
ST = TS/10000;
taskbase.freq_mat = cell2mat(taskbase.freq)

real_taskbase.start = taskbase.start(1:end-1)
taskbase.freq_mat = taskbase.freq_mat(1:length(real_taskbase.start(1:end-1)))

trial_numb = find(taskbase.freq_mat == [T]);

trial_start_times = (taskbase.start_nlx(trial_numb));
pulse_on1 = taskbase.pulse_on(trial_numb)

%% Turn pulse_on1 into a matrix
%Find max number of pulse ons (or the length) within the cells: pulse_on1
for trial_ind = 1:length(trial_numb)
    num_on(trial_ind) = length(pulse_on1{trial_ind});
end

max_length = max(num_on);
pulse_mat = nan(length(trial_numb), max_length); %Initialize the matrix

for trial_ind = 1:length(trial_numb)
    pulse_mat(trial_ind, 1:num_on(trial_ind)) = pulse_on1{trial_ind};
end

%% Index relative times in terms of start times and trial number. 
%Since taskbase.pulse_on times are within cells, getting 1st element per cell requires:
for i = 1:length(trial_numb);
    trial_numb_ind = trial_numb(i);
    pulse1(i) = taskbase.pulse_on{trial_numb_ind}(1);
    pulse_one = pulse1';
end

event_times = pulse_one;

return;

