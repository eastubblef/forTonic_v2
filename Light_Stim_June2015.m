function Light_Stim_June2015(taskbase_file, t_file)
%% Written by Andrew Wolf - given to BS 6/23/15
%This is simple code that I used to look at the pattern of neuron firing following
%the initiation of the light pulses (which in my case where at 8 hz. The
%neurons that were successfully optogenetically stimulated displayed a
%characteristic 8 peak response in the 1 second following the start of the
%light pulses

%Load necessary files 
load(taskbase_file);
x = t_to_mat(t_file);

%Convert imported times
x = x / 10000;
spikes = x;

%Pull out relative times of first light pulse for each trial
vector = cell2mat(taskbase.pulse_on);
firstpulse = vector(:,1);

%Transpose vectors to fit function input requirements
spikes = spikes';
firstpulse = firstpulse';
taskbase.start_nlx = taskbase.start_nlx';

%Produce raster based on GF function (+/- 1 second)
[ref_spike_times, trial_inds, plot_handle] = rasterGF(spikes, taskbase.start_nlx(1:end-1), firstpulse, [-1 1]);

%Produce PSTH based on GF function (+/- 1 second)
[ref_psth, plot_handle, psth_info] = psthGF(ref_spike_times, trial_inds, [-1 1], 10);

end

