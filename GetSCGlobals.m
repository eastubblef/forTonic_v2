%% GetSCGlobals
%% Central location where global variables are defined for coordination across
%% analyses.
%% 1/8/08

function GetSCGlobals

global SELECTIVITY_BATCHFILE; % the file containing all of the selectivity calculations for each cell (should be the most recent one calculated)
global SELECTIVITY_BATCHFILE_GT; % as above, but for go tone data
global MAX_MOVEMENT_TIME_TO_WATER; % max time from odorpokeout to waterpokein (discard trials where time is longer than this)
global MAX_MOVEMENT_TIME_TO_ODOR; % max time from waterpokeout to nextodorpokein
global MAX_TONERESPONSE_TIME; % max time from gotoneon to odorpokeout
global INTERPOKE_DIST; % in cm; distance between center and side ports 

global TRIAL_EVENTS; % numbers in saved raster_info corresponding to the trial events

% for epoch definitions (for selectivity analyses)
global PREMOVEMENT_TIME;
global PREGOTONE_TIME;
global POSTWATER_TIME;
global PREODORPOKE_TIME;
global epochs; % define epochs for selectivity analyses


global sc_layer_depths; % the estimated depths of SC layers (in mm)
global MM_PER_TURN; % mm per 1/8th turn of tetrode screw

global RESOLUTION; % set to 1000 to plot spike times to 1 msec resolution

global FAMILIAR_IND; %which mix_ind is the familiar and which is the novel odor
global NOVEL_IND; %which mix_ind is the familiar and which is the novel odor
global MIXTURE_IND; %which mix_ind is the mixture and which is the pure odor
global PURE_IND; %which mix_ind is the mixture and which is the pure odor

global FIRING_RATE_MIN; % minimum firing rate (per "category"); cells below this will be excluded
global NUM_TRIALS_MIN; % minimum number of trials (per "category"); cells below this will be excluded

global INVALID_VALUE_FLAG;  % flag for trials to be marked invalid (e.g., b/c of early waterpokeout (before epoch start))
global MAX_WVO_DELAY; % the maximum delay between waterpokein and watervalveon


global GT_CUTOFF; % go tone delays shorter than this are "short" trials; longer are "long" trials

global JOURNAL_NAME; % the name of the journal to format panel labels for


% base these values on data in Movement_times_by_session.ppt
MAX_MOVEMENT_TIME_TO_WATER = 1; % in sec
MAX_MOVEMENT_TIME_TO_ODOR = 1.5; % in sec
MAX_TONERESPONSE_TIME = 0.3; % in sec

%INTERPOKE_DIST = 5.4; % in cm; distance between center and side ports
INTERPOKE_DIST = 5.7; % in cm; distance between center and side ports - updated 3/24/08

% for epoch definitions (for selectivity analyses)
PREMOVEMENT_TIME = 0.1; % in sec
PREGOTONE_TIME = 0.5; % in sec
POSTWATER_TIME = 1; % in sec
PREODORPOKE_TIME = 0.5; % in sec

FIRING_RATE_MIN = 2; % spks/sec
NUM_TRIALS_MIN = 4;
%NUM_TRIALS_MIN = 10;

INVALID_VALUE_FLAG = 99999;
MAX_WVO_DELAY = 0.5;

% Define the trial events
TRIAL_EVENTS.OdorPokeIn = 1;
TRIAL_EVENTS.OdorValveOn = 2;
TRIAL_EVENTS.OdorPokeOut = 3;
TRIAL_EVENTS.WaterPokeIn = 4;
TRIAL_EVENTS.WaterValveOn = 5;
TRIAL_EVENTS.WaterPokeOut = 6;
TRIAL_EVENTS.NextOdorPokeIn = 7;
TRIAL_EVENTS.GoToneOn = 8; % NOTE: not chronological!


%% Define the epochs here - note that they are no longer chronological, due
%% to GoToneOn addition in TRIAL_EVENTS
epochs{1} = [888 -PREODORPOKE_TIME TRIAL_EVENTS.OdorPokeIn]; %Pre-OdorPoke (to test for bias)
epochs{2} = [TRIAL_EVENTS.OdorValveOn TRIAL_EVENTS.OdorPokeOut]; %Stimulus: OdorValveOn to OdorPokeOut
epochs{3} = [888 -PREMOVEMENT_TIME TRIAL_EVENTS.OdorPokeOut]; %Pre-movement: The PREMOVEMENT_TIME msec preceding OdorPokeOut
epochs{4} = [TRIAL_EVENTS.OdorPokeOut TRIAL_EVENTS.WaterPokeIn]; % Movement: OdorPokeOut to WaterPokeIn
epochs{5} = [999 TRIAL_EVENTS.WaterPokeIn POSTWATER_TIME]; % Reward: WaterPokeIn to POSTWATER_TIME
epochs{6} = [TRIAL_EVENTS.WaterPokeIn TRIAL_EVENTS.WaterPokeOut]; % Reward (alternative): WaterPokeIn to WaterPokeOut
epochs{7} = [888 -PREMOVEMENT_TIME TRIAL_EVENTS.WaterPokeOut]; % Pre-reinitiation movement: The PREMOVEMENT_TIME msec preceding WaterPokeOut
epochs{8} = [TRIAL_EVENTS.WaterPokeOut TRIAL_EVENTS.NextOdorPokeIn]; % Reinitiation movement: WaterPokeOut to NextOdorPokeIn

epochs{9} = [TRIAL_EVENTS.OdorValveOn TRIAL_EVENTS.GoToneOn]; % mvmt planning activity: odorvalveon to gotoneon
epochs{10} = [888 -PREGOTONE_TIME TRIAL_EVENTS.GoToneOn]; % mvmt planning activity: PREGOTONE_TIME to gotoneon
epochs{11} = [TRIAL_EVENTS.GoToneOn TRIAL_EVENTS.OdorPokeOut]; % mvmt initiation activity: gotoneon to odorpokeout
epochs{12} = [888 -1 TRIAL_EVENTS.GoToneOn 888 -0.5 TRIAL_EVENTS.GoToneOn]; % 1000-500 ms before gotone
    % (specific for short vs. long delay comparisons)

% % specifically for 'at water' epoch - different atwater epoch lengths
% epochs{1} = [999 4 0.25]; % Reward: WaterPokeIn for xx s
% epochs{2} = [999 4 0.5]; % Reward: WaterPokeIn for xx s
% epochs{3} = [999 4 0.75]; % Reward: WaterPokeIn for xx s
% epochs{4} = [999 4 1]; % Reward: WaterPokeIn for xx s

% Define the estimated depths of SC layers (in mm)
sc_layer_depths.zonal_top = 2.6;
sc_layer_depths.superficial_top = 2.7;
sc_layer_depths.intermediate_top = 3.3;
sc_layer_depths.deep_top = 4;
sc_layer_depths.deep_bottom = 5.4;

MM_PER_TURN = 0.04;

RESOLUTION = 1000;

FAMILIAR_IND = 1; % in ratbase structures, familiar mix_ind is 1 and novel mix_ind is 2
NOVEL_IND = 2; % in ratbase structures, familiar mix_ind is 1 and novel mix_ind is 2
MIXTURE_IND = 1; % in ratbase structures, mixture mix_ind is 1 and pure mix_ind is 2
PURE_IND = 2; % in ratbase structures, mixture mix_ind is 1 and pure mix_ind is 2

GT_CUTOFF = 0.75;

%JOURNAL_NAME = 'nature_neuroscience';
JOURNAL_NAME = 'neuron';


%all included, AND 'SALVAGED', nonpoor G017, JE7, JC3, and G042 cells, 1000
%iterations, invalid trials excluded; some more new analyses included; activity
%info for each cell included; atwater and outcome (including L/R) selectivity info included (different
%lengths of 'outcome' epoch), with restricted epochs (epoch ends w/ waterpokeout, instead of 'natural' end); 
% additional bin sizes for sliding selectivity added to address reviewer
% 1's comments - This was the batchfile used for the final
% submission of Felsen & Mainen, 2008
%SELECTIVITY_BATCHFILE = 'C:\Gidon\mainen_projects\superior_colliculus\batch_results\CompileSelectivities_20080716T092846';

% Also includes leftchoice_vs_rightchoice_easy, leftchoice_vs_rightchoice_hard, leftchoice_vs_rightchoice_pair1, and leftchoice_vs_rightchoice_pair2 analyses
%SELECTIVITY_BATCHFILE = 'C:\Gidon\mainen_projects\superior_colliculus\batch_results\CompileSelectivities_20090415T152604';
% And include some other new selectivity analyses
SELECTIVITY_BATCHFILE = 'C:\Gidon\mainen_projects\superior_colliculus\batch_results\CompileSelectivities_20091216T114600';

%% load a separate batchfile for the gotone data
%SELECTIVITY_BATCHFILE_GT = 'C:\Gidon\mainen_projects\superior_colliculus\batch_results\CompileSelectivities_20090108T094209';

%all included, salvaged, and recut nonpoor G046, G047, G057, and G059 cells
%SELECTIVITY_BATCHFILE_GT = 'C:\Gidon\mainen_projects\superior_colliculus\batch_results\CompileSelectivities_20090313T154643';
%SELECTIVITY_BATCHFILE_GT = 'C:\Gidon\mainen_projects\superior_colliculus\batch_results\CompileSelectivities_20090320T142218';
%SELECTIVITY_BATCHFILE_GT = 'C:\Gidon\mainen_projects\superior_colliculus\batch_results\CompileSelectivities_gt_20090410T105133';
%SELECTIVITY_BATCHFILE_GT = 'C:\Gidon\mainen_projects\superior_colliculus\batch_results\CompileSelectivities_gt_20090421T100355';
%SELECTIVITY_BATCHFILE_GT = 'C:\Gidon\mainen_projects\superior_colliculus\batch_results\CompileSelectivities_gt_20090506T130425';
%SELECTIVITY_BATCHFILE_GT = 'C:\Gidon\mainen_projects\superior_colliculus\batch_results\CompileSelectivities_gt_20090824T231747';
%SELECTIVITY_BATCHFILE_GT = 'C:\Gidon\mainen_projects\superior_colliculus\batch_results\CompileSelectivities_gt_20091104T101201';
%SELECTIVITY_BATCHFILE_GT = 'C:\Gidon\mainen_projects\superior_colliculus\batch_results\CompileSelectivities_gt_20091214T134911';
SELECTIVITY_BATCHFILE_GT = 'C:\Gidon\mainen_projects\superior_colliculus\batch_results\CompileSelectivities_gt_20100512T120333';

