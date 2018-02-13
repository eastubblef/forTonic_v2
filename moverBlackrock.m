
%% From Katie 3/15
% Most helpful file is TNC_MoverBehaviorExtract. 

% How to use: 
% ContData = TNC_MoverBehaviorExtract(filenamestr, targetName, dataRate, chan)
% 
% filenamestr: path to the file, for example, '/Users/dudmanj/Desktop/m42/m42-20110427-wide-plustim-008.ns5'
% targetName: will automatically save the data as this, can be whatever you want 
% dataRate: data type, needs to be 'ns4', 'ns3' or 'ns2' depending on file
% chan: For each channel, there will be a number associated with it. They'll be different depending on your setup. You need to find chan.rew, chan.thr, chan.lick, chan.x, chan.y. 
% You should be able to figure this numbers out by running HeadFixedBehaviorLoadingScript- specifically the Load the continuous behavior monitor channels section
% Within the dataAinp data file, look under dataAinp.MetaTags.ChannelID. Those numbers are the different numbers associated with specific blackrock channels.
% 
% Summary:
% Files from TONIC you need to have: TNC_MoverBehaviorExtract, HeadFixedBehaviorLoadingScript, openNSx, openNEV
% Use HeadFixedBehaviorLoadingScript to get information for the chan variable in TNC_MoverBehaviorExtract
% Create chan variable to be used in TNC_MoverBehaviorExtract- you'll need chan.rew, chan.thr, chan.lick, chan.x, chan.y (those are the variables currently in TNC_MoverBehaviorExtract)


