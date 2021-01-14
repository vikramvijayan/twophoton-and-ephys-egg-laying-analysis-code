% THIS IS ONLY STABLE CLASS 1 whole cell that are resting at 4+
% spikes/sec
% little to no current injection
% >15 minutes of recoridng

% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list, fly_ID, cell_ID] = import_ROI_names_SS53469_patching(ROIsSS53469, whatdoyouwant)

%function [recording] = make_small_recording_group()


% recording_list = {'2019_09_17_0003_spikes.mat','2020_09_09_0000_spikes.mat','2020_09_17_0001_spikes.mat'};
% ROI_list = [-1,-1,-1];

% incldudes a 1 spike/sec recording
% use this from now on
recording_list = {'2019_09_17_0003_spikes.mat','2020_09_09_0000_spikes.mat','2020_09_17_0001_spikes.mat','2020_09_15_0000_spikes.mat','2020_10_15_0001_spikes.mat','2020_10_15_0002_spikes.mat','2020_10_15_0003_spikes.mat'};
ROI_list = [-1,-1,-1,-1,-1,-1,-1];
fly_ID = [1,10,18,13,39,39,39];
cell_ID = [1,10,18,13,39,39,39];


% imported ROIsSS53469 as string

rec_list = recording_list;


end