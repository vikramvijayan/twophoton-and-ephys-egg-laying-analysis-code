% THIS IS ONLY STABLE CLASS 2 whole cell 
% little to no current injection
% there is often some drift in these recordings...
% >5 minutes of recoridng

% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list, fly_ID, cell_ID] = import_ROI_names_SS53469_patching_classII(ROIsSS53469, whatdoyouwant)

%function [recording] = make_small_recording_group()


recording_list = {'2019_09_17_0006_spikes.mat','2020_09_14_0002_spikes.mat','2020_09_14_0003_spikes.mat','2020_09_14_0004_spikes.mat'};
ROI_list = [-1,-1,-1,-1];
fly_ID = [3,12,12,12];
cell_ID = [3,12,12,12];

% imported ROIsSS53469 as string

rec_list = recording_list;


end