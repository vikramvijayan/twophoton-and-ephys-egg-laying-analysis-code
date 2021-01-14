
% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list, fly_ID, cell_ID] = import_ROI_names_SS53469_patching_VTmutant(ROIsSS53469, whatdoyouwant)

%function [recording] = make_small_recording_group()

recording_list = {'2020_10_13_0004_spikes.mat','2020_10_13_0008_spikes.mat','2020_10_13_0009_spikes.mat','2020_10_13_0013_spikes.mat'};
ROI_list = [-1,-1,-1,-1];
fly_ID = [35,36,36,37];
cell_ID = [35,36,36,37];

% imported ROIsSS53469 as string

rec_list = recording_list;


end