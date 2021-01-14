
% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list, fly_ID, cell_ID] = import_ROI_names_SS53469_patching_VTkir(ROIsSS53469, whatdoyouwant)

%function [recording] = make_small_recording_group()


recording_list = {'2020_10_08_0007_spikes.mat','2020_10_08_0008_spikes.mat','2020_10_08_0011_spikes.mat','2020_10_08_0018_spikes.mat','2020_10_08_0024_spikes.mat'};
ROI_list = [-1,-1,-1,-1,-1];

fly_ID = [26,28,28,30,30];
cell_ID = [27,28,28,30,30];
% imported ROIsSS53469 as string

rec_list = recording_list;


end