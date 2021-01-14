function [out, len] = group_get_frame_rates(recordings_list)


%% plotting correlations or time locked average
for rec_index = 1:1:length(recordings_list)
    tf = strfind([char(recordings_list(rec_index))],'spikes');
  modifiedStr2 = [char(recordings_list(rec_index))];

if(isempty(tf))
    modifiedStr = strrep([char(recordings_list(rec_index))], '.mat', '_stripped.mat');
    modifiedStr2 = strrep([char(recordings_list(rec_index))], '.mat', '_stripped_new_dlc_track.mat');
    end

    if(exist([modifiedStr2]))
        modifiedStr = modifiedStr2;
    end
    % modifiedStr = [char(recordings_list(rec_index))];
    recording = loadsinglerecording(modifiedStr);
    
    out(rec_index) = recording.tseries(1).median_imaging_interval;
    len(rec_index) = recording.tseries(1).Time_s(end)-recording.tseries(1).Time_s(1);
end