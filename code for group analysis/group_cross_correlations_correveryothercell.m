function [corr_structure] = group_cross_correlations_correveryothercell(recordings_list, ROI_num, fly_ID, cell_ID, backsub)

%% plotting correlations or time locked average

for rec_index = 1:2:length(recordings_list)
    % it is often faster to use the stripped files without the associated
    % images
    modifiedStr = strrep([char(recordings_list(rec_index))], '.mat', '_stripped.mat');
    modifiedStr2 = strrep([char(recordings_list(rec_index))], '.mat', '_stripped_new_dlc_track.mat');
    
    if(exist([modifiedStr2]))
        modifiedStr = modifiedStr2;
    end
    %modifiedStr = [char(recordings_list(rec_index))];
    recording = loadsinglerecording(modifiedStr);
    
    [recording] = processes_more_DLC_variables(recording);
    if(backsub == 1)
        [recording] = replace_df_over_f_withbackgroundsubtracted(recording);
    end
    
          if(backsub == 2)
        [recording] = replace_df_over_f_withbackgroundsubtracted_runningmean(recording);
        end
        
          if(backsub == 3)
        [recording] = replace_df_over_f_withbackgroundsubtracted_runningmean2(recording);
          end
          
    
    patch = ~ROI_num(rec_index);
    
    ROI_num1 = ROI_num(rec_index);
    ROI_num2 = ROI_num(rec_index+1);
    
    corr_of_signals = [];
    
    for m = 1:1:length(recording.tseries)
        signal1    = recording.tseries(m).df_over_f(ROI_num1,:)./nanmean(recording.tseries(m).df_over_f(ROI_num1,:));
        signal2   = recording.tseries(m).df_over_f(ROI_num2,:)./nanmean(recording.tseries(m).df_over_f(ROI_num2,:));
        
%          signal1    = recording.tseries(m).mean_in_ROI(ROI_num1,:)./nanmean(recording.tseries(m).mean_in_ROI(ROI_num1,:))-1;
%          signal2   = recording.tseries(m).mean_in_ROI(ROI_num2,:)./nanmean(recording.tseries(m).mean_in_ROI(ROI_num2,:))-1;
        signaltime = recording.tseries(m).Time_s;
        
        
        % look at signal correlated to body position
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal1, signal2, signaltime, signaltime, [-120:.1:120]);
        corr_of_signals = [corr_of_signals; out_corr];
        time_base_corr = out_corr_shift_ofvec1;
        
    end
    
    %%%% the loop above can be replaced with new vweaions of the
    %%%% plot_cross_correlation function. It is ahared

    corr_structure(rec_index).recording = recordings_list(rec_index);
    corr_structure(rec_index).recording_index = rec_index;
    corr_structure(rec_index).ROI_used = ROI_num(rec_index);
    corr_structure(rec_index).corr_to_signal = nanmean(corr_of_signals,1)';
    corr_structure(rec_index).time_base_corr = time_base_corr;
        corr_structure(rec_index).fly_ID = fly_ID(rec_index);

            corr_structure(rec_index).cell_ID = cell_ID(rec_index);

    
end
end