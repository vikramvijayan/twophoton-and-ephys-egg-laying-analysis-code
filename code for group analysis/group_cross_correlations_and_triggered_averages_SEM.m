function [corr_structure,egg_structure] = group_cross_correlations_and_triggered_averages_SEM(recordings_list, ROI_num, fly_ID, cell_ID, filter_out_chrimson_data, filter_chrimson_eggs, filter_out_egg_times, backsub)

%% note that signal over mean is being currently used to calculate cross correlations

%% make sure to set the time interval for the calculations (triggered averages)
% time_interval = [-240:.1:240];
time_interval = [-1200:.1:1200];
time_interval_mid = (length(time_interval)-1)./2;

%% make sure to set the time interval for cross correaltions
time_interval_corr = [-120:.1:120];

%% description of inputs
% set filter chrimson data to 1 if you want to remove all data during pulses
% (as well as 100 sec after the pulse) check code for exact number

% set filter chrimson eggs to remove all eggs in the same interval as above
% 0 = no filter
% 1 - filter out chrimson eggs
% 2 - only use chrimson eggs

% set ROI num to -1 if you want tu use CH1 patch
% if patch is 1 then use CH1
% if patch is 0 then use CH1 but remove pulses (and 100 seconds after)

% set filter_out_egg_times to 1 if you don't want to have any of the egg
% laying times inclued (currently 300 sec on either side of egg)


%% plotting correlations or time locked average
for rec_index = 1:1:length(recordings_list)
    
    % it is often faster to use the stripped files without the associated
    % images
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
    
    if(ROI_num(rec_index) ==-1)
        patch = 1;
    else
        patch = 0;
    end
    if(patch == 1)
        recording.tseries = 1;
    end
           
    %% initialize arrays for triggered averages
    
    data_to_plain   = [];
    data_to_sucrose = [];
    
    data_to_sucrose_control = [];
    data_to_plain_control = [];
    
    data_around_egg_vel = [];
    data_around_egg_vel_noave = [];
    data_around_egg_speed_2hz_hold = [];
    
    data_around_egg_dfoverf = [];
    data_around_egg_dfoverf_mean = [];
    data_around_egg_dfoverf_fo_is_mean = [];
    
    data_around_egg_sub = [];
    
    data_around_egg_body = [];
    data_around_egg_body_x = [];
    data_around_egg_body_y = [];
    
    data_around_egg_prob = [];
    data_around_egg_prob_x = [];
    data_around_egg_prob_y = [];
    
    data_around_egg_path = [];
    data_around_egg_body_only = [];
    data_around_egg_body_only_x = [];
    data_around_egg_body_only_y = [];
    
    data_around_egg_angle = [];
    data_around_egg_angle_only = [];
    
    data_around_egg_PWM = [];
    
    %% initialize arrays for correlations
    
    corr_to_vel = [];
    corr_to_vel_noave = [];
    corr_to_speed_2hz_hold = [];
    
    corr_to_body = [];
    corr_to_body_x = [];
    corr_to_body_y = [];
    
    corr_to_prob = [];
    corr_to_prob_x = [];
    corr_to_prob_y = [];
    
    corr_to_angle = [];
    corr_to_body_only = [];
    corr_to_body_only_x = [];
    corr_to_body_only_y = [];
    corr_to_body_path = [];   
    corr_to_angle_only = [];
    
    corr_to_sucrose = [];   
    corr_to_sucrose_t = [];
    corr_to_trans = []; 
    corr_to_absvaltrans = []; 
    corr_to_egg = [];
    corr_to_plainegg = []; 
    corr_to_sucroseegg = [];
    
    %% process more of the DLC movie (performing median normalization of lengths)
    
    abd_only_length_x = recording.movie1.abd_x_L3tip./nanmedian(recording.movie1.abd_x_L3tip);
    bodylength_x =  recording.movie1.abd_x_neck_tip./nanmedian(recording.movie1.abd_x_neck_tip);
    abd_only_length_y = recording.movie1.abd_y_L3tip./nanmedian(recording.movie1.abd_y_L3tip);
    bodylength_y =  recording.movie1.abd_y_neck_tip./nanmedian(recording.movie1.abd_y_neck_tip);
    
    prob_x =  recording.movie1.prob_x./nanmedian(recording.movie1.prob_x);
    prob_y =  recording.movie1.prob_y./nanmedian(recording.movie1.prob_y);
    
    movtime  = recording.movie1.time_stamps;
    sucrose = recording.movie1.sucrose;
    bodylength =  recording.movie1.abd_length./nanmedian(recording.movie1.abd_length);
    bodyangle =  recording.movie1.abd_angle;
    problength = recording.movie1.prob_length./nanmedian(recording.movie1.prob_length);
    abd_path_length = recording.movie1.abd_path_length./nanmedian(recording.movie1.abd_path_length);
    abd_only_length = recording.movie1.abd_only_length./nanmedian(recording.movie1.abd_only_length);
    abd_only_angle = recording.movie1.abd_only_angle;
    wheelvelocity = (pi*7*-25).*smooth([0; diff(unwrap(recording.movie1.filtered_wheel))],25) ./ (2*pi); %mm/sec 7mm radium
    wheelvelocity_noave = (pi*7*-25).*[0; diff(unwrap(recording.movie1.filtered_wheel))]./ (2*pi); %mm/sec 7mm radium
    
    %% this is new code to get a 2hz speed
    tmp_unwrapped_wheel = unwrap(recording.movie1.filtered_wheel);
    tmp_unwrapped_wheel_2hz = [];
    tmp_valx =tmp_unwrapped_wheel(1);
    
    % remaking a 2 hz array
    for i = 1:1:length(tmp_unwrapped_wheel)
        if(mod(i,25) == 13)
            tmp_valx = (tmp_unwrapped_wheel(i)+tmp_unwrapped_wheel(i-1))./2;
        end
        if(mod(i,25) == 0)
            tmp_valx = (tmp_unwrapped_wheel(i));
        end
        
        tmp_unwrapped_wheel_2hz(i) = tmp_valx;
    end
    tmp_unwrapped_wheel_2hz_speed = [0, (pi*7*-25).*abs(diff(tmp_unwrapped_wheel_2hz))./ (2*pi)];
    

    
    % this holds the value of the speed (so its not jumping between 0 and a
    % value). Since we are holding we are dividing by 12.5
    tmp_valsp = tmp_unwrapped_wheel_2hz_speed(1)./12.5;
    speed_2hz_hold = [];
    
    for i = 1:1:length(tmp_unwrapped_wheel_2hz_speed)
        if(tmp_unwrapped_wheel_2hz_speed(i) ~=0)
            tmp_valsp = tmp_unwrapped_wheel_2hz_speed(i)./12.5;
        end
        speed_2hz_hold(i) = tmp_valsp;
    end
    speed_2hz_hold = -1.*speed_2hz_hold';
    
        % prevent errant issues
    speed_2hz_hold(speed_2hz_hold>5) = 0;
    
    if(filter_chrimson_eggs == 1)
        recording.movie1.eggs = recording.movie1.eggs(recording.movie1.eggs_pulse == 0);
    end
    if(filter_chrimson_eggs == 2)
        recording.movie1.eggs = recording.movie1.eggs(recording.movie1.eggs_pulse == 1);
    end
    
    [a b] = find(recording.movie1.sucrose(recording.movie1.eggs) == 0);
    plain_egg = recording.movie1.eggs(b);
    plain_eggt = recording.movie1.time_stamps(recording.movie1.eggs(b));
    
    [a b] = find(recording.movie1.sucrose(recording.movie1.eggs) == 200 | recording.movie1.sucrose(recording.movie1.eggs) == 500 );
    sucrose_egg = recording.movie1.eggs(b);
    sucrose_eggt = recording.movie1.time_stamps(recording.movie1.eggs(b));
    
    [a b] = find(recording.movie1.sucrose(recording.movie1.eggs) >= 0);
    all_egg = recording.movie1.eggs(b);
    all_eggt = recording.movie1.time_stamps(recording.movie1.eggs(b));
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,all_eggt, time_interval);
    data_around_egg_sub = [data_around_egg_sub; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(wheelvelocity,movtime,all_eggt, time_interval);
    data_around_egg_vel = [data_around_egg_vel; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(wheelvelocity_noave,movtime,all_eggt, time_interval);
    data_around_egg_vel_noave = [data_around_egg_vel_noave; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(speed_2hz_hold,movtime,all_eggt, time_interval);
    data_around_egg_speed_2hz_hold = [data_around_egg_speed_2hz_hold; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(problength,movtime,all_eggt, time_interval);
    data_around_egg_prob = [data_around_egg_prob; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(prob_x,movtime,all_eggt, time_interval);
    data_around_egg_prob_x = [data_around_egg_prob_x; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(prob_y,movtime,all_eggt, time_interval);
    data_around_egg_prob_y = [data_around_egg_prob_y; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(bodylength,movtime,all_eggt, time_interval);
    data_around_egg_body = [data_around_egg_body; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(bodylength_x,movtime,all_eggt, time_interval);
    data_around_egg_body_x = [data_around_egg_body_x; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(bodylength_y,movtime,all_eggt, time_interval);
    data_around_egg_body_y = [data_around_egg_body_y; data_to_average_interp];
    
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_path_length,movtime,all_eggt, time_interval);
    data_around_egg_path = [data_around_egg_path; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_length,movtime,all_eggt, time_interval);
    data_around_egg_body_only = [data_around_egg_body_only; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_length_x,movtime,all_eggt, time_interval);
    data_around_egg_body_only_x = [data_around_egg_body_only_x; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_length_y,movtime,all_eggt, time_interval);
    data_around_egg_body_only_y = [data_around_egg_body_only_y; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(bodyangle,movtime,all_eggt, time_interval);
    data_around_egg_angle = [data_around_egg_angle; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_angle,movtime,all_eggt, time_interval);
    data_around_egg_angle_only = [data_around_egg_angle_only; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(recording.abf.PWMlaser_uWpermm2,recording.abf.Time_s,all_eggt, time_interval);
    data_around_egg_PWM = [data_around_egg_PWM; data_to_average_interp];
    
    data_around_egg_sub = repmat(data_around_egg_sub,length(recording.tseries),1);
    
    all_recordings_list = repmat(recordings_list(rec_index),length(recording.tseries),1);
    all_eggs_listindex = repmat(all_egg,length(recording.tseries),1);
    all_eggs_listtime = repmat(all_eggt,length(recording.tseries),1);
    all_recordings_list_index = repmat(rec_index,length(recording.tseries),1);
    all_ROI_index = repmat(ROI_num(rec_index),length(recording.tseries),1);
    all_cellID_list = repmat(cell_ID(rec_index),length(recording.tseries),1);
    all_flyID_list = repmat(fly_ID(rec_index),length(recording.tseries),1);

    for m = 1:1:length(recording.tseries)
        if(patch == 0)
            signal    = recording.tseries(m).df_over_f(ROI_num(rec_index),:);
            signal_over_mean    = recording.tseries(m).df_over_f(ROI_num(rec_index),:)./nanmean(recording.tseries(m).df_over_f(ROI_num(rec_index),:));
            signal_fo_is_mean  = recording.tseries(m).mean_in_ROI(ROI_num(rec_index),:)./nanmean(recording.tseries(m).mean_in_ROI(ROI_num(rec_index),:)) -1;
            
            
            signaltime = recording.tseries(m).Time_s;
            
            if(filter_out_chrimson_data)
                tseries_no_laser = interp1(recording.abf.Time_s, recording.abf.no_laser, recording.tseries(m).Time_s,'previous');
                signal = signal.*tseries_no_laser';
                signal_over_mean = signal_over_mean.*tseries_no_laser';
                signal_fo_is_mean = signal_over_mean.*tseries_no_laser';
            end
            
            if(filter_out_egg_times)
                tseries_no_egg = interp1(recording.movie1.time_stamps, recording.movie1.no_egg, recording.tseries(m).Time_s,'previous');
                signal = signal.*tseries_no_egg';
                signal_over_mean = signal_over_mean.*tseries_no_egg';
                signal_fo_is_mean = signal_over_mean.*tseries_no_egg';
            end
        end
        
        %% for patching data, the signal changes
        if(patch == 1)
            %signal    = recording.abf.CH1_patch(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000))-13;
            %signal    = 10000.*recording.abf.CH1_patch_spikes_conv_area_rect(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
            signal    = recording.abf.CH1_patch_spikes_removed(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000))-13;
            signaltime = recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
            
            if(filter_out_chrimson_data)
                signal    = recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)).*recording.abf.no_laser(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
            end
            
            if(filter_out_egg_times)
                disp 'Error this code has not been written'
            end
            signal_over_mean = signal;
            signal_fo_is_mean = signal;
        end
        
        %% compute egg triggered signal averages
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,all_eggt, time_interval);
        data_around_egg_dfoverf = [data_around_egg_dfoverf; data_to_average_interp];
        
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal_over_mean,signaltime,all_eggt, time_interval);
        data_around_egg_dfoverf_mean = [data_around_egg_dfoverf_mean; data_to_average_interp];
        
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal_fo_is_mean,signaltime,all_eggt, time_interval);
        data_around_egg_dfoverf_fo_is_mean = [data_around_egg_dfoverf_fo_is_mean; data_to_average_interp];
        
        %% compute correlations to signal
        
        % look at signal correlated to velocity (smoothed)
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, wheelvelocity, signaltime, movtime, time_interval_corr);
        corr_to_vel = [corr_to_vel; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, wheelvelocity_noave, signaltime, movtime, time_interval_corr);
        corr_to_vel_noave = [corr_to_vel_noave; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, speed_2hz_hold, signaltime, movtime, time_interval_corr);
        corr_to_speed_2hz_hold = [corr_to_speed_2hz_hold; out_corr];
        
        % look at signal correlated to body position
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, problength, signaltime, movtime, time_interval_corr);
        corr_to_prob = [corr_to_prob; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, prob_x, signaltime, movtime, time_interval_corr);
        corr_to_prob_x = [corr_to_prob_x; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, prob_y, signaltime, movtime, time_interval_corr);
        corr_to_prob_y = [corr_to_prob_y; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodylength, signaltime, movtime, time_interval_corr);
        corr_to_body = [corr_to_body; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodylength_x, signaltime, movtime, time_interval_corr);
        corr_to_body_x = [corr_to_body_x; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodylength_y, signaltime, movtime, time_interval_corr);
        corr_to_body_y = [corr_to_body_y; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodyangle, signaltime, movtime, time_interval_corr);
        corr_to_angle = [corr_to_angle; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_path_length, signaltime, movtime, time_interval_corr);
        corr_to_body_path = [corr_to_body_path; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_only_length, signaltime, movtime, time_interval_corr);
        corr_to_body_only = [corr_to_body_only; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_only_length_x, signaltime, movtime, time_interval_corr);
        corr_to_body_only_x = [corr_to_body_only_x; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_only_length_y, signaltime, movtime, time_interval_corr);
        corr_to_body_only_y = [corr_to_body_only_y; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_only_angle, signaltime, movtime, time_interval_corr);
        corr_to_angle_only = [corr_to_angle_only; out_corr];
        
        % look at signal correlated to substrate
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, sucrose, signaltime, movtime, time_interval_corr);
        corr_to_sucrose = [corr_to_sucrose; out_corr];
        
        % look at signal correlated to substrate threshold
        suc_thr = sucrose;
        suc_thr(suc_thr>0) = 1;
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, suc_thr, signaltime, movtime, time_interval_corr);
        corr_to_sucrose_t = [corr_to_sucrose_t; out_corr];
        
        % look at signal correlated to egg
        tmp_ar = zeros(1,length(movtime));
        tmp_ar(all_egg) = 1;
        tmp_ar = smooth(tmp_ar,25*10);
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, time_interval_corr);
        corr_to_egg = [corr_to_egg; out_corr];
        
        % look at signal correlated to plain egg
        tmp_ar = zeros(1,length(movtime));
        tmp_ar(plain_egg) = 1;
        tmp_ar = smooth(tmp_ar,25*10);
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, time_interval_corr);
        corr_to_plainegg = [corr_to_plainegg; out_corr];
        
        % look at signal correlated to sucrose egg
        tmp_ar = zeros(1,length(movtime));
        tmp_ar(sucrose_egg) = 1;
        tmp_ar = smooth(tmp_ar,25*10);
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, time_interval_corr);
        corr_to_sucroseegg = [corr_to_sucroseegg; out_corr];
        
        % look at signal correlated to transition (since behavior cam is 25 fps, transition is defined as 5 seconds
        transitions = diff(sucrose);
        [a b] = find(abs(transitions) > 0);
        transitionsabsval = zeros(1,length(transitions));
        transitionsabsval(b) = 1;
        [a b] = find(abs(transitions) < 0);
        transitionsabsval(b) = 1;
        transitionsabsval = smooth(transitionsabsval,10*25);
        [a b] = find((transitions) > 0);
        transitionssigned = zeros(1,length(transitions));
        transitionssigned(b) = 1;
        [a b] = find((transitions) < 0);
        transitionssigned(b) = -1;
        transitionssigned = smooth(transitionssigned,10*25);
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal,[0;transitionssigned], signaltime, movtime, time_interval_corr);
        corr_to_trans = [corr_to_trans; out_corr];
        
        [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal,[0;transitionsabsval], signaltime, movtime, time_interval_corr);
        corr_to_absvaltrans = [corr_to_absvaltrans; out_corr];
    end
    
    egg_structure(rec_index).time_interval = time_interval;
    egg_structure(rec_index).data_around_egg_sub = data_around_egg_sub;
    egg_structure(rec_index).data_around_egg_vel = repmat(data_around_egg_vel,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_vel_noave = repmat(data_around_egg_vel_noave,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_speed_2hz_hold = repmat(data_around_egg_speed_2hz_hold,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_body = repmat(data_around_egg_body,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_body_x = repmat(data_around_egg_body_x,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_body_y = repmat(data_around_egg_body_y,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_prob = repmat(data_around_egg_prob,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_prob_x = repmat(data_around_egg_prob_x,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_prob_y = repmat(data_around_egg_prob_y,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_path = repmat(data_around_egg_path,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_body_only = repmat(data_around_egg_body_only,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_body_only_x = repmat(data_around_egg_body_only_x,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_body_only_y = repmat(data_around_egg_body_only_y,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_angle = repmat(data_around_egg_angle,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_angle_only = repmat(data_around_egg_angle_only,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_PWM = repmat(data_around_egg_PWM,length(recording.tseries),1)';
   
    
    
    egg_structure(rec_index).data_around_egg = data_around_egg_dfoverf';
    egg_structure(rec_index).data_around_egg_mean = data_around_egg_dfoverf_mean';
    egg_structure(rec_index).data_around_egg_fo_is_mean = data_around_egg_dfoverf_fo_is_mean';

    
    
    egg_structure(rec_index).recording = all_recordings_list;
    egg_structure(rec_index).recording_index = all_recordings_list_index;
    egg_structure(rec_index).egg_time = all_eggs_listtime;
    egg_structure(rec_index).egg_index = all_eggs_listindex;
    egg_structure(rec_index).ROI_used = all_ROI_index;
    
    egg_structure(rec_index).cell_ID = all_cellID_list;
    egg_structure(rec_index).fly_ID = all_flyID_list;

    corr_structure(rec_index).fly_ID = fly_ID(rec_index);
    corr_structure(rec_index).cell_ID = cell_ID(rec_index);

    corr_structure(rec_index).recording = recordings_list(rec_index);
    corr_structure(rec_index).recording_index = rec_index;
    corr_structure(rec_index).ROI_used = ROI_num(rec_index);
    corr_structure(rec_index).time_base_corr = time_interval_corr;
    corr_structure(rec_index).corr_to_vel = nanmean(corr_to_vel,1)';
    corr_structure(rec_index).corr_to_vel_noave = nanmean(corr_to_vel_noave,1)';
    corr_structure(rec_index).corr_to_speed_2hz_hold = nanmean(corr_to_speed_2hz_hold,1)';
    corr_structure(rec_index).corr_to_body = nanmean(corr_to_body,1)';
    corr_structure(rec_index).corr_to_body_x = nanmean(corr_to_body_x,1)';
    corr_structure(rec_index).corr_to_body_y = nanmean(corr_to_body_y,1)';
    corr_structure(rec_index).corr_to_body_only = nanmean(corr_to_body_only,1)';
    corr_structure(rec_index).corr_to_body_only_x = nanmean(corr_to_body_only_x,1)';
    corr_structure(rec_index).corr_to_body_only_y = nanmean(corr_to_body_only_y,1)';
    corr_structure(rec_index).corr_to_prob = nanmean(corr_to_prob,1)';
    corr_structure(rec_index).corr_to_prob_x = nanmean(corr_to_prob_x,1)';
    corr_structure(rec_index).corr_to_prob_y = nanmean(corr_to_prob_y,1)';
    corr_structure(rec_index).corr_to_body_path = nanmean(corr_to_body_path,1)';
    corr_structure(rec_index).corr_to_angle = nanmean(corr_to_angle,1)';
    corr_structure(rec_index).corr_to_angle_only = nanmean(corr_to_angle_only,1)';
    corr_structure(rec_index).corr_to_sucrose = nanmean(corr_to_sucrose,1)';
    corr_structure(rec_index).corr_to_sucrose_t = nanmean(corr_to_sucrose_t,1)';
    corr_structure(rec_index).corr_to_trans = nanmean(corr_to_trans,1)';
    corr_structure(rec_index).corr_to_absvaltrans = nanmean(corr_to_absvaltrans,1)';
    corr_structure(rec_index).corr_to_egg = nanmean(corr_to_egg,1)';
    corr_structure(rec_index).corr_to_plainegg = nanmean(corr_to_plainegg,1)';
    corr_structure(rec_index).corr_to_sucroseegg = nanmean(corr_to_sucroseegg,1)';
    corr_structure(rec_index).recording = recordings_list(rec_index);
      
end
end


%% example plotting code for correlations and triggered averages
% 
% time_base_corr_to_body = corr_structure(1).time_base_corr(2:end-1);
% % %% example plotting code for correlations
% % if peak is neg, patch comes first. however...
% corr_to_body_x = [corr_structure(1:end).corr_to_body_x];
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_body_x')', nansem(corr_to_body_x')', time_base_corr_to_body', [.5 .5 .5]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_x,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen length, X only'); 
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3]);
% 
% 
% corr_to_body_x = [corr_structure(1:end).corr_to_angle];
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_body_x')', nansem(corr_to_body_x')', time_base_corr_to_body', [.5 .5 .5]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_x,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen angle'); 
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3]);
% 
% 
% corr_to_body_x = [corr_structure(1:end).corr_to_vel_noave];
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_body_x')', nansem(corr_to_body_x')', time_base_corr_to_body', [.5 .5 .5]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_x,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate velocity no ave'); 
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3]);
% 
% 
% corr_to_body_x = [corr_structure(1:end).corr_to_speed_2hz_hold];
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_body_x')', nansem(corr_to_body_x')', time_base_corr_to_body', [.5 .5 .5]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_x,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to speed 2hz hold'); 
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3]);
% 
% 
% corr_to_body_x = [corr_structure(1:end).corr_to_sucrose];
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_body_x')', nansem(corr_to_body_x')', time_base_corr_to_body', [.5 .5 .5]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_x,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to sucrose concentration'); 
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3]);
% 
% 
% corr_to_body_x=[corr_structure(1:end).corr_to_prob];
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_body_x')', nansem(corr_to_body_x')', time_base_corr_to_body', [.5 .5 .5]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_x,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to prob length'); 
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3]);
% % 
% % % 
% % % %% example plotting code for triggered averages
% % % 
% % % data_around_egg_vel = [egg_structure(1:end).data_around_egg_vel];
% % % data_around_egg_vel_noave = [egg_structure(1:end).data_around_egg_vel_noave];
% % % data_around_egg_speed_2hz_hold = [egg_structure(1:end).data_around_egg_speed_2hz_hold];
% % % data_around_egg_body = [egg_structure(1:end).data_around_egg_body];
% % % data_around_egg_body_x = [egg_structure(1:end).data_around_egg_body_x];
% % % data_around_egg_body_y = [egg_structure(1:end).data_around_egg_body_y];
% % % data_around_egg_prob = [egg_structure(1:end).data_around_egg_prob];
% % % data_around_egg_prob_x = [egg_structure(1:end).data_around_egg_prob_x];
% % % data_around_egg_prob_y = [egg_structure(1:end).data_around_egg_prob_y];
% % % data_around_egg_path = [egg_structure(1:end).data_around_egg_path];
% % % data_around_egg_body_only = [egg_structure(1:end).data_around_egg_body_only];
% % % data_around_egg_body_only_x = [egg_structure(1:end).data_around_egg_body_only_x];
% % % data_around_egg_body_only_y = [egg_structure(1:end).data_around_egg_body_only_y];
% % % data_around_egg_angle = [egg_structure(1:end).data_around_egg_angle];
% % % data_around_egg_angle_only = [egg_structure(1:end).data_around_egg_angle_only];
% % % data_around_egg = [egg_structure(1:end).data_around_egg];
% % % 
% % % % remove double plotting of behavior
% % % if(~isempty(data_around_egg))
% % %     [a b] = find(~isnan(data_around_egg(time_interval_mid,:)));
% % %     
% % %     data_around_egg_vel = data_around_egg_vel(:,b);  
% % %     data_around_egg_vel_noave = data_around_egg_vel_noave(:,b);
% % %     data_around_egg_speed_2hz_hold = data_around_egg_speed_2hz_hold(:,b);
% % %     data_around_egg_body = data_around_egg_body(:,b);
% % %     data_around_egg_body_x = data_around_egg_body_x(:,b);
% % %     data_around_egg_body_y = data_around_egg_body_y(:,b);  
% % %     data_around_egg_prob = data_around_egg_prob(:,b);
% % %     data_around_egg_prob_x = data_around_egg_prob_x(:,b);
% % %     data_around_egg_prob_y = data_around_egg_prob_y(:,b);   
% % %     data_around_egg_path = data_around_egg_path(:,b);
% % %     data_around_egg_body_only = data_around_egg_body_only(:,b);
% % %     data_around_egg_body_only_x = data_around_egg_body_only_x(:,b);
% % %     data_around_egg_body_only_y =data_around_egg_body_only_y(:,b);   
% % %     data_around_egg_angle = data_around_egg_angle(:,b);
% % %     data_around_egg_angle_only = data_around_egg_angle_only(:,b);
% % %     data_around_egg = data_around_egg(:,b);
% % % end
% % % 
% % % [pp,ppp] = size(data_around_egg_body);
% % % 
% % % figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% % % set(gca,'TickDir','out');
% % % if(~isempty(data_around_egg_body))
% % %     if(ppp > 1)
% % %         
% % %         patch_errorbar(nanmean(data_around_egg_vel_noave')', nansem(data_around_egg_vel_noave')', time_base_egg', [.5 .5 .5])
% % %     end
% % %     %plot(time_base_egg,data_around_egg_vel,'color',[.75,.75,.75]);
% % %     plot(time_base_egg,nanmean(data_around_egg_vel_noave,2),'-k');
% % % end
% % % ylabel('vel, mm/sec');
% % % xlabel('time (sec)'); %since peak is neg, patch comes first. however...
% % % axis manual;
% % % line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% % % 
% % % figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% % % set(gca,'TickDir','out');
% % % if(~isempty(data_around_egg_body))
% % %     if(ppp > 1)
% % %         
% % %         patch_errorbar(nanmean(data_around_egg_angle')', nansem(data_around_egg_angle')', time_base_egg', [.5 .5 .5])
% % %     end
% % %     %plot(time_base_egg,data_around_egg_vel,'color',[.75,.75,.75]);
% % %     plot(time_base_egg,nanmean(data_around_egg_angle,2),'-k');
% % % end
% % % ylabel('angle');
% % % xlabel('time (sec)'); %since peak is neg, patch comes first. however...
% % % axis manual;
% % % line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% % % 
% % % figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% % % set(gca,'TickDir','out');
% % % if(~isempty(data_around_egg_body))
% % %     if(ppp > 1)
% % %         
% % %         patch_errorbar(nanmean(data_around_egg_body_x')', nansem(data_around_egg_body_x')', time_base_egg', [.5 .5 .5])
% % %     end
% % %     %plot(time_base_egg,data_around_egg_vel,'color',[.75,.75,.75]);
% % %     plot(time_base_egg,nanmean(data_around_egg_body_x,2),'-k');
% % % end
% % % ylabel('body x');
% % % xlabel('time (sec)'); %since peak is neg, patch comes first. however...
% % % axis manual;
% % % line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% % % 
% % % 
% % % figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% % % set(gca,'TickDir','out');
% % % if(~isempty(data_around_egg_body))
% % %     if(ppp > 1)
% % %         
% % %         patch_errorbar(nanmean(data_around_egg_prob')', nansem(data_around_egg_prob')', time_base_egg', [.5 .5 .5])
% % %     end
% % %     %plot(time_base_egg,data_around_egg_vel,'color',[.75,.75,.75]);
% % %     plot(time_base_egg,nanmean(data_around_egg_prob,2),'-k');
% % % end
% % % ylabel('prob length');
% % % xlabel('time (sec)'); %since peak is neg, patch comes first. however...
% % % axis manual;
% % % line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% % % 
% % % 
% % 

