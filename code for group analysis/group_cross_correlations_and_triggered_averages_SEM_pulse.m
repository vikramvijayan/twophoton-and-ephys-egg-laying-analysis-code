function [corr_structure,pulse_structure] = group_cross_correlations_and_triggered_averages_SEM_pulse(recordings_list, ROI_num,  fly_ID, cell_ID,min_power, max_power, filter_out_pulse, backsub, pulse_duration_min,pulse_duration_max)

%% there is no corr_structure outputted

%% make sure to set the time interval for the calculations (triggered averages)
% time_interval = [-240:.1:240];
time_interval = [-1200:.1:1200];
time_interval_mid = (length(time_interval)-1)./2;

%% description of inputs
% set ROI num to -1 if you want tu use CH1 patch
% if patch is 1 then use CH1
% if patch is 0 then use CH1 but remove pulses (and 100 seconds after)

% set filter_out_pulse to 0 is you want all pulses

% set filter_out_pulse to 1 if you want pulses that didnt cause eggs (not
% with 60 sec before)

% set filter_out_pulse to 2 if you want chrimson pulse that activated egg (currently 60 sec before egg)


%% plotting time locked pulse average
for rec_index = 1:1:length(recordings_list)

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
        [recording] = replace_df_over_f_withbackgroundsubtracted_runningmeanexcludep(recording);
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
    
    data_around_pulse_vel = [];
    data_around_pulse_vel_noave = [];
    data_around_pulse_speed_2hz_hold = [];
    
    data_around_pulse_raw = [];
    data_around_pulse_raw_bck = [];
    
    data_around_pulse_dfoverf = [];
    data_around_pulse_dfoverf_mean = [];
    data_around_pulse_dfoverf_fo_is_mean = [];

    data_around_plain_pulse_dfoverf = [];
    data_around_sucrose_pulse_dfoverf = [];
    
    data_around_plain_pulse_dfoverf_mean = [];
    data_around_sucrose_pulse_dfoverf_mean = [];
    
    data_around_plain_pulse_dfoverf_fo_is_mean = [];
    data_around_sucrose_pulse_dfoverf_fo_is_mean = [];

    data_around_pulse_sub = [];
    data_around_plain_pulse_sub = [];
    data_around_sucrose_pulse_sub = [];
    
    data_around_pulse_temp = [];
    data_around_pulse_PWMlaseruW = [];
    data_around_pulse_puff = [];
    
    data_around_pulse_body = [];
    data_around_pulse_body_x = [];
    data_around_pulse_body_y = [];
    
    data_around_pulse_prob = [];
    data_around_pulse_prob_x = [];
    data_around_pulse_prob_y = [];
    
    data_around_pulse_path = [];
    data_around_pulse_body_only = [];
    data_around_pulse_body_only_x = [];
    data_around_pulse_body_only_y = [];
    
    data_around_pulse_angle = [];
    data_around_pulse_angle_only = [];
    
    data_around_pulse_egglaid = [];
    
    %% process more of the DLC movie (median normalize)
    
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
    
    % this is new code to get a 2hz speed
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
    
    % find all pulse times this is a bit complicted since the pulses have pulses (the stimulation is at a specific Hz)
    % this code will find the each beginning and end of a pulse train (5 pulses of 1
    % sec that are spaced 10 sec apart). That is it will find all 5 pulses
    
    temp_diffs = diff(recording.abf.PWMlaser);
    [a b] = find(temp_diffs > 2);
    a=a+1;
    temp_diffs2 = diff(a);
    [a1 b1] = find(temp_diffs2 > 10000);
    if(~isempty(a))
        tmp_times = [a(1); a(a1+1)];
    else
        tmp_times = [];
    end
    recording.abf.pulse_ON_times = [];
    
    for interate_pulses = 1:1:length(tmp_times)
        if(recording.abf.PWMlaser_uWpermm2(tmp_times(interate_pulses)) > min_power && recording.abf.PWMlaser_uWpermm2(tmp_times(interate_pulses)) < max_power)
            recording.abf.pulse_ON_times = [recording.abf.pulse_ON_times ;tmp_times(interate_pulses)] ;
        end
    end
    
    % pulse off is the last time it was ON
    temp_diffs = diff(recording.abf.PWMlaser);
    [a b] = find(temp_diffs < -2);
    temp_diffs2 = diff(a);
    [a1 b1] = find(temp_diffs2 > 10000);
    if(~isempty(a))
        tmp_times = [a(a1); a(end)];
    else
        tmp_times = [];
    end
    recording.abf.pulse_OFF_times = [];
    
    for interate_pulses = 1:1:length(tmp_times)
        if(recording.abf.PWMlaser_uWpermm2(tmp_times(interate_pulses)) > min_power && recording.abf.PWMlaser_uWpermm2(tmp_times(interate_pulses)) < max_power)
            recording.abf.pulse_OFF_times = [recording.abf.pulse_OFF_times ;tmp_times(interate_pulses)] ;
        end
    end
    
    tmp_ON= [];
    tmp_OFF = [];
    for interate_pulses = 1:1:length(recording.abf.pulse_OFF_times)
        if((recording.abf.pulse_OFF_times(interate_pulses)-recording.abf.pulse_ON_times(interate_pulses)) > pulse_duration_min*10000 && (recording.abf.pulse_OFF_times(interate_pulses)-recording.abf.pulse_ON_times(interate_pulses)) < pulse_duration_max*10000)
            tmp_ON = [tmp_ON; recording.abf.pulse_ON_times(interate_pulses)];
            tmp_OFF = [tmp_OFF; recording.abf.pulse_OFF_times(interate_pulses)];
        end
    end
    
    recording.abf.pulse_OFF_times = tmp_OFF;
    recording.abf.pulse_ON_times = tmp_ON;

    % convert abf times to movie index, either the first frame in movie after
    % pulse is ON, or last frame before it is OFF
    recording.movie1.pulse_ON_times_movie = [];
    pulse_ON_times = recording.abf.Time_s(recording.abf.pulse_ON_times);
    for iterate_pulse =1:1:length(pulse_ON_times)
        [a b] = find(recording.movie1.time_stamps-pulse_ON_times(iterate_pulse) >0,1,'first');
        recording.movie1.pulse_ON_times_movie(iterate_pulse) = a;
    end
    
    recording.movie1.pulse_OFF_times_movie = [];
    pulse_OFF_times = recording.abf.Time_s(recording.abf.pulse_OFF_times);
    for iterate_pulse =1:1:length(pulse_OFF_times)
        [a b] = find(recording.movie1.time_stamps-pulse_OFF_times(iterate_pulse) <0,1,'last');
        recording.movie1.pulse_OFF_times_movie(iterate_pulse) = a;
        
    end
    
    % this will get you 60 sec before the egg is laid as greater than 0,
    % any pulsees in this window will be causing egg-laying
    %recording.movie1.no_egg3_smooth = fliplr(causal_filter(2500, fliplr(recording.movie1.no_egg3)));
     recording.movie1.no_egg3_smooth = fliplr(causal_filter(25*60, fliplr(recording.movie1.no_egg3)));

    % this is code that filters the eggs depending on whether it was caused by a chrimson pulse or not
    if(filter_out_pulse == 2)
        [a b] = find(recording.movie1.no_egg3_smooth(recording.movie1.pulse_ON_times_movie) > 0);
        recording.movie1.pulse_ON_times_movie = recording.movie1.pulse_ON_times_movie(b);
        recording.movie1.pulse_OFF_times_movie = recording.movie1.pulse_OFF_times_movie(b);
    end
    
    if(filter_out_pulse == 1)
        [a b] = find(recording.movie1.no_egg3_smooth(recording.movie1.pulse_ON_times_movie) == 0);
        recording.movie1.pulse_ON_times_movie = recording.movie1.pulse_ON_times_movie(b);
        recording.movie1.pulse_OFF_times_movie = recording.movie1.pulse_OFF_times_movie(b);
    end
    
    [a b] = find(recording.movie1.sucrose(recording.movie1.pulse_ON_times_movie) == 0);
    plain_pulse = recording.movie1.pulse_ON_times_movie(b);
    plain_pulset = recording.movie1.time_stamps(recording.movie1.pulse_ON_times_movie(b));
    
    [a b] = find(recording.movie1.sucrose(recording.movie1.pulse_ON_times_movie) == 200 | recording.movie1.sucrose(recording.movie1.pulse_ON_times_movie) == 500 );
    sucrose_pulse = recording.movie1.pulse_ON_times_movie(b);
    sucrose_pulset = recording.movie1.time_stamps(recording.movie1.pulse_ON_times_movie(b));
    
    [a b] = find(recording.movie1.sucrose(recording.movie1.pulse_ON_times_movie) >= 0);
    all_pulse = recording.movie1.pulse_ON_times_movie(b);
    all_pulset = recording.movie1.time_stamps(recording.movie1.pulse_ON_times_movie(b));
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,all_pulset, time_interval);
    data_around_pulse_sub = [data_around_pulse_sub; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,plain_pulset, time_interval);
    data_around_plain_pulse_sub = [data_around_plain_pulse_sub; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,sucrose_pulset, time_interval);
    data_around_sucrose_pulse_sub = [data_around_sucrose_pulse_sub; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(recording.abf.Puff(1:1000:end),recording.abf.Time_s(1:1000:end),all_pulset, time_interval);
    data_around_pulse_puff = [data_around_pulse_puff; data_to_average_interp];
    
    %     if(~exist('recording.abf.PWMlaser_uWpermm2'))
    %         recording.abf.PWMlaser_uWpermm2 = zeros(length(recording.abf.PWMlaser),1);
    %     end
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(recording.abf.PWMlaser_uWpermm2(1:1000:end),recording.abf.Time_s(1:1000:end),all_pulset, time_interval);
    data_around_pulse_PWMlaseruW = [data_around_pulse_PWMlaseruW; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(recording.abf.Temp(1:1000:end),recording.abf.Time_s(1:1000:end),all_pulset, time_interval);
    data_around_pulse_temp = [data_around_pulse_temp; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(wheelvelocity,movtime,all_pulset, time_interval);
    data_around_pulse_vel = [data_around_pulse_vel; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(wheelvelocity_noave,movtime,all_pulset, time_interval);
    data_around_pulse_vel_noave = [data_around_pulse_vel_noave; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(speed_2hz_hold,movtime,all_pulset, time_interval);
    data_around_pulse_speed_2hz_hold = [data_around_pulse_speed_2hz_hold; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(problength,movtime,all_pulset, time_interval);
    data_around_pulse_prob = [data_around_pulse_prob; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(prob_x,movtime,all_pulset, time_interval);
    data_around_pulse_prob_x = [data_around_pulse_prob_x; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(prob_y,movtime,all_pulset, time_interval);
    data_around_pulse_prob_y = [data_around_pulse_prob_y; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(bodylength,movtime,all_pulset, time_interval);
    data_around_pulse_body = [data_around_pulse_body; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(bodylength_x,movtime,all_pulset, time_interval);
    data_around_pulse_body_x = [data_around_pulse_body_x; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(bodylength_y,movtime,all_pulset, time_interval);
    data_around_pulse_body_y = [data_around_pulse_body_y; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_path_length,movtime,all_pulset, time_interval);
    data_around_pulse_path = [data_around_pulse_path; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_length,movtime,all_pulset, time_interval);
    data_around_pulse_body_only = [data_around_pulse_body_only; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_length_x,movtime,all_pulset, time_interval);
    data_around_pulse_body_only_x = [data_around_pulse_body_only_x; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_length_y,movtime,all_pulset, time_interval);
    data_around_pulse_body_only_y = [data_around_pulse_body_only_y; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(bodyangle,movtime,all_pulset, time_interval);
    data_around_pulse_angle = [data_around_pulse_angle; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(causal_filter(25,recording.movie1.egglaid),movtime,all_pulset, time_interval);
    data_around_pulse_egglaid = [data_around_pulse_egglaid; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_angle,movtime,all_pulset, time_interval);
    data_around_pulse_angle_only = [data_around_pulse_angle_only; data_to_average_interp];
    

    data_around_sucrose_pulse_sub = repmat(data_around_sucrose_pulse_sub,length(recording.tseries),1);
    data_around_pulse_sub = repmat(data_around_pulse_sub,length(recording.tseries),1);
    data_around_plain_pulse_sub = repmat(data_around_plain_pulse_sub,length(recording.tseries),1);

    all_recordings_list = repmat(recordings_list(rec_index),length(recording.tseries),1);
    all_pulse_listindex = repmat(all_pulse,length(recording.tseries),1);
    all_pulse_listtime = repmat(all_pulset,length(recording.tseries),1);
    all_recordings_list_index = repmat(rec_index,length(recording.tseries),1);
    all_ROI_index = repmat(ROI_num(rec_index),length(recording.tseries),1);
    
     all_cellID_list = repmat(cell_ID(rec_index),length(recording.tseries),1);
    all_flyID_list = repmat(fly_ID(rec_index),length(recording.tseries),1);
    
    for m = 1:1:length(recording.tseries)
        if(patch == 0)
            signal = recording.tseries(m).df_over_f(ROI_num(rec_index),:);
            signal_over_mean = recording.tseries(m).df_over_f(ROI_num(rec_index),:)./nanmean(recording.tseries(m).df_over_f(ROI_num(rec_index),:));
            signal_fo_is_mean = recording.tseries(m).mean_in_ROI(ROI_num(rec_index),:)./nanmean(recording.tseries(m).mean_in_ROI(ROI_num(rec_index),:)) -1;

            signaltime = recording.tseries(m).Time_s;
            signal_raw = recording.tseries(m).mean_in_ROI(ROI_num(rec_index),:);
            signal_raw_putative_bck = recording.tseries(m).mean_in_ROI(end,:);  
        end
        
        %% for patching data, the signal changes
        if(patch ==1)
            signal    = recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
            signaltime = recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
            
            signal_raw    = recording.abf.CH1_patch_spikes(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
            signal_raw_putative_bck    = [];
            signal_over_mean = signal;
            signal_fo_is_mean = signal;

        end
        
        %% compute signal triggered signal averages
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,all_pulset, time_interval);
        data_around_pulse_dfoverf = [data_around_pulse_dfoverf; data_to_average_interp];
        
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal_over_mean,signaltime,all_pulset, time_interval);
        data_around_pulse_dfoverf_mean = [data_around_pulse_dfoverf_mean; data_to_average_interp];
        
                
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal_fo_is_mean,signaltime,all_pulset, time_interval);
        data_around_pulse_dfoverf_fo_is_mean = [data_around_pulse_dfoverf_fo_is_mean; data_to_average_interp];
        
        
        

        [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,plain_pulset, time_interval);
        data_around_plain_pulse_dfoverf = [data_around_plain_pulse_dfoverf; data_to_average_interp];
        
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,sucrose_pulset, time_interval);
        data_around_sucrose_pulse_dfoverf = [data_around_sucrose_pulse_dfoverf; data_to_average_interp];
        
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal_over_mean,signaltime,plain_pulset, time_interval);
        data_around_plain_pulse_dfoverf_mean = [data_around_plain_pulse_dfoverf_mean; data_to_average_interp];
        
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal_over_mean,signaltime,sucrose_pulset, time_interval);
        data_around_sucrose_pulse_dfoverf_mean = [data_around_sucrose_pulse_dfoverf_mean; data_to_average_interp];
        
           [time_base_to_return, data_to_average_interp]  = average_around_event(signal_fo_is_mean,signaltime,plain_pulset, time_interval);
        data_around_plain_pulse_dfoverf_fo_is_mean = [data_around_plain_pulse_dfoverf_fo_is_mean; data_to_average_interp];
        
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal_fo_is_mean,signaltime,sucrose_pulset, time_interval);
        data_around_sucrose_pulse_dfoverf_fo_is_mean = [data_around_sucrose_pulse_dfoverf_fo_is_mean; data_to_average_interp];
        
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal_raw,signaltime,all_pulset, time_interval);
        data_around_pulse_raw = [data_around_pulse_raw; data_to_average_interp];
        
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal_raw_putative_bck,signaltime,all_pulset, time_interval);
        data_around_pulse_raw_bck = [data_around_pulse_raw_bck; data_to_average_interp];
        
    end

    pulse_structure(rec_index).data_around_pulse_sub = data_around_pulse_sub;
    pulse_structure(rec_index).data_around_pulse_vel = repmat(data_around_pulse_vel,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_vel_noave = repmat(data_around_pulse_vel_noave,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_speed_2hz_hold = repmat(data_around_pulse_speed_2hz_hold,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_body = repmat(data_around_pulse_body,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_body_x = repmat(data_around_pulse_body_x,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_body_y = repmat(data_around_pulse_body_y,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_temp = repmat(data_around_pulse_temp,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_puff = repmat(data_around_pulse_puff,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_PWMlaseruW = repmat(data_around_pulse_PWMlaseruW,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_path = repmat(data_around_pulse_path,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_body_only = repmat(data_around_pulse_body_only,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_body_only_x = repmat(data_around_pulse_body_only_x,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_body_only_y = repmat(data_around_pulse_body_only_y,length(recording.tseries),1)';   
    pulse_structure(rec_index).data_around_pulse_angle = repmat(data_around_pulse_angle,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_angle_only = repmat(data_around_pulse_angle_only,length(recording.tseries),1)';

   

    
    
    pulse_structure(rec_index).data_around_pulse = data_around_pulse_dfoverf';
    pulse_structure(rec_index).data_around_pulse_mean = data_around_pulse_dfoverf_mean';
    pulse_structure(rec_index).data_around_pulse_fo_is_mean = data_around_pulse_dfoverf_fo_is_mean';

    pulse_structure(rec_index).data_around_pulse_plain = data_around_plain_pulse_dfoverf';
    pulse_structure(rec_index).data_around_pulse_plain_mean = data_around_plain_pulse_dfoverf_mean';
    pulse_structure(rec_index).data_around_pulse_plain_fo_is_mean = data_around_plain_pulse_dfoverf_fo_is_mean';

        
    pulse_structure(rec_index).data_around_pulse_sucrose = data_around_sucrose_pulse_dfoverf';
    pulse_structure(rec_index).data_around_pulse_sucrose_mean = data_around_sucrose_pulse_dfoverf_mean';
    pulse_structure(rec_index).data_around_pulse_sucrose_fo_is_mean = data_around_sucrose_pulse_dfoverf_fo_is_mean';

    pulse_structure(rec_index).data_around_pulse_raw = data_around_pulse_raw';
    pulse_structure(rec_index).data_around_pulse_raw_bck = data_around_pulse_raw_bck';   
    pulse_structure(rec_index).data_around_pulse_egglaid = repmat(data_around_pulse_egglaid,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_prob_x = repmat(data_around_pulse_prob_x,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_prob_y = repmat(data_around_pulse_prob_y,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_prob = repmat(data_around_pulse_prob,length(recording.tseries),1)';  
    pulse_structure(rec_index).recording = all_recordings_list;
    pulse_structure(rec_index).recording_index = all_recordings_list_index;
    pulse_structure(rec_index).pulse_time = all_pulse_listtime;
    pulse_structure(rec_index).pulse_index = all_pulse_listindex;
    pulse_structure(rec_index).ROI_used = all_ROI_index;
    
         pulse_structure(rec_index).cell_ID = all_cellID_list;
    pulse_structure(rec_index).fly_ID = all_flyID_list;

    corr_structure(rec_index).recording = recordings_list(rec_index);
    
    
end

end