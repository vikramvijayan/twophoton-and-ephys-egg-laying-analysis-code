function [corr_structure,pulse_structure] = group_cross_correlations_pulse_SEM_puff(recordings_list, ROI_num, min_power, max_power, filter_out_pulse, backsub, pulse_duration_min,pulse_duration_max)




% set ROI num to -1 if you want tu use CH1 patch
% if patch is 1 then use CH1
% if patch is 0 then use CH1 but remove pulses (and 100 seconds after)

% set filter_out_pulse to 0 is you want all pulses

% set filter_out_pulse to 1 if you want pulses that didnt cause eggs

% set filter_out_pulse to 2 if you owant chrimson pulse that activated egg (currently 100 sec before egg)


%%%%%%%%%%%
%% plotting correlations or time locked average

for rec_index = 1:1:length(recordings_list)
    
 
    
    %for rec_index = 10
    
    % it is often faster to use the stripped files without the associated
    % images
    modifiedStr = strrep([char(recordings_list(rec_index))], '.mat', '_stripped.mat');
    %modifiedStr = [char(recordings_list(rec_index))];
    recording = loadsinglerecording(modifiedStr);
    [recording] = processes_more_DLC_variables(recording);
    
       % comment this out, it is to have the puffs be used instead f PWM
     recording.abf.PWMlaser_uWpermm2 =  2.*recording.abf.Puff;
         recording.abf.PWMlaser =  2.*recording.abf.Puff;

    if(backsub == 1)
        [recording] = replace_df_over_f_withbackgroundsubtracted(recording);
    end
    if(ROI_num(rec_index) ==-1)
        patch = 1;
    else
        patch = 0;
    end
    if(patch == 1)
        recording.tseries = 1;
    end
    
    % recording.abf.CH1_patch_spikes_conv = 50000.*recording.abf.CH1_patch_spikes_conv_area_rect'./5;
    
    %%% this loop below is common code with plot cross correlations %%
    %% with one exception that is noted (a few lines added here, not in plot correlations)
    
    data_to_plain   = [];
    data_to_sucrose = [];
    
    data_to_sucrose_control = [];
    data_to_plain_control = [];
    time_base_trans = [];
    
    data_around_pulse_vel = [];
            data_around_pulse_vel_noave = [];

    data_around_pulse = [];
    data_around_plain_pulse = [];
    data_around_sucrose_pulse = [];
    time_base_pulse = [];
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
    
    data_around_pulse2 = [];
    data_around_plain_pulse2 = [];
    data_around_sucrose_pulse2 = [];
    time_base_pulse2 = [];
    data_around_pulse_sub2 = [];
    data_around_plain_pulse_sub2 = [];
    data_around_sucrose_pulse_sub2 = [];
    
    corr_to_vel = [];
    time_base_corr_to_vel = [];
    
    corr_to_body = [];
    corr_to_body_x = [];
    corr_to_body_y = [];
    
    corr_to_prob = [];
    corr_to_prob_x = [];
    corr_to_prob_y = [];
    
    
    time_base_corr_to_body = [];
    
    corr_to_angle = [];
    time_base_corr_to_angle = [];
    
    corr_to_body_only = [];
    corr_to_body_only_x = [];
    corr_to_body_only_y = [];
    
    corr_to_body_path = [];
    
    corr_to_angle_only = [];
    
    corr_to_sucrose = [];
    time_base_corr_to_sucrose = [];
    
    corr_to_sucrose_t = [];
    time_base_corr_to_sucrose_t = [];
    
    corr_to_trans = [];
    time_base_corr_to_trans = [];
    
    corr_to_absvaltrans = [];
    time_base_corr_to_absvaltrans = [];
    
    corr_to_pulse = [];
    time_base_corr_to_pulse = [];
    
    corr_to_plainpulse = [];
    time_base_corr_to_plainpulse = [];
    
    corr_to_sucrosepulse = [];
    time_base_corr_to_sucrosepulse = [];
    
    
    store_sucrose_mean = [];
    store_plain_mean   = [];
    store_sucrose_std  = [];
    store_plain_std    = [];
    
    movtime  = recording.movie1.time_stamps;
    sucrose = recording.movie1.sucrose;
    bodylength =  recording.movie1.abd_length;
    bodyangle =  recording.movie1.abd_angle;
    problength = recording.movie1.prob_length;
    abd_path_length = recording.movie1.abd_path_length;
    abd_only_length = recording.movie1.abd_only_length;
    abd_only_angle = recording.movie1.abd_only_angle;
wheelvelocity = (pi*7*-25).*smooth([0; diff(unwrap(recording.movie1.filtered_wheel))],25) ./ (2*pi); %mm/sec 7mm radium 
    wheelvelocity_noave = (pi*7*-25).*[0; diff(unwrap(recording.movie1.filtered_wheel))]./ (2*pi); %mm/sec 7mm radium

    % these are variables that were not processed earlier
    % process more of the DLC movie
    %[recording] = processes_more_DLC_variables(recording);
    abd_only_length_x = recording.movie1.abd_x_L3tip;
    bodylength_x =  recording.movie1.abd_x_neck_tip;
    abd_only_length_y = recording.movie1.abd_y_L3tip;
    bodylength_y =  recording.movie1.abd_y_neck_tip;
    
    prob_x =  recording.movie1.prob_x;
    prob_y =  recording.movie1.prob_y;
    
    
    % find all pulse times this is a bit complicted since the pulses have pulses (the stimulation is at a specific Hz)
    % this code will find the each beginning and end of a pulse train (5 pulses of 1
    % sec that are spaced 10 sec apart). That is it will find all 5 pulses
    
    temp_diffs = diff(recording.abf.PWMlaser);
    [a b] = find(temp_diffs > 2);
    a=a+1;
    temp_diffs2 = diff(a);
    [a1 b1] = find(temp_diffs2 > 1000);
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
    [a1 b1] = find(temp_diffs2 > 1000);
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

    
    recording.movie1.no_egg3 = causal_filter(25000, fliplr(recording.movie1.no_egg3));
    
    %convert abf times to movie index, either the first frame in movie after
    %pulse is ON, or last frame before it is OFF
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
    
    % this will get you 100 sec before the gg is laid as greater than 0,
    % any pulsees in this window will be causing egg-laying
    recording.movie1.no_egg3_smooth = flipud(causal_filter(25000, flipud(recording.movie1.no_egg3)));
    
    % this is code that filters the eggs depending onwhether it was caused by a chrimson pulse or not
    
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
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_sub = [data_around_pulse_sub; data_to_average_interp];
    time_base_pulse = time_base_to_return;
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,plain_pulset, [-240:.1:240]);
    data_around_plain_pulse_sub = [data_around_plain_pulse_sub; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,sucrose_pulset, [-240:.1:240]);
    data_around_sucrose_pulse_sub = [data_around_sucrose_pulse_sub; data_to_average_interp];
    
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(recording.abf.Puff(1:1000:end),recording.abf.Time_s(1:1000:end),all_pulset, [-240:.1:240]);
    data_around_pulse_puff = [data_around_pulse_puff; data_to_average_interp];
    
%     if(~exist('recording.abf.PWMlaser_uWpermm2'))
%         recording.abf.PWMlaser_uWpermm2 = zeros(length(recording.abf.PWMlaser),1);
%     end
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(recording.abf.PWMlaser_uWpermm2(1:1000:end),recording.abf.Time_s(1:1000:end),all_pulset, [-240:.1:240]);
    data_around_pulse_PWMlaseruW = [data_around_pulse_PWMlaseruW; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(recording.abf.Temp(1:1000:end),recording.abf.Time_s(1:1000:end),all_pulset, [-240:.1:240]);
    data_around_pulse_temp = [data_around_pulse_temp; data_to_average_interp];
    
   
    [time_base_to_return, data_to_average_interp]  = average_around_event(wheelvelocity,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_vel = [data_around_pulse_vel; data_to_average_interp];
    
      [time_base_to_return, data_to_average_interp]  = average_around_event(wheelvelocity_noave,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_vel_noave = [data_around_pulse_vel_noave; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(problength,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_prob = [data_around_pulse_prob; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(prob_x,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_prob_x = [data_around_pulse_prob_x; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(prob_y,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_prob_y = [data_around_pulse_prob_y; data_to_average_interp];
    
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(bodylength,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_body = [data_around_pulse_body; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(bodylength_x,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_body_x = [data_around_pulse_body_x; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(bodylength_y,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_body_y = [data_around_pulse_body_y; data_to_average_interp];
    
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_path_length,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_path = [data_around_pulse_path; data_to_average_interp];
    
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_length,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_body_only = [data_around_pulse_body_only; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_length_x,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_body_only_x = [data_around_pulse_body_only_x; data_to_average_interp];
    
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_length_y,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_body_only_y = [data_around_pulse_body_only_y; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(bodyangle,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_angle = [data_around_pulse_angle; data_to_average_interp];
    
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_angle,movtime,all_pulset, [-240:.1:240]);
    data_around_pulse_angle_only = [data_around_pulse_angle_only; data_to_average_interp];
    
    [time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,all_pulset, [-120:.1:120]);
    data_around_pulse_sub2 = [data_around_pulse_sub2; data_to_average_interp];
    time_base_pulse2 = time_base_to_return;
    [time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,plain_pulset, [-120:.1:120]);
    data_around_plain_pulse_sub2 = [data_around_plain_pulse_sub2; data_to_average_interp];
    [time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,sucrose_pulset, [-120:.1:120]);
    data_around_sucrose_pulse_sub2 = [data_around_sucrose_pulse_sub2; data_to_average_interp];
    
    data_around_sucrose_pulse_sub = repmat(data_around_sucrose_pulse_sub,length(recording.tseries),1);
    data_around_pulse_sub = repmat(data_around_pulse_sub,length(recording.tseries),1);
    data_around_plain_pulse_sub = repmat(data_around_plain_pulse_sub,length(recording.tseries),1);
    
    data_around_sucrose_pulse_sub2 = repmat(data_around_sucrose_pulse_sub2,length(recording.tseries),1);
    data_around_pulse_sub2 = repmat(data_around_pulse_sub2,length(recording.tseries),1);
    data_around_plain_pulse_sub2 = repmat(data_around_plain_pulse_sub2,length(recording.tseries),1);
    
    %% the lines below need to be added to newer versions of plots cross correlations
    all_recordings_list = repmat(recordings_list(rec_index),length(recording.tseries),1);
    all_pulse_listindex = repmat(all_pulse,length(recording.tseries),1);
    all_pulse_listtime = repmat(all_pulset,length(recording.tseries),1);
    all_recordings_list_index = repmat(rec_index,length(recording.tseries),1);
    all_ROI_index = repmat(ROI_num(rec_index),length(recording.tseries),1);
    %%
    
    
    for m = 1:1:length(recording.tseries)
        if(patch ==0)
            
            signal    = recording.tseries(m).df_over_f(ROI_num(rec_index),:);
            signaltime = recording.tseries(m).Time_s;
            
            
        end
        
        % for patching data, the signal changes
        if(patch ==1)
            signal    = recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
            signaltime = recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
            
            
        end
        
        
        wheelvelocity = (pi*7*-25).*smooth([0; diff(unwrap(recording.movie1.filtered_wheel))],25) ./ (2*pi); %mm/sec 7mm radium
    wheelvelocity_noave = (pi*7*-25).*[0; diff(unwrap(recording.movie1.filtered_wheel))]./ (2*pi); %mm/sec 7mm radium

        bodylength =  recording.movie1.abd_length;
        bodyangle =  recording.movie1.abd_angle;
        abd_path_length = recording.movie1.abd_path_length;
        abd_only_length = recording.movie1.abd_only_length;
        abd_only_angle = recording.movie1.abd_only_angle;
        abd_only_length_x = recording.movie1.abd_x_L3tip;
        bodylength_x =  recording.movie1.abd_x_neck_tip;
        abd_only_length_y = recording.movie1.abd_y_L3tip;
        bodylength_y =  recording.movie1.abd_y_neck_tip;
        prob_x =  recording.movie1.prob_x;
        prob_y =  recording.movie1.prob_y;
        problength = recording.movie1.prob_length;
        
        % look at signal at pulses
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,all_pulset, [-240:.1:240]);
        data_around_pulse = [data_around_pulse; data_to_average_interp];
        time_base_pulse = time_base_to_return;
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,plain_pulset, [-240:.1:240]);
        data_around_plain_pulse = [data_around_plain_pulse; data_to_average_interp];
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,sucrose_pulset, [-240:.1:240]);
        data_around_sucrose_pulse = [data_around_sucrose_pulse; data_to_average_interp];
        
        % look at signal at pulses
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,all_pulset, [-120:.1:120]);
        data_around_pulse2 = [data_around_pulse2; data_to_average_interp];
        time_base_pulse2 = time_base_to_return;
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,plain_pulset, [-120:.1:120]);
        data_around_plain_pulse2 = [data_around_plain_pulse2; data_to_average_interp];
        [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,sucrose_pulset, [-120:.1:120]);
        data_around_sucrose_pulse2 = [data_around_sucrose_pulse2; data_to_average_interp];
        %
        %         % look at signal correlated to velocity (smoothed)
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, wheelvelocity, signaltime, movtime, [-120:.1:120]);
        %         corr_to_vel = [corr_to_vel; out_corr];
        %         time_base_corr_to_vel = out_corr_shift_ofvec1;
        %
        %         % look at signal correlated to body position
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, problength, signaltime, movtime, [-120:.1:120]);
        %         corr_to_prob = [corr_to_prob; out_corr];
        %         time_base_corr_to_prob = out_corr_shift_ofvec1;
        %
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, prob_x, signaltime, movtime, [-120:.1:120]);
        %         corr_to_prob_x = [corr_to_prob_x; out_corr];
        %         time_base_corr_to_prob = out_corr_shift_ofvec1;
        %
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, prob_y, signaltime, movtime, [-120:.1:120]);
        %         corr_to_prob_y = [corr_to_prob_y; out_corr];
        %         time_base_corr_to_prob = out_corr_shift_ofvec1;
        %
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodylength, signaltime, movtime, [-120:.1:120]);
        %         corr_to_body = [corr_to_body; out_corr];
        %         time_base_corr_to_body = out_corr_shift_ofvec1;
        %
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodylength_x, signaltime, movtime, [-120:.1:120]);
        %         corr_to_body_x = [corr_to_body_x; out_corr];
        %         time_base_corr_to_body = out_corr_shift_ofvec1;
        %
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodylength_y, signaltime, movtime, [-120:.1:120]);
        %         corr_to_body_y = [corr_to_body_y; out_corr];
        %         time_base_corr_to_body = out_corr_shift_ofvec1;
        %
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodyangle, signaltime, movtime, [-120:.1:120]);
        %         corr_to_angle = [corr_to_angle; out_corr];
        %         time_base_corr_to_angle = out_corr_shift_ofvec1;
        %
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_path_length, signaltime, movtime, [-120:.1:120]);
        %         corr_to_body_path = [corr_to_body_path; out_corr];
        %
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_only_length, signaltime, movtime, [-120:.1:120]);
        %         corr_to_body_only = [corr_to_body_only; out_corr];
        %
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_only_length_x, signaltime, movtime, [-120:.1:120]);
        %         corr_to_body_only_x = [corr_to_body_only_x; out_corr];
        %
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_only_length_y, signaltime, movtime, [-120:.1:120]);
        %         corr_to_body_only_y = [corr_to_body_only_y; out_corr];
        %
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_only_angle, signaltime, movtime, [-120:.1:120]);
        %         corr_to_angle_only = [corr_to_angle_only; out_corr];
        %
        %         % look at signal correlated to substrate
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, sucrose, signaltime, movtime, [-120:.1:120]);
        %         corr_to_sucrose = [corr_to_sucrose; out_corr];
        %         time_base_corr_to_sucrose = out_corr_shift_ofvec1;
        %
        %         % look at signal correlated to substrate threshold
        %         suc_thr = sucrose;
        %         suc_thr(suc_thr>0) = 1;
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, suc_thr, signaltime, movtime, [-120:.1:120]);
        %         corr_to_sucrose_t = [corr_to_sucrose_t; out_corr];
        %         time_base_corr_to_sucrose_t = out_corr_shift_ofvec1;
        %
        %         % look at signal correlated to pulse [using 2 seconds after pulse
        %         % starts -- for a 1 sec pulse, this includes aother second)
        %         tmp_ar = zeros(1,length(movtime));
        %         tmp_ar(all_pulse:(all_pulse+50)) = 1;
        %         %tmp_ar = smooth(tmp_ar,25*10);
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, [-120:.1:120]);
        %         corr_to_pulse = [corr_to_pulse; out_corr];
        %         time_base_corr_to_pulse = out_corr_shift_ofvec1;
        %
        %         % look at signal correlated to plain pulse
        %         tmp_ar = zeros(1,length(movtime));
        %         tmp_ar(plain_pulse:(plain_pulse+50)) = 1;
        %         %tmp_ar = smooth(tmp_ar,25*10);
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, [-120:.1:120]);
        %         corr_to_plainpulse = [corr_to_plainpulse; out_corr];
        %         time_base_corr_to_plainpulse = out_corr_shift_ofvec1;
        %
        %         % look at signal correlated to sucrose pulse
        %         tmp_ar = zeros(1,length(movtime));
        %         tmp_ar(sucrose_pulse:(sucrose_pulse+50)) = 1;
        %         %tmp_ar = smooth(tmp_ar,25*10);
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, [-120:.1:120]);
        %         corr_to_sucrosepulse = [corr_to_sucrosepulse; out_corr];
        %         time_base_corr_to_sucrosepulse = out_corr_shift_ofvec1;
        %
        %         % look at signal correlated to transition (since behavior cam is 25 fps, transition is defined as 5 seconds
        %         transitions = diff(sucrose);
        %         [a b] = find(abs(transitions) > 0);
        %         transitionsabsval = zeros(1,length(transitions));
        %         transitionsabsval(b) = 1;
        %         [a b] = find(abs(transitions) < 0);
        %         transitionsabsval(b) = 1;
        %         transitionsabsval = smooth(transitionsabsval,10*25);
        %         [a b] = find((transitions) > 0);
        %         transitionssigned = zeros(1,length(transitions));
        %         transitionssigned(b) = 1;
        %         [a b] = find((transitions) < 0);
        %         transitionssigned(b) = -1;
        %         transitionssigned = smooth(transitionssigned,10*25);
        %
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal,[0;transitionssigned], signaltime, movtime, [-120:.1:120]);
        %         corr_to_trans = [corr_to_trans; out_corr];
        %         time_base_corr_to_trans = out_corr_shift_ofvec1;
        %
        %         [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal,[0;transitionsabsval], signaltime, movtime, [-120:.1:120]);
        %         corr_to_absvaltrans = [corr_to_absvaltrans; out_corr];
        %         time_base_corr_to_absvaltrans = out_corr_shift_ofvec1;
    end
    %%%% the loop above can be replaced with new vweaions of the
    %%%% plot_cross_correlation function. It is ahared
    
    pulse_structure(rec_index).data_around_pulse_sub = data_around_pulse_sub;
    pulse_structure(rec_index).data_around_pulse_vel = repmat(data_around_pulse_vel,length(recording.tseries),1)';
        pulse_structure(rec_index).data_around_pulse_vel_noave = repmat(data_around_pulse_vel_noave,length(recording.tseries),1)';

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
    pulse_structure(rec_index).data_around_pulse = data_around_pulse';
    
    pulse_structure(rec_index).data_around_pulse_prob_x = repmat(data_around_pulse_prob_x,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_prob_y = repmat(data_around_pulse_prob_y,length(recording.tseries),1)';
    pulse_structure(rec_index).data_around_pulse_prob = repmat(data_around_pulse_prob,length(recording.tseries),1)';
    
    pulse_structure(rec_index).recording = all_recordings_list;
    pulse_structure(rec_index).recording_index = all_recordings_list_index;
    pulse_structure(rec_index).pulse_time = all_pulse_listtime;
    pulse_structure(rec_index).pulse_index = all_pulse_listindex;
    pulse_structure(rec_index).ROI_used = all_ROI_index;
    
    %
    %
    %
    corr_structure(rec_index).recording = recordings_list(rec_index);
    %     corr_structure(rec_index).recording_index = rec_index;
    %     corr_structure(rec_index).ROI_used = ROI_num(rec_index);
    %
    %
    %     corr_structure(rec_index).corr_to_vel = nanmean(corr_to_vel,1)';
    %     corr_structure(rec_index).time_base_corr_to_vel = time_base_corr_to_vel;
    %
    %
    %     % figure; hold on;  title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    %     % set(gca,'TickDir','out');
    %     % plot(time_base_corr_to_vel,corr_to_vel,'color',[.75,.75,.75]);
    %     % plot(time_base_corr_to_vel,nanmean(corr_to_vel),'-k');
    %     % ylabel('corr coeff');
    %     % xlabel('shift signal, correlate to smoothed velocity'); %since peak is neg, patch comes first. however...
    %     % axis manual;
    %     % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    %
    %     corr_structure(rec_index).corr_to_body = nanmean(corr_to_body,1)';
    %     corr_structure(rec_index).time_base_corr_to_body = time_base_corr_to_body;
    %
    %     corr_structure(rec_index).corr_to_body_x = nanmean(corr_to_body_x,1)';
    %     corr_structure(rec_index).corr_to_body_y = nanmean(corr_to_body_y,1)';
    %
    %     % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    %     % set(gca,'TickDir','out');
    %     % plot(time_base_corr_to_body,corr_to_body,'color',[.75,.75,.75]);
    %     % plot(time_base_corr_to_body,nanmean(corr_to_body),'-k');
    %     % ylabel('corr coeff');
    %     % xlabel('shift signal, correlate to abdomen length'); %since peak is neg, patch comes first. however...
    %     % axis manual;
    %     % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    %
    %     corr_structure(rec_index).corr_to_body_only = nanmean(corr_to_body_only,1)';
    %     corr_structure(rec_index).corr_to_body_only_x = nanmean(corr_to_body_only_x,1)';
    %     corr_structure(rec_index).corr_to_body_only_y = nanmean(corr_to_body_only_y,1)';
    %
    %     % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    %     % set(gca,'TickDir','out');
    %     % plot(time_base_corr_to_body,corr_to_body_only,'color',[.75,.75,.75]);
    %     % plot(time_base_corr_to_body,nanmean(corr_to_body_only),'-k');
    %     % ylabel('corr coeff');
    %     % xlabel('shift signal, correlate to abdomen length L3 to tip'); %since peak is neg, patch comes first. however...
    %     % axis manual;
    %     % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    %
    %     corr_structure(rec_index).corr_to_body_path = nanmean(corr_to_body_path,1)';
    %     % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    %     % set(gca,'TickDir','out');
    %     % plot(time_base_corr_to_body,corr_to_body_path,'color',[.75,.75,.75]);
    %     % plot(time_base_corr_to_body,nanmean(corr_to_body_path),'-k');
    %     % ylabel('corr coeff');
    %     % xlabel('shift signal, correlate to abdomen path length L3 to tip'); %since peak is neg, patch comes first. however...
    %     % axis manual;
    %     % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    %
    %     corr_structure(rec_index).corr_to_angle = nanmean(corr_to_angle,1)';
    %     corr_structure(rec_index).time_base_corr_to_angle = time_base_corr_to_angle;
    %
    %     % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    %     % set(gca,'TickDir','out');
    %     % plot(time_base_corr_to_angle,corr_to_angle,'color',[.75,.75,.75]);
    %     % plot(time_base_corr_to_angle,nanmean(corr_to_angle),'-k');
    %     % ylabel('corr coeff');
    %     % xlabel('shift signal, correlate to abdomen angle'); %since peak is neg, patch comes first. however...
    %     % axis manual;
    %     % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    %
    %     corr_structure(rec_index).corr_to_angle_only = nanmean(corr_to_angle_only,1)';
    %     % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    %     % set(gca,'TickDir','out');
    %     % plot(time_base_corr_to_angle,corr_to_angle_only,'color',[.75,.75,.75]);
    %     % plot(time_base_corr_to_angle,nanmean(corr_to_angle_only),'-k');
    %     % ylabel('corr coeff');
    %     % xlabel('shift signal, correlate to abdomen angle L3 to tip'); %since peak is neg, patch comes first. however...
    %     % axis manual;
    %     % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    %
    %     corr_structure(rec_index).corr_to_sucrose = nanmean(corr_to_sucrose,1)';
    %     corr_structure(rec_index).time_base_corr_to_sucrose = time_base_corr_to_sucrose;
    %
    %     % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    %     % set(gca,'TickDir','out');
    %     % plot(time_base_corr_to_sucrose,corr_to_sucrose,'color',[.75,.75,.75]);
    %     % plot(time_base_corr_to_sucrose,nanmean(corr_to_sucrose),'-k');
    %     % ylabel('corr coeff');
    %     % xlabel('shift signal, correlate to sucrose concentration'); %since peak is neg, patch comes first. however...
    %     % axis manual;
    %     % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    %
    %     corr_structure(rec_index).corr_to_sucrose_t = nanmean(corr_to_sucrose_t,1)';
    %     corr_structure(rec_index).time_base_corr_to_sucrose_t = time_base_corr_to_sucrose_t;
    %
    %
    %     % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    %     % set(gca,'TickDir','out');
    %     % plot(time_base_corr_to_sucrose_t,corr_to_sucrose_t,'color',[.75,.75,.75]);
    %     % plot(time_base_corr_to_sucrose_t,nanmean(corr_to_sucrose_t),'-k');
    %     % ylabel('corr coeff');
    %     % xlabel('shift signal, correlate to sucrose thresholded concentration'); %since peak is neg, patch comes first. however...
    %     % axis manual;
    %     % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    %
    %     corr_structure(rec_index).corr_to_trans = nanmean(corr_to_trans,1)';
    %     corr_structure(rec_index).time_base_corr_to_trans = time_base_corr_to_trans;
    %
    %
    %     % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    %     % set(gca,'TickDir','out');
    %     % plot(time_base_corr_to_trans,corr_to_trans,'color',[.75,.75,.75]);
    %     % plot(time_base_corr_to_trans,nanmean(corr_to_trans),'-k');
    %     % ylabel('corr coeff');
    %     % xlabel('shift signal, correlate to signed transitions (10 sec around trans is 1 or -1)'); %since peak is neg, patch comes first. however...
    %     % axis manual;
    %     % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    %
    %     corr_structure(rec_index).corr_to_absvaltrans = nanmean(corr_to_absvaltrans,1)';
    %     corr_structure(rec_index).time_base_corr_to_absvaltrans = time_base_corr_to_absvaltrans;
    %
    %
    %     % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    %     % set(gca,'TickDir','out');
    %     % plot(time_base_corr_to_absvaltrans,corr_to_absvaltrans,'color',[.75,.75,.75]);
    %     % plot(time_base_corr_to_absvaltrans,nanmean(corr_to_absvaltrans),'-k');
    %     % ylabel('corr coeff');
    %     % xlabel('shift signal, correlate to abs val transitions (10 sec around trans is )'); %since peak is neg, patch comes first. however...
    %     % axis manual;
    %     % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    %
    %     corr_structure(rec_index).corr_to_pulse = nanmean(corr_to_pulse,1)';
    %     corr_structure(rec_index).time_base_corr_to_pulse = time_base_corr_to_pulse;
    %
    %     % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    %     % set(gca,'TickDir','out');
    %     % plot(time_base_corr_to_pulse,corr_to_pulse,'color',[.75,.75,.75]);
    %     % plot(time_base_corr_to_pulse,nanmean(corr_to_pulse),'-k');
    %     % ylabel('corr coeff');
    %     % xlabel('shift signal, correlate to pulses'); %since peak is neg, patch comes first. however...
    %     % axis manual;
    %     % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    %
    %     corr_structure(rec_index).corr_to_plainpulse = nanmean(corr_to_plainpulse,1)';
    %     corr_structure(rec_index).time_base_corr_to_plainpulse = time_base_corr_to_plainpulse;
    %
    %     % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    %     % set(gca,'TickDir','out');
    %     % plot(time_base_corr_to_plainpulse,corr_to_plainpulse,'color',[.75,.75,.75]);
    %     % plot(time_base_corr_to_plainpulse,nanmean(corr_to_plainpulse),'-k');
    %     % ylabel('corr coeff');
    %     % xlabel('shift signal, correlate to plain pulses'); %since peak is neg, patch comes first. however...
    %     % axis manual;
    %     % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    %
    %     corr_structure(rec_index).corr_to_sucrosepulse = nanmean(corr_to_sucrosepulse,1)';
    %     corr_structure(rec_index).time_base_corr_to_sucrosepulse = time_base_corr_to_sucrosepulse;
    %
    %     % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    %     % set(gca,'TickDir','out');
    %     % plot(time_base_corr_to_sucrosepulse,corr_to_sucrosepulse,'color',[.75,.75,.75]);
    %     % plot(time_base_corr_to_sucrosepulse,nanmean(corr_to_sucrosepulse),'-k');
    %     % ylabel('corr coeff');
    %     % xlabel('shift signal, correlate to sucrose pulses'); %since peak is neg, patch comes first. however...
    %     % axis manual;
    %     % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
end

%
% corr_to_vel = [corr_structure(1:end).corr_to_vel];
%
% figure; hold on;  %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_vel')', nansem(corr_to_vel')', time_base_corr_to_vel', [.5 .5 .5])
% %plot(time_base_corr_to_vel,corr_to_vel,'color',[.75,.75,.75]);
% plot(time_base_corr_to_vel,nanmean(corr_to_vel,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to smoothed velocity'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3])
%
%
% corr_to_body = [corr_structure(1:end).corr_to_body];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_body')', nansem(corr_to_body')', time_base_corr_to_body', [.5 .5 .5])
% %plot(time_base_corr_to_body,corr_to_body,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen length'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3])
%
% corr_to_body_x = [corr_structure(1:end).corr_to_body_x];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_body_x')', nansem(corr_to_body_x')', time_base_corr_to_body', [.5 .5 .5])
% %plot(time_base_corr_to_body,corr_to_body_x,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_x,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen length, X only'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3])
%
%
% corr_to_body_y = [corr_structure(1:end).corr_to_body_y];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_body_y')', nansem(corr_to_body_y')', time_base_corr_to_body', [.5 .5 .5])
% %plot(time_base_corr_to_body,corr_to_body_y,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_y,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen length, Y only'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3])
%
% corr_to_body_only = [corr_structure(1:end).corr_to_body_only];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_body_only')', nansem(corr_to_body_only')', time_base_corr_to_body', [.5 .5 .5])
% %plot(time_base_corr_to_body,corr_to_body_only,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_only,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen length L3 to tip'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3])
%
%
%
% corr_to_body_only_x = [corr_structure(1:end).corr_to_body_only_x];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_body_only_x')', nansem(corr_to_body_only_x')', time_base_corr_to_body', [.5 .5 .5])
% %plot(time_base_corr_to_body,corr_to_body_only_x,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_only_x,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen length L3 to tip, X only'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3])
%
%
% corr_to_body_only_y = [corr_structure(1:end).corr_to_body_only_y];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_body_only_y')', nansem(corr_to_body_only_y')', time_base_corr_to_body', [.5 .5 .5])
% %plot(time_base_corr_to_body,corr_to_body_only_y,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_only_y,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen length L3 to tip, Y only'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3])
%
%
% corr_to_body_path = [corr_structure(1:end).corr_to_body_path];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_body_path')', nansem(corr_to_body_path')', time_base_corr_to_body', [.5 .5 .5])
% %plot(time_base_corr_to_body,corr_to_body_path,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_path,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen path length L3 to tip'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3])
%
% corr_to_angle = [corr_structure(1:end).corr_to_angle];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_angle')', nansem(corr_to_angle')', time_base_corr_to_angle', [.5 .5 .5])
% %plot(time_base_corr_to_angle,corr_to_angle,'color',[.75,.75,.75]);
% plot(time_base_corr_to_angle,nanmean(corr_to_angle,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen angle'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3])
%
% corr_to_angle_only = [corr_structure(1:end).corr_to_angle_only];
%
% % figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% % set(gca,'TickDir','out');
% % plot(time_base_corr_to_angle,corr_to_angle_only,'color',[.75,.75,.75]);
% % plot(time_base_corr_to_angle,nanmean(corr_to_angle_only),'-k');
% % ylabel('corr coeff');
% % xlabel('shift signal, correlate to abdomen angle L3 to tip'); %since peak is neg, patch comes first. however...
% % axis manual;
% %line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% %
%
% corr_to_sucrose = [corr_structure(1:end).corr_to_sucrose];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_sucrose')', nansem(corr_to_sucrose')', time_base_corr_to_sucrose', [.5 .5 .5])
% %plot(time_base_corr_to_sucrose,corr_to_sucrose,'color',[.75,.75,.75]);
% plot(time_base_corr_to_sucrose,nanmean(corr_to_sucrose,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to sucrose concentration'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3])
%
% corr_to_sucrose_t = [corr_structure(1:end).corr_to_sucrose_t];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_sucrose_t')', nansem(corr_to_sucrose_t')', time_base_corr_to_sucrose_t', [.5 .5 .5])
% %plot(time_base_corr_to_sucrose_t,corr_to_sucrose_t,'color',[.75,.75,.75]);
% plot(time_base_corr_to_sucrose_t,nanmean(corr_to_sucrose_t,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to sucrose thresholded concentration'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3])
%
% corr_to_trans = [corr_structure(1:end).corr_to_trans];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_trans')', nansem(corr_to_trans')', time_base_corr_to_trans', [.5 .5 .5])
% %plot(time_base_corr_to_trans,corr_to_trans,'color',[.75,.75,.75]);
% plot(time_base_corr_to_trans,nanmean(corr_to_trans,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to signed transitions (10 sec around trans is 1 or -1)'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3])
%
% corr_to_absvaltrans = [corr_structure(1:end).corr_to_absvaltrans];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_absvaltrans')', nansem(corr_to_absvaltrans')', time_base_corr_to_absvaltrans', [.5 .5 .5])
% %plot(time_base_corr_to_absvaltrans,corr_to_absvaltrans,'color',[.75,.75,.75]);
% plot(time_base_corr_to_absvaltrans,nanmean(corr_to_absvaltrans,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abs val transitions (10 sec around trans is )'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.3,.3])
%
% corr_to_pulse = [corr_structure(1:end).corr_to_pulse];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_pulse')', nansem(corr_to_pulse')', time_base_corr_to_pulse', [.5 .5 .5])
% %plot(time_base_corr_to_pulse,corr_to_pulse,'color',[.75,.75,.75]);
% plot(time_base_corr_to_pulse,nanmean(corr_to_pulse,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to pulses'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.4,.4])
%
% corr_to_plainpulse = [corr_structure(1:end).corr_to_plainpulse];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_plainpulse')', nansem(corr_to_plainpulse')', time_base_corr_to_plainpulse', [.5 .5 .5])
% %plot(time_base_corr_to_plainpulse,corr_to_plainpulse,'color',[.75,.75,.75]);
% plot(time_base_corr_to_plainpulse,nanmean(corr_to_plainpulse,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to plain pulses'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.4,.4])
%
%
% corr_to_sucrosepulse = [corr_structure(1:end).corr_to_sucrosepulse];
%
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% patch_errorbar(nanmean(corr_to_sucrosepulse')', nansem(corr_to_sucrosepulse')', time_base_corr_to_sucrosepulse', [.5 .5 .5])
% %plot(time_base_corr_to_sucrosepulse,corr_to_sucrosepulse,'color',[.75,.75,.75]);
% plot(time_base_corr_to_sucrosepulse,nanmean(corr_to_sucrosepulse,2),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to sucrose pulses'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% set(gca,'ylim',[-.4,.4])
%

% plotting pulse triggered averages
data_around_pulse_vel = [pulse_structure(1:end).data_around_pulse_vel];
data_around_pulse_vel_noave = [pulse_structure(1:end).data_around_pulse_vel_noave];


data_around_pulse_body = [pulse_structure(1:end).data_around_pulse_body];
data_around_pulse_body_x = [pulse_structure(1:end).data_around_pulse_body_x];
data_around_pulse_body_y = [pulse_structure(1:end).data_around_pulse_body_y];

data_around_pulse_path = [pulse_structure(1:end).data_around_pulse_path];
data_around_pulse_body_only = [pulse_structure(1:end).data_around_pulse_body_only];
data_around_pulse_body_only_x = [pulse_structure(1:end).data_around_pulse_body_only_x];
data_around_pulse_body_only_y = [pulse_structure(1:end).data_around_pulse_body_only_y];

data_around_pulse_angle = [pulse_structure(1:end).data_around_pulse_angle];
data_around_pulse_angle_only = [pulse_structure(1:end).data_around_pulse_angle_only];
data_around_pulse = [pulse_structure(1:end).data_around_pulse];

data_around_pulse_prob = [pulse_structure(1:end).data_around_pulse_prob];
data_around_pulse_prob_x = [pulse_structure(1:end).data_around_pulse_prob_x];
data_around_pulse_prob_y = [pulse_structure(1:end).data_around_pulse_prob_y];

data_around_pulse_temp = [pulse_structure(1:end).data_around_pulse_temp];
data_around_pulse_PWMlaseruW = [pulse_structure(1:end).data_around_pulse_PWMlaseruW];
data_around_pulse_puff = [pulse_structure(1:end).data_around_pulse_puff];

% remove double plotting of behavior
if(~isempty(data_around_pulse))
[a b] = find(~isnan(data_around_pulse(2400,:)));

data_around_pulse_vel = data_around_pulse_vel(:,b);
data_around_pulse_vel_noave = data_around_pulse_vel_noave(:,b);

data_around_pulse_body = data_around_pulse_body(:,b);
data_around_pulse_body_x = data_around_pulse_body_x(:,b);
data_around_pulse_body_y = data_around_pulse_body_y(:,b);
data_around_pulse_path = data_around_pulse_path(:,b);
data_around_pulse_body_only = data_around_pulse_body_only(:,b);
data_around_pulse_body_only_x = data_around_pulse_body_only_x(:,b);
data_around_pulse_body_only_y = data_around_pulse_body_only_y(:,b);
data_around_pulse_angle = data_around_pulse_angle(:,b);
data_around_pulse_angle_only = data_around_pulse_angle_only(:,b);
data_around_pulse = data_around_pulse(:,b);
data_around_pulse_prob = data_around_pulse_prob(:,b);
data_around_pulse_prob_x = data_around_pulse_prob_x(:,b);
data_around_pulse_prob_y = data_around_pulse_prob_y(:,b);
data_around_pulse_temp = data_around_pulse_temp(:,b);
data_around_pulse_PWMlaseruW = data_around_pulse_PWMlaseruW(:,b);
data_around_pulse_puff = data_around_pulse_puff(:,b);
end

[pp,ppp] = size(data_around_pulse_body);


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
    if(ppp > 1)
    patch_errorbar(nanmean(data_around_pulse_vel_noave')', nansem(data_around_pulse_vel_noave')', time_base_pulse', [.5 .5 .5]);
    end
    %plot(time_base_pulse,data_around_pulse_vel,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_vel_noave,2),'-k');
end
ylabel('vel, mm/sec');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_vel')', nansem(data_around_pulse_vel')', time_base_pulse', [.5 .5 .5])
        end
    %plot(time_base_pulse,data_around_pulse_vel,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_vel,2),'-k');
end
ylabel('vel, mm/sec');
xlabel('time (sec) 1s smooth'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);

% d = [];
% for i =1:1:size(data_around_pulse_vel,2)
%     d(:,i) = causal_filter(200,data_around_pulse_vel(:,i));
% end
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% if(~isempty(data_around_pulse_body))
%     patch_errorbar(nanmean(d')', nansem(d')', time_base_pulse', [.5 .5 .5])
%     %plot(time_base_pulse,data_around_pulse_vel,'color',[.75,.75,.75]);
%     plot(time_base_pulse,nanmean(d,2),'-k');
% end
% ylabel('vel, mm/sec (smth 1s), then 20s causal smth');
% xlabel('time (sec)'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% snapnow;
% set(gca,'xlim',[-60,60]);



figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_body')', nansem(data_around_pulse_body')', time_base_pulse', [.5 .5 .5])
        end
    % plot(time_base_pulse,data_around_pulse_body,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_body,2),'-k');
end
ylabel('neck to tip length');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_body_x')', nansem(data_around_pulse_body_x')', time_base_pulse', [.5 .5 .5])
        end
    % plot(time_base_pulse,data_around_pulse_body_x,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_body_x,2),'-k');
end
ylabel('neck to tip length, X only');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_body_y')', nansem(data_around_pulse_body_y')', time_base_pulse', [.5 .5 .5])
        end
    %  plot(time_base_pulse,data_around_pulse_body_y,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_body_y,2),'-k');
end
ylabel('neck to tip length, Y only');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_path')', nansem(data_around_pulse_path')', time_base_pulse', [.5 .5 .5])
        end
    %  plot(time_base_pulse,data_around_pulse_path,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_path,2),'-k');
end
ylabel('neck to tip path length');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_body_only')', nansem(data_around_pulse_body_only')', time_base_pulse', [.5 .5 .5])
        end
    % plot(time_base_pulse,data_around_pulse_body_only,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_body_only,2),'-k');
end
ylabel('L3 to tip length');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_body_only_x')', nansem(data_around_pulse_body_only_x')', time_base_pulse', [.5 .5 .5])
        end
    %   plot(time_base_pulse,data_around_pulse_body_only_x,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_body_only_x,2),'-k');
end
ylabel('L3 to tip length, X only');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_body_only_y')', nansem(data_around_pulse_body_only_y')', time_base_pulse', [.5 .5 .5])
        end
    % plot(time_base_pulse,data_around_pulse_body_only_y,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_body_only_y,2),'-k');
end
ylabel('L3 to tip length, Y only');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_angle')', nansem(data_around_pulse_angle')', time_base_pulse', [.5 .5 .5])
        end
    %  plot(time_base_pulse,data_around_pulse_angle,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_angle,2),'-k');
end
ylabel('neck to tip angle');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);



figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_prob')', nansem(data_around_pulse_prob')', time_base_pulse', [.5 .5 .5])
        end
    %  plot(time_base_pulse,data_around_pulse_angle,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_prob,2),'-k');
end
ylabel('proboscis length');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_prob_x')', nansem(data_around_pulse_prob_x')', time_base_pulse', [.5 .5 .5])
        end
    %  plot(time_base_pulse,data_around_pulse_angle,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_prob_x,2),'-k');
end
ylabel('proboscis length x');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_prob_y')', nansem(data_around_pulse_prob_y')', time_base_pulse', [.5 .5 .5])
        end
    %  plot(time_base_pulse,data_around_pulse_angle,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_prob_y,2),'-k');
end
ylabel('proboscis length y');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);



figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse')', nansem(data_around_pulse')', time_base_pulse', [.5 .5 .5])
        end
    % plot(time_base_pulse,data_around_pulse,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse,2),'-k');
end
ylabel('signal around pulse');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
if(~patch)
    set(gca,'ylim',[-.5,1.5])
    %set(gca,'ylim',[-.5,3.5])
end
snapnow;
set(gca,'xlim',[-60,60]);




figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_PWMlaseruW')', nansem(data_around_pulse_PWMlaseruW')', time_base_pulse', [.5 .5 .5])
        end
    %  plot(time_base_pulse,data_around_pulse_angle,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_PWMlaseruW,2),'-k');
end
ylabel('PWM laser uW');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);



figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_temp')', nansem(data_around_pulse_temp')', time_base_pulse', [.5 .5 .5])
        end
    %  plot(time_base_pulse,data_around_pulse_angle,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_temp,2),'-k');
end
ylabel('Temp of batch (C)');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);



figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_puff')', nansem(data_around_pulse_puff')', time_base_pulse', [.5 .5 .5])
        end
    %  plot(time_base_pulse,data_around_pulse_angle,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_puff,2),'-k');
end
ylabel('Air puff');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);


end