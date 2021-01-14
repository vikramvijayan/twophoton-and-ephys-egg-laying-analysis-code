function [filter_structure] = group_pass_filter(recordings_list, ROI_num, fly_ID, cell_ID, filter_out_chrimson_data, filter_out_egg_times, backsub)
% this function just returns all the data back with .1 interpolation
% can be used to filter signal to find egg-laying like events etc

% set filter chrimson data to 1 if yu want to remove all data during pulses
% (as well as 100 sec after the pulse) check code for exact number

% set ROI num to -1 if you want tu use CH1 patch
% if patch is 1 then use CH1
% if patch is 0 then use CH1 but remove pulses (and 100 seconds after)

% set filter_out_egg_times to 1 if you don't want to have any of the egg
% laying times inclued (currently 300 sec on either side of egg)
%%%%%%%%%%%
%% plotting correlations or time locked average

%for rec_index = 1:1:2
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
    
    %%% this loop below is common code with plot cross correlations %%
    %% with one exception that is noted (a few lines added here, not in plot correlations)
    
    for m = 1:1:length(recording.tseries)
        if(patch == 0)
            signal    = recording.tseries(m).df_over_f(ROI_num(rec_index),:);
            signal_over_mean    = recording.tseries(m).df_over_f(ROI_num(rec_index),:)./nanmean(recording.tseries(m).df_over_f(ROI_num(rec_index),:));
            signal_fo_is_mean    = (recording.tseries(m).mean_in_ROI(ROI_num(rec_index),:)./nanmean(recording.tseries(m).mean_in_ROI(ROI_num(rec_index),:))-1);

            signaltime = recording.tseries(m).Time_s;
            
            if(filter_out_chrimson_data)
                tseries_no_laser = interp1(recording.abf.Time_s, recording.abf.no_laser, recording.tseries(m).Time_s,'previous');
                signal = signal.*tseries_no_laser';
                signal_over_mean = signal_over_mean.*tseries_no_laser';
                signal_fo_is_mean = signal_fo_is_mean.*tseries_no_laser';

            end
            if(filter_out_egg_times == 1)
                tseries_no_egg = interp1(recording.movie1.time_stamps, recording.movie1.no_egg, recording.tseries(m).Time_s,'previous');
                signal = signal.*tseries_no_egg';
                signal_over_mean = signal_over_mean.*tseries_no_egg';
                signal_fo_is_mean = signal_fo_is_mean.*tseries_no_egg';
            end
%             if(filter_out_egg_times == 2)
%                 not yet written
%                 tseries_no_egg = interp1(recording.movie1.time_stamps, recording.movie1.no_egg, recording.tseries(m).Time_s,'previous');
%                 signal = signal.*tseries_no_egg';
%                 signal_over_mean = signal_over_mean.*tseries_no_egg';
%             end
        end
        
        % for patching data, the signal changes
        if(patch == 1)
            signal    = recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
            signaltime = recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
            signal_over_mean = [];
            signal_fo_is_mean = [];
            if(filter_out_chrimson_data)
                signal    = recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)).*recording.abf.no_laser(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
            end
            if(filter_out_egg_times)
                disp 'Error this code has not been written'
            end
            
        end
        
        
        signal_interp    = interp1(signaltime, signal, signaltime(1):.1:signaltime(end));
        signal_over_mean_interp    = interp1(signaltime, signal_over_mean, signaltime(1):.1:signaltime(end));
        signal_fo_is_mean_interp    = interp1(signaltime, signal_fo_is_mean, signaltime(1):.1:signaltime(end));

        signaltime_interp = signaltime(1):.1:signaltime(end);
        
        abd_only_length_x = recording.movie1.abd_x_L3tip./nanmedian(recording.movie1.abd_x_L3tip);
        bodylength_x =  recording.movie1.abd_x_neck_tip./nanmedian(recording.movie1.abd_x_neck_tip);
        abd_only_length_y = recording.movie1.abd_y_L3tip./nanmedian(recording.movie1.abd_y_L3tip);
        bodylength_y =  recording.movie1.abd_y_neck_tip./nanmedian(recording.movie1.abd_y_neck_tip);
        prob_x =  recording.movie1.prob_x./nanmedian(recording.movie1.prob_x);
        prob_y =  recording.movie1.prob_y./nanmedian(recording.movie1.prob_y);
        movtime  = recording.movie1.time_stamps';
        sucrose = recording.movie1.sucrose;
        bodylength =  recording.movie1.abd_length./nanmedian(recording.movie1.abd_length);
        bodyangle =  recording.movie1.abd_angle;
        problength = recording.movie1.prob_length./nanmedian(recording.movie1.prob_length);
        abd_path_length = recording.movie1.abd_path_length./nanmedian(recording.movie1.abd_path_length);
        abd_only_length = recording.movie1.abd_only_length./nanmedian(recording.movie1.abd_only_length);
        abd_only_angle = recording.movie1.abd_only_angle;
        wheelvelocity = (pi*7*-25).*smooth([0; diff(unwrap(recording.movie1.filtered_wheel))],25) ./ (2*pi); %mm/sec 7mm radium
        wheelvelocity = wheelvelocity';
        wheelvelocity_noave = (pi*7*-25).*[0; diff(unwrap(recording.movie1.filtered_wheel))]./ (2*pi); %mm/sec 7mm radium
        wheelvelocity_noave = wheelvelocity_noave';
        eggs_vector = zeros(1,length(recording.movie1.sucrose));
        eggs_vector(recording.movie1.eggs) = 1;
        filtered_wheel = recording.movie1.filtered_wheel';
        
        
        
        
        
        
        
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
        
        
        
        
        
        
        filter_structure(rec_index).tseries(m).signal_interp = signal_interp;
        filter_structure(rec_index).tseries(m).signal_over_mean_interp = signal_over_mean_interp;
        filter_structure(rec_index).tseries(m).signal_fo_is_mean_interp = signal_fo_is_mean_interp;

        filter_structure(rec_index).tseries(m).signaltime_interp = signaltime_interp;
        filter_structure(rec_index).tseries(m).wheelvelocity_interp = interp1(movtime, wheelvelocity, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).wheelvelocity_noave_interp = interp1(movtime, wheelvelocity_noave, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).bodylength_interp =  interp1(movtime, bodylength, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).bodyangle_interp =  interp1(movtime, bodyangle, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).abd_path_length_interp = interp1(movtime, abd_path_length, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).abd_only_length_interp = interp1(movtime, abd_only_length, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).abd_only_angle_interp = interp1(movtime, abd_only_angle, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).abd_only_length_x_interp = interp1(movtime, abd_only_length_x, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).bodylength_x_interp =  interp1(movtime, bodylength_x, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).abd_only_length_y_interp = interp1(movtime, abd_only_length_y, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).bodylength_y_interp =  interp1(movtime, bodylength_y, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).prob_x_interp =  interp1(movtime, prob_x, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).prob_y_interp =  interp1(movtime, prob_y, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).problength_interp = interp1(movtime, problength, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).eggs_vector = interp1(movtime, eggs_vector, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).filtered_wheel = interp1(movtime, filtered_wheel, signaltime(1):.1:signaltime(end));
        filter_structure(rec_index).tseries(m).speed_2hz_hold = interp1(movtime, speed_2hz_hold, signaltime(1):.1:signaltime(end));
        
        filter_structure(rec_index).fly_ID = fly_ID(rec_index);
        filter_structure(rec_index).cell_ID = cell_ID(rec_index);
        
        filter_structure(rec_index).recording = recordings_list(rec_index);
        filter_structure(rec_index).recording_index = rec_index;
        filter_structure(rec_index).ROI_used = ROI_num(rec_index);
    end
    
    
    
    
end