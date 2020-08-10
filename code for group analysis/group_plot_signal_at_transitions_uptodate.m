function [transition_structure, transition_structure_notperfly] = group_plot_signal_at_transitions_uptodate(recordings_list, ROI_num, filter_out_chrimson_data, filter_out_egg_times, backsub)
% ROI_num = -1 use patch
%% plot signal after transitions

% this is all transitions and never plots a second transitions
% data is stopped when there is another transition in the positive or
% negative portions of this plot

% transition structure has a few thigns added to it (data)
% col 1,2,3,4,5,6 are:
% 1 time of transition (time onto the new substrate)
% 2 current sbstrate
% 3 previous substrate
% 4 time of last 0
% 5 time of last 200
% 6 time of last 500
% 7 is the speed in the 10 secounds around transition

transition_structure_notperfly.s200_to_p    = [];
transition_structure_notperfly.p_to_s200    = [];
transition_structure_notperfly.s500_to_p    = [];
transition_structure_notperfly.p_to_s500    = [];
transition_structure_notperfly.s500_to_s200 = [];
transition_structure_notperfly.s200_to_s500 = [];

transition_structure_notperfly.s200_to_p_vel    = [];
transition_structure_notperfly.p_to_s200_vel    = [];
transition_structure_notperfly.s500_to_p_vel    = [];
transition_structure_notperfly.p_to_s500_vel    = [];
transition_structure_notperfly.s500_to_s200_vel = [];
transition_structure_notperfly.s200_to_s500_vel = [];

transition_structure_notperfly.s200_to_p_vel_noave    = [];
transition_structure_notperfly.p_to_s200_vel_noave     = [];
transition_structure_notperfly.s500_to_p_vel_noave     = [];
transition_structure_notperfly.p_to_s500_vel_noave     = [];
transition_structure_notperfly.s500_to_s200_vel_noave  = [];
transition_structure_notperfly.s200_to_s500_vel_noave  = [];

transition_structure_notperfly.s200_to_p_speed_2hz_hold    = [];
transition_structure_notperfly.p_to_s200_speed_2hz_hold    = [];
transition_structure_notperfly.s500_to_p_speed_2hz_hold    = [];
transition_structure_notperfly.p_to_s500_speed_2hz_hold    = [];
transition_structure_notperfly.s500_to_s200_speed_2hz_hold = [];
transition_structure_notperfly.s200_to_s500_speed_2hz_hold = [];

transition_structure_notperfly.s200_to_p_mean    = [];
transition_structure_notperfly.p_to_s200_mean    = [];
transition_structure_notperfly.s500_to_p_mean    = [];
transition_structure_notperfly.p_to_s500_mean    = [];
transition_structure_notperfly.s500_to_s200_mean = [];
transition_structure_notperfly.s200_to_s500_mean = [];

transition_structure_notperfly.s200_to_p_fo_is_mean    = [];
transition_structure_notperfly.p_to_s200_fo_is_mean     = [];
transition_structure_notperfly.s500_to_p_fo_is_mean     = [];
transition_structure_notperfly.p_to_s500_fo_is_mean     = [];
transition_structure_notperfly.s500_to_s200_fo_is_mean  = [];
transition_structure_notperfly.s200_to_s500_fo_is_mean  = [];

transition_structure_notperfly.s200_to_p_filtwheel    = [];
transition_structure_notperfly.p_to_s200_filtwheel    = [];
transition_structure_notperfly.s500_to_p_filtwheel    = [];
transition_structure_notperfly.p_to_s500_filtwheel    = [];
transition_structure_notperfly.s500_to_s200_filtwheel = [];
transition_structure_notperfly.s200_to_s500_filtwheel = [];

transition_structure_notperfly.s200_to_p_data    = [];
transition_structure_notperfly.p_to_s200_data    = [];
transition_structure_notperfly.s500_to_p_data    = [];
transition_structure_notperfly.p_to_s500_data    = [];
transition_structure_notperfly.s500_to_s200_data = [];
transition_structure_notperfly.s200_to_s500_data = [];

transition_structure_notperfly.s200_to_p_median_proboscis    = [];
transition_structure_notperfly.p_to_s200_median_proboscis     = [];
transition_structure_notperfly.s500_to_p_median_proboscis      = [];
transition_structure_notperfly.p_to_s500_median_proboscis      = [];
transition_structure_notperfly.s500_to_s200_median_proboscis   = [];
transition_structure_notperfly.s200_to_s500_median_proboscis   = [];


for rec_index = 1:1:length(recordings_list)
    modifiedStr = strrep([char(recordings_list(rec_index))], '.mat', '_stripped.mat');
    modifiedStr2 = strrep([char(recordings_list(rec_index))], '.mat', '_stripped_new_dlc_track.mat');
    
    if(exist([modifiedStr2]))
        modifiedStr = modifiedStr2;
    end
    %   modifiedStr = [char(recordings_list(rec_index))];
    
    recording = loadsinglerecording(modifiedStr);
    
    if(ROI_num(rec_index) ==-1)
        patch = 1;
    else
        patch = 0;
    end
    if(patch ==1)
        recording.tseries = 1;
    end
    
    [recording] = processes_more_DLC_variables(recording);
    if(backsub == 1)
        [recording] = replace_df_over_f_withbackgroundsubtracted(recording);
    end
        
    cntr = 0;
    cell_store_trans = {};
    cell_store_trans_mean = {};
    cell_store_trans_fo_is_mean = {};

    cell_store_trans_timebase = {};
    cell_store_filtered_wheel = {};
    cell_store_median_proboscis = {};
    cell_store_sucrose = {};
    cell_store_vel = [];
    cell_store_vel_noave = [];
    cell_store_speed_2hz_hold = [];

    store_trans = [];
    store_trans_moredetails = [];
    
    out_signal = [];
    out_signal_mean = [];
    out_signal_fo_is_mean = [];

    out_sucrose = [];
    out_filtered_wheel = [];
    out_timebase = [];
    out_median_proboscis = [];
    out_vel = [];
    out_vel_noave = [];
    out_speed_2hz_hold = [];

    %% calculate speed/velocity
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
    
    
    %% start main loop
    for ts = 1:1:length(recording.tseries)
        ds = diff(recording.movie1.sucrose);
        ds = [0,ds];
        [a b] = find(ds);
        b = [1,b, length(recording.movie1.sucrose)];
        
        for i =2:1:(length(b)-1)
            
            startt =recording.movie1.time_stamps(b(i-1));
            stopt = recording.movie1.time_stamps(b(i+1)-1);
            
            if(patch == 0)
                [time_base_to_return, data_to_average_interp] = average_around_event(recording.tseries(ts).df_over_f(ROI_num(rec_index),:), recording.tseries(ts).Time_s, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
                [time_base_to_return, data_to_average_interp_mean] = average_around_event(recording.tseries(ts).df_over_f(ROI_num(rec_index),:)./nanmean(recording.tseries(ts).df_over_f(ROI_num(rec_index),:)), recording.tseries(ts).Time_s, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));        
                [time_base_to_return, data_to_average_interp_fo_is_mean] = average_around_event( (recording.tseries(ts).mean_in_ROI(ROI_num(rec_index),:)./nanmean(recording.tseries(ts).mean_in_ROI(ROI_num(rec_index),:))-1), recording.tseries(ts).Time_s, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));

            end
            
            if(patch == 1)
                [time_base_to_return, data_to_average_interp] = average_around_event(recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
                data_to_average_interp_fo_is_mean = data_to_average_interp;
                data_to_average_interp_mean = data_to_average_interp;
                data_to_average_interp = data_to_average_interp;
            end
            
            if(patch == 0 && filter_out_chrimson_data == 1)
                [time_base_to_return_laser, data_to_average_interp_laser] = average_around_event(recording.abf.no_laser2, recording.abf.Time_s, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
                data_to_average_interp_laser(data_to_average_interp_laser>0) = NaN;
                data_to_average_interp_laser(data_to_average_interp_laser==0) = 1;
                data_to_average_interp = data_to_average_interp.*data_to_average_interp_laser;
                data_to_average_interp_mean = data_to_average_interp_mean.*data_to_average_interp_laser;
                data_to_average_interp_fo_is_mean = data_to_average_interp_fo_is_mean.*data_to_average_interp_laser;               
            end
            
            if(patch == 0 && filter_out_egg_times == 1)
                [time_base_to_return_laser, data_to_average_interp_noegg] = average_around_event(recording.movie1.no_egg2, recording.movie1.time_stamps, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
                data_to_average_interp_noegg(data_to_average_interp_noegg>0) = NaN;
                data_to_average_interp_noegg(data_to_average_interp_noegg==0) = 1;
                data_to_average_interp = data_to_average_interp.*data_to_average_interp_noegg;
                data_to_average_interp_mean = data_to_average_interp_mean.*data_to_average_interp_noegg;
                data_to_average_interp_fo_is_mean = data_to_average_interp_fo_is_mean.*data_to_average_interp_noegg;
            end
            
            if(patch == 1  && filter_out_chrimson_data == 1)
                [time_base_to_return_laser, data_to_average_interp_laser] = average_around_event(recording.abf.no_laser2(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
                data_to_average_interp_laser(data_to_average_interp_laser>0) = NaN;
                data_to_average_interp_laser(data_to_average_interp_laser==0) = 1;
                data_to_average_interp = data_to_average_interp.*data_to_average_interp_laser;
                data_to_average_interp_mean = data_to_average_interp_mean.*data_to_average_interp_laser;
                data_to_average_interp_fo_is_mean = data_to_average_interp_fo_is_mean.*data_to_average_interp_laser;

            end
            
            if(patch == 1 && filter_out_egg_times == 1)
                disp 'Code ot yet written'
            end
            
            [time_base_to_return1, data_to_average_interp1] = average_around_event(recording.movie1.sucrose, recording.movie1.time_stamps, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
            [time_base_to_return2, data_to_average_interp2] = average_around_event(unwrap(recording.movie1.filtered_wheel), recording.movie1.time_stamps, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
            [time_base_to_return3, data_to_average_interp3] = average_around_event(recording.movie1.prob_length./nanmedian(recording.movie1.prob_length), recording.movie1.time_stamps, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
            [time_base_to_return4, data_to_average_interp4] = average_around_event(wheelvelocity, recording.movie1.time_stamps, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
            [time_base_to_return5, data_to_average_interp5] = average_around_event(wheelvelocity_noave, recording.movie1.time_stamps, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
            [time_base_to_return6, data_to_average_interp6] = average_around_event(speed_2hz_hold, recording.movie1.time_stamps, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));

            % time between transition (i-1) and i+1 has to be greater
            % than 4 seconds and has to be in transition zone (removes
            % errant)
            
            if( (b(i+1)-b(i-1)) > 100)  
                cntr = cntr+1; 
                cell_store_trans_timebase{cntr} = time_base_to_return;
                cell_store_trans{cntr} =  data_to_average_interp;
                cell_store_trans_mean{cntr} =  data_to_average_interp_mean;
                cell_store_trans_fo_is_mean{cntr} =  data_to_average_interp_fo_is_mean;
                cell_store_sucrose{cntr} = data_to_average_interp1;
                cell_store_filtered_wheel{cntr} = data_to_average_interp2;
                cell_store_median_proboscis{cntr} = data_to_average_interp3;
                cell_store_vel{cntr} = data_to_average_interp4;
                cell_store_vel_noave{cntr} = data_to_average_interp5;
                cell_store_speed_2hz_hold{cntr} = data_to_average_interp6;
    
                store_trans(cntr,:) = [b(i-1), b(i), b(i+1)-1];
            end
        end
    end
    
    for i = 1:1:cntr
        out_signal(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_trans{i}, -240:.1:900,'previous');
        out_signal_mean(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_trans_mean{i}, -240:.1:900,'previous');
        out_signal_fo_is_mean(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_trans_fo_is_mean{i}, -240:.1:900,'previous');

        out_sucrose(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_sucrose{i}, -240:.1:900,'previous');
        out_filtered_wheel(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_filtered_wheel{i}, -240:.1:900,'previous');
        out_timebase(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_trans_timebase{i}, -240:.1:900,'previous');
        out_median_proboscis(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_median_proboscis{i}, -240:.1:900,'previous');
        out_vel(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_vel{i}, -240:.1:900,'previous');
        out_vel_noave(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_vel_noave{i}, -240:.1:900,'previous');
        out_speed_2hz_hold(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_speed_2hz_hold{i}, -240:.1:900,'previous');
    end
    
    out_signal_copy = out_signal;
    out_signal(~any(~isnan(out_signal), 2),:)=[];
    out_signal_mean(~any(~isnan(out_signal_copy), 2),:)=[];
    out_signal_fo_is_mean(~any(~isnan(out_signal_copy), 2),:)=[];

    out_sucrose(~any(~isnan(out_signal_copy), 2),:)=[];
    out_filtered_wheel(~any(~isnan(out_signal_copy), 2),:)=[];
    out_timebase(~any(~isnan(out_signal_copy), 2),:)=[];
    out_median_proboscis(~any(~isnan(out_signal_copy), 2),:)=[];
    out_vel(~any(~isnan(out_signal_copy), 2),:)=[];
    out_vel_noave(~any(~isnan(out_signal_copy), 2),:)=[];
    out_speed_2hz_hold(~any(~isnan(out_signal_copy), 2),:)=[];
    store_trans(~any(~isnan(out_signal_copy), 2),:)=[];
    
    % plotting after calculations are made
    % no filtering of data here (except transition time)
    % first just plot where each of the transitions are
    
    %%%%%%
    
    store_trans_moredetails = store_trans(:,2);
    
    [rr,c] = size(out_signal);
    s200_to_p = [];
    p_to_s200 = []; 
    s200_to_s500 = [];
    s500_to_s200 = [];
    s500_to_p = [];
    p_to_s500 = [];

    for r = 1:1:rr
        store_trans_moredetails(r,2) = recording.movie1.sucrose(store_trans(r,2));
        store_trans_moredetails(r,3) = recording.movie1.sucrose(store_trans(r,2)-1);
        
        [last0 last0p] = find(recording.movie1.sucrose(1:(store_trans(r,2)-1)) == 0,1,'last');
        [last500 last500p] = find(recording.movie1.sucrose(1:(store_trans(r,2)-1)) == 500,1,'last');
        [last200 last200p] = find(recording.movie1.sucrose(1:(store_trans(r,2)-1)) == 200,1,'last');
        
        if(isempty(last0p))
            last0p = inf;
        end
        if(isempty(last500p))
            last500p = inf;
        end
        if(isempty(last200p))
            last200p = inf;
        end
        store_trans_moredetails(r,4) =  store_trans(r,2)-last0p;
        store_trans_moredetails(r,5) =  store_trans(r,2)-last200p;
        store_trans_moredetails(r,6) =  store_trans(r,2)-last500p;
        
        if(recording.movie1.sucrose(store_trans(r,2)) == 200 && recording.movie1.sucrose(store_trans(r,2)-1) == 0)
            p_to_s200 = [p_to_s200, r];
        end
        
        if(recording.movie1.sucrose(store_trans(r,2)) == 0 && recording.movie1.sucrose(store_trans(r,2)-1) == 200)
            s200_to_p = [s200_to_p, r];    
        end
        
        if(recording.movie1.sucrose(store_trans(r,2)) == 500 && recording.movie1.sucrose(store_trans(r,2)-1) == 0)
            p_to_s500 = [p_to_s500, r];
        end
        
        if(recording.movie1.sucrose(store_trans(r,2)) == 0 && recording.movie1.sucrose(store_trans(r,2)-1) == 500)
            s500_to_p = [s500_to_p, r];    
        end
        
        if(recording.movie1.sucrose(store_trans(r,2)) == 500 && recording.movie1.sucrose(store_trans(r,2)-1) == 200)
            s200_to_s500 = [s200_to_s500, r];
        end
        
        if(recording.movie1.sucrose(store_trans(r,2)) == 200 && recording.movie1.sucrose(store_trans(r,2)-1) == 500)
            s500_to_s200 = [s500_to_s200, r];    
        end  
    end
    
    transition_structure(rec_index).s200_to_p = nanmean(out_signal(s200_to_p,:),1)';
    transition_structure(rec_index).p_to_s200 = nanmean(out_signal(p_to_s200,:),1)';
    transition_structure(rec_index).s500_to_p = nanmean(out_signal(s500_to_p,:),1)';
    transition_structure(rec_index).p_to_s500 = nanmean(out_signal(p_to_s500,:),1)';
    transition_structure(rec_index).s500_to_s200 = nanmean(out_signal(s500_to_s200,:),1)';
    transition_structure(rec_index).s200_to_s500 = nanmean(out_signal(s200_to_s500,:),1)';
    
    transition_structure(rec_index).s200_to_p_mean = nanmean(out_signal_mean(s200_to_p,:),1)';
    transition_structure(rec_index).p_to_s200_mean = nanmean(out_signal_mean(p_to_s200,:),1)';
    transition_structure(rec_index).s500_to_p_mean = nanmean(out_signal_mean(s500_to_p,:),1)';
    transition_structure(rec_index).p_to_s500_mean = nanmean(out_signal_mean(p_to_s500,:),1)';
    transition_structure(rec_index).s500_to_s200_mean = nanmean(out_signal_mean(s500_to_s200,:),1)';
    transition_structure(rec_index).s200_to_s500_mean = nanmean(out_signal_mean(s200_to_s500,:),1)';
    
        transition_structure(rec_index).s200_to_p_fo_is_mean = nanmean(out_signal_fo_is_mean(s200_to_p,:),1)';
    transition_structure(rec_index).p_to_s200_fo_is_mean = nanmean(out_signal_fo_is_mean(p_to_s200,:),1)';
    transition_structure(rec_index).s500_to_p_fo_is_mean = nanmean(out_signal_fo_is_mean(s500_to_p,:),1)';
    transition_structure(rec_index).p_to_s500_fo_is_mean = nanmean(out_signal_fo_is_mean(p_to_s500,:),1)';
    transition_structure(rec_index).s500_to_s200_fo_is_mean = nanmean(out_signal_fo_is_mean(s500_to_s200,:),1)';
    transition_structure(rec_index).s200_to_s500_fo_is_mean = nanmean(out_signal_fo_is_mean(s200_to_s500,:),1)';
    
    transition_structure_notperfly.s200_to_p    = [transition_structure_notperfly.s200_to_p ; (out_signal(s200_to_p,:))];
    transition_structure_notperfly.p_to_s200    = [transition_structure_notperfly.p_to_s200 ;(out_signal(p_to_s200,:))];
    transition_structure_notperfly.s500_to_p    = [transition_structure_notperfly.s500_to_p ;(out_signal(s500_to_p,:))];
    transition_structure_notperfly.p_to_s500    = [transition_structure_notperfly.p_to_s500 ;(out_signal(p_to_s500,:))];
    transition_structure_notperfly.s500_to_s200 = [transition_structure_notperfly.s500_to_s200 ;(out_signal(s500_to_s200,:))];
    transition_structure_notperfly.s200_to_s500 = [transition_structure_notperfly.s200_to_s500 ;(out_signal(s200_to_s500,:))];
    
    transition_structure_notperfly.s200_to_p_mean    = [transition_structure_notperfly.s200_to_p_mean ; (out_signal_mean(s200_to_p,:))];
    transition_structure_notperfly.p_to_s200_mean    = [transition_structure_notperfly.p_to_s200_mean ;(out_signal_mean(p_to_s200,:))];
    transition_structure_notperfly.s500_to_p_mean    = [transition_structure_notperfly.s500_to_p_mean ;(out_signal_mean(s500_to_p,:))];
    transition_structure_notperfly.p_to_s500_mean    = [transition_structure_notperfly.p_to_s500_mean ;(out_signal_mean(p_to_s500,:))];
    transition_structure_notperfly.s500_to_s200_mean = [transition_structure_notperfly.s500_to_s200_mean ;(out_signal_mean(s500_to_s200,:))];
    transition_structure_notperfly.s200_to_s500_mean = [transition_structure_notperfly.s200_to_s500_mean ;(out_signal_mean(s200_to_s500,:))];
    
        transition_structure_notperfly.s200_to_p_fo_is_mean    = [transition_structure_notperfly.s200_to_p_fo_is_mean ; (out_signal_fo_is_mean(s200_to_p,:))];
    transition_structure_notperfly.p_to_s200_fo_is_mean    = [transition_structure_notperfly.p_to_s200_fo_is_mean ;(out_signal_fo_is_mean(p_to_s200,:))];
    transition_structure_notperfly.s500_to_p_fo_is_mean    = [transition_structure_notperfly.s500_to_p_fo_is_mean ;(out_signal_fo_is_mean(s500_to_p,:))];
    transition_structure_notperfly.p_to_s500_fo_is_mean    = [transition_structure_notperfly.p_to_s500_fo_is_mean ;(out_signal_fo_is_mean(p_to_s500,:))];
    transition_structure_notperfly.s500_to_s200_fo_is_mean = [transition_structure_notperfly.s500_to_s200_fo_is_mean ;(out_signal_fo_is_mean(s500_to_s200,:))];
    transition_structure_notperfly.s200_to_s500_fo_is_mean = [transition_structure_notperfly.s200_to_s500_fo_is_mean ;(out_signal_fo_is_mean(s200_to_s500,:))];
    
    transition_structure_notperfly.s200_to_p_filtwheel    = [transition_structure_notperfly.s200_to_p_filtwheel ; (out_filtered_wheel(s200_to_p,:))];
    transition_structure_notperfly.p_to_s200_filtwheel   = [transition_structure_notperfly.p_to_s200_filtwheel ;(out_filtered_wheel(p_to_s200,:))];
    transition_structure_notperfly.s500_to_p_filtwheel    = [transition_structure_notperfly.s500_to_p_filtwheel ;(out_filtered_wheel(s500_to_p,:))];
    transition_structure_notperfly.p_to_s500_filtwheel    = [transition_structure_notperfly.p_to_s500_filtwheel ;(out_filtered_wheel(p_to_s500,:))];
    transition_structure_notperfly.s500_to_s200_filtwheel = [transition_structure_notperfly.s500_to_s200_filtwheel ;(out_filtered_wheel(s500_to_s200,:))];
    transition_structure_notperfly.s200_to_s500_filtwheel = [transition_structure_notperfly.s200_to_s500_filtwheel ;(out_filtered_wheel(s200_to_s500,:))];
    
    transition_structure_notperfly.s200_to_p_vel    = [transition_structure_notperfly.s200_to_p_vel ; (out_vel(s200_to_p,:))];
    transition_structure_notperfly.p_to_s200_vel   = [transition_structure_notperfly.p_to_s200_vel ;(out_vel(p_to_s200,:))];
    transition_structure_notperfly.s500_to_p_vel    = [transition_structure_notperfly.s500_to_p_vel ;(out_vel(s500_to_p,:))];
    transition_structure_notperfly.p_to_s500_vel    = [transition_structure_notperfly.p_to_s500_vel ;(out_vel(p_to_s500,:))];
    transition_structure_notperfly.s500_to_s200_vel = [transition_structure_notperfly.s500_to_s200_vel ;(out_vel(s500_to_s200,:))];
    transition_structure_notperfly.s200_to_s500_vel = [transition_structure_notperfly.s200_to_s500_vel ;(out_vel(s200_to_s500,:))];
    
    transition_structure_notperfly.s200_to_p_vel_noave    = [transition_structure_notperfly.s200_to_p_vel_noave ; (out_vel_noave(s200_to_p,:))];
    transition_structure_notperfly.p_to_s200_vel_noave   = [transition_structure_notperfly.p_to_s200_vel_noave ;(out_vel_noave(p_to_s200,:))];
    transition_structure_notperfly.s500_to_p_vel_noave    = [transition_structure_notperfly.s500_to_p_vel_noave ;(out_vel_noave(s500_to_p,:))];
    transition_structure_notperfly.p_to_s500_vel_noave    = [transition_structure_notperfly.p_to_s500_vel_noave ;(out_vel_noave(p_to_s500,:))];
    transition_structure_notperfly.s500_to_s200_vel_noave = [transition_structure_notperfly.s500_to_s200_vel_noave ;(out_vel_noave(s500_to_s200,:))];
    transition_structure_notperfly.s200_to_s500_vel_noave = [transition_structure_notperfly.s200_to_s500_vel_noave ;(out_vel_noave(s200_to_s500,:))];
    
    transition_structure_notperfly.s200_to_p_speed_2hz_hold    = [transition_structure_notperfly.s200_to_p_speed_2hz_hold ; (out_speed_2hz_hold(s200_to_p,:))];
    transition_structure_notperfly.p_to_s200_speed_2hz_hold   = [transition_structure_notperfly.p_to_s200_speed_2hz_hold ;(out_speed_2hz_hold(p_to_s200,:))];
    transition_structure_notperfly.s500_to_p_speed_2hz_hold   = [transition_structure_notperfly.s500_to_p_speed_2hz_hold ;(out_speed_2hz_hold(s500_to_p,:))];
    transition_structure_notperfly.p_to_s500_speed_2hz_hold    = [transition_structure_notperfly.p_to_s500_speed_2hz_hold ;(out_speed_2hz_hold(p_to_s500,:))];
    transition_structure_notperfly.s500_to_s200_speed_2hz_hold = [transition_structure_notperfly.s500_to_s200_speed_2hz_hold ;(out_speed_2hz_hold(s500_to_s200,:))];
    transition_structure_notperfly.s200_to_s500_speed_2hz_hold = [transition_structure_notperfly.s200_to_s500_speed_2hz_hold ;(out_speed_2hz_hold(s200_to_s500,:))];
    
    transition_structure_notperfly.s200_to_p_median_proboscis    = [transition_structure_notperfly.s200_to_p_median_proboscis   ; (out_median_proboscis(s200_to_p,:))];
    transition_structure_notperfly.p_to_s200_median_proboscis     = [transition_structure_notperfly.p_to_s200_median_proboscis   ;(out_median_proboscis(p_to_s200,:))];
    transition_structure_notperfly.s500_to_p_median_proboscis      = [transition_structure_notperfly.s500_to_p_median_proboscis   ;(out_median_proboscis(s500_to_p,:))];
    transition_structure_notperfly.p_to_s500_median_proboscis      = [transition_structure_notperfly.p_to_s500_median_proboscis   ;(out_median_proboscis(p_to_s500,:))];
    transition_structure_notperfly.s500_to_s200_median_proboscis   = [transition_structure_notperfly.s500_to_s200_median_proboscis   ;(out_median_proboscis(s500_to_s200,:))];
    transition_structure_notperfly.s200_to_s500_median_proboscis   = [transition_structure_notperfly.s200_to_s500_median_proboscis   ;(out_median_proboscis(s200_to_s500,:))];
    
    transition_structure_notperfly.s200_to_p_data    = [transition_structure_notperfly.s200_to_p_data ; (store_trans_moredetails(s200_to_p,:))];
    transition_structure_notperfly.p_to_s200_data    = [transition_structure_notperfly.p_to_s200_data ;(store_trans_moredetails(p_to_s200,:))];
    transition_structure_notperfly.s500_to_p_data    = [transition_structure_notperfly.s500_to_p_data ;(store_trans_moredetails(s500_to_p,:))];
    transition_structure_notperfly.p_to_s500_data    = [transition_structure_notperfly.p_to_s500_data ;(store_trans_moredetails(p_to_s500,:))];
    transition_structure_notperfly.s500_to_s200_data = [transition_structure_notperfly.s500_to_s200_data ;(store_trans_moredetails(s500_to_s200,:))];
    transition_structure_notperfly.s200_to_s500_data = [transition_structure_notperfly.s200_to_s500_data ;(store_trans_moredetails(s200_to_s500,:))];
    
end
% 
% %% plottig code for individual flys and mean
% 
% f_200_to_0 = figure;
% f_0_to_200 = figure;
% f_200_to_500 = figure;
% f_500_to_200 = figure;
% f_500_to_0 = figure;
% f_0_to_500 = figure;
% 
% 
% figure(f_200_to_0); hold on; title(['200 mM sucrose to plain transitions']);
% %ylim([-0.5 3])
% set(gca,'TickDir','out');
% 
% figure(f_0_to_200); hold on; title(['plain to 200 mM sucrose transitions' ]);
% %ylim([-0.5 3])
% set(gca,'TickDir','out');
% 
% figure(f_500_to_0); hold on; title(['500 mM sucrose to plain transitions' ]);
% %ylim([-0.5 3])
% set(gca,'TickDir','out');
% 
% figure(f_0_to_500); hold on; title(['plain to 500 mM sucrose transitions' ]);
% %ylim([-0.5 3])
% set(gca,'TickDir','out');
% 
% figure(f_200_to_500); hold on; title(['200 mM sucrose to 500 mM sucrose transitions' ]);
% %ylim([-0.5 3])
% set(gca,'TickDir','out');
% 
% figure(f_500_to_200); hold on; title(['500 mM sucrose to 200 mM sucrose transitions' ]);
% %ylim([-0.5 3])
% set(gca,'TickDir','out');
% 
% out_signal = [transition_structure(1:end).s200_to_p];
% 
% figure(f_200_to_0); hold on;
% if(~isempty(out_signal))
%     qq = nanmean(out_signal,2);
%     plot(-240:.1:900,(out_signal),'color',[220 220 220]./255,'linewidth',1);
%     plot(-240:.1:0,qq(1:2401),'color',[33 113 181]./255,'linewidth',2);
%     plot(0:.1:900,qq((end-9000):end),'color',[107 174 214]./255,'linewidth',2);
%     
%     xlim([-45,45]);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out_signal = [transition_structure(1:end).p_to_s200];
% 
% if(~isempty(out_signal))
%     
%     figure(f_0_to_200); hold on;
%     qq = nanmean(out_signal,2);
%     plot(-240:.1:900,(out_signal),'color',[220 220 220]./255,'linewidth',1);
%     plot(-240:.1:0,qq(1:2401),'color',[107 174 214]./255,'linewidth',2);
%     plot(0:.1:900,qq((end-9000):end),'color',[33 113 181]./255,'linewidth',2);
%     
%     xlim([-45,45]);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out_signal = [transition_structure(1:end).s500_to_p];
% if(~isempty(out_signal))
%     
%     figure(f_500_to_0); hold on;
%     qq = nanmean(out_signal,2);
%     plot(-240:.1:900,(out_signal),'color',[220 220 220]./255,'linewidth',1);
%     plot(-240:.1:0,qq(1:2401),'color',[8 48 107]./255,'linewidth',2);
%     plot(0:.1:900,qq((end-9000):end),'color',[107 174 214]./255,'linewidth',2);
%     
%     xlim([-45,45]);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out_signal = [transition_structure(1:end).p_to_s500];
% if(~isempty(out_signal))
%     
%     figure(f_0_to_500); hold on;
%     qq = nanmean(out_signal,2);
%     plot(-240:.1:900,(out_signal),'color',[220 220 220]./255,'linewidth',1);
%     plot(-240:.1:0,qq(1:2401),'color',[107 174 214]./255,'linewidth',2);
%     plot(0:.1:900,qq((end-9000):end),'color',[8 48 107]./255,'linewidth',2);
%     
%     xlim([-45,45]);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out_signal = [transition_structure(1:end).s500_to_s200];
% if(~isempty(out_signal))
%     
%     figure(f_500_to_200); hold on;
%     qq = nanmean(out_signal,2);
%     plot(-240:.1:900,(out_signal),'color',[220 220 220]./255,'linewidth',1);
%     plot(-240:.1:0,qq(1:2401),'color',[8 48 107]./255,'linewidth',2);
%     plot(0:.1:900,qq((end-9000):end),'color',[33 113 181]./255,'linewidth',2);
%     
%     xlim([-45,45]);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out_signal = [transition_structure(1:end).s200_to_s500];
% if(~isempty(out_signal))
%     
%     figure(f_200_to_500); hold on;
%     qq = nanmean(out_signal,2);
%     plot(-240:.1:900,(out_signal),'color',[220 220 220]./255,'linewidth',1);
%     plot(-240:.1:0,qq(1:2401),'color',[33 113 181]./255,'linewidth',2);
%     plot(0:.1:900,qq((end-9000):end),'color',[8 48 107]./255,'linewidth',2);
%     
%     xlim([-45,45]);
% end
% 
% figure(f_200_to_0);  set(gca,'ylim',[-.5,1]); pause(1); snapnow;
% 
% figure(f_0_to_200); set(gca,'ylim',[-.5,1]);  pause(1);snapnow;
% 
% figure(f_500_to_0); set(gca,'ylim',[-.5,1]);  pause(1); snapnow;
% 
% figure(f_0_to_500); set(gca,'ylim',[-.5,1]); pause(1); snapnow;
% 
% figure(f_200_to_500); set(gca,'ylim',[-.5,1]); pause(1); snapnow;
% 
% figure(f_500_to_200); set(gca,'ylim',[-.5,1]); pause(1); snapnow;
% 
% figure(f_200_to_0);  set(gca,'ylim',[0,.5]); pause(1); snapnow;
% 
% figure(f_0_to_200); set(gca,'ylim',[0,.5]); pause(1);snapnow;
% 
% figure(f_500_to_0);  set(gca,'ylim',[0,.5]);  pause(1); snapnow;
% 
% figure(f_0_to_500);  set(gca,'ylim',[0,.5]); pause(1); snapnow;
% 
% figure(f_200_to_500);  set(gca,'ylim',[0,.5]); pause(1); snapnow;
% 
% figure(f_500_to_200);  set(gca,'ylim',[0,.5]); pause(1); snapnow;
% 
% 
% figure(f_200_to_0); xlim([-120,120]);   pause(1); snapnow;
% figure(f_0_to_200); xlim([-120,120]);  pause(1); snapnow;
% figure(f_500_to_0); xlim([-120,120]);  pause(1);  snapnow;
% figure(f_0_to_500); xlim([-120,120]);  pause(1); snapnow;
% figure(f_200_to_500); xlim([-120,120]); pause(1);  snapnow;
% figure(f_500_to_200); xlim([-120,120]);   pause(1); snapnow;

out=1;

end

