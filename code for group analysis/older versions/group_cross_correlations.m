

function [corr_structure,egg_structure] = group_cross_correlations(recordings_list, ROI_num, filter_out_chrimson_data, filter_chrimson_eggs, filter_out_egg_times, backsub)


% set filter chrimson data to 1 if yu want to remove all data during pulses
% (as well as 100 sec after the pulse) check code for exact number

% set filter chrimson eggs to remove all eggs in the same interval as above
% 0 = no filter
% 1 - filter out chrimson eggs
% 2 - only use chrimson eggs

% set ROI num to 0 if you want tu use CH1 patch
% if patch is 1 then use CH1
% if patch is 0 then use CH1 but remove pulses (and 100 seconds after)

% set filter_out_egg_times to 1 if you don't want to have any of the egg
% laying times inclued (currently 300 sec on either side of egg)

%%%%%%%%%%%
%% plotting correlations or time locked average

for rec_index = 1:1:length(recordings_list)
    %for rec_index = 10
    %  rec_index
    
    % it is often faster to use the stripped files without the associated
    % images
    %modifiedStr = strrep([char(recordings_list(rec_index))], '.mat', '_stripped.mat');
    modifiedStr = [char(recordings_list(rec_index))];
    recording = loadsinglerecording(modifiedStr);
    
    [recording] = processes_more_DLC_variables(recording);
    if(backsub == 1)
        [recording] = replace_df_over_f_withbackgroundsubtracted(recording);
    end
    
    patch = ~ROI_num(rec_index);
if(patch == 1)
    recording.tseries = 1;
end

%%% this loop below is common code with plot cross correlations %%
%% with one exception that is noted (a few lines added here, not in plot correlations)

data_to_plain   = [];
data_to_sucrose = [];

data_to_sucrose_control = [];
data_to_plain_control = [];
time_base_trans = [];

data_around_egg_vel = [];

data_around_egg = [];
data_around_plain_egg = [];
data_around_sucrose_egg = [];
time_base_egg = [];
data_around_egg_sub = [];
data_around_plain_egg_sub = [];
data_around_sucrose_egg_sub = [];

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

data_around_egg2 = [];
data_around_plain_egg2 = [];
data_around_sucrose_egg2 = [];
time_base_egg2 = [];
data_around_egg_sub2 = [];
data_around_plain_egg_sub2 = [];
data_around_sucrose_egg_sub2 = [];

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

corr_to_egg = [];
time_base_corr_to_egg = [];

corr_to_plainegg = [];
time_base_corr_to_plainegg = [];

corr_to_sucroseegg = [];
time_base_corr_to_sucroseegg = [];


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

% process more of the DLC movie

abd_only_length_x = recording.movie1.abd_x_L3tip;
bodylength_x =  recording.movie1.abd_x_neck_tip;
abd_only_length_y = recording.movie1.abd_y_L3tip;
bodylength_y =  recording.movie1.abd_y_neck_tip;

prob_x =  recording.movie1.prob_x;
prob_y =  recording.movie1.prob_y;


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


[time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,all_eggt, [-240:.1:240]);
data_around_egg_sub = [data_around_egg_sub; data_to_average_interp];
time_base_egg = time_base_to_return;

[time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,plain_eggt, [-240:.1:240]);
data_around_plain_egg_sub = [data_around_plain_egg_sub; data_to_average_interp];

[time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,sucrose_eggt, [-240:.1:240]);
data_around_sucrose_egg_sub = [data_around_sucrose_egg_sub; data_to_average_interp];

[time_base_to_return, data_to_average_interp]  = average_around_event(wheelvelocity,movtime,all_eggt, [-240:.1:240]);
data_around_egg_vel = [data_around_egg_vel; data_to_average_interp];


[time_base_to_return, data_to_average_interp]  = average_around_event(problength,movtime,all_eggt, [-240:.1:240]);
data_around_egg_prob = [data_around_egg_prob; data_to_average_interp];

[time_base_to_return, data_to_average_interp]  = average_around_event(prob_x,movtime,all_eggt, [-240:.1:240]);
data_around_egg_prob_x = [data_around_egg_prob_x; data_to_average_interp];

[time_base_to_return, data_to_average_interp]  = average_around_event(prob_y,movtime,all_eggt, [-240:.1:240]);
data_around_egg_prob_y = [data_around_egg_prob_y; data_to_average_interp];


[time_base_to_return, data_to_average_interp]  = average_around_event(bodylength,movtime,all_eggt, [-240:.1:240]);
data_around_egg_body = [data_around_egg_body; data_to_average_interp];

[time_base_to_return, data_to_average_interp]  = average_around_event(bodylength_x,movtime,all_eggt, [-240:.1:240]);
data_around_egg_body_x = [data_around_egg_body_x; data_to_average_interp];

[time_base_to_return, data_to_average_interp]  = average_around_event(bodylength_y,movtime,all_eggt, [-240:.1:240]);
data_around_egg_body_y = [data_around_egg_body_y; data_to_average_interp];


[time_base_to_return, data_to_average_interp]  = average_around_event(abd_path_length,movtime,all_eggt, [-240:.1:240]);
data_around_egg_path = [data_around_egg_path; data_to_average_interp];


[time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_length,movtime,all_eggt, [-240:.1:240]);
data_around_egg_body_only = [data_around_egg_body_only; data_to_average_interp];

[time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_length_x,movtime,all_eggt, [-240:.1:240]);
data_around_egg_body_only_x = [data_around_egg_body_only_x; data_to_average_interp];


[time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_length_y,movtime,all_eggt, [-240:.1:240]);
data_around_egg_body_only_y = [data_around_egg_body_only_y; data_to_average_interp];

[time_base_to_return, data_to_average_interp]  = average_around_event(bodyangle,movtime,all_eggt, [-240:.1:240]);
data_around_egg_angle = [data_around_egg_angle; data_to_average_interp];


[time_base_to_return, data_to_average_interp]  = average_around_event(abd_only_angle,movtime,all_eggt, [-240:.1:240]);
data_around_egg_angle_only = [data_around_egg_angle_only; data_to_average_interp];

[time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,all_eggt, [-120:.1:120]);
data_around_egg_sub2 = [data_around_egg_sub2; data_to_average_interp];
time_base_egg = time_base_to_return;
[time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,plain_eggt, [-120:.1:120]);
data_around_plain_egg_sub2 = [data_around_plain_egg_sub2; data_to_average_interp];
[time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,sucrose_eggt, [-120:.1:120]);
data_around_sucrose_egg_sub2 = [data_around_sucrose_egg_sub2; data_to_average_interp];

data_around_sucrose_egg_sub = repmat(data_around_sucrose_egg_sub,length(recording.tseries),1);
data_around_egg_sub = repmat(data_around_egg_sub,length(recording.tseries),1);
data_around_plain_egg_sub = repmat(data_around_plain_egg_sub,length(recording.tseries),1);

data_around_sucrose_egg_sub2 = repmat(data_around_sucrose_egg_sub2,length(recording.tseries),1);
data_around_egg_sub2 = repmat(data_around_egg_sub2,length(recording.tseries),1);
data_around_plain_egg_sub2 = repmat(data_around_plain_egg_sub2,length(recording.tseries),1);

%% the lines below need to be added to newer versions of plots cross correlations
  all_recordings_list = repmat(recordings_list(rec_index),length(recording.tseries),1);
    all_eggs_listindex = repmat(all_egg,length(recording.tseries),1);
    all_eggs_listtime = repmat(all_eggt,length(recording.tseries),1);
    all_recordings_list_index = repmat(rec_index,length(recording.tseries),1);
    all_ROI_index = repmat(ROI_num(rec_index),length(recording.tseries),1);
    %%
    
for m = 1:1:length(recording.tseries)
    if(patch == 0)
        signal    = recording.tseries(m).df_over_f(ROI_num,:);
        signaltime = recording.tseries(m).Time_s;
        
        if(filter_out_chrimson_data)
            tseries_no_laser = interp1(recording.abf.Time_s, recording.abf.no_laser, recording.tseries(m).Time_s,'previous');
            signal = signal.*tseries_no_laser';
        end
        
        if(filter_out_egg_times)
            tseries_no_egg = interp1(recording.movie1.time_stamps, recording.movie1.no_egg, recording.tseries(m).Time_s,'previous');
            signal = signal.*tseries_no_egg';
        end
    end
    
    % for patching data, the signal changes
    if(patch == 1)
        signal    = recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
        signaltime = recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
        
        if(filter_out_chrimson_data)
            signal    = recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)).*recording.abf.no_laser(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
        end
        
        if(filter_out_egg_times)
            disp 'Error this code has not been written'
        end
        
    end
    
    
    
    wheelvelocity = (pi*7*-25).*smooth([0; diff(unwrap(recording.movie1.filtered_wheel))],25) ./ (2*pi); %mm/sec 7mm radium
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
    
    % look at signal at eggs
    [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,all_eggt, [-240:.1:240]);
    data_around_egg = [data_around_egg; data_to_average_interp];
    time_base_egg = time_base_to_return;
    [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,plain_eggt, [-240:.1:240]);
    data_around_plain_egg = [data_around_plain_egg; data_to_average_interp];
    [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,sucrose_eggt, [-240:.1:240]);
    data_around_sucrose_egg = [data_around_sucrose_egg; data_to_average_interp];
    
    % look at signal at eggs
    [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,all_eggt, [-120:.1:120]);
    data_around_egg2 = [data_around_egg2; data_to_average_interp];
    time_base_egg2 = time_base_to_return;
    [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,plain_eggt, [-120:.1:120]);
    data_around_plain_egg2 = [data_around_plain_egg2; data_to_average_interp];
    [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,sucrose_eggt, [-120:.1:120]);
    data_around_sucrose_egg2 = [data_around_sucrose_egg2; data_to_average_interp];
    
    % look at signal correlated to velocity (smoothed)
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, wheelvelocity, signaltime, movtime, [-120:.1:120]);
    corr_to_vel = [corr_to_vel; out_corr];
    time_base_corr_to_vel = out_corr_shift_ofvec1;
    
    % look at signal correlated to body position
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, problength, signaltime, movtime, [-120:.1:120]);
    corr_to_prob = [corr_to_prob; out_corr];
    time_base_corr_to_prob = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, prob_x, signaltime, movtime, [-120:.1:120]);
    corr_to_prob_x = [corr_to_prob_x; out_corr];
    time_base_corr_to_prob = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, prob_y, signaltime, movtime, [-120:.1:120]);
    corr_to_prob_y = [corr_to_prob_y; out_corr];
    time_base_corr_to_prob = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodylength, signaltime, movtime, [-120:.1:120]);
    corr_to_body = [corr_to_body; out_corr];
    time_base_corr_to_body = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodylength_x, signaltime, movtime, [-120:.1:120]);
    corr_to_body_x = [corr_to_body_x; out_corr];
    time_base_corr_to_body = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodylength_y, signaltime, movtime, [-120:.1:120]);
    corr_to_body_y = [corr_to_body_y; out_corr];
    time_base_corr_to_body = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodyangle, signaltime, movtime, [-120:.1:120]);
    corr_to_angle = [corr_to_angle; out_corr];
    time_base_corr_to_angle = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_path_length, signaltime, movtime, [-120:.1:120]);
    corr_to_body_path = [corr_to_body_path; out_corr];
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_only_length, signaltime, movtime, [-120:.1:120]);
    corr_to_body_only = [corr_to_body_only; out_corr];
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_only_length_x, signaltime, movtime, [-120:.1:120]);
    corr_to_body_only_x = [corr_to_body_only_x; out_corr];
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_only_length_y, signaltime, movtime, [-120:.1:120]);
    corr_to_body_only_y = [corr_to_body_only_y; out_corr];
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, abd_only_angle, signaltime, movtime, [-120:.1:120]);
    corr_to_angle_only = [corr_to_angle_only; out_corr];
    
    % look at signal correlated to substrate
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, sucrose, signaltime, movtime, [-120:.1:120]);
    corr_to_sucrose = [corr_to_sucrose; out_corr];
    time_base_corr_to_sucrose = out_corr_shift_ofvec1;
    
    % look at signal correlated to substrate threshold
    suc_thr = sucrose;
    suc_thr(suc_thr>0) = 1;
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, suc_thr, signaltime, movtime, [-120:.1:120]);
    corr_to_sucrose_t = [corr_to_sucrose_t; out_corr];
    time_base_corr_to_sucrose_t = out_corr_shift_ofvec1;
    
    % look at signal correlated to egg
    tmp_ar = zeros(1,length(movtime));
    tmp_ar(all_egg) = 1;
    tmp_ar = smooth(tmp_ar,25*10);
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, [-120:.1:120]);
    corr_to_egg = [corr_to_egg; out_corr];
    time_base_corr_to_egg = out_corr_shift_ofvec1;
    
    % look at signal correlated to plain egg
    tmp_ar = zeros(1,length(movtime));
    tmp_ar(plain_egg) = 1;
    tmp_ar = smooth(tmp_ar,25*10);
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, [-120:.1:120]);
    corr_to_plainegg = [corr_to_plainegg; out_corr];
    time_base_corr_to_plainegg = out_corr_shift_ofvec1;
    
    % look at signal correlated to sucrose egg
    tmp_ar = zeros(1,length(movtime));
    tmp_ar(sucrose_egg) = 1;
    tmp_ar = smooth(tmp_ar,25*10);
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, [-120:.1:120]);
    corr_to_sucroseegg = [corr_to_sucroseegg; out_corr];
    time_base_corr_to_sucroseegg = out_corr_shift_ofvec1;
    
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
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal,[0;transitionssigned], signaltime, movtime, [-120:.1:120]);
    corr_to_trans = [corr_to_trans; out_corr];
    time_base_corr_to_trans = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal,[0;transitionsabsval], signaltime, movtime, [-120:.1:120]);
    corr_to_absvaltrans = [corr_to_absvaltrans; out_corr];
    time_base_corr_to_absvaltrans = out_corr_shift_ofvec1;
end
   
%%%% the loop above can be replaced with new vweaions of the
%%%% plot_cross_correlation function. It is ahared

    egg_structure(rec_index).data_around_egg_sub = data_around_egg_sub;
    egg_structure(rec_index).data_around_egg_vel = repmat(data_around_egg_vel,length(recording.tseries),1)';
    
    egg_structure(rec_index).data_around_egg_body = repmat(data_around_egg_body,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_body_x = repmat(data_around_egg_body_x,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_body_y = repmat(data_around_egg_body_y,length(recording.tseries),1)';
    
    egg_structure(rec_index).data_around_egg_path = repmat(data_around_egg_path,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_body_only = repmat(data_around_egg_body_only,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_body_only_x = repmat(data_around_egg_body_only_x,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_body_only_y = repmat(data_around_egg_body_only_y,length(recording.tseries),1)';
    
    egg_structure(rec_index).data_around_egg_angle = repmat(data_around_egg_angle,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg_angle_only = repmat(data_around_egg_angle_only,length(recording.tseries),1)';
    egg_structure(rec_index).data_around_egg = data_around_egg';
    egg_structure(rec_index).recording = all_recordings_list;
    egg_structure(rec_index).recording_index = all_recordings_list_index;
    egg_structure(rec_index).egg_time = all_eggs_listtime;
    egg_structure(rec_index).egg_index = all_eggs_listindex;
    egg_structure(rec_index).ROI_used = all_ROI_index;
    
    
    
    
    corr_structure(rec_index).recording = recordings_list(rec_index);
    corr_structure(rec_index).recording_index = rec_index;
    corr_structure(rec_index).ROI_used = ROI_num(rec_index);
    
    
    corr_structure(rec_index).corr_to_vel = nanmean(corr_to_vel,1)';
    corr_structure(rec_index).time_base_corr_to_vel = time_base_corr_to_vel;
    
    
    % figure; hold on;  title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    % set(gca,'TickDir','out');
    % plot(time_base_corr_to_vel,corr_to_vel,'color',[.75,.75,.75]);
    % plot(time_base_corr_to_vel,nanmean(corr_to_vel),'-k');
    % ylabel('corr coeff');
    % xlabel('shift signal, correlate to smoothed velocity'); %since peak is neg, patch comes first. however...
    % axis manual;
    % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
    corr_structure(rec_index).corr_to_body = nanmean(corr_to_body,1)';
    corr_structure(rec_index).time_base_corr_to_body = time_base_corr_to_body;
    
    corr_structure(rec_index).corr_to_body_x = nanmean(corr_to_body_x,1)';
    corr_structure(rec_index).corr_to_body_y = nanmean(corr_to_body_y,1)';
    
    % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    % set(gca,'TickDir','out');
    % plot(time_base_corr_to_body,corr_to_body,'color',[.75,.75,.75]);
    % plot(time_base_corr_to_body,nanmean(corr_to_body),'-k');
    % ylabel('corr coeff');
    % xlabel('shift signal, correlate to abdomen length'); %since peak is neg, patch comes first. however...
    % axis manual;
    % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
    corr_structure(rec_index).corr_to_body_only = nanmean(corr_to_body_only,1)';
    corr_structure(rec_index).corr_to_body_only_x = nanmean(corr_to_body_only_x,1)';
    corr_structure(rec_index).corr_to_body_only_y = nanmean(corr_to_body_only_y,1)';
    
    % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    % set(gca,'TickDir','out');
    % plot(time_base_corr_to_body,corr_to_body_only,'color',[.75,.75,.75]);
    % plot(time_base_corr_to_body,nanmean(corr_to_body_only),'-k');
    % ylabel('corr coeff');
    % xlabel('shift signal, correlate to abdomen length L3 to tip'); %since peak is neg, patch comes first. however...
    % axis manual;
    % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
    corr_structure(rec_index).corr_to_body_path = nanmean(corr_to_body_path,1)';
    % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    % set(gca,'TickDir','out');
    % plot(time_base_corr_to_body,corr_to_body_path,'color',[.75,.75,.75]);
    % plot(time_base_corr_to_body,nanmean(corr_to_body_path),'-k');
    % ylabel('corr coeff');
    % xlabel('shift signal, correlate to abdomen path length L3 to tip'); %since peak is neg, patch comes first. however...
    % axis manual;
    % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
    corr_structure(rec_index).corr_to_angle = nanmean(corr_to_angle,1)';
    corr_structure(rec_index).time_base_corr_to_angle = time_base_corr_to_angle;
    
    % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    % set(gca,'TickDir','out');
    % plot(time_base_corr_to_angle,corr_to_angle,'color',[.75,.75,.75]);
    % plot(time_base_corr_to_angle,nanmean(corr_to_angle),'-k');
    % ylabel('corr coeff');
    % xlabel('shift signal, correlate to abdomen angle'); %since peak is neg, patch comes first. however...
    % axis manual;
    % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
    corr_structure(rec_index).corr_to_angle_only = nanmean(corr_to_angle_only,1)';
    % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    % set(gca,'TickDir','out');
    % plot(time_base_corr_to_angle,corr_to_angle_only,'color',[.75,.75,.75]);
    % plot(time_base_corr_to_angle,nanmean(corr_to_angle_only),'-k');
    % ylabel('corr coeff');
    % xlabel('shift signal, correlate to abdomen angle L3 to tip'); %since peak is neg, patch comes first. however...
    % axis manual;
    % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
    corr_structure(rec_index).corr_to_sucrose = nanmean(corr_to_sucrose,1)';
    corr_structure(rec_index).time_base_corr_to_sucrose = time_base_corr_to_sucrose;
    
    % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    % set(gca,'TickDir','out');
    % plot(time_base_corr_to_sucrose,corr_to_sucrose,'color',[.75,.75,.75]);
    % plot(time_base_corr_to_sucrose,nanmean(corr_to_sucrose),'-k');
    % ylabel('corr coeff');
    % xlabel('shift signal, correlate to sucrose concentration'); %since peak is neg, patch comes first. however...
    % axis manual;
    % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
    corr_structure(rec_index).corr_to_sucrose_t = nanmean(corr_to_sucrose_t,1)';
    corr_structure(rec_index).time_base_corr_to_sucrose_t = time_base_corr_to_sucrose_t;
    
    
    % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    % set(gca,'TickDir','out');
    % plot(time_base_corr_to_sucrose_t,corr_to_sucrose_t,'color',[.75,.75,.75]);
    % plot(time_base_corr_to_sucrose_t,nanmean(corr_to_sucrose_t),'-k');
    % ylabel('corr coeff');
    % xlabel('shift signal, correlate to sucrose thresholded concentration'); %since peak is neg, patch comes first. however...
    % axis manual;
    % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
    corr_structure(rec_index).corr_to_trans = nanmean(corr_to_trans,1)';
    corr_structure(rec_index).time_base_corr_to_trans = time_base_corr_to_trans;
    
    
    % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    % set(gca,'TickDir','out');
    % plot(time_base_corr_to_trans,corr_to_trans,'color',[.75,.75,.75]);
    % plot(time_base_corr_to_trans,nanmean(corr_to_trans),'-k');
    % ylabel('corr coeff');
    % xlabel('shift signal, correlate to signed transitions (10 sec around trans is 1 or -1)'); %since peak is neg, patch comes first. however...
    % axis manual;
    % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
    corr_structure(rec_index).corr_to_absvaltrans = nanmean(corr_to_absvaltrans,1)';
    corr_structure(rec_index).time_base_corr_to_absvaltrans = time_base_corr_to_absvaltrans;
    
    
    % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    % set(gca,'TickDir','out');
    % plot(time_base_corr_to_absvaltrans,corr_to_absvaltrans,'color',[.75,.75,.75]);
    % plot(time_base_corr_to_absvaltrans,nanmean(corr_to_absvaltrans),'-k');
    % ylabel('corr coeff');
    % xlabel('shift signal, correlate to abs val transitions (10 sec around trans is )'); %since peak is neg, patch comes first. however...
    % axis manual;
    % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
    corr_structure(rec_index).corr_to_egg = nanmean(corr_to_egg,1)';
    corr_structure(rec_index).time_base_corr_to_egg = time_base_corr_to_egg;
    
    % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    % set(gca,'TickDir','out');
    % plot(time_base_corr_to_egg,corr_to_egg,'color',[.75,.75,.75]);
    % plot(time_base_corr_to_egg,nanmean(corr_to_egg),'-k');
    % ylabel('corr coeff');
    % xlabel('shift signal, correlate to eggs'); %since peak is neg, patch comes first. however...
    % axis manual;
    % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
    corr_structure(rec_index).corr_to_plainegg = nanmean(corr_to_plainegg,1)';
    corr_structure(rec_index).time_base_corr_to_plainegg = time_base_corr_to_plainegg;
    
    % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    % set(gca,'TickDir','out');
    % plot(time_base_corr_to_plainegg,corr_to_plainegg,'color',[.75,.75,.75]);
    % plot(time_base_corr_to_plainegg,nanmean(corr_to_plainegg),'-k');
    % ylabel('corr coeff');
    % xlabel('shift signal, correlate to plain eggs'); %since peak is neg, patch comes first. however...
    % axis manual;
    % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
    corr_structure(rec_index).corr_to_sucroseegg = nanmean(corr_to_sucroseegg,1)';
    corr_structure(rec_index).time_base_corr_to_sucroseegg = time_base_corr_to_sucroseegg;
    
    % figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
    % set(gca,'TickDir','out');
    % plot(time_base_corr_to_sucroseegg,corr_to_sucroseegg,'color',[.75,.75,.75]);
    % plot(time_base_corr_to_sucroseegg,nanmean(corr_to_sucroseegg),'-k');
    % ylabel('corr coeff');
    % xlabel('shift signal, correlate to sucrose eggs'); %since peak is neg, patch comes first. however...
    % axis manual;
    % line('Xdata',[0,0],'YData',[-100,100],'Color','c');
    
end
corr_to_vel = [corr_structure(1:end).corr_to_vel];

figure; hold on;  %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_vel')', nanstd(corr_to_vel')', time_base_corr_to_vel', [.5 .5 .5])
plot(time_base_corr_to_vel,corr_to_vel,'color',[.75,.75,.75]);
plot(time_base_corr_to_vel,nanmean(corr_to_vel,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to smoothed velocity'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.3,.3])


corr_to_body = [corr_structure(1:end).corr_to_body];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_body')', nanstd(corr_to_body')', time_base_corr_to_body', [.5 .5 .5])
plot(time_base_corr_to_body,corr_to_body,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen length'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.3,.3])

corr_to_body_x = [corr_structure(1:end).corr_to_body_x];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_body_x')', nanstd(corr_to_body_x')', time_base_corr_to_body', [.5 .5 .5])
plot(time_base_corr_to_body,corr_to_body_x,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body_x,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen length, X only'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.3,.3])


corr_to_body_y = [corr_structure(1:end).corr_to_body_y];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_body_y')', nanstd(corr_to_body_y')', time_base_corr_to_body', [.5 .5 .5])
plot(time_base_corr_to_body,corr_to_body_y,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body_y,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen length, Y only'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.3,.3])

corr_to_body_only = [corr_structure(1:end).corr_to_body_only];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_body_only')', nanstd(corr_to_body_only')', time_base_corr_to_body', [.5 .5 .5])
plot(time_base_corr_to_body,corr_to_body_only,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body_only,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen length L3 to tip'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.3,.3])



corr_to_body_only_x = [corr_structure(1:end).corr_to_body_only_x];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_body_only_x')', nanstd(corr_to_body_only_x')', time_base_corr_to_body', [.5 .5 .5])
plot(time_base_corr_to_body,corr_to_body_only_x,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body_only_x,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen length L3 to tip, X only'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.3,.3])


corr_to_body_only_y = [corr_structure(1:end).corr_to_body_only_y];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_body_only_y')', nanstd(corr_to_body_only_y')', time_base_corr_to_body', [.5 .5 .5])
plot(time_base_corr_to_body,corr_to_body_only_y,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body_only_y,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen length L3 to tip, Y only'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.3,.3])


corr_to_body_path = [corr_structure(1:end).corr_to_body_path];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_body_path')', nanstd(corr_to_body_path')', time_base_corr_to_body', [.5 .5 .5])
plot(time_base_corr_to_body,corr_to_body_path,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body_path,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen path length L3 to tip'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.3,.3])

corr_to_angle = [corr_structure(1:end).corr_to_angle];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_angle')', nanstd(corr_to_angle')', time_base_corr_to_angle', [.5 .5 .5])
plot(time_base_corr_to_angle,corr_to_angle,'color',[.75,.75,.75]);
plot(time_base_corr_to_angle,nanmean(corr_to_angle,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen angle'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.3,.3])

corr_to_angle_only = [corr_structure(1:end).corr_to_angle_only];

% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_angle,corr_to_angle_only,'color',[.75,.75,.75]);
% plot(time_base_corr_to_angle,nanmean(corr_to_angle_only),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen angle L3 to tip'); %since peak is neg, patch comes first. however...
% axis manual;
%line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
%

corr_to_sucrose = [corr_structure(1:end).corr_to_sucrose];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_sucrose')', nanstd(corr_to_sucrose')', time_base_corr_to_sucrose', [.5 .5 .5])
plot(time_base_corr_to_sucrose,corr_to_sucrose,'color',[.75,.75,.75]);
plot(time_base_corr_to_sucrose,nanmean(corr_to_sucrose,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to sucrose concentration'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.3,.3])

corr_to_sucrose_t = [corr_structure(1:end).corr_to_sucrose_t];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_sucrose_t')', nanstd(corr_to_sucrose_t')', time_base_corr_to_sucrose_t', [.5 .5 .5])
plot(time_base_corr_to_sucrose_t,corr_to_sucrose_t,'color',[.75,.75,.75]);
plot(time_base_corr_to_sucrose_t,nanmean(corr_to_sucrose_t,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to sucrose thresholded concentration'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.3,.3])

corr_to_trans = [corr_structure(1:end).corr_to_trans];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_trans')', nanstd(corr_to_trans')', time_base_corr_to_trans', [.5 .5 .5])
plot(time_base_corr_to_trans,corr_to_trans,'color',[.75,.75,.75]);
plot(time_base_corr_to_trans,nanmean(corr_to_trans,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to signed transitions (10 sec around trans is 1 or -1)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.3,.3])

corr_to_absvaltrans = [corr_structure(1:end).corr_to_absvaltrans];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_absvaltrans')', nanstd(corr_to_absvaltrans')', time_base_corr_to_absvaltrans', [.5 .5 .5])
plot(time_base_corr_to_absvaltrans,corr_to_absvaltrans,'color',[.75,.75,.75]);
plot(time_base_corr_to_absvaltrans,nanmean(corr_to_absvaltrans,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abs val transitions (10 sec around trans is )'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.3,.3])

corr_to_egg = [corr_structure(1:end).corr_to_egg];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_egg')', nanstd(corr_to_egg')', time_base_corr_to_egg', [.5 .5 .5])
plot(time_base_corr_to_egg,corr_to_egg,'color',[.75,.75,.75]);
plot(time_base_corr_to_egg,nanmean(corr_to_egg,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to eggs'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.4,.4])

corr_to_plainegg = [corr_structure(1:end).corr_to_plainegg];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_plainegg')', nanstd(corr_to_plainegg')', time_base_corr_to_plainegg', [.5 .5 .5])
plot(time_base_corr_to_plainegg,corr_to_plainegg,'color',[.75,.75,.75]);
plot(time_base_corr_to_plainegg,nanmean(corr_to_plainegg,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to plain eggs'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.4,.4])


corr_to_sucroseegg = [corr_structure(1:end).corr_to_sucroseegg];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
patch_errorbar(nanmean(corr_to_sucroseegg')', nanstd(corr_to_sucroseegg')', time_base_corr_to_sucroseegg', [.5 .5 .5])
plot(time_base_corr_to_sucroseegg,corr_to_sucroseegg,'color',[.75,.75,.75]);
plot(time_base_corr_to_sucroseegg,nanmean(corr_to_sucroseegg,2),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to sucrose eggs'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.4,.4])


% plotting egg triggered averages
data_around_egg_vel = [egg_structure(1:end).data_around_egg_vel];


data_around_egg_body = [egg_structure(1:end).data_around_egg_body];
data_around_egg_body_x = [egg_structure(1:end).data_around_egg_body_x];
data_around_egg_body_y = [egg_structure(1:end).data_around_egg_body_y];

data_around_egg_path = [egg_structure(1:end).data_around_egg_path];
data_around_egg_body_only = [egg_structure(1:end).data_around_egg_body_only];
data_around_egg_body_only_x = [egg_structure(1:end).data_around_egg_body_only_x];
data_around_egg_body_only_y = [egg_structure(1:end).data_around_egg_body_only_y];

data_around_egg_angle = [egg_structure(1:end).data_around_egg_angle];
data_around_egg_angle_only = [egg_structure(1:end).data_around_egg_angle_only];
data_around_egg = [egg_structure(1:end).data_around_egg];

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    patch_errorbar(nanmean(data_around_egg_vel')', nanstd(data_around_egg_vel')', time_base_egg', [.5 .5 .5])
    plot(time_base_egg,data_around_egg_vel,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_vel,2),'-k');
end
ylabel('vel, mm/sec (smth 1s)');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    patch_errorbar(nanmean(data_around_egg_body')', nanstd(data_around_egg_body')', time_base_egg', [.5 .5 .5])
    plot(time_base_egg,data_around_egg_body,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_body,2),'-k');
end
ylabel('neck to tip length');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    patch_errorbar(nanmean(data_around_egg_body_x')', nanstd(data_around_egg_body_x')', time_base_egg', [.5 .5 .5])
    plot(time_base_egg,data_around_egg_body_x,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_body_x,2),'-k');
end
ylabel('neck to tip length, X only');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');



figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    patch_errorbar(nanmean(data_around_egg_body_y')', nanstd(data_around_egg_body_y')', time_base_egg', [.5 .5 .5])
    plot(time_base_egg,data_around_egg_body_y,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_body_y,2),'-k');
end
ylabel('neck to tip length, Y only');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');



figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    patch_errorbar(nanmean(data_around_egg_path')', nanstd(data_around_egg_path')', time_base_egg', [.5 .5 .5])
    plot(time_base_egg,data_around_egg_path,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_path,2),'-k');
end
ylabel('neck to tip path length');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    patch_errorbar(nanmean(data_around_egg_body_only')', nanstd(data_around_egg_body_only')', time_base_egg', [.5 .5 .5])
    plot(time_base_egg,data_around_egg_body_only,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_body_only,2),'-k');
end
ylabel('L3 to tip length');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    patch_errorbar(nanmean(data_around_egg_body_only_x')', nanstd(data_around_egg_body_only_x')', time_base_egg', [.5 .5 .5])
    plot(time_base_egg,data_around_egg_body_only_x,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_body_only_x,2),'-k');
end
ylabel('L3 to tip length, X only');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    patch_errorbar(nanmean(data_around_egg_body_only_y')', nanstd(data_around_egg_body_only_y')', time_base_egg', [.5 .5 .5])
    plot(time_base_egg,data_around_egg_body_only_y,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_body_only_y,2),'-k');
end
ylabel('L3 to tip length, Y only');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');



figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    patch_errorbar(nanmean(data_around_egg_angle')', nanstd(data_around_egg_angle')', time_base_egg', [.5 .5 .5])
    plot(time_base_egg,data_around_egg_angle,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_angle,2),'-k');
end
ylabel('neck to tip angle');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    patch_errorbar(nanmean(data_around_egg')', nanstd(data_around_egg')', time_base_egg', [.5 .5 .5])
    plot(time_base_egg,data_around_egg,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg,2),'-k');
end
ylabel('signal around egg');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'ylim',[-.5,1.5])

end