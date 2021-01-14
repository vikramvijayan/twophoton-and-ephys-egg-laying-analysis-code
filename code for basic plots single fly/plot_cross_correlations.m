function [out] = plot_cross_correlations(recording, ROI_num, filter_out_chrimson_data, filter_chrimson_eggs)

%% make sure to set the time interval for the calculations (triggered averages)
time_interval = [-240:.1:240];

%% make sure to set the time interval for cross correaltions
time_interval_corr = [-120:.1:120];


% set filter chrimson data to 1 if yu want to remove all data during pulses
% (as well as 100 sec after the pulse) check code for exact number

% set filter chrimson eggs to remove all eggs in the same interval as above
% 0 = no filter
% 1 - filter out chrimson eggs
% 2 - only use chrimson eggs

% set ROI num to -1 if you want tu use CH1 patch
% if patch is 1 then use CH1
% if patch is 0 then use CH1 but remove pulses (and 100 seconds after)

  if(ROI_num ==-1)
        patch = 1;
    else
        patch = 0;
  end
    
if(patch == 1)
    recording.tseries = 1;
end
%%%%%%%%%%%
%% plotting correlations or time locked average


data_to_plain   = [];
data_to_sucrose = [];

data_to_sucrose_control = [];
data_to_plain_control = [];
time_base_trans = [];

data_around_egg_vel = [];
data_around_egg_vel_noave = [];
data_around_egg_speed_2hz_hold = [];

    

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
corr_to_vel_noave = [];
corr_to_speed_2hz_hold = [];

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
  wheelvelocity_noave = (pi*7*-25).*[0; diff(unwrap(recording.movie1.filtered_wheel))]./ (2*pi); %mm/sec 7mm radium
  
% process more of the DLC movie
%[recording] = processes_more_DLC_variables(recording);
abd_only_length_x = recording.movie1.abd_x_L3tip;
bodylength_x =  recording.movie1.abd_x_neck_tip;
abd_only_length_y = recording.movie1.abd_y_L3tip;
bodylength_y =  recording.movie1.abd_y_neck_tip;

prob_x =  recording.movie1.prob_x;
prob_y =  recording.movie1.prob_y;

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
time_base_egg = time_base_to_return;

[time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,plain_eggt, time_interval);
data_around_plain_egg_sub = [data_around_plain_egg_sub; data_to_average_interp];

[time_base_to_return, data_to_average_interp]  = average_around_event(sucrose,movtime,sucrose_eggt, time_interval);
data_around_sucrose_egg_sub = [data_around_sucrose_egg_sub; data_to_average_interp];

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


data_around_sucrose_egg_sub = repmat(data_around_sucrose_egg_sub,length(recording.tseries),1);
data_around_egg_sub = repmat(data_around_egg_sub,length(recording.tseries),1);
data_around_plain_egg_sub = repmat(data_around_plain_egg_sub,length(recording.tseries),1);

data_around_sucrose_egg_sub2 = repmat(data_around_sucrose_egg_sub2,length(recording.tseries),1);
data_around_egg_sub2 = repmat(data_around_egg_sub2,length(recording.tseries),1);
data_around_plain_egg_sub2 = repmat(data_around_plain_egg_sub2,length(recording.tseries),1);

for m = 1:1:length(recording.tseries)
    if(patch == 0)
        signal    = recording.tseries(m).df_over_f(ROI_num,:);
        signaltime = recording.tseries(m).Time_s;
        
        if(filter_out_chrimson_data)
            tseries_no_laser = interp1(recording.abf.Time_s, recording.abf.no_laser, recording.tseries(m).Time_s,'previous');
            signal = signal.*tseries_no_laser';
        end
    end
    
    % for patching data, the signal changes
    if(patch == 1)
        %signal    = recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
        signaltime = recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
        
         signal    = recording.abf.CH1_patch(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000))-13;
         signal    = recording.abf.CH1_patch_spikes_conv_area_rect(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));

        if(filter_out_chrimson_data)
            signal    = recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)).*recording.abf.no_laser(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
        end
        
    end
    
    

    
    % look at signal at eggs
    [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,all_eggt, time_interval);
    data_around_egg = [data_around_egg; data_to_average_interp];
    time_base_egg = time_base_to_return;
    [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,plain_eggt, time_interval);
    data_around_plain_egg = [data_around_plain_egg; data_to_average_interp];
    [time_base_to_return, data_to_average_interp]  = average_around_event(signal,signaltime,sucrose_eggt, time_interval);
    data_around_sucrose_egg = [data_around_sucrose_egg; data_to_average_interp];
    

    
    % look at signal correlated to velocity (smoothed)
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, wheelvelocity, signaltime, movtime, time_interval_corr);
    corr_to_vel = [corr_to_vel; out_corr];
    time_base_corr_to_vel = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, wheelvelocity_noave, signaltime, movtime, time_interval_corr);
    corr_to_vel_noave = [corr_to_vel_noave; out_corr];
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, speed_2hz_hold, signaltime, movtime, time_interval_corr);
    corr_to_speed_2hz_hold = [corr_to_speed_2hz_hold; out_corr];
    
    % look at signal correlated to body position
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, problength, signaltime, movtime, time_interval_corr);
    corr_to_prob = [corr_to_prob; out_corr];
    time_base_corr_to_prob = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, prob_x, signaltime, movtime, time_interval_corr);
    corr_to_prob_x = [corr_to_prob_x; out_corr];
    time_base_corr_to_prob = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, prob_y, signaltime, movtime, time_interval_corr);
    corr_to_prob_y = [corr_to_prob_y; out_corr];
    time_base_corr_to_prob = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodylength, signaltime, movtime, time_interval_corr);
    corr_to_body = [corr_to_body; out_corr];
    time_base_corr_to_body = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodylength_x, signaltime, movtime, time_interval_corr);
    corr_to_body_x = [corr_to_body_x; out_corr];
    time_base_corr_to_body = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodylength_y, signaltime, movtime, time_interval_corr);
    corr_to_body_y = [corr_to_body_y; out_corr];
    time_base_corr_to_body = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, bodyangle, signaltime, movtime, time_interval_corr);
    corr_to_angle = [corr_to_angle; out_corr];
    time_base_corr_to_angle = out_corr_shift_ofvec1;
    
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
    time_base_corr_to_sucrose = out_corr_shift_ofvec1;
    
    % look at signal correlated to substrate threshold
    suc_thr = sucrose;
    suc_thr(suc_thr>0) = 1;
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, suc_thr, signaltime, movtime, time_interval_corr);
    corr_to_sucrose_t = [corr_to_sucrose_t; out_corr];
    time_base_corr_to_sucrose_t = out_corr_shift_ofvec1;
    
    % look at signal correlated to egg
    tmp_ar = zeros(1,length(movtime));
    tmp_ar(all_egg) = 1;
    tmp_ar = smooth(tmp_ar,25*10);
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, time_interval_corr);
    corr_to_egg = [corr_to_egg; out_corr];
    time_base_corr_to_egg = out_corr_shift_ofvec1;
    
    % look at signal correlated to plain egg
    tmp_ar = zeros(1,length(movtime));
    tmp_ar(plain_egg) = 1;
    tmp_ar = smooth(tmp_ar,25*10);
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, time_interval_corr);
    corr_to_plainegg = [corr_to_plainegg; out_corr];
    time_base_corr_to_plainegg = out_corr_shift_ofvec1;
    
    % look at signal correlated to sucrose egg
    tmp_ar = zeros(1,length(movtime));
    tmp_ar(sucrose_egg) = 1;
    tmp_ar = smooth(tmp_ar,25*10);
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, time_interval_corr);
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
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal,[0;transitionssigned], signaltime, movtime, time_interval_corr);
    corr_to_trans = [corr_to_trans; out_corr];
    time_base_corr_to_trans = out_corr_shift_ofvec1;
    
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal,[0;transitionsabsval], signaltime, movtime, time_interval_corr);
    corr_to_absvaltrans = [corr_to_absvaltrans; out_corr];
    time_base_corr_to_absvaltrans = out_corr_shift_ofvec1;
end

figure; hold on;  title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_vel,corr_to_vel,'color',[.75,.75,.75]);
plot(time_base_corr_to_vel,nanmean(corr_to_vel),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to smoothed velocity'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on;  title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_vel,corr_to_speed_2hz_hold,'color',[.75,.75,.75]);
plot(time_base_corr_to_vel,nanmean(corr_to_speed_2hz_hold),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to 2hz hold speed'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on;  title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_vel,corr_to_vel_noave,'color',[.75,.75,.75]);
plot(time_base_corr_to_vel,nanmean(corr_to_vel_noave),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to velocity'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_body,corr_to_prob,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_prob),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to proboscis length'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_body,corr_to_prob_x,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_prob_x),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to proboscis length, X only'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');


figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_body,corr_to_prob_y,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_prob_y),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to proboscis length, Y only'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');


figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_body,corr_to_body,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen length'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_body,corr_to_body_x,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body_x),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen length, X only'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');


figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_body,corr_to_body_y,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body_y),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen length, Y only'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_body,corr_to_body_only,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body_only),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen length L3 to tip'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_body,corr_to_body_only_x,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body_only_x),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen length L3 to tip, X only'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');


figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_body,corr_to_body_only_y,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body_only_y),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen length L3 to tip, Y only'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_body,corr_to_body_path,'color',[.75,.75,.75]);
plot(time_base_corr_to_body,nanmean(corr_to_body_path),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen path length L3 to tip'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_angle,corr_to_angle,'color',[.75,.75,.75]);
plot(time_base_corr_to_angle,nanmean(corr_to_angle),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen angle'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_angle,corr_to_angle_only,'color',[.75,.75,.75]);
plot(time_base_corr_to_angle,nanmean(corr_to_angle_only),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abdomen angle L3 to tip'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_sucrose,corr_to_sucrose,'color',[.75,.75,.75]);
plot(time_base_corr_to_sucrose,nanmean(corr_to_sucrose),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to sucrose concentration'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_sucrose_t,corr_to_sucrose_t,'color',[.75,.75,.75]);
plot(time_base_corr_to_sucrose_t,nanmean(corr_to_sucrose_t),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to sucrose thresholded concentration'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_trans,corr_to_trans,'color',[.75,.75,.75]);
plot(time_base_corr_to_trans,nanmean(corr_to_trans),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to signed transitions (10 sec around trans is 1 or -1)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_absvaltrans,corr_to_absvaltrans,'color',[.75,.75,.75]);
plot(time_base_corr_to_absvaltrans,nanmean(corr_to_absvaltrans),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to abs val transitions (10 sec around trans is )'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_egg,corr_to_egg,'color',[.75,.75,.75]);
plot(time_base_corr_to_egg,nanmean(corr_to_egg),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to eggs'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_plainegg,corr_to_plainegg,'color',[.75,.75,.75]);
plot(time_base_corr_to_plainegg,nanmean(corr_to_plainegg),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to plain eggs'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_sucroseegg,corr_to_sucroseegg,'color',[.75,.75,.75]);
plot(time_base_corr_to_sucroseegg,nanmean(corr_to_sucroseegg),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to sucrose eggs'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

% plot signal around egg

average1 = nanmean(data_around_plain_egg,1);
average2 = nanmean(data_around_sucrose_egg,1);
average3 = nanmean(data_around_egg,1);

r1 = 0;
r2 = 0;
r3 = 0;
clrzz = [ 107./255,174./255,214./255; [33 113 181]./255; [33 113 181]./255; 8./255,48./255,107./255; 8./255,48./255,107./255; 8./255,48./255,107./255];

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average signal around eggs laid on plain');

if(~isempty(data_around_plain_egg))
    [r,c] = size(data_around_plain_egg);
    for sizes = 1:1:r
        for lineseg = 1 : (length(data_around_plain_egg)-1)
            if(~isnan(data_around_plain_egg_sub(sizes,lineseg)))
                line('XData',time_base_egg(lineseg:lineseg+1), 'YData',data_around_plain_egg(sizes,lineseg:lineseg+1), 'Color',clrzz(data_around_plain_egg_sub(sizes,lineseg)./100+1,:));
            end
        end
    end
end

if(~isempty(data_around_plain_egg))
    hold on;
    plot(time_base_egg,average1,'k'); % egg plain
    data_around_plain_egg(~any(~isnan(data_around_plain_egg), 2),:)=[];
    [r1] = sum(~isnan(data_around_plain_egg(:,2401)));
end

ylabel([num2str(r1) ' eggs laid on plain']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average signal around eggs laid on sucrose');

if(~isempty(data_around_sucrose_egg))
    [r,c] = size(data_around_sucrose_egg);
    for sizes = 1:1:r
        for lineseg = 1 : (length(data_around_sucrose_egg)-1)
            if(~isnan(data_around_sucrose_egg_sub(sizes,lineseg)))
                line('XData',time_base_egg(lineseg:lineseg+1), 'YData',data_around_sucrose_egg(sizes,lineseg:lineseg+1), 'Color',clrzz(data_around_sucrose_egg_sub(sizes,lineseg)./100+1,:));
            end
        end
    end
end

if(~isempty(data_around_sucrose_egg))
    hold on;
    plot(time_base_egg,average2,'k');% egg sucrose
    data_around_sucrose_egg(~any(~isnan(data_around_sucrose_egg), 2),:)=[];
    [r2] = sum(~isnan(data_around_sucrose_egg(:,2401)));
end

ylabel([num2str(r2) ' eggs laid on sucrose']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average signal around eggs');

if(~isempty(data_around_egg))
    hold on;
    
    [r,c] = size(data_around_egg);
    for sizes = 1:1:r
        for lineseg = 1 : (length(data_around_egg)-1)
            line('XData',time_base_egg(lineseg:lineseg+1), 'YData',data_around_egg(sizes,lineseg:lineseg+1), 'Color',[.75 .75 .75]);
        end
    end
    
    plot(time_base_egg,average3,'Color','k'); % egg
    %data_around_egg(any(isnan(data_around_egg), 2),:)=[];
    %[r3,c] = size(data_around_egg);
end

ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

% plotting average of various body length around egg
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average vel mm/sec (1s smoth)');
average = nanmean(data_around_egg_vel,1);

if(~isempty(data_around_egg_vel))
    hold on;
    
    [r,c] = size(data_around_egg_vel);
    for sizes = 1:1:r
        plot(time_base_egg, data_around_egg_vel(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_egg,average,'Color','k'); % egg
end
ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

% plotting average of various prob length around egg
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average proboscis length around eggs');
average = nanmean(data_around_egg_prob,1);

if(~isempty(data_around_egg_prob))
    hold on;
    
    [r,c] = size(data_around_egg_prob);
    for sizes = 1:1:r
        plot(time_base_egg, data_around_egg_prob(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_egg,average,'Color','k'); % egg
end
ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

% plotting average of various prob length around egg
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average proboscis length around eggs, X only');
average = nanmean(data_around_egg_prob_x,1);

if(~isempty(data_around_egg_prob_x))
    hold on;
    
    [r,c] = size(data_around_egg_prob_x);
    for sizes = 1:1:r
        plot(time_base_egg, data_around_egg_prob_x(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_egg,average,'Color','k'); % egg
end
ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

% plotting average of various prob length around egg
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average proboscis length around eggs, Y only');
average = nanmean(data_around_egg_prob_y,1);

if(~isempty(data_around_egg_prob_y))
    hold on;
    
    [r,c] = size(data_around_egg_prob_y);
    for sizes = 1:1:r
        plot(time_base_egg, data_around_egg_prob_y(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_egg,average,'Color','k'); % egg
end
ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

% plotting average of various body length around egg
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body length around eggs');
average = nanmean(data_around_egg_body,1);

if(~isempty(data_around_egg_body))
    hold on;
    
    [r,c] = size(data_around_egg_body);
    for sizes = 1:1:r
        plot(time_base_egg, data_around_egg_body(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_egg,average,'Color','k'); % egg
end
ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

% plotting average of various body length around egg
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body length around eggs, X only');
average = nanmean(data_around_egg_body_x,1);

if(~isempty(data_around_egg_body_x))
    hold on;
    
    [r,c] = size(data_around_egg_body_x);
    for sizes = 1:1:r
        plot(time_base_egg, data_around_egg_body_x(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_egg,average,'Color','k'); % egg
end
ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');



% plotting average of various body length around egg
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body length around eggs, Y only');
average = nanmean(data_around_egg_body_y,1);

if(~isempty(data_around_egg_body_y))
    hold on;
    
    [r,c] = size(data_around_egg_body_y);
    for sizes = 1:1:r
        plot(time_base_egg, data_around_egg_body_y(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_egg,average,'Color','k'); % egg
end
ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

% plotting average of various body path metrics around egg
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body path length around eggs');
average = nanmean(data_around_egg_path,1);

if(~isempty(data_around_egg_path))
    hold on;
    
    [r,c] = size(data_around_egg_path);
    for sizes = 1:1:r
        plot(time_base_egg, data_around_egg_path(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_egg,average,'Color','k'); % egg
end
ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

% plotting average of various body path metrics around egg
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body length L3 to tip around eggs');
average = nanmean(data_around_egg_body_only,1);

if(~isempty(data_around_egg_body_only))
    hold on;
    
    [r,c] = size(data_around_egg_body_only);
    for sizes = 1:1:r
        plot(time_base_egg, data_around_egg_body_only(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_egg,average,'Color','k'); % egg
end
ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

% plotting average of various body path metrics around egg
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body length L3 to tip around eggs, X only');
average = nanmean(data_around_egg_body_only_x,1);

if(~isempty(data_around_egg_body_only_x))
    hold on;
    
    [r,c] = size(data_around_egg_body_only_x);
    for sizes = 1:1:r
        plot(time_base_egg, data_around_egg_body_only_x(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_egg,average,'Color','k'); % egg
end
ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

% plotting average of various body path metrics around egg
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body length L3 to tip around eggs, Y only');
average = nanmean(data_around_egg_body_only_y,1);

if(~isempty(data_around_egg_body_only_y))
    hold on;
    
    [r,c] = size(data_around_egg_body_only_y);
    for sizes = 1:1:r
        plot(time_base_egg, data_around_egg_body_only_y(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_egg,average,'Color','k'); % egg
end
ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');



% plotting average of various body path metrics around egg
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body angle around eggs');
average = nanmean(data_around_egg_angle,1);

if(~isempty(data_around_egg_angle))
    hold on;
    
    [r,c] = size(data_around_egg_angle);
    for sizes = 1:1:r
        plot(time_base_egg, data_around_egg_angle(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_egg,average,'Color','k'); % egg
end
ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

% plotting average of various body path metrics around egg
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body angle L3 to tip around eggs');
average = nanmean(data_around_egg_angle_only,1);

if(~isempty(data_around_egg_angle_only))
    hold on;
    
    [r,c] = size(data_around_egg_angle_only);
    for sizes = 1:1:r
        plot(time_base_egg, data_around_egg_angle_only(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_egg,average,'Color','k'); % egg
end
ylabel([num2str(r1+r2) ' total eggs']);
axis manual;
line('Xdata',[0,0],'YData',[-100,100],'Color','c');

%%%%%%%%%%

end

