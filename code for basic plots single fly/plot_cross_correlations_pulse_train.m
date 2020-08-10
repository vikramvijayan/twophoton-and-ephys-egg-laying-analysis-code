function [out] = plot_cross_correlations_pulse_train(recording, ROI_num)

% set ROI num to -1 if you want tu use CH1 patch
% if patch is 1 then use CH1
% if patch is 0 then use CH1 but remove pulses (and 100 seconds after)

  if(ROI_num ==-1)
        patch = 1;
    else
        patch = 0;
  end
%%%%%%%%%%%
%% plotting correlations or time locked average


data_to_plain   = [];
data_to_sucrose = [];

data_to_sucrose_control = [];
data_to_plain_control = [];
time_base_trans = [];

data_around_pulse_vel = [];

data_around_pulse = [];
data_around_plain_pulse = [];
data_around_sucrose_pulse = [];
time_base_pulse = [];
data_around_pulse_sub = [];
data_around_plain_pulse_sub = [];
data_around_sucrose_pulse_sub = [];

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


% process more of the DLC movie
%[recording] = processes_more_DLC_variables(recording);
abd_only_length_x = recording.movie1.abd_x_L3tip;
bodylength_x =  recording.movie1.abd_x_neck_tip;
abd_only_length_y = recording.movie1.abd_y_L3tip;
bodylength_y =  recording.movie1.abd_y_neck_tip;

prob_x =  recording.movie1.prob_x;
prob_y =  recording.movie1.prob_y;

% find all pulse times this is a bit complicted since the pulses have pulses (the stimulation is at a specific Hz)
% this code will find the single beginning and end of a full pulse train (5 pulses of 1
% sec that are spaced 10 sec apart)

temp_diffs = diff(recording.abf.PWMlaser);
[a b] = find(temp_diffs > 2);
a=a+1;
temp_diffs2 = diff(a);
[a1 b1] = find(temp_diffs2 > 1000000);
recording.abf.pulse_ON_times = [a(1); a(a1+1)];

% pulse off is the last time it was ON
temp_diffs = diff(recording.abf.PWMlaser);
[a b] = find(temp_diffs < -2);
temp_diffs2 = diff(a);
[a1 b1] = find(temp_diffs2 > 1000000);
recording.abf.pulse_OFF_times = [a(a1); a(end)];

recording.movie1.pulse_ON_times_movie = [];
%convert abf times to movie index, either the first frame in movie after
%pulse is ON, or last frame before it is OFF
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

[time_base_to_return, data_to_average_interp]  = average_around_event(wheelvelocity,movtime,all_pulset, [-240:.1:240]);
data_around_pulse_vel = [data_around_pulse_vel; data_to_average_interp];


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
time_base_pulse = time_base_to_return;
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

for m = 1:1:length(recording.tseries)
    
  if(patch ==0)
        
        signal    = recording.tseries(m).df_over_f(ROI_num,:);
        signaltime = recording.tseries(m).Time_s;
        
    end
    
    wheelvelocity = (pi*7*-25).*smooth([0; diff(unwrap(recording.movie1.filtered_wheel))],25) ./ (2*pi); %mm/sec 7mm radium
    
      
    % for patching data, the signal changes
    if(patch ==1)
        signal    = recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
        signaltime = recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
    end
    
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
    
    % look at signal correlated to pulse [using the full pulse train of 45seconds])
    tmp_ar = zeros(1,length(movtime));
    tmp_ar(all_pulse:(all_pulse+45*25)) = 1;
    %tmp_ar = smooth(tmp_ar,25*10);
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, [-120:.1:120]);
    corr_to_pulse = [corr_to_pulse; out_corr];
    time_base_corr_to_pulse = out_corr_shift_ofvec1;
    
    % look at signal correlated to plain pulse
    tmp_ar = zeros(1,length(movtime));
    tmp_ar(plain_pulse:(plain_pulse+45*25)) = 1;
    %tmp_ar = smooth(tmp_ar,25*10);
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, [-120:.1:120]);
    corr_to_plainpulse = [corr_to_plainpulse; out_corr];
    time_base_corr_to_plainpulse = out_corr_shift_ofvec1;
    
    % look at signal correlated to sucrose pulse
    tmp_ar = zeros(1,length(movtime));
    tmp_ar(sucrose_pulse:(sucrose_pulse+45*25)) = 1;
    %tmp_ar = smooth(tmp_ar,25*10);
    [out_corr,out_corr_shift_ofvec1] = regular_xcorr(signal, tmp_ar, signaltime, movtime, [-120:.1:120]);
    corr_to_sucrosepulse = [corr_to_sucrosepulse; out_corr];
    time_base_corr_to_sucrosepulse = out_corr_shift_ofvec1;
    
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
% 
% figure; hold on;  title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_vel,corr_to_vel,'color',[.75,.75,.75]);
% plot(time_base_corr_to_vel,nanmean(corr_to_vel),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to smoothed velocity'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_body,corr_to_prob,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_prob),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to proboscis length'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_body,corr_to_prob_x,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_prob_x),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to proboscis length, X only'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_body,corr_to_prob_y,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_prob_y),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to proboscis length, Y only'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_body,corr_to_body,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen length'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_body,corr_to_body_x,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_x),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen length, X only'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_body,corr_to_body_y,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_y),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen length, Y only'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_body,corr_to_body_only,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_only),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen length L3 to tip'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_body,corr_to_body_only_x,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_only_x),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen length L3 to tip, X only'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_body,corr_to_body_only_y,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_only_y),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen length L3 to tip, Y only'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_body,corr_to_body_path,'color',[.75,.75,.75]);
% plot(time_base_corr_to_body,nanmean(corr_to_body_path),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen path length L3 to tip'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_angle,corr_to_angle,'color',[.75,.75,.75]);
% plot(time_base_corr_to_angle,nanmean(corr_to_angle),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen angle'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_angle,corr_to_angle_only,'color',[.75,.75,.75]);
% plot(time_base_corr_to_angle,nanmean(corr_to_angle_only),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abdomen angle L3 to tip'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_sucrose,corr_to_sucrose,'color',[.75,.75,.75]);
% plot(time_base_corr_to_sucrose,nanmean(corr_to_sucrose),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to sucrose concentration'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_sucrose_t,corr_to_sucrose_t,'color',[.75,.75,.75]);
% plot(time_base_corr_to_sucrose_t,nanmean(corr_to_sucrose_t),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to sucrose thresholded concentration'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_trans,corr_to_trans,'color',[.75,.75,.75]);
% plot(time_base_corr_to_trans,nanmean(corr_to_trans),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to signed transitions (10 sec around trans is 1 or -1)'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% 
% figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
% set(gca,'TickDir','out');
% plot(time_base_corr_to_absvaltrans,corr_to_absvaltrans,'color',[.75,.75,.75]);
% plot(time_base_corr_to_absvaltrans,nanmean(corr_to_absvaltrans),'-k');
% ylabel('corr coeff');
% xlabel('shift signal, correlate to abs val transitions (10 sec around trans is )'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_pulse,corr_to_pulse,'color',[.75,.75,.75]);
plot(time_base_corr_to_pulse,nanmean(corr_to_pulse),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to pulses'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_plainpulse,corr_to_plainpulse,'color',[.75,.75,.75]);
plot(time_base_corr_to_plainpulse,nanmean(corr_to_plainpulse),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to plain pulses'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
plot(time_base_corr_to_sucrosepulse,corr_to_sucrosepulse,'color',[.75,.75,.75]);
plot(time_base_corr_to_sucrosepulse,nanmean(corr_to_sucrosepulse),'-k');
ylabel('corr coeff');
xlabel('shift signal, correlate to sucrose pulses'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');

% plot signal around pulse

average1 = nanmean(data_around_plain_pulse,1);
average2 = nanmean(data_around_sucrose_pulse,1);
average3 = nanmean(data_around_pulse,1);

r1 = 0;
r2 = 0;
r3 = 0;
clrzz = [ 107./255,174./255,214./255; [33 113 181]./255; [33 113 181]./255; 8./255,48./255,107./255; 8./255,48./255,107./255; 8./255,48./255,107./255];

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average signal around pulses occuring on plain');

if(~isempty(data_around_plain_pulse))
    [r,c] = size(data_around_plain_pulse);
    for sizes = 1:1:r
        for lineseg = 1 : (length(data_around_plain_pulse)-1)
            if(~isnan(data_around_plain_pulse_sub(sizes,lineseg)))
                line('XData',time_base_pulse(lineseg:lineseg+1), 'YData',data_around_plain_pulse(sizes,lineseg:lineseg+1), 'Color',clrzz(data_around_plain_pulse_sub(sizes,lineseg)./100+1,:));
            end
        end
    end
end

if(~isempty(data_around_plain_pulse))
    hold on;
    plot(time_base_pulse,average1,'k'); % pulse plain
    data_around_plain_pulse(~any(~isnan(data_around_plain_pulse), 2),:)=[];
    [r1] = sum(~isnan(data_around_plain_pulse(:,2401)));
end

ylabel([num2str(r1) ' pulses occuring on plain']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average signal around pulses occuring on sucrose');

if(~isempty(data_around_sucrose_pulse))
    [r,c] = size(data_around_sucrose_pulse);
    for sizes = 1:1:r
        for lineseg = 1 : (length(data_around_sucrose_pulse)-1)
            if(~isnan(data_around_sucrose_pulse_sub(sizes,lineseg)))
                line('XData',time_base_pulse(lineseg:lineseg+1), 'YData',data_around_sucrose_pulse(sizes,lineseg:lineseg+1), 'Color',clrzz(data_around_sucrose_pulse_sub(sizes,lineseg)./100+1,:));
            end
        end
    end
end

if(~isempty(data_around_sucrose_pulse))
    hold on;
    plot(time_base_pulse,average2,'k');% pulse sucrose
    data_around_sucrose_pulse(~any(~isnan(data_around_sucrose_pulse), 2),:)=[];
    [r2] = sum(~isnan(data_around_sucrose_pulse(:,2401)));
end

ylabel([num2str(r2) ' pulses occuring on sucrose']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average signal around pulses');

if(~isempty(data_around_pulse))
    hold on;
    
    [r,c] = size(data_around_pulse);
    for sizes = 1:1:r
        for lineseg = 1 : (length(data_around_pulse)-1)
            line('XData',time_base_pulse(lineseg:lineseg+1), 'YData',data_around_pulse(sizes,lineseg:lineseg+1), 'Color',[.75 .75 .75]);
        end
    end
    
    plot(time_base_pulse,average3,'Color','k'); % pulse
    %data_around_pulse(any(isnan(data_around_pulse), 2),:)=[];
    %[r3,c] = size(data_around_pulse);
end

ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

% plotting average of various body length around pulse
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average vel mm/sec (1s smoth)');
average = nanmean(data_around_pulse_vel,1);

if(~isempty(data_around_pulse_vel))
    hold on;
    
    [r,c] = size(data_around_pulse_vel);
    for sizes = 1:1:r
        plot(time_base_pulse, data_around_pulse_vel(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_pulse,average,'Color','k'); % pulse
end
ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

% plotting average of various prob length around pulse
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average proboscis length around pulses');
average = nanmean(data_around_pulse_prob,1);

if(~isempty(data_around_pulse_prob))
    hold on;
    
    [r,c] = size(data_around_pulse_prob);
    for sizes = 1:1:r
        plot(time_base_pulse, data_around_pulse_prob(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_pulse,average,'Color','k'); % pulse
end
ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

% plotting average of various prob length around pulse
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average proboscis length around pulses, X only');
average = nanmean(data_around_pulse_prob_x,1);

if(~isempty(data_around_pulse_prob_x))
    hold on;
    
    [r,c] = size(data_around_pulse_prob_x);
    for sizes = 1:1:r
        plot(time_base_pulse, data_around_pulse_prob_x(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_pulse,average,'Color','k'); % pulse
end
ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

% plotting average of various prob length around pulse
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average proboscis length around pulses, Y only');
average = nanmean(data_around_pulse_prob_y,1);

if(~isempty(data_around_pulse_prob_y))
    hold on;
    
    [r,c] = size(data_around_pulse_prob_y);
    for sizes = 1:1:r
        plot(time_base_pulse, data_around_pulse_prob_y(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_pulse,average,'Color','k'); % pulse
end
ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

% plotting average of various body length around pulse
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body length around pulses');
average = nanmean(data_around_pulse_body,1);

if(~isempty(data_around_pulse_body))
    hold on;
    
    [r,c] = size(data_around_pulse_body);
    for sizes = 1:1:r
        plot(time_base_pulse, data_around_pulse_body(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_pulse,average,'Color','k'); % pulse
end
ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

% plotting average of various body length around pulse
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body length around pulses, X only');
average = nanmean(data_around_pulse_body_x,1);

if(~isempty(data_around_pulse_body_x))
    hold on;
    
    [r,c] = size(data_around_pulse_body_x);
    for sizes = 1:1:r
        plot(time_base_pulse, data_around_pulse_body_x(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_pulse,average,'Color','k'); % pulse
end
ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);



% plotting average of various body length around pulse
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body length around pulses, Y only');
average = nanmean(data_around_pulse_body_y,1);

if(~isempty(data_around_pulse_body_y))
    hold on;
    
    [r,c] = size(data_around_pulse_body_y);
    for sizes = 1:1:r
        plot(time_base_pulse, data_around_pulse_body_y(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_pulse,average,'Color','k'); % pulse
end
ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

% plotting average of various body path metrics around pulse
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body path length around pulses');
average = nanmean(data_around_pulse_path,1);

if(~isempty(data_around_pulse_path))
    hold on;
    
    [r,c] = size(data_around_pulse_path);
    for sizes = 1:1:r
        plot(time_base_pulse, data_around_pulse_path(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_pulse,average,'Color','k'); % pulse
end
ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

% plotting average of various body path metrics around pulse
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body length L3 to tip around pulses');
average = nanmean(data_around_pulse_body_only,1);

if(~isempty(data_around_pulse_body_only))
    hold on;
    
    [r,c] = size(data_around_pulse_body_only);
    for sizes = 1:1:r
        plot(time_base_pulse, data_around_pulse_body_only(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_pulse,average,'Color','k'); % pulse
end
ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

% plotting average of various body path metrics around pulse
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body length L3 to tip around pulses, X only');
average = nanmean(data_around_pulse_body_only_x,1);

if(~isempty(data_around_pulse_body_only_x))
    hold on;
    
    [r,c] = size(data_around_pulse_body_only_x);
    for sizes = 1:1:r
        plot(time_base_pulse, data_around_pulse_body_only_x(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_pulse,average,'Color','k'); % pulse
end
ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

% plotting average of various body path metrics around pulse
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body length L3 to tip around pulses, Y only');
average = nanmean(data_around_pulse_body_only_y,1);

if(~isempty(data_around_pulse_body_only_y))
    hold on;
    
    [r,c] = size(data_around_pulse_body_only_y);
    for sizes = 1:1:r
        plot(time_base_pulse, data_around_pulse_body_only_y(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_pulse,average,'Color','k'); % pulse
end
ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);



% plotting average of various body path metrics around pulse
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body angle around pulses');
average = nanmean(data_around_pulse_angle,1);

if(~isempty(data_around_pulse_angle))
    hold on;
    
    [r,c] = size(data_around_pulse_angle);
    for sizes = 1:1:r
        plot(time_base_pulse, data_around_pulse_angle(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_pulse,average,'Color','k'); % pulse
end
ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

% plotting average of various body path metrics around pulse
figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average body angle L3 to tip around pulses');
average = nanmean(data_around_pulse_angle_only,1);

if(~isempty(data_around_pulse_angle_only))
    hold on;
    
    [r,c] = size(data_around_pulse_angle_only);
    for sizes = 1:1:r
        plot(time_base_pulse, data_around_pulse_angle_only(sizes,:), 'Color',[.75 .75 .75]);
    end
    
    plot(time_base_pulse,average,'Color','k'); % pulse
end
ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

%%%%%%%%%%

%plot signal around pulse

average1 = nanmean(data_around_plain_pulse2,1);
average2 = nanmean(data_around_sucrose_pulse2,1);
average3 = nanmean(data_around_pulse2,1);

r1 = 0;
r2 = 0;
r3 = 0;

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average signal around pulses occuring on plain');
set(gca,'xlim',[-60,60]);

if(~isempty(data_around_plain_pulse2))
    [r,c] = size(data_around_plain_pulse2);
    for sizes = 1:1:r
        for lineseg = 1 : (length(data_around_plain_pulse2)-1)
            if(~isnan(data_around_plain_pulse_sub2(sizes,lineseg)))
                line('XData',time_base_pulse2(lineseg:lineseg+1), 'YData',data_around_plain_pulse2(sizes,lineseg:lineseg+1), 'Color',clrzz(data_around_plain_pulse_sub2(sizes,lineseg)./100+1,:));
            end
        end
    end
end

if(~isempty(data_around_plain_pulse2))
    hold on;
    plot(time_base_pulse2,average1,'k'); % pulse plain
    data_around_plain_pulse2(~any(~isnan(data_around_plain_pulse2), 2),:)=[];
    [r1] = sum(~isnan(data_around_plain_pulse2(:,1201)));
end

ylabel([num2str(r1) ' pulses occuring on plain']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average signal around pulses occuring on sucrose');
set(gca,'xlim',[-60,60]);

if(~isempty(data_around_sucrose_pulse2))
    [r,c] = size(data_around_sucrose_pulse2);
    for sizes = 1:1:r
        for lineseg = 1 : (length(data_around_sucrose_pulse2)-1)
            if(~isnan(data_around_sucrose_pulse_sub2(sizes,lineseg)))
                line('XData',time_base_pulse2(lineseg:lineseg+1), 'YData',data_around_sucrose_pulse2(sizes,lineseg:lineseg+1), 'Color',clrzz(data_around_sucrose_pulse_sub2(sizes,lineseg)./100+1,:));
            end
        end
    end
end

if(~isempty(data_around_sucrose_pulse2))
    hold on;
    plot(time_base_pulse2,average2,'k');% pulse sucrose
    data_around_sucrose_pulse2(~any(~isnan(data_around_sucrose_pulse2), 2),:)=[];
    [r2] = sum(~isnan(data_around_sucrose_pulse2(:,1201)));
end

ylabel([num2str(r2) ' pulses occuring on sucrose']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);

figure; hold on; title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
xlabel('average signal around pulses');

if(~isempty(data_around_pulse2))
    hold on;
    [r,c] = size(data_around_pulse2);
    for sizes = 1:1:r
        for lineseg = 1 : (length(data_around_pulse2)-1)
            line('XData',time_base_pulse2(lineseg:lineseg+1), 'YData',data_around_pulse2(sizes,lineseg:lineseg+1), 'Color',[.5 .5 .5]);
        end
    end
    plot(time_base_pulse2,average3,'Color','k'); % pulse
    data_around_pulse(any(isnan(data_around_pulse), 2),:)=[];
    [r3,c] = size(data_around_pulse2);
end

ylabel([num2str(r1+r2) ' total pulses']);
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
set(gca,'xlim',[-60,60]);
end

