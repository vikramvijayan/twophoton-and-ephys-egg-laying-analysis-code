
function [out] = correlate_images_to_traces(recording)

warning('off');
%[recording] = processes_more_DLC_variables(recording);

tsstacksconcat = [];
suc_concat = [];
suc_threshold_concat = [];

wheelvel_concat = [];
body_concat = [];
body_angle_concat = [];
body_concat_x = [];
prob_length_concat = [];

trans_concat = [];
transsig_concat = [];
p_concat = [];
s_concat = [];
signal_concat = [];
pulse_concat = [];

tmp_ar_plainegg = zeros(1,length(recording.movie1.sucrose));
tmp_ar_sucroseegg = zeros(1,length(recording.movie1.sucrose));

[a b] = find(recording.movie1.sucrose(recording.movie1.eggs) == 0);
tmp_ar_plainegg(recording.movie1.eggs(b)) = 1;

[a b] = find(recording.movie1.sucrose(recording.movie1.eggs) > 0);
tmp_ar_sucroseegg(recording.movie1.eggs(b)) = 1;

tmp_ar_plainegg = smooth(tmp_ar_plainegg,10*25);
tmp_ar_sucroseegg = smooth(tmp_ar_sucroseegg,10*25);

wheelvelocity = -1.*smooth([0; diff(unwrap(recording.movie1.filtered_wheel))],25);
bodynorm =  recording.movie1.abd_length;
bodyangle =  recording.movie1.abd_angle;
sucrose =      recording.movie1.sucrose;
bodynorm_x =    recording.movie1.abd_x_neck_tip;
prob_length = recording.movie1.prob_length;

wheelvelocity = interp1(recording.movie1.time_stamps, wheelvelocity, recording.abf.Time_s,'previous');
bodynorm = interp1(recording.movie1.time_stamps, bodynorm, recording.abf.Time_s,'previous');
bodyangle = interp1(recording.movie1.time_stamps, bodyangle, recording.abf.Time_s,'previous');
bodynorm_x = interp1(recording.movie1.time_stamps, bodynorm_x, recording.abf.Time_s,'previous');
prob_length = interp1(recording.movie1.time_stamps, prob_length, recording.abf.Time_s,'previous');

array_eggs = interp1(recording.movie1.time_stamps, tmp_ar_sucroseegg, recording.abf.Time_s,'previous');
array_eggp = interp1(recording.movie1.time_stamps, tmp_ar_plainegg, recording.abf.Time_s,'previous');


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

transitionssigned = interp1(recording.movie1.time_stamps, [0; transitionssigned], recording.abf.Time_s,'previous');
transitionsabsval = interp1(recording.movie1.time_stamps, [0; transitionsabsval], recording.abf.Time_s,'previous');


% % % find all pulse times this is a bit complicted since the pulses have pulses (the stimulation is at a specific Hz)
% % % this code will find the single beginning and end of a full pulse train (5 pulses of 1
% % % sec that are spaced 10 sec apart)
% % temp_diffs = diff(recording.abf.PWMlaser);
% % [a b] = find(temp_diffs > 3);
% % a=a+1;
% % temp_diffs2 = diff(a);
% % [a1 b1] = find(temp_diffs2 > 1000000);
% % recording.abf.pulse_ON_times = [a(1); a(a1+1)];
% % 
% % % pulse off is the last time it was ON
% % temp_diffs = diff(recording.abf.PWMlaser);
% % [a b] = find(temp_diffs < -3);
% % temp_diffs2 = diff(a);
% % [a1 b1] = find(temp_diffs2 > 1000000);
% % recording.abf.pulse_OFF_times = [a(a1); a(end)];
% % 
% % recording.movie1.pulse_ON_times_movie = [];
% % %convert abf times to movie index, either the first frame in movie after
% % %pulse is ON, or last frame before it is OFF
% % pulse_ON_times = recording.abf.Time_s(recording.abf.pulse_ON_times);
% % for iterate_pulse =1:1:length(pulse_ON_times)
% %     [a b] = find(recording.movie1.time_stamps-pulse_ON_times(iterate_pulse) >0,1,'first');
% %     recording.movie1.pulse_ON_times_movie(iterate_pulse) = a;
% % end
% % 
% % recording.movie1.pulse_OFF_times_movie = [];
% % pulse_OFF_times = recording.abf.Time_s(recording.abf.pulse_OFF_times);
% % for iterate_pulse =1:1:length(pulse_OFF_times)
% %     [a b] = find(recording.movie1.time_stamps-pulse_OFF_times(iterate_pulse) <0,1,'last');
% %     recording.movie1.pulse_OFF_times_movie(iterate_pulse) = a;
% %     
% % end
% % 
% % pulse = zeros(1,length(recording.movie1.time_stamps));
% % for j =1:1:length(recording.movie1.pulse_ON_times_movie)
% %     for i =1:1:length(pulse)
% %         if ( (i >= recording.movie1.pulse_ON_times_movie(j)) & (i <= recording.movie1.pulse_OFF_times_movie(j)))
% %             pulse(i) = 1;
% %         end
% %     end
% % end

% % pulse = interp1(recording.movie1.time_stamps, pulse, recording.abf.Time_s,'previous');

for tsall = 1:1:length(recording.tseries)
    
    % tsstacksconcat = [recording.tseries(tsall).tsStack(:,:,1,:,:)];
    mntsstacksconcat = [mean(double(recording.tseries(tsall).tsStack(:,:,1,:,:)),5)];
    
    suc_concat = [suc_concat; recording.abf.sucrose(recording.tseries(tsall).middle_index)];
    wheelvel_concat = [wheelvel_concat; wheelvelocity(recording.tseries(tsall).middle_index)];
    body_concat = [body_concat; bodynorm(recording.tseries(tsall).middle_index)];
    body_angle_concat = [body_angle_concat; bodyangle(recording.tseries(tsall).middle_index)];
%    pulse_concat = [pulse_concat; pulse(recording.tseries(tsall).middle_index)];
    body_concat_x = [body_concat_x; bodynorm_x(recording.tseries(tsall).middle_index)];
    prob_length_concat = [prob_length_concat; prob_length(recording.tseries(tsall).middle_index)];

    transsig_concat = [transsig_concat; transitionssigned(recording.tseries(tsall).middle_index)];
    trans_concat = [trans_concat; transitionsabsval(recording.tseries(tsall).middle_index)];
    s_concat = [s_concat; array_eggs(recording.tseries(tsall).middle_index)];
    p_concat = [p_concat; array_eggp(recording.tseries(tsall).middle_index)];
    
    suc_threshold_concat= suc_concat;
    suc_threshold_concat(suc_threshold_concat>0) = 1;
    
    [corr_images1b] = corrimage_tobehavior(mntsstacksconcat, double(suc_concat), 'correlates average z planes to sucrose trace');
    
    %[corr_images1b_2] = corrimage_tobehavior(mntsstacksconcat, double(suc_threshold_concat), 'correlates average z planes to sucrose threshold trace');
    [corr_images2b] = corrimage_tobehavior(mntsstacksconcat, double(wheelvel_concat), 'correlates average z planes to wheel vel trace');
        [corr_images2c] = corrimage_tobehavior(mntsstacksconcat, double(prob_length_concat), 'correlates average z planes to proboscis length');

    [corr_images3b] = corrimage_tobehavior(mntsstacksconcat, double(body_concat), 'correlates average z planes to body length');
    [corr_images3b_2] = corrimage_tobehavior(mntsstacksconcat, double(body_concat_x), 'correlates average z planes to body length x');
    [corr_images3b_3] = corrimage_tobehavior(mntsstacksconcat, double(body_angle_concat), 'correlates average z planes to body angle');
    
    [corr_images4b] = corrimage_tobehavior(mntsstacksconcat, double(transsig_concat), 'correlates average z planes to transitions signed');
    [corr_images5b] = corrimage_tobehavior(mntsstacksconcat, double(trans_concat), 'correlates average z planes to transitions abs');
    
    [corr_images6b] = corrimage_tobehavior(mntsstacksconcat, double(s_concat+p_concat), 'correlates average z planes to eggs');
    %[corr_images6b_2] = corrimage_tobehavior(mntsstacksconcat, double(p_concat), 'correlates average z planes to plain eggs');
    %[corr_images6b_3] = corrimage_tobehavior(mntsstacksconcat, double(s_concat), 'correlates average z planes to sucrose eggs');
    %[corr_images7] = corrimage_tobehavior(mntsstacksconcat, double(pulse_concat), 'correlates average z planes to pulses');
    
    for roi_iterate = 1:1:recording.tseries(tsall).ROI_number
        signal_concat = recording.tseries(tsall).df_over_f(roi_iterate,:);
        [corr_images0b] = corrimage_tobehavior(mntsstacksconcat, double(signal_concat)', ['correlates average z planes to signal in ROI ' num2str(roi_iterate)]);
    end
    
    
    tsstacksconcat = [];
    suc_concat = [];
    suc_threshold_concat = [];
    wheelvel_concat = [];
    body_concat = [];
    body_concat_x = [];
    body_angle_concat = [];
    trans_concat = [];
    transsig_concat = [];
    p_concat = [];
    s_concat = [];
    signal_concat = [];
    pulse_concat = [];
    prob_length_concat = [];
    
end

warning('on');

end



