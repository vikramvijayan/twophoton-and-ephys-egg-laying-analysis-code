function [recording] = convert_PWMlaser_to_power(recording)

% calibration table for chrimson 655nm in uW/mm2
% 11/13/2019 values
PWM_level =   [1,   .5,  .25, .2,  .1,  .05, .025, .01];
uWmm2_level = [650, 380, 250, 230, 157, 135, 130, 115];

% % calivration for 595 nm on 11/27 [this is a bit off from what I wrote in
% % the file but it is fine]
%  PWM_level =   [1,   .5, .01];
%  uWmm2_level = [640, 370, 115];
% 
% % calivration for 617 nm on 11/27
% PWM_level =   [1,   .5, .01];
% uWmm2_level = [614, 400, 115];
% 
% % calivration for 565 nm on 11/27
% PWM_level =   [1,   .5, .01];
% uWmm2_level = [640, 370, 115];

% find all pulse times this is a bit complicted since the pulses have pulses (the stimulation is at a specific Hz or PWM)
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
recording.abf.pulse_ON_times = tmp_times;



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
recording.abf.pulse_OFF_times = tmp_times;


smoothed_PWMlaser = smooth(recording.abf.PWMlaser,1000);

recording.abf.PWMlaser_uWpermm2 = zeros(length(recording.abf.PWMlaser),1);
for i = 1:1:length(recording.abf.pulse_ON_times)
    val_of_laser = smoothed_PWMlaser(ceil((recording.abf.pulse_ON_times(i)+recording.abf.pulse_OFF_times(i))./2));
    vq = interp1(PWM_level,uWmm2_level,val_of_laser/5,'linear');
    recording.abf.PWMlaser_uWpermm2((recording.abf.pulse_ON_times(i)):(recording.abf.pulse_OFF_times(i))) = vq;
end

figure; hold on; plot(recording.abf.PWMlaser_uWpermm2)
hold on; plot(recording.abf.PWMlaser)

end

