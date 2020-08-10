
function [out] = correlate_images_to_traces_eachz(recording)

warning('off');
tsstacksconcat = [];
suc_concat = [];
suc_threshold_concat = [];

wheelvel_concat = [];
body_concat = [];
body_angle_concat = [];

trans_concat = [];
transsig_concat = [];
p_concat = [];
s_concat = [];
signal_concat = [];

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

wheelvelocity = interp1(recording.movie1.time_stamps, wheelvelocity, recording.abf.Time_s,'previous');
bodynorm = interp1(recording.movie1.time_stamps, bodynorm, recording.abf.Time_s,'previous');
bodyangle = interp1(recording.movie1.time_stamps, bodyangle, recording.abf.Time_s,'previous');
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

for tsall = 1:1:length(recording.tseries)
    
    tsstacksconcat = [recording.tseries(tsall).tsStack(:,:,1,:,:)];
    %mntsstacksconcat = [mean(double(recording.tseries(tsall).tsStack(:,:,1,:,:)),5)];
    
    suc_concat = [suc_concat; recording.abf.sucrose(recording.tseries(tsall).middle_index)];
    wheelvel_concat = [wheelvel_concat; wheelvelocity(recording.tseries(tsall).middle_index)];
    body_concat = [body_concat; bodynorm(recording.tseries(tsall).middle_index)];
    body_angle_concat = [body_angle_concat; bodyangle(recording.tseries(tsall).middle_index)];
    
    transsig_concat = [transsig_concat; transitionssigned(recording.tseries(tsall).middle_index)];
    trans_concat = [trans_concat; transitionsabsval(recording.tseries(tsall).middle_index)];
    s_concat = [s_concat; array_eggs(recording.tseries(tsall).middle_index)];
    p_concat = [p_concat; array_eggp(recording.tseries(tsall).middle_index)];
    %signal_concat = [signal_concat; recording.tseries(tsall).df_over_f(roiss,:)];
    
    suc_threshold_concat= suc_concat;
    suc_threshold_concat(suc_threshold_concat>0) = 1;
    
   % [corr_images0b] = corrimage_tobehavior(double(tsstacksconcat), double(signal_concat)', 'correlates each z planes to signal in designated ROI');
    [corr_images1b] = corrimage_tobehavior(double(tsstacksconcat), double(suc_concat), 'correlates each z planes to sucrose trace');
    
    [corr_images1b_2] = corrimage_tobehavior(double(tsstacksconcat), double(suc_threshold_concat), 'correlates each z planes to sucrose threshold trace');
    [corr_images2b] = corrimage_tobehavior(double(tsstacksconcat), double(wheelvel_concat), 'correlates each z planes to wheel vel trace');
    
    [corr_images3b] = corrimage_tobehavior(double(tsstacksconcat), double(body_concat), 'correlates each z planes to body length');
    [corr_images3b_2] = corrimage_tobehavior(double(tsstacksconcat), double(body_angle_concat), 'correlates each z planes to body angle');
    
    [corr_images4b] = corrimage_tobehavior(double(tsstacksconcat), double(transsig_concat), 'correlates each z planes to transitions signed');
    [corr_images5b] = corrimage_tobehavior(double(tsstacksconcat), double(trans_concat), 'correlates each z planes to transitions abs');
    
%     [corr_images6b] = corrimage_tobehavior(double(tsstacksconcat), double(s_concat+p_concat), 'correlates each z planes to eggs');
%     [corr_images6b_2] = corrimage_tobehavior(double(tsstacksconcat), double(p_concat), 'correlates each z planes to plain eggs');
%     [corr_images6b_3] = corrimage_tobehavior(double(tsstacksconcat), double(s_concat), 'correlates each z planes to sucrose eggs');
    
    
    
    tsstacksconcat = [];
    suc_concat = [];
    suc_threshold_concat = [];
    wheelvel_concat = [];
    body_concat = [];
    body_angle_concat = [];
    trans_concat = [];
    transsig_concat = [];
    p_concat = [];
    s_concat = [];
    signal_concat = [];
    
    
end

warning('on');

end



