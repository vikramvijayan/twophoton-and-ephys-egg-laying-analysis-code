% parameters used[recording]=add_spike_to_recording(recording,500,-30);
% for 2017 11 17 0002N_truncated

% for 9/17 0003 use 500, .25 with highpess filtering recording.abf.CH1_patch = highpass(recording.abf.CH1_patch,100,10000);
% for 9/17 0004 use 500, .1 with highpess filtering recording.abf.CH1_patch = highpass(recording.abf.CH1_patch,100,10000);
% for 9/17 0006 use 500, .35 with highpess filtering recording.abf.CH1_patch = highpass(recording.abf.CH1_patch,10,10000);
% for 9/18 0007 use 250, 2.25 with no highpass filtering;
% for 9/18 0009 use 250, 1.2 with no highpass filtering;
% for 9/18 0010 use 100, .05 with no highpass filtering;
% for 9/18 0011 use 100, .05 with no highpass filtering;
% for 9/18 0011 redofilter use 100, .08 with no highpass filtering;

% for 9/19 0000 use 100, .2 with no highpess filtering
% for 9/19 0002 use 100, .4 with highpess filtering recording.abf.CH1_patch = highpass(recording.abf.CH1_patch,100,10000);

function [recording_out]=add_spike_to_recording(recording,minpeakdist,minpeakh)


% filter then find spikes
% added for 9/17/2019 0003
%recording.abf.CH1_patch_filtered = highpass(recording.abf.CH1_patch,100,10000);
recording.abf.CH1_patch_filtered = recording.abf.CH1_patch;

[locs pk] = peakseek(recording.abf.CH1_patch_filtered,minpeakdist,minpeakh);

q = zeros(1,length(recording.abf.CH1_patch_filtered));
q(locs)=1;

recording.abf.CH1_patch_spikes = q;

% this convolution does not preserve the area under the curve just the max
window = 10000; %1 sec
w = gausswin(window); 
y = filter(w,1,recording.abf.CH1_patch_spikes);

y = [ y((window/2)+1:end), NaN(1,window/2)];

recording.abf.CH1_patch_spikes_conv = y;




% repeat with mor generous parameters so taht spikes removed is very robust
% commented out at the moment
%[locs pk] = peakseek(recording.abf.CH1_patch,1,-35);




% q = zeros(1,length(recording.abf.CH1_patch));
% q(locs)=1;
% 
% w = gausswin(4*50000); %20 seconds
% y = filter(w,1,q);
% 
% y = [ y(100000:end), NaN(1,100000)];

qq = recording.abf.CH1_patch;
[a b] = find(y >0);
qq(b) = NaN;

recording.abf.CH1_patch_spikes_removed = qq;




% this convolution preserves the area
window = 50000; %5 sec
w = gausswin(window); 
y = filter(w./(sum(w)),1,recording.abf.CH1_patch_spikes);

y = [ y((window/2)+1:end), NaN(1,window/2)];

recording.abf.CH1_patch_spikes_conv_area = y;
recording.abf.CH1_patch_spikes_conv_area = recording.abf.CH1_patch_spikes_conv_area';


window = 50000; %5 sec

y = filter((1/window)*ones(1,window),1,recording.abf.CH1_patch_spikes);
y = [ y((window/2)+1:end), NaN(1,window/2)];

recording.abf.CH1_patch_spikes_conv_area_rect = y;

% this line of code was in the pre 2020 code and probably made the rect
% convolution equal to the gaussian one
%recording.abf.CH1_patch_spikes_conv_area_rect = recording.abf.CH1_patch_spikes_conv_area';
recording.abf.CH1_patch_spikes_conv_area_rect = recording.abf.CH1_patch_spikes_conv_area_rect';





recording.abf.CH1_patch_spikes = recording.abf.CH1_patch_spikes';
recording.abf.CH1_patch_spikes_conv = recording.abf.CH1_patch_spikes_conv';






% recording.abf.CH1_patch= recording.abf.CH1_patch(1:end-100000);
% recording.abf.CH1_patch_spikes_removed= recording.abf.CH1_patch_spikes_removed(1:end-100000);
% recording.abf.CH1_patch_spikes= recording.abf.CH1_patch_spikes(1:end-100000);
% recording.abf.CH1_patch_spikes_conv= recording.abf.CH1_patch_spikes_conv(1:end-100000);
% recording.abf.Time_s = recording.abf.Time_s(1:end-100000);
% recording.abf.sucrose = recording.abf.sucrose(1:end-100000);

recording_out = recording;

end
