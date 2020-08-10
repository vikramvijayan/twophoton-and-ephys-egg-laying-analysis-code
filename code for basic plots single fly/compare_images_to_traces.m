
function [out] = compare_images_to_traces(recording)

warning('off');
tsstacksconcat = [];
pulse_concat = [];

% find all pulse times this is a bit complicted since the pulses have pulses (the stimulation is at a specific Hz)
% this code will find the single beginning and end of a full pulse train (5 pulses of 1
% sec that are spaced 10 sec apart)
temp_diffs = diff(recording.abf.PWMlaser);
[a b] = find(temp_diffs > 3);
a=a+1;
temp_diffs2 = diff(a);
[a1 b1] = find(temp_diffs2 > 1000000);
recording.abf.pulse_ON_times = [a(1); a(a1+1)];

% pulse off is the last time it was ON
temp_diffs = diff(recording.abf.PWMlaser);
[a b] = find(temp_diffs < -3);
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

pulse = zeros(1,length(recording.movie1.time_stamps));
for j =1:1:length(recording.movie1.pulse_ON_times_movie)
    for i =1:1:length(pulse)
        if ( (i >= recording.movie1.pulse_ON_times_movie(j)) & (i <= recording.movie1.pulse_OFF_times_movie(j)))
            pulse(i) = 1;
        end
    end
end

pulse = interp1(recording.movie1.time_stamps, pulse, recording.abf.Time_s,'previous');

for tsall = 1:1:length(recording.tseries)
    
    % tsstacksconcat = [recording.tseries(tsall).tsStack(:,:,1,:,:)];
    mntsstacksconcat = [mean(double(recording.tseries(tsall).tsStack(:,:,1,:,:)),5)];
    pulse_concat = [pulse_concat; pulse(recording.tseries(tsall).middle_index)];
    
    [a b] = find(pulse_concat == 1);
    [a1 b1] = find(pulse_concat == 0);
    
    duringpulse = mean(mntsstacksconcat(:,:,1,a),4);
    duringnopulse = mean(mntsstacksconcat(:,:,1,a1),4);
    figure; hold on;
    caxis('auto')
    %subplot(1,4,1);hold on;
    title(['no pulse ' num2str(tsall)]);
    imagesc(flipud(duringnopulse)); colorbar; colormap('jet'); axis tight; axis equal; cl = caxis;
    
        figure; hold on;
    %subplot(1,4,2);hold on;
    title(['pulse train ' num2str(tsall)]);
    imagesc(flipud(duringpulse)); colorbar; axis tight; axis equal; colormap('jet'); caxis(cl);
    
        figure; hold on;
    %subplot(1,4,3);hold on;
    title(['log pulse train / no pulse' num2str(tsall)]);
    imagesc(flipud(log2(duringpulse./duringnopulse))); colorbar; colormap('jet'); axis tight; axis equal; 
     
    figure; hold on;
    %subplot(1,4,4); hold on;
    title(['pulse train - no pulse' num2str(tsall)]);
    imagesc(flipud(duringpulse-duringnopulse)); colorbar; colormap('jet'); axis tight; axis equal; 
    
    
    tsstacksconcat = [];
    pulse_concat = [];
    
    
end

warning('on');

end



