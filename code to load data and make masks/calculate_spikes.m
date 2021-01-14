% similar to add_spike_to_recording, but it dosent need a recording as
% input
function [spikes_out]=calculate_spikes(recording,minpeakdist,minpeakh)


% filter then find spikes
CH1_patch_filtered = highpass(recording,10,10000);

%recording.abf.CH1_patch_filtered = recording.abf.CH1_patch;

[locs pk] = peakseek(CH1_patch_filtered,minpeakdist,minpeakh);

q = zeros(1,length(CH1_patch_filtered));
q(locs)=1;

spikes_out = q;


end

