% parameters used[recording]=add_spike_to_recording(recording,500,-30);
% for 2017 11 17 0002N_truncated

% [trace_out]=add_proboscis_ext_event(head_to_prob,25,30);
function [trace_out]=add_proboscis_ext_event(trace,minpeakdist,minpeakh)



[locs pk] = peakseek(trace,minpeakdist,minpeakh);

trace_out = zeros(1,length(trace));
trace_out(locs)=1;





end
