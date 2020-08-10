function [ filt_out ] = causal_filter(window_in_abf, vector )

% a window of 10000 would be 1 sec
filt_out = filter((1/window_in_abf)*ones(1,window_in_abf),1,vector);
end

