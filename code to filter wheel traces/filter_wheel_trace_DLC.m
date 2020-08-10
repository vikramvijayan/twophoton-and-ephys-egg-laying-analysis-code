
function trace_out = filter_wheel_trace_DLC(trace)



% now make sure that the difference between values is never larger than a
% theshold
% if it is, then just use (hold) the previous value
for i = 1:1:(length(trace)-1)
    if(abs(diff(unwrap([trace(i), trace(i+1)]))) > .2) %was 1
        trace(i+1) = trace(i);
    end 
end



trace_out = trace;
end
