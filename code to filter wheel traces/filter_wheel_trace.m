
function trace_out = filter_wheel_trace(trace)

%figure; scatter(1:length(trace), trace,'m'); hold on;


% fill in all zeros (identified by values less than .5)  with the last legitimate value
[a b] = find(trace < .5);

for i = 1:1:length(a)
    [a1 b1] = find(trace(1:a(i)) > .9,1,'last');
    if(~isempty(a1))
        trace(a(i)) = trace(a1(1));
    else
        trace(a(i)) = 1.5;
    end
end

% now scale the data between 0 and 2*pi
trace = trace-min(trace);
trace = trace./max(trace);
trace = trace.*2*pi;


% now make sure that the difference between values is never larger than a
% theshold
% if it is, then just use (hold) the previous value
for i = 1:1:(length(trace)-1)
    if(abs(diff(unwrap([trace(i), trace(i+1)]))) > 1.5) %was 1
        trace(i+1) = trace(i);
    end 
end

% this smoothing was added 9/4/2018
trace = smooth(unwrap(trace),51);
trace = mod(trace,2*pi);

% % now make sure that the difference between values is never larger than a
% % theshold
% % if it is, then just use (hold) the previous value
% for i = (length(trace)):-1:2
%     if(abs(diff(unwrap([trace(i), trace(i-1)]))) > 1.5) %was 1
%         trace(i-1) = trace(i);
%     end 
% end


%scatter(1:length(trace),trace,'k');

trace_out = trace;
end
