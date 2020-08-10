function [time_base_to_return, data_to_average_interp] = average_around_event_nearest(data_to_average, time_base, time_of_events, time_base_to_return)

data_to_average_interp = [];
for i = 1:1:length(time_of_events)
        data_to_average_interp(i,:) = interp1(time_base-time_of_events(i), data_to_average, time_base_to_return,'nearest');   
end
average = nanmean(data_to_average_interp,1);
end