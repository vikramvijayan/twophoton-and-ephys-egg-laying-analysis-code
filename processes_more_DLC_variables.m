
function [recording] = processes_more_DLC_variables(recording)

% recording.abf.CH1_patch_spikes = recording.abf.CH1_patch_spikes';
% recording.abf.CH1_patch_spikes_conv = recording.abf.CH1_patch_spikes_conv';

if(~isfield(recording.abf,'PWMlaser_uWpermm2'))
    recording.abf.PWMlaser_uWpermm2 = zeros(1,length(recording.abf.Time_s));
end

% filter out 300 sec before and after an egg
recording.movie1.no_egg = zeros(1,length(recording.movie1.filtered_wheel));
recording.movie1.egglaid = zeros(1,length(recording.movie1.filtered_wheel));

recording.movie1.no_egg(recording.movie1.eggs)=1;
recording.movie1.egglaid(recording.movie1.eggs)=1;

% no_egg3 is ones where there is egg, 0 otherwise (with no smoothing)
recording.movie1.no_egg3 =recording.movie1.no_egg;

% 600 seconds total, at 25fps
recording.movie1.no_egg = smooth(recording.movie1.no_egg, 600*25);
recording.movie1.no_egg(recording.movie1.no_egg > 0) = NaN;
recording.movie1.no_egg(recording.movie1.no_egg == 0) = 1;

recording.movie1.no_egg = recording.movie1.no_egg';

% no_egg is NaN whree these is egg, 1 otherwise
% no_egg2 is ones where there is egg, 0 otherwise

recording.movie1.no_egg2 = recording.movie1.no_egg;
recording.movie1.no_egg2(recording.movie1.no_egg2 == 1) = 0;
recording.movie1.no_egg2(isnan(recording.movie1.no_egg2)) = 1;


% not calculating the below for non chimson rcordings
% process the PWL laser to find times when there are no pulses
[a b] = find(recording.abf.PWMlaser > 0.1);
recording.abf.no_laser = zeros(length(recording.abf.PWMlaser),1);
recording.abf.no_laser2 = zeros(length(recording.abf.PWMlaser),1);

recording.abf.no_laser(a) = 1;

%10000 is 1 second
%recording.abf.no_laser = causal_filter(200000, recording.abf.no_laser );

% 100 seconds after pulse ends is removed
% this looks complicated because I need to causual filter a subsampled
% signal otherwise it takes too long
recording.abf.no_laser = causal_filter(10000, recording.abf.no_laser(1:100:end));
recording.abf.no_laser = interp1(recording.abf.Time_s(1:100:end), recording.abf.no_laser, recording.abf.Time_s,'previous');

% no laser 2 is ones when there is laser, 0 otherwise
recording.abf.no_laser2(recording.abf.no_laser > 0) = 1;

% no laser is Nan when there is laser 1 otherwise
recording.abf.no_laser(recording.abf.no_laser > 0) = NaN;
recording.abf.no_laser(recording.abf.no_laser == 0) = 1;


% find eggs without pulses close by (60 seconds before)
 recording.movie1.eggs_with_pulse = [];
 eggcntr=0;
laser_movie = interp1(recording.abf.Time_s,recording.abf.no_laser2,recording.movie1.time_stamps);
laser_movie = fliplr(causal_filter(25*60, fliplr(laser_movie)));
for i = 1:1:length(recording.movie1.eggs)
    if(laser_movie(recording.movie1.eggs(i)) > 0)
        recording.movie1.eggs_pulse(i) = 1;
        eggcntr = eggcntr+1;
        recording.movie1.eggs_with_pulse(eggcntr) = recording.movie1.eggs(i);
    else
        recording.movie1.eggs_pulse(i) = 0;
    end
end

if(isempty(recording.movie1.eggs))
    recording.movie1.eggs_pulse = [];
end




% these are variables that were not processed earlier
% like length of abdomen in ony X coordinates etc
recording.movie1.abd_x_L3tip = (recording.movie1.DLC(10,:)-recording.movie1.DLC(1,:));
recording.movie1.abd_x_neck_tip =  (recording.movie1.DLC(13,:)-recording.movie1.DLC(1,:));
recording.movie1.prob_x =  -1.*(recording.movie1.DLC(13,:)-recording.movie1.DLC(16,:));

recording.movie1.abd_y_L3tip = (recording.movie1.DLC(2,:)-recording.movie1.DLC(11,:));
recording.movie1.abd_y_neck_tip =  (recording.movie1.DLC(2,:)-recording.movie1.DLC(14,:));
recording.movie1.prob_y =  (recording.movie1.DLC(17,:)-recording.movie1.DLC(14,:));


for iterate_DLC = 2:1:length(recording.movie1.abd_x_L3tip)
    if((recording.movie1.DLC(3,iterate_DLC) <.95) | (recording.movie1.DLC(12,iterate_DLC) <.95))
        recording.movie1.abd_x_L3tip(iterate_DLC) = NaN;
        recording.movie1.abd_y_L3tip(iterate_DLC) = NaN;
    end
end

for iterate_DLC = 2:1:length(recording.movie1.prob_y)
    if((recording.movie1.DLC(15,iterate_DLC) <.95) | (recording.movie1.DLC(18,iterate_DLC) <.95))
        recording.movie1.prob_x(iterate_DLC) = NaN;
        recording.movie1.prob_y(iterate_DLC) = NaN;
    end
end


for iterate_DLC = 2:1:length(recording.movie1.abd_x_neck_tip)
    if((recording.movie1.DLC(3,iterate_DLC) <.95) | (recording.movie1.DLC(15,iterate_DLC) <.95))
        recording.movie1.abd_x_neck_tip(iterate_DLC) = NaN;
        recording.movie1.abd_y_neck_tip(iterate_DLC) = NaN;
    end
end

% cant really use this because of the risk of negative values for the
% other 3 values calculated here
mean_b = nanmedian(recording.movie1.abd_x_neck_tip);
tmp = recording.movie1.abd_x_neck_tip;
for iterate_DLC = 2:1:length(recording.movie1.abd_x_neck_tip)
    if(abs(tmp(iterate_DLC-1) -  tmp(iterate_DLC)) > mean_b/4)
        %recording.movie1.abd_path_length(iterate_DLC) = recording.movie1.abd_path_length(iterate_DLC-1);
        recording.movie1.abd_x_neck_tip(iterate_DLC) = NaN;
        tmp(iterate_DLC) = tmp(iterate_DLC-1);
    end
end

recording.movie1.prob_length = sqrt( (recording.movie1.prob_x).^2 +(recording.movie1.prob_y).^2);

mean_b = nanmedian(recording.movie1.prob_length);
tmp = recording.movie1.prob_length;
for iterate_DLC = 2:1:length(recording.movie1.abd_length)
    if(abs(tmp(iterate_DLC-1) -  tmp(iterate_DLC)) > mean_b/4)
        % recording.movie1.abd_only_length(iterate_DLC) = recording.movie1.abd_only_length(iterate_DLC-1);
        recording.movie1.prob_length(iterate_DLC) = NaN;
        tmp(iterate_DLC) = tmp(iterate_DLC-1);
    end
end
