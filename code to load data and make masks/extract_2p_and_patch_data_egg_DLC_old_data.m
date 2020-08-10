function out = extract_2p_and_patch_data_egg_DLC_old_data(data_matfile, movie1)

substrate_vector = [0,0,500];

%load([data_matfile '.mat']);
load([data_matfile]);
% this is optional if the tsStack needs to be redone
%recording = reload_tsStack(recording);
if(~isempty(movie1))
    
    %change the name of the movie
    recording.movie1.filename = movie1;
    
    [a b] = find(diff(recording.cam.Time_s) > .25,1,'first');
    
    %  modifiedStr = strrep([recording.movie1.filename(1:1:end-6) 'DeepCut_resnet50_GPU_05062019May6shuffle1_1030000.csv'], '.avi', '');
    modifiedStr = strrep([recording.movie1.filename 'DeepCut_resnet50_GPU_05062019May6shuffle1_1030000.csv'], '.avi', '');
    
    
    [bodypart]= importDeepLabCutfile_VV(modifiedStr, 4, inf);
    recording.movie1.DLC = bodypart(2:19,:);
end

recording.movie1.abd_length = sqrt( (recording.movie1.DLC(1,:)-recording.movie1.DLC(13,:)).^2 +(recording.movie1.DLC(2,:)-recording.movie1.DLC(14,:)).^2);
recording.movie1.abd_angle = atan( ( (recording.movie1.DLC(2,:)-recording.movie1.DLC(14,:))) ./ ((recording.movie1.DLC(13,:)-recording.movie1.DLC(1,:))));
recording.movie1.abd_likehood = recording.movie1.DLC(3,:)+recording.movie1.DLC(15,:);

% adding three line segments
recording.movie1.abd_path_length = sqrt( (recording.movie1.DLC(1,:)-recording.movie1.DLC(4,:)).^2 +(recording.movie1.DLC(2,:)-recording.movie1.DLC(5,:)).^2)+ sqrt( (recording.movie1.DLC(4,:)-recording.movie1.DLC(7,:)).^2 +(recording.movie1.DLC(5,:)-recording.movie1.DLC(8,:)).^2)+sqrt( (recording.movie1.DLC(7,:)-recording.movie1.DLC(10,:)).^2 +(recording.movie1.DLC(8,:)-recording.movie1.DLC(11,:)).^2);
recording.movie1.abd_only_length = sqrt( (recording.movie1.DLC(1,:)-recording.movie1.DLC(10,:)).^2 +(recording.movie1.DLC(2,:)-recording.movie1.DLC(11,:)).^2);
recording.movie1.abd_only_angle = atan( ( (recording.movie1.DLC(2,:)-recording.movie1.DLC(11,:))) ./ ((recording.movie1.DLC(10,:)-recording.movie1.DLC(1,:))));
recording.movie1.abd_only_angle = unwrap(2*recording.movie1.abd_only_angle)./2+pi;

for iterate_DLC = 2:1:length(recording.movie1.abd_likehood)
    if((recording.movie1.DLC(3,iterate_DLC) <.95) | (recording.movie1.DLC(15,iterate_DLC) <.95))
        recording.movie1.abd_length(iterate_DLC) = recording.movie1.abd_length(iterate_DLC-1);
        recording.movie1.abd_angle(iterate_DLC) = recording.movie1.abd_angle(iterate_DLC-1);
        recording.movie1.abd_angle(iterate_DLC) = NaN;
        recording.movie1.abd_length(iterate_DLC) = NaN;
        
    end
end

for iterate_DLC = 2:1:length(recording.movie1.abd_likehood)
    if( (recording.movie1.DLC(3,iterate_DLC) <.95) | (recording.movie1.DLC(6,iterate_DLC) <.95) | (recording.movie1.DLC(9,iterate_DLC) <.95) | (recording.movie1.DLC(12,iterate_DLC) <.95))
        recording.movie1.abd_path_length(iterate_DLC) = recording.movie1.abd_path_length(iterate_DLC-1);
        recording.movie1.abd_path_length(iterate_DLC) = NaN;
    end
end

for iterate_DLC = 2:1:length(recording.movie1.abd_likehood)
    if((recording.movie1.DLC(3,iterate_DLC) <.95) | (recording.movie1.DLC(12,iterate_DLC) <.95))
        recording.movie1.abd_only_length(iterate_DLC) = recording.movie1.abd_only_length(iterate_DLC-1);
        recording.movie1.abd_only_angle(iterate_DLC) = recording.movie1.abd_only_angle(iterate_DLC-1);
        recording.movie1.abd_only_length(iterate_DLC) = NaN;
        recording.movie1.abd_only_angle(iterate_DLC) = NaN;
        
    end
end

tmp = recording.movie1.abd_angle;
for iterate_DLC = 2:1:length(recording.movie1.abd_likehood)
    if(abs(circ_dist(tmp(iterate_DLC-1) ,  tmp(iterate_DLC))) > pi/6)
        %recording.movie1.abd_angle(iterate_DLC) = recording.movie1.abd_angle(iterate_DLC-1);
        recording.movie1.abd_angle(iterate_DLC) = NaN;
        tmp(iterate_DLC) = tmp(iterate_DLC-1);
        
    end
end

tmp = recording.movie1.abd_only_angle;
for iterate_DLC = 2:1:length(recording.movie1.abd_likehood)
    if(abs(circ_dist(tmp(iterate_DLC-1) ,  tmp(iterate_DLC))) > pi/6)
        %recording.movie1.abd_only_angle(iterate_DLC) = recording.movie1.abd_only_angle(iterate_DLC-1);
        recording.movie1.abd_only_angle(iterate_DLC) =  NaN;
        tmp(iterate_DLC) = tmp(iterate_DLC-1);
        
    end
end


mean_b = nanmedian(recording.movie1.abd_path_length);
tmp = recording.movie1.abd_path_length;
for iterate_DLC = 2:1:length(recording.movie1.abd_likehood)
    if(abs(tmp(iterate_DLC-1) -  tmp(iterate_DLC)) > mean_b/4)
        %recording.movie1.abd_path_length(iterate_DLC) = recording.movie1.abd_path_length(iterate_DLC-1);
        recording.movie1.abd_path_length(iterate_DLC) = NaN;
        tmp(iterate_DLC) = tmp(iterate_DLC-1);
    end
end

mean_b = nanmedian(recording.movie1.abd_length);
tmp = recording.movie1.abd_length;
for iterate_DLC = 2:1:length(recording.movie1.abd_likehood)
    if(abs(tmp(iterate_DLC-1) -  tmp(iterate_DLC)) > mean_b/4)
        %recording.movie1.abd_length(iterate_DLC) = recording.movie1.abd_length(iterate_DLC-1);
        recording.movie1.abd_length(iterate_DLC) = NaN;
        tmp(iterate_DLC) = tmp(iterate_DLC-1);
    end
end

mean_b = nanmedian(recording.movie1.abd_only_length);
tmp = recording.movie1.abd_only_length;
for iterate_DLC = 2:1:length(recording.movie1.abd_likehood)
    if(abs(tmp(iterate_DLC-1) -  tmp(iterate_DLC)) > mean_b/4)
        % recording.movie1.abd_only_length(iterate_DLC) = recording.movie1.abd_only_length(iterate_DLC-1);
        recording.movie1.abd_only_length(iterate_DLC) = NaN;
        tmp(iterate_DLC) = tmp(iterate_DLC-1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filter_wheel = filter_wheel_trace( recording.abf.Wheel(abs(round(recording.movie2.time_stamps.*10000))-1));
filter_wheel = unwrap(filter_wheel) - filter_wheel(recording.movie2.bound);
filter_wheel = mod(filter_wheel,2*pi);
sucrose = zeros(1,length(filter_wheel));

recording.movie2.circle_x = recording.movie2.DLC(1,:);
recording.movie2.circle_y = recording.movie2.DLC(2,:);

circle_parameters = CircleFitByPratt([recording.movie2.circle_x',recording.movie2.circle_y']);

filter_wheel_DLC = atan2(recording.movie2.circle_y-circle_parameters(2),recording.movie2.circle_x-circle_parameters(1));

filter_wheel_DLC = unwrap(filter_wheel_DLC) - filter_wheel_DLC(recording.movie2.bound);
filter_wheel_DLC = mod(filter_wheel_DLC,2*pi);
filter_wheel_DLC = 2*pi-filter_wheel_DLC;
filter_wheel_DLC = filter_wheel_DLC';





if(length(substrate_vector) == 2)
    if(substrate_vector == [0,500])
        [a b] = find(filter_wheel_DLC > pi);
        sucrose(a) = 500;
        [a b] = find(filter_wheel_DLC < pi);
        sucrose(a) = 0;
    end
    
else
    
    if(substrate_vector == [0,500,200])
        [a b] = find(filter_wheel_DLC >  0);
        sucrose(a) = 200;
        [a b] = find(filter_wheel_DLC > 2*pi/3);
        sucrose(a) = 500;
        [a b] = find(filter_wheel_DLC > 4*pi/3);
        sucrose(a) = 0;
    end
    
    if(substrate_vector == [0,200,500])
        [a b] = find(filter_wheel_DLC >  0);
        sucrose(a) = 500;
        [a b] = find(filter_wheel_DLC > 2*pi/3);
        sucrose(a) = 200;
        [a b] = find(filter_wheel_DLC > 4*pi/3);
        sucrose(a) = 0;
    end
    
    if(substrate_vector == [200,500,0])
        [a b] = find(filter_wheel_DLC >  0);
        sucrose(a) = 0;
        [a b] = find(filter_wheel_DLC > 2*pi/3);
        sucrose(a) = 500;
        [a b] = find(filter_wheel_DLC > 4*pi/3);
        sucrose(a) = 200;
    end
    
      if(substrate_vector == [0,0,500])
        [a b] = find(filter_wheel_DLC >  0);
        sucrose(a) = 500;
        [a b] = find(filter_wheel_DLC > 2*pi/3);
        sucrose(a) = 0;
        [a b] = find(filter_wheel_DLC > 4*pi/3);
        sucrose(a) = 0;
      end
    
end

recording.movie2.sucrose = sucrose';
recording.movie2.filtered_wheel = filter_wheel;
recording.movie2.filtered_wheel_DLC = filter_wheel_DLC;
recording.abf.sucrose = interp1(recording.movie2.time_stamps, recording.movie2.sucrose, recording.abf.Time_s,'previous');

recording.movie1.filtered_wheel = interp1(recording.movie2.time_stamps, recording.movie2.filtered_wheel_DLC, recording.movie1.time_stamps,'previous');
recording.movie1.sucrose = interp1(recording.movie2.time_stamps, recording.movie2.sucrose, recording.movie1.time_stamps,'previous');
recording.movie1.sucrose = recording.movie1.sucrose';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a b] = find(recording.abf.sucrose == 1);
recording.abf.sucrose(a) = 500;

[a b] = find(recording.movie1.sucrose == 1);
recording.movie1.sucrose(b) = 500;

%save([data_matfile '_DLC.mat'],'recording');
%save([data_matfile '.mat'],'recording');
save(['C:\Users\vikra\Documents\Postdoc\working folder\tmp\' data_matfile],'recording');
warning('on');
out = 1;
end