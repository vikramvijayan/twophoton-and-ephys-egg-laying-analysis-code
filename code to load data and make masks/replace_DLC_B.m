function recording = replace_DLC_B(recording ,DLC_csv_B)

substrate_vector = [0,0,0];

%currently the alst 2 csvs are hard coded
warning('off');

% this is new code that is aliging to the end of the movie. can't use the
% first pause anymore because of the fact we dont have triggers

%modifiedStr = DLC_csv_B;
        modifiedStr = strrep([recording.movie2.filename(1:1:end-6) 'DeepCut_resnet50_GPU_05062019BMay6shuffle1_1030000.csv'], '.avi', '');

        
[bodypart]= importDeepLabCutfile_VV(modifiedStr, 4, inf);

if(length(recording.movie2.DLC) == length(bodypart(2:13,:)))
    
          [bodypart]= importDeepLabCutfile_VV(modifiedStr, 4, inf);
        recording.movie2.DLC = bodypart(2:13,:);

        
        
        
        
filter_wheel = filter_wheel_trace( recording.abf.Wheel(abs(round(recording.movie2.time_stamps.*10000))-1));
filter_wheel = unwrap(filter_wheel) - filter_wheel(recording.movie2.bound);
filter_wheel = mod(filter_wheel,2*pi);
sucrose = zeros(1,length(filter_wheel));

recording.movie2.circle_x = recording.movie2.DLC(1,:);
recording.movie2.circle_y = recording.movie2.DLC(2,:);

%circle_parameters = CircleFitByPratt([recording.movie2.circle_x(2000:end)',recording.movie2.circle_y(2000:end)']);
circle_parameters = CircleFitByPratt([recording.movie2.circle_x',recording.movie2.circle_y']);

filter_wheel_DLC = atan2(recording.movie2.circle_y-circle_parameters(2),recording.movie2.circle_x-circle_parameters(1));

filter_wheel_DLC = unwrap(filter_wheel_DLC) - filter_wheel_DLC(recording.movie2.bound);
filter_wheel_DLC = mod(filter_wheel_DLC,2*pi);
filter_wheel_DLC = 2*pi-filter_wheel_DLC;
filter_wheel_DLC = filter_wheel_DLC';

%filter_wheel_DLC = filter_wheel_trace_DLC(filter_wheel_DLC);



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
    
        
    if(substrate_vector == [200,200,500])
        [a b] = find(filter_wheel_DLC >  0);
        sucrose(a) = 500;
        [a b] = find(filter_wheel_DLC > 2*pi/3);
        sucrose(a) = 200;
        [a b] = find(filter_wheel_DLC > 4*pi/3);
        sucrose(a) = 200;
    end
    
    if(substrate_vector == [0,500,500])
        [a b] = find(filter_wheel_DLC >  0);
        sucrose(a) = 500;
        [a b] = find(filter_wheel_DLC > 2*pi/3);
        sucrose(a) = 500;
        [a b] = find(filter_wheel_DLC > 4*pi/3);
        sucrose(a) = 0;
    end
    
        
    if(substrate_vector == [0,200,200])
        [a b] = find(filter_wheel_DLC >  0);
        sucrose(a) = 200;
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
    
        if(substrate_vector == [200,200,500])
        [a b] = find(filter_wheel_DLC >  0);
        sucrose(a) = 500;
        [a b] = find(filter_wheel_DLC > 2*pi/3);
        sucrose(a) = 200;
        [a b] = find(filter_wheel_DLC > 4*pi/3);
        sucrose(a) = 200;
        end
    
    if(substrate_vector == [0,500,0])
        [a b] = find(filter_wheel_DLC >  0);
        sucrose(a) = 0;
        [a b] = find(filter_wheel_DLC > 2*pi/3);
        sucrose(a) = 500;
        [a b] = find(filter_wheel_DLC > 4*pi/3);
        sucrose(a) = 0;
    end
    
    
    if(substrate_vector == [0,0,0])
        [a b] = find(filter_wheel_DLC >  0);
        sucrose(a) = 0;
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
%recording.movie1.filtered_wheel = recording.movie1.filtered_wheel';
recording.movie1.sucrose = interp1(recording.movie2.time_stamps, recording.movie2.sucrose, recording.movie1.time_stamps,'previous');
recording.movie1.sucrose = recording.movie1.sucrose';
        
    
else
    disp 'mismatch in lengths, nothing was done'
    
end
end
