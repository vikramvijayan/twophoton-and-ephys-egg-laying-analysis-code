function recording = replace_DLC_A(recording ,DLC_csv_A)

%currently the alst 2 csvs are hard coded
warning('off');

% this is new code that is aliging to the end of the movie. can't use the
% first pause anymore because of the fact we dont have triggers

modifiedStr = DLC_csv_A;
[bodypart]= importDeepLabCutfile_VV(modifiedStr, 4, inf);

if(length(recording.movie1.DLC) == length(bodypart(2:19,:)))
    
    recording.movie1.DLC = bodypart(2:19,:);
    
    
    recording.movie1.abd_length = sqrt( (recording.movie1.DLC(1,:)-recording.movie1.DLC(13,:)).^2 +(recording.movie1.DLC(2,:)-recording.movie1.DLC(14,:)).^2);
    recording.movie1.abd_angle = atan( ( (recording.movie1.DLC(2,:)-recording.movie1.DLC(14,:))) ./ ((recording.movie1.DLC(13,:)-recording.movie1.DLC(1,:))));
    
    recording.movie1.prob_length = sqrt( (recording.movie1.DLC(16,:)-recording.movie1.DLC(13,:)).^2 +(recording.movie1.DLC(17,:)-recording.movie1.DLC(14,:)).^2);
    
    recording.movie1.abd_path_length = sqrt( (recording.movie1.DLC(1,:)-recording.movie1.DLC(4,:)).^2 +(recording.movie1.DLC(2,:)-recording.movie1.DLC(5,:)).^2)+ sqrt( (recording.movie1.DLC(4,:)-recording.movie1.DLC(7,:)).^2 +(recording.movie1.DLC(5,:)-recording.movie1.DLC(8,:)).^2)+sqrt( (recording.movie1.DLC(7,:)-recording.movie1.DLC(10,:)).^2 +(recording.movie1.DLC(8,:)-recording.movie1.DLC(11,:)).^2);
    recording.movie1.abd_only_length = sqrt( (recording.movie1.DLC(1,:)-recording.movie1.DLC(10,:)).^2 +(recording.movie1.DLC(2,:)-recording.movie1.DLC(11,:)).^2);
    recording.movie1.abd_only_angle = atan( ( (recording.movie1.DLC(2,:)-recording.movie1.DLC(11,:))) ./ ((recording.movie1.DLC(10,:)-recording.movie1.DLC(1,:))));
    recording.movie1.abd_only_angle = unwrap(2*recording.movie1.abd_only_angle)./2+pi;
    
    for iterate_DLC = 2:1:length(recording.movie1.abd_length)
        if((recording.movie1.DLC(15,iterate_DLC) <.95) | (recording.movie1.DLC(18,iterate_DLC) <.95))
            recording.movie1.prob_length(iterate_DLC) = NaN;
        end
    end
    
    for iterate_DLC = 2:1:length(recording.movie1.abd_length)
        if((recording.movie1.DLC(3,iterate_DLC) <.95) | (recording.movie1.DLC(15,iterate_DLC) <.95))
            recording.movie1.abd_length(iterate_DLC) = recording.movie1.abd_length(iterate_DLC-1);
            recording.movie1.abd_angle(iterate_DLC) = recording.movie1.abd_angle(iterate_DLC-1);
            recording.movie1.abd_angle(iterate_DLC) = NaN;
            recording.movie1.abd_length(iterate_DLC) = NaN;
            
        end
    end
    
    for iterate_DLC = 2:1:length(recording.movie1.abd_length)
        if( (recording.movie1.DLC(3,iterate_DLC) <.95) | (recording.movie1.DLC(6,iterate_DLC) <.95) | (recording.movie1.DLC(9,iterate_DLC) <.95) | (recording.movie1.DLC(12,iterate_DLC) <.95))
            recording.movie1.abd_path_length(iterate_DLC) = recording.movie1.abd_path_length(iterate_DLC-1);
            recording.movie1.abd_path_length(iterate_DLC) = NaN;
        end
    end
    
    for iterate_DLC = 2:1:length(recording.movie1.abd_length)
        if((recording.movie1.DLC(3,iterate_DLC) <.95) | (recording.movie1.DLC(12,iterate_DLC) <.95))
            recording.movie1.abd_only_length(iterate_DLC) = recording.movie1.abd_only_length(iterate_DLC-1);
            recording.movie1.abd_only_angle(iterate_DLC) = recording.movie1.abd_only_angle(iterate_DLC-1);
            recording.movie1.abd_only_length(iterate_DLC) = NaN;
            recording.movie1.abd_only_angle(iterate_DLC) = NaN;
            
        end
    end
    
    tmp = recording.movie1.abd_angle;
    for iterate_DLC = 2:1:length(recording.movie1.abd_length)
        if(abs(circ_dist(tmp(iterate_DLC-1) ,  tmp(iterate_DLC))) > pi/6)
            %recording.movie1.abd_angle(iterate_DLC) = recording.movie1.abd_angle(iterate_DLC-1);
            recording.movie1.abd_angle(iterate_DLC) = NaN;
            tmp(iterate_DLC) = tmp(iterate_DLC-1);
            
        end
    end
    
    tmp = recording.movie1.abd_only_angle;
    for iterate_DLC = 2:1:length(recording.movie1.abd_length)
        if(abs(circ_dist(tmp(iterate_DLC-1) ,  tmp(iterate_DLC))) > pi/6)
            %recording.movie1.abd_only_angle(iterate_DLC) = recording.movie1.abd_only_angle(iterate_DLC-1);
            recording.movie1.abd_only_angle(iterate_DLC) =  NaN;
            tmp(iterate_DLC) = tmp(iterate_DLC-1);
            
        end
    end
    
    mean_b = nanmedian(recording.movie1.abd_path_length);
    tmp = recording.movie1.abd_path_length;
    for iterate_DLC = 2:1:length(recording.movie1.abd_length)
        if(abs(tmp(iterate_DLC-1) -  tmp(iterate_DLC)) > mean_b/4)
            %recording.movie1.abd_path_length(iterate_DLC) = recording.movie1.abd_path_length(iterate_DLC-1);
            recording.movie1.abd_path_length(iterate_DLC) = NaN;
            tmp(iterate_DLC) = tmp(iterate_DLC-1);
        end
    end
    
    mean_b = nanmedian(recording.movie1.abd_length);
    tmp = recording.movie1.abd_length;
    for iterate_DLC = 2:1:length(recording.movie1.abd_length)
        if(abs(tmp(iterate_DLC-1) -  tmp(iterate_DLC)) > mean_b/4)
            %recording.movie1.abd_length(iterate_DLC) = recording.movie1.abd_length(iterate_DLC-1);
            recording.movie1.abd_length(iterate_DLC) = NaN;
            tmp(iterate_DLC) = tmp(iterate_DLC-1);
        end
    end
    
    mean_b = nanmedian(recording.movie1.abd_only_length);
    tmp = recording.movie1.abd_only_length;
    for iterate_DLC = 2:1:length(recording.movie1.abd_length)
        if(abs(tmp(iterate_DLC-1) -  tmp(iterate_DLC)) > mean_b/4)
            % recording.movie1.abd_only_length(iterate_DLC) = recording.movie1.abd_only_length(iterate_DLC-1);
            recording.movie1.abd_only_length(iterate_DLC) = NaN;
            tmp(iterate_DLC) = tmp(iterate_DLC-1);
        end
    end
    
else
    disp 'mismatch in lengths, nothing was done'
    
end
end
