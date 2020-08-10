% sample run if this code:
% [ k ] = split_images_by_heading(recording, 2, [0,pi-pi/8; pi+pi/8, 2*pi],1,recording.tseries(2).Time_s(1),500,0);
%query_klusterimage_zproj(recording.tseries(2).tsStack,[k],'max','mean')



% split frames based on current heading
% onbtl frames that pass the fly speed lower bound in cm/sec are used
% the last kluster is all the images that didnt pass the threshold

function [ k ] = split_images_by_heading(recording, t_num, split_angles, speed_lower_bound, ts, te, offset_head, exclude_dark )

[a b]   = find(recording.tseries(t_num).Time_s >= ts,1, 'first');
[a1 b1] = find(recording.tseries(t_num).Time_s <= te ,1,'last');



[r, c] = size(split_angles);

k = (r+1).*ones(1,length(recording.tseries(t_num).actual_head));


if(~exclude_dark)
    for i = 1:1:r
        current = split_angles(i,1);
        next = split_angles(i,2);
        
        [a2 b2] = find( ( mod((offset_head+recording.tseries(t_num).actual_head(a:a1)),2*pi) >= current) & (mod((offset_head+recording.tseries(t_num).actual_head(a:a1)),2*pi) < (next)) & (recording.tseries(t_num).fvel_inst(a:a1) > speed_lower_bound));
        
        k(a2+a-1) = i;
        
    end
end



if(exclude_dark)
    for i = 1:1:r
        current = split_angles(i,1);
        next = split_angles(i,2);
        
                [a2 b2] = find(  ((recording.tseries(t_num).actual_head(a:a1) < (pi-pi/8)) | (recording.tseries(t_num).actual_head(a:a1) > (pi+pi/8))) & ( mod((offset_head+recording.tseries(t_num).actual_head(a:a1)),2*pi) >= current) & (mod((offset_head+recording.tseries(t_num).actual_head(a:a1)),2*pi) < (next)) & (recording.tseries(t_num).fvel_inst(a:a1) > speed_lower_bound));

        
        k(a2+a-1) = i;
        
    end
end




end

