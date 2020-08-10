% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list] = import_ROI_names_SS53469imageR65H10chrimson(ROIsSS53469, whatdoyouwant)

%function [recording] = make_small_recording_group()

recording_list = {'2019_11_13_0000.mat', '2019_11_13_0001.mat', '2019_11_13_0002.mat', '2019_11_13_0003.mat'};

% imported ROIsSS53469 as string

for i = 1:1:length(recording_list)
    % brighteset to dimmest
    ROI(i).cell1 = str2num(ROIsSS53469(i,1));
    ROI(i).cell2 = str2num(ROIsSS53469(i,2));
    ROI(i).cell3 = str2num(ROIsSS53469(i,3));
    ROI(i).cell4 = str2num(ROIsSS53469(i,4));
    ROI(i).cell5 = str2num(ROIsSS53469(i,5));
    ROI(i).cell6 = str2num(ROIsSS53469(i,6));
    
    
end


ROI_list = [];
rec_list = [];
cnt = 0;

if(strcmp(whatdoyouwant,'cell'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).cell1)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).cell1(j);
        end
        for j = 1:1:length(ROI(i).cell2)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).cell2(j);
        end
        for j = 1:1:length(ROI(i).cell3)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).cell3(j);
        end
        for j = 1:1:length(ROI(i).cell4)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).cell4(j);
        end
        for j = 1:1:length(ROI(i).cell5)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).cell5(j);
        end
        for j = 1:1:length(ROI(i).cell6)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).cell6(j);
        end
    end
end


if(strcmp(whatdoyouwant,'only_brightest_cell'))
    for i = 1:1:length(recording_list)
        cell_taken = 0;
        for j = 1:1:length(ROI(i).cell1)
            if(cell_taken == 0)
                cnt = cnt+1;
                rec_list{cnt} = char(recording_list(i));
                ROI_list(cnt) = ROI(i).cell1(j);
                cell_taken=1;
            end
        end
        
    end
end





ROI_names =ROI;

end