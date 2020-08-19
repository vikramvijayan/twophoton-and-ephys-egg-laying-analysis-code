function recording = extract_2p_and_patch_data_egg(tsStack_cell, ROI_cell, abf_filename, interactive_mask_for_2p_trig, movie1, movie2, alignmovie1, alignmovie2, substrate_vector, DLC_csv_abdomen,DLC_csv)

%currently the alst 2 csvs are hard coded
warning('off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% abf is a structure within 'recording'
% It keeps all the data (or interpolates it) at the abf framerate which is
% typically 100 us (this is automatically determined from the abf metadata)
[d,si,h] = abfload(abf_filename);

[rws,cls] = size(d);
recording.abfname = abf_filename;

% this is an older version of the abf imported without the two Patch
% channels

if(cls == 13)
    recording.abf.CH1_patch = zeros(length(d(:,1)),1);
    recording.abf.CH2_patch = zeros(length(d(:,1)),1);
    recording.abf.FrameSt2p = d(:,3-2);
    recording.abf.Pockels = d(:,14-2);
    recording.abf.PWMlaser = d(:,13-2);
    recording.abf.Time_s = (0:si/1000000:(length(recording.abf.FrameSt2p)-1)*(si/1000000))';
    recording.abf.Puff = d(:,12-2);
    recording.abf.Temp = d(:,4-2);
    recording.abf.WalkCamTri = d(:,11-2);
    recording.abf.Wheel = d(:,15-2);
    recording.abf.xstim = d(:,9-2);
    recording.abf.ystim = d(:,10-2);
    recording.abf.current_mode = d(:,5-2);
else
    
    
    recording.abf.CH1_patch = d(:,1);
    recording.abf.CH2_patch = d(:,2);
    recording.abf.FrameSt2p = d(:,3);
    recording.abf.Pockels = d(:,14);
    recording.abf.PWMlaser = d(:,13);
    recording.abf.Time_s = (0:si/1000000:(length(recording.abf.FrameSt2p)-1)*(si/1000000))';
    recording.abf.Puff = d(:,12);
    recording.abf.Temp = d(:,4);
    recording.abf.WalkCamTri = d(:,11);
    recording.abf.Wheel = d(:,15);
    recording.abf.xstim = d(:,9);
    recording.abf.ystim = d(:,10);
    recording.abf.current_mode = d(:,5);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find when the walk cam is trigerred (assuming down triggered, CHECK!)

% note there is delay between the end of the trigger and the output
% calculation (from fictrac) that isnt taken care of (apperas to be 5 to 10 ms)
% also I do not place the behavior camera in between the triggers. That is,
% there could be another 10 to 20 ms (at 50Hz) delay in behavior

% [a b] = find(diff(recording.abf.WalkCamTri)> 1.5);
% [a1 b1] = find(diff(recording.abf.WalkCamTri) < -1.5);

% using a different Arduino trigger
[a b] = find(diff(recording.abf.WalkCamTri) > 1);
[a1 b1] = find(diff(recording.abf.WalkCamTri) < -1);

a1_new = a1(1);
for i = 1:1:(length(a1)-1)
    if((a1(i+1)-a1(i)) > 5)
        a1_new = [a1_new a1(i+1)];
    end
end

a1 = a1_new;

recording.cam.trigger_start = a1;
recording.cam.median_diff = median(diff(a1));
recording.cam.Time_s = recording.abf.Time_s(a1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% makes an interactive mask for selecting 2p triggers
% (select before and then after the triggers of interest)

tmp_cnt = 1;
mask_tmp = [];
mask = [];
orient_fig =  figure; plot(recording.abf.FrameSt2p,'r');
if(interactive_mask_for_2p_trig)
    while(tmp_cnt < length(recording.abf.FrameSt2p))
        tmp_fig = figure; plot(recording.abf.FrameSt2p(tmp_cnt:min(tmp_cnt+10000000,length(recording.abf.FrameSt2p))));
        [mask_tmp,~]=ginput;
        mask = [mask; floor(mask_tmp)+tmp_cnt-1];
        tmp_cnt = tmp_cnt+10000000;
        close(tmp_fig);
    end
    close(orient_fig);
    
    mask_applied = zeros(length(recording.abf.FrameSt2p),1);
    for i =1:2:length(mask)
        mask_applied(mask(i):mask(i+1)) = 1;
    end
else
    mask_applied = ones(length(recording.abf.FrameSt2p),1);
end


% finding the time of all 2 photon frame triggers
[a b] = find(diff(recording.abf.FrameSt2p.*mask_applied)> 2.4);
[a1 b1] = find(diff(recording.abf.FrameSt2p.*mask_applied) < -2.4);

if(length(a1) > 0)
    a1_new = a1(1);
    for i = 1:1:(length(a1)-1)
        if((a1(i+1)-a1(i)) > 5)
            a1_new = [a1_new a1(i+1)];
        end
    end
    a1 = a1_new;
    
    
    a_new = a(1);
    for i = 1:1:(length(a)-1)
        if((a(i+1)-a(i)) > 5)
            a_new = [a_new a(i+1)];
        end
    end
    a = a_new;
end

% whenever the 2ptrigger goes up and then down and there are pockels, there is a frame
% this used to be looking between down triggers, now it looks between up
% and down to get better image placement

trig = [];
for i =1:1:(length(a1))
    tmp_mean = mean(recording.abf.Pockels(a(i):a1(i)));
    if(tmp_mean > 0.005) % need to use this smaller number if the exposure is fast
        trig = [trig (recording.abf.Time_s(a(i)) + recording.abf.Time_s(a1(i)))/2];
    end
end

% may need to adjust the trig array if there were aborted tseries etc

%trig = [trig(1:(1250*15)),trig((end-15*2000+1):(end))];
%trig = [trig(1:(17*4042))];

recording.abf.middle_of_trig = trig';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go through each of the 2p time series made in this session and compute
% ROI values for etc.

number_of_tseries = length(ROI_cell);

for i = 1:1:number_of_tseries
    ROI = cell2mat(ROI_cell(i));
    tsStack = cell2mat(tsStack_cell(i));
    size_ROI = size(ROI);
    size_tsStack = size(tsStack);
    
    recording.tseries(i).name = [recording.abfname ' tseries ' num2str(i)];
    
    recording.tseries(i).ROI_number = size_ROI(4);
    recording.tseries(i).ROI = ROI;
    recording.tseries(i).tsStack = tsStack;
    
    if(length(size_tsStack) == 5)
        recording.tseries(i).Total2pFrames = size_tsStack(4);
        recording.tseries(i).zslices = size_tsStack(5);
        recording.tseries(i).mean_in_ROI = zeros(size_ROI(4),size_tsStack(4));
    else
        recording.tseries(i).Total2pFrames = size_tsStack(3);
        recording.tseries(i).zslices = 1;
        recording.tseries(i).mean_in_ROI = zeros(size_ROI(4),size_tsStack(3));
    end
    
    % calculate df over f etc in different ROIs
    for k =1:1:recording.tseries(i).ROI_number
        for j = 1:1:recording.tseries(i).Total2pFrames
            if(length(size_tsStack) == 5)
                tmp = tsStack(:,:,1,j,:);
                tmp = squeeze(tmp);
                tmp_masked_data = tmp(ROI(:,:,:,k) == 1);
                
                [a b]  = find(ROI(:,:,:,k) == 1);
                recording.tseries(i).mean_in_ROI(k,j) = sum(tmp_masked_data)./length(a);
                if(recording.tseries(i).mean_in_ROI(k,j) == 0)
                    recording.tseries(i).mean_in_ROI(k,j) = NaN;
                end
                recording.tseries(i).sum_in_ROI(k,j) = sum(tmp_masked_data);
                recording.tseries(i).median_in_ROI(k,j) = median(tmp_masked_data);
                
                %go to each z slice for a given ROI and computes stats
                %seperately
                for m = 1:1:size_ROI(3)
                    tmpm = tmp(:,:,m);
                    tmp_masked_data = tmpm(ROI(:,:,m,k) == 1);
                    [a b] = find(ROI(:,:,m,k) == 1);
                    recording.tseries(i).mean_in_ROI_indep_z(k,j,m) = sum(tmp_masked_data)./length(a);
                    recording.tseries(i).sum_in_ROI_indep_z(k,j,m) = sum(tmp_masked_data);
                    recording.tseries(i).median_in_ROI_indep_z(k,j,m) = median(tmp_masked_data);
                    
                end
                
                %fit a gaussian to the means replacing NAN with 0
                tofit = squeeze(recording.tseries(i).mean_in_ROI_indep_z(k,j,:));
                tofit(isnan(tofit)) = 0;
                %f = fit((1:1:recording.tseries(i).zslices)',tofit, 'gauss1');
                recording.tseries(i).mean_in_ROI_gauss(k,j) = 0;
                
                % this is for 2 channels
                if( (size_tsStack(3)) == 2)
                    tmp = tsStack(:,:,2,j,:);
                    tmp = squeeze(tmp);
                    tmp_masked_data = tmp(ROI(:,:,:,k) == 1);
                    [a b]  = find(ROI(:,:,:,k) == 1);
                    recording.tseries(i).mean_in_ROI2(k,j) = sum(tmp_masked_data)./length(a);
                    if(recording.tseries(i).mean_in_ROI2(k,j) == 0)
                        recording.tseries(i).mean_in_ROI2(k,j) = NaN;
                    end
                    recording.tseries(i).sum_in_ROI2(k,j) = sum(tmp_masked_data);
                    recording.tseries(i).median_in_ROI2(k,j) = median(tmp_masked_data);
                    
                    for m = 1:1:size_ROI(3)
                        tmpm = tmp(:,:,m);
                        tmp_masked_data = tmpm(ROI(:,:,m,k) == 1);
                        [a b] = find(ROI(:,:,m,k) == 1);
                        recording.tseries(i).mean_in_ROI_indep_z2(k,j,m) = sum(tmp_masked_data)./length(a);
                        recording.tseries(i).sum_in_ROI_indep_z2(k,j,m) = sum(tmp_masked_data);
                        recording.tseries(i).median_in_ROI_indep_z2(k,j,m) = median(tmp_masked_data);
                        
                    end
                    
                    %fit a gaussian to the means replacing NAN with 0
                    tofit = squeeze(recording.tseries(i).mean_in_ROI_indep_z2(k,j,:));
                    tofit(isnan(tofit)) = 0;
                    %f = fit((1:1:recording.tseries(i).zslices)',tofit, 'gauss1');
                    recording.tseries(i).mean_in_ROI_gauss2(k,j) = 0;
                    
                end
            else
                tmp = tsStack(:,:,j);
                tmp = squeeze(tmp);
                tmp_masked_data = tmp(ROI(:,:,1,k) == 1);
                [a b]  = find(ROI(:,:,1,k) == 1);
                recording.tseries(i).mean_in_ROI(k,j) = sum(tmp_masked_data)./length(a);
                recording.tseries(i).sum_in_ROI(k,j) = sum(tmp_masked_data);
                recording.tseries(i).median_in_ROI(k,j) = median(tmp_masked_data);
                
            end
        end
    end
    
    for  k =1:1:recording.tseries(i).ROI_number
        recording.tseries(i).df_over_f_usingmean(k,:) = (recording.tseries(i).mean_in_ROI(k,:)-mean(recording.tseries(i).mean_in_ROI(k,:)))./mean(recording.tseries(i).mean_in_ROI(k,:));
        low_val = prctile(recording.tseries(i).mean_in_ROI(k,:),5);
        recording.tseries(i).df_over_f(k,:) = (recording.tseries(i).mean_in_ROI(k,:)-low_val)./low_val;
        
        low_val = prctile(recording.tseries(i).mean_in_ROI_gauss(k,:),5);
        recording.tseries(i).df_over_f_gauss(k,:) = (recording.tseries(i).mean_in_ROI_gauss(k,:)-low_val)./low_val;
        
        
        %for each z slice indep
        
        % take max sum of the mean z slice
        low_val = prctile(max(recording.tseries(i).mean_in_ROI_indep_z(k,:,:),[],3),5);
        recording.tseries(i).df_over_f_maxofmean(k,:) = (max(recording.tseries(i).mean_in_ROI_indep_z(k,:,:),[],3)-low_val)./low_val;
        
        % take mean of max (and one on each side) of the mean z slice
        [~, max_slice] = max(recording.tseries(i).mean_in_ROI_indep_z(k,:,:),[],3);
        
        [a1 b1] = find(~isnan(recording.tseries(i).mean_in_ROI_indep_z(k,1,:)),1,'first');
        [a2 b2] = find(~isnan(recording.tseries(i).mean_in_ROI_indep_z(k,1,:)),1,'last');
        
        % make sure you dont go beyond bounds
        bottom_slice = max(b1.*ones(1,length(max_slice)), max_slice-1);
        top_slice = min(b2.*ones(1,length(max_slice)), max_slice+1);
        
        low_val = prctile(mean(recording.tseries(i).mean_in_ROI_indep_z(k,:,(bottom_slice):(top_slice)),3),5);
        recording.tseries(i).df_over_f_3maxofmean(k,:) = (mean(recording.tseries(i).mean_in_ROI_indep_z(k,:,(bottom_slice):(top_slice)),3)-low_val)./low_val;
        
        % take max of the median z slice
        low_val = prctile(max(recording.tseries(i).median_in_ROI_indep_z(k,:,:),[],3),5);
        recording.tseries(i).df_over_f_maxofmed(k,:) = double((max(recording.tseries(i).median_in_ROI_indep_z(k,:,:),[],3)-low_val))./double(low_val);
        
        if( (size_tsStack(3)) == 2)
            recording.tseries(i).df_over_f_usingmean2(k,:) = (recording.tseries(i).mean_in_ROI2(k,:)-mean(recording.tseries(i).mean_in_ROI2(k,:)))./mean(recording.tseries(i).mean_in_ROI2(k,:));
            low_val = prctile(recording.tseries(i).mean_in_ROI2(k,:),5);
            recording.tseries(i).df_over_f2(k,:) = (recording.tseries(i).mean_in_ROI2(k,:)-low_val)./low_val;
            
            low_val = prctile(recording.tseries(i).mean_in_ROI_gauss(k,:),5);
            recording.tseries(i).df_over_f_gauss2(k,:) = (recording.tseries(i).mean_in_ROI_gauss2(k,:)-low_val)./low_val;
            
            %for each z slice indep
            
            % take max of the mean z slice
            low_val = prctile(max(recording.tseries(i).mean_in_ROI_indep_z2(k,:,:),[],3),5);
            recording.tseries(i).df_over_f_maxofmean2(k,:) = (max(recording.tseries(i).mean_in_ROI_indep_z2(k,:,:),[],3)-low_val)./low_val;
            
            % take mean of max (and one on each side) of the mean z slice
            [~, max_slice] = max(recording.tseries(i).mean_in_ROI_indep_z2(k,:,:),[],3);
            
            % make sure you dont go beyond bounds
            [a1 b1] = find(~isnan(recording.tseries(i).mean_in_ROI_indep_z(k,1,:)),1,'first');
            [a2 b2] = find(~isnan(recording.tseries(i).mean_in_ROI_indep_z(k,1,:)),1,'last');
            
            bottom_slice = max(b1.*ones(1,length(max_slice)), max_slice-1);
            top_slice = min(b2.*ones(1,length(max_slice)), max_slice+1);
            
            low_val = prctile(mean(recording.tseries(i).mean_in_ROI_indep_z2(k,:,(bottom_slice):(top_slice)),3),5);
            recording.tseries(i).df_over_f_3maxofmean2(k,:) = (mean(recording.tseries(i).mean_in_ROI_indep_z2(k,:,(bottom_slice):(top_slice)),3)-low_val)./low_val;
            
            % take max of the median z slice
            low_val = prctile(max(recording.tseries(i).median_in_ROI_indep_z2(k,:,:),[],3),5);
            recording.tseries(i).df_over_f_maxofmed2(k,:) = double((max(recording.tseries(i).median_in_ROI_indep_z2(k,:,:),[],3)-low_val))./double(low_val);
        end
    end
    
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% go through each of the tseries and make sure the images are assigned a
% time using abf triggers.

% we will use the average time for the full zseries as the time that
% zseries was acquired. essentially we are placing the full z projection at
% the middle timepoint
total_2p_frames_tmp = 0;
cnt = 1;
cnt2 = 1;
recording.abf.middle_of_trig_ave_for_z=[];
recording.abf.start_of_trig_ave_for_z=[];
recording.abf.end_of_trig_ave_for_z=[];

for i = 1:1:number_of_tseries
    
    total_2p_frames_tmp = total_2p_frames_tmp + recording.tseries(i).zslices.*recording.tseries(i).Total2pFrames;
    while( ((cnt+recording.tseries(i).zslices-1) <= length(recording.abf.middle_of_trig)) && ((cnt+recording.tseries(i).zslices-1) <= total_2p_frames_tmp))
        recording.abf.middle_of_trig_ave_for_z(cnt2) = mean(recording.abf.middle_of_trig(cnt:(cnt+recording.tseries(i).zslices-1)));
        recording.abf.start_of_trig_ave_for_z(cnt2) = recording.abf.middle_of_trig(cnt);
        recording.abf.end_of_trig_ave_for_z(cnt2)   = recording.abf.middle_of_trig(cnt+recording.tseries(i).zslices-1);
        cnt = cnt + recording.tseries(i).zslices;
        cnt2 = cnt2+1;
    end
end

if(length(recording.abf.middle_of_trig) ~= total_2p_frames_tmp)
    error( 'There is inconsistency in the number of triggers and number of frames');
end


% the closest index to an abf sampling that matches all of the two photon
% frames (that is, what abf index was the 2p image captured at)

for i =1:1:length(recording.abf.middle_of_trig_ave_for_z)
    recording.abf.twoPimages_index(i) = ceil(recording.abf.middle_of_trig_ave_for_z(i)*(1000000/si));
    recording.abf.twoPimages_index_start(i) = ceil(recording.abf.start_of_trig_ave_for_z(i)*(1000000/si));
    recording.abf.twoPimages_index_end(i) = ceil(recording.abf.end_of_trig_ave_for_z(i)*(1000000/si));
end

% find the median imaging interval for each tseries in abf units
for i = 1:1:number_of_tseries
    recording.tseries(i).median_imaging_interval = median( recording.abf.twoPimages_index_end -  recording.abf.twoPimages_index_start);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for each tseries, save values from abf at the time that 2p images were taken
offset = 1;

for i = 1:1:number_of_tseries
    offset_end = offset + recording.tseries(i).Total2pFrames -1;
    
    % these are higher res arrays that are at abf timestep
    first_index = recording.abf.twoPimages_index(offset);
    last_index = recording.abf.twoPimages_index(offset_end);
    
    recording.tseries(i).high_res_index = first_index:last_index;
    
    % these are the arrays that save data from the abf at the middle of the
    % 2 photon image
    
    recording.tseries(i).middle_index = recording.abf.twoPimages_index(offset:offset_end);
    
    recording.tseries(i).Time_s = recording.abf.Time_s(recording.abf.twoPimages_index(offset:offset_end));
    recording.tseries(i).xstim = recording.abf.xstim(recording.abf.twoPimages_index(offset:offset_end));
    recording.tseries(i).ystim = recording.abf.ystim(recording.abf.twoPimages_index(offset:offset_end));
    recording.tseries(i).Temp = recording.abf.Temp(recording.abf.twoPimages_index(offset:offset_end));
    recording.tseries(i).current_mode = recording.abf.current_mode(recording.abf.twoPimages_index(offset:offset_end));
    recording.tseries(i).Wheel = recording.abf.Wheel(recording.abf.twoPimages_index(offset:offset_end));
    recording.tseries(i).CH1_patch = recording.abf.CH1_patch(recording.abf.twoPimages_index(offset:offset_end));
    recording.tseries(i).CH2_patch = recording.abf.CH2_patch(recording.abf.twoPimages_index(offset:offset_end));
    
    offset = offset_end+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load movies
if(~isempty(movie1))
    recording.movie1.filename = movie1;
    modifiedStr = strrep([recording.movie1.filename '_processed'], '.', '_');
    
    % if(strcmp(movie1(end-2:end),'fmf'))
    % add the processed data to the movie structure
    
    if exist([modifiedStr '.mat'])
        load([modifiedStr '.mat']);
        recording.movie1.time_stamps = stamp;
        recording.movie1.eggs = egg_array_num;
    else
        disp 'No movie 1'
    end
    %  end
end


if(~isempty(movie2))
    recording.movie2.filename = movie2;
    modifiedStr = strrep([recording.movie2.filename '_processed'], '.', '_');
    
    %  if(strcmp(movie2(end-2:end),'fmf'))
    % add the processed data to the movie structure
    
    if exist([modifiedStr '.mat'])
        load([modifiedStr '.mat']);
        recording.movie2.time_stamps = stamp;
        recording.movie2.bound = bound;
    else
        disp 'No movie 2'
    end
    %  end
end



% if there is a movie find the alignment points
% find first 1 sec pause in triggers and use this to align

% this is new code that is aliging to the end of the movie. can't use the
% first pause anymore because of the fact we dont have triggers
if(~isempty(movie1))
    [a b] = find(diff(recording.cam.Time_s) > .25,1,'first');
    if(alignmovie1 > 0)
        recording.movie1.abf_time = recording.cam.Time_s((a+1-alignmovie1+1) : ((a+1-alignmovie1+1)+length(recording.movie1.time_stamps)-1));
        
    else
        modifiedStr = strrep([recording.movie1.filename(1:1:end-6) 'DeepCut_resnet50_GPU_05062019May6shuffle1_1030000.csv'], '.avi', '');
        
        
        [bodypart]= importDeepLabCutfile_VV(modifiedStr, 4, inf);
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
        
%         mean_b = nanmedian(recording.movie1.prob_length);
%         tmp = recording.movie1.prob_length;
%         for iterate_DLC = 2:1:length(recording.movie1.abd_length)
%             if(abs(tmp(iterate_DLC-1) -  tmp(iterate_DLC)) > mean_b/4)
%                 % recording.movie1.abd_only_length(iterate_DLC) = recording.movie1.abd_only_length(iterate_DLC-1);
%                 recording.movie1.prob_length(iterate_DLC) = NaN;
%                 tmp(iterate_DLC) = tmp(iterate_DLC-1);
%             end
%         end
        
        recording.movie1.time_stamps = recording.cam.Time_s((end-length(recording.movie1.time_stamps)+1):1:end);
        disp(['movie1 extra DLC ' num2str(length(bodypart(2,:)) - length(recording.movie1.time_stamps))]);
        disp(['movie1 cam extra triggers ' num2str(length(recording.cam.Time_s) - length(recording.movie1.time_stamps))]);
        
    end
end


if(~isempty(movie2))
    [a b] = find(diff(recording.cam.Time_s) > .25,1,'first');
    if(alignmovie2 > 0)
        recording.movie2.abf_time = recording.cam.Time_s((a+1-alignmovie2+1) : ((a+1-alignmovie2+1)+length(recording.movie2.time_stamps)-1));
    else
        modifiedStr = strrep([recording.movie2.filename(1:1:end-6) 'DeepCut_resnet50_GPU_05062019BMay6shuffle1_1030000.csv'], '.avi', '');
        [bodypart]= importDeepLabCutfile_VV(modifiedStr, 4, inf);
        recording.movie2.DLC = bodypart(2:13,:);
        recording.movie2.time_stamps = recording.cam.Time_s((end-length(recording.movie2.time_stamps)+1):1:end);
    end
    disp(['movie2 extra DLC ' num2str(length(bodypart(2,:)) - length(recording.movie2.time_stamps))]);
    disp(['movie2 extra triggers ' num2str(length(recording.cam.Time_s) - length(recording.movie2.time_stamps))]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter the wheel trace and find substrate transitions

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

warning('on');
end