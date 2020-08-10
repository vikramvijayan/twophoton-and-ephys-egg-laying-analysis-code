function recording = extract_2p_and_patch_data_egg_redo_tsStack(tsStack_cell, ROI_cell, recording)

%currently the alst 2 csvs are hard coded
warning('off');

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

warning('on');
end