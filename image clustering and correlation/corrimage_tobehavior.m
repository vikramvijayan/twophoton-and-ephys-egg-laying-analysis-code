function [corr_images] = corrimage_tobehavior(tsStack,behavior, titlefff)
warning('off');

% clusters the pixel by pixel correlations from tsStack and returns the cluster
% use query_klusterimage to plit the clusters

ts_size = size(tsStack);


if(length(ts_size) == 5)
    for i = 1:1:ts_size(5)
        
        current_plane = tsStack(:,:,1,:,i);
        current_plane = squeeze(current_plane);
        
        current_plane_size = size(current_plane);
        
        reshape_current_plane = reshape(current_plane,[current_plane_size(1)*current_plane_size(2),current_plane_size(3)]);
        
        
        concant_beh = [behavior'; reshape_current_plane];
        out1 = corrcoef(double(concant_beh)','rows','complete');
        out1 = out1(2:end,1);
        
        out2 = out1(:);
        
        out3 = reshape(out2,[current_plane_size(1),current_plane_size(2)]);
        
        corr_images(:,:,i) = out3;
        
        figure; hold on;
        title([titlefff ' slice ' num2str(i)]);
        imagesc(flipud(corr_images(:,:,i))); colorbar; colormap('jet'); axis tight; axis equal;
        
        warning('on');
        
    end
    
else
    
    current_plane = tsStack(:,:,1,:);
    current_plane = squeeze(current_plane);
    
    current_plane_size = size(current_plane);
    
    reshape_current_plane = reshape(current_plane,[current_plane_size(1)*current_plane_size(2),current_plane_size(3)]);
    
    
    concant_beh = [behavior'; reshape_current_plane];
    out1 = corrcoef(double(concant_beh)','rows','complete');
    out1 = out1(2:end,1);
    
    out2 = out1(:);
    
    out3 = reshape(out2,[current_plane_size(1),current_plane_size(2)]);
    
    corr_images(:,:) = out3;
    
    figure; hold on;
    title([titlefff]);
    imagesc(flipud(corr_images(:,:))); colorbar; colormap('jet'); axis equal; axis tight;
    if(ispublishing())
        snapnow
        close;
    end
    warning('on');
end
