function [km] = klusterimage_zproj(tsStack, transform,k)

% clusters images from tsStack (the full zp proj using either mean or max transform) and returns the cluster
% use query_klusterimage_zproj to plit the clusters

if(strcmp(transform,'max'))
    current_plane = max(tsStack(:,:,1,:,:),[],5);
end
if(strcmp(transform,'mean'))
    current_plane = mean(tsStack(:,:,1,:,:),5);
end
current_plane = squeeze(current_plane);
current_plane_size = size(current_plane);
reshape_current_plane = reshape(current_plane,[current_plane_size(1)*current_plane_size(2),current_plane_size(3)]);

mean_for_df = mean(double(reshape_current_plane)')';

df = double(reshape_current_plane) - double(repmat(mean_for_df,1,size(reshape_current_plane,2)));
 df_over_f = df ./ double(repmat(mean_for_df,1,size(reshape_current_plane,2)));
 
km = kmeans(double(df_over_f)',k);
end