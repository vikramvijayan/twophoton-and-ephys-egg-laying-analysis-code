function [km] = klusterimage(tsStack,k)

% clusters images from tsStack and returns the cluster
% use query_klusterimage to plit the clusters

ts_size = size(tsStack);

for i = 1:1:ts_size(5)
    current_plane = tsStack(:,:,1,:,i);
    current_plane = squeeze(current_plane);
    current_plane_size = size(current_plane);

    reshape_current_plane = reshape(current_plane,[current_plane_size(1)*current_plane_size(2),current_plane_size(3)]);
    
    km(i,:) = kmeans(double(reshape_current_plane)',k);
end