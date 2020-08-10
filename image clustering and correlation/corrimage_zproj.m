function [corr_images, km] = corrimage_zproj(tsStack,k,transform)

if(strcmp(transform,'max'))
    current_plane = max(tsStack(:,:,1,:,:),[],5);
end
if(strcmp(transform,'mean'))
    current_plane = mean(tsStack(:,:,1,:,:),5);
end


ts_size = size(tsStack);

i =1;

current_plane = tsStack(:,:,1,:,i);
current_plane = squeeze(current_plane);

current_plane_size = size(current_plane);

%     for j = 1:1:current_plane_size(3);
%         current_plane(:,:,j) = imgaussfilt(current_plane(:,:,j),4);
%     end

reshape_current_plane = reshape(current_plane,[current_plane_size(1)*current_plane_size(2),current_plane_size(3)]);

out1 = corrcoef(double(reshape_current_plane'));

out2 = out1(:);

out3 = reshape(out2,[current_plane_size(1),current_plane_size(2),current_plane_size(1)*current_plane_size(2)]);

out4 = reshape(out2,[current_plane_size(1)*current_plane_size(2),current_plane_size(1)*current_plane_size(2)]);

km(i,:) = kmeans(out4,k);
%km(i,:) = clusterdata(out4,'maxclust',10);

corr_images(:,:,1,:,i) = out3;

%
% for i = 1:1:qsize(1)
%     for j = 1:1:qsize(2)
%         out(i,j) = xcorr(squeeze(query), squeeze(q(i,j,:)));
%     end
% end



