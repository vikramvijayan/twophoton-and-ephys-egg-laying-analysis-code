function [corr_images, km] = corrimage(tsStack,k, plot_clusters)

% clusters the pixel by pixel correlations from tsStack and returns the cluster
% use query_klusterimage to plit the clusters

ts_size = size(tsStack);

    
    
for i = 1:1:ts_size(5)
    
    i
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
    
    
    if(plot_clusters)
        figure; hold on;

        for num_clus = 1:1:k
            [a b] = find(km(i,:) == num_clus);
            hold on;
            subplot(k, 1,num_clus);
            imagesc(mean(out3(:,:,b),3)); colorbar; colormap('jet'); axis equal; axis tight;
        end
        
        a = axes;
        t1 = title(['klusters total ' num2str(k) ' zslice ' num2str(i)]);
        a.Visible = 'off'; % set(a,'Visible','off');
        t1.Visible = 'on'; % set(t1,'Visible','on');

    end
    
    
    %
    % for i = 1:1:qsize(1)
    %     for j = 1:1:qsize(2)
    %         out(i,j) = xcorr(squeeze(query), squeeze(q(i,j,:)));
    %     end
    % end
end


