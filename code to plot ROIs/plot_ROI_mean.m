function out = plot_ROI_mean(tsStack, ROI)

%figure; hold on; axis off;
%suptitle(tseries.name);


warning('off');

size_ROI = size(ROI);
%size_tsStack = size(tsStack);

ROI_number = size_ROI(4);
zslices = size_ROI(3);



c = 1;
for i = 1:1:ROI_number
    %subplot(ROI_number,1,c); hold on; axis off;
    figure; hold on; axis off; %title([titlef ' ROI ' num2str(i)]);
    tmp = mean(tsStack(:,:,1,:,:),5);
    
    imagesc(flipud(mean(tmp,4)));
    set(gca,'YDir','normal')
    colormap('bone');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    for j = 1:1:zslices
        
        b = bwboundaries(flipud(ROI(:,:,j,i)));
        for k = 1:numel(b)
            plot(b{k}(:,2), b{k}(:,1), 'r', 'Linewidth', 2)
        end
        %   BW2 = bwperim(ROI(:,:,j,i));
        axis equal
    end
    c = c+1;
    
    if(ispublishing())
        snapnow;
        close;
    end
end
warning('on');
