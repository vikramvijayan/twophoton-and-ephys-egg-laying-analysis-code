function out = plot_ROI(tsStack, ROI)
warning('off');

size_ROI = size(ROI);
size_tsStack = size(tsStack);

ROI_number = size_ROI(4);
zslices = size_tsStack(5);

figure;

c = 1;
for i = 1:1:ROI_number
    for j = 1:1:zslices
        smplot(ROI_number,zslices,c); hold on; axis off;
       % imagesc(flipud(max(tsStack(:,:,1,:,j),[],4)));
                imagesc(flipud(mean(tsStack(:,:,1,:,j),4)));
        set(gca,'YDir','normal')
        %colormap('bone');
                colormap('jet');

        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        
        b = bwboundaries(flipud(ROI(:,:,j,i)));
        for k = 1:numel(b)
            plot(b{k}(:,2), b{k}(:,1), 'r', 'Linewidth', 2)
        end
        %   BW2 = bwperim(ROI(:,:,j,i));
        axis equal
        c = c+1;
    end
end
warning('on');
