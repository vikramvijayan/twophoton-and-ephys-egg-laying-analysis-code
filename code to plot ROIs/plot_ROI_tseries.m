function out = plot_ROI_tseries(tseries)

figure; hold on; axis off;
%suptitle(tseries.name);

ROI = tseries.ROI;
tsStack = tseries.tsStack;

warning('off');

size_ROI = size(ROI);
%size_tsStack = size(tsStack);

ROI_number = size_ROI(4);
zslices = tseries.zslices;


c = 1;
for i = 1:1:ROI_number
    for j = 1:1:zslices
        subplot(ROI_number,zslices,c); hold on; axis off;

        imagesc(flipud(max(tsStack(:,:,1,:,j),[],4)));
        set(gca,'YDir','normal')
        colormap('bone');
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
