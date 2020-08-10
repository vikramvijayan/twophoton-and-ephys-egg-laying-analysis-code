function out = query_klusterimage_zproj(tsStack, query_frames, transform1, transform2)

warning('off');
figure;
c = 1;
uniq_kl = max(query_frames(:));
for i =1:1:(uniq_kl)
    smplot(1,(uniq_kl),c); hold on;
    [a b] = find(query_frames(:) == i);
    
    if(~isempty(a))
        
        if(strcmp(transform1,'max'))
            tmp = max(tsStack(:,:,1,a,:),[],5);
            tmp_all = max(tsStack(:,:,1,:,:),[],5);
        end
        if(strcmp(transform1,'mean'))
            tmp = mean(tsStack(:,:,1,a,:),5);
            tmp_all = mean(tsStack(:,:,1,:,:),5);
        end
        
        tmp = squeeze(tmp);
        tmp_all = squeeze(tmp_all);
        
        if(strcmp(transform2,'max'))
            data = max(tmp,[],3);
            data_all = max(tmp_all,[],3);
        end
        if(strcmp(transform2,'mean'))
            data = mean(tmp,3);
            data_all = mean(tmp_all,3);
        end
        
        %       imshow((data),[0 2^16]); colormap('jet');
        %imagesc((data-data_all)); colormap('jet');
        % imagesc((data)); colormap('jet');
        xlabel(num2str(    length(a)));
        imagesc(flipud(log2(data./data_all))); colormap('jet'); caxis([-.05,.05]);
    end
        axis equal;   axis off;
    
    c = c+1;
end

out =1;
warning('on');