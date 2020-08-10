function out = query_klusterimage(tsStack, query_frames,transform)

[zs ~] = size(query_frames);
warning('off');
figure;
c = 1;
for j =1:1:zs
    uniq_kl = unique(query_frames(j,:));
    for i =1:1:length(uniq_kl)
        smplot(zs,length(uniq_kl),c);
        [a b] = find(query_frames(j,:) == i);
        tmp = tsStack(:,:,1,b,j);
        tmp = squeeze(tmp);
        if(strcmp(transform,'max'))
            data = max(tmp,[],3);
        end
        if(strcmp(transform,'mean'))
            data = mean(tmp,3);
        end
        
        imshow((data),[0 2^16]); colormap('jet');
        axis equal;    axis off;
        c = c+1;
    end
end
out =1;
warning('on');