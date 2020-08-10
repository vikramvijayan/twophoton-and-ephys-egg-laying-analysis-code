function saveallfigs_2p(dir)

h = get(0,'children');
for i=1:length(h)
    %set(h(i),'Visible','on');
    figure(h(i));
    %   title = get(gca,'Title');
    %   tit_edit = regexprep(title.String,'/',' per ');
    %   tit_edit = regexprep(tit_edit,'>',' gt ');
    
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, [dir 'figure' num2str(i) '.fig']);
    saveas(gcf, [dir 'figure' num2str(i) '.png']);
    
    %saveas(h(i), [dir  tit_edit '.png']);
    close(h(i));
end