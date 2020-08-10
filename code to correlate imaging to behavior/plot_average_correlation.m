function [] = plot_average_correlation( out_corr, out_corr_shift_ofvec1, shift_vec1_forward_sec,shift_vec1_backward_sec )

dt = .01;
figure; hold on;
plot(out_corr_shift_ofvec1,out_corr,'Color',[.5,.5,.5]);
plot(out_corr_shift_ofvec1,mean(out_corr),'b');
set(gca,'xtick',[-1.*(1/dt)*shift_vec1_backward_sec:2/dt:shift_vec1_forward_sec*(1/dt)]);
set(gca,'xticklabel',-1.*shift_vec1_backward_sec:2:shift_vec1_forward_sec);
line([0 0], [-1, 1],'color','k');
xlabel('Time (sec)');


end

