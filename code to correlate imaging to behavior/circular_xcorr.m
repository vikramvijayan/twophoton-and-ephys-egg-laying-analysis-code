function [out_corr,out_corr_shift_ofvec1] = circular_xcorr(vect1, vect2, time_vec1, time_vec2, time_base_to_return)

shift_vec1_forward_sec = time_base_to_return(end);
shift_vec1_backward_sec = -1.*time_base_to_return(1);

dt = mean(diff(time_base_to_return));

vect1_interp = interp1(time_vec1, vect1, [max(time_vec1(1),time_vec2(1)):dt:min(time_vec1(end),time_vec2(end))],'previous'); 
vect2_interp = interp1(time_vec2, vect2, [max(time_vec1(1),time_vec2(1)):dt:min(time_vec1(end),time_vec2(end))],'previous'); 

% Vector 1 is shifted forward and backward while circularly correlating.
% Input must be in radians
% That is if the max corr is at shift greater than 0, then vec1 comes
% after vec2
out_corr1 = [];
out_corr1_shift_ofvec1 = [];

for i = 1:1:(shift_vec1_forward_sec*(1/dt))
    out_corr1(i) = circ_corrcc(vect1_interp(i:end),vect2_interp(1:(end+1-i)));
    out_corr1_shift_ofvec1(i) = i-1;
end

out_corr2 = [];
out_corr2_shift_ofvec1 = [];

for i = 1:1:(shift_vec1_backward_sec*(1/dt))
    out_corr2(i) = circ_corrcc(vect2_interp(i:end),vect1_interp(1:(end+1-i)));
    out_corr2_shift_ofvec1(i) = -1*(i-1);
    
end

out_corr = [fliplr(out_corr2), out_corr1(2:end)];
out_corr_shift_ofvec1 = [fliplr(out_corr2_shift_ofvec1), out_corr1_shift_ofvec1(2:end)];

our_corr_shift_ofvec1 = [-1.*shift_vec1_backward_sec:dt:shift_vec1_forward_sec];

% figure; hold on;
% plot(out_corr_shift_ofvec1,out_corr,'k');
% set(gca,'xtick',[-1.*(1/dt)*shift_vec1_backward_sec:1/dt:shift_vec1_forward_sec*(1/dt)]);
% set(gca,'xticklabel',-1.*shift_vec1_backward_sec:1:shift_vec1_forward_sec);
% line([0 0], [-1, 1],'color','k');
% xlabel('Time (sec)');