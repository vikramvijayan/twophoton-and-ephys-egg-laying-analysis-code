function [no_puff] = find_no_puff(puff, buffer_in_frames )
[a b] = find(puff > .65);

q = zeros(length(puff),1);
q(a) = 1;
q2 = smooth(q,buffer_in_frames*2+1);

[no_puff b] = find(q2 == 0);
