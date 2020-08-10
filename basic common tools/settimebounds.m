function [ outtime, outvector ] = settimebounds( intime, invector, timelimit1, timelimit2  )

% for things like eggs or transitions that may not be sorted in time
[a b] = sort(intime, 'ascend');
invector = invector(b);
intime = intime(b);

[a b] = find(intime >= timelimit1,1,'first');
[a1 b1] = find(intime <= timelimit2,1,'last');

outtime = intime(a:a1);
outvector = invector(a:a1);

end

