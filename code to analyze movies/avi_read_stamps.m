function [ stamps] = avi_read_stamps( filename )

v = VideoReader(filename);

% goes through each frame of the avi and calculates the three different
% parts of the timestamp (cyclecount, secondcount, and cycle offset)

len = 0;
while(hasFrame(v))
    len = len+1;
    q(:,:,:) = readFrame(v);
    m(:) = q(1,1:4,1);
    mm(:) = convert_d(double(m(:)));
    cyclecount(len) = bin2dec(mm(8:20));
    secondcount(len) = bin2dec(mm(1:7));
    cycleoff(len) = bin2dec(mm(21:32));
end

% assembles the three different parts of the timestamp appropriately
stamps = secondcount + (cyclecount+cycleoff./3072)./8000;
stampsnotwrapped = stamps;


% very inefficiently corrects for wrap around of the timestamp. That is,
% after 128 seconds the timestamp resets and we obviously don't want this.
add = 0;
addnum = 0;
for i = 1:1:(length(stamps)-1)
    if(stamps(i+1)-stamps(i) < 0)
        addnum = addnum+128;
    end
    add = [add addnum];
end

stamps = stamps+add;

end

