function recording = reload_tsStack(recording)
for i = 0:1:(length(recording.tseries)-1)
   % recording.tseries(i+1).tsStack = TIFFStack([recording.abfname(1:end-4) '-00' num2str(i) '_reg_interleaved.tif'],[],[1,recording.tseries(i+1).zslices, recording.tseries(i+1).Total2pFrames],1);
   %   recording.tseries(i+1).tsStack = TIFFStack([recording.abfname(1:end-4) '-00' num2str(i) '_reg.tif'],[],[1,recording.tseries(i+1).zslices, recording.tseries(i+1).Total2pFrames],1);
  % recording.tseries(i+1).tsStack = TIFFStack([recording.abfname(1:end-4) '-00' num2str(i) '.tif'],[],[1,recording.tseries(i+1).zslices, recording.tseries(i+1).Total2pFrames],1);
   recording.tseries(i+1).tsStack = TIFFStack([recording.abfname(1:end-4) '-00' num2str(i) '_new_reg_script2.tif'],[],[1,recording.tseries(i+1).zslices, recording.tseries(i+1).Total2pFrames],1);

   recording.tseries(i+1).tsStack = permute( recording.tseries(i+1).tsStack,[1 2 3 5 4]);

end