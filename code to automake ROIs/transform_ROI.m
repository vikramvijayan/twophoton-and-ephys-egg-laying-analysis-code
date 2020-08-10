function [ROI2] = transform_ROI(tsStack1, tsStack2, ROI1)
warning('off');
q2 = mean(mean(tsStack2,5),4);
q1 = mean(mean(tsStack1,5),4);
[optimizer, metric] = imregconfig('multimodal');

optimizer.MaximumIterations = 5000;
optimizer.Epsilon = 1e-10;

tform = imregtform(q1,q2,'affine',optimizer, metric);

ROI2 = imwarp(ROI1,tform,'OutputView',imref2d(size(q2)));
warning('on');
