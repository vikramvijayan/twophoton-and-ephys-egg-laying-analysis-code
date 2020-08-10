function [tsStack] = load2pstack(fn, c, z, t, newROI, TakeMax)

if (z ~= 1)
        tsStack = TIFFStack(fn,[],[c z t],1);

    %tsStack = TIFFStack(fn,[]);
    tsStack = permute(tsStack,[1 2 3 5 4]);
else
    tsStack = TIFFStack(fn,[],[c t z],1);
end

size(tsStack)
if(newROI)
    ROI_select(tsStack,fn, TakeMax);
end
end