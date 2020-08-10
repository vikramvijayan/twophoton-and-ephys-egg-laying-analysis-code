function recording = reload_tsStack_and_use_NoRM_registration(recording)
for i = 0:1:(length(recording.tseries)-1)
    
    
    
    name = [recording.abfname(1:end-4) '-00' num2str(i) '.tif'];
    % this part is commented but need to uncomment torun
    tempfile = TIFFStack(name,[],[1,recording.tseries(i+1).zslices, recording.tseries(i+1).Total2pFrames],1);
    tempfile = permute(tempfile,[1 2 3 5 4]);
    tempfile = mean(tempfile,5);
    tempfile = squeeze(tempfile);
    
    M_final = new_registration_script_singleZ(tempfile,name);


    t = [];
    for j=1:1:recording.tseries(i+1).zslices
        t = cat(3,t,M_final);
    end
    res = saveastiff(t, [name(1:(end-4)) '_new_reg_script.tif']);

    recording.tseries(i+1).tsStack = t;
    
    % for some reason the tiff written by the above function cant be opened
    % with TiffStack. need to resave them in ImageJ and then run the below.
    
    %     file_name = [name(1:(end-4)) '_new_reg_script.tif'];
    %
    %     tempfile = TIFFStack(file_name,[],[1,recording.tseries(i+1).zslices, recording.tseries(i+1).Total2pFrames],1);
    %     tempfile = permute( recording.tseries(i+1).tsStack,[1 2 3 5 4]);
    %         recording.tseries(i+1).tsStack = tempfile;
    
    
end