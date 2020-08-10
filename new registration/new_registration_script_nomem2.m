% based on NoRM demo in CalmAN
function M_final = new_registration_script_nomem2(name,path)

name_small  = name(1:(end-4));
Path4Fiji='C:\Fiji.app\ImageJ-win64.exe';
Path4postreg_Macro='C:\Users\vikra\Documents\Postdoc\common_commands_scripts_macros\postreg_mergeZC_to_hyperstack.ijm';

% 
% disp 'loading tsstack'
% 
% Y = TIFFStack([name]);
% disp 'done loading tsstack'
% 
% Y = permute(Y,[1 2 4 3]);
% 
% size(Y)
% %at this point you have a 4d array x,y,z,t
% 
% %% now try non-rigid motion correction (also in parallel)
% disp 'starting registation'
% 
% % It is critical that the grid size (first two numbers) be less than two times the dimension of the image
% 
% options=NoRMCorreSetParms('d1             ',size(Y,1),...          %% dataset info % number of rows
%     'd2             ',size(Y,2),...           % number of cols
%     'd3             ',size(Y,3),...           % number of planes (for 3d imaging, default: 1) 
%     'grid_size      ',[24,24,2],...             %% patches % size of non-overlapping regions (default: [d1,d2,d3])
%     'overlap_pre    ',[32,32,2],...              % size of overlapping region (default: [32,32,16])
%     'min_patch_size ',[32,32,2],...              % minimum size of patch (default: [32,32,16])
%     'min_diff       ',[16,16,2],...               % minimum difference between patches (default: [16,16,5])
%     'us_fac         ',50,...                      % upsampling factor for subpixel registration (default: 20)
%     'mot_uf         ',[4,4,1],...                 % degree of patches upsampling (default: [4,4,1])
%     'max_dev        ',[3,3,1],...                 % maximum deviation of patch shift from rigid shift (default: [3,3,1])
%     'overlap_post   ',[32,32,24],...              % size of overlapping region after upsampling (default: [32,32,16])
%     'max_shift      ',[25,25,1],...               % maximum rigid shift in each direction (default: [15,15,5])
%     'phase_flag     ',false,...                   % flag for using phase correlation (default: false)
%     'shifts_method  ','fft',...                   % method to apply shifts ('FFT','cubic','linear')
%     'upd_template   ',true,...                   %% template updating % flag for online template updating (default: true)
%     'init_batch     ',200,...                     % length of initial batch (default: 100)
%     'bin_width      ',200,...                     % width of each bin (default: 200)
%     'buffer_width   ',50,...                      % number of local means to keep in memory (default: 50)
%     'method         ',{'median';'mean'},...       % method for averaging the template (default: {'median';'mean})
%     'iter           ',1,...                       % number of data passes (default: 1)
%     'boundary       ','copy',...                  % method of boundary treatment 'NaN','copy','zero','template' (default: 'copy')
%     'add_value      ',0,...                      %% misc % add dc value to data (default: 0)
%     'use_parallel   ',false,...                   % for each frame, update patches in parallel (default: false)
%     'memmap         ',false,...                   % flag for saving memory mapped motion corrected file (default: false)
%     'mem_filename   ','motion_corrected.mat',...  % name for memory mapped file (default: 'motion_corrected.mat')
%     'mem_batch_size ',5000,...                    % batch size during memory mapping for speed (default: 5000)
%     'print_msg      ',true,...                    % flag for printing progress to command line (default: true)
%     'plot_flag      ',false,...                  %% plotting % flag for plotting results in real time (default: false)
%     'make_avi       ',false,...                   % flag for making movie (default: false)
%     'name           ','motion_corrected.avi',...  % name for movie (default: 'motion_corrected.avi')
%     'fr             ',30,...                      % frame rate for movie (default: 30)
%     'output_type    ','mat',...                  %% output type % 'mat' (load in memory), 'memmap', 'tiff', 'hdf5', 'bin' (default:mat)
%     'h5_groupname   ','mov',...                   % name for hdf5 dataset (default: 'mov')
%     'h5_filename    ','motion_corrected.h5',...   % name for hdf5 saved file (default: 'motion_corrected.h5')
%     'tiff_filename  ','motion_corrected.tif',...  % name for saved tiff stack (default: 'motion_corrected.tif')
%     'output_filename','',...                      % name for saved file will be used if `h5_,tiff_filename` are not specified
%     'use_windowing  ',false,...                  %% use windowing % flag for windowing data before fft (default: false)
%     'window_length  ',0.5,...                     % length of window on each side of the signal as a fraction of signal length | total length = length(signal)(1 + 2*window_length). (default: 0.5)
%     'bitsize        ',2,...                      %% bitsize for reading .raw files % (default: 2 (uint16). other choices 1 (uint8), 4 (single), 8 (double))
%     'correct_bidir  ',false,...                  %% offset from bidirectional sampling % check for offset due to bidirectional scanning (default: true)
%     'nFrames        ',50,...                      % number of frames to average (default: 50)
%     'bidir_us       ',10,...                      % upsampling factor for bidirectional sampling (default: 10)
%     'col_shift      ',[]);                        % known bi-directional offset provided by the user (default: [])
% 
% % % 
% tic; [M2,shifts2,template2,options_nonrigid] = normcorre_batch(single(Y),options); toc
% %tic; [M2,shifts2,template2,options_nonrigid] = normcorre(single(Y),options); toc
% 
% disp 'done registration'
% 
% 
% disp 'saving z*c .tif file(s)';

tic;
cd(path);
% 
% M2=uint16(permute(M2,[1,2,4,3])); % yxtzc
% for k=1:size(M2,4) % z
%     saveastiff(M2(:,:,:,k),[path, '\', num2str(k),'_1','.tif']); % FileName=z_c.tif
% end
toc

 disp 'merging z*c .tif file(s) to a xyczt-hyperstack';

    % find small FileNames.tif in the folder
    FileName=dir(['*.tif']);
    FileName_size=zeros(length(FileName),1);
    for j=1:length(FileName)
    FileName_size(j,1)=length(FileName(j,1).name);
    end
    FileName(FileName_size>8,:)=[];
    FileName={FileName.name}';
    PathFileName=cell(length(FileName),1);
    FileName2Fiji_c1='title=1 open';
    FileName2Fiji_c2='title=2 open';
    for j=1:length(FileName)
        if ispc; PathFileName{j,1}=[path,'\',FileName{j,1}];
        elseif ismac; PathFileName{j,1}=[path,'/',FileName{j,1}]; end
        if j==1
           PathFileName2Fiji=PathFileName{j,1};
        else PathFileName2Fiji=[PathFileName2Fiji,' ',PathFileName{j,1}];
        end
        if strcmp(FileName{j,1}(end-4),'1')
           FileName2Fiji_c1=[FileName2Fiji_c1,' image',num2str(FileName{j,1}(1:end-6)),'=',FileName{j,1}];
        else FileName2Fiji_c2=[FileName2Fiji_c2,' image',num2str(FileName{j,1}(1:end-6)),'=',FileName{j,1}];
        end
    end
    PathFileName2Fiji=[PathFileName2Fiji,' ',path,path(end-6:end),'_reg.tif'];
    if ispc; PathFileName2Fiji=strrep(PathFileName2Fiji,'\','\\'); end
    if ispc; system([Path4Fiji,' -batch ',Path4postreg_Macro,' "',PathFileName2Fiji,'/',FileName2Fiji_c1,'/',FileName2Fiji_c2,'"']);
    elseif ismac; system([Path4Fiji,' -batch ',Path4postreg_Macro,' "',PathFileName2Fiji,'\\',FileName2Fiji_c1,'\\',FileName2Fiji_c2,'"']); end
%    delete(FileName{:});
%    FileName=[Path(1).folder(end-5:end),'.tif'];
%    delete(FileName);

