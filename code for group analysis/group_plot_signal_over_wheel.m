% need to run this on a test group

% currently edititing thos 03/19/2020
% For each data point need to output the last time that the fly was on all
% other substrates (as well as the time since the trasition) -- which is
% just the min of above


function [wheel_structure] = group_plot_signal_over_wheel(recordings_list, ROI_num)

wheel_structure.median_wheel = [];
wheel_structure.mean_wheel = [];
wheel_structure.position_wheel = [];
wheel_structure.s0vals = [];
wheel_structure.s0vals_interp = [];

wheel_structure.s0vals_since0trans = [];
wheel_structure.s0vals_since200 = [];
wheel_structure.s0vals_since500 = [];
wheel_structure.s200vals = [];
wheel_structure.s200vals_interp = [];

wheel_structure.s200vals_since200trans = [];
wheel_structure.s200vals_since0 = [];
wheel_structure.s200vals_since500 = [];

wheel_structure.s500vals = [];
wheel_structure.s500vals_interp = [];
wheel_structure.s500vals_since500trans = [];
wheel_structure.s500vals_since0 = [];
wheel_structure.s500vals_since200 = [];


for rec_index = 1:1:length(recordings_list)
    
    
    
    %for rec_index = 10
    
    % it is often faster to use the stripped files without the associated
    % images
    modifiedStr = strrep([char(recordings_list(rec_index))], '.mat', '_stripped.mat');
    modifiedStr2 = strrep([char(recordings_list(rec_index))], '.mat', '_stripped_new_dlc_track.mat');
    
    if(exist([modifiedStr2]))
        modifiedStr = modifiedStr2;
    end
    recording = loadsinglerecording(modifiedStr);
    [recording] = processes_more_DLC_variables(recording);
    
    
    
    
    
    signal_store_interp = [];
    signal_store = [];
    signal_time_store = [];
    filtered_wheel_store = [];
    sucrose_store = [];
    sucrose_store_interp = [];
    
    for i = 1:1:length(recording.tseries)
        
        signal      = recording.tseries(i).df_over_f(ROI_num(rec_index),:)./nanmean(recording.tseries(i).df_over_f(ROI_num(rec_index),:));
        signal_time = recording.tseries(i).Time_s;
        interp_filt_wheel = interp1(recording.movie1.time_stamps, recording.movie1.filtered_wheel, signal_time,'previous');
        interp_sucrose = interp1(recording.movie1.time_stamps, recording.movie1.sucrose, signal_time,'previous');
        
        signal_interp = interp1(signal_time, signal, signal_time(1):.1:signal_time(end),'previous');
        sucrose_interp = interp1(recording.movie1.time_stamps, recording.movie1.sucrose,  signal_time(1):.1:signal_time(end),'previous');
        
        
        sucrose_store_interp = [sucrose_store_interp; sucrose_interp'];
        signal_store_interp  = [signal_store_interp; signal_interp'];
        signal_store = [signal_store; signal'];
        signal_time_store = [signal_time_store; signal_time];
        filtered_wheel_store = [filtered_wheel_store; interp_filt_wheel];
        sucrose_store = [sucrose_store, interp_sucrose'];
    end
    
    
    [a b] = find(sucrose_store_interp' ==200);
    wheel_structure(rec_index).s200vals_interp = signal_store_interp(b);
    [a b] = find(sucrose_store_interp' ==500);
    wheel_structure(rec_index).s500vals_interp = signal_store_interp(b);
    [a b] = find(sucrose_store_interp' ==0);
    wheel_structure(rec_index).s0vals_interp = signal_store_interp(b);
    
    
    
    [a b] = find(sucrose_store' ==200);
    %find last time sucrose was a particular value
    
    for i = 1:1:length(a)
        [a1 b1] = (min(abs(signal_time_store(a(i))-recording.movie1.time_stamps)));
        [last0 last0p] = find(recording.movie1.sucrose(1:b1) == 0,1,'last');
        [last500 last500p] = find(recording.movie1.sucrose(1:b1) == 500,1,'last');
        [last200t last200tp] = find(diff(recording.movie1.sucrose(1:b1)) ~= 0,1,'last');
        if(isempty(last0p))
            last0p = inf;
        end
        if(isempty(last500p))
            last500p = inf;
        end
        if(isempty(last200tp))
            last200tp = inf;
        end
        wheel_structure(rec_index).s200vals_since200trans(i) = b1-last200tp;
        wheel_structure(rec_index).s200vals_since0(i) = b1-last0p;
        wheel_structure(rec_index).s200vals_since500(i) = b1-last500p;
    end
    
    wheel_structure(rec_index).s200vals = signal_store(a);
    wheel_structure(rec_index).s200mean = nanmean(signal_store(a));
    wheel_structure(rec_index).s200time = signal_time_store(a);
    


        
        
    
    [a b] = find(sucrose_store' ==500);
    
    for i = 1:1:length(a)
        [a1 b1] = (min(abs(signal_time_store(a(i))-recording.movie1.time_stamps)));
        [last0 last0p] = find(recording.movie1.sucrose(1:b1) == 0,1,'last');
        [last200 last200p] = find(recording.movie1.sucrose(1:b1) == 200,1,'last');
        [last500t last500tp] = find(diff(recording.movie1.sucrose(1:b1)) ~= 0,1,'last');
        if(isempty(last0p))
            last0p = inf;
        end
        if(isempty(last500tp))
            last500tp = inf;
        end
        if(isempty(last200p))
            last200p = inf;
        end
        wheel_structure(rec_index).s500vals_since500trans(i) = b1-last500tp;
        wheel_structure(rec_index).s500vals_since0(i) = b1-last0p;
        wheel_structure(rec_index).s500vals_since200(i) = b1-last200p;
    end
    
    wheel_structure(rec_index).s500vals = signal_store(a);
    wheel_structure(rec_index).s500mean = nanmean(signal_store(a));
    wheel_structure(rec_index).s500time = signal_time_store(a);
    
    [a b] = find(sucrose_store' == 0);
    
    for i = 1:1:length(a)
        [a1 b1] = (min(abs(signal_time_store(a(i))-recording.movie1.time_stamps)));
        [last500 last500p] = find(recording.movie1.sucrose(1:b1) == 500,1,'last');
        [last200 last200p] = find(recording.movie1.sucrose(1:b1) == 200,1,'last');
        [last0t last0tp] = find(diff(recording.movie1.sucrose(1:b1)) ~= 0,1,'last');
        if(isempty(last0tp))
            last0tp = inf;
        end
        if(isempty(last500p))
            last500p = inf;
        end
        if(isempty(last200p))
            last200p = inf;
        end
        wheel_structure(rec_index).s0vals_since0trans(i) = b1-last0tp;
        wheel_structure(rec_index).s0vals_since500(i) = b1-last500p;
        wheel_structure(rec_index).s0vals_since200(i) = b1-last200p;
    end
    
    
    wheel_structure(rec_index).s0vals = signal_store(a);
    wheel_structure(rec_index).s0mean = nanmean(signal_store(a));
    wheel_structure(rec_index).s0time = signal_time_store(a);
    
    
    cnt = 1;
    m = [];
    medin = [];
    
    tmp_filtered_wheel_store = [filtered_wheel_store; filtered_wheel_store+2*pi; filtered_wheel_store+4*pi; filtered_wheel_store+6*pi];
    tmp_signal_store = [signal_store; signal_store; signal_store; signal_store];
    
    for iterat = 0:pi/50:(6*pi)
        [a b] = find(tmp_filtered_wheel_store >= iterat & tmp_filtered_wheel_store < (iterat+pi/4));
        m(cnt) = nanmean(tmp_signal_store(a));
        medin(cnt) = nanmedian(tmp_signal_store(a));
        
        cnt = cnt+1;
    end
    
    [a b] = sort(mod(pi/8+(0:pi/50:6*pi),2*pi),'ascend');
    
    wheel_structure(rec_index) .median_wheel= medin(b);
    wheel_structure(rec_index).mean_wheel = m(b);
    wheel_structure(rec_index).position_wheel  = a;
    
    
end


% figure; hold on;
% for i =1:1:length(wheel_structure)
%     plot(wheel_structure(i).position_wheel(:),wheel_structure(i).mean_wheel(:),'color',[.75,.75,.75]);
%     m(i,:) = wheel_structure(i).mean_wheel(:);
% end
% plot(wheel_structure(1).position_wheel(:),mean(m),'k','linewidth',3)
% set(gca,'xlim',[0 2*pi]);
% snapnow;
% 
% figure; hold on; scatter(rand(length([wheel_structure.s0mean]),1)',[wheel_structure.s0mean], 10,[107 174 214]./255,'filled');
% scatter(3+rand(length([wheel_structure.s200mean]),1)',[wheel_structure.s200mean], 10,[33 113 181]./255,'filled');
% scatter(6+rand(length([wheel_structure.s500mean]),1)',[wheel_structure.s500mean], 10,[8 48 107]./255,'filled');
% 
% 
% scatter(.5,mean([wheel_structure.s0mean]), 30,[107 174 214]./255,'filled');
% scatter(3.5,mean([wheel_structure.s200mean]), 30,[33 113 181]./255,'filled');
% scatter(6.5,mean([wheel_structure.s500mean]), 30,[8 48 107]./255,'filled');
% errorbar(.5,mean([wheel_structure.s0mean]),std([wheel_structure.s0mean])./sqrt(length([wheel_structure.s0mean])),'Color',[107 174 214]./255);
% errorbar(3.5,mean([wheel_structure.s200mean]),std([wheel_structure.s200mean])./sqrt(length([wheel_structure.s200mean])),'Color',[33 113 181]./255);
% errorbar(6.5,mean([wheel_structure.s500mean]),std([wheel_structure.s500mean])./sqrt(length([wheel_structure.s500mean])),'Color',[8 48 107]./255);
% 
% text(1,mean([wheel_structure.s0mean]),num2str(mean([wheel_structure.s0mean])));
% text(4.5,mean([wheel_structure.s200mean]),num2str(mean([wheel_structure.s200mean])));
% text(7.5,mean([wheel_structure.s500mean]),num2str(mean([wheel_structure.s500mean])));
% 
% set(gca,'xlim',[0 9]); snapnow;

end

