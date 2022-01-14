%% Master File for Processing Data
clear; clc; close all;
stimulus = 'Odorants ("o"), cell culture ("c") or honeybee ("h") files?  ';
stimulus = input(stimulus,'s');
file_type = 'Spike sorted (ss) or RMS (r) data?  ';
file_type = input(file_type,'s');
input_files = 'Process multiple ("m") or individual ("s") file(s)? ';
input_files = input(input_files,'s');

filepath = '/Users/Xander/Documents/MATLAB/Neural_Recordings';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit variables below this line
time = [0 4];          % time relative to stimulus onset; must be within [-2:8]
new_bin_size = 50;
PSTH_smooth = 1;
PCA_smooth = 1;
PCA_line_mult = 500;      % line from origin every x msecs, MUST be set to multiple of new bin_size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit variables above this line
if input_files == 'm'
    filepath2 = 'Master_files';
    if stimulus == 'c'
        include_odors = 'Include odors? (y / n) ';
        include_odors = input(include_odors,'s');
        if file_type == 'r'
            filename = 'Cell_Culture_RMS_Master.mat';
            load([filepath '/' filepath2 '/' filename]);
            var_ext = '_RMS_data_filt';
            files = who([filepath '/' filepath2 '/' filename], ['*' var_ext]);
            if include_odors == 'n'
                files = files([1:3,5:6]);
            end
            bin_mult = new_bin_size/bin_size;
            bins_per_sec = 1000/(bin_size*bin_mult);
            for cycle_classes = 1:numel(files)
                Odorants{cycle_classes} = files{cycle_classes}(1:end-numel(var_ext));
                data_temp = squeeze(mean(eval(files{cycle_classes}),2));
                data(:,:,:,cycle_classes) = data_temp(:,:,(time(1)+2)*bins_per_sec+1:(time(2)+2)*bins_per_sec);
            end
        elseif file_type == 'ss'  
            filename = 'Cell_Culture_SS_Master.mat';
            load([filepath '/' filepath2 '/' filename]);
            var_ext = '_spikes_bin';
            files = who([filepath '/' filepath2 '/' filename], ['*' var_ext]);
            if include_odors == 'n'
                files = files([1:3,5:6]);
            end
            bin_mult = new_bin_size/bin_size;
            bins_per_sec = 1000/(bin_size*bin_mult);
            for cycle_classes = 1:numel(files)
                Odorants{cycle_classes} = files{cycle_classes}(1:end-numel(var_ext));
                data_temp = eval(files{cycle_classes});
                data(:,:,:,cycle_classes) = data_temp(:,:,(time(1)+2)*bins_per_sec+1:(time(2)+2)*bins_per_sec);
            end         
        end    
        if include_odors == 'y'
            Colors = [255,0,0 ; 0,255,0 ; 150,150,150 ; 255,255,0 ; 0,0,0 ; 0,0,255 ; 0,255,255]./255; 
        elseif include_odors == 'n'
            Colors = [255,0,0 ; 0,255,0 ; 150,150,150 ; 0,0,0 ; 0,0,255]./255; 
        end
    elseif stimulus == 'o'
        if file_type == 'r'
            filename = 'Odor_RMS_Master.mat';
            load([filepath '/' filepath2 '/' filename]);
            var_ext = '_RMS_data_filt';
            files = who([filepath '/' filepath2 '/' filename], ['*' var_ext]);
            bin_mult = new_bin_size/bin_size;
            bins_per_sec = 1000/(bin_size*bin_mult);
            for cycle_classes = 1:numel(files)
                Odorants{cycle_classes} = files{cycle_classes}(1:end-numel(var_ext));
                data_temp = squeeze(mean(eval(files{cycle_classes}),2));
                data(:,:,:,cycle_classes) = data_temp(:,:,(time(1)+2)*bins_per_sec+1:(time(2)+2)*bins_per_sec);
            end
        elseif file_type == 'ss'
            filename = 'Odor_SS_Master.mat';
            load([filepath '/' filepath2 '/' filename]);
            var_ext = '_spikes_bin';
            files = who([filepath '/' filepath2 '/' filename], ['*' var_ext]);
            bin_mult = new_bin_size/bin_size;
            bins_per_sec = 1000/(bin_size*bin_mult);
            for cycle_classes = 1:numel(files)
                Odorants{cycle_classes} = files{cycle_classes}(1:end-numel(var_ext));
                data_temp = eval(files{cycle_classes});
                data(:,:,:,cycle_classes) = data_temp(:,:,(time(1)+2)*bins_per_sec+1:(time(2)+2)*bins_per_sec);
            end
        end
        Colors = [180,0,0 ; 0,220,0 ; 102,178,255 ; 255,10,10 ; 192,192,192 ; 255,200,200 ; 0,225,225 ; 204,255,103 ; 153,0,153 ; 50,50,50]./255;       
%     elseif stimulus == 'h'
%         if file_type == 'r'
%         elseif file_type == 'ss'
%         end
    end
    
    mu = mean(data,[3 4]);
    data_norm = data - mu;
    odor_data_norm = squeeze(mean(data_norm,2));
    mu = mean(odor_data_norm,[2 3]);
    class_mu = squeeze(mean(odor_data_norm,2));   

    set(0, 'DefaultTextInterpreter', 'none');    % format underscores in title to avoid subscripting
    time = [stim_on_trim+time(1) stim_on_trim+time(2)];
    PCA_line_mult = PCA_line_mult/new_bin_size;
    titles = [num2str(time(1)-stim_on_trim) ' to ' num2str(time(2)-stim_on_trim) ' seconds after stimulus onset'];

    
    % pred_counts = zeros(numel(files),numel(files),size(data_norm,2));       % pre allocate array with zeros
    pred_counts_zeros = zeros(1,numel(files));                              % allocate array with zeros for padding
    for cycle_trials = 1:size(data,2)
        train_trials = find(1:size(data,2) ~= cycle_trials);           % training trials index
        train_data = squeeze(mean(data(:,train_trials,:,:),2));        % training data- trials average 

        test_trials = setdiff(1:size(data,2),train_trials);            % testing trials index
        test_data = squeeze(data(:,cycle_trials,:,:));                  % testing data
        test_data_temp = permute(test_data,[3 2 1]);

        for cycle_classes = 1:numel(files)
            train_data_temp(cycle_classes,:) = mean(train_data(:,:,cycle_classes),2);   % mean for training data with respect to each class
        end  
        for cycle_classes=1:numel(files)
            pred_data = squeeze(vecnorm(test_data_temp(cycle_classes,:,:)-permute(train_data_temp,[1 3 2]),2,3))';    % compute L2 norm for all time bins
            [~,pred_class] = min(pred_data,[],2);
            pred_counts_temp = accumarray(pred_class,1)';
            if numel(pred_counts_temp) ~= numel(files)
                pred_counts_temp(numel(pred_counts_zeros)) = 0;     % pad vector with 0 events if necessary
            end
            pred_counts(cycle_classes,:,cycle_trials) = pred_counts_temp;
        end 
    end
    pred_counts = sum(pred_counts,3)/(sum(pred_counts,'all')/numel(files));


    PCA_v2(odor_data_norm, Colors, Odorants, PCA_smooth, ['PCA: ' titles], time, stim_on_trim, bins_per_sec, PCA_line_mult)
    LDA_v2(odor_data_norm, class_mu, mu, Colors, Odorants, files, ['LDA: ' titles])
	Confusion(pred_counts,Odorants,['Confusion: ' titles])


%%%%%%%%%%%%%%%%%%% Single day analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif input_files == 's'
    Date = 'Specify date (MM_DD_YYYY): ';
    Date = input(Date,'s');
    Position = 'Specify position: ';
    Position = input(Position,'s');
    Tetrode = 'Specify tetrode: ';
    Tetrode = input(Tetrode,'s');
        
    if file_type == 'r'
        if stimulus == 'o'
            filepath2 = '/MAT_files/Odors';    
            Colors = [180,0,0 ; 0,220,0 ; 102,178,255 ; 255,10,10 ; 192,192,192 ; 255,200,200 ; 0,225,225 ; 204,255,103 ; 153,0,153 ; 50,50,50]./255;
%         elseif stimulus == 'cc'
%             filepath2 = 'MAT_files/Cell_Culture';
        end
        Odorant = 'Specify odorant name or all: ';
        Odorant = input(Odorant,'s');
        files = dir([filepath '/' filepath2 '/' Date '/Position_' Position]);
        if strcmp(Odorant, 'all')
            files = files(~ismember({files.name},{'.','..',['Position_' Position '.mat']}));
        else
            files = files(ismember({files.name},{[Odorant '.mat']})); 
        end
        for cycle_odors = 1:numel(files)
            Odorant = files(cycle_odors).name(1:end-4);
            if cycle_odors == 1
                load([filepath '/' filepath2 '/' Date '/Position_' Position '/' files(cycle_odors).name],[Odorant '_data_filt'],'stim_on','stim_off','sample_rate','total_time');
                load([filepath '/' filepath2 '/' Date '/Position_' Position '/Position_' num2str(Position) '.mat'],[Odorant '_RMS_data_filt'],'stim_on_trim','stim_off_trim','bin_size','trim_time');
                odor_RMS_data = eval([files(cycle_odors).name(1:end-4) '_RMS_data_filt']);

                trim_time = [stim_on_trim+time(1) stim_on_trim+time(2)];
                time = [stim_on+time(1) stim_on+time(2)];
                
                bin_mult = new_bin_size/bin_size;
                bins_per_sec = 1000/(bin_size*bin_mult);
            else
                load([filepath '/' filepath2 '/' Date '/Position_' Position '/' files(cycle_odors).name],[Odorant '_data_filt']);
                load([filepath '/' filepath2 '/' Date '/Position_' Position '/Position_' num2str(Position) '.mat'],[Odorant '_RMS_data_filt']);

            end
            odor_data = eval([Odorant '_data_filt']);        
            odor_data = odor_data(str2double(Tetrode),:,:,time(1)*sample_rate+1:(time(2)*sample_rate));            
            
            odor_RMS_data_base = mean(odor_RMS_data(str2double(Tetrode),:,:,(stim_on_trim-2)*bins_per_sec+1:stim_on_trim*bins_per_sec),4);
            odor_RMS_data_temp = odor_RMS_data(str2double(Tetrode),:,:,trim_time(1)*bins_per_sec+1:trim_time(2)*bins_per_sec);
            odor_RMS_data_temp = odor_RMS_data_temp-odor_RMS_data_base;

            exp_pars = struct('Date',Date,'Odorant',Odorant,'Position',Position,'stim_on',stim_on,'stim_on_trim',stim_on_trim,'stim_off',stim_off,'stim_off_trim',stim_off_trim,'sample_rate',sample_rate,'bins_per_sec',bins_per_sec,'total_time',total_time,'trim_time',trim_time);
            for cycle_chans = 1:size(odor_RMS_data_temp,2)            
                PSVT(exp_pars,odor_data(:,cycle_chans,:,:),str2double(Tetrode),cycle_chans,time)
                RMS_PSVT(exp_pars,odor_RMS_data_temp(:,cycle_chans,:,:),str2double(Tetrode),cycle_chans,trim_time)
    %             RMS_PSTH(exp_pars,file,cycle_tetrodes,cycle_channels,time)
    %             RMS_Raster(exp_pars,file,cycle_tetrodes,cycle_channels,time)
            end
        end

    elseif file_type == 'ss'
        Cell = 'Specify cell: ';
        Cell = input(Cell,'s');
        filepath2 = 'Spike_Sorted';
        if stimulus == 'o'
            filepath2 = [filepath2 '/Odors/Trim_Txt_files'];
            Colors = [180,0,0 ; 0,220,0 ; 102,178,255 ; 255,10,10 ; 192,192,192 ; 255,200,200 ; 0,225,225 ; 204,255,103 ; 153,0,153 ; 50,50,50]./255;
        elseif stimulus == 'c'
            filepath2 = [filepath2 '/Cell_Culture/Trim_Txt_files'];
                Colors = [180,0,0 ; 0,220,0 ; 102,178,255 ; 255,10,10 ; 192,192,192]./255;           
        end
        
        files = dir([filepath '/' filepath2 '/' Date '/Position_' Position '/Tetrode_' num2str(Tetrode) '/Cell' num2str(Cell) '/*.txt']);        
        load([filepath '/MAT_files/Odors/' Date '/Position_' Position '/Position_' num2str(Position) '.mat'],'stim_on','stim_on_trim','stim_off','stim_off_trim','trim_time','bin_size');
        set(0, 'DefaultTextInterpreter', 'none');    % format underscores in title to avoid subscripting

        trim_time = [stim_on_trim+time(1) stim_on_trim+time(2)];
%         time = [stim_on+time(1) stim_on+time(2)];
              
        Raster_title = ['Raster: ' num2str(trim_time(1)-stim_on_trim) ' to ' num2str(trim_time(2)-stim_on_trim) ' seconds after stimulus onset'];
        PSTH_title = ['PSTH: ' num2str(trim_time(1)-stim_on_trim) ' to ' num2str(trim_time(2)-stim_on_trim) ' seconds after stimulus onset'];
  
        raster_fig = figure('Position', [40 0 600 1000]);
        PSTH_fig = figure('Position', [40 0 600 1000]);
        for cycle_odors = 1:numel(files)
            odor = files(cycle_odors).name(1:end-4);
            spike_times = readmatrix([files(cycle_odors).folder '/' files(cycle_odors).name]);  % reads in raw spike times from txt files
            spikes = [];
            for cycle_trials = 1:size(spike_times,2)        % bins spike times and counts number of spikes per bin
                spikes_temp = nonzeros(spike_times(:,cycle_trials));
                spikes_temp = spikes_temp(spikes_temp >= trim_time(1) & spikes_temp <= trim_time(2));
                spikes_bin(cycle_trials,:) = histcounts(spikes_temp,'BinLimits',[trim_time(1) trim_time(2)],'BinWidth',bin_size/1000);

                if cycle_trials ~= 1 && size(spikes,1) < size(spikes_temp,1)
                    spikes = padarray(spikes,size(spikes_temp,1)-size(spikes,1),'post');
                elseif cycle_trials ~= 1 && size(spikes,1) > size(spikes_temp,1)
                    spikes_temp = padarray(spikes_temp,size(spikes,1)-size(spikes_temp,1),'post');
                end
                spikes = [spikes spikes_temp];
            end
            figure(raster_fig)
            Raster(cycle_odors) = subplot(numel(files),1,cycle_odors);   
            plot_Raster_Master(spikes,time,stim_on_trim,stim_off_trim,odor,Cell,cycle_odors,numel(files),Colors,date,Raster_title);

            figure(PSTH_fig)
            PSTH(cycle_odors) = subplot(numel(files),1,cycle_odors); 
            [PSTH_max(cycle_odors), PSTH_patch(cycle_odors)] = plot_PSTH_Master(spikes_bin,time,stim_on_trim,stim_off_trim,odor,Cell,cycle_odors,numel(files),Colors,date,PSTH_smooth,1000/new_bin_size,PSTH_title);  

            var_name = [odor '_spike_times'];
            eval([var_name '= spikes;']);
        end
        linkaxes(PSTH,'y');
        PSTH(1,1).YLim = [0 max(PSTH_max)];
        linkprop(PSTH_patch, 'YData');
        PSTH_patch(1,1).YData = [0 max(PSTH_max) max(PSTH_max) 0];
    end
end