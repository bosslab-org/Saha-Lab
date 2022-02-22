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
PCA_line_mult = 250;      % line from origin every x msecs, MUST be set to multiple of new bin_size
lambda = 0;               % Vanilla LDA when lambda is zero; spherical LDA when lambda is one; must be within [0:1]
LDA_dims = 106;             % Reduces number of LDA-transformed data dimensions to specified value- only applies to confusion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit variables above this line
if input_files == 'm'
    filepath2 = 'Master_files';
    if stimulus == 'c'
        include_odors = 'Include odors? (y / n) ';
        include_odors = input(include_odors,'s');
        daywise_comp = 'Daywise comparison? (y /n) ';
        daywise_comp = input(daywise_comp,'s');
        if file_type == 'r'
            if daywise_comp == 'n'
                filename = 'Cell_Culture_RMS_Master.mat';
            elseif daywise_comp == 'y'
                select_day = 'Select day (1, 2, 3, 4): ';
                select_day = input(select_day,'s');
                switch select_day
                    case '1', filename = 'Cell_Culture_RMS_Day1.mat';
                    case '2', filename = 'Cell_Culture_RMS_Day2.mat';
                    case '3', filename = 'Cell_Culture_RMS_Day3.mat';
                    case '4', filename = 'Cell_Culture_RMS_Day4.mat';
                end
            end
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
            if daywise_comp == 'n'
                filename = 'Cell_Culture_SS_Master.mat';
            elseif daywise_comp == 'y'
                select_day = 'Select day (1, 2, 3, 4): ';
                select_day = input(select_day,'s');
                switch select_day
                    case '1', filename = 'Cell_Culture_SS_Day1.mat';
                    case '2', filename = 'Cell_Culture_SS_Day2.mat';
                    case '3', filename = 'Cell_Culture_SS_Day3.mat';
                    case '4', filename = 'Cell_Culture_SS_Day4.mat';
                end
            end
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
%             filename = 'Odor_RMS_Master.mat';
%             load([filepath '/' filepath2 '/' filename]);
%             var_ext = '_RMS_data_filt';
%             files = who([filepath '/' filepath2 '/' filename], ['*' var_ext]);
%             bin_mult = new_bin_size/bin_size;
%             bins_per_sec = 1000/(bin_size*bin_mult);
%             for cycle_classes = 1:numel(files)
%                 Odorants{cycle_classes} = files{cycle_classes}(1:end-numel(var_ext));
%                 data_temp = squeeze(mean(eval(files{cycle_classes}),2));
%                 data(:,:,:,cycle_classes) = data_temp(:,:,(time(1)+2)*bins_per_sec+1:(time(2)+2)*bins_per_sec);
%             end
%         elseif file_type == 'ss'
%             filename = 'Odor_SS_Master.mat';
%             load([filepath '/' filepath2 '/' filename]);
%             var_ext = '_spikes_bin';
%             files = who([filepath '/' filepath2 '/' filename], ['*' var_ext]);
%             bin_mult = new_bin_size/bin_size;
%             bins_per_sec = 1000/(bin_size*bin_mult);
%             for cycle_classes = 1:numel(files)
%                 Odorants{cycle_classes} = files{cycle_classes}(1:end-numel(var_ext));
%                 data_temp = eval(files{cycle_classes});
%                 data(:,:,:,cycle_classes) = data_temp(:,:,(time(1)+2)*bins_per_sec+1:(time(2)+2)*bins_per_sec);
%             end
%         end
    end
    
    mu = mean(data,[3 4]);
    data_norm = data - mu;                      % need to normalize training and testing data separately for confusion matrix
    odor_data_norm = squeeze(mean(data_norm,2));  

    set(0, 'DefaultTextInterpreter', 'none');    % format underscores in title to avoid subscripting
    time = [stim_on_trim+time(1) stim_on_trim+time(2)];
    PCA_line_mult = PCA_line_mult/new_bin_size;
    titles = [num2str(time(1)-stim_on_trim) ' to ' num2str(time(2)-stim_on_trim) ' seconds after stimulus onset'];
    
    pred_counts_zeros = zeros(1,numel(files));                         % preallocate array with zeros for padding
    for cycle_trials = 1:size(data,2)
        train_trials = find(1:size(data,2) ~= cycle_trials);           % training trials index
        train_data = data(:,train_trials,:,:);                         % training data- trials average         
        HD_train_data_mean = permute((mean(train_data,[2 3])),[4 1 2 3]);        % Training data class averages in high dimensional feature space
        [LDA_train_data,LDA_train_data_mean] = LDA_transform(train_data, files, lambda, LDA_dims);   % Training data class averages in LDA feature space            
        LDA_train_data_mean = squeeze(LDA_train_data_mean);

        test_trials = setdiff(1:size(data,2),train_trials);            % testing trials index
        test_data = data(:,cycle_trials,:,:);                          % testing data        
        HD_test_data = permute(test_data,[4 3 1 2]);
        [LDA_test_data,~] = LDA_transform(test_data, files, lambda, LDA_dims);

        for cycle_classes=1:numel(files)
            [HD_trial_pred_counts_temp(cycle_classes,1), HD_bin_pred_counts(cycle_classes,:,cycle_trials)] = Test_data(HD_test_data,HD_train_data_mean, files, pred_counts_zeros, cycle_classes);                  % Hi-D feature space testing
            [LDA_trial_pred_counts_temp(cycle_classes,1), LDA_bin_pred_counts(cycle_classes,:,cycle_trials)] = Test_data(LDA_test_data,LDA_train_data_mean, files, pred_counts_zeros, cycle_classes);           % LDA feature space testing  
        end
        HD_trial_pred_counts(:,:,cycle_trials) = accumarray([HD_trial_pred_counts_temp,(1:numel(files))'],ones(1,numel(files)));
        LDA_trial_pred_counts(:,:,cycle_trials) = accumarray([LDA_trial_pred_counts_temp,(1:numel(files))'],ones(1,numel(files)));
    end
    HD_bin_pred_counts = sum(HD_bin_pred_counts,3)/(sum(HD_bin_pred_counts,'all')/numel(files));
    HD_trial_pred_counts = sum(HD_trial_pred_counts,3)/(sum(HD_trial_pred_counts,'all')/numel(files));
    LDA_bin_pred_counts = sum(LDA_bin_pred_counts,3)/(sum(LDA_bin_pred_counts,'all')/numel(files));
    LDA_trial_pred_counts = sum(LDA_trial_pred_counts,3)/(sum(LDA_trial_pred_counts,'all')/numel(files));

    PCA_fig = PCA_v2(odor_data_norm, Colors, Odorants, PCA_smooth, ['PCA (' file_type '): ' titles], time, stim_on_trim, bins_per_sec, PCA_line_mult);
    LDA_fig = LDA_v2(odor_data_norm, Colors, Odorants, files, lambda, ['LDA (' file_type '): ' titles]);
	Conf_fig_HD_bin = Confusion(HD_bin_pred_counts,Odorants,['HD Binwise Confusion (' file_type '): ' titles]);
    Conf_fig_HD_trial = Confusion(HD_trial_pred_counts,Odorants,['HD Trialwise Confusion (' file_type '): ' titles]);
	Conf_fig_LDA_bin = Confusion(LDA_bin_pred_counts,Odorants,['LDA Binwise Confusion (' file_type '): ' titles]);
	Conf_fig_LDA_trial = Confusion(LDA_trial_pred_counts,Odorants,['LDA Trialwise Confusion (' file_type '): ' titles]);

    
    
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