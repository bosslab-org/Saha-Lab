function Read_Intan_v2
tic;
filepath = '/Users/Xander/Documents/MATLAB/Neural_Recordings';
master_file_path = 'Master_files';
filepath2 = 'RHD_files';
filepath3 = 'Odors';
date = '08_19_2021';
Position = '2';

trim_time = [-5 14];    % trim time relative to stimulus onset
total_time = numel(trim_time(1):trim_time(2)-1);

%%%%%%%%%%%%%%%%Adjust master_matrix to appropriate variable
master_data = 'Oor_RMS_Data.mat';
master_data_filt = 'Oor_RMS_Data_Filt.mat';
write_path = 'MAT_files';
trial = 5;
trial_selection = zeros(5,4);                   % pre-allocate 5x4 matrix of zeros
% trial_selection(1,:) = [4 8 11 220000];         % 11 second trial
trial_selection(2,:) = [5 9 19 400000];         % 19 second trial, 20 second recording (1 sec. pre-buffer + 19 sec. recording + 1 sec. post-buffer)
trial_selection(3,:) = [10 14 36 740000];       % 36 second trial, 37 second recording
trial_selection(4,:) = [10 14 33 680000];       % 33 second trial, 36 second recording
trial_selection(5,:) = [15 19 40 820000];       % 40 second trial, 43 second recording (1 sec. pre-buffer + 41 sec. recording + 1 sec. post-buffer)  
trial_selection = trial_selection(trial,:);

trial_time = trial_selection(1,3);
trimmed_trial_length = trial_selection(1,4);

date_temp = datestr(datetime(date, 'InputFormat','MM_dd_yyyy'));            % return date of experiment
date_temp = [date_temp(4:6) , '_' , date_temp(1:2) , '_' , date_temp(10:11)];                   % reformat date for saving purposes 
files = dir([filepath '/' filepath2 '/' filepath3 '/' date '/Position_' Position]);
files = files(~ismember({files.name},{'.','..','.DS_Store'}));

fcutlow = 300;                  % low cut off frequency in Hz
bin_size = 50;                  % bin size in milliseconds

bins_per_sec = 1000/bin_size;   % convert milliseconds to seconds
for cycle_odors = 4 %1:length(files)        
    Odorant{cycle_odors} = [files(cycle_odors).name];
    Odor = Odorant{cycle_odors}(1:end-3);
    experiment = [date, '/', 'Position_', Position, '/' , Odorant{cycle_odors}];
    trials = dir([filepath '/' filepath2 '/' filepath3 '/' experiment '/' '*.rhd']);
    trials = {trials.name};
    trial_idx = 1:length(trials);
    fprintf(['Reading ' Odorant{cycle_odors} '...']); fprintf(1, '\n');
    data_bin = [];
    data_filt_bin = [];
%%%%%%%%%%%%%%%%Above is custom code added to Intan read_file%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Below is read_file code provided by Intan%%%%%%%%%%%%%%%%
    for cycle_trials = trial_idx
        filename = [filepath '/' filepath2  '/' filepath3 '/' experiment '/' trials{cycle_trials}];
        fid = fopen(filename, 'r');
        s = dir(filename);
        filesize = s.bytes;
        magic_number = fread(fid, 1, 'uint32');
        data_file_main_version_number = fread(fid, 1, 'int16');
        data_file_secondary_version_number = fread(fid, 1, 'int16');
        if (data_file_main_version_number == 1)
            num_samples_per_data_block = 60;
        else
            num_samples_per_data_block = 128;
        end    
        sample_rate = fread(fid, 1, 'single');
        dsp_enabled = fread(fid, 1, 'int16');
        actual_dsp_cutoff_frequency = fread(fid, 1, 'single');
        actual_lower_bandwidth = fread(fid, 1, 'single');
        actual_upper_bandwidth = fread(fid, 1, 'single');
        desired_dsp_cutoff_frequency = fread(fid, 1, 'single');
        desired_lower_bandwidth = fread(fid, 1, 'single');
        desired_upper_bandwidth = fread(fid, 1, 'single');
        notch_filter_mode = fread(fid, 1, 'int16');
        notch_filter_frequency = 60;
        desired_impedance_test_frequency = fread(fid, 1, 'single');
        actual_impedance_test_frequency = fread(fid, 1, 'single');
        notes = struct('note1', fread_QString(fid),'note2', fread_QString(fid),'note3', fread_QString(fid) );
        num_temp_sensor_channels = fread(fid, 1, 'int16');
        board_mode = fread(fid, 1, 'int16');
        if (data_file_main_version_number > 1)
            reference_channel = fread_QString(fid);
        end
        frequency_parameters = struct('amplifier_sample_rate', sample_rate,'aux_input_sample_rate', sample_rate / 4,'supply_voltage_sample_rate', sample_rate / num_samples_per_data_block,'board_adc_sample_rate', sample_rate,'board_dig_in_sample_rate', sample_rate,'desired_dsp_cutoff_frequency', desired_dsp_cutoff_frequency,'actual_dsp_cutoff_frequency', actual_dsp_cutoff_frequency,'dsp_enabled', dsp_enabled,'desired_lower_bandwidth', desired_lower_bandwidth,'actual_lower_bandwidth', actual_lower_bandwidth,'desired_upper_bandwidth', desired_upper_bandwidth,'actual_upper_bandwidth', actual_upper_bandwidth,'notch_filter_frequency', notch_filter_frequency,'desired_impedance_test_frequency', desired_impedance_test_frequency,'actual_impedance_test_frequency', actual_impedance_test_frequency );
        spike_trigger_struct = struct('voltage_trigger_mode', {},'voltage_threshold', {},'digital_trigger_channel', {},'digital_edge_polarity', {} );
        new_trigger_channel = struct(spike_trigger_struct);
        spike_triggers = struct(spike_trigger_struct);
        channel_struct = struct('native_channel_name', {},'custom_channel_name', {},'native_order', {},'custom_order', {},'board_stream', {},'chip_channel', {},'port_name', {},'port_prefix', {},'port_number', {},'electrode_impedance_magnitude', {},'electrode_impedance_phase', {} );
        new_channel = struct(channel_struct);
        amplifier_channels = struct(channel_struct);
        aux_input_channels = struct(channel_struct);
        supply_voltage_channels = struct(channel_struct);
        board_adc_channels = struct(channel_struct);
        board_dig_in_channels = struct(channel_struct);
        board_dig_out_channels = struct(channel_struct);
        amplifier_index = 1;
        aux_input_index = 1;
        supply_voltage_index = 1;
        board_adc_index = 1;
        board_dig_in_index = 1;
        board_dig_out_index = 1;
        number_of_signal_groups = fread(fid, 1, 'int16');
        for signal_group = 1:number_of_signal_groups
            signal_group_name = fread_QString(fid);
            signal_group_prefix = fread_QString(fid);
            signal_group_enabled = fread(fid, 1, 'int16');
            signal_group_num_channels = fread(fid, 1, 'int16');
            signal_group_num_amp_channels = fread(fid, 1, 'int16');
            if (signal_group_num_channels > 0 && signal_group_enabled > 0)
                new_channel(1).port_name = signal_group_name;
                new_channel(1).port_prefix = signal_group_prefix;
                new_channel(1).port_number = signal_group;
                for signal_channel = 1:signal_group_num_channels
                    new_channel(1).native_channel_name = fread_QString(fid);
                    new_channel(1).custom_channel_name = fread_QString(fid);
                    new_channel(1).native_order = fread(fid, 1, 'int16');
                    new_channel(1).custom_order = fread(fid, 1, 'int16');
                    signal_type = fread(fid, 1, 'int16');
                    channel_enabled = fread(fid, 1, 'int16');
                    new_channel(1).chip_channel = fread(fid, 1, 'int16');
                    new_channel(1).board_stream = fread(fid, 1, 'int16');
                    new_trigger_channel(1).voltage_trigger_mode = fread(fid, 1, 'int16');
                    new_trigger_channel(1).voltage_threshold = fread(fid, 1, 'int16');
                    new_trigger_channel(1).digital_trigger_channel = fread(fid, 1, 'int16');
                    new_trigger_channel(1).digital_edge_polarity = fread(fid, 1, 'int16');
                    new_channel(1).electrode_impedance_magnitude = fread(fid, 1, 'single');
                    new_channel(1).electrode_impedance_phase = fread(fid, 1, 'single');
                    if (channel_enabled)
                        switch (signal_type)
                            case 0
                                amplifier_channels(amplifier_index) = new_channel;
                                spike_triggers(amplifier_index) = new_trigger_channel;
                                amplifier_index = amplifier_index + 1;
                            case 1
                                aux_input_channels(aux_input_index) = new_channel;
                                aux_input_index = aux_input_index + 1;
                            case 2
                                supply_voltage_channels(supply_voltage_index) = new_channel;
                                supply_voltage_index = supply_voltage_index + 1;
                            case 3
                                board_adc_channels(board_adc_index) = new_channel;
                                board_adc_index = board_adc_index + 1;
                            case 4
                                board_dig_in_channels(board_dig_in_index) = new_channel;
                                board_dig_in_index = board_dig_in_index + 1;
                            case 5
                                board_dig_out_channels(board_dig_out_index) = new_channel;
                                board_dig_out_index = board_dig_out_index + 1;
                            otherwise
                                error('Unknown channel type');
                        end
                    end
                end
            end
        end
        num_amplifier_channels = amplifier_index - 1;
        num_aux_input_channels = aux_input_index - 1;
        num_supply_voltage_channels = supply_voltage_index - 1;
        num_board_adc_channels = board_adc_index - 1;
        num_board_dig_in_channels = board_dig_in_index - 1;
        num_board_dig_out_channels = board_dig_out_index - 1;
        bytes_per_block = num_samples_per_data_block * 4;  % timestamp data
        bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_amplifier_channels;
        bytes_per_block = bytes_per_block + (num_samples_per_data_block / 4) * 2 * num_aux_input_channels;
        bytes_per_block = bytes_per_block + 1 * 2 * num_supply_voltage_channels;
        bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_board_adc_channels;
        if (num_board_dig_in_channels > 0)
            bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
        end
        if (num_board_dig_out_channels > 0)
            bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
        end
        data_present = 0;
        bytes_remaining = filesize - ftell(fid);
        if (bytes_remaining > 0)
            data_present = 1;
        end
        num_data_blocks = bytes_remaining / bytes_per_block;
        num_amplifier_samples = num_samples_per_data_block * num_data_blocks;
        num_aux_input_samples = (num_samples_per_data_block / 4) * num_data_blocks;
        num_supply_voltage_samples = 1 * num_data_blocks;
        num_board_adc_samples = num_samples_per_data_block * num_data_blocks;
        num_board_dig_in_samples = num_samples_per_data_block * num_data_blocks;
        num_board_dig_out_samples = num_samples_per_data_block * num_data_blocks;
        record_time = num_amplifier_samples / sample_rate;
        t_amplifier = zeros(1, num_amplifier_samples);
        amplifier_data = zeros(num_amplifier_channels, num_amplifier_samples);
        aux_input_data = zeros(num_aux_input_channels, num_aux_input_samples);
        supply_voltage_data = zeros(num_supply_voltage_channels, num_supply_voltage_samples);
        board_adc_data = zeros(num_board_adc_channels, num_board_adc_samples);
        board_dig_in_data = zeros(num_board_dig_in_channels, num_board_dig_in_samples);
        board_dig_in_raw = zeros(1, num_board_dig_in_samples);
        board_dig_out_data = zeros(num_board_dig_out_channels, num_board_dig_out_samples);
        board_dig_out_raw = zeros(1, num_board_dig_out_samples);
        amplifier_index = 1;
        aux_input_index = 1;
        supply_voltage_index = 1;
        board_adc_index = 1;
        board_dig_in_index = 1;
        board_dig_out_index = 1;
        for i=1:num_data_blocks
            if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 2) || (data_file_main_version_number > 1))
                t_amplifier(amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'int32');
            else
                t_amplifier(amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint32');
            end
            if (num_amplifier_channels > 0)
                amplifier_data(:, amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_amplifier_channels], 'uint16')';
            end
            if (num_aux_input_channels > 0)
                aux_input_data(:, aux_input_index:(aux_input_index + (num_samples_per_data_block / 4) - 1)) = fread(fid, [(num_samples_per_data_block / 4), num_aux_input_channels], 'uint16')';
            end
            if (num_supply_voltage_channels > 0)
                supply_voltage_data(:, supply_voltage_index) = fread(fid, [1, num_supply_voltage_channels], 'uint16')';
            end
            if (num_board_adc_channels > 0)
                board_adc_data(:, board_adc_index:(board_adc_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_board_adc_channels], 'uint16')';
            end
            if (num_board_dig_in_channels > 0)
                board_dig_in_raw(board_dig_in_index:(board_dig_in_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint16');
            end
            if (num_board_dig_out_channels > 0)
                board_dig_out_raw(board_dig_out_index:(board_dig_out_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint16');
            end

            amplifier_index = amplifier_index + num_samples_per_data_block;
            aux_input_index = aux_input_index + (num_samples_per_data_block / 4);
            supply_voltage_index = supply_voltage_index + 1;
            board_adc_index = board_adc_index + num_samples_per_data_block;
            board_dig_in_index = board_dig_in_index + num_samples_per_data_block;
            board_dig_out_index = board_dig_out_index + num_samples_per_data_block;
        end
        fclose(fid);
        for i=1:num_board_dig_in_channels
           mask = 2^(board_dig_in_channels(i).native_order) * ones(size(board_dig_in_raw));
           board_dig_in_data(i, :) = (bitand(board_dig_in_raw, mask) > 0);
        end
        amplifier_data = 0.195 * (amplifier_data - 32768); % units = microvolts        
        num_gaps = sum(diff(t_amplifier) ~= 1);
        if (num_gaps ~= 0)
            fprintf(1, 'Warning: %d gaps in timestamp data found.  Time scale will not be uniform!\n', num_gaps);
        end
        
%%%%%%%%%%%%%%%%Above is read_file code provided by Intan%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Below is custom code added to Intan read_file%%%%%%%%%%%%

        channels_idx = [amplifier_channels.custom_order];
        [~,channels_idx] = sort(channels_idx, 'ascend');
        Ny_lim = sample_rate/2;                         % Nyquist limit to avoid aliasing
        [b,a] = butter(6,fcutlow/Ny_lim,'high');        % High pass butterworth filter transfer function coefficients

        if num_board_dig_in_channels == 1               % Manual stimulus time selection
            stimulus_onset = trial_selection(1,1);
            stimulus_offset = trial_selection(1,2);
            fprintf('Manual stimulus selection enabled. ');
        elseif num_board_dig_in_channels == 2           % Read external signal for stimulus time
            stimulus_onset = round((find(board_dig_in_data(2,:), 1, 'first')-sample_rate)/sample_rate);
            stimulus_offset = round((find(board_dig_in_data(2,:), 1, 'last')-sample_rate)/sample_rate);
        end
        ch_per_tetrode = 4;
        samples_per_bin = sample_rate/bins_per_sec;
        num_bins = bins_per_sec*total_time;
        for cycle_tetrodes = 1:(num_amplifier_channels/ch_per_tetrode)
            tetrode_channels = ch_per_tetrode*(cycle_tetrodes-1)+1:ch_per_tetrode*cycle_tetrodes;
            for cycle_channels=1:ch_per_tetrode
                if data_file_main_version_number < 3
                   data(cycle_tetrodes,cycle_channels,cycle_trials,:) =  notch_filter(amplifier_data(channels_idx(tetrode_channels(cycle_channels)),(sample_rate+1:trimmed_trial_length)), sample_rate, notch_filter_frequency, 10);
                else
                   data(cycle_tetrodes,cycle_channels,cycle_trials,:) =  amplifier_data(channels_idx(tetrode_channels(cycle_channels)),(sample_rate+1:trimmed_trial_length));  
                end
                data_filt(cycle_tetrodes,cycle_channels,cycle_trials,:) = filtfilt(b,a,squeeze(data(cycle_tetrodes,cycle_channels,cycle_trials,:)));
                
                data_bin(cycle_tetrodes,cycle_channels,cycle_trials,:) = rms(reshape(data(cycle_tetrodes,cycle_channels,cycle_trials,(stimulus_onset+trim_time(1))*sample_rate+1:(stimulus_onset+trim_time(2))*sample_rate),samples_per_bin,num_bins));                
                data_filt_bin(cycle_tetrodes,cycle_channels,cycle_trials,:) = rms(reshape(data_filt(cycle_tetrodes,cycle_channels,cycle_trials,(stimulus_onset+trim_time(1))*sample_rate+1:(stimulus_onset+trim_time(2))*sample_rate),samples_per_bin,num_bins));                
                
                Date_temp{cycle_tetrodes,1} = date_temp;
%                 data_bin(cycle_tetrodes,cycle_channels,cycle_trials,:) = rms(reshape(data(cycle_tetrodes,cycle_channels,cycle_trials,:),samples_per_bin,num_bins));
%                 data_filt_bin(cycle_tetrodes,cycle_channels,cycle_trials,:) = rms(reshape(data_filt(cycle_tetrodes,cycle_channels,cycle_trials,:),samples_per_bin,num_bins));
            end
        end
            
        fprintf('Finished file %.0f of %.0f. ', cycle_trials, length(trials))
        fprintf(1, 'Channels: %.0f. %0.3f seconds of data. ', num_amplifier_channels, record_time);
        fprintf(1, 'Sampling rate: %0.2f kHz. \n',sample_rate / 1000);

    end
    
    var_name_data = [Odor '_data'];
    eval([var_name_data '= data;'])
    
    var_name_data_filt = [Odor '_data_filt'];
    eval([var_name_data_filt '= data_filt;'])
    
    if ~exist([filepath '/' write_path '/' filepath3 '/' date '/Position_' Position])
        mkdir([filepath '/' write_path '/' filepath3 '/' date '/Position_' Position]);
    end
    
    if ~exist([filepath '/' master_file_path])
        mkdir([filepath '/' master_file_path]);
    end
%%%%%%%%%%%%%%Save data within odor experiment folder
    fprintf(1, '\nSaving %s data... \n', Odor);    
    save([filepath '/' write_path '/' filepath3 '/' date '/Position_' Position '/' Odor],...
        [Odor '_data'], [Odor '_data_filt'], 'stimulus_onset','stimulus_offset','trial_time','sample_rate');
%%%%%%%%%%%%%%Day-wise and master matrix variables (below)
%     var_name_data_bin = [Odor '_data'];
%     var_name_data_filt_bin = [Odor '_data_filt'];
%     eval([var_name_data_bin '= data_bin;'])
%     eval([var_name_data_filt_bin '= data_filt_bin;'])
%%%%%%%%%%%%%%Save data to day-wise file
%     if cycle_odors == 1
%         fprintf(1, 'Creating new RMS-based experiment file... \n');
%         stimulus_on = abs(trim_time(1));
%         stimulus_off = abs(trim_time(1))+4;
% %         total_time = abs(trim_time(1))+trim_time(2);        
%         save([filepath '/' write_path '/' filepath3 '/' date '/Position_' Position '/Position_' Position],[Odor '_data'], [Odor '_data_filt'], 'stimulus_on','stimulus_off','total_time','sample_rate', 'bin_size');
%     else
%         fprintf(1, 'Appending %s RMS data to experiment file... \n', Odor);
%         save([filepath '/' write_path '/' filepath3 '/' date '/Position_' Position '/Position_' Position],[Odor '_data'], [Odor '_data_filt'],'-append');
%     end
%%%%%%%%%%%%%%Add data to master matrix
%     if ~(exist([filepath '/' master_file_path '/' master_data]))
%         fprintf(1, 'Creating new RMS-based master matrix... \n\n');
%         Date = Date_temp;
%         stimulus_on = abs(trim_time(1));
%         stimulus_off = abs(trim_time(1))+4;
% %         total_time = abs(trim_time(1))+trim_time(2);
%         save([filepath '/' master_file_path '/' master_data], [Odor '_data'], 'Date', 'stimulus_on','stimulus_off','total_time','sample_rate', 'bin_size','-v7.3');
%         save([filepath '/' master_file_path '/' master_data_filt], [Odor '_data_filt'], 'Date', 'stimulus_on','stimulus_off','total_time','sample_rate', 'bin_size','-v7.3');
%     elseif isempty(who('-file', [filepath '/' master_file_path '/' master_data], [Odor '_data']))
%         fprintf(1, 'Adding %s RMS data to master matrix... \n\n', Odor);
%         save([filepath '/' master_file_path '/' master_data],[Odor '_data'],'-append');
%         save([filepath '/' master_file_path '/' master_data_filt], [Odor '_data_filt'],'-append');
%     elseif ~isempty(who('-file', [filepath '/' master_file_path '/' master_data], [Odor '_data']))
%         fprintf(1, 'Appending %s RMS data to master matrix... \n\n', Odor);
% %%%%%%%%%%%%%%Add unfiltered data to master matrix
%         load([filepath '/' master_file_path '/' master_data],'Date',[Odor '_data']);
%         if cycle_odors == 1
%             Date = vertcat(Date, Date_temp);
%             data_bin = cat(1,eval([Odor '_data']),data_bin);
%             eval([var_name_data_bin '= data_bin;'])
%             save([filepath '/' master_file_path '/' master_data],[Odor '_data'], 'Date','-append');
%         else
%             data_bin = cat(1,eval([Odor '_data']),data_bin);
%             eval([var_name_data_bin '= data_bin;'])
%             save([filepath '/' master_file_path '/' master_data],[Odor '_data'],'-append');
%         end
% %%%%%%%%%%%%%%Add filtered data to master matrix
%         load([filepath '/' master_file_path '/' master_data_filt],'Date',[Odor '_data_filt']);
%         if cycle_odors == 1
%             Date = vertcat(Date, Date_temp);
%             data_filt_bin = cat(1,eval([Odor '_data_filt']),data_filt_bin);
%             eval([var_name_data_filt_bin '= data_filt_bin;'])
%             save([filepath '/' master_file_path '/' master_data_filt],[Odor '_data_filt'], 'Date','-append');
%         else
%             data_filt_bin = cat(1,eval([Odor '_data_filt']),data_filt_bin);
%             eval([var_name_data_filt_bin '= data_filt_bin;'])
%             save([filepath '/' master_file_path '/' master_data_filt], [Odor '_data_filt'],'-append');
%         end
%     end
fprintf(1, '%d%% complete \n\n',round((cycle_odors/length(files))*100));
end
toc;
fprintf(['Finished ' date ': Position ' Position '\n\n']);
return

function a = fread_QString(fid)
a = '';
length = fread(fid, 1, 'uint32');
if length == hex2num('ffffffff')
    return;
end
length = length / 2;
for i=1:length
    a(i) = fread(fid, 1, 'uint16');
end
return
function out = notch_filter(in, fSample, fNotch, Bandwidth)

tstep = 1/fSample;
Fc = fNotch*tstep;

L = length(in);

d = exp(-2*pi*(Bandwidth/2)*tstep);
b = (1 + d*d)*cos(2*pi*Fc);
a0 = 1;
a1 = -b;
a2 = d*d;
a = (1 + d*d)/2;
b0 = 1;
b1 = -2*cos(2*pi*Fc);
b2 = 1;

out = zeros(size(in));
out(1) = in(1);  
out(2) = in(2);

for i=3:L
    out(i) = (a*b2*in(i-2) + a*b1*in(i-1) + a*b0*in(i) - a2*out(i-2) - a1*out(i-1))/a0;
end

return