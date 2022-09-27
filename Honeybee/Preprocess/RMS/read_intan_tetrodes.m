function [data_filt, exp_pars] = read_intan_tetrodes(stimulus, filepath, readpath, exp_path)

write_path = 'MAT_files';
fcutlow = 300;                  % low cut off frequency in Hz

stim = stimulus(1:end-3);
trials = dir([filepath '/' readpath '/' exp_path '/' stimulus '/' '*.rhd']);
trials = {trials.name};
fprintf(['Reading ' stimulus '...']); fprintf(1, '\n');

for cycle_trials = 1:length(trials)
    filename = [filepath '/' readpath '/' exp_path '/' stimulus '/' trials{cycle_trials}];
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
%     fprintf(['record time = ' num2str(record_time)])
    
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

    channels_idx = [amplifier_channels.custom_order];
    [~,channels_idx] = sort(channels_idx, 'ascend');
    
    stim_on_temp(cycle_trials) = round((find(board_dig_in_data(2,:), 1, 'first')-sample_rate)/sample_rate,3);
    stim_off_temp(cycle_trials) = round((find(board_dig_in_data(2,:), 1, 'last')-sample_rate)/sample_rate,3);
    
    %%%%%%%%
    if contains(exp_path,'11_02_2021') && contains(stimulus,'PMH(1)') && cycle_trials==5
        zero_pad = zeros(1,940320-size(amplifier_data,2));
        amplifier_data = cat(2,amplifier_data,repmat(zero_pad,4,1));
        board_dig_in_data = cat(2,board_dig_in_data,repmat(zero_pad,2,1));
        record_time = record_time + size(zero_pad,2)/sample_rate;
    end
    %%%%%%%%
    
    trial_end_temp = stim_on_temp(cycle_trials)+(round(record_time)-stim_on_temp(cycle_trials)-0.9);

    fprintf(['Stimulus duration: ' num2str(stim_on_temp(cycle_trials)) ' to ' num2str(stim_off_temp(cycle_trials)) ' seconds \n']);
    
    [b,a] = butter(6,fcutlow/(sample_rate/2),'high');           % High pass butterworth filter transfer function coefficients
    
    ch_per_tetrode = 4;
    for cycle_tetrodes = 1:(num_amplifier_channels/ch_per_tetrode)
        tetrode_channels = ch_per_tetrode*(cycle_tetrodes-1)+1:ch_per_tetrode*cycle_tetrodes;
        for cycle_channels=1:ch_per_tetrode
            if data_file_main_version_number < 3
                data_filt(cycle_tetrodes,cycle_channels,cycle_trials,:) =  filtfilt(b,a,notch_filter(amplifier_data(channels_idx(tetrode_channels(cycle_channels)),(1:trial_end_temp*sample_rate)), sample_rate, notch_filter_frequency, 10));
            else
                data_filt(cycle_tetrodes,cycle_channels,cycle_trials,:) =  filtfilt(b,a,amplifier_data(channels_idx(tetrode_channels(cycle_channels)),(1:trial_end*sample_rate)));  
            end
            [data_processed(cycle_tetrodes,cycle_channels,cycle_trials,:), ~] = process_artifacts(squeeze(data_filt(cycle_tetrodes,cycle_channels,cycle_trials,:)),...
                sample_rate/1000, 8, stim_on_temp(cycle_trials), stim_off_temp(cycle_trials), sample_rate, 0);
        end
    end      
    fprintf('Finished processing trial %.0f of %.0f. ', cycle_trials, length(trials))
    fprintf(1, 'Channels: %.0f. %0.3f seconds of data. ', num_amplifier_channels, record_time);
    fprintf(1, 'Sampling rate: %0.2f kHz. \n',sample_rate / 1000);
end

% Align stimulus onsets to nearest 1 msec
start_time = round((1.1-(stim_on_temp-floor(stim_on_temp)))*sample_rate);
end_trim_time = round((stim_on_temp-floor(stim_on_temp))*sample_rate);
data_filt = data_filt(:,:,:,start_time+1:end-end_trim_time);
data_processed = data_processed(:,:,:,start_time+1:end-end_trim_time);

stim_on = floor(stim_on_temp(1));
stim_off = floor(stim_off_temp(1));
total_time = size(data_filt,4)/sample_rate;
exp_pars = [stim_on, stim_off, sample_rate];

%%%%%%%%%%%%%%Save data_filt%%%%%%%%%%%%%%       
var_name_data_filt = [stim '_data_filt'];
eval([var_name_data_filt '= data_filt;'])  
if ~exist([filepath '/' write_path '/' exp_path],'dir')
    mkdir([filepath '/' write_path '/' exp_path]);
end

% if ~exist([filepath '/' write_path '/' exp_path '/' stim '.mat'],'file')
    fprintf(1, 'Saving %s data... \n', stim);    
    save([filepath '/' write_path '/' exp_path '/' stim], [stim '_data_filt'], 'stim_on','stim_off','total_time','sample_rate');
    fprintf(['Finished ' exp_path ': ' stim '\n\n']);
% else
%     fprintf('data_filt previously saved. MAT file not overwritten.\n')
% end

%%%%%%%%%%%%%%Save data_processed%%%%%%%%%%%%%%       
var_name_data_filt = [stim '_data_processed'];
eval([var_name_data_filt '= data_processed;'])  
if ~exist([filepath '/' write_path '/Processed/' exp_path],'dir')
    mkdir([filepath '/' write_path '/Processed/' exp_path]);
end

% if ~exist([filepath '/' write_path '/Processed/' exp_path '/' stim '.mat'],'file')
    fprintf(1, 'Saving %s data... \n', stim);    
    save([filepath '/' write_path '/Processed/' exp_path '/' stim], [stim '_data_filt'], 'stim_on','stim_off','total_time','sample_rate');
    fprintf(['Finished ' exp_path ': ' stim '\n\n']);
% else
%     fprintf('data_processed previously saved. MAT file not overwritten.\n')
% end
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