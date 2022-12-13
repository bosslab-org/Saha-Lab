function mat_to_igor(filepath, stimulus, day, position, data, sample_rate)

basepath = [filepath '/MAT_to_IGOR'];

for cycle_tetrodes = 1:size(data,1)
    write_path = [basepath '/' day '/Position_' num2str(position) '/Tetrode_' num2str(cycle_tetrodes)];
    if ~(exist(write_path, 'dir'))
        mkdir(write_path)
    end
    for cycle_trials = 1:size(data,3)     % cycle through trials
        for cycle_channels = 1:size(data,2)        % cycle through channels
            cycle_data = squeeze(data(cycle_tetrodes,cycle_channels,cycle_trials,:));
            cycle_data = resample(cycle_data,15000,sample_rate);
            cycle_data = int16(cycle_data);
            fid=fopen([write_path '/' stimulus '_t' num2str(cycle_trials,'%02.2d') '.' num2str(cycle_channels,'%02.2d')],'w','ieee-be');
            fwrite(fid,squeeze(cycle_data),'int16');
            fclose(fid);
        end
    end
end
return