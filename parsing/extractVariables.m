function s = extractVariables(variables_to_process,directory,seagliderSN)
%EXTRACTVARIABLES
% INPUTS: files to process, variables to process
% OUTPUT: struct s that contains each variable formatted as cell array
    % Initialize output struct with cell arrays (column)

    % Get list of all matching '.nc' files
    tag = sprintf('p%d*.nc',seagliderSN);
    file_pattern = fullfile(directory, tag);
    files = dir(file_pattern);
    
    if isempty(files)
        error('No p686*.nc files found in %s', directory);
    end

    s = struct();
    s.profile_numbers = [];
    for v = 1:length(variables_to_process)
        var_name = variables_to_process{v};
        s.(var_name) = cell(length(files), 1);  % Pre-allocate nx1 cell array
    end
    
    % Initialize counter for valid profiles
    profile_count = 0;
    
    % Extract data from each file
    fprintf('\nExtracting data from files...\n');
    for i = 1:length(files)
        file_path = fullfile(directory, files(i).name);
        file_name = files(i).name;
        
        % Extract profile number from filename (p686xxxx.nc)
        profile_num_str = regexp(file_name, 'p686(\d+)\.nc', 'tokens');
        
        if ~isempty(profile_num_str)
            profile_num = str2double(profile_num_str{1}{1});
            profile_count = profile_count + 1;
            s.profile_numbers(profile_count, 1) = profile_num;
            
            fprintf('  Profile %04d:\n', profile_num);
            
            for v = 1:length(variables_to_process)
                var_name = variables_to_process{v};
                
                try
                    % Read the data
                    data_read = ncread(file_path, var_name);
                    
                    % Convert QC flags from char to double
                    if contains(var_name, '_qc') || strcmp(var_name, 'CTD_qc')
                        data_read = str2double(data_read);
                    end
                    
                    % Store as column vector in cell array
                    s.(var_name){profile_count} = data_read(:);
                    
                    fprintf('    %s: %d values\n', var_name, length(data_read));
                    
                catch ME
                    % Store empty array if variable not found
                    s.(var_name){profile_count} = [];
                    fprintf('    %s: NOT FOUND or ERROR\n', var_name);
                end
            end
        end
    end
    
    % Trim cell arrays to actual number of profiles found
    s.profile_numbers = s.profile_numbers(1:profile_count);
    for v = 1:length(variables_to_process)
        var_name = variables_to_process{v};
        s.(var_name) = s.(var_name)(1:profile_count);
    end
    
    fprintf('\nProcessing complete! Struct contains %d variables across %d profiles.\n', ...
        length(variables_to_process), length(s.profile_numbers));
    
end

