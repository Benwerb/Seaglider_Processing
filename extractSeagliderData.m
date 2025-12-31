function [s, data] = extractSeagliderData(dir_path, missionID)
    % Extract variables from multiple Seaglider NetCDF files and return as struct
    % Input: dir_path - directory containing p686*.nc files
    % Output: s - struct with all variables as cell arrays (one cell per profile)
    %         data - processed data structure with CTD and GPS data
    % ***Notes***
    % The ascent and descent profiles are flipped and I can't figure
    % out why. Temp fix is to switch the logical so grab the ascent as time
    % < start of climb. Should figure this issue out in the future BW
    % 12/10/2025

    % Variables to extract from .nc files
    variables_to_process = {...
        'log_gps_time',... % first fix post-dive, last fix pre-dive, first fix post-dive
        'log_gps_lat',... % first fix post-dive, last fix pre-dive, first fix post-dive
        'log_gps_lon',... % first fix post-dive, last fix pre-dive, first fix post-dive
        'log_dive',... % dive number
        'fet_time',... % time of fet sample
        'fet_Vrse',... % sensor value from fet
        'fet_Vrsestd',... % sensor value from fet
        'fet_Vk',... % sensor value from fet
        'fet_Vkstd',... % sensor value from fet
        'fet_Ik',... % sensor value from fet
        'fet_Ib',... % sensor value from fet
        'ctd_pressure', ... % sensor value from CTD
        'ctd_depth',... % sensor value from CTD
        'ctd_time',... % time of CTD sample
        'temperature',... % temperature (in situ) corrected for thermistor first-order lag
        'salinity',... % salinity corrected for thermal-inertia effects (PSU)
        'conductivity',... % Conductivity corrected for anomalies
        'aanderaa4831_dissolved_oxygen',... % Dissolved oxygen concentration
        'aanderaa4831_results_time',... % time of aanderaa optode sample
        'temperature_qc',... % quality control flags
        'salinity_qc',... % quality control flags
        'aanderaa4831_qc',... % quality control flags
        'aanderaa4831_drift_gain',... % gain correction for aanderaa optode
        'latlong_qc',... % quality control flags
        'CTD_qc',... % quality control flags
        'sigma_theta',... % sigma theta sea water
        'theta',... % theta sea water
        'density',... % rho sea water
        'start_of_climb_time' % elapsed seconds after dive start
    };
    % Get list of all pxxx*.nc files, i.e. xxx: 686
    tag = sprintf('p%s*.nc',missionID(4:6));
    file_pattern = fullfile(dir_path, tag);
    
    % file_pattern = fullfile(dir_path, 'p686*.nc');
    files = dir(file_pattern);
    
    if isempty(files)
        error('No matching files found in %s', dir_path);
    end
    
    % fprintf('Found %d NetCDF files\n', length(files));
    
    % Initialize output struct with cell arrays (column)
    s = struct();
    s.profile_numbers = [];
    for v = 1:length(variables_to_process)
        var_name = variables_to_process{v};
        s.(var_name) = cell(length(files), 1);  % Pre-allocate nx1 cell array
    end
    
    % Initialize counter for valid profiles
    profile_count = 0;
    
    % Extract data from each file
    % fprintf('\nExtracting data from files...\n');
    for i = 1:length(files)
        file_path = fullfile(dir_path, files(i).name);
        file_name = files(i).name;
        
        % Extract profile number from filename (p686xxxx.nc)
        profile_num_str = regexp(file_name, 'p686(\d+)\.nc', 'tokens');
        
        if ~isempty(profile_num_str)
            profile_num = str2double(profile_num_str{1}{1});
            profile_count = profile_count + 1;
            s.profile_numbers(profile_count, 1) = profile_num;
            
            % fprintf('  Profile %04d:\n', profile_num);
            
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
                    
                    % fprintf('    %s: %d values\n', var_name, length(data_read));
                    
                catch ME
                    % Store empty array if variable not found
                    s.(var_name){profile_count} = [];
                    % fprintf('    %s: NOT FOUND or ERROR\n', var_name);
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
    
    % fprintf('\nProcessing complete! Struct contains %d variables across %d profiles.\n', ...
        % length(variables_to_process), length(s.profile_numbers));
    
    % Create standardized data structure
    % fprintf('\nCreating standardized data structure...\n');
    data = struct();

    % GPS data structure
    data.gps = struct();
    
    % Extract GPS data (2nd value = dive start, 3rd value = dive end)
    n_profiles = length(s.log_gps_time);
    
    data.gps.time = struct();
    data.gps.time.divestart = zeros(n_profiles, 1);
    data.gps.time.diveend = zeros(n_profiles, 1);
    
    data.gps.lat = struct();
    data.gps.lat.divestart = zeros(n_profiles, 1);
    data.gps.lat.diveend = zeros(n_profiles, 1);
    
    data.gps.lon = struct();
    data.gps.lon.divestart = zeros(n_profiles, 1);
    data.gps.lon.diveend = zeros(n_profiles, 1);
    
    data.gps.ndive = struct();
    data.gps.ndive.divestart = nan(n_profiles, 1);
    data.gps.ndive.diveend = nan(n_profiles, 1);
    
    % Loop through profiles to extract GPS values
    for i = 1:n_profiles
        % Time
        if length(s.log_gps_time{i}) >= 3
            data.gps.time.divestart(i) = s.log_gps_time{i}(2);
            data.gps.time.diveend(i) = s.log_gps_time{i}(3);
        else
            data.gps.time.divestart(i) = NaN;
            data.gps.time.diveend(i) = NaN;
        end
        
        % Latitude
        if length(s.log_gps_lat{i}) >= 3
            data.gps.lat.divestart(i) = s.log_gps_lat{i}(2);
            data.gps.lat.diveend(i) = s.log_gps_lat{i}(3);
        else
            data.gps.lat.divestart(i) = NaN;
            data.gps.lat.diveend(i) = NaN;
        end
        
        % Longitude
        if length(s.log_gps_lon{i}) >= 3
            data.gps.lon.divestart(i) = s.log_gps_lon{i}(2);
            data.gps.lon.diveend(i) = s.log_gps_lon{i}(3);
        else
            data.gps.lon.divestart(i) = NaN;
            data.gps.lon.diveend(i) = NaN;
        end
    end
    
    % Empty GPS fields (not in NetCDF files)
    data.gps.tfix = struct();
    data.gps.nsat = struct();
    data.gps.minsnr = struct();
    data.gps.meansnr = struct();
    data.gps.maxsnr = struct();
    data.gps.hdop = struct();
    data.gps.gpsstat = struct();
    data.gps.wingstat = struct();
    data.gps.istat = struct();
    
    % CTD data structure
    data.ctd = struct();
    data.ctd.info = struct();
    data.ctd.p = s.ctd_pressure;
    data.ctd.t = s.temperature;
    data.ctd.c = s.conductivity;
    data.ctd.n = nan;
    data.ctd.ndive = cell2mat(s.log_dive);
    data.ctd.depth = s.ctd_depth;
    data.ctd.s = s.salinity;
    data.ctd.theta = s.theta;
    data.ctd.sigma = s.sigma_theta;
    data.ctd.rho = s.density; % double check this is correct with OSU. Needs to match ctdvars_sp2.m: data.ctd.rho{idive}=sw_dens(data.ctd.s{idive},data.ctd.t{idive},data.ctd.p{idive})-1000;
    % CTD time is relative to divestart (for bin_sp2.m compatibility)
    data.ctd.time = cell(size(s.ctd_time));
    data.ctd.start_of_climb_time = cell(size(s.start_of_climb_time));
    for i = 1:length(s.ctd_time)
        if ~isempty(s.ctd_time{i})
            data.ctd.time{i} = s.ctd_time{i} - data.gps.time.divestart(i);
            data.ctd.start_of_climb_time{i} = s.start_of_climb_time{i};
        else
            data.ctd.time{i} = [];
            data.ctd.start_of_climb_time{i} = [];
        end
    end

    % Determine cast phase (upcast/downcast) for each profile
    % Phase = 1 for downcast (descent), 0 for upcast (ascent)
    data.ctd.phase = cell(size(data.ctd.p));
    for i = 1:length(data.ctd.p)
        if ~isempty(data.ctd.p{i})
            data.ctd.phase{i} = data.ctd.time{i} < data.ctd.start_of_climb_time{i};
        else
            data.ctd.phase{i} = [];
        end
    end
    
    % pH data structure
    data.ph = struct();
    data.ph.info = struct();
    data.ph.Vrse = s.fet_Vrse;
    data.ph.vrse_std = s.fet_Vrsestd;
    data.ph.vk = s.fet_Vk;
    data.ph.vk_std = s.fet_Vkstd;
    data.ph.ik = s.fet_Ik;
    data.ph.ib = s.fet_Ib;
    data.ph.ndive = cell2mat(s.log_dive);
    data.ph.p = cell(size(s.ctd_pressure));  % Removed duplicate line    
    data.ph.depth = cell(size(s.ctd_pressure));  % Removed duplicate line    
    % fet time relative to divestart (for bin_sp2 compatibility)
    data.ph.time = cell(size(s.fet_time));
    for i = 1:length(s.fet_time)
        if ~isempty(s.fet_time{i})
            data.ph.time{i} = s.fet_time{i} - data.gps.time.divestart(i);
        else
            data.ph.time{i} = [];
        end
    end
    
    % No pressure variable for fet. Interp p from ctd onto fet time grid
    for i = 1:length(s.ctd_pressure)
        if ~isempty(s.ctd_pressure{i}) && ~isempty(data.ph.time{i})
            % Remove duplicate times while preserving order
            [~, ia] = unique(data.ctd.time{i}, 'stable');
            
            % Only interpolate if we have valid data
            if length(ia) > 1
                data.ph.p{i} = interp1(data.ctd.time{i}(ia), data.ctd.p{i}(ia), ...
                                       data.ph.time{i}, 'linear', 'extrap');
                data.ph.depth{i} = interp1(data.ctd.time{i}(ia), data.ctd.depth{i}(ia), ...
                                       data.ph.time{i}, 'linear', 'extrap');
                % data.ph.depth{i} = sw_dpth(data.ph.p{i},data.gps.lat.divestart(i));
            else
                data.ph.p{i} = [];
                data.ph.depth{i} = [];
            end
        else
            data.ph.p{i} = [];
            data.ph.depth{i} = [];
        end
    end
    % Determine cast phase (upcast/downcast) for each profile
    % Phase = 1 for downcast (ascent), 0 for upcast (ascent)
    data.ph.phase = cell(size(data.ph.p));
    for i = 1:length(data.ph.p)
        if ~isempty(data.ph.p{i})
            data.ph.phase{i} = data.ph.time{i} < data.ctd.start_of_climb_time{i};
        else
            data.ph.phase{i} = [];
        end
    end
    %% Get pH cal
    try
        databaseName = 'glidata';
        username = 'glidata_user';
        password = 'IoFUTBeaQDppSYcmBebA4rV8SJOEMCFI';
        server = 'dpg-d2jobg3e5dus738ce5vg-a.oregon-postgres.render.com';
        port = '5432';
        conn = connect_glidata(databaseName, ...
                        username,password,server,port);
        query = sprintf("SELECT * FROM sensor_deployment_view " + ...
           "WHERE mission_id = '%s' " + ...
           "ORDER BY date DESC LIMIT 1", missionID);
        pH_cal = fetch(conn, query);
        if isopen(conn)
            close(conn)
        end
        disp('got pH calibration')
        k0 = pH_cal.k0_seawater;
        k2 = pH_cal.k2_fp_c0;
        Pcoefs = [pH_cal.fp_k1,pH_cal.fp_k2,pH_cal.fp_k3,...
                    pH_cal.fp_k4,pH_cal.fp_k5,pH_cal.fp_k6]';
    catch
        k0 = NaN;
        k2 = NaN;
        Pcoefs = [NaN, NaN, NaN, NaN, NaN, NaN];
    end
    % Calculate pH
    % interpolate t,s to p of pressure and compute ph
    ndive = length(data.ph.p);
    data.ph.ph = cell(ndive,1);
    data.ph.s = cell(ndive, 1);
    data.ph.t = cell(ndive, 1);
    for n = 1:ndive
        if ~isempty(data.ph.Vrse{n})
            [~,iuse] = unique(data.ctd.time{n});
            if sum(iuse)>1
                ss = interp1(data.ctd.time{n}(iuse),data.ctd.s{n}(iuse),data.ph.time{n},'linear','extrap'); % interpolate in time rather than pressure since time is monotonic
                tt = interp1(data.ctd.time{n}(iuse),data.ctd.t{n}(iuse),data.ph.time{n},'linear','extrap');
                [~,data.ph.ph{n}] = phcalc_jp(data.ph.Vrse{n},data.ph.p{n},tt,ss,k0,k2,Pcoefs);
                data.ph.s{n} = ss(:);
                data.ph.t{n} = tt(:);
            else
                data.ph.ph{n} = nan(size(data.ph.Vrse{n}));
                data.ph.s{n} = nan(size(data.ph.Vrse{n}));
                data.ph.t{n} = nan(size(data.ph.Vrse{n}));
            end
        end
    end

    % Dissolved Oxygen data structure
    data.dox = struct();
    data.dox.info = struct();
    data.dox.oxphase = cell(size(s.ctd_pressure));
    data.dox.thermv = cell(size(s.ctd_pressure));
    data.dox.oxconc = s.aanderaa4831_dissolved_oxygen;
    data.dox.thermt = cell(size(s.ctd_pressure));
    data.dox.ndive = cell2mat(s.log_dive);
    data.dox.p = cell(size(s.ctd_pressure));  
    data.dox.depth = cell(size(s.ctd_pressure));
    data.dox.ox = s.aanderaa4831_dissolved_oxygen;
    data.dox.oxumolkg = s.aanderaa4831_dissolved_oxygen;
    % fet time relative to divestart (for bin_sp2 compatibility)
    data.dox.time = cell(size(s.aanderaa4831_results_time));
    for i = 1:length(s.aanderaa4831_results_time)
        if ~isempty(s.aanderaa4831_results_time{i})
            data.dox.time{i} = s.aanderaa4831_results_time{i} - data.gps.time.divestart(i);
        else
            data.dox.time{i} = [];
        end
    end
    
    % No pressure variable for fet. Interp p from ctd onto fet time grid
    for i = 1:length(s.ctd_pressure)
        if ~isempty(s.ctd_pressure{i}) && ~isempty(data.dox.time{i})
            % Remove duplicate times while preserving order
            [~, ia] = unique(data.ctd.time{i}, 'stable');
            
            % Only interpolate if we have valid data
            if length(ia) > 1
                data.dox.p{i} = interp1(data.ctd.time{i}(ia), data.ctd.p{i}(ia), ...
                                       data.dox.time{i}, 'linear', 'extrap');
                data.dox.depth{i} = interp1(data.ctd.time{i}(ia), data.ctd.depth{i}(ia), ...
                                       data.dox.time{i}, 'linear', 'extrap');
                % data.ph.depth{i} = sw_dpth(data.ph.p{i},data.gps.lat.divestart(i));
            else
                data.dox.p{i} = [];
                data.dox.depth{i} = [];
            end
        else
            data.dox.p{i} = [];
            data.dox.depth{i} = [];
        end
    end
    % Determine cast phase (upcast/downcast) for each profile
    % Phase = 1 for downcast (ascent), 0 for upcast (ascent)
    data.dox.phase = cell(size(data.dox.p));
    for i = 1:length(data.dox.p)
        if ~isempty(data.dox.p{i})
            data.dox.phase{i} = data.dox.time{i} < data.ctd.start_of_climb_time{i};
        else
            data.dox.phase{i} = [];
        end
    end


    % Quality control structure
    data.qual = struct();
    
    % GPS quality control
    data.qual.gps = struct();
    data.qual.gps.divestart = zeros(n_profiles, 1);
    data.qual.gps.diveend = zeros(n_profiles, 1);
    
    % Extract QC flags and convert (1 or less -> 0, otherwise keep value)
    % I am spoofing QC here for now by assuming everything is 0
    for i = 1:n_profiles
        if ~isempty(s.latlong_qc{i}) && length(s.latlong_qc{i}) >= 1
            qc_val = s.latlong_qc{i}(1);
            data.qual.gps.divestart(i) = (qc_val <= 1) * 0 + (qc_val > 1) * qc_val;
            data.qual.gps.diveend(i) = (qc_val <= 1) * 0 + (qc_val > 1) * qc_val;
        else
            data.qual.gps.divestart(i) = NaN;
            data.qual.gps.diveend(i) = NaN;
        end
    end
    
    % CTD quality control (all zeros, same structure as data.ctd)
    data.qual.ctd = struct();
    data.qual.ctd.t = cell(size(data.ctd.t));
    data.qual.ctd.c = cell(size(data.ctd.c));
    data.qual.ctd.p = cell(size(data.ctd.p));
    data.qual.ctd.s = cell(size(data.ctd.s));
    
    for i = 1:length(data.ctd.t)
        data.qual.ctd.t{i} = zeros(size(data.ctd.t{i}));
        data.qual.ctd.c{i} = zeros(size(data.ctd.c{i}));
        data.qual.ctd.p{i} = zeros(size(data.ctd.p{i}));
        data.qual.ctd.s{i} = zeros(size(data.ctd.s{i}));
    end
    
    % Depth average QC
    data.qual.u = zeros(size(data.ctd.t));
    data.qual.v = zeros(size(data.ctd.t));

    % pH QC
    data.qual.ph = struct();
    data.qual.ph.ph = cell(size(data.ctd.t));
    for i = 1:length(data.ph.Vrse)
        data.qual.ph.ph{i} = zeros(size(data.ph.Vrse{i}));
    end

    % Dissolved Oxygen QC
    data.qual.dox = struct();
    data.qual.dox.ox = cell(size(data.ctd.t));
    for i = 1:length(data.dox.oxumolkg)
        data.qual.dox.ox{i} = zeros(size(data.dox.oxumolkg{i}));
    end

    % data.opt = struct();

    % Top-level convenience arrays (nx2: column 1 = divestart, column 2 = diveend)
    data.lat = double([data.gps.lat.divestart, data.gps.lat.diveend]);
    data.lon = double([data.gps.lon.divestart, data.gps.lon.diveend]);
    data.time = double([data.gps.time.divestart, data.gps.time.diveend]);
    
    data.u = NaN(size(data.ctd.t));
    data.v = NaN(size(data.ctd.t));

    fprintf('Phase determination complete.\n');
end