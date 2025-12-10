clear all;
%% Notes
% QC:
% 0: QC_NO_CHANGE
% 1: QC_GOOD
% 2: QC_PROBABLY_GOOD
% 3: QC_PROBABLY_BAD
% 4: QC_BAD
% 5: QC_CHANGED
% 6: QC_UNSAMPLED
% 8: QC_INTERPOLATED
% 9: QC_MISSING
%% Goal
% Parse netCDF so that it can be passed into binsat_sp2
%% Processing step 1
% Directory containing NetCDF files
% dir_path = '\\atlas\ftp\incoming\Fred\'; % live path
dir_path = 'C:\Users\bwerb\Documents\GitHub\Seaglider_Processing\data'; % local test path
missionID = '25B68601';
% Variables to extract from .nc files
% t,s already have qc flags applied and bad data is nan
variables_to_process = {...
    'log_gps_time',... % dive end
    'log_gps_lat',... % dive end
    'log_gps_lon',... % dive end
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
    'aanderaa4831_dissolved_oxygen',... % Dissolved oxygen concentration, calculated from optode tcphase corrected for salininty and depth
    'aanderaa4831_results_time',... % time of aanderaa optode sample
    'temperature_qc',... % quality control flags, see top comments
    'salinity_qc',... % quality control flags, see top comments
    'aanderaa4831_qc',... % quality control flags, see top comments
    'aanderaa4831_drift_gain',... % gain correction for aanderaa optode. Check if already applied
    'latlong_qc',... % quality control flags, see top comments
    'CTD_qc',... % quality control flags, see top comments
    'sigma_theta',... % sigma theta sea water
    'start_of_climb_time',... % elapsed seconds after dive start when second (positive) apogee pump starts
    'log_GPS1',... % string reported in logfile for GPS1 fix (first surface position before dive)
    'log_GPS2',... % string reported in logfile for GPS2 fix (last surface position before dive)
    'log_GPS' % String reported in logfile for GPS fix (first surface position after dive)
    };
% 'log_dive'
% 'log_start'
% Get list of all p686*.nc files
file_pattern = fullfile(dir_path, 'p686*.nc');
files = dir(file_pattern);
% files = files(1:2);
% Initialize storage for all variables
variable_data = struct();
for v = 1:length(variables_to_process)
    var_name = variables_to_process{v};
    variable_data.(var_name).max_length = 0;
    variable_data.(var_name).profile_data = {};
    variable_data.(var_name).profile_numbers = [];
end

% Determine maximum array length for each variable
fprintf('Scanning files to determine maximum length for each variable...\n');
for i = 1:length(files)
    file_path = fullfile(dir_path, files(i).name);
    
    for v = 1:length(variables_to_process)
        var_name = variables_to_process{v};
        try
            data = ncread(file_path, var_name);
            variable_data.(var_name).max_length = max(variable_data.(var_name).max_length, length(data));
        catch
            % Variable might not exist in this file, skip silently
        end
    end
end

% Extract data and pad to max_length for each variable
fprintf('\nExtracting data from %d files...\n', length(files));
for i = 1:length(files)
    file_path = fullfile(dir_path, files(i).name);
    file_name = files(i).name;
    
    % Extract profile number from filename (p686xxxx.nc)
    profile_num_str = regexp(file_name, 'p686(\d+)\.nc', 'tokens');
    
    if ~isempty(profile_num_str)
        profile_num = str2double(profile_num_str{1}{1});
        
        fprintf('  Profile %04d:\n', profile_num);
        
        for v = 1:length(variables_to_process)
            var_name = variables_to_process{v};
            max_length = variable_data.(var_name).max_length; % could make this max length of all data
            
            try
                % Read the data
                data = ncread(file_path, var_name);
                
                % Convert QC flags from char to double
                if contains(var_name, '_qc') || strcmp(var_name, 'CTD_qc')
                    data = str2double(data);
                end
            
            % Pad with NaN if shorter than max_length
            if length(data) < max_length
                data = [data(:); NaN(max_length - length(data), 1)];
            end
                
                % Store data and profile number
                variable_data.(var_name).profile_data{end+1} = data;
                variable_data.(var_name).profile_numbers(end+1) = profile_num;
                
                fprintf('    %s: %d values (padded to %d)\n', ...
                        var_name, sum(~isnan(data)), max_length);
            catch ME
                fprintf('    %s: NOT FOUND or ERROR\n', var_name);
            end
        end
    end
end

% Convert to arrays and create output variables
fprintf('\nCreating output arrays...\n');
for v = 1:length(variables_to_process)
    var_name = variables_to_process{v};
    
    if ~isempty(variable_data.(var_name).profile_data)
        % Create matrix from cell array
        data_matrix = cell2mat(variable_data.(var_name).profile_data);
        
        % Assign to workspace with variable name
        assignin('base', var_name, data_matrix);
        
        fprintf('  %s: %d rows x %d columns (assigned to workspace)\n', ...
                var_name, size(data_matrix, 1), size(data_matrix, 2));
    else
        fprintf('  %s: NO DATA EXTRACTED\n', var_name);
    end
end

fprintf('\nProcessing complete!\n');
%%
s_raw = struct;
s_raw.tc = temperature;
s_raw.psal = salinity;
s_raw.ctd_time = ctd_time;
s_raw.depth = ctd_depth;
s_raw.pres = ctd_pressure;
s_raw.fet_time = fet_time;
s_raw.vrse = fet_Vrse;
s_raw.vrse_std = fet_Vrsestd;
s_raw.vk = fet_Vk;
s_raw.vk_std = fet_Vkstd;
s_raw.ik = fet_Ik;
s_raw.ib = fet_Ib;
s_raw.doxy = aanderaa4831_dissolved_oxygen;
s_raw.doxy_time = aanderaa4831_results_time;
s_raw.sigma = sigma_theta;

s_raw.tc_QC = temperature_qc;
s_raw.psal_QC = salinity_qc;
s_raw.doxy_QC = aanderaa4831_qc;
s_raw.aanderaa4831_drift_gain = aanderaa4831_drift_gain;
s_raw.latlong_QC = latlong_qc;
s_raw.CTD_QC = CTD_qc;

s_raw.unixtime = log_gps_time(1,:); % verify this top row is correct
s_raw.sdn = datenum(datetime(s_raw.unixtime, 'ConvertFrom', 'posixtime'));
s_raw.lat = log_gps_lat(1,:);
s_raw.lon = log_gps_lon(1,:);
s_raw.divenum = log_dive;
s_raw.missionID = missionID;

s_raw.start_of_climb_time = start_of_climb_time + s_raw.unixtime(1,:); % need to make sure I add start of dive time
s_raw.castdir = [zeros(1,length(s_raw.depth(1,:))); diff(runmean(s_raw.depth,5))];
% 
% s_raw.castdir(s_raw.castdir < 0) = 1;
% s_raw.castdir(s_raw.castdir >= 0) = 0;
%% Interpolate all variables onto fet_time grid
fprintf('\nInterpolating variables onto fet_time grid...\n');

% List of CTD variables to interpolate (on ctd_time grid)
ctd_vars = {'tc', 'psal', 'depth', 'pres'};
ctd_qc_vars = {'tc_QC','psal_QC','CTD_QC'};
% Initialize interpolated structure
s = struct;

% Copy FET variables (already on fet_time grid)
s.fet_time = s_raw.fet_time;
s.vrse = s_raw.vrse;
s.vrse_std = s_raw.vrse_std;
s.vk = s_raw.vk;
s.vk_std = s_raw.vk_std;
s.ik = s_raw.ik;
s.ib = s_raw.ib;
s.start_of_climb_time = s_raw.start_of_climb_time;
% Interpolate CTD variables from ctd_time to fet_time
for v = 1:length(ctd_vars)
    var_name = ctd_vars{v};
    fprintf('  Interpolating %s (from ctd_time)...\n', var_name);
    
    % Get dimensions
    [n_points_fet, n_profiles] = size(s_raw.fet_time);
    s.(var_name) = NaN(n_points_fet, n_profiles);
    
    % Interpolate each profile individually
    for i = 1:n_profiles
        % Get valid (non-NaN) data for this profile
        valid_ctd = ~isnan(s_raw.ctd_time(:,i)) & ~isnan(s_raw.(var_name)(:,i));
        valid_fet = ~isnan(s_raw.fet_time(:,i));
        
        if sum(valid_ctd) > 1 && sum(valid_fet) > 0
            try
                % Interpolate using linear interpolation
                s.(var_name)(valid_fet, i) = interp1(...
                    s_raw.ctd_time(valid_ctd, i), ...
                    s_raw.(var_name)(valid_ctd, i), ...
                    s_raw.fet_time(valid_fet, i), ...
                    'linear', NaN);
            catch ME
                warning('Could not interpolate %s for profile %d: %s', var_name, i, ME.message);
            end
        end
    end
end
% % % % Interpolate CTD variables from ctd_time to fet_time
% % % for v = 1:length(ctd_qc_vars)
% % %     var_name = ctd_qc_vars{v};
% % %     fprintf('  Interpolating %s (from ctd_time)...\n', var_name);
% % % 
% % %     % Get dimensions
% % %     [n_points_fet, n_profiles] = size(s_raw.fet_time);
% % %     s.(var_name) = NaN(n_points_fet, n_profiles);
% % % 
% % %     % Interpolate each profile individually
% % %     for i = 1:n_profiles
% % %         % Get valid (non-NaN) data for this profile
% % %         valid_ctd = ~isnan(s_raw.ctd_time(:,i)) & ~isnan(s_raw.(var_name)(:,i));
% % %         valid_fet = ~isnan(s_raw.fet_time(:,i));
% % % 
% % %         if sum(valid_ctd) > 1 && sum(valid_fet) > 0
% % %             try
% % %                 % Interpolate using linear interpolation
% % %                 s.(var_name)(valid_fet, i) = interp1(...
% % %                     s_raw.ctd_time(valid_ctd, i), ...
% % %                     s_raw.(var_name)(valid_ctd, i), ...
% % %                     s_raw.fet_time(valid_fet, i), ...
% % %                     'linear', NaN);
% % %             catch ME
% % %                 warning('Could not interpolate %s for profile %d: %s', var_name, i, ME.message);
% % %             end
% % %         end
% % %     end
% % % end
% Handle dissolved oxygen separately (on its own time grid)
fprintf('  Interpolating doxy (from doxy_time to fet_time)...\n');
[n_points_fet, n_profiles] = size(s_raw.fet_time);
s.doxy = NaN(n_points_fet, n_profiles);

% Might not have to do this. (Interp CTD onto fet time grid to compute pH
% then put everything back onto the CTD time grid)
% For oxygen, interpolate from doxy_time to fet_time
for i = 1:n_profiles
    % Get valid data for oxygen on doxy_time grid
    valid_doxy = ~isnan(s_raw.doxy_time(:, i)) & ~isnan(s_raw.doxy(:, i));
    valid_fet = ~isnan(s_raw.fet_time(:, i));
    
    if sum(valid_doxy) > 1 && sum(valid_fet) > 0
        try
            % Interpolate oxygen from doxy_time to fet_time
            s.doxy(valid_fet, i) = interp1(...
                s_raw.doxy_time(valid_doxy, i), ...
                s_raw.doxy(valid_doxy, i), ...
                s_raw.fet_time(valid_fet, i), ...
                'linear', NaN);
        catch ME
            warning('Could not interpolate doxy for profile %d: %s', i, ME.message);
        end
    end
end

fprintf('\nInterpolation complete! All variables now on fet_time grid.\n');
fprintf('Interpolated data stored in structure: s_interp\n');
%% Split upcast/downcast
s.upcast = s.fet_time > s_raw.start_of_climb_time; % Make sure this is the correct time
%% Add other vars
s.unixtime = log_gps_time(1,:); % verify this top row is correct
s.sdn = datenum(datetime(s_raw.unixtime, 'ConvertFrom', 'posixtime'));
s.lat = log_gps_lat(1,:);
s.lon = log_gps_lon(1,:);
s.divenum = log_dive;
s.missionID = missionID;
%% Get pH cal and calc pH
try
    missionID = '25B68601';
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
    Pcoefs = [pH_cal.fp_k1,pH_cal.fp_k2,pH_cal.fp_k3,...
                pH_cal.fp_k4,pH_cal.fp_k5,pH_cal.fp_k6]';
    [~, s.pHin] = phcalc_jp(s.vrse,s.pres,s.tc,s.psal,pH_cal.k0_seawater,pH_cal.k2_fp_c0,Pcoefs); % BW 9/25/2025 change to pH_total
catch
    pH_cal = readtable('parsing\ph_cal.csv'); % local path
    pH_cal = pH_cal(end,:);
    Pcoefs = [pH_cal.fp_k1,pH_cal.fp_k2,pH_cal.fp_k3,...
                pH_cal.fp_k4,pH_cal.fp_k5,pH_cal.fp_k6]';
    [~, s.pHin] = phcalc_jp(s.vrse,s.pres,s.tc,s.psal,pH_cal.k0,pH_cal.k2,Pcoefs); % BW 9/25/2025 change to pH_total
end
%% carb
s = calc_carb(s);
%% Quick table maker
[n, m] = size(s.pHin);
nm = n*m;
sdn = repmat(s.sdn,n,1);
% sdnu = repmat(s.sdnu,n,1);
unixtime = repmat(s.unixtime,n,1);
% unixtimeu = repmat(s.unixtimeu,n,1);
lat = repmat(s.lat,n,1);
% latu = repmat(s.latu,n,1);
lon = repmat(s.lon,n,1);
% lonu = repmat(s.lonu,n,1);
% u = repmat(s.u,n,1);
% v = repmat(s.v,n,1);
divenum = repmat(1:m,n,1);
T = table();
T.mission_id = string(repmat(s.missionID,length(sdn(:)),1));
T.divenumber = int64(divenum(:));
T.sdn = sdn(:);
T.unixtime = unixtime(:);
T.lat = lat(:);
T.lon = lon(:);
% T.sdnu = sdnu(:);
% T.unixtimeu = unixtimeu(:);
% T.latu = latu(:);
% T.lonu = lonu(:);
T.depth = s.depth(:);
T.pres = s.pres(:);
T.tc = s.tc(:);
T.psal = s.psal(:);
% T.chla = s.chla(:);
T.doxy = s.doxy(:);
T.phin = s.pHin(:);
% T.theta = s.theta(:);
% T.sigma = s.sigma(:);
% T.rho = s.rho(:);
% Carbonate parameters
fields = {'pH25atm', 'pHin_canb', 'dic_canb','ta_canb','po4_canb',...
    'sil_canb','pHin_canbtadic','pco2in','fco2in','satarin','satarin'};
for i = 1:length(fields)
    field = fields{i};
    if isfield(s, field)
        T.(lower(field)) = s.(field)(:);
    else
        T.(lower(field)) = NaN(nm, 1);
    end
end
% T.u = u(:);
% T.v = v(:);
% % T.udop = s.udop(:);
% T.vdop = s.vdop(:);
% T.abs = s.abs(:);
% pH diagnostics
T.vrse = s.vrse(:);
T.vrse_std = s.vrse_std(:);
T.vk = s.vk(:);
T.vk_std = s.vk_std(:);
T.ik = s.ik(:);
T.ib = s.ib(:);
T.castdir = s.upcast(:); % 0: downcast, 1: upcast
%%
iup = T.castdir == 1;
T_up = T(iup,:);
%% Save file
% fname = ['\\atlas.shore.mbari.org\ProjectLibrary\901805_Coastal_Biogeochemical_Sensing\Spray_Data\',missionID,'\sat_files\',[missionID,'sat.mat']];
% save(fname,"s");
% %% Push database
% inan = isnan(T.depth);
% T(inan,:) = []; % remove all rows where depth is nan
% %%
% missionID = '25B68601';
% databaseName = 'glidata';
% username = 'glidata_user';
% password = 'IoFUTBeaQDppSYcmBebA4rV8SJOEMCFI';
% server = 'dpg-d2jobg3e5dus738ce5vg-a.oregon-postgres.render.com';
% port = '5432';
% conn = connect_glidata(databaseName, ...
%                 username,password,server,port);
% 
% query = sprintf([...
%         'SELECT mission_id, divenumber, depth ' ...
%         'FROM public.real_time_binned_id_%s ' ...
%         'ORDER BY mission_id ASC, divenumber ASC, depth ASC'], ...
%         missionID);
%     existingKeys = fetch(conn, query);
%     % Create logical index of existing rows
%     if ~isempty(existingKeys)
%         % Compare using setdiff with 'rows'
%         [~, newRowIdx] = setdiff(...
%             [T.mission_id, T.divenumber, T.depth], ...
%             [string(existingKeys.mission_id), ...
%             existingKeys.divenumber, existingKeys.depth], ...
%             'rows');
%         newRows = T(newRowIdx, :);
%     else
%         newRows = T;
%     end
%     % Only upload if there are new rows
%     if height(newRows) > 0
%         % Write to the partitioned table
%         tableName = sprintf('real_time_binned_id_%s',...
%             string(missionID));
%         sqlwrite(conn, tableName, newRows);
%         fprintf('Uploaded %d new rows successfully\n',...
%             height(newRows));
%     else
%         disp('No new rows to upload');
%     end
% 
%     if isopen(conn)
%         close(conn)
%     end