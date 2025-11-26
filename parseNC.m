% Directory containing NetCDF files
dir_path = '\\atlas\ftp\incoming\Fred\';

% ========== CONFIGURATION: Add variables to process here ==========
variables_to_process = {'log_GPS', 'log_gps_time', 'log_gps_lat', 'log_gps_lon', 'log_dive', 'fet_time','fet_Vrse', 'fet_Vrsestd', 'fet_Vk', 'fet_Vkstd', 'fet_Ik', 'fet_Ib', 'ctd_pressure', 'ctd_depth', 'ctd_time', 'temperature', 'salinity', 'aanderaa4831_dissolved_oxygen','aanderaa4831_results_time'};
% ==================================================================

% Get list of all p686*.nc files
file_pattern = fullfile(dir_path, 'p686*.nc');
files = dir(file_pattern);

% Initialize storage for all variables
variable_data = struct();
for v = 1:length(variables_to_process)
    var_name = variables_to_process{v};
    variable_data.(var_name).max_length = 0;
    variable_data.(var_name).profile_data = {};
    variable_data.(var_name).profile_numbers = [];
end

% First pass: determine maximum array length for each variable
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

% Display max lengths
for v = 1:length(variables_to_process)
    var_name = variables_to_process{v};
    fprintf('  %s: maximum length = %d\n', var_name, variable_data.(var_name).max_length);
end

% Second pass: extract data and pad to max_length for each variable
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
            max_length = variable_data.(var_name).max_length;
            
            try
                % Read the data
                data = ncread(file_path, var_name);
                
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
s = struct;
s.tc = temperature;
s.psal = salinity;
s.ctd_time = ctd_time;
s.depth = ctd_depth;
s.pres = ctd_pressure;
s.fet_time = fet_time;
s.vrse = fet_Vrse;
s.vrse_std = fet_Vrsestd;
s.vk = fet_Vk;
s.vk_std = fet_Vkstd;
s.ik = fet_Ik;
s.ib = fet_Ib;
s.doxy = aanderaa4831_dissolved_oxygen;
s.doxy_time = aanderaa4831_results_time;
%% Interpolate all variables onto fet_time grid
fprintf('\nInterpolating variables onto fet_time grid...\n');

% List of CTD variables to interpolate (on ctd_time grid)
ctd_vars = {'tc', 'psal', 'depth', 'pres'};

% Initialize interpolated structure
s_interp = struct;

% Copy FET variables (already on fet_time grid)
s_interp.fet_time = s.fet_time;
s_interp.vrse = s.vrse;
s_interp.vrse_std = s.vrse_std;
s_interp.vk = s.vk;
s_interp.vk_std = s.vk_std;
s_interp.ik = s.ik;
s_interp.ib = s.ib;

% Interpolate CTD variables from ctd_time to fet_time
for v = 1:length(ctd_vars)
    var_name = ctd_vars{v};
    fprintf('  Interpolating %s (from ctd_time)...\n', var_name);
    
    % Get dimensions
    [n_points_fet, n_profiles] = size(s.fet_time);
    s_interp.(var_name) = NaN(n_points_fet, n_profiles);
    
    % Interpolate each profile individually
    for i = 1:n_profiles
        % Get valid (non-NaN) data for this profile
        valid_ctd = ~isnan(s.ctd_time(:,i)) & ~isnan(s.(var_name)(:,i));
        valid_fet = ~isnan(s.fet_time(:,i));
        
        if sum(valid_ctd) > 1 && sum(valid_fet) > 0
            try
                % Interpolate using linear interpolation
                s_interp.(var_name)(valid_fet, i) = interp1(...
                    s.ctd_time(valid_ctd, i), ...
                    s.(var_name)(valid_ctd, i), ...
                    s.fet_time(valid_fet, i), ...
                    'linear', NaN);
            catch ME
                warning('Could not interpolate %s for profile %d: %s', var_name, i, ME.message);
            end
        end
    end
end

% Handle dissolved oxygen separately (on its own time grid)
fprintf('  Interpolating doxy (from doxy_time to fet_time)...\n');
[n_points_fet, n_profiles] = size(s.fet_time);
s_interp.doxy = NaN(n_points_fet, n_profiles);

% For oxygen, interpolate from doxy_time to fet_time
for i = 1:n_profiles
    % Get valid data for oxygen on doxy_time grid
    valid_doxy = ~isnan(s.doxy_time(:, i)) & ~isnan(s.doxy(:, i));
    valid_fet = ~isnan(s.fet_time(:, i));
    
    if sum(valid_doxy) > 1 && sum(valid_fet) > 0
        try
            % Interpolate oxygen from doxy_time to fet_time
            s_interp.doxy(valid_fet, i) = interp1(...
                s.doxy_time(valid_doxy, i), ...
                s.doxy(valid_doxy, i), ...
                s.fet_time(valid_fet, i), ...
                'linear', NaN);
        catch ME
            warning('Could not interpolate doxy for profile %d: %s', i, ME.message);
        end
    end
end

fprintf('\nInterpolation complete! All variables now on fet_time grid.\n');
fprintf('Interpolated data stored in structure: s_interp\n');
%% Get pH cal and calc pH
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

[~, s_interp.pHin] = phcalc_jp(s_interp.vrse,s_interp.pres,s_interp.tc,s_interp.psal,pH_cal.k0_seawater,pH_cal.k2_fp_c0,Pcoefs); % BW 9/25/2025 change to pH_total
%% s_interp -> s
s = s_interp;
s.unixtime = log_gps_time(1,:); % verify this top row is correct
s.sdn = datenum(datetime(s.unixtime, 'ConvertFrom', 'posixtime'));
s.lat = log_gps_lat(1,:);
s.lon = log_gps_lon(1,:);
s.divenum = log_dive;
s.missionID = missionID;
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
%% Save file
fname = ['\\atlas.shore.mbari.org\ProjectLibrary\901805_Coastal_Biogeochemical_Sensing\Spray_Data\',missionID,'\sat_files\',[missionID,'sat.mat']];
save(fname,"s");
%% Push database
inan = isnan(T.depth);
T(inan,:) = []; % remove all rows where depth is nan
%%
missionID = '25B68601';
databaseName = 'glidata';
username = 'glidata_user';
password = 'IoFUTBeaQDppSYcmBebA4rV8SJOEMCFI';
server = 'dpg-d2jobg3e5dus738ce5vg-a.oregon-postgres.render.com';
port = '5432';
conn = connect_glidata(databaseName, ...
                username,password,server,port);

query = sprintf([...
        'SELECT mission_id, divenumber, depth ' ...
        'FROM public.real_time_binned_id_%s ' ...
        'ORDER BY mission_id ASC, divenumber ASC, depth ASC'], ...
        missionID);
    existingKeys = fetch(conn, query);
    % Create logical index of existing rows
    if ~isempty(existingKeys)
        % Compare using setdiff with 'rows'
        [~, newRowIdx] = setdiff(...
            [T.mission_id, T.divenumber, T.depth], ...
            [string(existingKeys.mission_id), ...
            existingKeys.divenumber, existingKeys.depth], ...
            'rows');
        newRows = T(newRowIdx, :);
    else
        newRows = T;
    end
    % Only upload if there are new rows
    if height(newRows) > 0
        % Write to the partitioned table
        tableName = sprintf('real_time_binned_id_%s',...
            string(missionID));
        sqlwrite(conn, tableName, newRows);
        fprintf('Uploaded %d new rows successfully\n',...
            height(newRows));
    else
        disp('No new rows to upload');
    end

    if isopen(conn)
        close(conn)
    end