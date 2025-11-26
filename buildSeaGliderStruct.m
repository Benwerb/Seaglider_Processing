function s = buildSeaGliderStruct(data)
%BUILDSPRAY2STRUCT Create s struct from the bindata struct
[n, m] = size(data.pHin);
% Build s
s = struct;
s.unixtime = data.time';
s.unixtimeu = data.timeu';
s.sdn = datenum(datetime(data.time', 'ConvertFrom', 'posixtime'));
s.sdnu = datenum(datetime(data.timeu', 'ConvertFrom', 'posixtime'));
s.lat = data.lat';
s.lon = data.lon';
s.latu = data.latu';
s.lonu = data.lonu';
s.u = data.u';
s.v = data.v';
s.unixtime_ = [data.time'; data.timeu'];
s.sdn_ = [s.sdn; s.sdnu];
s.lat_ = [data.lat'; data.latu'];
s.lon_ = [data.lon'; data.lonu'];
s.depth = repmat(data.depth,1,m);
s.pres = sw_pres(s.depth,s.lat);

% Assign fields with NaN fallback
s.tc = getFieldOrNaN(data, 't', n, m);
s.psal = getFieldOrNaN(data, 's', n, m);
s.chla = getFieldOrNaN(data, 'fl', n, m);
s.doxy = getFieldOrNaN(data, 'oxumolkg', n, m);
s.oxconc = getFieldOrNaN(data, 'oxconc', n, m);
s.ox = getFieldOrNaN(data, 'ox', n, m);
s.theta = getFieldOrNaN(data, 'theta', n, m);
s.sigma = getFieldOrNaN(data, 'sigma', n, m);
s.rho = getFieldOrNaN(data, 'rho', n, m);
s.vrse = getFieldOrNaN(data, 'Vrse', n, m);
s.vrse_std = getFieldOrNaN(data, 'Vrse_std', n, m);
s.vk = getFieldOrNaN(data, 'Vk', n, m);
s.vk_std = getFieldOrNaN(data, 'Vk_std', n, m);
s.ik = getFieldOrNaN(data, 'Ik', n, m);
s.ib = getFieldOrNaN(data, 'Ib', n, m);
s.pHin = getFieldOrNaN(data, 'pHin', n, m);
s.udop = getFieldOrNaN(data, 'udop', n, m);
s.vdop = getFieldOrNaN(data, 'vdop', n, m);
s.abs = getFieldOrNaN(data, 'abs', n, m);
end

function value = getFieldOrNaN(struct_data, fieldname, n, m)
    % Return field value if it exists and is non-empty, otherwise return NaN array
    if isfield(struct_data, fieldname) && ~isempty(struct_data.(fieldname))
        value = struct_data.(fieldname);
    else
        value = nan(n, m);
    end
end