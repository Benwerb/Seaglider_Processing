function T = buildSpray2Table(s)
%BUILDSPRAY2TABLE Creates a table from the s struct that can be uploaded to
%the spray database
[n, m] = size(s.pHin);
nm = n*m;
sdn = repmat(s.sdn,n,1);
sdnu = repmat(s.sdnu,n,1);
unixtime = repmat(s.unixtime,n,1);
unixtimeu = repmat(s.unixtimeu,n,1);
lat = repmat(s.lat,n,1);
latu = repmat(s.latu,n,1);
lon = repmat(s.lon,n,1);
lonu = repmat(s.lonu,n,1);
u = repmat(s.u,n,1);
v = repmat(s.v,n,1);
divenum = repmat(1:m,n,1);
T = table();
T.mission_id = string(repmat(s.missionID,length(sdn(:)),1));
T.divenumber = int64(divenum(:));
T.sdn = sdn(:);
T.unixtime = unixtime(:);
T.lat = lat(:);
T.lon = lon(:);
T.sdnu = sdnu(:);
T.unixtimeu = unixtimeu(:);
T.latu = latu(:);
T.lonu = lonu(:);
T.depth = s.depth(:);
T.pres = s.pres(:);
T.tc = s.tc(:);
T.psal = s.psal(:);
T.chla = s.chla(:);
T.doxy = s.doxy(:);
T.phin = s.pHin(:);
T.theta = s.theta(:);
T.sigma = s.sigma(:);
T.rho = s.rho(:);
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
T.u = u(:);
T.v = v(:);
T.udop = s.udop(:);
T.vdop = s.vdop(:);
T.abs = s.abs(:);
% pH diagnostics
T.vrse = s.vrse(:);
T.vrse_std = s.vrse_std(:);
T.vk = s.vk(:);
T.vk_std = s.vk_std(:);
T.ik = s.ik(:);
T.ib = s.ib(:);
end

