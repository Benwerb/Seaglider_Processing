dir_path = '\\atlas\ftp\incoming\Fred\';
missionID = '25B68601';

[s, data] = extractSeagliderData(dir_path, missionID);
bindata = bin_sp2(data,0,5,1000,'p','ascend',1,'descend',0);

bindata.dt = datetime(bindata.time, 'ConvertFrom', 'posixtime');
bindata.dn = datenum(bindata.dt);
pgrid = repmat(bindata.p,1,length(bindata.time));
timegrid = repmat(bindata.dn',length(bindata.p),1);
phgrid = bindata.ph;
tgrid = bindata.t;
sgrid = bindata.s;
fig = figure();

tl = tiledlayout(3,1,"TileSpacing","tight");
title(tl,'Trinidad Head Line OSU')
subtitle(tl,'Seaglider 686')

nexttile
contourf(timegrid,pgrid,tgrid)
datetick('x',2)
set(gca,'Ydir','rev')
colorbar
title('Temperature [^oC]')

nexttile
contourf(timegrid,pgrid,sgrid)
datetick('x',2)
set(gca,'Ydir','rev')
colorbar
title('Salinity [PSU]')


nexttile
contourf(timegrid,pgrid,phgrid)
datetick('x',2)
set(gca,'Ydir','rev')
colorbar
title('pHin: GliderFET009')