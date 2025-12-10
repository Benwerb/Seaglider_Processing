% dir_path = '\\atlas\ftp\incoming\Fred\';
dir_path = 'C:\Users\bwerb\Documents\GitHub\Seaglider_Processing\data';
missionID = '25B68601';
SN = 686;
basepath = ['\\atlas.shore.mbari.org\ProjectLibrary\' ...
            '901805_Coastal_Biogeochemical_Sensing\Spray_Data'];
missionPath = fullfile(basepath, missionID);
matPath = fullfile(missionPath,'real_time_binned',...
                [missionID,'sat.mat']);

[s1, data] = extractSeagliderData(dir_path, missionID);
%%
bindata = bin_sp2(data,0,5,1000,'d','ascend',1,'descend',0,'both',0);
s = buildSpray2Struct(bindata);
s = calc_carb(s);
s.missionID = missionID;
s.SN = SN;
T = buildSpray2Table(s);
% save(matPath,"s",'-mat')



profileFigure(T)
%%
for prof = 5
    ctd_p = data.ctd.p{prof,:}; 
    ctd_t = data.ctd.time{prof,:};
    ctd_phase = data.ctd.phase{prof,:};

    ph = data.ph.ph{prof,:};
    ph_p = data.ph.p{prof,:};
    ph_phase = data.ph.phase{prof,:};
   
    %
    figure(1)
    clf
    subplot(121)
    % plot(t_test,p_test)
    hold on
    plot(ctd_t(ctd_phase),ctd_p(ctd_phase),'Color','red')
    plot(ctd_t(~ctd_phase),ctd_p(~ctd_phase),'Color','black')
    xline(data.ctd.start_of_climb_time{prof},'--')
    set(gca,"YDir","reverse")
    pause(.1)
    legend('descent phase','ascent phase','start of climb time',Location='southeast')
    xlabel('seconds since dive start')
    ylabel('pressure')
    title(['Profile ', num2str(data.ph.ndive(prof))])

    subplot(122)
    hold on
    plot(ph(ph_phase),ph_p(ph_phase),'Color','red')
    plot(ph(~ph_phase),ph_p(~ph_phase),'Color','black')
    set(gca,"YDir","reverse")
    % ylabel('pressure')
    xlabel('pH')
    % legend('descent phase','ascent phase',Location='best')
end
%%
g = Spray2Data('25B68601');
g.ascend = 0;
g.descend = 1; % something is reversed here. Seems like descent is ascent
% but also I think they turned on ascent only so the last 2 profiles are
% marked as all 0's for ascend
g.forceProcess = 1;
g.savefigure = 1;
g.sendEmails = 1;
g.updateODSS = 1;
g.recalccarb = 0;
g.applyphcorr = 0;
g.process_full_pipeline;
%%
g = Spray2Data('25B20901');
g.forceProcess = 1;
g.savefigure = 1;
g.sendEmails = 1;
g.updateODSS = 1;
g.recalccarb = 1;
g.applyphcorr = 1;
g.process_full_pipeline;