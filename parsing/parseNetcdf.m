function [data, s] = parseNetcdf(directory,seagliderSN)
%PARSE Summary of this function goes here
%   INPUTS: directory containing the files, seagliderSN: glider SN 
% OUTPUT: data, a struct with all parameters organized into sensor structs.
% The format matches the processing functions for Spray2 so that data can
% be processed and displayed with the method.

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
    
    s = extractVariables(variables_to_process,directory,seagliderSN);
    data = [];
end

