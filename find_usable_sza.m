function [vect_sza_approved] = find_usable_sza(datetime_sza_est,min_angle,max_angle)     
camera_rotational_angle = 320; %clockwise from north
                DST = false; % daylight savings time has been removed

        % calculate solar position
        lat = 41.325; % [arc-degrees] latitude
        long = -70.5667; % [arc-degrees] longitude
        TZ = -4; % [hrs] offset from UTC, during standard time
        rot = 0; % [arc-degrees] rotation clockwise from north
        [angles ,projection] = solarPosition(datetime_sza_est,lat,long,TZ,rot,DST); % [arc-degrees], time in UTC
        solar_azimuth_angle = abs(angles(:,2)-camera_rotational_angle);
       % figure
       % plot(datetime_sza_est,solar_azimuth_angle.*(solar_azimuth_angle>115).*(solar_azimuth_angle<155))
       % datetickzoom
        vect_sza_approved = (solar_azimuth_angle>min_angle).*(solar_azimuth_angle<max_angle);
        
end