%load the Eco and CDOm data
load('/home/flagg/Dropbox (MIT)/PikaL_Dev/chl_data_mvco/code/code/ecochlsmooth2018.mat')
load('/home/flagg/Dropbox (MIT)/PikaL_Dev/chl_data_mvco/code/code/CDOMsmooth2018.mat')
load('/home/flagg/Dropbox (MIT)/PikaL_Dev/downloadedDatasets/MetDat_s-2.mat')


startNumber = 20000;
startNumberCDOM = 1; %212 and 213 are trouble, 2467-2470

chl_tss_time = ecosmooth(startNumber:end,1);
cdom_time = CDOMsmooth(startNumberCDOM:end,1);
cleaned_CDOM_time = cdom_time
cleaned_CDOM_time(2467:2468,:) = cleaned_CDOM_time(2466,:)
cleaned_CDOM_time(2469:2470,:) = cleaned_CDOM_time(2471,:)
cdom_time = cleaned_CDOM_time

close all; subplot(414); plot((chl_tss_time),ecosmooth(startNumber:end,7)); datetickzoom('x','mmm dd hh','keepticks'); ylabel('ug/l Chl a'); grid on;
subplot(411); plot(chl_tss_time,ecosmooth(startNumber:end,8)); title('\beta_7_0_0 (1/m)'); datetickzoom('x','mmm dd hh','keepticks'); ylabel('\beta_7_0_0 (1/m)'); grid on;
subplot(412); scatter(cdom_time,CDOMsmooth(startNumberCDOM:end,4),'.'); title('CDOM'); datetickzoom('x','mmm dd hh','keepticks'); ylabel('CDOM'); grid on;

%get the date data from the MVCO data
mvco_range_start = 432320
mvco_range=mvco_range_start:size(Yday,2);
for(i= mvco_range)
    date_number_mvco(i-432319) = datenum(year(i),MO(i),day(i),HH(i),MM(i),SS(i));
end
% plot wind speed data
subplot(413); scatter(date_number_mvco,wspd_son3d_mean(mvco_range),'.'); title('Mean Windspeed'); datetickzoom('x','mmm dd hh','keepticks'); ylabel('Windspeed m/s'); grid on;
%subplot(413); plot(date_number_mvco,wspd_son3d_std(mvco_range)); title('Windspeed STD'); datetickzoom('x','mmm dd hh','keepticks'); ylabel('Windspeed STD m/s'); grid on;

figure; subplot(411); plot(date_number_mvco); subplot(412); plot( HH(mvco_range))
%%
X=load('2018_MetDat_s.C99','ascii');

yd_met=X(:,2);
Year = X(:,1);
Month = X(:,3);
Day = X(:,4);
Hour = X(:,5);
Minute = X(:,6);
Second = X(:,7);

Wspd_son3D_mean=X(:,8);

for(i= 1:length(Year))
    date_number_mvco(i) = datenum(Year(i),Month(i),Day(i),Hour(i),Minute(i),Second(i));
end
% Wdir_son3D_mean=X(:,9);
% Tair_vai_mean=X(:,10);
% RH_vai_mean=X(:,11);
% Pressure_vai_mean=X(:,12);
% IR_campmt_median=X(:,13);
% Solar_campmt_median=X(:,14);
% Rain_campmt=X(:,15);
% Tson_son3D_mean=X(:,16);
% 
% Wspd_son3D_std=X(:,17);
% v_son3D_std=X(:,18);
% w_son3D_std=X(:,19);
% 
% uv_son3D_var=X(:,20);
% uw_son3D_var=X(:,21);
% vw_son3D_var=X(:,22);
% 
% wTv_son3D_var=X(:,23);
% wTv_son3DVaiPTU_var=X(:,24);
% wT_son3DVaiPTU_var=X(:,25);
% wrhov_son3DVaiPTULicor_var=X(:,26);
% wco2_son3DVaiPTULicor_var=X(:,27);
% 
% lrec_3d = X(:,26);
% lrec_vai=X(:,29);
% lrec_li=X(:,30);
% lrec_camp=X(:,31);

clear X ans n2

%%
%pull out data points for each of the given matlab times. 
% given a column of times, it provides a column of the data sections that
% match those times from the MVCO dataset (OG times are from the MVCO
% dataset) 
%close all
clear mvco_matched_camera_times_TSS
load('surfaceAndDepthData.mat')

cameraTimes = dates';

%plot(date_number_mvco); hold on; plot(cameraTimes); legend('date_num_mvco','camera_times')
for(i = 1:4612)%length(cameraTimes))
    [min_val,mvco_matched_camera_times_BETA_depth(i)] = min(abs(depth_datenums_b'-cameraTimes(i)));
    [min_val,mvco_matched_camera_times_CHL_depth(i)] = min(abs(depth_datenums_chl'-cameraTimes(i)));
    [min_val,mvco_matched_camera_times_CDOM_depth(i)] = min(abs(depth_datenums_cdom'-cameraTimes(i)));
    
    
    [min_val,mvco_matched_camera_times_BETA_surface(i)] = min(abs(surface_datenums'-cameraTimes(i)));
    [min_val,mvco_matched_camera_times_CHL_surface(i)] = min(abs(surface_datenums'-cameraTimes(i)));
    [min_val,mvco_matched_camera_times_CDOM_surface(i)] = min(abs(surface_datenums'-cameraTimes(i)));
    [min_val,mvco_matched_camera_times_wind(i)] = min(abs(date_number_mvco'-cameraTimes(i)));

end



%TSS_values = depth_b( mvco_matched_camera_times_TSS_CHL_BETA);
%TSS_times = depth_datenums_chl(mvco_matched_camera_times_TSS_CHL_BETA);
b_values_depth = depth_b( mvco_matched_camera_times_BETA_depth);
b_times_depth = depth_datenums_chl(mvco_matched_camera_times_BETA_depth);
chl_values_depth = depth_chl( mvco_matched_camera_times_CHL_depth);
chl_times_depth = depth_datenums_chl(mvco_matched_camera_times_CHL_depth);
CDOM_values_depth = depth_cdom( mvco_matched_camera_times_CDOM_depth);
CDOM_times_depth = depth_datenums_cdom(mvco_matched_camera_times_CDOM_depth);

b_values_surface = surface_median_backscatter( mvco_matched_camera_times_BETA_surface)';
b_times_surface = surface_datenums(mvco_matched_camera_times_BETA_surface);
chl_values_surface = surface_median_chl( mvco_matched_camera_times_CHL_surface)';
chl_times_surface = surface_datenums(mvco_matched_camera_times_CHL_surface);
CDOM_values_surface = surface_median_cdom( mvco_matched_camera_times_CDOM_surface)';
CDOM_times_surface = surface_datenums(mvco_matched_camera_times_CDOM_surface);
wind_values = Wspd_son3D_mean(mvco_matched_camera_times_wind)
wind_times = date_number_mvco(mvco_matched_camera_times_wind)




chl_values_manual = [5.009 6.352 5.814 5.627 8.379 9.825 11.02 11.81 7.847 12.25 12.31 11.88 12.37 10.39 4.749 5.288 4.586 5.808 9.789 8.99 6.437 4.556 3.588 11.82 13.79 9.196 11.05 21.96 9.965 7.079 8.567 10.19 5.548 3.134 11.89 8.022 8.561 11.54 6.104 3.86 2.529 2.299 1.785 1.537 1.119 1.022 .8954 .7986 .7683 .8651 .6171 .8046 .7865 .8772 1.089 1.047 1.972 1.24 1.24 1.258 1.016 1.059 1.077 .9801 1.174 1.252 1.331 1.307 1.283 1.004 .7865 .9377 .7502 1.035 1.228 1.174 1.041 1.107 .8409 .7441   ]
chl_values_start_day = datenum(2018,03,23,08,00,00)
chl_times_manual = (0:size(chl_values_manual,2)-1) + chl_values_start_day
interped_chl_values = interp1(chl_times_manual,chl_values_manual,chl_times_depth)

figure
plot(chl_times_depth, interped_chl_values); hold on; plot(chl_times_depth, chl_values_depth);
hold on
plot(chl_times_surface, chl_values_surface); 
title('Chl'); 
 datetickzoom('x', 'HH:MM mm/dd ');
 
figure
plot(b_times_depth, b_values_depth);
hold on
plot(b_times_surface, b_values_surface); 
title('backscatter'); 
 datetickzoom('x', 'HH:MM mm/dd ');
 
figure
plot(CDOM_times_depth, CDOM_values_depth); 
hold on
plot(CDOM_times_surface, CDOM_values_surface); 
title('CDOM'); 
 datetickzoom('x', 'HH:MM mm/dd ');
 
figure
plot(wind_times, wind_values); 
title('wind'); 
 datetickzoom('x', 'HH:MM mm/dd ');

save MVCO_DATA_VALUES_surface_depth_camera_aligned.mat b_values_depth b_times_depth chl_values_depth chl_times_depth CDOM_values_depth CDOM_times_depth b_values_surface b_times_surface chl_values_surface chl_times_surface CDOM_values_surface CDOM_times_surface wind_values wind_times
% start = 23 
%chl values = 5.009 6.352 5.814 5.627 8.379 9.825 11.02 11.81 7.847 12.25
%12.31 11.88 12.37 10.39 4.749 5.288 4.586 5.808 9.789 8.99 6.437 4.556
%3.588 11.82 13.79 9.196 11.05 21.96 9.965 7.079 8.567 10.19 5.548 3.134
% 11.89 8.022 8.561 11.54 6.104 3.86 2.529 2.299 