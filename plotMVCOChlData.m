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
%pull out data points for each of the given matlab times. 
% given a column of times, it provides a column of the data sections that
% match those times from the MVCO dataset (OG times are from the MVCO
% dataset) 
%close all
clear mvco_matched_camera_times_TSS
figure
cameraTimes = dates';

plot(date_number_mvco); hold on; plot(cameraTimes); legend('date_num_mvco','camera_times')
for(i = 1:2200)%length(cameraTimes))
    [min_val,mvco_matched_camera_times_TSS_CHL_BETA(i)] = min(abs(chl_tss_time'-cameraTimes(i)));
    [min_val,mvco_matched_camera_times_CDOM(i)] = min(abs(cdom_time'-cameraTimes(i)));
    [min_val,mvco_matched_camera_times_wind(i)] = min(abs(date_number_mvco'-cameraTimes(i)));

end
figure
plot(mvco_matched_camera_times_TSS_CHL_BETA)
%mvco_matched_camera_times_TSS=mvco_matched_camera_times_TSS+startNumber;
figure


TSS_values = ecosmooth( mvco_matched_camera_times_TSS_CHL_BETA+startNumber,5);
TSS_times = chl_tss_time(mvco_matched_camera_times_TSS_CHL_BETA);
b_values = ecosmooth( mvco_matched_camera_times_TSS_CHL_BETA+startNumber,8);
b_times = chl_tss_time(mvco_matched_camera_times_TSS_CHL_BETA);
chl_values = ecosmooth( mvco_matched_camera_times_TSS_CHL_BETA+startNumber,7);
chl_times = chl_tss_time(mvco_matched_camera_times_TSS_CHL_BETA);
CDOM_values = CDOMsmooth( mvco_matched_camera_times_CDOM+startNumberCDOM,4);
CDOM_times = chl_tss_time(mvco_matched_camera_times_CDOM);
wind_values = wspd_son3d_mean(mvco_matched_camera_times_wind + mvco_range_start)
wind_times = date_number_mvco(mvco_matched_camera_times_wind)


plot(chl_times,chl_values); datetickzoom('x','mmm dd hh','keepticks');



chl_values_manual = [5.009 6.352 5.814 5.627 8.379 9.825 11.02 11.81 7.847 12.25 12.31 11.88 12.37 10.39 4.749 5.288 4.586 5.808 9.789 8.99 6.437 4.556 3.588 11.82 13.79 9.196 11.05 21.96 9.965 7.079 8.567 10.19 5.548 3.134 11.89 8.022 8.561 11.54 6.104 3.86 2.529 2.299]
chl_values_start_day = datenum(2018,03,23,08,00,00)
chl_times_manual = (0:size(chl_values_manual,2)-1) + chl_values_start_day
interped_chl_values = interp1(chl_times_manual,chl_values_manual,chl_times)

figure
plot(chl_times, interped_chl_values); hold on; plot(chl_times, chl_values); 

%save MVCO_DATA_VALUES.mat TSS_values TSS_times b_values b_times chl_values chl_times CDOM_values CDOM_times  interped_chl_values wind_values wind_times
% start = 23 
%chl values = 5.009 6.352 5.814 5.627 8.379 9.825 11.02 11.81 7.847 12.25
%12.31 11.88 12.37 10.39 4.749 5.288 4.586 5.808 9.789 8.99 6.437 4.556
%3.588 11.82 13.79 9.196 11.05 21.96 9.965 7.079 8.567 10.19 5.548 3.134
% 11.89 8.022 8.561 11.54 6.104 3.86 2.529 2.299 