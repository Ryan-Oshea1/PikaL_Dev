%THis script loads data from the Eco triplet, and plots it
clc
close all
clear all
load('mvco_triplet_data.mat')
mvco_triplet_dates = string(table2array(MVCOdata062218full(:,1)));
mvco_triplet_times = string(table2array(MVCOdata062218full(:,2))) ; 
mvco_triplet_650 = table2array(MVCOdata062218full(:,3)) ;
mvco_triplet_650_msmnt = table2array(MVCOdata062218full(:,4));
mvco_triplet_695 = table2array(MVCOdata062218full(:,5));
mvco_triplet_695_msmnt = table2array(MVCOdata062218full(:,6));
mvco_triplet_460 = table2array(MVCOdata062218full(:,7));
mvco_triplet_460_msmnt = table2array(MVCOdata062218full(:,8));
mvco_triplet_column_9 = table2array(MVCOdata062218full(:,9));
mvco_date_and_time = mvco_triplet_dates + ' ' + mvco_triplet_times;
mvco_datetimes = datetime(mvco_date_and_time,'InputFormat','MM/dd/yy HH:mm:ss');
mvco_datenums = datenum(mvco_datetimes);
%Convert from EST to UTC time 
mvco_datenums = mvco_datenums + (60*4)/1440 % this adds 1/1440 * 4 hours of time to convert to UTC



% from dev file
scale_factor_650 = 3.274e-06;
offset_650 = 50;

scale_factor_695 = .0121;
offset_695 = 50;

scale_factor_460 = .0904;
offset_460 = 47;

%useable_range = [1:2114 2116:28230];
%useable_range = [1:2114 2116:7325 7327:(7336) (7338):27230];

deployment_time_surface = '04/11/18 15:00:00'; % left at 7 am, left tower by 11 am
deployment_time_datetime_surface = datetime(deployment_time_surface, 'InputFormat','MM/dd/yy HH:mm:ss');
deployment_time_datenum_surface = datenum(deployment_time_datetime_surface);

after_deployment_surface = mvco_datenums > deployment_time_datenum_surface;

after_deployment_found_surface = find(after_deployment_surface);

time_usable = (mvco_datenums(after_deployment_found_surface(1)) < mvco_datenums ) .* (mvco_datenums< mvco_datenums(28700));
usable_650 = mvco_triplet_650 == 650;
usable_695 = mvco_triplet_695 == 695;
usable_460 = mvco_triplet_460 == 460;

usable_650_msmnt = mvco_triplet_650_msmnt < 4160; %4130 +-60 is max output
usable_695_msmnt = mvco_triplet_695_msmnt < 4160;
usable_460_msmnt = mvco_triplet_460_msmnt < 4160;



useable_range = find(time_usable.*usable_650.*usable_460.*usable_695.*usable_650_msmnt.*usable_695_msmnt.*usable_460_msmnt);
%remove the out of water data 
%useable_range = useable_range

%Begin plotting!!!
figure(101);
plot(mvco_datenums(useable_range)); % dat5etickzoom('x', 'HH:MM mm/dd ');
title('datenums'); 
%title('650');   datetickzoom('x', 'HH:MM mm/dd '); ylabel('Wavelength (nm)');% ylim([500 700])

figure(102);

%650
subplot(331);
scatter(mvco_datenums(useable_range), mvco_triplet_650(useable_range)); 
title('650');   datetickzoom('x', 'HH:MM mm/dd ');   ylabel('Wavelength (nm)');% ylim([500 700])

%measurements
subplot(334);
plot(mvco_datenums(useable_range),mvco_triplet_650_msmnt(useable_range)); % datetickzoom('x', 'HH:MM mm/dd ');
title('650 msmnt'); datetickzoom('x', 'HH:MM mm/dd ');  ylabel('Intensity ');% ylim([500 700])

%650
subplot(332);
scatter(mvco_datenums(useable_range), mvco_triplet_695(useable_range)); 
title('695');   datetickzoom('x', 'HH:MM mm/dd ');   ylabel('Wavelength (nm)');% ylim([500 700])

%measurements
subplot(335);
plot(mvco_datenums(useable_range),mvco_triplet_695_msmnt(useable_range)); % datetickzoom('x', 'HH:MM mm/dd ');
title('695 msmnt'); datetickzoom('x', 'HH:MM mm/dd ');  ylabel('Intensity ');% ylim([500 700])

%650
subplot(333);
scatter(mvco_datenums(useable_range), mvco_triplet_460(useable_range)); 
title('460');   datetickzoom('x', 'HH:MM mm/dd ');   ylabel('Wavelength (nm)');% ylim([500 700])

%measurements
subplot(336);
plot(mvco_datenums(useable_range),mvco_triplet_460_msmnt(useable_range)); % datetickzoom('x', 'HH:MM mm/dd ');
title('460 msmnt'); datetickzoom('x', 'HH:MM mm/dd ');  ylabel('Intensity ');% ylim([500 700])



backscatter_triplet = (mvco_triplet_650_msmnt-offset_650)  * scale_factor_650;

chl_triplet = (mvco_triplet_695_msmnt-offset_695)  * scale_factor_695;

cdom_triplet =(mvco_triplet_460_msmnt-offset_460)  * scale_factor_460;

usable_backscatter = backscatter_triplet > 0; %4130 +-60 is max output
usable_chl = chl_triplet > 0;
usable_cdom = cdom_triplet > 0;

useable_range = find(time_usable.*usable_650.*usable_460.*usable_695.*usable_650_msmnt.*usable_695_msmnt.*usable_460_msmnt.*usable_backscatter.*usable_cdom.*usable_chl)

% filter each individual time section
unique_times = unique(mvco_datenums(useable_range));
i=min(useable_range);
while( i < (max(useable_range)+1))
    if( sum(i == useable_range) > 0) %if i is in usable range
        same_times = abs((mvco_datenums(useable_range) - mvco_datenums(i)))< (4/1440); % 3 minutes is 3/1440 in datenum
        median_chl(i) = median(chl_triplet(useable_range(find(same_times))));
        median_backscatter(i) = median(backscatter_triplet(useable_range(find(same_times))));
        median_cdom(i) = median(cdom_triplet(useable_range(find(same_times))));
        i = i+1; %max(find(same_times)) +10
        %    disp('stuck 1')
   
    else
        i = i+1;
   %     disp('stuck 2')
    end
    
end
%load in CDOM and chl data from tower for comparison
load('/home/flagg/Dropbox (MIT)/PikaL_Dev/chl_data_mvco/06_23_18_data/mvco_data_06_23_18/CDOMsmooth2018.mat')
load('/home/flagg/Dropbox (MIT)/PikaL_Dev/chl_data_mvco/06_23_18_data/mvco_data_06_23_18/ecochlsmooth2018.mat')


deployment_time_depth= '03/24/18 10:00:00'; % left at 7 am, left tower by 11 am
deployment_time_datetime_depth = datetime(deployment_time_depth, 'InputFormat','MM/dd/yy HH:mm:ss');
deployment_time_datenum_depth = datenum(deployment_time_datetime_depth);




b_times = ecosmooth(:,1)
b_values = ecosmooth(:,8)

chl_times = ecosmooth(:,1)
chl_values = ecosmooth(:,7)

CDOM_times = CDOMsmooth(:,1)
CDOM_values = CDOMsmooth(:,4)

after_deployment_depth_b = b_times > deployment_time_datenum_depth;
after_deployment_depth_chl = chl_times > deployment_time_datenum_depth;
after_deployment_depth_CDOM = CDOM_times > deployment_time_datenum_depth;

after_deployment_found_depth_b = find(after_deployment_depth_b);
after_deployment_found_depth_chl = find(after_deployment_depth_chl);
after_deployment_found_depth_CDOM = find(after_deployment_depth_CDOM);


surface_datenums = mvco_datenums(useable_range);
surface_median_chl =  median_chl(useable_range);
surface_median_backscatter =  median_backscatter(useable_range);
surface_median_cdom =  median_cdom(useable_range);

depth_datenums_cdom = CDOM_times(after_deployment_found_depth_CDOM);
depth_datenums_chl = chl_times(after_deployment_found_depth_chl);
depth_datenums_b = b_times(after_deployment_found_depth_b);

depth_chl = chl_values(after_deployment_found_depth_chl);
depth_b = b_values(after_deployment_found_depth_b);
depth_cdom = CDOM_values(after_deployment_found_depth_CDOM);


%calibrated real world units 
%measurements
subplot(337);
plot(surface_datenums,surface_median_backscatter); % datetickzoom('x', 'HH:MM mm/dd ');
hold on;
plot(depth_datenums_b,depth_b)
hold off
title('backscatter'); datetickzoom('x', 'HH:MM mm/dd ');  ylabel('Backscatter 650 1/m '); legend('surface', 'depth')% ylim([500 700])

%650
subplot(338);
plot(surface_datenums,surface_median_chl); 
hold on;
plot(depth_datenums_chl,depth_chl)
hold off
title('chl');   datetickzoom('x', 'HH:MM mm/dd ');   ylabel('chl \mug/L'); legend('surface', 'depth')% ylim([500 700])

%measurements
subplot(339);
plot(surface_datenums,surface_median_cdom); % datetickzoom('x', 'HH:MM mm/dd ');
hold on;
plot(depth_datenums_cdom,depth_cdom)
hold off
title('cdom'); datetickzoom('x', 'HH:MM mm/dd ');  ylabel('CDOM ppb ');  legend('surface', 'depth') % ylim([500 700])
 




figure(103)
%calibrated real world units 
%measurements
subplot(311);
plot(surface_datenums,surface_median_backscatter); % datetickzoom('x', 'HH:MM mm/dd ');
hold on;
plot(depth_datenums_b,depth_b)
hold off
title('backscatter'); datetickzoom('x', 'HH:MM mm/dd ');  ylabel('Backscatter 650 1/m '); legend('surface', 'depth')% ylim([500 700])

%650
subplot(312);
plot(surface_datenums, surface_median_chl); 
hold on;
plot(depth_datenums_chl,depth_chl)
hold off
title('chl');   datetickzoom('x', 'HH:MM mm/dd ');   ylabel('chl \mug/L'); legend('surface', 'depth')% ylim([500 700])

%measurements
subplot(313);
plot(surface_datenums,surface_median_cdom); % datetickzoom('x', 'HH:MM mm/dd ');
hold on;
plot(depth_datenums_cdom,depth_cdom)
hold off
title('cdom'); datetickzoom('x', 'HH:MM mm/dd ');  ylabel('CDOM ppb (QSDE)');  legend('surface', 'depth') % ylim([500 700])

save('surfaceAndDepthData.mat','surface_datenums','surface_median_backscatter','surface_median_chl','surface_median_cdom','depth_datenums_b','depth_datenums_chl','depth_datenums_cdom','depth_b','depth_chl','depth_cdom')
