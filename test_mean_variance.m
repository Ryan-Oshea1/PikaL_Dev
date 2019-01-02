

clear mean_rows variance_rows
for(j=1:12)
    exposure_time=(j)*20000
%[mean_rows(:,j),variance_rows(:,j)] = mean_variance_calculator(exposure_time,0,'/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/flats_no_polarizer','/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/flats_polarizer','/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/darks/',0,0,CalFolder,camera_wavelengths)
[mean_rows(:,j),variance_rows(:,j)] = mean_variance_calculator(exposure_time,0,'/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS','/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS_pol','/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS_darks/',0,0,CalFolder,camera_wavelengths)

end

for(j=14:23)
    exposure_time=(j)*20000
%[mean_rows(:,j),variance_rows(:,j)] = mean_variance_calculator(exposure_time,0,'/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/flats_no_polarizer','/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/flats_polarizer','/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/darks/',0,0,CalFolder,camera_wavelengths)
[mean_rows(:,j-1),variance_rows(:,j-1)] = mean_variance_calculator(exposure_time,0,'/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS','/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS_pol','/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS_darks/',0,0,CalFolder,camera_wavelengths)

end
%%
close all
row_location = 150
column_location = 1:22
variance_data = variance_rows(row_location,column_location);
mean_data = mean_rows(row_location,column_location);
scatter(mean_data,variance_data)
xlabel('Variance (counts squared)')
ylabel('Mean (counts)')
title('Retreiving ADC gain from mean and variance')
coefficients = polyfit(mean_data,variance_data,1)

gain = coefficients(1)
hold on
plot(polyval(coefficients,1:700))
hold off


representable_electrons = gain*4095

%all rows
variance_data_all = variance_rows(:)./2;
mean_data_all = mean_rows(:);
figure
scatter(mean_data_all,variance_data_all,[],[.6 .6 .6])
xlabel('Mean (counts)')
ylabel('Variance (counts^2)')
title({['Calculating ADC Gain from'] ['Variance vs. Mean']})

coefficients = polyfit(mean_data_all,variance_data_all,1)
hold on
plot(polyval(coefficients,1:1900),'k','LineWidth',3)
hold off

gain = coefficients(1)
offset = coefficients(2)
representable_electrons = gain*4095
txt = "Y=" + num2str(gain) + "*X + " +  num2str(offset) 
txt_gain = "1/(ADC Gain)=" + num2str(gain)
ADC_gain =1/gain
txt_adc_gain = "ADC Gain=" + num2str(ADC_gain)
text(25,230,txt,'FontSize',25)
text(25,210,txt_gain,'FontSize',25)
text(25,190,txt_adc_gain,'FontSize',25)
legend('Camera Measured','First Order Regression')
set(gca,'FontSize',25)

% plot the residuals
figure(3);
scatter(mean_data_all,variance_data_all - polyval(coefficients,mean_data_all))
axis([0 2000 -40 40 ])
grid on
title('Residuals of Variance Vs. Mean Fit')
%% 
% the next graph I need to make should showcase the radiance, calculate the
% electronics expected noise, and the std of the radiance over the 900
% spatial pixels, compare it to the on board radiometer.
%close all
[snr_image,] = snr_calculator(400000,0,'/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS','/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS_pol','/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS_darks/',0,0,CalFolder,camera_wavelengths,7.9346,gain);


%%
%calculate the SNR using the two different methods on each image from the
%real-wrold dataset
%This will show how well each method can remove the environmental
%variations on a pixel by pixel basis


[snr_image,] = snr_calculator(400000,0,'/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS','/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS_pol','/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS_darks/',0,0,CalFolder,camera_wavelengths,7.9346,ADC_gain);


%then we can see how much that propagates to the different estimation
%algorithms

%%
figure(11)
i= 60
subplot(211);imagesc(all_images(i).set_num_9_cal(:,:)); colorbar; caxis([0 10])
subplot(212);imagesc(all_images(i).set_num_7_cal(:,:)); colorbar; caxis([0 10])

figure(12)
subplot(211);imagesc(all_images(i).set_num_9_rrs(:,:)); colorbar;caxis([0 .07])
subplot(212);imagesc(all_images(i).set_num_7_rrs(:,:)); colorbar; caxis([0 .07])

figure(13)
i= 60
subplot(211);imagesc(all_images(i).set_num_9_cal_pol(:,:)); colorbar; caxis([0 10])
subplot(212);imagesc(all_images(i).set_num_7_cal_pol(:,:)); colorbar; caxis([0 10])

figure(14)
subplot(211);imagesc(all_images(i).set_num_9_rrs_pol(:,:)); colorbar;caxis([0 .07])
subplot(212);imagesc(all_images(i).set_num_7_rrs_pol(:,:)); colorbar; caxis([0 .07])
%%
close all
figure(10)
i=2016%2016%1905%1849%1512%1457%1344%58%337%55
spatial_pixel_number = 100
plot(camera_wavelengths(3,1:150),median(all_images(i).set_num_9_rrs(:,spatial_pixel_number),2))
hold on
plot(camera_wavelengths(1,:),median(all_images(i).set_num_7_rrs_pol(:,spatial_pixel_number),2))
plot(camera_wavelengths(3,1:150),median(all_images(i).set_num_9_rrs_pol(:,spatial_pixel_number),2))
%plot(camera_wavelengths(1,:),median(all_images(i).set_num_7_rrs_pol_2(:,spatial_pixel_number),2))
%plot(camera_wavelengths(1,:),median(all_images(i).set_num_7_rrs_pol_3(:,spatial_pixel_number),2))


legend('uncorrected','rrs_lee','pol_1','pol_2','pol_3')


%%
% now we calculate the theoretical signal-to-noise ratio
% of the real-world signals
%close all
close all
set_num = 9;
camera_wavelength_row = 3 ; 
camera_wavelength_columns = 1:150;
     xrange = 79:228;
     yrange = 132:356;
     spatial_range = 200:225;

%start_point = 1289
% 506/561 are dark high chl, light high chl, 4027,4083 are dark/light of
% lw, (i=[506,561,4027,4083])
% chl
count = 0
for(j= 1:length(find(find_usable_sza(water_quality_parameter_times_chl,115,155))))
  count = count+1;
    j
    k=find(find_usable_sza(water_quality_parameter_times_chl,115,155));
    i = k(j)
  %  if(size(all_images(i).set_num_9_cal_pol,1) == 300)

  if((size(all_images(i).set_num_9_cal_pol,2) >1) && (size(all_images(i).set_num_9_cal,2) >1) && (size(all_images(i).set_num_9_rrs_lee,2) >1) )
%172%171%170169%112
%uncorr '/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS_darks/'
[SNR_image_theoretical_uncorrected, theoretical_snr_uncorrected(:,count),spatial_SNR_uncorrected(:,count), spatial_rel_practical_uncorrected(:,count),environmental_SNR(:,count)] = snr_calculator_real_world_2(all_images(i).set_num_9_cal,all_images(i).set_num_9,exposure_Im_All(i,set_num),0,'/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/darks/',camera_wavelengths,ADC_gain,spatial_range);
%
%pol
[SNR_image_theoretical_polarized, theoretical_snr_polarized(:,count),spatial_SNR_polarized(:,count), spatial_rel_practical_polarized(:,count),environmental_SNR_polarized(:,count)] = snr_calculator_real_world_2(all_images(i).set_num_9_cal_pol,all_images(i).set_num_9_pol,exposure_Im_All_pol(i,set_num),0,'/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/darks/',camera_wavelengths,ADC_gain,spatial_range);
%
%srg NEED to backsolve for the radiance
[SNR_image_theoretical_lee, theoretical_snr_lee(:,count),spatial_SNR_lee(:,count), spatial_rel_practical_lee(:,count),environmental_SNR_lee(:,count)] = snr_calculator_real_world_2(all_images(i).set_num_9_rrs_lee_2_t,all_images(i).set_num_9,exposure_Im_All(i,set_num),0,'/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/darks/',camera_wavelengths,ADC_gain,spatial_range);

% 
%          figure(5);
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),median(all_images(i).set_num_9(xrange,yrange),2)); 
%       hold on;
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),median(all_images(i).set_num_9_pol(xrange,yrange),2)); 
      
      theoretical_snr_holder(count,:) = nanmedian(theoretical_snr_uncorrected(:,count),2);
      spatial_snr_holder(count,:) = nanmedian(spatial_SNR_uncorrected(:,count),2);
      radiance_holder(count,:) = nanmedian(all_images(i).set_num_9_cal,2);
      %plot typical radiance above typical SNR
%        figure(9);
%        subplot(211)
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),median(all_images(i).set_num_9_cal,2)); 
%       hold on;
%  %     plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),median(all_images(i).set_num_9_cal_pol,2)); 
%       %legend('Uncorrected')%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')
%       title('Spectral Radiance')
%       xlabel('Wavelength (nm)')
%       ylabel('Radiance (W/m^2/um/sr)')
%       set(gca,'FontSize',20)
%       subplot(212)
%       hold on
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),nanmedian(theoretical_snr_uncorrected(:,count),2)); 
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),median(spatial_SNR_uncorrected(:,count),2)); 
% 
%       legend('Shot-noise limited SNR', 'Spatially Derived SNR')%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')
%       title('Theoretical Spectral SNR')
%       xlabel('Wavelength (nm)')
%       ylabel('SNR')
%       set(gca,'FontSize',20)
% 
% 
%      % hold on ; plot(camera_wavelengths(1,:),median_SNR); hold off; 
%      % plot(camera_wavelengths(1,:),mean_rad./STD_rad); hold off;
%       legend('uncorrected/Lee','Polarized')%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')
%       title('Raw Digital Numbers')
%       xlabel('Wavelength (nm)')
%       ylabel('Digital Numbers')
%       set(gca,'FontSize',20)
%       hold off
% 
% 
%          figure(6);
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_rel_practical_uncorrected); 
%       hold on;
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_rel_practical_polarized); 
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_rel_practical_lee); 
% 
% 
%      % hold on ; plot(camera_wavelengths(1,:),median_SNR); hold off; 
%      % plot(camera_wavelengths(1,:),mean_rad./STD_rad); hold off;
%       legend('uncorrected','Polarized','Lee')%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')
%       title('Spatial Rel. to Theoretical')
%       xlabel('Wavelength (nm)')
%       ylabel('rel. SNR')
%       set(gca,'FontSize',20)
%       hold off
%       axis([400 1000 0 1 ])
% 
%       figure(7);
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_SNR_uncorrected); 
%       hold on;
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_SNR_polarized); 
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_SNR_lee); 
%        
% 
%      % hold on ; plot(camera_wavelengths(1,:),median_SNR); hold off; 
%      % plot(camera_wavelengths(1,:),mean_rad./STD_rad); hold off;
%       legend('uncorrected','Polarized','Lee')%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')
%       title('Spatial SNR')
%       xlabel('Wavelength (nm)')
%       ylabel('rel. SNR')
%       set(gca,'FontSize',20)
%       hold off
%       axis([400 1000 0 60 ])
%       
%             figure(8);
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),environmental_SNR); 
%       hold on;
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),environmental_SNR_polarized); 
%       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),environmental_SNR_lee); 
%        
% 
%      % hold on ; plot(camera_wavelengths(1,:),median_SNR); hold off; 
%      % plot(camera_wavelengths(1,:),mean_rad./STD_rad); hold off;
%       legend('uncorrected','Polarized','Lee')%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')
%       title('Environmental SNR')
%       xlabel('Wavelength (nm)')
%       ylabel('rel. SNR')
%       set(gca,'FontSize',20)
%       hold off
%       axis([400 1000 0 250 ])
%     
% hold the data in an array so that I can take the spectral median of
% unpolarized/polarized/srg
%      waitforbuttonpress
%pause(5)
end
end
%%
 theoretical_snr_holder(count,:) = nanmedian(theoretical_snr_uncorrected,2);
      spatial_snr_holder(count,:) = nanmedian(spatial_SNR_uncorrected,2);
      radiance_holder(count,:) = nanmedian(all_images(i).set_num_9_cal,2);
      %%
      %plot typical radiance above typical SNR
       figure(9);
       subplot(211)
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),radiance_holder(1,:),'k','LineWidth',2); 
      hold on;
            plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),radiance_holder(2,:),'m','LineWidth',2); 
                  plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),radiance_holder(3,:),'b','LineWidth',2); 
                        plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),radiance_holder(4,:),'r','LineWidth',2); 



      
      
 %     plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),median(all_images(i).set_num_9_cal_pol,2)); 
      %legend('Uncorrected')%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')
      title('Spectral Radiance')
      xlabel('Wavelength (nm)')
      ylabel({'Radiance'  '(W/m^2/um/sr)'})
            legend('Dark High Chl.', 'Bright High Chl.','Dark Low Chl.', 'Bright Low Chl.')%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')

      set(gca,'FontSize',20)
      axis([380 1000 0 12])
      subplot(212)
     %       axis([380 1000 0 80])

      hold on
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),theoretical_snr_holder(1,:)/sqrt(2),'k','LineWidth',2); 
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_snr_holder(1,:),'k-.','LineWidth',2); 

      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),theoretical_snr_holder(2,:)/sqrt(2),'m','LineWidth',2); 
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_snr_holder(2,:),'m-.','LineWidth',2); 
      
      
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),theoretical_snr_holder(3,:)/sqrt(2),'b','LineWidth',2); 
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_snr_holder(3,:),'b-.','LineWidth',2); 
      
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),theoretical_snr_holder(4,:)/sqrt(2),'r','LineWidth',2); 
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_snr_holder(4,:),'r-.','LineWidth',2); 
      legend('Shot-noise limited SNR', 'Spatially derived SNR')%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')
      title('Spectral SNR')
      xlabel('Wavelength (nm)')
      ylabel('SNR')
      set(gca,'FontSize',20)
      
      %%
      %plot median SNR's
      %       figure(7);
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),nanmedian(spatial_SNR_uncorrected,2),'k-.','LineWidth',2);
      hold on;
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),nanmedian(spatial_SNR_polarized,2),'b-.','LineWidth',2);
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),nanmedian(spatial_SNR_lee,2),'m-.','LineWidth',2);
      
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),nanmedian(theoretical_snr_uncorrected,2),'k','LineWidth',2); 
      hold on;
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),nanmedian(theoretical_snr_polarized,2),'b','LineWidth',2);
      %plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),nanmedian(theoretical_snr_lee,2)); 
      


     % hold on ; plot(camera_wavelengths(1,:),median_SNR); hold off; 
     % plot(camera_wavelengths(1,:),mean_rad./STD_rad); hold off;
      legend('Uncorr. spatial','Pol. spatial','Lee spatial','Uncorr. shot-noise','Pol. shot-noise')%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')
      title('Median SNR of Timeseries')
      xlabel('Wavelength (nm)')
      ylabel('SNR')
      set(gca,'FontSize',20)
      hold off
      axis([400 1000 0 100 ])