%removes the glinty and dark images from the dataset, by looking at the SZA

for(i=1:4612)%4612)% goes through all of the images
      if(1 == find_usable_sza(water_quality_parameter_times_chl(i),105,165))
          % calculate the software glint corrected image
          disp('software')
          all_images(i).set_num_7_rrs_lee =  software_correct_all_data(  all_images, camera_wavelengths, water_quality_parameter_times_chl, water_quality_parameter_values_chl, wavelength_range_numer_chl,wavelength_range_denom_chl,order_chl,set_num,set_num_name_rrs_array,noise_cutoff,noise_min,tss_range_end,verbose,spatial_angle_start,spatial_angle_end,angle_spacing,indices_training_days_chl,indices_validation_days_chl,1,wavelengths_rrs_opt,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,wavelength835B,srs_spectrum,spatial_pixel_angles,i);

      else
           % remove the glinty or dark image
          disp('remove')
          all_images(i).set_num_7_rrs = []
          all_images(i).set_num_7_rrs_pol = []
      end
      
end
%%
%removes the glinty and dark images from the dataset, by looking at the SZA
set_num = 9
for(i=3470:4612)%4612)% goes through all of the images
      if(1 == find_usable_sza(water_quality_parameter_times_chl(i),105,165))
          % calculate the software glint corrected image
          disp('software')
          %all_images(i).set_num_7_pol = all_images_pol_105_165(i).set_num_7_pol;
          %all_images_pol_105_165(i).set_num_7_pol = [];
          
          %all_images(i).set_num_7 = all_images_unpol_105_165(i).set_num_7;
          %all_images_unpol_105_165(i).set_num_7 = [];
          all_images(i).set_num_9_rrs_lee =  software_correct_all_data2(  all_images, camera_wavelengths, water_quality_parameter_times_chl, water_quality_parameter_values_chl, wavelength_range_numer_chl,wavelength_range_denom_chl,order_chl,set_num,set_num_name_rrs_array,noise_cutoff,noise_min,tss_range_end,verbose,spatial_angle_start,spatial_angle_end,angle_spacing,indices_training_days_chl,indices_validation_days_chl,1,wavelengths_rrs_opt,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,wavelength835B,srs_spectrum,spatial_pixel_angles,i);

      else
           % remove the glinty or dark image
          disp('remove')
         % all_images(i).set_num_7_pol = []
          %all_images(i).set_num_7_rrs_pol = []
      end
      
end
%%
%removes the glinty and dark images from the dataset, by looking at the SZA
%TESTER
set_num = 9
for(i=1)%4612)% goes through all of the images
      if(1 == find_usable_sza(water_quality_parameter_times_chl(i),105,165))
          % calculate the software glint corrected image
          disp('software')
          %all_images(i).set_num_7_pol = all_images_pol_105_165(i).set_num_7_pol;
          %all_images_pol_105_165(i).set_num_7_pol = [];
          
          %all_images(i).set_num_7 = all_images_unpol_105_165(i).set_num_7;
          %all_images_unpol_105_165(i).set_num_7 = [];
          [rrs_lee_image,trs_image,srs_image,fresnel_array,initial_values_array] =  software_correct_all_data_test(  all_images, camera_wavelengths, water_quality_parameter_times_chl, water_quality_parameter_values_chl, wavelength_range_numer_chl,wavelength_range_denom_chl,order_chl,set_num,set_num_name_rrs_array,noise_cutoff,noise_min,tss_range_end,verbose,spatial_angle_start,spatial_angle_end,angle_spacing,indices_training_days_chl,indices_validation_days_chl,1,wavelengths_rrs_opt,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,wavelength835B,srs_spectrum,spatial_pixel_angles,i);

      else
           % remove the glinty or dark image
          disp('remove')
         % all_images(i).set_num_7_pol = []
          %all_images(i).set_num_7_rrs_pol = []
      end
      
end
%%
xlswrite('/home/flagg/Dropbox (MIT)/PikaL_Dev/excel/test_rrs_data',rrs_lee_image)
xlswrite('/home/flagg/Dropbox (MIT)/PikaL_Dev/excel/test_trs_data',trs_image)
xlswrite('/home/flagg/Dropbox (MIT)/PikaL_Dev/excel/test_srs_data',srs_image)
xlswrite('/home/flagg/Dropbox (MIT)/PikaL_Dev/excel/test_fresnel_data',fresnel_array)
xlswrite('/home/flagg/Dropbox (MIT)/PikaL_Dev/excel/initial_values_array',initial_values_array)


%%
set_num = 9
          software_correct_all_data2(  all_images, camera_wavelengths, water_quality_parameter_times_chl, water_quality_parameter_values_chl, wavelength_range_numer_chl,wavelength_range_denom_chl,order_chl,set_num,set_num_name_rrs_array,noise_cutoff,noise_min,tss_range_end,verbose,spatial_angle_start,spatial_angle_end,angle_spacing,indices_training_days_chl,indices_validation_days_chl,1,wavelengths_rrs_opt,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,wavelength835B,srs_spectrum,spatial_pixel_angles,1);

          %%
          close all
          image_number = 1;
          spatial_pixel_loc = 21%109%76%414;
          input_image_rrs = all_images(image_number).set_num_7_rrs;
          input_image_rrs_lee = all_images(image_number).set_num_7_rrs_lee;
          input_image_rrs_pol = all_images(image_number).set_num_7_rrs_pol;

          figure(1)
          subplot(211)
          imagesc(all_images(image_number).set_num_7_rrs);colorbar;caxis([-.001 .0055])
          subplot(212)
          imagesc(all_images(image_number).set_num_7_rrs_lee);colorbar;caxis([0 .005])
              
          figure(2)
          subplot(211)
          plot(nanmean(input_image_rrs(:,spatial_pixel_loc),2))
          subplot(212)
          hold on
         % plot(nanmean(input_image_rrs_lee(:,spatial_pixel_loc),2)+.9906)
          plot(camera_wavelengths(1,:),nanmean(input_image_rrs_lee(:,spatial_pixel_loc-1),2))
          plot(camera_wavelengths(1,:),nanmean(input_image_rrs_pol(:,spatial_pixel_loc),2)*2.5)
          plot(camera_wavelengths(1,:),nanmean(input_image_rrs(:,spatial_pixel_loc),2))

          hold off
          legend('lee','pol','uncorr')
%%
           figure(3)
           hold on
          subplot(211)
          plot(nanmean(input_image_rrs(:,spatial_pixel_loc),2))
          subplot(212)
          plot(input_image_rrs_lee(:,spatial_pixel_loc)-input_image_rrs_lee(:,spatial_pixel_loc-1))
          nanmean(input_image_rrs_lee(:,spatial_pixel_loc)-input_image_rrs_lee(:,spatial_pixel_loc-1))