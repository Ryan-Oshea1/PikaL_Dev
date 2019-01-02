%Ryan O'Shea
%04/12/18
%This function reads in the nearest calibration file that does not have
%saturated pixels and the nearest dark value. It provides the dark
%corrected calibration image and the 
% i want to pass in the exposure time that I am interested in, the location
% of the darks, the location of the 100 cal images, the folder for the
% radiances, and have this function spit out the mean and variance for each
% line. I can then use that data to make a graph

function [SNR_image_theoretical, theoretical_snr,spatial_SNR, spatial_rel_practical,environmental_SNR] = snr_calculator_real_world(original_image_radiance,original_image_raw,exposure_time,gain,path_to_darks,camera_wavelengths,ADC_gain,spatial_range)

%find the image with the nearest exposure time
image_exposure_time = exposure_time;
image_exposure_time_str = num2str(image_exposure_time);

%find the exposure time
if(exposure_time >480000)
     rounded_exposure_time = round(480000/20000) *20000;
else
 rounded_exposure_time = round(exposure_time/20000) *20000;
end
 final_exposure_time = rounded_exposure_time +21;
 final_exposure_time_str = num2str(final_exposure_time);

 %define variables for a single image
     path_set = "1";
     bin = 1;
     xrange = 316:915;
     yrange = 528:1427;
     size(original_image_raw)
    %resize image
    original_image_raw_resized = double(original_image_raw(xrange,yrange));
     
    %read the image
  %  figure(1);imagesc(original_image_raw_resized);

% if the currently selected exposure is not divisibile, then we round to
% find nearest dark
  rounded_exposure_time_image = round(image_exposure_time/20000) *20000;
   if(rounded_exposure_time_image ~= 500000)
 final_exposure_time_image = rounded_exposure_time_image +21;
   else
 final_exposure_time_image = 480021;
   end
 
 final_exposure_time_str_image = num2str(final_exposure_time_image);
     
 %find the paths to the dark images, and the dark image and read noise
 %values
      full_path_darks_image = path_to_darks + "set_" + path_set + "_" + final_exposure_time_str_image + "_" + num2str(gain) + "/AEG/autoExpAndGain" + num2str(1) + ".png"; % '/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/darks/set_'
      dark_img_sub = png_value_convertor(char(full_path_darks_image ));
      dark_img = double(dark_img_sub(xrange,yrange));
      
      full_path_read_noise = path_to_darks + "set_" + path_set + "_" + num2str(21) + "_" + num2str(gain) + "/AEG/autoExpAndGain" + num2str(1) + ".png";
      read_noise_image = png_value_convertor(char(full_path_read_noise ));
      read_noise = double(read_noise_image(xrange,yrange));
      %%%%%%%%%CURRENTLY I AM REMOVING THE IMAGINARY AND NEGATIVE DATA
      %remove dark noise from radiance im
      original_image_raw_unbiased = double(original_image_raw_resized - dark_img);
      original_image_raw_unbiased_electrons = ADC_gain*original_image_raw_unbiased;
      original_image_raw_unbiased_electrons(abs(imag(original_image_raw_unbiased_electrons))>0) = NaN;
      original_image_raw_unbiased_electrons(original_image_raw_unbiased_electrons<0) = NaN;
      
      
      full_signal_electrons = double(original_image_raw_resized) *ADC_gain;

      SNR_image_theoretical = full_signal_electrons./sqrt(original_image_raw_unbiased_electrons + (dark_img - read_noise).^2 + (read_noise).^2);
        SNR_image_theoretical(SNR_image_theoretical<0) = NaN;
        SNR_image_theoretical(abs(imag(SNR_image_theoretical))>0) = NaN;
        
      std_raw_data = nanstd(original_image_raw_unbiased_electrons,0,2);
      median_SNR = nanmedian(SNR_image_theoretical,2);

      size(median_SNR)
     
      %spatial_range = 860:900
      STD_rad = std(original_image_radiance(:,spatial_range),0,2);
      mean_rad = nanmean(original_image_radiance(:,spatial_range),2);
           
      
      
      %plot the raw read/dark/og raw image for clarity
%       figure(1); subplot(311); imagesc(read_noise); colorbar; title('read noise electrons'); caxis([0 80])
%       subplot(312); imagesc(dark_img); colorbar; title('dark noise electrons');  caxis([0 80]);
%       subplot(313); imagesc(original_image_raw_unbiased_electrons); colorbar; caxis([-1 0]); title('signal  electrons'); 
%       
%       %plpot SNR image, and the og raw image
%       figure(2); subplot(211); imagesc(imag(SNR_image_theoretical)); colorbar; ;title('pixel to pixel SNR calculated');
%       subplot(212); imagesc(original_image_raw_unbiased_electrons.^(.5)); colorbar; caxis([0 80]);title('pixel to pixel SNR photon shot-noise estimated');

      %plot the difference between shot nosie limited?
       %     figure(3); imagesc(SNR_image_theoretical - original_image_raw_unbiased_electrons.^(.5)); colorbar; caxis([0 1]);title('SNR Difference');
      
      theoretical_snr = median(SNR_image_theoretical,2);
      practical_theoretical_snr_limit = theoretical_snr/sqrt(2);
      spatial_SNR = mean_rad./STD_rad;
      radiance_noise_theoretical = mean_rad./theoretical_snr;
      1./theoretical_snr;
      1./spatial_SNR;
      environmental_SNR = 1./(1./spatial_SNR - 1./theoretical_snr);
      
      spatial_rel_practical = spatial_SNR./theoretical_snr; 
      
%       figure(4);
%       plot(camera_wavelengths(1,:),theoretical_snr); 
%       hold on
%       plot(camera_wavelengths(1,:),practical_theoretical_snr_limit); 
% 
%      % hold on ; plot(camera_wavelengths(1,:),median_SNR); hold off; 
%       plot(camera_wavelengths(1,:),spatial_SNR); hold off;
%       legend('Theor. SNR','Cal. Lim. Theor. SNR','Spatial SNR')%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')
%       title('Theoretical SNR vs. Spatial SNR')
%       xlabel('Wavelength (nm)')
%       ylabel('SNR')
%       set(gca,'FontSize',20)
%       
%       
%             figure(5);
%       plot(camera_wavelengths(1,:),spatial_rel_practical); 
% 
% 
%      % hold on ; plot(camera_wavelengths(1,:),median_SNR); hold off; 
%      % plot(camera_wavelengths(1,:),mean_rad./STD_rad); hold off;
%       legend('spatial rel. practical')%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')
%       title('Spatial Rel. to Practical')
%       xlabel('Wavelength (nm)')
%       ylabel('rel. SNR')
%       set(gca,'FontSize',20)
      
%       figure(5)
%       imagesc(full_image_cal_3_median_100)
%       colorbar
%       
%       
%       figure(5); imagesc(dark_img)t
%       
%  %Subtract darks from each image
%  hyperspectralImageAE_cal_unbiased = double(hyperspectralImageAE_cal - dark_img);
%  hyperspectralImageAE_cal_2_unbiased = double(hyperspectralImageAE_cal_2 - dark_img);
%   figure(1); subplot(211); imagesc(hyperspectralImageAE_cal_unbiased);colorbar; subplot(212); imagesc(hyperspectralImageAE_cal_2_unbiased);colorbar;  
%   
%  size(median(hyperspectralImageAE_cal_unbiased,2))
%  med_hyp = median(hyperspectralImageAE_cal_unbiased,2)
%  med_hyp_2 = median(hyperspectralImageAE_cal_2_unbiased,2)
% 
%  ratio_a_b = double(med_hyp) ./ double(med_hyp_2);
%  figure(2); plot(camera_wavelengths(1,:),med_hyp);  hold on; plot(camera_wavelengths(1,:),med_hyp_2);
%  figure(3); plot( ratio_a_b)  
%  
%  pixel_to_pixel_variations =  hyperspectralImageAE_cal_unbiased - ratio_a_b .* hyperspectralImageAE_cal_2_unbiased;
%    figure(4); imagesc(pixel_to_pixel_variations);colorbar; 
% 
%  variance_rows =                std(pixel_to_pixel_variations,0,2).^2;
%  mean_rows = mean(hyperspectralImageAE_cal_unbiased,2);
               
               
%interpolate radiance values to camera wavelengths

%interpolated_radiance_camera_wavelengths(:) = interp1(wavelength835D,spectrum835D(1:255),camera_wavelengths(bin,:));
%camera_image_of_radiometer = ones(size(image)) .* (squeeze(interpolated_radiance_camera_wavelengths(1,1:size(image))))';

%need to subtract darks from hyperspectralImageAE_cal and from image, but
%need to select the closest exposure time darks
%calibration_image = camera_image_of_radiometer ./ double(hyperspectralImageAE_cal-dark_cal) * double(str2num(final_exposure_time_str));    


%calibrated_image = double(image-dark_img)/double(image_exposure_time).*double(calibration_image);
        
        
% output: camera_image_of_radiometer ./ double(all_images(set_num).spectrum_image_cal) * double(str2num(all_images(set_num).cal_exposure));
%calibrated_image=calibrated_image;
end