%Ryan O'Shea
%04/12/18
%This function reads in the nearest calibration file that does not have
%saturated pixels and the nearest dark value. It provides the dark
%corrected calibration image and the 
% i want to pass in the exposure time that I am interested in, the location
% of the darks, the location of the 100 cal images, the folder for the
% radiances, and have this function spit out the mean and variance for each
% line. I can then use that data to make a graph

function [SNR_image] = snr_calculator(exposure_time,gain,path_to_unpolarized_cal_images,path_to_polarized_cal_images,path_to_darks,select_polarized_images,select_gain_images,CalFolder,camera_wavelengths,ADC_gain,ADC_gain2)
%find the image with the nearest exposure time
figureNumber=1
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

     path_set = "1";
     bin = 1;
     xrange = 316:915;
     yrange = 528:1427;
  

  
%if the calibration image is saturated, decrement the exposure time image
%that is read, and the exposure time, do this iteratively until the
%calibration image is not saturated at all   
  %make sure to calibrate on the closest non-saturated image
  max_image_value = 4095;
  saturation_cutoff_value = 4000;
  while ( (max_image_value > saturation_cutoff_value) )
      
    %create the path for the image and the radiometer
    if select_gain_images == 0 && select_polarized_images == 0
    full_image_path = path_to_unpolarized_cal_images  + "/set_" + path_set + "_" + final_exposure_time_str + "_" + num2str(gain) + "/AEG/autoExpAndGain1.png"; %pol_set replaced by set
    full_image_path_2 = path_to_unpolarized_cal_images  + "/set_" + path_set + "_" + final_exposure_time_str + "_" + num2str(gain) + "/AEG/autoExpAndGain2.png"; %pol_set replaced by set

    full_image_path_rad = path_to_unpolarized_cal_images  + "/rad_" + final_exposure_time_str + "_" + num2str(gain) ;
    else
            if select_gain_images == 0 && select_polarized_images == 1
                 full_image_path = path_to_unpolarized_cal_images  + "/set_" + path_set + "_" + final_exposure_time_str + "_" + num2str(gain) + "/AEG/autoExpAndGain1.png";
                 full_image_path_rad = path_to_unpolarized_cal_images  + "/rad_" + final_exposure_time_str + "_" + num2str(gain) ;
            end %full_image_path = path_to_polarized_cal_images  + "/pol_set_" + path_set + "_" + final_exposure_time_str + "_" + num2str(gain)
    end

    %read the image
    hyperspectralImageAE_cal=png_value_convertor(char(full_image_path ));

   % imagesc(hyperspectralImageAE_cal);
    max_image_value = max(max(hyperspectralImageAE_cal));
    hyperspectralImageAE_cal = hyperspectralImageAE_cal(xrange,yrange);
    
    hyperspectralImageAE_cal_2 = png_value_convertor(char(full_image_path_2));
    hyperspectralImageAE_cal_2 = hyperspectralImageAE_cal_2(xrange,yrange);
   
    % if there is no radiometer data, read a different file name
        datafile_835D = [char(full_image_path_rad)];
        datafile_835D_ramses =  dir((datafile_835D));
        if (length(datafile_835D_ramses) <3)
            disp('nan wavelength 835D');
        else
%             datafile_835D_ramses = extractfield(datafile_835D_ramses,'name');
%             datafile_835D_ramses = datafile_835D_ramses(3);
%             datafile_835D = [datafile_835D '/' datafile_835D_ramses{1}];
% %             [wavelength835D, spectrum835D(1:255)] = ramses_get_835B([CalFolder '835D_deployedRadiometer/'],'835D', datafile_835D);
%             if isnan(wavelength835D)
%                 disp('nan wavelength 835D');
%             end
        end    
        
  end
  
  

% if the currently selected exposure is not divisibile, then we round to
% find nearest dark
 rounded_exposure_time_dark = round(exposure_time/20000) *20000;
 if(rounded_exposure_time_dark ~= 500000)
 final_exposure_time_dark = rounded_exposure_time_dark +21;
 else
  final_exposure_time_dark = 480021;
 end
 final_exposure_time_str_dark = num2str(final_exposure_time_dark);
 
  rounded_exposure_time_image = round(image_exposure_time/20000) *20000;
   if(rounded_exposure_time_image ~= 500000)
 final_exposure_time_image = rounded_exposure_time_image +21;
   else
 final_exposure_time_image = 480021;
   end
 
 final_exposure_time_str_image = num2str(final_exposure_time_image);
 
%read in the dark values for both the calibration image and the original
%image
full_path_darks_calibration = path_to_darks + "set_" + path_set + "_" + final_exposure_time_str_dark + "_" + num2str(gain) + "/AEG/autoExpAndGain1.png"
full_path_darks_image = path_to_darks + "set_" + path_set + "_" + final_exposure_time_str_image + "_" + num2str(gain) + "/AEG/autoExpAndGain1.png"; % '/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/darks/set_'

                dark_cal = imread(char(full_path_darks_calibration ));
                dark_cal_high = bitshift(dark_cal, 16-12);
                dark_cal_low= bitshift(dark_cal, 16-28);
                dark_cal=dark_cal_high+dark_cal_low;
                dark_cal = dark_cal(xrange,yrange);
              %  figure; imagesc(dark_cal);

                dark_img = imread(char(full_path_darks_image ));
                dark_img_high = bitshift(dark_img, 16-12);
                dark_img_low= bitshift(dark_img, 16-28);
                dark_img=dark_img_high+dark_img_low;
                dark_img = dark_img(xrange,yrange);

               % figure; imagesc(dark_img);
                
               
  for(dark_index=   1:100)
      full_path_darks_image = path_to_darks + "set_" + path_set + "_" + final_exposure_time_str_image + "_" + num2str(gain) + "/AEG/autoExpAndGain" + num2str(dark_index) + ".png"; % '/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/darks/set_'
      dark_img_sub = png_value_convertor(char(full_path_darks_image ));
      dark_img_holder{dark_index} = dark_img_sub(xrange,yrange);
      dark_index
      
      full_path_read_noise = path_to_darks + "set_" + path_set + "_" + num2str(21) + "_" + num2str(gain) + "/AEG/autoExpAndGain" + num2str(dark_index) + ".png";
      read_noise_image = png_value_convertor(char(full_path_read_noise ));
      read_noise_image_holder{dark_index} = read_noise_image(xrange,yrange);
      
      full_image_path_cal = path_to_unpolarized_cal_images  + "/set_" + path_set + "_" + num2str(440021) + "_" + num2str(gain) + "/AEG/autoExpAndGain" + num2str(dark_index) + ".png" %pol_set replaced by set
      full_image_cal = png_value_convertor(char(full_image_path_cal ));
      full_image_cal_holder{dark_index} = full_image_cal(xrange,yrange);

      
      
      
  end

      full_image_cal_3_dim = cat(3,full_image_cal_holder{:});
      full_image_cal_3_median_100 = median(full_image_cal_3_dim,3);
      full_image_cal_3_median_1 = full_image_cal_3_dim(:,:,1)
      full_image_cal_3_median_30 = median(full_image_cal_3_dim(:,:,1:30),3);

  
      read_noise_images_3_dim = cat(3,read_noise_image_holder{:});
      read_noise_3_median = median(read_noise_images_3_dim,3);
      read_noise_median_electrons = ADC_gain*double(read_noise_3_median);
            read_noise_median_electrons2 = ADC_gain2*double(read_noise_3_median);

      dark_images_3_dim = cat(3,dark_img_holder{:});         
      dark_images_3_median = median(dark_images_3_dim,3);         
      dark_img = dark_images_3_median;
      dark_image_median_electrons = ADC_gain*double(dark_img);
            dark_image_median_electrons2 = ADC_gain2*double(dark_img);

      hyperspectralImageAE_cal_unbiased = double(hyperspectralImageAE_cal - dark_img);
      hyperspectralImageAE_cal_unbiased_electrons = ADC_gain*hyperspectralImageAE_cal_unbiased;
            hyperspectralImageAE_cal_unbiased_electrons2 = ADC_gain2*hyperspectralImageAE_cal_unbiased;

      full_signal_electrons = double(hyperspectralImageAE_cal) *ADC_gain;
            full_signal_electrons2 = double(hyperspectralImageAE_cal) *ADC_gain2;

      SNR_image = full_signal_electrons./sqrt(hyperspectralImageAE_cal_unbiased_electrons + (dark_image_median_electrons - read_noise_median_electrons).^2 + (read_noise_median_electrons).^2);
            SNR_image2 = full_signal_electrons2./sqrt(hyperspectralImageAE_cal_unbiased_electrons2 + (dark_image_median_electrons2 - read_noise_median_electrons2).^2 + (read_noise_median_electrons2).^2);

      STD_plot = std(hyperspectralImageAE_cal_unbiased_electrons,0,2);
      size(STD_plot)
      median_SNR = median(SNR_image,2);
            median_SNR2 = median(SNR_image2,2);

      size(median_SNR)
      
      radiometrically_corrected_1 = calibrated_camera_image_dark( png_value_convertor(char(full_image_path )),exposure_time,0,'/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/integrating_sphere_calibrations/no_polarizer_MVCO_unscratched_08_27_18','/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/integrating_sphere_calibrations/polarizer_MVCO_unscratched_08_27_18',0,0,CalFolder,camera_wavelengths,(dark_img),full_image_cal_3_median_1);      
      radiometrically_corrected_30 = calibrated_camera_image_dark( png_value_convertor(char(full_image_path )),exposure_time,0,'/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/integrating_sphere_calibrations/no_polarizer_MVCO_unscratched_08_27_18','/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/integrating_sphere_calibrations/polarizer_MVCO_unscratched_08_27_18',0,0,CalFolder,camera_wavelengths,(dark_img),full_image_cal_3_median_30);      
      radiometrically_corrected_100= calibrated_camera_image_dark( png_value_convertor(char(full_image_path )),exposure_time,0,'/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/integrating_sphere_calibrations/no_polarizer_MVCO_unscratched_08_27_18','/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/integrating_sphere_calibrations/polarizer_MVCO_unscratched_08_27_18',0,0,CalFolder,camera_wavelengths,(dark_img),full_image_cal_3_median_100);      

      
      
      
      STD_rad_1 = std(radiometrically_corrected_1,0,2);
      mean_rad_1 = mean(radiometrically_corrected_1,2);
            
      STD_rad_30 = std(radiometrically_corrected_30,0,2);
      mean_rad_30 = mean(radiometrically_corrected_30,2);    
     
      STD_rad_100 = std(radiometrically_corrected_100,0,2);
      mean_rad_100 = mean(radiometrically_corrected_100,2);
      
      figure(1); subplot(311); imagesc(read_noise_median_electrons); colorbar; title('read noise electrons'); caxis([0 80])
      subplot(312); imagesc(dark_image_median_electrons); colorbar; title('dark noise electrons');  caxis([0 80]);
      subplot(313); imagesc(hyperspectralImageAE_cal_unbiased_electrons); colorbar; title('signal  electrons'); 
      
      
      figure(2); subplot(211); imagesc(SNR_image); colorbar; caxis([0 80]);title('pixel to pixel SNR calculated');
      subplot(212); imagesc(hyperspectralImageAE_cal_unbiased_electrons.^(.5)); colorbar; caxis([0 80]);title('pixel to pixel SNR photon shot-noise estimated');

            figure(3); imagesc(SNR_image - hyperspectralImageAE_cal_unbiased_electrons.^(.5)); colorbar; caxis([0 1]);title('SNR Difference');
      
      
      figure(4);
                 plot(camera_wavelengths(1,:),median_SNR2,'LineWidth',2,'color',[0 0 0]); hold off; 
            hold on ; plot(camera_wavelengths(1,:),median_SNR2./sqrt(2),'LineWidth',2,'color',[.5 .5 .5  ]); hold off; 

       hold on ;    plot(camera_wavelengths(1,:),median(full_signal_electrons,2)./STD_plot,'LineWidth',2,'color',[.75 .75 .75 ]); 
     % hold on ; plot(camera_wavelengths(1,:),median_SNR); hold off; 
      hold on; plot(camera_wavelengths(1,:),mean_rad_100./STD_rad_100,'--','LineWidth',2,'color',[0 0 0]); hold off
      hold on; plot(camera_wavelengths(1,:),mean_rad_1./STD_rad_1,'-.','LineWidth',2,'color',[.5 .5 .5  ]); hold off;
  %          hold on; plot(camera_wavelengths(1,:),mean_rad_30./STD_rad_30); hold off;

      legend_1 = legend('Theor. SNR',['Theor. SNR/$$\sqrt{2}$$'], 'Spatial SNR','Rad. Spatial SNR 100','Rad. Spatial SNR 1')%'Rad. Spatial SNR 30'
      set(legend_1,'Interpreter', 'latex')
      title('Theoretical SNR vs. Spatial SNR')
      xlabel('Wavelength (nm)')
      ylabel('SNR')
      set(gca,'FontSize',20)
      xlim([380 1050])
      figure(5)
      imagesc(full_image_cal_3_median_100)
      colorbar
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