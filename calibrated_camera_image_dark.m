%Ryan O'Shea
%04/12/18
%This function reads in the nearest calibration file that does not have
%saturated pixels and the nearest dark value. It provides the dark
%corrected calibration image and the 

function calibrated_image = calibrated_camera_image_dark(image,exposure_time,gain,path_to_unpolarized_cal_images,path_to_polarized_cal_images,select_polarized_images,select_gain_images,CalFolder,camera_wavelengths,dark_image_median,calibration_image_median)
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
 %break out
if ((size(image,1) ~= 300) &&  (size(image,1) ~=  600) &  (size(image,1) ~= 1200))
    disp('incorrect size, will not calibrate')
    calibrated_image = nan
    return
end
% determine which level of binning
 if size(image,1) == 300
     path_set = "3";
     bin =3;
     xrange = 79:228;
     yrange = 132:356;
 end
 
 if size(image,1) == 600
     path_set = "2";
     bin = 2;
     xrange = 158:457;
     yrange = 264:713;
 end
 
  if size(image,1) == 1200
     path_set = "1";
     bin = 1;
     xrange = 316:915;
     yrange = 528:1427;
  end
  
  %resize image
  image = image(xrange,yrange);

  
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
    full_image_path_rad = path_to_unpolarized_cal_images  + "/rad_" + final_exposure_time_str + "_" + num2str(gain) ;
    else
            if select_gain_images == 0 && select_polarized_images == 1
                 full_image_path = path_to_polarized_cal_images  + "/set_" + path_set + "_" + final_exposure_time_str + "_" + num2str(gain) + "/AEG/autoExpAndGain1.png";
                 full_image_path_rad = path_to_polarized_cal_images  + "/rad_" + final_exposure_time_str + "_" + num2str(gain) ;
            end %full_image_path = path_to_polarized_cal_images  + "/pol_set_" + path_set + "_" + final_exposure_time_str + "_" + num2str(gain)
    end

    %read the image
    %fix its bits, so that it displays properly
    hyperspectralImageAE_cal = imread(char(full_image_path ));
    hyperspectralImageAEHigh_cal = bitshift(hyperspectralImageAE_cal, 16-12);
    hyperspectralImageAELow_cal= bitshift(hyperspectralImageAE_cal, 16-28);
    hyperspectralImageAE_cal=hyperspectralImageAEHigh_cal+hyperspectralImageAELow_cal;

   % imagesc(hyperspectralImageAE_cal);
    max_image_value = max(max(hyperspectralImageAE_cal));
    hyperspectralImageAE_cal = hyperspectralImageAE_cal(xrange,yrange);
    
    
    if max_image_value > saturation_cutoff_value
    
         exposure_time = exposure_time - 20000;
         rounded_exposure_time = round(exposure_time/20000) *20000;
         final_exposure_time = rounded_exposure_time +21;
         final_exposure_time_str = num2str(final_exposure_time);
 
    else
    % if there is no radiometer data, read a different file 
   %     full_image_path_rad = '/media/flagg/15AFC52042FBF9AC/Pika_L/Calibration Files/flats_100_images_IS/rad_440021_0/'
        datafile_835D = [char(full_image_path_rad)];
        datafile_835D_ramses =  dir((datafile_835D));
        if (length(datafile_835D_ramses) <3)
%             if exposure_time > 21000
%          exposure_time = exposure_time - 20000;
%          rounded_exposure_time = round(exposure_time/20000) *20000;
%          final_exposure_time = rounded_exposure_time +21;
%          final_exposure_time_str = num2str(final_exposure_time)
%             else
%           exposure_time = exposure_time + 20000;
%          rounded_exposure_time = round(exposure_time/20000) *20000;
%          final_exposure_time = rounded_exposure_time +21;
%          final_exposure_time_str = num2str(final_exposure_time)
%             end
             if bin == 3
            full_image_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/downloadedDatasets/calibration/hyperspectralImages/set_9/AE/autoExp1.png';
            full_exposure_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/downloadedDatasets/calibration/hyperspectralImages/set_9/AE/AutoExposure.txt';

             end

             if bin ==2 
            full_image_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/downloadedDatasets/calibration/hyperspectralImages/set_8/AE/autoExp1.png';
            full_exposure_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/downloadedDatasets/calibration/hyperspectralImages/set_8/AE/AutoExposure.txt';

             end

              if bin ==1
            full_image_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/downloadedDatasets/calibration/hyperspectralImages/set_7/AE/autoExp1.png';
            full_exposure_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/downloadedDatasets/calibration/hyperspectralImages/set_7/AE/AutoExposure.txt';
              end

              rad_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/downloadedDatasets/calibration/waterRadiometer/';
                hyperspectralImageAE_cal = imread(char(full_image_path ));
                hyperspectralImageAEHigh_cal = bitshift(hyperspectralImageAE_cal, 16-12);
                hyperspectralImageAELow_cal= bitshift(hyperspectralImageAE_cal, 16-28);
                hyperspectralImageAE_cal=hyperspectralImageAEHigh_cal+hyperspectralImageAELow_cal;
                hyperspectralImageAE_cal = hyperspectralImageAE_cal(xrange,yrange);
                
                exposure_time = str2num(fileread(full_exposure_path));
                     final_exposure_time = (exposure_time);
                     final_exposure_time_str = num2str(final_exposure_time);
                     
                datafile_835D = [char(rad_path)];
                datafile_835D_ramses =  dir((datafile_835D));    
                datafile_835D_ramses = extractfield(datafile_835D_ramses,'name');
                datafile_835D_ramses = datafile_835D_ramses(3);
                datafile_835D = [datafile_835D '/' datafile_835D_ramses{1}];
                [wavelength835D, spectrum835D(1:255)] = ramses_get_835B([CalFolder '835D_deployedRadiometer/'],'835D', datafile_835D);
                if isnan(wavelength835D)
                disp('nan wavelength 835D');
                end
                     
                     
                     
        else
            
            
            
    datafile_835D_ramses = extractfield(datafile_835D_ramses,'name');
    datafile_835D_ramses = datafile_835D_ramses(3);
    datafile_835D = [datafile_835D '/' datafile_835D_ramses{1}];
    [wavelength835D, spectrum835D(1:255)] = ramses_get_835B([CalFolder '835D_deployedRadiometer/'],'835D', datafile_835D);
    if isnan(wavelength835D)
        disp('nan wavelength 835D');
    end
            
            

        end
        
    end

        
        
  end
  
  
  
%   %try to get a radiance reading near in time from the radiometer
%     if select_gain_images == 0 && select_polarized_images == 0
%     full_image_path_rad = path_to_unpolarized_cal_images  + "/pol_rad_" + final_exposure_time_str + "_" + num2str(gain) 
%     else
%     %full_image_path = path_to_polarized_cal_images  + "/pol_set_" + path_set + "_" + final_exposure_time_str + "_" + num2str(gain)
%     end
% 
%         datafile_835D = [char(full_image_path_rad)];
%         datafile_835D_ramses =  dir((datafile_835D));
%         
% if length(datafile_835D_ramses) < 3;
%     disp('radiometer read error for calibration')
%     spectrum835D(1:255) = nan;
%     wavelength835D = [[305.422398662160,308.757055000290,312.092334241280,315.428225058750,318.764716126320,322.101796117610,325.439453706240,328.777677565830,332.116456370000,335.455778792370,338.795633506560,342.136009186190,345.476894504880,348.818278136250,352.160148753920,355.502495031510,358.845305642640,362.188569260930,365.532274560000,368.876410213470,372.220964894960,375.565927278090,378.911286036480,382.257029843750,385.603147373520,388.949627299410,392.296458295040,395.643629034030,398.991128190000,402.338944436570,405.687066447360,409.035482895990,412.384182456080,415.733153801250,419.082385605120,422.431866541310,425.781585283440,429.131530505130,432.481690880000,435.832055081670,439.182611783760,442.533349659890,445.884257383680,449.235323628750,452.586537068720,455.937886377210,459.289360227840,462.640947294230,465.992636250000,469.344415768770,472.696274524160,476.048201189790,479.400184439280,482.752212946250,486.104275384320,489.456360427110,492.808456748240,496.160553021330,499.512637920000,502.864700117870,506.216728288560,509.568711105690,512.920637242880,516.272495373750,519.624274171920,522.975962311010,526.327548464640,529.679021306430,533.030369510000,536.381581748970,539.732646696960,543.083553027590,546.434289414480,549.784844531250,553.135207051520,556.485365648910,559.835308997040,563.185025769530,566.534504640000,569.883734282070,573.232703369360,576.581400575490,579.929814574080,583.277934038750,586.625747643120,589.973244060810,593.320411965440,596.667240030630,600.013716930000,603.359831337170,606.705571925760,610.050927369390,613.395886341680,616.740437516250,620.084569566720,623.428271166710,626.771530989840,630.114337709730,633.456680000000,636.798546534270,640.139925986160,643.480807029290,646.821178337280,650.161028583750,653.500346442320,656.839120586610,660.177339690240,663.514992426830,666.852067470000,670.188553493370,673.524439170560,676.859713175190,680.194364180880,683.528380861250,686.861751889920,690.194465940510,693.526511686640,696.857877801930,700.188552960000,703.518525834470,706.847785098960,710.176319427090,713.504117492480,716.831167968750,720.157459529520,723.482980848410,726.807720599040,730.131667455030,733.454810090000,736.777137177570,740.098637391360,743.419299404990,746.739111892080,750.058063526250,753.376142981120,756.693338930310,760.009640047440,763.325035006130,766.639512480000,769.953061142670,773.265669667760,776.577326728890,779.888020999680,783.197741153750,786.506475864720,789.814213806210,793.120943651840,796.426654075230,799.731333750000,803.034971349770,806.337555548160,809.639075018790,812.939518435280,816.238874471250,819.537131800320,822.834279096110,826.130305032240,829.425198282330,832.718947520000,836.011541418870,839.302968652560,842.593217894690,845.882277818880,849.170137098750,852.456784407920,855.742208420010,859.026397808640,862.309341247430,865.591027410000,868.871444969970,872.150582600960,875.428428976590,878.704972770480,881.980202656250,885.254107307520,888.526675397910,891.797895601040,895.067756590530,898.336247040000,901.603355623070,904.869071013360,908.133381884490,911.396276910080,914.657744763750,917.917774119120,921.176353649810,924.433472029440,927.689117931630,930.943280030000,934.195946998170,937.447107509760,940.696750238390,943.944863857680,947.191437041250,950.436458462720,953.679916795710,956.921800713840,960.162098890730,963.400800000000,966.637892715270,969.873365710160,973.107207658290,976.339407233280,979.569953108750,982.798833958320,986.026038455610,989.251555274240,992.475373087830,995.697480570000,998.917866394370,1002.13651923456,1005.35342776419,1008.56858065688,1011.78196658625,1014.99357422592,1018.20339224951,1021.41140933064,1024.61761414293,1027.82199536000,1031.02454165547,1034.22524170296,1037.42408417609,1040.62105774848,1043.81615109375,1047.00935288552,1050.20065179741,1053.39003650304,1056.57749567603,1059.76301799000,1062.94659211857,1066.12820673536,1069.30785051399,1072.48551212808,1075.66118025125,1078.83484355712,1082.00649071931,1085.17611041144,1088.34369130713,1091.50922208000,1094.67269140367,1097.83408795176,1100.99340039789,1104.15061741568,1107.30572767875,1110.45871986072,1113.60958263521,1116.75830467584,1119.90487465623,1123.04928125000,1126.19151313077,1129.33155897216,1132.46940744779,1135.60504723128,1138.73846699625,1141.86965541632]];;
%     calibrated_image = zeros(size(image));
% else
%     datafile_835D_ramses = extractfield(datafile_835D_ramses,'name');
%     datafile_835D_ramses = datafile_835D_ramses(3);
%     datafile_835D = [datafile_835D '/' datafile_835D_ramses{1}];
%     [wavelength835D, spectrum835D(1:255)] = ramses_get_835B([CalFolder '835D_deployedRadiometer/'],'835D', datafile_835D);
%     if isnan(wavelength835D)
%         disp('nan wavelength 835D');
% 
% end
  
%        exposure_time_rad = exposure_time;
        
%read radiance values from the files
% while(length(datafile_835D_ramses) < 3)
%         if exposure_time < 250000
%          exposure_time_rad = exposure_time_rad + 20000;
%         else
%          exposure_time_rad = exposure_time_rad - 20000;
%         end
%         
%          rounded_exposure_time_rad = round(exposure_time_rad/20000) *20000;
%          final_exposure_time_rad = rounded_exposure_time_rad +21;
%          final_exposure_time_str_rad = num2str(final_exposure_time_rad);
%          
%     if select_gain_images == 0 && select_polarized_images == 0
%     full_image_path_rad = path_to_unpolarized_cal_images  + "/pol_rad_" + final_exposure_time_str_rad + "_" + num2str(gain) 
%     else
%     %full_image_path = path_to_polarized_cal_images  + "/pol_set_" + path_set + "_" + final_exposure_time_str + "_" + num2str(gain)
%     end
% 
%         datafile_835D = [char(full_image_path_rad)];
%         datafile_835D_ramses =  dir((datafile_835D));   
% 
% end


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
full_path_darks_calibration = '/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/darks/set_' + path_set + "_" + final_exposure_time_str_dark + "_" + num2str(gain) + "/AEG/autoExpAndGain1.png";
full_path_darks_image = '/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/darks/set_' + path_set + "_" + final_exposure_time_str_image + "_" + num2str(gain) + "/AEG/autoExpAndGain1.png";

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
                dark_img = dark_image_median;
               % figure; imagesc(dark_img);
                
%interpolate radiance values to camera wavelengths

interpolated_radiance_camera_wavelengths(:) = interp1(wavelength835D,spectrum835D(1:255),camera_wavelengths(bin,:));
camera_image_of_radiometer = ones(size(image)) .* (squeeze(interpolated_radiance_camera_wavelengths(1,1:size(image))))';

%need to subtract darks from hyperspectralImageAE_cal and from image, but
%need to select the closest exposure time darks
hyperspectralImageAE_cal = calibration_image_median;
calibration_image = camera_image_of_radiometer ./ double(hyperspectralImageAE_cal-dark_img) *   double(str2num(final_exposure_time_str));  %  double(440021);  double(str2num(final_exposure_time_str));
figure(9); imagesc(hyperspectralImageAE_cal - image); colorbar;
% size(image)
% size(dark_img)
% size(image_exposure_time)
% size(calibration_image)

calibrated_image = double(image-dark_img)/double(image_exposure_time).*double(calibration_image);
        
        
% output: camera_image_of_radiometer ./ double(all_images(set_num).spectrum_image_cal) * double(str2num(all_images(set_num).cal_exposure));
calibrated_image=calibrated_image;
end