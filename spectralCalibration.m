%read an image in
load('camera_wavelengths.mat')
load('spectralCal_1_1.mat')
path = '/home/flagg/Desktop/downloadedDatasets/spectralCalibrationImages/mercuryArgonSpectralCal.png'
%path = '/home/flagg/Desktop/downloadedDatasets/spectralCalibrationImages/Image__2016-02-13__10-53-23.png'
hyperspectralImageAE = imread(path);
hyperspectralImageAE_OG = hyperspectralImageAE;
hyperspectralImageAEHigh=bitshift(hyperspectralImageAE, 16-12);
hyperspectralImageAElow = bitshift(hyperspectralImageAE, 16-28);
hyperspectralImageAE =hyperspectralImageAEHigh+ hyperspectralImageAElow;
hyperspectralImageAE = hyperspectralImageAE;


%output_image_comparison = png_value_convertor(path);
%figure(29); subplot(211); ( imagesc(hyperspectralImageAE)); subplot(212); imagesc( (output_image_comparison));


figure(30); subplot(411); imagesc(hyperspectralImageAE)

range_spatial = 250:356;
subplot(412); imagesc(hyperspectralImageAE(316:915,range_spatial))

subplot(413); plot(mean(hyperspectralImageAE(316:915,range_spatial),2)/(max(mean(hyperspectralImageAE(316:915,range_spatial),2))))

subplot(414); plot(camera_wavelengths(1,:),mean(hyperspectralImageAE(316:915,range_spatial),2)/(max(mean(hyperspectralImageAE(316:915,range_spatial),2))));% xlim([400 1000])

figure(300); imagesc(hyperspectralImageAE(1010:1200,300:340))



load('camera_wavelengths.mat')
radiometerPath =  '/home/flagg/Desktop/downloadedDatasets/spectralCalibrationImages/radiometericArgon_835D/rad__'
datafile_835B_cal = radiometerPath;
datafile_835B_ramses_cal =  dir(datafile_835B_cal)
datafile_835B_ramses_cal = extractfield(datafile_835B_ramses_cal,'name')
datafile_835B_ramses_cal = datafile_835B_ramses_cal(3);
set_num_name_calibrated_array = ["set_num_1_cal";"set_num_2_cal";"set_num_3_cal";"set_num_4_cal";"set_num_5_cal";"set_num_6_cal";"set_num_7_cal";"set_num_8_cal";"set_num_9_cal";"set_num_10_cal";"set_num_11_cal";"set_num_12_cal";]


%create spatial mean of the usefule spectrum (for the first image) of

datafile_835B_cal = [datafile_835B_cal '/' datafile_835B_ramses_cal{1}]

CalFolder = '/home/flagg/Desktop/downloadedDatasets/code/Calibration Files/' %'/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/spectralCalibrationImages/radiometericArgon_835D/rad__'
[wavelength835B_cal, spectrum835B_cal] = ramses_get_835B([CalFolder '835D_deployedRadiometer/'],'835D', datafile_835B_cal);

for(i = 316:915)
%OG: camera_wavelengths(i-315) = 87.25381 + .92718 * i + 1.009192E-4*i^2;
%camera_wavelengths(i-315) = 87.25381-10 + .94018 * i + 1.96E-4*i^2;
%camera_wavelengths(1,i-315) = 87.25381 + .92718 * (i+8) + 1.009192E-4*(i+8)^2; % their camera calibration with the 8 pixel offset input
camera_wavelengths(1,i-315) = 87.25381 + .92718 * (i+8) + 1.009192E-4*(i+8)^2; % my camera calibration with the xenon mercury lamp
%camera_wavelengths(1,i-315) = p(3) + p(2) * (i-315) + p(1)*(i-315)^2; % my camera calibration with the xenon mercury lamp
%camera_wavelengths(1,i-315) = camera_wavelengths(1,i-315) + offset + i*(12.5-offset)/(315+364);
%camera_wavelengths(i-315) = 297.9 + .05086 * i + 9.433E-4*i^2;
end

figure(31); subplot(211);plot(camera_wavelengths(1,:),mean(hyperspectralImageAE(316:915,132:356),2)/(max(mean(hyperspectralImageAE(316:915,132:356),2)))); xlim([400 1000])
title('pre spectral calibration')
figure(31); subplot(211); hold on;  plot(wavelength835B_cal,spectrum835B_cal/max(spectrum835B_cal)); xlim([400 1000]);

p= polyfit(spectralCal_1_1(:,1),spectralCal_1_1(:,2),2)
camera_wavelengths_recalibrated = polyval(p,1:600)

figure(31); subplot(212);plot(camera_wavelengths_recalibrated,mean(hyperspectralImageAE(316:915,132:356),2)/(max(mean(hyperspectralImageAE(316:915,132:356),2)))); xlim([400 1000])
title('post spectral calibration')
figure(31); subplot(212);hold on;  plot(wavelength835B_cal,spectrum835B_cal/max(spectrum835B_cal)); xlim([400 1000]);


%%
figure(33);plot(mean(hyperspectralImageAE(316:1200,132:356),2)/(max(mean(hyperspectralImageAE(316:915,132:356),2))))

%%
figure(34); plot((hyperspectralImageAE(1:1200,371)))



%%
%create a spectral calibration for the second order data
red_image_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/spectralCalibrationImages/Image__2016-02-13__10-53-23_RED.png'
red_image = png_value_convertor(red_image_path);
figure(35); imagesc(red_image); colorbar

green_image_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/spectralCalibrationImages/Image__2016-02-13__11-58-32_GREEN.png'
green_image = imread(green_image_path);
figure(36); imagesc(green_image); colorbar

green_image_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/spectralCalibrationImages/Image__2016-02-13__11-58-32_GREEN.png'
green_image = imread(green_image_path);

image_450_path = '/home/flagg/Desktop/downloadedDatasets/spectralCalibrationImages/Image__2016-02-11__16-47-44_450_nm_filter.png'
image_450 = imread(image_450_path);

image_500_path = '/home/flagg/Desktop/downloadedDatasets/spectralCalibrationImages/Image__2016-02-11__16-51-37_500nm_filter.png'
image_500 = imread(image_500_path);

image_550_path = '/home/flagg/Desktop/downloadedDatasets/spectralCalibrationImages/Image__2016-02-11__16-53-43_550_nm_filter.png'
image_550 = imread(image_550_path);

image_580_path = '/home/flagg/Desktop/downloadedDatasets/spectralCalibrationImages/Image__2016-02-11__16-58-45_580_nm_filter.png'
image_580 = imread(image_580_path);

image_600_path = '/home/flagg/Desktop/downloadedDatasets/spectralCalibrationImages/Image__2016-02-11__17-04-28_600nm_filter.png'
image_600 = imread(image_600_path);

image_630_path = '/home/flagg/Desktop/downloadedDatasets/spectralCalibrationImages/Image__2016-02-11__17-07-42_630nm_filter.png'
image_630 = imread(image_630_path);

image_light_source_path = '/home/flagg/Desktop/downloadedDatasets/spectralCalibrationImages/Image__2016-02-11__17-12-36_camera_light_source.png'
image_source = imread(image_light_source_path);

image_670_path = '/home/flagg/Desktop/downloadedDatasets/spectralCalibrationImages/Image__2016-02-11__16-44-35_670_nm_filter.png'
image_670 = imread(image_670_path);

image_700_path = '/home/flagg/Desktop/downloadedDatasets/spectralCalibrationImages/Image__2016-02-11__16-33-55_700nm_filter.png'
image_700 = imread(image_700_path);


figure(37); subplot(331); imagesc( image_500); title('500 nm Image')
figure(37); subplot(332); imagesc( green_image); title('Green laser Image (532 nm)')
figure(37); subplot(333); imagesc(  image_550); title('550 nm Image')
figure(37); subplot(334); imagesc( image_580); title('580 nm Image')
figure(37); subplot(335); imagesc(  image_600); title('600 nm Image')
figure(37); subplot(336); imagesc(  image_630); title('630 nm Image')
figure(37); subplot(338); imagesc( image_670); title('670 nm Image')
figure(37); subplot(339); imagesc( image_700); title('700 nm Image')
figure(37); subplot(337); imagesc(  red_image); title('Red Laser Image (likely 650 nm)')

summed_image = mean(image_500,2) + mean(image_550,2) + mean(image_580,2) + mean(image_600,2) + mean(image_630,2) + mean(image_670,2) %+ mean(image_700,2)% + mean(red_image,2) + mean(green_image,2)

figure(41); subplot(211); imagesc(image_source); subplot(212); plot(camera_wavelengths(1,:),image_source(315:914,1424))
%%

%1 is wavelength, 2 is the pixel number?
load('second_order_spectral_calibration.mat')
p= polyfit(second_order_spectral_calibration(:,2),second_order_spectral_calibration(:,1),2)
camera_wavelengths_recalibrated_second_order_1216= polyval(p,1:1216)
camera_wavelengths_recalibrated_second_order_1200= camera_wavelengths_recalibrated_second_order_1216(9:1208)


figure(38); plot(camera_wavelengths_recalibrated_second_order_1200, mean(green_image(:,144:152),2))
figure(38); plot(camera_wavelengths_recalibrated_second_order_1200, mean(green_image(:,144:152),2))







figure(39); subplot(331); imagesc( 1:size(image_500,2),camera_wavelengths_recalibrated_second_order_1216, image_500); title('500 nm Image')
figure(39); subplot(332); imagesc( 1:size(green_image,2),camera_wavelengths_recalibrated_second_order_1200,green_image); title('Green laser Image (532 nm)')
figure(39); subplot(333); imagesc( 1:size(image_500,2),camera_wavelengths_recalibrated_second_order_1216, image_550); title('550 nm Image')
figure(39); subplot(334); imagesc( 1:size(image_500,2),camera_wavelengths_recalibrated_second_order_1216, image_580); title('580 nm Image')
figure(39); subplot(335); imagesc( 1:size(image_500,2),camera_wavelengths_recalibrated_second_order_1216, image_600); title('600 nm Image')
figure(39); subplot(336); imagesc( 1:size(image_500,2),camera_wavelengths_recalibrated_second_order_1216, image_630); title('630 nm Image')
figure(39); subplot(338); imagesc( 1:size(image_500,2),camera_wavelengths_recalibrated_second_order_1216, image_670); title('670 nm Image')
figure(39); subplot(339); imagesc( 1:size(image_500,2),camera_wavelengths_recalibrated_second_order_1216, image_700); title('700 nm Image')
figure(39); subplot(337); imagesc( 1:size(red_image,2),camera_wavelengths_recalibrated_second_order_1200, red_image); title('Red Laser Image (likely 650 nm)')

figure(40); subplot(211); plot(camera_wavelengths_recalibrated_second_order_1216)
figure(40); subplot(212); plot(camera_wavelengths_recalibrated_second_order_1200)




%%
%apply it to a real image 
base_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/MVCO_datasets/2016-03-21T22:30:04+00:00' % non-saturated
%base_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/MVCO_datasets/2016-03-21T23:15:15+00:00' %saturated
base_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/MVCO_datasets/2016-04-19T23:45:06+00:00'
%base_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/MVCO_datasets/2016-03-20T23:00:13+00:00'
base_path = '/home/flagg/Dropbox (MIT)/PikaL_Dev/MVCO_datasets/2016-03-16T17:00:08+00:00'
real_image_path = [base_path  '/hyperspectralImages/set_1/AE/autoExp1.png']%'/home/flagg/Dropbox (MIT)/PikaL_Dev/MVCO_datasets/2016-03-21T21:45:08+00:00/hyperspectralImages/set_1/AE/autoExp1.png'
real_world_mvco_image = png_value_convertor(real_image_path);
spatial_pixel = 850:880
spatial_pixel = 500:570%620:670 ;%850:880
spatial_pixel_2 = 1410:1470
xlim_val = [575 700];

figure(410); imagesc(real_world_mvco_image)
figure(42); subplot(211); plot(camera_wavelengths(1,100:300),    mean(real_world_mvco_image(415:615,spatial_pixel),2)); xlim(xlim_val); title('Camera Image MVCO: First Order'); ylabel('digital number'); grid on; hold on; plot(camera_wavelengths_recalibrated_second_order_1200(900:1200),    mean(real_world_mvco_image(900:1200,spatial_pixel),2).*first_order_divided_by_second_order'); legend('First Order Camera', 'Second Order Scaled');subplot(212); plot(camera_wavelengths_recalibrated_second_order_1200(900:1200),    mean(real_world_mvco_image(900:1200,spatial_pixel),2)'); xlim(xlim_val);title('Camera Image MVCO: Second Order'); grid on
xlabel('wavelength (nm)'); ylabel('digital number')

radiometerPath =  [base_path  '/skyRadiometer']%'/home/flagg/Dropbox (MIT)/PikaL_Dev/MVCO_datasets/2016-03-21T22:30:04+00:00/irradiometer'
datafile_835B_cal = radiometerPath;
datafile_835B_ramses_cal =  dir(datafile_835B_cal)
datafile_835B_ramses_cal = extractfield(datafile_835B_ramses_cal,'name')
datafile_835B_ramses_cal = datafile_835B_ramses_cal(3);
set_num_name_calibrated_array = ["set_num_1_cal";"set_num_2_cal";"set_num_3_cal";"set_num_4_cal";"set_num_5_cal";"set_num_6_cal";"set_num_7_cal";"set_num_8_cal";"set_num_9_cal";"set_num_10_cal";"set_num_11_cal";"set_num_12_cal";]


%create spatial mean of the usefule spectrum (for the first image) of

datafile_835B_cal = [datafile_835B_cal '/' datafile_835B_ramses_cal{1}]

CalFolder = '/home/flagg/Desktop/downloadedDatasets/code/Calibration Files/' %'/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/spectralCalibrationImages/radiometericArgon_835D/rad__'
[wavelength835B_cal, spectrum835B_cal] = ramses_get_835B([CalFolder '835B_deployedRadiometer/'],'835B', datafile_835B_cal);



radiometerPath= [base_path  '/irradiometer']
datafile_8396_cal = radiometerPath;
datafile_8396_ramses_cal =  dir(datafile_8396_cal)
datafile_8396_ramses_cal = extractfield(datafile_8396_ramses_cal,'name')
datafile_8396_ramses_cal = datafile_8396_ramses_cal(3);
datafile_8396_cal = [datafile_8396_cal '/' datafile_8396_ramses_cal{1}]
[wavelength8396_cal, spectrum8396_cal] = ramses_get_irradiometer([CalFolder '8396_deployedIrradiometer/'],'8396', datafile_8396_cal);

figure(420); subplot(311);plot(wavelength8396_cal, spectrum8396_cal);xlim([xlim_val]); title('Irradiance')
figure(420); subplot(312);plot(wavelength835B_cal, spectrum835B_cal);xlim([xlim_val]); title('sky radiance')

radiometerPath= [base_path  '/waterRadiometer']
datafile_835D_cal = radiometerPath;
datafile_835D_ramses_cal =  dir(datafile_835D_cal)
datafile_835D_ramses_cal = extractfield(datafile_835D_ramses_cal,'name')
datafile_835D_ramses_cal = datafile_835D_ramses_cal(3);
datafile_835D_cal = [datafile_835D_cal '/' datafile_835D_ramses_cal{1}]
[wavelength835D_cal, spectrum835D_cal] = ramses_get_835B([CalFolder '835D_deployedRadiometer/'],'835D', datafile_835D_cal);
figure(420); subplot(313);plot(wavelength835D_cal, spectrum835D_cal);xlim([xlim_val]); title('water radiance')

figure(421);  plot(camera_wavelengths_recalibrated_second_order_1200(900:1200), second_order_rad_correction.* mean(real_world_mvco_image(900:1200,spatial_pixel),2)); hold on; plot(camera_wavelengths(1,100:300),    mean(real_world_mvco_image(415:615,spatial_pixel),2).*first_order_rad_correction); hold on; plot(wavelength835D_cal, spectrum835D_cal); 
figure(422);   plot(camera_wavelengths(1,100:300),    mean(real_world_mvco_image(415:615,spatial_pixel_2),2).*first_order_rad_correction);  hold on; plot(camera_wavelengths(1,100:300),    mean(real_world_mvco_image(415:615,spatial_pixel),2).*first_order_rad_correction); hold on; plot(wavelength835D_cal, spectrum835D_cal); 



hold on; plot(wavelength835B_cal, spectrum835B_cal); xlim([xlim_val]); title('Radiometrically corrected first and second orders:saturated first order'); legend('Second order', 'First Order', 'Radiometer (WLR)', 'Sky Radiance'); xlabel('wavelength (nm)'); ylabel('Radiance')
   sky_rad_etimate = (mean(real_world_mvco_image(315:914,spatial_pixel),2).*first_order_rad_correction_full -  mean(real_world_mvco_image(315:914,spatial_pixel_2),2).*first_order_rad_correction_full)/(.0408 - .025);
figure(425);   plot(camera_wavelengths(1,:),    sky_rad_etimate);  
hold on; plot(wavelength835B_cal, spectrum835B_cal);  title('Radiometrically corrected first and second orders:saturated first order'); legend('Sky rad estimate ', 'Radiometer (WLR)'); xlabel('wavelength (nm)'); ylabel('Radiance')


%%
%here I try to pull 2 calibration images first/second orders to see how
%they compare
figure(44)
imagesc(png_value_convertor('/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/spectralCalibrationImages/Image__2016-02-13__09-01-09.png'))

figure(45)
second_order = png_value_convertor('/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/flats_835D_cleaned_0322_polarizerFlipped_all_polarized/pol_set_1_480021_0/AEG/autoExpAndGain1.png');
first_order = png_value_convertor('/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/flats_835D_cleaned_0322_polarizerFlipped_all_polarized/pol_set_1_480021_0/AEG/autoExpAndGain1.png');

subplot(211)
imagesc(second_order)
subplot(212)
imagesc(first_order)

figure(46)
plot(camera_wavelengths(1,100:300),    mean(first_order(415:615,spatial_pixel),2)/1200); hold on; plot(camera_wavelengths_recalibrated_second_order_1200(900:1200),    mean(second_order(900:1200,spatial_pixel),2)/40); xlim(xlim_val)

first_order_interpolated_to_second_order_wavelengths = interp1(camera_wavelengths(1,100:300),mean(first_order(415:615,spatial_pixel),2),camera_wavelengths_recalibrated_second_order_1200(900:1200))

figure(47)
plot(camera_wavelengths_recalibrated_second_order_1200(900:1200),   first_order_interpolated_to_second_order_wavelengths/1200); hold on; plot(camera_wavelengths_recalibrated_second_order_1200(900:1200),    mean(second_order(900:1200,spatial_pixel),2)/40); xlim(xlim_val)

first_order_divided_by_second_order = first_order_interpolated_to_second_order_wavelengths ./  mean(second_order(900:1200,spatial_pixel),2)'

radiometerPath= '/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/flats_835D_cleaned_0322_polarizerFlipped_all_polarized/pol_rad_480021_0'
datafile_835D_cal = radiometerPath;
datafile_835D_ramses_cal =  dir(datafile_835D_cal)
datafile_835D_ramses_cal = extractfield(datafile_835D_ramses_cal,'name')
datafile_835D_ramses_cal = datafile_835D_ramses_cal(3);
datafile_835D_cal = [datafile_835D_cal '/' datafile_835D_ramses_cal{1}]
[wavelength835D_cal, spectrum835D_cal] = ramses_get_835B([CalFolder '835D_deployedRadiometer/'],'835D', datafile_835D_cal);

rad_to_second_order_wavelengths =  interp1(wavelength835D_cal,spectrum835D_cal,camera_wavelengths_recalibrated_second_order_1200(900:1200))
rad_to_first_order_wavelengths = interp1(wavelength835D_cal,spectrum835D_cal,camera_wavelengths(1,100:300))

first_order_rad_correction = rad_to_first_order_wavelengths'./ mean(first_order(415:615,spatial_pixel),2) 
second_order_rad_correction = rad_to_second_order_wavelengths'./ mean(second_order(900:1200,spatial_pixel),2)

%%

rad_to_first_order_wavelengths_full = interp1(wavelength835D_cal,spectrum835D_cal,camera_wavelengths(1,:))

first_order_rad_correction_full = rad_to_first_order_wavelengths_full'./ mean(first_order(315:914,spatial_pixel),2) 
