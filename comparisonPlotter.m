%cd /home/flagg/Desktop/downloadedDatasets/code
%Fri April 6th 2018 at 15:34:23 = Wed Feb 24 20:40:23 (UTC) 2016
select_polarized_images = 0
select_gain_images = 0


%folders that can be changed
folderOfInterest = '/home/flagg/Dropbox (MIT)/PikaL_Dev/MVCO_datasets/';
CalFolder = '/home/flagg/Desktop/downloadedDatasets/code/Calibration Files/'; 
folderOfInterest_cal = '/home/flagg/Dropbox (MIT)/PikaL_Dev/downloadedDatasets/calibration';

%set folders
folders = dir(folderOfInterest); 
folders(1:2) = []; 


%%

%Comparison Plotter and calibration application for 1 file
% char array of setnum values
set_num_name_array = ["set_num_1";"set_num_2";"set_num_3";"set_num_4";"set_num_5";"set_num_6";"set_num_7";"set_num_8";"set_num_9";"set_num_10";"set_num_11";"set_num_12";]
set_num_name_array_numbers = ["set_1";"set_2";"set_3";"set_4";"set_5";"set_6";"set_7";"set_8";"set_9";"set_10";"set_11";"set_12";]

set_num_name_calibrated_array = ["set_num_1_cal";"set_num_2_cal";"set_num_3_cal";"set_num_4_cal";"set_num_5_cal";"set_num_6_cal";"set_num_7_cal";"set_num_8_cal";"set_num_9_cal";"set_num_10_cal";"set_num_11_cal";"set_num_12_cal";]

%%
%
% for(i=1:(length(folders)-2)) 
%     date{i} = folders(i).date
% end
%%
%create spatial mean of the usefule spectrum (for the first image) of
%camera data
%spectrum=mean(ImAll(316:915,528:1427,:),2);


%%
%create subplot of camera data
%figure(1); subplot(212) ; pcolor(datenum(date)-5/24-datenum([2018 2 26 0 0 0]),1:600,squeeze(spectrum)); shading flat; colorbar; caxis([0 20]); title('Camera Water (D) RAW')

%create subplot of water radiometer data for comparison
%figure(1); subplot(211); pcolor(datenum(date)-5/24-datenum([2018 2 26 0 0 0]),1:255,spectrum835D); shading flat; colorbar; caxis([0 20]);title('Water (D)  RAD')


%%
%Read in the radiometer data from the calibration file

datafile_835B_cal = [folderOfInterest_cal '/waterRadiometer'];
datafile_835B_ramses_cal =  dir(datafile_835B_cal);
datafile_835B_ramses_cal = extractfield(datafile_835B_ramses_cal,'name');
datafile_835B_ramses_cal = datafile_835B_ramses_cal(3);
datafile_835B_cal = [datafile_835B_cal '/' datafile_835B_ramses_cal{1}];

[wavelength835B_cal, spectrum835B_cal] = ramses_get_835B([CalFolder '835D_deployedRadiometer/'],'835D', datafile_835B_cal);
%spectrum_mean_cal = mean(hyperspectralImageAE_cal(316:915,:),2);
%figure(10); plot(wavelength835B_cal,spectrum835B_cal/max(spectrum835B_cal));hold on; plot(camera_wavelengths,spectrum_mean_cal/max(spectrum_mean_cal))
%%
% we have to do this for all 12 sets
%read in the calibration data for the first Image, pull out the useful
%section, scale by the exposure time (ignore dark reads for now)
if select_polarized_images == 1
    HyperImage = '/hyperspectralImagesPolarized/';
    pol = '_pol';
else
    HyperImage = '/hyperspectralImages/';
    pol= '';
end

if select_gain_images == 1
    gain = '_gain';
else
    gain = '';
end

%cal_image_path = [folderOfInterest_cal HyperImage char(set_num_name_array_numbers(set_num)) '/AE/autoExp1.png']


%reads in the calibration information for all 12 sets, depending on
%polarized or not polarized. in our case only set 7:9 matter, if they have
%poor calibrations this will be an issue. I should instead write a script
%that reads the correct calibration image and pulls out the correct data
%based on the exposure time (ignore gains for now)
for(set_num = 1:12)

hyperspectralImageAE_cal = imread([folderOfInterest_cal HyperImage char(set_num_name_array_numbers(set_num)) '/AE/autoExp1.png']);
hyperspectralImageAEHigh_cal = bitshift(hyperspectralImageAE_cal, 16-12);
hyperspectralImageAELow_cal= bitshift(hyperspectralImageAE_cal, 16-28);
hyperspectralImageAE_cal=hyperspectralImageAEHigh_cal+hyperspectralImageAELow_cal;
all_images(set_num).calibration = hyperspectralImageAE_cal;
all_images(set_num).cal_exposure = fileread([folderOfInterest_cal HyperImage char(set_num_name_array_numbers(set_num)) '/AE/AutoExposure.txt']);

if(size(hyperspectralImageAE_cal,2) == 1920 && size(hyperspectralImageAE_cal,1) == 1200)
   all_images(set_num).spectrum_image_cal  = hyperspectralImageAE_cal(316:915,528:1427); 
end

if(size(hyperspectralImageAE_cal,2) == 960 && size(hyperspectralImageAE_cal,1) == 600)
   all_images(set_num).spectrum_image_cal  = hyperspectralImageAE_cal(158:457,264:713); 
end

if(size(hyperspectralImageAE_cal,2) == 480 && size(hyperspectralImageAE_cal,1) == 300)
   all_images(set_num).spectrum_image_cal  = hyperspectralImageAE_cal(79:228,132:356); 
end

%exposure_cal(set_num,:) =  fileread([folderOfInterest_cal HyperImage char(set_num_name_array_numbers(set_num))  '/AE/AutoExposure.txt']);
figure(2); subplot(211); imagesc(hyperspectralImageAE_cal); colorbar; title(all_images(set_num).cal_exposure);
figure(2); subplot(212); plot(median(hyperspectralImageAE_cal,2));title(all_images(set_num).cal_exposure);
end
%read in the calibration image from folder

%pull out useful data section by looking at width of images (also figure
%out binning here)
%spectrum_image_cal = hyperspectralImageAE_cal(316:915,528:1427);
%figure(11); imagesc(spectrum_image_cal); colorbar;

%assign a wavelength to each pixel
%camera_wavelengths = [387.07, 388.06, 389.06, 390.05, 391.05, 392.04, 393.04, 394.03, 395.03, 396.03, 397.02, 398.02, 399.02, 400.02, 401.01, 402.01, 403.01, 404.01, 405.01, 406.0, 407.0, 408.0, 409.0, 410.0, 411.0, 412.0, 413.0, 414.0, 415.0, 416.0, 417.0, 418.01, 419.01, 420.01, 421.01, 422.01, 423.01, 424.02, 425.02, 426.02, 427.03, 428.03, 429.03, 430.04, 431.04, 432.05, 433.05, 434.06, 435.06, 436.07, 437.07, 438.08, 439.08, 440.09, 441.09, 442.1, 443.11, 444.11, 445.12, 446.13, 447.14, 448.14, 449.15, 450.16, 451.17, 452.18, 453.19, 454.2, 455.21, 456.21, 457.22, 458.23, 459.24, 460.25, 461.27, 462.28, 463.29, 464.3, 465.31, 466.32, 467.33, 468.35, 469.36, 470.37, 471.38, 472.4, 473.41, 474.42, 475.44, 476.45, 477.46, 478.48, 479.49, 480.51, 481.52, 482.54, 483.55, 484.57, 485.58, 486.6, 487.62, 488.63, 489.65, 490.67, 491.68, 492.7, 493.72, 494.74, 495.75, 496.77, 497.79, 498.81, 499.83, 500.85, 501.87, 502.89, 503.91, 504.93, 505.95, 506.97, 507.99, 509.01, 510.03, 511.05, 512.07, 513.09, 514.11, 515.14, 516.16, 517.18, 518.2, 519.23, 520.25, 521.27, 522.3, 523.32, 524.35, 525.37, 526.39, 527.42, 528.44, 529.47, 530.49, 531.52, 532.54, 533.57, 534.6, 535.62, 536.65, 537.68, 538.7, 539.73, 540.76, 541.79, 542.81, 543.84, 544.87, 545.9, 546.93, 547.96, 548.99, 550.02, 551.05, 552.08, 553.11, 554.14, 555.17, 556.2, 557.23, 558.26, 559.29, 560.32, 561.35, 562.39, 563.42, 564.45, 565.48, 566.52, 567.55, 568.58, 569.62, 570.65, 571.68, 572.72, 573.75, 574.79, 575.82, 576.86, 577.89, 578.93, 579.96, 581.0, 582.04, 583.07, 584.11, 585.14, 586.18, 587.22, 588.26, 589.29, 590.33, 591.37, 592.41, 593.45, 594.49, 595.52, 596.56, 597.6, 598.64, 599.68, 600.72, 601.76, 602.8, 603.84, 604.88, 605.93, 606.97, 608.01, 609.05, 610.09, 611.13, 612.18, 613.22, 614.26, 615.31, 616.35, 617.39, 618.44, 619.48, 620.52, 621.57, 622.61, 623.66, 624.7, 625.75, 626.79, 627.84, 628.88, 629.93, 630.98, 632.02, 633.07, 634.12, 635.16, 636.21, 637.26, 638.31, 639.36, 640.4, 641.45, 642.5, 643.55, 644.6, 645.65, 646.7, 647.75, 648.8, 649.85, 650.9, 651.95, 653.0, 654.05, 655.1, 656.15, 657.2, 658.26, 659.31, 660.36, 661.41, 662.47, 663.52, 664.57, 665.63, 666.68, 667.73, 668.79, 669.84, 670.9, 671.95, 673.01, 674.06, 675.12, 676.17, 677.23, 678.28, 679.34, 680.4, 681.45, 682.51, 683.57, 684.63, 685.68, 686.74, 687.8, 688.86, 689.92, 690.97, 692.03, 693.09, 694.15, 695.21, 696.27, 697.33, 698.39, 699.45, 700.51, 701.57, 702.63, 703.69, 704.76, 705.82, 706.88, 707.94, 709.0, 710.07, 711.13, 712.19, 713.26, 714.32, 715.38, 716.45, 717.51, 718.57, 719.64, 720.7, 721.77, 722.83, 723.9, 724.96, 726.03, 727.1, 728.16, 729.23, 730.3, 731.36, 732.43, 733.5, 734.56, 735.63, 736.7, 737.77, 738.84, 739.9, 740.97, 742.04, 743.11, 744.18, 745.25, 746.32, 747.39, 748.46, 749.53, 750.6, 751.67, 752.74, 753.82, 754.89, 755.96, 757.03, 758.1, 759.18, 760.25, 761.32, 762.4, 763.47, 764.54, 765.62, 766.69, 767.76, 768.84, 769.91, 770.99, 772.06, 773.14, 774.21, 775.29, 776.37, 777.44, 778.52, 779.59, 780.67, 781.75, 782.83, 783.9, 784.98, 786.06, 787.14, 788.22, 789.29, 790.37, 791.45, 792.53, 793.61, 794.69, 795.77, 796.85, 797.93, 799.01, 800.09, 801.17, 802.25, 803.33, 804.42, 805.5, 806.58, 807.66, 808.74, 809.83, 810.91, 811.99, 813.08, 814.16, 815.24, 816.33, 817.41, 818.5, 819.58, 820.67, 821.75, 822.84, 823.92, 825.01, 826.09, 827.18, 828.27, 829.35, 830.44, 831.53, 832.61, 833.7, 834.79, 835.88, 836.96, 838.05, 839.14, 840.23, 841.32, 842.41, 843.5, 844.59, 845.68, 846.77, 847.86, 848.95, 850.04, 851.13, 852.22, 853.31, 854.4, 855.49, 856.59, 857.68, 858.77, 859.86, 860.96, 862.05, 863.14, 864.24, 865.33, 866.42, 867.52, 868.61, 869.71, 870.8, 871.9, 872.99, 874.09, 875.18, 876.28, 877.37, 878.47, 879.57, 880.66, 881.76, 882.86, 883.96, 885.05, 886.15, 887.25, 888.35, 889.45, 890.54, 891.64, 892.74, 893.84, 894.94, 896.04, 897.14, 898.24, 899.34, 900.44, 901.54, 902.64, 903.74, 904.85, 905.95, 907.05, 908.15, 909.25, 910.36, 911.46, 912.56, 913.67, 914.77, 915.87, 916.98, 918.08, 919.19, 920.29, 921.39, 922.5, 923.6, 924.71, 925.82, 926.92, 928.03, 929.13, 930.24, 931.35, 932.45, 933.56, 934.67, 935.78, 936.88, 937.99, 939.1, 940.21, 941.32, 942.43, 943.53, 944.64, 945.75, 946.86, 947.97, 949.08, 950.19, 951.3, 952.41, 953.53, 954.64, 955.75, 956.86, 957.97, 959.08, 960.2, 961.31, 962.42, 963.53, 964.65, 965.76, 966.87, 967.99, 969.1, 970.22, 971.33, 972.45, 973.56, 974.68, 975.79, 976.91, 978.02, 979.14, 980.25, 981.37, 982.49, 983.6, 984.72, 985.84, 986.96, 988.07, 989.19, 990.31, 991.43, 992.55, 993.67, 994.78, 995.9, 997.02, 998.14, 999.26, 1000.38, 1001.5, 1002.62, 1003.74, 1004.87, 1005.99, 1007.11, 1008.23, 1009.35, 1010.47, 1011.6, 1012.72, 1013.84, 1014.96, 1016.09, 1017.21, 1018.33, 1019.46, 1020.58, 1\021.71];
%%
%coefficients = polyfit(

%%set the camera wavelengths, the radiometer wavelengths, and then
%%interpolate the radiometer wavelengths to the camera wavelengths
offset = 3.3;
for(i = 316:915)
%OG: camera_wavelengths(i-315) = 87.25381 + .92718 * i + 1.009192E-4*i^2;
%camera_wavelengths(i-315) = 87.25381-10 + .94018 * i + 1.96E-4*i^2;
camera_wavelengths(1,i-315) = 87.25381 + .92718 * (i+8) + 1.009192E-4*(i+8)^2; % their camera calibration with the 8 pixel offset input
%camera_wavelengths(1,i-315) = 87.25381 + .92718 * (i+8) + 1.009192E-4*(i+8)^2; % my camera calibration with the xenon mercury lamp
%camera_wavelengths(1,i-315) = p(3) + p(2) * (i-315) + p(1)*(i-315)^2; % my camera calibration with the xenon mercury lamp
%camera_wavelengths(1,i-315) = camera_wavelengths(1,i-315) + offset + i*(12.5-offset)/(315+364);
%camera_wavelengths(i-315) = 297.9 + .05086 * i + 9.433E-4*i^2;
end
for(i = 158:457)
camera_wavelengths(2,i-157) = (camera_wavelengths(1,(i*2)-315) +  camera_wavelengths(1,(i*2) - 315 +1))/2 -1 ;
end
for(i = 79:228)
camera_wavelengths(3,i-78) = (camera_wavelengths(1,(i*4)-315) +  camera_wavelengths(1,(i*4) - 315 +1) + camera_wavelengths(1,(i*4)-315 +2) +  camera_wavelengths(1,(i*4) - 315 +3))/4 -3 ;
end

% goood offset for blue : camera_wavelengths = camera_wavelengths +12.7-5;%
%camera_wavelengths = camera_wavelengths +12.7;
%camera_wavelengths = camera_wavelengths +12.7-5;
for(i=1:255)
D_wavelength_check(i) = 2.987550000E+02  + (3.333050000E+00)* (i+1) + (3.284410000E-04)* (i+1)^2 + (-1.887730000E-06)*(i+1)^3;
end
%D_wavelength_check = D_wavelength_check + 15.9;
% Interpolate radiometer radiance to pika L wavelengths
for( i=1:3)
interpolated_radiance_camera_wavelengths(i,:) = interp1(D_wavelength_check,spectrum835B_cal,camera_wavelengths(i,:));
%     for ( j = 1:size(spectrum8396,2))
%     interpolated_irrad_camera_wavelengths(i,j,:) = interp1(wavelength8396,spectrum8396(:,j),camera_wavelengths(i,:));
%     end
end
%%
figure(3); subplot(211); plot(camera_wavelengths(1,:),interpolated_radiance_camera_wavelengths(1,:)); xlim([300 1000]);
figure(3); subplot(212); plot(camera_wavelengths(2,1:300),mean(all_images(8).spectrum_image_cal,2));  xlim([300 1000]);



%
%check D wavelengths

%%
% divide the radiometer radiance by the digital numbers of the camera
for(set_num = [1, 4, 7, 10])
   
%creates a camera image sized copy of the radiance
camera_image_of_radiometer = ones(size(all_images(set_num).spectrum_image_cal)) .* (squeeze(interpolated_radiance_camera_wavelengths(1,:)))';

%the calibrated image is the radiometer divided by the image calibration
%multiplied by the exposure time  of the calibration file
all_images(set_num).calibrated_camera_image = camera_image_of_radiometer ./ double(all_images(set_num).spectrum_image_cal) * double(str2num(all_images(set_num).cal_exposure));
figure(4); subplot(211);imagesc(all_images(set_num).calibrated_camera_image); colorbar; title('calibrated camera image gains');caxis([0 .1E-6]);


end
%
%binned_cal_2 =softwareBin(all_images(1).spectrum_image_cal,2);
%binned_cal_4 =softwareBin(all_images(1).spectrum_image_cal,4);
%figure(20)
%imagesc(binned_cal_2)
%


for(set_num = [2, 5, 8, 11])
   
camera_image_of_radiometer = ones(size(all_images(set_num).spectrum_image_cal)) .* (squeeze(interpolated_radiance_camera_wavelengths(2,1:size(all_images(set_num).spectrum_image_cal,1))))';
all_images(set_num).calibrated_camera_image = camera_image_of_radiometer ./ double(all_images(set_num).spectrum_image_cal) * double(str2num(all_images(set_num).cal_exposure));
%all_images(set_num).calibrated_camera_image = camera_image_of_radiometer ./ double(binned_cal_2) * double(str2num(all_images(1).cal_exposure));
figure(4); subplot(211);imagesc(all_images(set_num).calibrated_camera_image); colorbar; title('calibrated camera image gains');caxis([0 .1E-6]);


end

for(set_num = [3, 6, 9, 12])
   
camera_image_of_radiometer = ones(size(all_images(set_num).spectrum_image_cal)) .* (squeeze(interpolated_radiance_camera_wavelengths(3,1:size(all_images(set_num).spectrum_image_cal,1))))';
all_images(set_num).calibrated_camera_image = camera_image_of_radiometer ./ double(all_images(set_num).spectrum_image_cal) * double(str2num(all_images(set_num).cal_exposure));
figure(4); subplot(211);imagesc(all_images(set_num).calibrated_camera_image); colorbar; title('calibrated camera image gains');caxis([0 .1E-6]);


end
%radiance_from_camera = calibrated_camera_image .* double(spectrum_image_cal) *500000;
%figure(4); subplot(212); plot(camera_wavelengths, radiance_from_camera); hold on; plot(camera_wavelengths, interpolated_radiance_camera_wavelengths)
% finally divide by exposure time (in microseconds)

%
%multiply every ImAll by the 2D radiance converted image
%calibratableImage = double(ImAll(316:915,528:1427,:));
%%
% clc
% for(set_num = [7])%[1,4,7,10])
%         clear calibratableImage;
% 
%     for (i =1:size(all_images,2))
%         inBetweenImage = double(all_images(i).(set_num_name_array(set_num) + pol+gain));%(316:915,528:1427,:));
%             if(size(inBetweenImage,2) >10)
%             calibratableImage(:,:,i) = inBetweenImage(316:915,528:1427);
%             if select_polarized_images == 1
%             calibratableImage(:,:,i) = calibratableImage(:,:,i) / exposure_Im_All_pol(i,set_num) .* all_images(set_num).calibrated_camera_image;
%             else
%             calibratableImage(:,:,i) = calibratableImage(:,:,i) / exposure_Im_All(i,set_num) .* all_images(set_num).calibrated_camera_image;
%             end
%             if select_gain_images == 1
%                 if select_polarized_images ==1
%                 calibratableImage(:,:,i) = calibratableImage(:,:,i)/ (10^(gain_Im_All_pol_gain(i,set_num)/20));       
%                 else
%                 calibratableImage(:,:,i) = calibratableImage(:,:,i)/ (10^(gain_Im_All(i,set_num)/20));       
%                 end
%             end
%             all_images(i).(set_num_name_calibrated_array(set_num)) =  calibratableImage(:,:,i);
%             end
%             
%     end
% 
% end
% 
% for(set_num = [8])%[2,5,8,11])
%     clear calibratableImage;
%     for (i =1:size(all_images,2))
%         inBetweenImage = double(all_images(i).(set_num_name_array(set_num)+pol+gain));%(316:915,528:1427,:));
%             if(size(inBetweenImage,2) >10)
%             calibratableImage(:,:,i) = inBetweenImage(158:457,264:713);
%             if select_polarized_images == 1
%             calibratableImage(:,:,i) = calibratableImage(:,:,i) / exposure_Im_All_pol(i,set_num) .* all_images(set_num).calibrated_camera_image;
%             else
%             calibratableImage(:,:,i) = calibratableImage(:,:,i) / exposure_Im_All(i,set_num) .* all_images(set_num).calibrated_camera_image;
%             end
%            if select_gain_images == 1
%                 if select_polarized_images ==1
%                 calibratableImage(:,:,i) = calibratableImage(:,:,i)/ (10^(gain_Im_All_pol_gain(i,set_num)/20));       
%                 else
%                 calibratableImage(:,:,i) = calibratableImage(:,:,i)/ (10^(gain_Im_All(i,set_num)/20));       
%                 end
%             end
%             all_images(i).(set_num_name_calibrated_array(set_num)) =  calibratableImage(:,:,i);
%             end
%             
%     end
% 
% end
% 
% for(set_num = 9)%[3,6,9,12])
%     clear calibratableImage;
%     for (i =1:size(all_images,2))
%         inBetweenImage = double(all_images(i).(set_num_name_array(set_num)+pol+gain));%(316:915,528:1427,:));
%             if(size(inBetweenImage,2) >10)
%             calibratableImage(:,:,i) = inBetweenImage(79:228,132:356);
%             if select_polarized_images == 1
%             calibratableImage(:,:,i) = calibratableImage(:,:,i) / exposure_Im_All_pol(i,set_num) .* all_images(set_num).calibrated_camera_image;
%             else
%             calibratableImage(:,:,i) = calibratableImage(:,:,i) / exposure_Im_All(i,set_num) .* all_images(set_num).calibrated_camera_image;
%             end
%             if select_gain_images == 1
%                 if select_polarized_images ==1
%                 calibratableImage(:,:,i) = calibratableImage(:,:,i)/ (10^(gain_Im_All_pol_gain(i,set_num)/20));       
%                 else
%                 calibratableImage(:,:,i) = calibratableImage(:,:,i)/ (10^(gain_Im_All(i,set_num)/20));       
%                 end
%             end
%             all_images(i).(set_num_name_calibrated_array(set_num)) =  calibratableImage(:,:,i);
%             end
%             
%     end
% 
% end
% %then take the mean spectrum and plot vs. water
% %177
% 

for(set_num = 7:9)%[3,6,9,12])
    for (i =1:size(all_images,2))
        set_num
        i
        all_images(i).(set_num_name_calibrated_array(set_num)) = calibrated_camera_image( (all_images(i).(set_num_name_array(set_num))),exposure_Im_All(i,set_num),0,'/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/rad_cal_long_0322_no_pol','/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/flats_835D_cleaned_0322_polarizerFlipped_all_polarized',select_polarized_images,select_gain_images,CalFolder,camera_wavelengths);
    end
end

figure(5); imagesc(calibratableImage(:,:,6)); colorbar
%%
%i=1; set_num = 7;
%calibrated_camera_image(calibratableImage(:,:,i),exposure_Im_All(i,set_num),0,'/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/rad_cal_long_0322_no_pol','/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/flats_835D_cleaned_0322_polarizerFlipped_all_polarized',select_polarized_images,select_gain_images)
%%
spectrum=nanmedian(calibratableImage(:,:,:),2);


%%
%create subplot of camera data
number =  10% 148 good, 84 good , 65ish, 177
set_num = 7%close(7)
clear Rrs
if(set_num == 2 || set_num == 5 || set_num == 8 ||set_num == 11)
figure(7); subplot(211); plot(D_wavelength_check, spectrum835D(:,number)); hold on;  plot(camera_wavelengths(2,1:size(median(all_images(number).(set_num_name_calibrated_array(set_num)),2),1)), median(all_images(number).(set_num_name_calibrated_array(set_num)),2)); 

            if select_polarized_images == 1
                if select_gain_images == 1
                title(['Exposure: ' num2str(exposure_Im_All_pol_gain(number,set_num)) ' gain: '  num2str(gain_Im_All_pol_gain(number,set_num)) ]) %*(500000/exposure_Im_All(number,set_num))^2)

                else    
                title(exposure_Im_All_pol(number,set_num)) %*(500000/exposure_Im_All(number,set_num))^2)
                end
            else
                if select_gain_images == 1
                title(['Exposure: ' num2str(exposure_Im_All_pol_gain(number,set_num)) ' gain: '  num2str(gain_Im_All(number,set_num)) ])
                else   
                title(exposure_Im_All(number,set_num)) %*(500000/exposure_Im_All(number,set_num))^2)
                end
            end
            for(i =1:size(all_images,2))
                if(size(median(all_images(i).(set_num_name_calibrated_array(set_num)),2)',2) >10)

Rrs(i,:) = squeeze(median(all_images(i).(set_num_name_calibrated_array(set_num)),2)') ./ interp1(wavelength8396,spectrum8396(:,i),camera_wavelengths(2,1:300));
                end
            end
            
figure(14); plot(camera_wavelengths(2,1:300),Rrs(number,:)); title('Rrs')

figure(7); subplot(212); plot(D_wavelength_check, spectrum835D(:,number)); hold on; plot(squeeze(camera_wavelengths(2,1:size(median(all_images(number).(set_num_name_calibrated_array(set_num)),2),1))), all_images(number).(set_num_name_calibrated_array(set_num)))
end

if(set_num == 3 || set_num == 6 || set_num == 9 ||set_num == 12)
figure(7); subplot(211); plot(D_wavelength_check, spectrum835D(:,number)); hold on;  plot(camera_wavelengths(3,1:size(median(all_images(number).(set_num_name_calibrated_array(set_num)),2),1)), median(all_images(number).(set_num_name_calibrated_array(set_num)),2));  %*(double(str2num(all_images(set_num).cal_exposure))/exposure_Im_All(number,set_num))^2
            if select_polarized_images == 1
                if select_gain_images == 1
                title(['Exposure: ' num2str(exposure_Im_All_pol_gain(number,set_num)) ' gain: '  num2str(gain_Im_All_pol_gain(number,set_num)) ]) %*(500000/exposure_Im_All(number,set_num))^2)

                else    
                title(exposure_Im_All_pol(number,set_num)) %*(500000/exposure_Im_All(number,set_num))^2)
                end
            else
                if select_gain_images == 1
                title(['Exposure: ' num2str(exposure_Im_All_pol_gain(number,set_num)) ' gain: '  num2str(gain_Im_All(number,set_num)) ])
                else   
                title(exposure_Im_All(number,set_num)) %*(500000/exposure_Im_All(number,set_num))^2)
                end
            end
            for(i =1:size(all_images,2))
                if(size(median(all_images(i).(set_num_name_calibrated_array(set_num)),2)',2) >10)
Rrs(i,:) = squeeze(median(all_images(i).(set_num_name_calibrated_array(set_num)),2)') ./ interp1(wavelength8396,spectrum8396(:,i),camera_wavelengths(3,1:150));
                end
            end
             figure(14); plot(camera_wavelengths(3,1:150),Rrs(number,:)); title('Rrs')
figure(7); subplot(212); plot(D_wavelength_check, spectrum835D(:,number)); hold on; plot(squeeze(camera_wavelengths(3,1:size(median(all_images(number).(set_num_name_calibrated_array(set_num)),2),1))), all_images(number).(set_num_name_calibrated_array(set_num)));grid on;
end

if(set_num == 1 || set_num == 4 || set_num == 7 ||set_num == 10)
figure(7); subplot(211); plot(D_wavelength_check, spectrum835D(:,number)); hold on;  plot(camera_wavelengths(1,1:size(median(all_images(number).(set_num_name_calibrated_array(set_num)),2),1)), median(all_images(number).(set_num_name_calibrated_array(set_num)),2)); 
            if select_polarized_images == 1
                if select_gain_images == 1
                title(['Exposure: ' num2str(exposure_Im_All_pol_gain(number,set_num)) ' gain: '  num2str(gain_Im_All_pol_gain(number,set_num)) ]) %*(500000/exposure_Im_All(number,set_num))^2)

                else    
                title(exposure_Im_All_pol(number,set_num)) %*(500000/exposure_Im_All(number,set_num))^2)
                end
            else
                if select_gain_images == 1
                title(['Exposure: ' num2str(exposure_Im_All_pol_gain(number,set_num)) ' gain: '  num2str(gain_Im_All(number,set_num)) ])
                else   
                title(exposure_Im_All(number,set_num)) %*(500000/exposure_Im_All(number,set_num))^2)
                end
            end
            
            for(i =1:size(all_images,2))
                if(size(median(all_images(i).(set_num_name_calibrated_array(set_num)),2)',2) >10)
Rrs(i,:) = squeeze(median(all_images(i).(set_num_name_calibrated_array(set_num)),2)') ./ interp1(wavelength8396,spectrum8396(:,i),camera_wavelengths(1,:));
                end
            end
             figure(14); plot(camera_wavelengths(1,1:600),Rrs(number,:)); title('Rrs')
figure(7); subplot(212); plot(D_wavelength_check, spectrum835D(:,number)); hold on; plot(squeeze(camera_wavelengths(1,1:size(median(all_images(number).(set_num_name_calibrated_array(set_num)),2),1))), all_images(number).(set_num_name_calibrated_array(set_num)))
end

clear spectrum
for(i =1:size(all_images,2))
                if(size(median(all_images(i).(set_num_name_calibrated_array(set_num)),2)',2) >10)
spectrum(i,:)=squeeze(median(all_images(i).(set_num_name_calibrated_array(set_num)),2)');


%Rrs(i,:) = squeeze(median(all_images(i).(set_num_name_calibrated_array(set_num)),2)') ./ interp1(wavelength8396,spectrum8396(:,i),camera_wavelengths(1,:));
                end
end

%%
%       figure(14)
%plot(camera_wavelengths(1,1:600),Rrs(number,:))

%plot(camera_wavelengths(2,1:300),Rrs(number,:))
%plot(camera_wavelengths(3,1:150),Rrs(number,:))

%%

%[folders.month] = [folders.name(1:10)]
offset_time_seconds = posixtime(datetime('2018-04-11 15:20:00')) - posixtime(datetime('2016-02-29 20:26:05'))
for(i=1:length(folders)-3)
[folders(i).month] = [folders(i).name(1:10)];
[folders(i).day] = [folders(i).name(12:19)];
month = folders(i).month ;
day = folders(i).day ;
date=[month, ' ',day];
[folders(i).timeS] = posixtime(datetime(date)) + offset_time_seconds;
[folders(i).date_corrected] = datetime(folders(i).timeS,'ConvertFrom','posixtime');
date_timestamp(i) = folders(i).date_corrected;
end
%%
clc
dates = datenum(date_timestamp(1:size(all_images,2)));
figure(6); subplot(411); pcolor(dates,wavelength835D,spectrum835D); shading flat; colorbar; title('Water (D) RAD'); axis ij; caxis([0 20]); yticks([0:50:200]); datetickzoom('x', 'HH:MM mm/dd '); ylabel('Wavelength (nm)');  yticks([ 400 600 800 1000]);
figure(6); subplot(412) ; pcolor(dates,camera_wavelengths(1,:),squeeze(spectrum)'); shading flat; colorbar; axis ij; title('Camera Water (D) RAD');  caxis([0 20]); datetickzoom('x', 'HH:MM mm/dd '); ylabel('Wavelength (nm)'); 
figure(6); subplot(413) ; pcolor(dates,wavelength8396,squeeze(spectrum8396)); shading flat; colorbar; axis ij; title('Ed');  caxis([0 1300   ]); datetickzoom('x', 'HH:MM mm/dd '); ylabel('Wavelength (nm)'); 
figure(6); subplot(414) ; pcolor(dates,camera_wavelengths(1,:),squeeze(Rrs)'); shading flat; colorbar; axis ij; title('Rrs');  caxis([0 .06   ]) ; datetickzoom('x', 'HH:MM mm/dd '); ylabel('Wavelength (nm)'); 


%%
close all
%794 no rad, 905, 906 
%test the calibrate camera code, to make sure it works in the event of a
%lack of radiometer code
i=720
camera_image_1 = calibrated_camera_image(all_images(i).set_num_7,exposure_Im_All(i,7),0,'/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/rad_cal_long_0322_no_pol','/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/flats_835D_cleaned_0322_polarizerFlipped_all_polarized',select_polarized_images,select_gain_images,CalFolder,camera_wavelengths);
figure(58); subplot(311); imagesc(camera_image_1); colorbar; caxis([0 20])
camera_image_2 = calibrated_camera_image(all_images(i).set_num_8,exposure_Im_All(i,8),0,'/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/rad_cal_long_0322_no_pol','/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/flats_835D_cleaned_0322_polarizerFlipped_all_polarized',select_polarized_images,select_gain_images,CalFolder,camera_wavelengths);
figure(58); subplot(312); imagesc(camera_image_2); colorbar; caxis([0 20])
camera_image_3 = calibrated_camera_image(all_images(i).set_num_9,exposure_Im_All(i,9),0,'/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/rad_cal_long_0322_no_pol','/home/flagg/Dropbox (MIT)/PikaL_Dev/Calibration Files/flats_835D_cleaned_0322_polarizerFlipped_all_polarized',select_polarized_images,select_gain_images,CalFolder,camera_wavelengths);
figure(58); subplot(313); imagesc(camera_image_3); colorbar; caxis([0 20])



figure(59); subplot(311); plot(camera_wavelengths(1,:), (camera_image_1)); xlim([400 900]); ylim([0 12])
figure(59); subplot(312); plot(camera_wavelengths(2,1:size(camera_image_2,1)), (camera_image_2)); xlim([400 900]) ; ylim([0 12])
figure(59); subplot(313); plot(camera_wavelengths(3,1:size(camera_image_3,1)), (camera_image_3)); xlim([400 900]); ylim([0 12])

%%
%save the requisite data for Jennie

%Jennie's path in Pika L Dev
path_to_jennie_data = '/home/flagg/Dropbox (MIT)/PikaL_Dev/data_for_Jennie/MVCO_rad_data.mat'
save(path_to_jennie_data, 'all_images', 'camera_wavelengths', 'date_timestamp', 'exposure_Im_All', 'interpolated_radiance_camera_wavelengths', 'spectrum835B', 'spectrum835D', 'spectrum8396', '-v7.3'  )






