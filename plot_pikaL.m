%%%%%%%%%%%%%%%%%%%%%%
%current data compiler use from this base
%/home/flagg/Dropbox (MIT)/PikaL_Dev/CodeBase/matlab/currentDataCompiler_03_16_18
%/home/flagg/Desktop/downloadedDatasets/code
%
% This script is intended to be used to view all of the data simultaneously
% to get a holistic view of the data
% It can be used in the field for troubleshooting 
clc
cd '/home/flagg/Dropbox (MIT)/PikaL_Dev/CodeBase/matlab/current_data_compiler'

folderOfInterest = '/home/flagg/Dropbox (MIT)/PikaL_Dev/MVCO_datasets/';
HyperImage = '/hyperspectralImages/';
CalFolder = '/home/flagg/Desktop/downloadedDatasets/code/Calibration Files/'; 
folders = dir(folderOfInterest); 
folders(1:2) = []; 

%%

for k = 2200:length(folders)-2; 
    set_num = 1;
    fget = dir([folderOfInterest folders(k).name HyperImage 'set_' num2str(set_num) '/AE/']);
if length(fget) < 3;
    ImAll(:,:,k) = nan;
    continue
end


%read a single file 
dateAndTime = folders(k).name;
%folderOfInterest = [datasetLocalFolderDatasets dateAndTime '/']

%read contents from folder of interest, and graph it for interpretation
upwardFacingCamera = imread([folderOfInterest folders(k).name '/upwards.jpg']);
downwardFacingCamera = imread([folderOfInterest folders(k).name '/downwards.jpg']);

temperature = fileread([folderOfInterest folders(k).name '/temperatureLog']);
temperature = temperature(end-5:end-1)

figure(1)
subplot(331)
imagesc(upwardFacingCamera)
title(['Upward: temperature: ' temperature])

subplot(332)
imagesc(downwardFacingCamera)
title(['downWard: date and time: ' dateAndTime])


hyperspectralImageAE = imread([folderOfInterest folders(k).name HyperImage 'set_' num2str(set_num) '/AE/autoExp1.png']);
hyperspectralImageAE_OG = hyperspectralImageAE;
hyperspectralImageAEHigh=bitshift(hyperspectralImageAE, 16-12);
hyperspectralImageAElow = bitshift(hyperspectralImageAE, 16-28);

hyperspectralImageAE =hyperspectralImageAEHigh+ hyperspectralImageAElow;
hyperspectralImageAE = hyperspectralImageAE;

%ImAll(:,:,k) = hyperspectralImageAE; 
exposure = fileread([folderOfInterest folders(k).name HyperImage 'set_' num2str(set_num) '/AE/AutoExposure.txt']);

%hyperspectral camera images
figure(1)
subplot(334)
%hyperspectralImageAE_OG = bitshift(hyperspectralImageAE_OG,4);
imagesc(hyperspectralImageAE)
bitOfPoint4 = bitget(hyperspectralImageAE_OG(1,1025),16:-1:1,'uint16');
bitOfPoint7 = bitget(hyperspectralImageAE_OG(1,7),16:-1:1,'uint16');


title(['Set #: ' num2str(set_num) ' AutoExposure: ' exposure ])
colorbar

hyperspectralImageAEG = imread([folderOfInterest folders(k).name  HyperImage 'set_' num2str(set_num) '/AEG/autoExpAndGain1.png']);
hyperspectralImageAEGHigh = bitshift(hyperspectralImageAEG, 16-12);
hyperspectralImageAEGLow = bitshift(hyperspectralImageAEG, 16-28);
hyperspectralImageAEG=hyperspectralImageAEGHigh+hyperspectralImageAEGLow;
exposure = fileread([folderOfInterest folders(k).name  HyperImage 'set_' num2str(set_num) '/AEG/AutoExposure.txt']);
gain = fileread([folderOfInterest folders(k).name  HyperImage 'set_' num2str(set_num) '/AEG/AutoGain.txt']);

subplot(335)
imagesc(hyperspectralImageAEG)
title(['Set #: ' num2str(set_num) ' AutoExpAndGain 1, exposure: ' exposure 'gain: ' gain])
colorbar
    start=43;
spectrum=mean(hyperspectralImageAE(316:915,528:1427),2);
    
subplot(336)
plot(1:600,spectrum)
    title('Average spectrum from pixel 0 to 600')
    
imagefolders = dir([folderOfInterest folders(k).name  HyperImage]);
imagefolders(1:2) = []; 
dummy1 = 1; dummy2 = 1;

for i = 1:length(imagefolders); 
    d = dir([folderOfInterest folders(k).name  HyperImage imagefolders(i).name]);
    d(1:2) = []; 
    d1 = dir([folderOfInterest folders(k).name  HyperImage imagefolders(i).name '/' d(1).name]);
    d1(1:2) = []; 
    d1 = struct2cell(d1); 
    if sum(strcmp(d1(1,:),'autoExp1.png')) == 1;
        I1(dummy1) = i;
        dummy1 = dummy1+1;
    end
    d2 = dir([folderOfInterest folders(k).name  HyperImage imagefolders(i).name '/' d(2).name]);
    d2(1:2) = []; 
    d2 = struct2cell(d2); 
    if sum(strcmp(d2(1,:),'autoExpAndGain1.png')) == 1;
        I2(dummy2) = i;
        dummy2 = dummy2+1;
    end
end

Iplot = intersect(I1,I2); 
clear I1 I2

for i = 1:length(Iplot);
  set_num = i;
hyperspectralImageAE = imread([folderOfInterest folders(k).name  HyperImage imagefolders(Iplot(i)).name '/AE/autoExp1.png']);
hyperspectralImageAEHigh = bitshift(hyperspectralImageAE, 16-12);
hyperspectralImageAELow= bitshift(hyperspectralImageAE, 16-28);
hyperspectralImageAE=hyperspectralImageAEHigh+hyperspectralImageAELow;
exposure = fileread([folderOfInterest folders(k).name  HyperImage imagefolders(Iplot(i)).name '/AE/AutoExposure.txt']);

binning=floor(1936/(size(hyperspectralImageAE,2)));
    figure(2)
subplot(7,4,(i-1)*2+1)
xoffset=(528/binning);
xwidth=900/binning;
xendpoint=xoffset +xwidth;
xrange=xoffset:(xoffset+xwidth-1);


yoffset=(316/binning);
ywidth=600/binning;
yendpoint=yoffset +ywidth;
yrange=yoffset:(yoffset+ywidth-1);

imagesc(hyperspectralImageAE(yrange,xrange))
title(['Set #: ' imagefolders(Iplot(i)).name ' AutoExposure: ' exposure ])
colorbar

hyperspectralImageAEG = imread([folderOfInterest folders(k).name  HyperImage imagefolders(Iplot(i)).name '/AEG/autoExpAndGain1.png']);
hyperspectralImageAEG_High = bitshift(hyperspectralImageAEG, 16-12);

hyperspectralImageAEG_low = bitshift(hyperspectralImageAEG, -12);

hyperspectralImageAEG =hyperspectralImageAEG_High + hyperspectralImageAEG_low;
exposure = fileread([folderOfInterest folders(k).name  HyperImage imagefolders(Iplot(i)).name '/AEG/AutoExposure.txt']);
gain = fileread([folderOfInterest folders(k).name  HyperImage imagefolders(Iplot(i)).name '/AEG/AutoGain.txt']);
start+(i-1)+2;
subplot(7,4,(i-1)*2+2)
imagesc(hyperspectralImageAEG(yrange,xrange))
title(['Set #: ' imagefolders(Iplot(i)).name ' AutoExpAndGain 1, exposure: ' exposure 'gain: ' gain])
colorbar
end
%

%%%%%%%%%
% 835B
datafile_835B = [folderOfInterest folders(k).name '/skyRadiometer'];
datafile_835B_ramses =  dir(datafile_835B);
datafile_835B_ramses = extractfield(datafile_835B_ramses,'name');
datafile_835B_ramses = datafile_835B_ramses(3);
datafile_835B = [datafile_835B '/' datafile_835B_ramses{1}];

[wavelength835B, spectrum835B] = ramses_get_835B([CalFolder '835B_deployedRadiometer/'],'835B', datafile_835B);

% 8396
datafile_8396 = [folderOfInterest folders(k).name '/irradiometer'];
datafile_8396_ramses =  dir(datafile_8396);
datafile_8396_ramses = extractfield(datafile_8396_ramses,'name');
datafile_8396_ramses = datafile_8396_ramses(3);
datafile_8396 = [datafile_8396 '/' datafile_8396_ramses{1}];

[wavelength8396, spectrum8396] = ramses_get_irradiometer([CalFolder '8396_deployedIrradiometer/'],'8396', datafile_8396);

% 835D
datafile_835D = [folderOfInterest folders(k).name '/waterRadiometer'];
datafile_835D_ramses =  dir(datafile_835D);
datafile_835D_ramses = extractfield(datafile_835D_ramses,'name');
datafile_835D_ramses = datafile_835D_ramses(3);
datafile_835D = [datafile_835D '/' datafile_835D_ramses{1}];

[wavelength835D, spectrum835D] = ramses_get_835B([CalFolder '835D_deployedRadiometer/'],'835D', datafile_835D);

figure(1)
subplot(337)
plot(wavelength835B,spectrum835B)
xlabel('Wavelength (nm)'); 
ylabel('Radiance (mW/(sr*m^2*nm)')
title('835B skyRadiometer')

subplot(338)
plot(wavelength835D,spectrum835D)
xlabel('Wavelength (nm)'); 
ylabel('Radiance (mW/(sr*m^2*nm)')
title('835D waterRadiometer')

subplot(339)
plot(wavelength8396,spectrum8396)
xlabel('Wavelength (nm)'); 
ylabel('Irradioance (mW/(m^2*nm)')
title('8396 irradiometer')
input('want to see the next data section?');
close(figure(2))
end


