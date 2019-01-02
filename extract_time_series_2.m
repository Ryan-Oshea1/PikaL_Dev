%this code allows us to save specific sections of the datasets

%folders that can be changed
folderOfInterest = '/media/flagg/15AFC52042FBF9AC/Pika_L/MVCO_datasets/';
folderOfInterest_restart = '/media/flagg/15AFC52042FBF9AC/Pika_L/MVCO_datasets/MVCO_restart/';

CalFolder = '/home/flagg/Desktop/downloadedDatasets/code/Calibration Files/'; 

%dont change these
HyperImage = '/hyperspectralImages/';
HyperImagePol = '/hyperspectralImagesPolarized/';
folders_og = dir(folderOfInterest)
folders_restart = dir(folderOfInterest_restart)


%set these parameters to download different datasets
select_polarized_images =  1;
select_gain_images = 0;
set_num_range = 8;
read_radiometer_data = 0;
folder_range_og = 1:(length(folders_og)-4);
folder_range_restart = 3:(length(folders_restart));

%turn to tables for easy concatenation
folders_og_table = struct2table(folders_og(folder_range_og))
folders_restart_table = struct2table(folders_restart(folder_range_restart))

combined_folders_table =  [folders_og_table;folders_restart_table]


% for counter = 1:length(field_names)
%    folders_og.(field_names{counter}) = folders_restart.(field_names{counter})
% end

folders = table2struct(combined_folders_table)%[dir(folderOfInterest)]; 

% remove the bad data from the stack
folders(1:2) = []; 
folders(2828) = []; 
folders(3022) = []; 
folders(3106) = []; 
folders(3655) = []; 
folders(3947) = []; 


folder_range_og = 1:(length(folders));
        
    



 %-2 is to account for the two other terminal files that are in the same
 %folder

%%

%
clc
for k = folder_range_og;
        if(find_usable_sza(water_quality_parameter_times_chl(k),105,165) == 1)

set_num = 9
set_num_range = 9
k
for(set_num = set_num_range)
    
if select_polarized_images == 0 && select_gain_images == 0
    disp('unpol un gain')
set_num_path = [folders(k).folder '/' folders(k).name HyperImage 'set_' num2str(set_num)];
    fget = dir([set_num_path '/AE/']);
    %if the first camera image is empty, set everything to empty?
    if isempty(fget);
    %ImAll(1:1200,1:1920,k,set_num) = nan;
%     spectrum835B(:,k) = nan;
%     spectrum8396(:,k) = nan;
%     spectrum835D(:,k) = nan;
    continue;
    end
    
    fget(1:2) = [];
    fget = struct2cell(fget);
    if sum(strcmp(fget(1,:),'autoExp1.png')) < 1;
    %ImAll(1:1200,1:1920,k,set_num) = nan;
%     spectrum835B(:,k) = nan;
%     spectrum835D(:,k) = nan;
%     spectrum8396(:,k) = nan;
    continue;
    end
      
%    iterate through all of the data sets
    hyperspectralImageAE = imread([set_num_path '/AE/autoExp1.png']);
    hyperspectralImageAE_OG = hyperspectralImageAE;
    hyperspectralImageAEHigh=bitshift(hyperspectralImageAE, 16-12);
    hyperspectralImageAElow = bitshift(hyperspectralImageAE, 16-28);
    hyperspectralImageAE =hyperspectralImageAEHigh+ hyperspectralImageAElow;
    hyperspectralImageAE = hyperspectralImageAE;
    %ImAll(:,:,k) = hyperspectralImageAE;
    if(set_num == 1)
       all_images(k).set_num_1 =  hyperspectralImageAE;
    end
    if(set_num == 2)
       all_images(k).set_num_2 =  hyperspectralImageAE;
    end
    if(set_num == 3)
       all_images(k).set_num_3 =  hyperspectralImageAE;
    end
    if(set_num == 4)
       all_images(k).set_num_4 =  hyperspectralImageAE;
    end
    if(set_num == 5)
       all_images(k).set_num_5 =  hyperspectralImageAE;
    end
    if(set_num == 6)
       all_images(k).set_num_6 =  hyperspectralImageAE;
    end
    if(set_num == 7)
       all_images(k).set_num_7 =  hyperspectralImageAE;
    end
    if(set_num == 8)
       all_images(k).set_num_8 =  hyperspectralImageAE;
    end
    if(set_num == 9)
        disp('writing image 9')
       all_images(k).set_num_9 =  hyperspectralImageAE;
    end
    if(set_num == 10)
       all_images(k).set_num_10 =  hyperspectralImageAE;
    end
    if(set_num == 11)
       all_images(k).set_num_11 =  hyperspectralImageAE;
    end
    if(set_num == 12)
       all_images(k).set_num_12 =  hyperspectralImageAE;
    end
    
    if size(double(str2num(fileread([set_num_path '/AE/AutoExposure.txt']))),1) > 0
        exposure_Im_All(k,set_num) = double(str2num(fileread([set_num_path '/AE/AutoExposure.txt'])));
    end
end


if select_polarized_images == 1 && select_gain_images == 0
    set_num_path = [folders(k).folder '/' folders(k).name HyperImagePol 'set_' num2str(set_num)];
    
        fget = dir([set_num_path '/AE/']);
    %if the first camera image is empty, set everything to empty?
    if isempty(fget);
    %ImAll(1:1200,1:1920,k,set_num) = nan;
%     spectrum835B(:,k) = nan;
%     spectrum8396(:,k) = nan;
%     spectrum835D(:,k) = nan;
    continue;
    end
    
    fget(1:2) = [];
    fget = struct2cell(fget);
    if sum(strcmp(fget(1,:),'autoExp1.png')) < 1;
    %ImAll(1:1200,1:1920,k,set_num) = nan;
%     spectrum835B(:,k) = nan;
%     spectrum835D(:,k) = nan;
%     spectrum8396(:,k) = nan;
    continue;
    end
    
    
 %   iterate through all of the data sets
    hyperspectralImageAE = imread([set_num_path '/AE/autoExp1.png']);
    hyperspectralImageAE_OG = hyperspectralImageAE;
    hyperspectralImageAEHigh=bitshift(hyperspectralImageAE, 16-12);
    hyperspectralImageAElow = bitshift(hyperspectralImageAE, 16-28);
    hyperspectralImageAE =hyperspectralImageAEHigh+ hyperspectralImageAElow;
    hyperspectralImageAE = hyperspectralImageAE;
    %ImAll(:,:,k) = hyperspectralImageAE;
    if(set_num == 1)
       all_images(k).set_num_1_pol =  hyperspectralImageAE;
    end
    if(set_num == 2)
       all_images(k).set_num_2_pol =  hyperspectralImageAE;
    end
    if(set_num == 3)
       all_images(k).set_num_3_pol =  hyperspectralImageAE;
    end
    if(set_num == 4)
       all_images(k).set_num_4_pol =  hyperspectralImageAE;
    end
    if(set_num == 5)
       all_images(k).set_num_5_pol =  hyperspectralImageAE;
    end
    if(set_num == 6)
       all_images(k).set_num_6_pol =  hyperspectralImageAE;
    end
    if(set_num == 7)
       all_images(k).set_num_7_pol =  hyperspectralImageAE;
    end
    if(set_num == 8)
       all_images(k).set_num_8_pol =  hyperspectralImageAE;
    end
    if(set_num == 9)
       all_images(k).set_num_9_pol =  hyperspectralImageAE;
    end
    if(set_num == 10)
       all_images(k).set_num_10_pol =  hyperspectralImageAE;
    end
    if(set_num == 11)
       all_images(k).set_num_11_pol =  hyperspectralImageAE;
    end
    if(set_num == 12)
       all_images(k).set_num_12_pol =  hyperspectralImageAE;
    end
    
    if size(double(str2num(fileread([set_num_path '/AE/AutoExposure.txt']))),1) > 0
        exposure_Im_All_pol(k,set_num) = double(str2num(fileread([set_num_path '/AE/AutoExposure.txt'])));
    end

end

    %exposure_Im_All(k,set_num) = double(str2num(fileread([folderOfInterest folders(k).name HyperImage 'set_' num2str(set_num) '/AE/AutoExposure.txt'])));




if select_polarized_images == 1 && select_gain_images == 1
    
    set_num_path = [folders(k).folder '/' folders(k).name HyperImagePol 'set_' num2str(set_num)];
        fget = dir([set_num_path '/AEG/']);
    %if the first camera image is empty, set everything to empty?
    if isempty(fget);
        empty = 1
    %ImAll(1:1200,1:1920,k,set_num) = nan;
    spectrum835B(:,k) = nan;
    spectrum8396(:,k) = nan;
    spectrum835D(:,k) = nan;
    continue;
    end
    
    fget(1:2) = [];
    fget = struct2cell(fget);
    if sum(strcmp(fget(1,:),'autoExpAndGain1.png')) < 1;
    %ImAll(1:1200,1:1920,k,set_num) = nan;
    spectrum835B(:,k) = nan;
    spectrum835D(:,k) = nan;
    spectrum8396(:,k) = nan;
    continue;
    end
    
    
 %   iterate through all of the data sets
    hyperspectralImageAE = imread([set_num_path '/AEG/autoExpAndGain1.png']);
    hyperspectralImageAE_OG = hyperspectralImageAE;
    hyperspectralImageAEHigh=bitshift(hyperspectralImageAE, 16-12);
    hyperspectralImageAElow = bitshift(hyperspectralImageAE, 16-28);
    hyperspectralImageAE =hyperspectralImageAEHigh+ hyperspectralImageAElow;
    hyperspectralImageAE = hyperspectralImageAE;
    %ImAll(:,:,k) = hyperspectralImageAE;
    if(set_num == 1)
       all_images(k).set_num_1_pol_gain =  hyperspectralImageAE;
    end
    if(set_num == 2)
       all_images(k).set_num_2_pol_gain =  hyperspectralImageAE;
    end
    if(set_num == 3)
       all_images(k).set_num_3_pol_gain =  hyperspectralImageAE;
    end
    if(set_num == 4)
       all_images(k).set_num_4_pol_gain =  hyperspectralImageAE;
    end
    if(set_num == 5)
       all_images(k).set_num_5_pol_gain =  hyperspectralImageAE;
    end
    if(set_num == 6)
       all_images(k).set_num_6_pol_gain =  hyperspectralImageAE;
    end
    if(set_num == 7)
       all_images(k).set_num_7_pol_gain =  hyperspectralImageAE;
    end
    if(set_num == 8)
       all_images(k).set_num_8_pol_gain =  hyperspectralImageAE;
    end
    if(set_num == 9)
       all_images(k).set_num_9_pol_gain =  hyperspectralImageAE;
    end
    if(set_num == 10)
       all_images(k).set_num_10_pol_gain =  hyperspectralImageAE;
    end
    if(set_num == 11)
       all_images(k).set_num_11_pol_gain =  hyperspectralImageAE;
    end
    if(set_num == 12)
       all_images(k).set_num_12_pol_gain =  hyperspectralImageAE;
    end
    
    if size(double(str2num(fileread([set_num_path '/AEG/AutoExposure.txt']))),1) > 0 && size(double(str2num(fileread([set_num_path '/AEG/AutoGain.txt']))),1) > 0
        exposure_Im_All_pol_gain(k,set_num) = double(str2num(fileread([set_num_path '/AEG/AutoExposure.txt'])));
        gain_Im_All_pol_gain(k,set_num) = double(str2num(fileread([set_num_path '/AEG/AutoGain.txt'])));
    end

end

    %exposure_Im_All(k,set_num) = double(str2num(fileread([folderOfInterest folders(k).name HyperImage 'set_' num2str(set_num) '/AE/AutoExposure.txt'])));


if select_polarized_images == 0 && select_gain_images == 1
    
    set_num_path = [folders(k).folder '/' folders(k).name HyperImage 'set_' num2str(set_num)];
        fget = dir([set_num_path '/AEG/']);
    %if the first camera image is empty, set everything to empty?
    if isempty(fget);
        empty = 1;
    %ImAll(1:1200,1:1920,k,set_num) = nan;
    spectrum835B(:,k) = nan;
    spectrum8396(:,k) = nan;
    spectrum835D(:,k) = nan;
    continue;
    end
    
    fget(1:2) = [];
    fget = struct2cell(fget);
    if sum(strcmp(fget(1,:),'autoExpAndGain1.png')) < 1;
    %ImAll(1:1200,1:1920,k,set_num) = nan;
    spectrum835B(:,k) = nan;
    spectrum835D(:,k) = nan;
    spectrum8396(:,k) = nan;
    continue;
    end
    
    
 %   iterate through all of the data sets
    hyperspectralImageAE = imread([set_num_path '/AEG/autoExpAndGain1.png']);
    hyperspectralImageAE_OG = hyperspectralImageAE;
    hyperspectralImageAEHigh=bitshift(hyperspectralImageAE, 16-12);
    hyperspectralImageAElow = bitshift(hyperspectralImageAE, 16-28);
    hyperspectralImageAE =hyperspectralImageAEHigh+ hyperspectralImageAElow;
    hyperspectralImageAE = hyperspectralImageAE;
    %ImAll(:,:,k) = hyperspectralImageAE;
    if(set_num == 1)
       all_images(k).set_num_1_gain =  hyperspectralImageAE;
    end
    if(set_num == 2)
       all_images(k).set_num_2_gain =  hyperspectralImageAE;
    end
    if(set_num == 3)
       all_images(k).set_num_3_gain =  hyperspectralImageAE;
    end
    if(set_num == 4)
       all_images(k).set_num_4_gain =  hyperspectralImageAE;
    end
    if(set_num == 5)
       all_images(k).set_num_5_gain =  hyperspectralImageAE;
    end
    if(set_num == 6)
       all_images(k).set_num_6_gain =  hyperspectralImageAE;
    end
    if(set_num == 7)
       all_images(k).set_num_7_gain =  hyperspectralImageAE;
    end
    if(set_num == 8)
       all_images(k).set_num_8_gain =  hyperspectralImageAE;
    end
    if(set_num == 9)
       all_images(k).set_num_9_gain =  hyperspectralImageAE;
    end
    if(set_num == 10)
       all_images(k).set_num_10_gain =  hyperspectralImageAE;
    end
    if(set_num == 11)
       all_images(k).set_num_11_gain =  hyperspectralImageAE;
    end
    if(set_num == 12)
       all_images(k).set_num_12_gain =  hyperspectralImageAE;
    end
    
    if size(double(str2num(fileread([set_num_path '/AEG/AutoExposure.txt']))),1) > 0 && size(double(str2num(fileread([set_num_path '/AEG/AutoGain.txt']))),1) > 0
        exposure_Im_All(k,set_num) = double(str2num(fileread([set_num_path '/AEG/AutoExposure.txt'])));
        gain_Im_All(k,set_num) = double(str2num(fileread([set_num_path '/AEG/AutoGain.txt'])));
    end

end


% 
    if read_radiometer_data == 1
        datafile_835B = [folderOfInterest folders(k).name '/skyRadiometer'];
        datafile_835B_ramses =  dir(datafile_835B);
        if length(datafile_835B_ramses) < 3;
            spectrum835B(1:255,k) = nan;
        else
            datafile_835B_ramses = extractfield(datafile_835B_ramses,'name');
            datafile_835B_ramses = datafile_835B_ramses(3);
            datafile_835B = [datafile_835B '/' datafile_835B_ramses{1}];
            [wavelength835B, spectrum835B(1:255,k)] = ramses_get_835B([CalFolder '835B_deployedRadiometer/'],'835B', datafile_835B);
            %if there was no data in the first file, try the second file
            if isnan(wavelength835B)
                disp('Reading second radiometer file, since nan wavelength 835B');
                datafile_835B = [folderOfInterest folders(k).name '/skyRadiometer'];
                datafile_835B_ramses =  dir(datafile_835B);
                if length(datafile_835B_ramses) < 3;
                    spectrum835B(1:255,k) = nan;
                else
                    datafile_835B_ramses = extractfield(datafile_835B_ramses,'name');
                    datafile_835B_ramses = datafile_835B_ramses(4);
                    datafile_835B = [datafile_835B '/' datafile_835B_ramses{1}];
                    [wavelength835B, spectrum835B(1:255,k)] = ramses_get_835B([CalFolder '835B_deployedRadiometer/'],'835B', datafile_835B);
                        if isnan(wavelength835B)
                             disp('Second Radiometer Read Failed 835B');
                        end
                end
            end
        end
        
        % 8396
        datafile_8396 = [folderOfInterest folders(k).name '/irradiometer'];
        datafile_8396_ramses =  dir(datafile_8396);
        if length(datafile_8396_ramses) < 3;
            spectrum8396(1:255,k) = nan;
        else
            datafile_8396_ramses = extractfield(datafile_8396_ramses,'name');
            datafile_8396_ramses = datafile_8396_ramses(3);
            datafile_8396 = [datafile_8396 '/' datafile_8396_ramses{1}];
            [wavelength8396, spectrum8396(1:255,k)] = ramses_get_irradiometer([CalFolder '8396_deployedIrradiometer/'],'8396', datafile_8396);
            wavelength8396;
            if isnan(wavelength8396)
                disp('Reading second radiometer file, since nan wavelength 8396');
                datafile_8396 = [folderOfInterest folders(k).name '/irradiometer'];
                datafile_8396_ramses =  dir(datafile_8396);
                if length(datafile_8396_ramses) < 3;
                    spectrum8396(1:255,k) = nan;
                else
                    datafile_8396_ramses = extractfield(datafile_8396_ramses,'name');
                    datafile_8396_ramses = datafile_8396_ramses(4);
                    datafile_8396 = [datafile_8396 '/' datafile_8396_ramses{1}];
                    [wavelength8396, spectrum8396(1:255,k)] = ramses_get_irradiometer([CalFolder '8396_deployedIrradiometer/'],'8396', datafile_8396);
                        if isnan(wavelength8396)
                             disp('Second Radiometer Read Failed 8396');
                        end
                end
            end
        end
        % 835D
        datafile_835D = [folderOfInterest folders(k).name '/waterRadiometer'];
        datafile_835D_ramses =  dir(datafile_835D);
        if length(datafile_835D_ramses) < 3;
            spectrum835D(1:255,k) = nan;
        else
            datafile_835D_ramses = extractfield(datafile_835D_ramses,'name');
            datafile_835D_ramses = datafile_835D_ramses(3);
            datafile_835D = [datafile_835D '/' datafile_835D_ramses{1}];
            [wavelength835D, spectrum835D(1:255,k)] = ramses_get_835B([CalFolder '835D_deployedRadiometer/'],'835D', datafile_835D);
            if isnan(wavelength835D)
                disp('Reading second radiometer file, since nan wavelength 835D');
                datafile_835D = [folderOfInterest folders(k).name '/waterRadiometer'];
                datafile_835D_ramses =  dir(datafile_835D);
                if length(datafile_835D_ramses) < 3;
                    spectrum835D(1:255,k) = nan;
                else
                    datafile_835D_ramses = extractfield(datafile_835D_ramses,'name');
                    datafile_835D_ramses = datafile_835D_ramses(4);
                    datafile_835D = [datafile_835D '/' datafile_835D_ramses{1}];
                    [wavelength835D, spectrum835D(1:255,k)] = ramses_get_835B([CalFolder '835D_deployedRadiometer/'],'835D', datafile_835D);
                        if isnan(wavelength835D)
                             disp('Second Radiometer Read Failed 835D');
                        end
                end
            end
  
        end
    end
end
        end
end
%%
%clear the unusable data

%%
clear datafile_835B datafile_835D datafile_8396 datafile_835B_ramses datafile_835D_ramses datafile_8396_ramses
clear fget folder_range HyperImage HyperImagePol hyperspectralImageAE hyperspectralImageAE_OG hyperspectralImageAEHigh hyperspectralImageAElow
clear read_radiometer_data select_gain_images select_polarized_images set_num_path set_num_range
clear k

%%
%clear the unusable data
all