%interpolate spectrums
clc
ed_interpolated_ls = interp1(wavelength8396,spectrum8396,wavelength835B);

srs_spectrum = spectrum835B./ed_interpolated_ls;
%%
                figure(66);subplot(212); plot(wavelength835B,srs_spectrum(:,921));
                % srs = srs_spectrum(:,921);%ls(:,i)./ed(:,i);
                %figure(66); plot(wavelength835B,srs);
                
                %%
                %save the srg data
                index = 4193
                practice_image = all_images(index).(set_num_name_rrs_array(set_num));
                srs_practice = srs_spectrum(:,index)
                srs_wavelengths_practice = wavelength835B
                trs_wavelengths_practice = camera_wavelengths(1,:)
                save('practice_srg_correction', 'practice_image','srs_practice','srs_wavelengths_practice','trs_wavelengths_practice')