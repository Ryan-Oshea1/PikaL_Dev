function [rrs_lee_image] = software_correct_all_data2(  all_images,camera_wavelengths, water_quality_parameter_times, water_quality_parameter_values, wavelength_range_numer,wavelength_range_denom,order,set_num,set_num_name_rrs_array,noise_cutoff,noise_min,tss_range_end,verbose,spatial_pixel_start,spatial_pixel_end,angle_spacing,indices_training,indices_validation,srg_correction,wavelengths_rrs_opt,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,ls_wavelengths,srs_spectrum,spatial_pixel_angles,i)
%for(i = 1:2)

if(1 == find_usable_sza(water_quality_parameter_times(i),115,155))
if(set_num == 9)
    camera_row = 3;
    camera_columns = 1:150;
end
if(set_num == 7)
    camera_row = 1;
    camera_columns = 1:600;
end

for ( pixel = 1:size(all_images(i).(set_num_name_rrs_array(set_num)),2)) % for all of the spatial pixels
pixel
  %  close all
       % tic
      if(1 == find_usable_sza(water_quality_parameter_times(i),115,155))
            %    if((hour(water_quality_parameter_times(i)) == 12) ||   (hour(water_quality_parameter_times(i)) == 13))
           % training_attempts = training_attempts+1;

        if(size(all_images(i).(set_num_name_rrs_array(set_num)),1) >100 &&  size(all_images(i).(set_num_name_rrs_array(set_num)),2) >100)
              if(srg_correction == 1)

                i
                trs = nanmean((all_images(i).(set_num_name_rrs_array(set_num))(:,pixel)),2);
                trs_wavelengths = camera_wavelengths(camera_row,camera_columns);
                
                finite_trs = isfinite(trs);
                finite_trs_loc = find(finite_trs);
                size_trs = size(finite_trs_loc);
                
                if(size_trs <2)
                    rrs_lee_holder(:,pixel) = NaN(600,1);
                    disp('Setting to NaN due to trs not existing')
                    i
                else
                    %calculate srs of the indices we need
                    %ed_interped_to_ls = interp1(ed_wavelengths,ed,ls_wavelengths);
                    srs = srs_spectrum(:,i);%ls(:,i)./ed(:,i);
                  %  figure(66); plot(ls_wavelengths,srs); title(["number: " num2str(i)]);
                    %define srs_wavelengths
                    srs_wavelengths = ls_wavelengths;

                    %correct the average spatial data

                    [R_s, R_p, fresnel,theta_i_deg,theta_t_deg] = fresnelReflectanceCalculator(1,1.34, 90-spatial_pixel_angles(pixel));

                    %generate initial values for the min function
                    pgxd = [1,1,1,1]; % should not affect IV
                    %calculate IV
                    [err, pgxd_0] = rrs_optimization_func(pgxd,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);
                    Initial_values_pgxd0 = pgxd_0;
                    if(isreal(pgxd_0))
                                        %call the function here
                                        minimization_function_rrs_handle = @(pgxd)rrs_optimization_func(pgxd,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);

                                        lower_bounds = [.003,.001,.0001,0]; % ZHONGPING, I set the delta bound to be positive

                                        upper_bounds = [1,1,1,1];%[Inf, Inf, Inf, Inf] ZHONGPING, I set the bounds to be 1 for everything

                                        %optimize!
                                        minimization_parameters = fmincon(minimization_function_rrs_handle,pgxd_0,[],[],[],[],lower_bounds,upper_bounds );
                                        [err, pgxd_0, rrs_mod, rrs_mea, bbp_650] = rrs_optimization_func(minimization_parameters,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);
                                        disp('size rrs_mea')

                                        size(rrs_mea)
                                        rrs_lee_holder(:,pixel) = rrs_mea;
                    else
                        rrs_lee_holder(:,pixel) = nan(1,length(camera_columns))  ; 
                        disp('unreal pgxd_0')
                    end
                end

            end
            

        end
      end  
     % toc
    
end
%if(isreal(pgxd_0))
%all_images(i).set_num_7_rrs_lee = rrs_lee_holder
rrs_lee_image = rrs_lee_holder;
clear rrs_lee_holder
%else
    %    rrs_lee_image = []
%end
else
    rrs_lee_image = []
end
%end
end