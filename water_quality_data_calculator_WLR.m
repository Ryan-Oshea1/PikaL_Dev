%function [MAPE,RMSE,corr_coeff_ratio,not_nan_loc,predicted_beta,forty_deg_red_reflectance_numer] = water_quality_data_calculator_WLR(  WLR, Ed, WLR_lambda, Ed_lambda, water_quality_parameter_times, water_quality_parameter_values, wavelength_range_numer,wavelength_range_denom,order,set_num,set_num_name_rrs_array,noise_cutoff,noise_min,tss_range_end,verbose,spatial_pixel_start,spatial_pixel_end,angle_spacing,indices_training,indices_validation)
function [MAPE,RMSE,corr_coeff_ratio,not_nan_loc,predicted_beta,forty_deg_red_reflectance_numer] = water_quality_data_calculator_WLR(  WLR, Ed, WLR_lambda, Ed_lambda, water_quality_parameter_times_chl, water_quality_parameter_values_chl, wavelength_range_numer,wavelength_range_denom,order,noise_cutoff,noise_min,indices_training,indices_validation,srg_correction,wavelengths_rrs_opt,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,ls_wavelengths,srs_spectrum,all_images)

%clear forty_deg_red_reflectance_numer forty_deg_reflectance_denom   
%clear forty_deg_red_reflectance
% Ed_lambda = wavelength8396;
% WLR_lambda = wavelength835D;
% Ed = spectrum8396;
% WLR = spectrum835D;
%order = 2
%interpolate Ed to WLR wavelengths
Ed_wlr_lambda = interp1(Ed_lambda,Ed,WLR_lambda);

%calculate Rrs of WLR
Rrs_wlr = WLR./Ed_wlr_lambda;

%calculate the sample points based off of the input wavelength range
%wavelength_range_numer = 540:560;
%wavelength_range_denom = 480:500;

pixel_range_numer = ((WLR_lambda - min(wavelength_range_numer))>0).*((WLR_lambda - max(wavelength_range_numer))<0);
pixel_range_numer = find(pixel_range_numer);

wavelength_range_numer2 = 435:450;
pixel_range_numer2 = ((WLR_lambda - min(wavelength_range_numer2))>0).*((WLR_lambda - max(wavelength_range_numer2))<0);
pixel_range_numer2 = find(pixel_range_numer2);

wavelength_range_numer3 = 500:520;
pixel_range_numer3 = ((WLR_lambda - min(wavelength_range_numer3))>0).*((WLR_lambda - max(wavelength_range_numer3))<0);
pixel_range_numer3 = find(pixel_range_numer3);

pixel_range_denom = ((WLR_lambda - min(wavelength_range_denom))>0).*((WLR_lambda - max(wavelength_range_denom))<0);
pixel_range_denom = find(pixel_range_denom);

training_attempts = 0;
validation_attempts=0;
%use Rrs off WLR to calculate WQP



    for(i=indices_training)
     %   if( (hour(water_quality_parameter_times_chl(i)) == 12) || (hour(water_quality_parameter_times_chl(i)) == 13))
         if(1 == find_usable_sza(water_quality_parameter_times_chl(i),115,155))
            training_attempts = training_attempts+1;

     if(srg_correction == 1)

                i;
                trs = all_images(i).WLR_rrs_lee;%(Rrs_wlr(:,i));
                trs_wavelengths = WLR_lambda;
                
                finite_trs = isfinite(trs);
                finite_trs_loc = find(finite_trs);
                size_trs = size(finite_trs_loc);
                
                if(size_trs(1) <2)
                    value_numer = NaN;
                    value_denom = NaN;
                    disp('Setting to NaN due to trs not existing')
                    i
                else
%                     %calculate srs of the indices we need
%                     %ed_interped_to_ls = interp1(ed_wavelengths,ed,ls_wavelengths);
%                     srs = srs_spectrum(:,i);%ls(:,i)./ed(:,i);
%                   %  figure(66); plot(ls_wavelengths,srs); title(["number: " num2str(i)]);
%                     %define srs_wavelengths
%                     srs_wavelengths = ls_wavelengths;
% 
%                     %correct the average spatial data
% 
%                     [R_s, R_p, fresnel,theta_i_deg,theta_t_deg] = fresnelReflectanceCalculator(1,1.34, 90-45); %always at 45 degrees
% 
%                     %generate initial values for the min function
%                     pgxd = [1,1,1,1]; % should not affect IV
%                     %calculate IV
%                     [err, pgxd_0] = rrs_optimization_func(pgxd,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths)
%                     if(~isnan(err) && isreal(pgxd_0))
%                         %call the function here
%                         minimization_function_rrs_handle = @(pgxd)rrs_optimization_func(pgxd,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);
% 
%                         lower_bounds = [.003,.001,.0001,0]; % ZHONGPING, I set the delta bound to be positive
% 
%                         upper_bounds = [1,1,1,1];%[Inf, Inf, Inf, Inf] ZHONGPING, I set the bounds to be 1 for everything
% 
%                         %optimize!
%                         minimization_parameters = fmincon(minimization_function_rrs_handle,pgxd_0,[],[],[],[],lower_bounds,upper_bounds );
%                         [err, pgxd_0, rrs_mod, rrs_mea, bbp_650] = rrs_optimization_func(minimization_parameters,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);

                        value_numer = nanmean(nanmean(trs(pixel_range_numer)));
                        value_numer2 = nanmean(nanmean(trs(pixel_range_numer2)));
                        value_numer3 = nanmean(nanmean(trs(pixel_range_numer3)));

                        value_numer = max([value_numer,value_numer2,value_numer3]);

                        value_denom = nanmean(nanmean(trs(pixel_range_denom)));
%                     else
%                     value_numer = NaN;
%                     value_denom = NaN;
%                     disp('Setting to NaN due objective function error not existing fo IV')
%                     i
%                     end
                end
            else
            value_numer = nanmean(nanmean(Rrs_wlr(pixel_range_numer,i)));
            value_numer2 = nanmean(nanmean(Rrs_wlr(pixel_range_numer2,i)));
            value_numer3 = nanmean(nanmean(Rrs_wlr(pixel_range_numer3,i)));

            value_numer = max([value_numer,value_numer2,value_numer3]);

            value_denom = nanmean(nanmean(Rrs_wlr(pixel_range_denom,i)));
         end
            if( ((noise_min < value_numer) && (value_numer < noise_cutoff)) && ((noise_min < value_denom) && (value_denom < noise_cutoff)))
                forty_deg_red_reflectance_numer(i,:) = value_numer;
                forty_deg_reflectance_denom(i,:) = value_denom;
                not_nan_loc_finder_training(i) = 1;
             %   disp(['satisfy constraints value numer' num2str(value_numer) 'value denom' num2str(value_denom) ])

            else
                
                
               % disp(['Does not satisfy constraints value numer' num2str(value_numer) 'value denom' num2str(value_denom) ])
                forty_deg_red_reflectance_numer(i,:) = NaN;          
                forty_deg_reflectance_denom(i,:) = NaN;
                    not_nan_loc_finder_training(i) = 0;

            end
        end
    end


    forty_deg_red_reflectance_numer = forty_deg_red_reflectance_numer./forty_deg_reflectance_denom;

    not_nan_loc = (find(not_nan_loc_finder_training));%find(~isnan(forty_deg_red_reflectance_numer) );
     not_nan_loc = not_nan_loc';%.* ( ( (hour(water_quality_parameter_times_chl(not_nan_loc)) == 12) + (hour(water_quality_parameter_times_chl(not_nan_loc)) == 13)));
     not_nan_loc_loc = find(not_nan_loc);
     not_nan_loc = not_nan_loc(not_nan_loc_loc);
training_successes = length(not_nan_loc);

    y = log10((water_quality_parameter_values_chl(not_nan_loc))');  
    x = log10(forty_deg_red_reflectance_numer(not_nan_loc)');
    p = polyfit(x,y,order);
    y1= polyval(p,x);
  %    figure(7); scatter(x,y); hold on; scatter(x,y1); xlabel('Ratio'); ylabel('predicted and actual'); legend('Actual conc', 'Predicted Conc')


clear forty_deg_red_reflectance_numer forty_deg_reflectance_denom 
for(i=indices_validation)
  %if( (hour(water_quality_parameter_times_chl(i)) == 12) || (hour(water_quality_parameter_times_chl(i)) == 13))
   if(1 == find_usable_sza(water_quality_parameter_times_chl(i),115,155))
         validation_attempts=validation_attempts+1;

    if(srg_correction == 1)
    
                i;
                trs = all_images(i).WLR_rrs_lee;%(Rrs_wlr(:,i));
                trs_wavelengths = WLR_lambda;
                
                finite_trs = isfinite(trs);
                finite_trs_loc = find(finite_trs);
                size_trs = size(finite_trs_loc);
                
                if(size_trs(1) <2)
                    value_numer = NaN;
                    value_denom = NaN;
                    disp('Setting to NaN due to trs not existing')
                    i
                else
                    %calculate srs of the indices we need
                    %ed_interped_to_ls = interp1(ed_wavelengths,ed,ls_wavelengths);
%                     srs = srs_spectrum(:,i);%ls(:,i)./ed(:,i);
%                   %  figure(66); plot(ls_wavelengths,srs); title(["number: " num2str(i)]);
%                     %define srs_wavelengths
%                     srs_wavelengths = ls_wavelengths;
% 
%                     %correct the average spatial data
% 
%                     [R_s, R_p, fresnel,theta_i_deg,theta_t_deg] = fresnelReflectanceCalculator(1,1.34, 90-45); %always at 45 degrees
% 
%                     %generate initial values for the min function
%                     pgxd = [1,1,1,1]; % should not affect IV
%                     %calculate IV
%                     [err, pgxd_0] = rrs_optimization_func(pgxd,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths)
%                     if(~isnan(err) && isreal(pgxd_0))
%                         %call the function here
%                         minimization_function_rrs_handle = @(pgxd)rrs_optimization_func(pgxd,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);
% 
%                         lower_bounds = [.003,.001,.0001,0]; % ZHONGPING, I set the delta bound to be positive
% 
%                         upper_bounds = [1,1,1,1];%[Inf, Inf, Inf, Inf] ZHONGPING, I set the bounds to be 1 for everything
% 
%                         %optimize!
%                         minimization_parameters = fmincon(minimization_function_rrs_handle,pgxd_0,[],[],[],[],lower_bounds,upper_bounds );
%                         [err, pgxd_0, rrs_mod, rrs_mea, bbp_650] = rrs_optimization_func(minimization_parameters,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);
                        value_numer = nanmean(nanmean(trs(pixel_range_numer)));
                        value_numer2 = nanmean(nanmean(trs(pixel_range_numer2)));
                        value_numer3 = nanmean(nanmean(trs(pixel_range_numer3)));

                        value_numer = max([value_numer,value_numer2,value_numer3]);

                        value_denom = nanmean(nanmean(trs(pixel_range_denom)));
%                     else
%                                             value_numer = NaN;
%                     value_denom = NaN;
%                     disp('Setting to NaN due objective function error not existing fo IV')
%                     i
%                     end
                end
            else
            value_numer = nanmean(nanmean(Rrs_wlr(pixel_range_numer,i)));
            value_numer2 = nanmean(nanmean(Rrs_wlr(pixel_range_numer2,i)));
            value_numer3 = nanmean(nanmean(Rrs_wlr(pixel_range_numer3,i)));

            value_numer = max([value_numer,value_numer2,value_numer3]);
            value_denom = nanmean(nanmean(Rrs_wlr(pixel_range_denom,i)));
    end
        if( ((noise_min < value_numer) && (value_numer < noise_cutoff)) && ((noise_min < value_denom) && (value_denom < noise_cutoff)))
            forty_deg_red_reflectance_numer(i,:) = value_numer;
            forty_deg_reflectance_denom(i,:) = value_denom;
            not_nan_loc_finder(i) = 1;

        else
            
            % disp(['Does not satisfy constraints value numer' num2str(value_numer) 'value denom' num2str(value_denom) ])
            forty_deg_red_reflectance_numer(i,:) = NaN;          
            forty_deg_reflectance_denom(i,:) = NaN;
            not_nan_loc_finder(i) = 0;
        end
    
  end
end

forty_deg_red_reflectance_numer = forty_deg_red_reflectance_numer./forty_deg_reflectance_denom;

not_nan_loc = unique(find(not_nan_loc_finder)); %(~isnan(forty_deg_red_reflectance_numer(indices_validation)));    % +tss_range_end -1;
 not_nan_loc = not_nan_loc';%.* (   (hour(water_quality_parameter_times_chl(not_nan_loc)) == 12) + (hour(water_quality_parameter_times_chl(not_nan_loc)) == 13)   );
 not_nan_loc_loc = find(not_nan_loc);
 not_nan_loc = not_nan_loc(not_nan_loc_loc);
 validation_successes = length(not_nan_loc);

 %ratio
red_green_ratio_pred = log10(forty_deg_red_reflectance_numer(not_nan_loc)');  
predicted_beta = polyval(p,red_green_ratio_pred);

predicted_beta_holder(:) = polyval(p,red_green_ratio_pred);

  %  figure;scatter(red_green_ratio_pred,(water_quality_parameter_values_chl(not_nan_loc))); hold on; scatter(red_green_ratio_pred,predicted_beta); xlabel('Ratio'); ylabel('predicted and actual'); legend('Actual conc', 'Predicted Conc')
 %   figure(35);subplot(211);scatter(red_green_ratio_pred,log10(water_quality_parameter_values_chl(not_nan_loc)),10,'r'); hold on; scatter(red_green_ratio_pred,predicted_beta,'g'); xlabel('Ratio'); ylabel('predicted and actual'); legend('Actual conc', 'Predicted Conc')
  %  subplot(212);scatter(10.^red_green_ratio_pred,10.^log10(water_quality_parameter_values_chl(not_nan_loc)),'k'); hold on; scatter(10.^red_green_ratio_pred,10.^predicted_beta,'r'); xlabel('Ratio'); ylabel('predicted and actual'); legend('Actual conc', 'Predicted Conc')


% get correlation coefficient of Ratio and the NTU
% get statistics

ccr = corrcoef(10.^predicted_beta,(water_quality_parameter_values_chl(not_nan_loc)));
corr_coeff_ratio = ccr(1,2)
RMSE = sqrt(mean((water_quality_parameter_values_chl(not_nan_loc) - 10.^predicted_beta').^2));
MAPE = 100*mean(abs((water_quality_parameter_values_chl(not_nan_loc)- 10.^predicted_beta')./water_quality_parameter_values_chl(not_nan_loc)));%mean absolute percentage error


[min_MAPE_val,min_MAPE_loc] = min(MAPE);
%max value of beta
predicted_beta = predicted_beta_holder(min_MAPE_loc,:);
percent_error =(validation_successes+training_successes)/((training_attempts)+(validation_attempts)) *100;
forty_deg_red_reflectance_numer = percent_error;


end