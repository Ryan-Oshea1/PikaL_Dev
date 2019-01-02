function [angle_center,MAPE,RMSE,corr_coeff_ratio,not_nan_loc,predicted_beta,forty_deg_red_reflectance_numer] = water_quality_data_calculator(  all_images,camera_wavelengths, water_quality_parameter_times, water_quality_parameter_values, wavelength_range_numer,wavelength_range_denom,order,set_num,set_num_name_rrs_array,noise_cutoff,noise_min,tss_range_end,verbose,spatial_pixel_start,spatial_pixel_end,angle_spacing,indices_training,indices_validation,srg_correction,wavelengths_rrs_opt,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,ls_wavelengths,srs_spectrum)
index = 1;

%for ( angle_increment = spatial_pixel_start:angle_spacing:spatial_pixel_end)
%for ( angle_increment =linspace(0,13.5-angle_spacing,(13.5-angle_spacing)/angle_spacing))
for ( angle_increment =1:225)%linspace(0,13.5-angle_spacing,(13.5-angle_spacing)/angle_spacing))
angle_increment
%     clear forty_deg_red_reflectance_numer forty_deg_reflectance_denom 
% 
%     angle_start = 39.75+angle_increment;
%     angle_end = 39.75+ angle_increment+angle_spacing;
%     spatial_pixel_start = round((900/13.5)* (53.25-angle_end));
%     spatial_pixel_end = round((900/13.5)* (53.25-angle_start));
%     spatial_pixel_range = spatial_pixel_start:spatial_pixel_end;

    clear forty_deg_red_reflectance_numer forty_deg_reflectance_denom 

    angle_start = 39.75+angle_increment;
    angle_end = 39.75+ angle_increment+angle_spacing;
    spatial_pixel_start = ceil((224/13.5)* (53.25-angle_end)+.0001);
    spatial_pixel_end = ceil((224/13.5)* (53.25-angle_start)+.0001);
    %spatial_pixel_range = spatial_pixel_start:spatial_pixel_end;
    spatial_pixel_range = angle_increment;
%         angle_start = 36+angle_increment;
%     angle_end = 36.5 + angle_increment;
%     spatial_pixel_start = round((900/17.6)* (54-angle_end));
%     spatial_pixel_end = round((900/17.6)* (54-angle_start));
%     spatial_pixel_range = spatial_pixel_start:spatial_pixel_end;
    camera_bin_select = 3;
pixel_range_numer = ((camera_wavelengths(camera_bin_select,:) - min(wavelength_range_numer))>0).*((camera_wavelengths(camera_bin_select,:) - max(wavelength_range_numer))<0);
pixel_range_numer = find(pixel_range_numer);

wavelength_range_numer2 = 435:450;
pixel_range_numer2 = ((camera_wavelengths(camera_bin_select,:) - min(wavelength_range_numer2))>0).*((camera_wavelengths(camera_bin_select,:) - max(wavelength_range_numer2))<0);
pixel_range_numer2 = find(pixel_range_numer2);

wavelength_range_numer3 = 500:520;
pixel_range_numer3 = ((camera_wavelengths(camera_bin_select,:) - min(wavelength_range_numer3))>0).*((camera_wavelengths(camera_bin_select,:) - max(wavelength_range_numer3))<0);
pixel_range_numer3 = find(pixel_range_numer3);

pixel_range_denom = ((camera_wavelengths(camera_bin_select,:) - min(wavelength_range_denom))>0).*((camera_wavelengths(camera_bin_select,:) - max(wavelength_range_denom))<0);
pixel_range_denom = find(pixel_range_denom);

training_attempts = 0;
validation_attempts=0;
    clear forty_deg_red_reflectance
  %  close all
    for(i=indices_training)
       % tic
      if(1 == find_usable_sza(water_quality_parameter_times(i),115,155))
            %    if((hour(water_quality_parameter_times(i)) == 12) ||   (hour(water_quality_parameter_times(i)) == 13))
            training_attempts = training_attempts+1;

        if(size(all_images(i).(set_num_name_rrs_array(set_num)),1) >100 &&  size(all_images(i).(set_num_name_rrs_array(set_num)),2) >100)
              if(srg_correction == 1)

                i
                trs = nanmean((all_images(i).(set_num_name_rrs_array(set_num))(:,spatial_pixel_range)),2);
                trs_wavelengths = camera_wavelengths(camera_bin_select,:);
                
                finite_trs = isfinite(trs);
                finite_trs_loc = find(finite_trs);
                size_trs = size(finite_trs_loc);
                
                if(size_trs <2)
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
%                     [R_s, R_p, fresnel,theta_i_deg,theta_t_deg] = fresnelReflectanceCalculator(1,1.34, 90-(angle_start+angle_end)/2);
% 
%                     %generate initial values for the min function
%                     pgxd = [1,1,1,1]; % should not affect IV
%                     %calculate IV
%                     [err, pgxd_0] = rrs_optimization_func(pgxd,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);
%                     Initial_values_pgxd0 = pgxd_0;
%                     if(isreal(pgxd_0))
%                                         %call the function here
%                                         minimization_function_rrs_handle = @(pgxd)rrs_optimization_func(pgxd,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);
% 
%                                         lower_bounds = [.003,.001,.0001,0]; % ZHONGPING, I set the delta bound to be positive
% 
%                                         upper_bounds = [1,1,1,1];%[Inf, Inf, Inf, Inf] ZHONGPING, I set the bounds to be 1 for everything
% 
%                                         %optimize!
%                                         minimization_parameters = fmincon(minimization_function_rrs_handle,pgxd_0,[],[],[],[],lower_bounds,upper_bounds );
%                                         [err, pgxd_0, rrs_mod, rrs_mea, bbp_650] = rrs_optimization_func(minimization_parameters,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);
                                           rrs_mea=all_images(i).set_num_9_rrs_lee_2_t(:,spatial_pixel_range);
                                           size(rrs_mea)
                                        value_numer = nanmean(nanmean(rrs_mea(pixel_range_numer)));
                                        value_numer2 = nanmean(nanmean(rrs_mea(pixel_range_numer2)));
                                        value_numer3 = nanmean(nanmean(rrs_mea(pixel_range_numer3)));

                                        value_numer = max([value_numer,value_numer2,value_numer3]);

                                        value_denom = nanmean(nanmean(rrs_mea(pixel_range_denom)));
                 %   else
                 %       value_numer = NaN;
                 %       value_denom = NaN
                  %  end
                end
            else
                value_numer = nanmean(nanmean(all_images(i).(set_num_name_rrs_array(set_num))(pixel_range_numer,spatial_pixel_range)));
                value_numer2 = nanmean(nanmean(all_images(i).(set_num_name_rrs_array(set_num))(pixel_range_numer2,spatial_pixel_range)));
                value_numer3 = nanmean(nanmean(all_images(i).(set_num_name_rrs_array(set_num))(pixel_range_numer3,spatial_pixel_range)));

                value_numer = max([value_numer,value_numer2,value_numer3]);

                value_denom = nanmean(nanmean(all_images(i).(set_num_name_rrs_array(set_num))(pixel_range_denom,spatial_pixel_range)));
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
     % toc
     end


    forty_deg_red_reflectance_numer = forty_deg_red_reflectance_numer./forty_deg_reflectance_denom;

    not_nan_loc = (find(not_nan_loc_finder_training));%find(~isnan(forty_deg_red_reflectance_numer) );
     not_nan_loc = not_nan_loc';%.* (   (hour(water_quality_parameter_times(not_nan_loc)) == 12) +   (hour(water_quality_parameter_times(not_nan_loc)) == 13)  );
     not_nan_loc_loc = find(not_nan_loc);
     not_nan_loc = not_nan_loc(not_nan_loc_loc);
training_successes = length(not_nan_loc);
    y = log10((water_quality_parameter_values(not_nan_loc))');  
    x = log10(forty_deg_red_reflectance_numer(not_nan_loc)');
    p = polyfit(x,y,order);
    y1= polyval(p,x);
    p;
    if(verbose == 1)
     figure(30);  subplot(211); scatter(x,y,'k'); hold on; scatter(x,y1,'r'); subplot(212); scatter(10.^x,10.^y,'k'); hold on; scatter(10.^x,10.^y1,'r')
    end

% myfittype = fittype('a * (x)^2.7848 +b',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'a','b'});
% myfit = fit(x',y',myfittype);
% if(verbose == 1)
% figure(30);hold on; plot(myfit,x,y); xlabel('ratio'); ylabel('concentration');
% end
CORRELATION_COEFFICIENT = corrcoef(x,y);

if(verbose == 1)
    figure(31); 
    plot((water_quality_parameter_values(not_nan_loc))/max((water_quality_parameter_values(not_nan_loc))));  xlabel('Data timeseries'); ylabel('normalized concentration');
end


clear forty_deg_red_reflectance_numer forty_deg_reflectance_denom 
for(i=indices_validation)
 % if((hour(water_quality_parameter_times(i)) == 12) ||   (hour(water_quality_parameter_times(i)) == 13))
      if(1 == find_usable_sza(water_quality_parameter_times(i),115,155))
        validation_attempts=validation_attempts+1;
    if(size(all_images(i).(set_num_name_rrs_array(set_num)),1) >100 &&  size(all_images(i).(set_num_name_rrs_array(set_num)),2) >100)
        if(srg_correction == 1)

                i
                trs = nanmean((all_images(i).(set_num_name_rrs_array(set_num))(:,spatial_pixel_range)),2);
                trs_wavelengths = camera_wavelengths(camera_bin_select,:);
                
                finite_trs = isfinite(trs);
                finite_trs_loc = find(finite_trs);
                size_trs = size(finite_trs_loc);
                
                if(size_trs <2)
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
%                     [R_s, R_p, fresnel,theta_i_deg,theta_t_deg] = fresnelReflectanceCalculator(1,1.34, 90-(angle_start+angle_end)/2);
% 
%                     %generate initial values for the min function
%                     pgxd = [1,1,1,1]; % should not affect IV
%                     %calculate IV
%                     [err, pgxd_0] = rrs_optimization_func(pgxd,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);
%                        Initial_values_pgxd0 = pgxd_0;
%                         if(isreal(pgxd_0))
% 
%                     %call the function here
%                     minimization_function_rrs_handle = @(pgxd)rrs_optimization_func(pgxd,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);
% 
%                     lower_bounds = [.003,.001,.0001,0]; % ZHONGPING, I set the delta bound to be positive
% 
%                     upper_bounds = [1,1,1,1];%[Inf, Inf, Inf, Inf] ZHONGPING, I set the bounds to be 1 for everything
% 
%                     %optimize!
%                     minimization_parameters = fmincon(minimization_function_rrs_handle,pgxd_0,[],[],[],[],lower_bounds,upper_bounds );
    %                     [err, pgxd_0, rrs_mod, rrs_mea, bbp_650] = rrs_optimization_func(minimization_parameters,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);
                    rrs_mea=all_images(i).set_num_9_rrs_lee_2_t(:,spatial_pixel_range);
                    size(rrs_mea)
                    value_numer = nanmean(nanmean(rrs_mea(pixel_range_numer)));
                    value_numer2 = nanmean(nanmean(rrs_mea(pixel_range_numer2)));
                    value_numer3 = nanmean(nanmean(rrs_mea(pixel_range_numer3)));
                    
                    value_numer = max([value_numer,value_numer2,value_numer3]);
                    value_denom = nanmean(nanmean(rrs_mea(pixel_range_denom)));
                    
                  %      else
                %        value_numer = NaN;
               %         value_denom = NaN
                %    end
                end
        else
            value_numer = nanmean(nanmean(all_images(i).(set_num_name_rrs_array(set_num))(pixel_range_numer,spatial_pixel_range)));
            value_numer2 = nanmean(nanmean(all_images(i).(set_num_name_rrs_array(set_num))(pixel_range_numer2,spatial_pixel_range)));
            value_numer3 = nanmean(nanmean(all_images(i).(set_num_name_rrs_array(set_num))(pixel_range_numer3,spatial_pixel_range)));

            value_numer = max([value_numer,value_numer2,value_numer3]);
            value_denom = nanmean(nanmean(all_images(i).(set_num_name_rrs_array(set_num))(pixel_range_denom,spatial_pixel_range)));
        end
        
        if( (((noise_min < value_numer) && (value_numer < noise_cutoff)) && ((noise_min < value_denom) && (value_denom < noise_cutoff))))
            forty_deg_red_reflectance_numer(i,:) = value_numer;
            forty_deg_reflectance_denom(i,:) = value_denom;
            not_nan_loc_finder(i) = 1;
%             if(i<1000)
%                 i
%               %  not_nan_loc_finder(i)
%             end
        else
            
            % disp(['Does not satisfy constraints value numer' num2str(value_numer) 'value denom' num2str(value_denom) ])
            forty_deg_red_reflectance_numer(i,:) = NaN;          
            forty_deg_reflectance_denom(i,:) = NaN;
            not_nan_loc_finder(i) = 0;
        end
    end
  end
end

%not_nan_loc_finder
forty_deg_red_reflectance_numer = forty_deg_red_reflectance_numer./forty_deg_reflectance_denom;



not_nan_loc = unique(find(not_nan_loc_finder)); %(~isnan(forty_deg_red_reflectance_numer(indices_validation)));    % +tss_range_end -1;
%find(not_nan_loc_finder)
%not_nan_loc = find(~isnan(forty_deg_red_reflectance) );
 not_nan_loc = not_nan_loc';%.* (   (hour(water_quality_parameter_times(not_nan_loc)) == 12) +   (hour(water_quality_parameter_times(not_nan_loc)) == 13)   );
 not_nan_loc_loc = find(not_nan_loc);
 not_nan_loc = not_nan_loc(not_nan_loc_loc);
 validation_successes = length(not_nan_loc);

%not_nan_loc_holder (index,:) = not_nan_loc;
%p=[[-.1054,-3.0777,1.7627,-1.9369,.2389]]
red_green_ratio_pred = log10(forty_deg_red_reflectance_numer(not_nan_loc)');  
predicted_beta = polyval(p,red_green_ratio_pred);

predicted_beta_holder(index).value = polyval(p,red_green_ratio_pred);
not_nan_loc_holder(index).value = not_nan_loc;

if(verbose == 1)
    figure(32);subplot(211);scatter(red_green_ratio_pred,log10(water_quality_parameter_values(not_nan_loc)),10,[(1-angle_increment/13.5) 0 angle_increment/13.5]); hold on; scatter(red_green_ratio_pred,predicted_beta,'g'); xlabel('Ratio'); ylabel('predicted and actual'); legend('Actual conc', 'Predicted Conc')
    subplot(212);scatter(10.^red_green_ratio_pred,10.^log10(water_quality_parameter_values(not_nan_loc)),'k'); hold on; scatter(10.^red_green_ratio_pred,10.^predicted_beta,'r'); xlabel('Ratio'); ylabel('predicted and actual'); legend('Actual conc', 'Predicted Conc')
end
%waitforbuttonpress
if(verbose == 1)
    oc4_graph = imread('oc4_plot.png');
    figure(37);
    oc4_graph = flipdim(oc4_graph,2);

    % image('CData',oc4_graph,'XData',[logspace(.01, 25)],'YData',[logspace(.001, 100)])
     xlim([.1 25]);
     ylim([ .001 100])
     hold on
     scatter(10.^red_green_ratio_pred,10.^log10(water_quality_parameter_values(not_nan_loc)),'*g'); hold on; scatter(10.^red_green_ratio_pred,10.^predicted_beta,'b'); xlabel('Ratio'); ylabel('predicted and actual'); legend('Actual conc', 'Predicted Conc')
     set(gca,'xscale','log'); set(gca,'yscale','log');
waitforbuttonpress
end
% get correlation coefficient of Ratio and the NTU
% get statistics
%[ccr, p_values] = corrcoef(10.^forty_deg_red_reflectance_numer(not_nan_loc),(water_quality_parameter_values(not_nan_loc)));
[ccr, p_values] = corrcoef(10.^predicted_beta,(water_quality_parameter_values(not_nan_loc)));

corr_coeff_ratio(index) = ccr(1,2);
p_value = p_values(1,2);
RMSE(index) = sqrt(mean((water_quality_parameter_values(not_nan_loc) - 10.^predicted_beta').^2));
MAPE(index) = 100*mean(abs((water_quality_parameter_values(not_nan_loc)- 10.^predicted_beta'))./(water_quality_parameter_values(not_nan_loc)));%mean absolute percentage error
%angle_center(index) = 39.75+angle_increment;
angle_center(index) = 53.25 - 13.5*(angle_increment-1)/224;

percent_error(index) =(validation_successes+training_successes)/((training_attempts)+(validation_attempts)) *100;
    index=index+1;
end

if(verbose == 1)
    figure(33);scatter((water_quality_parameter_values(not_nan_loc)),10.^predicted_beta); xlabel('Fluorometer Measured'); ylabel('Optically Estimated'); set(gca,'xscale','log'); set(gca,'yscale','log');
    hold on
    plot(.5:15,.5:15)
    xlim([.4 15])
    ylim([.4 15])
end

[min_MAPE_val,min_MAPE_loc] = min(MAPE);
%max value of beta
predicted_beta = predicted_beta_holder;
not_nan_loc = not_nan_loc_holder;
%not_nan_loc = not_nan_loc_holder(min_MAPE_loc,:);

forty_deg_red_reflectance_numer = percent_error
corr_coeff_ratio
end