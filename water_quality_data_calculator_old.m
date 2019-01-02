function [angle_center,MAPE,RMSE,corr_coeff_ratio,not_nan_loc,predicted_beta,forty_deg_red_reflectance_numer] = water_quality_data_calculator(  all_images, water_quality_parameter_times, water_quality_parameter_values, wavelength_range_numer,wavelength_range_denom,order,set_num,set_num_name_rrs_array,noise_cutoff,noise_min,tss_range_end,verbose,spatial_pixel_start,spatial_pixel_end,angle_spacing,indices_training,indices_validation)
index = 1;

for ( angle_increment = spatial_pixel_start:angle_spacing:spatial_pixel_end)
    clear forty_deg_red_reflectance_numer forty_deg_reflectance_denom 

    angle_start = 36+angle_increment;
    angle_end = 36.5 + angle_increment;
    spatial_pixel_start = round((900/17.6)* (54-angle_end));
    spatial_pixel_end = round((900/17.6)* (54-angle_start));
    spatial_pixel_range = spatial_pixel_start:spatial_pixel_end;

    
    
    
    clear forty_deg_red_reflectance
    close all
    for(i=indices_training)
        if(size(all_images(i).(set_num_name_rrs_array(set_num)),1) >400 &&  size(all_images(i).(set_num_name_rrs_array(set_num)),2) >400)
            value_numer = nanmean(nanmean(all_images(i).(set_num_name_rrs_array(set_num))(wavelength_range_numer,spatial_pixel_range)));
            value_denom = nanmean(nanmean(all_images(i).(set_num_name_rrs_array(set_num))(wavelength_range_denom,spatial_pixel_range)));

            if( ((noise_min < value_numer) && (value_numer < noise_cutoff)) && ((noise_min < value_denom) && (value_denom < noise_cutoff)))
                forty_deg_red_reflectance_numer(i,:) = value_numer;
                forty_deg_reflectance_denom(i,:) = value_denom;
                            not_nan_loc_finder_training(i) = 1;

            else
                %disp('Does not satisfy constraints')
                forty_deg_red_reflectance_numer(i,:) = NaN;          
                forty_deg_reflectance_denom(i,:) = NaN;
                    not_nan_loc_finder_training(i) = 0;

            end
        end
    end


    forty_deg_red_reflectance_numer = forty_deg_red_reflectance_numer./forty_deg_reflectance_denom;

    not_nan_loc = (find(not_nan_loc_finder_training));%find(~isnan(forty_deg_red_reflectance_numer) );
     not_nan_loc = not_nan_loc'.* (  (hour(water_quality_parameter_times(not_nan_loc)) == 12));
     not_nan_loc_loc = find(not_nan_loc);
     not_nan_loc = not_nan_loc(not_nan_loc_loc);

    y = (water_quality_parameter_values(not_nan_loc))';  
    x = forty_deg_red_reflectance_numer(not_nan_loc)';
    p = polyfit(x,y,order);
    y1= polyval(p,x);
    
    if(verbose == 1)
     figure(30); subplot(311); scatter(x,y); hold on; scatter(x,y1)
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
    figure(30); subplot(312);
    plot((water_quality_parameter_values(not_nan_loc))/max((water_quality_parameter_values(not_nan_loc))));  xlabel('Data timeseries'); ylabel('normalized concentration');
end


clear forty_deg_red_reflectance_numer forty_deg_reflectance_denom 
for(i=indices_validation)
    if(size(all_images(i).(set_num_name_rrs_array(set_num)),1) >400 &&  size(all_images(i).(set_num_name_rrs_array(set_num)),2) >400)
        value_numer = nanmean(nanmean(all_images(i).(set_num_name_rrs_array(set_num))(wavelength_range_numer,spatial_pixel_range)));
        value_denom = nanmean(nanmean(all_images(i).(set_num_name_rrs_array(set_num))(wavelength_range_denom,spatial_pixel_range)));

        if( ((noise_min < value_numer) && (value_numer < noise_cutoff)) && ((noise_min < value_denom) && (value_denom < noise_cutoff)))
            forty_deg_red_reflectance_numer(i,:) = value_numer;
                        forty_deg_reflectance_denom(i,:) = value_denom;
            not_nan_loc_finder(i) = 1;
        else
            forty_deg_red_reflectance_numer(i,:) = NaN;          
            forty_deg_reflectance_denom(i,:) = NaN;
            not_nan_loc_finder(i) = 0;
        end
    end
end

forty_deg_red_reflectance_numer = forty_deg_red_reflectance_numer./forty_deg_reflectance_denom;

not_nan_loc = unique(find(not_nan_loc_finder)); %(~isnan(forty_deg_red_reflectance_numer(indices_validation)));    % +tss_range_end -1;

%not_nan_loc = find(~isnan(forty_deg_red_reflectance) );
 not_nan_loc = not_nan_loc'.* (   (hour(water_quality_parameter_times(not_nan_loc)) == 12)   );
 not_nan_loc_loc = find(not_nan_loc);
 not_nan_loc = not_nan_loc(not_nan_loc_loc);
 
%not_nan_loc_holder (index,:) = not_nan_loc;

red_green_ratio_pred = forty_deg_red_reflectance_numer(not_nan_loc)';  
predicted_beta = polyval(p,red_green_ratio_pred);

predicted_beta_holder(index,:) = polyval(p,red_green_ratio_pred);

if(verbose == 1)
    figure(30);subplot(313);scatter(red_green_ratio_pred,(water_quality_parameter_values(not_nan_loc))); hold on; scatter(red_green_ratio_pred,predicted_beta); xlabel('Ratio'); ylabel('predicted and actual'); legend('Actual conc', 'Predicted Conc')
end

% get correlation coefficient of Ratio and the NTU
% get statistics
ccr = corrcoef(forty_deg_red_reflectance_numer(not_nan_loc),(water_quality_parameter_values(not_nan_loc)));
corr_coeff_ratio(index) = ccr(1,2);
RMSE(index) = sqrt(mean((water_quality_parameter_values(not_nan_loc) - predicted_beta').^2));
MAPE(index) = 100*mean(abs((water_quality_parameter_values(not_nan_loc)- predicted_beta')./water_quality_parameter_values(not_nan_loc)));%mean absolute percentage error
angle_center(index) = 36.25+angle_increment;
    index=index+1;

end

[min_MAPE_val,min_MAPE_loc] = min(MAPE);
%max value of beta
predicted_beta = predicted_beta_holder(min_MAPE_loc,:);
%not_nan_loc = not_nan_loc_holder(min_MAPE_loc,:);


end