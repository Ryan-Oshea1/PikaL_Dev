for (j= 1:4612)
    j
    if(find_usable_sza(water_quality_parameter_times_chl(j),105,165) == 1)
    %all_images(j).set_num_7_sza = all_images(j).set_num_7;
    else
       all_images(j).set_num_7 = []; 
    end
    
end