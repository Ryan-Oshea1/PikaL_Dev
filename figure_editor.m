%figure editor
%thsi code is used to standardize alll 1 figure and 2 figure images in my 
% this code will plot all of the individual figures
close all
%standardization parameters
fontsize_standard = 20;
black =[0 0 0];
dark_grey = [.5 .5 .5];
light_grey = [.75 .75 .75];
white = [.99 .99 .99];
medium_grey = [.6 .6 .6];

darker_grey = [.4 .4 .4];

correlation_coefficient_view_angle = 1;
remotely_vs_in_situ = 2;
adc_figure = 3 ;
snr_vs_wavelength = 4;
chl_beta_timeseries = 5;
chl_timeseries = 6;
spectral_snr_fig = 7;
median_snr_fig = 8 ;

plot_figure_number = correlation_coefficient_view_angle;
%%
close all
%plot figure correlation_coefficieny_view_angle
if(plot_figure_number == correlation_coefficient_view_angle)
    %perhaps add load the necesary data here
current_figure = figure(44);
set(gcf, 'units','points','position', [ 0 100 500 850])
% rows, columns, gap, marginh lower upper, wargin w left right)
ha = tight_subplot(2,1,[.05 .4],[.1 .013],[.18 .01])
          for ii = 1:2; axes(ha(ii));  end
         set(ha(1:2),'XTickLabel',''); set(ha,'YTickLabel','')
          
          
figure_size = get(current_figure, 'position')
figure_width = figure_size(3);
figure_height = figure_size(4);
%unpolarized = k, polarized = green, srg is red, wlr has x's, pol srg is
%magenta
axes(ha(1));

chlorophyll_text_position = [.181 .947 .26 .04];
annotation('rectangle',chlorophyll_text_position,'FaceColor','white')
an = annotation('textbox',chlorophyll_text_position,'String','Clorophyll a','EdgeColor','black','fontsize',fontsize_standard,'FitBoxToText','off','FontWeight','normal')

%subplot_1 = subplot(211)
hold on; plot(angle_center_chl_pol,corr_coeff_ratio_chl_pol_median,'Color',dark_grey,'LineWidth',3);  hold off;
hold on; plot(angle_center_chl_srg,corr_coeff_ratio_chl_srg_median,'Color',black,'LineWidth',3); hold off;
hold on; plot(angle_center_chl,corr_coeff_ratio_chl_median,'Color',light_grey,'LineWidth',3);  hold off;
hold on; scatter(wlr_angle,corr_coeff_ratio_wlr_chl_median,3000,light_grey,'.'); hold off;
hold on; scatter(wlr_angle,corr_coeff_ratio_wlr_chl_srg_median,3000,black,'.'); hold off;
hold on; scatter(wlr_angle,corr_coeff_ratio_wlr_chl_median,230,black,'o'); hold off;
hold on; scatter(wlr_angle,corr_coeff_ratio_wlr_chl_srg_median,230,white,'o'); hold off;

%plotting labels
%ylabel('Correlation Coefficient (R^2)')
%title({'Correlation Coefficient Between Remotely-Estimated' ' Chl a & Fluorometer-Measured Chl a'})
%xlabel('View Angle \theta_v (degrees)')
set(gca,'fontsize',fontsize_standard)
set(gca,'XLim',[39.5 53.5],'XTick',40:2:53)
%set(gca, 'XTick',[])
grid on
ax = gca
ax.GridLineStyle = '-'
ax.GridColor = 'k'
ax.GridAlpha = .4
%set(subplot_1,'position', [l, b , h, w])

%%%Beta
%subplot(212)
axes(ha(2));
beta_text_position = [.181 .478 .09 .04];
annotation('rectangle',beta_text_position,'FaceColor','white')
annotation('textbox',beta_text_position,'interpreter','latex','String','$${\beta}_{650}$$','EdgeColor','black','fontsize',fontsize_standard,'FitBoxToText','off','FontWeight','bold')

hold on; plot(angle_center_chl_pol,corr_coeff_ratio_beta_pol_median,'Color',dark_grey,'LineWidth',3);  hold off;
hold on; plot(angle_center_chl_srg,corr_coeff_ratio_beta_srg_median,'Color',black,'LineWidth',3); hold off;
hold on; plot(angle_center_chl,corr_coeff_ratio_beta_median,'Color',light_grey,'LineWidth',3);   hold off;
hold on; scatter(wlr_angle,corr_coeff_ratio_wlr_beta_median,3000,light_grey,'.'); hold off;
hold on; scatter(wlr_angle,corr_coeff_ratio_wlr_beta_srg_median,3000,black,'.'); hold off;
hold on; scatter(wlr_angle,corr_coeff_ratio_wlr_beta_median,230,black,'o'); hold off;
hold on; scatter(wlr_angle,corr_coeff_ratio_wlr_beta_srg_median,230,white,'o'); hold off;
%plotting labels
%ylabel('Correlation Coefficient (R^2)')
%title({'Correlation Coefficient Between Remotely-Estimated' ' Beta_{650} & Fluorometer-Measured Beta_{650}'})
grid on
xlabel(['View Angle \theta_v (' char(0176) ')'])
set(gca,'fontsize',fontsize_standard)
set(gca,'XLim',[39.5 53.5],'XTick',40:2:53)
%for the legend
%legendmarkeradjust(50) 
legend({'Camera Uncorrected' ,'Camera Polarized', 'Camera SSGC',  'Radiometer Raw', 'Radiometer SSGC'},'Position', [.23 .13 .3 .1]) %'Chl a Pol. Srg', '\beta_7_0_0 Pol. Srg',
legendmarkeradjust(50)
hlegend = findobj(gcf,'Type','Legend')
set(hlegend,'Position', [.26 .15 .3 .1])


ax = gca;
ax.GridLineStyle = '-';
ax.GridColor = 'k';
ax.GridAlpha = .4;
pos1 = get(ha(1),'position');
pos2 = get(ha(2),'position');
pos2(1) = pos2(1)-.07;
pos= [pos2(1:3) pos1(2)+pos1(4)-pos2(2)]
hAOuter=axes('position',pos,'visible','off');
hL=ylabel(hAOuter,'Correlation Coefficient (R^2)','fontsize',fontsize_standard,'visible','on');

 yticklabels(ha(1),'auto')
 xticklabels(ha(1),'auto')
 
  yticklabels(ha(2),'auto')
 xticklabels(ha(2),'auto')
set(gca,'FontName','Times')
%suplabel('super label' , 'y')
end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the next figure which is the dual plots of beta and chlorophyll a
if(plot_figure_number == remotely_vs_in_situ)
current_figure = figure(47);
set(gcf, 'units','points','position', [ 0 100 500 850])

ha = tight_subplot(2,1,[.05 .4],[.1 .013],[.18 .05])
          for ii = 1:2; axes(ha(ii));  end
         %set(ha(1:2),'XTickLabel',''); set(ha,'YTickLabel','')
 
axes(ha(1));
chlorophyll_text_position = [.181 .947 .404 .04];
annotation('rectangle',chlorophyll_text_position,'FaceColor','white')
an = annotation('textbox',chlorophyll_text_position,'String','Chlorophyll a (\mug/l)','EdgeColor','black','fontsize',fontsize_standard,'FitBoxToText','off','FontWeight','normal')

hold on


scatter(nan_median_daily_chl(day_range,truth_val),nan_median_daily_chl(day_range,median_value),3000,light_grey,'.')

scatter(nan_median_daily_chl_pol(day_range,truth_val),nan_median_daily_chl_pol(day_range,median_value),2000,dark_grey,'.')

scatter(nan_median_daily_chl_srg(day_range,truth_val),nan_median_daily_chl_srg(day_range,median_value),800,'.k')

disp(['  Min. MAPE Unpol. CHL: ' num2str(min(MAPE_chl)) '   Min. MAPE Pol.: ' num2str(min(MAPE_chl_pol)) '   Min. MAPE Srg.: ' num2str(min(MAPE_chl_srg))  '   Min. MAPE WLR: ' num2str(min(MAPE_wlr_chl)) '   Min. MAPE WLR Srg: ' num2str(min(MAPE_wlr_chl_srg)) ]); 

set(gca,'xscale','log')
set(gca,'yscale','log')

hold on; plot(0.4:1:20,0.4:1:20,'k','LineWidth',2); hold off;
hold on; plot(.4:1:20,(.4:1:20)*1.35,'-.k','LineWidth',2); hold off;
hold on; plot(.6:1:20,(.6:1:20)*.65,'-.k','LineWidth',2); hold off;

%legend('Camera Uncorrected','Camera Polarized', 'Camera SSGC',   '1:1','1:1+/-35%'); %'Radiometer Uncorrected','Radiometer SSGC',
grid on; %'Pol. Srg'
set(gca,'FontSize',fontsize_standard)
axis([.55 14 .55 22])
ax = gca;
ax.GridLineStyle = '-';
ax.GridColor = 'k';
ax.GridAlpha = .8;
ax.MinorGridAlpha = 1;
box on
%xticks([.5 5 10 20])
%set(ha(1:2),'XTickLabel','auto')  
 yticklabels(ha(1),'auto')
 xticklabels(ha(1),'auto')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_with_spaces = ['$${\beta}_{650}{\space}(1/sr)$$']
axes(ha(2));
beta_text_position = [.181 .478 .188 .04];
annotation('rectangle',beta_text_position,'FaceColor','white')
annotation('textbox',beta_text_position,'interpreter','latex','String',beta_with_spaces,'EdgeColor','black','fontsize',fontsize_standard,'FitBoxToText','off','FontWeight','bold')

scatter(nan_median_daily_beta(day_range,truth_val),nan_median_daily_beta(day_range,median_value),3000,light_grey,'.')
hold on
scatter(nan_median_daily_beta_pol(day_range,truth_val),nan_median_daily_beta_pol(day_range,median_value),2000,dark_grey,'.')
scatter(nan_median_daily_beta_srg(day_range,truth_val),nan_median_daily_beta_srg(day_range,median_value),800,'.k')
hold off
hold on; plot(0.0006:.001:.01,.0006:.001:.01,'k','LineWidth',2); hold off;
hold on; plot(0.0006:.001:.01,(.0006:.001:.01)*1.35,'-.k','LineWidth',2); hold off;
hold on; plot(0.0006:.001:.011,(.0006:.001:.011)*.65,'-.k','LineWidth',2); hold off;

set(gca,'xscale','log')
set(gca,'yscale','log')
axis([0.0006 .01 0.0004 .01])
grid on; %'Pol. Srg'
%title({['Remotely-Estimated vs. In-situ Measured \Beta_{650} (1/sr)']})%

%ylabel('Remotely-Estimated'); 
xlabel('In-situ Measured'); 
set(gca,'FontSize',fontsize_standard)
legend('Camera Uncorrected','Camera Polarized', 'Camera SSGC',   '1:1','1:1+/-35%'); %'Radiometer Uncorrected','Radiometer SSGC',
legendmarkeradjust(50)
hlegend = findobj(gcf,'Type','Legend')
set(hlegend,'Position', [.57 .13 .3 .1])

disp(['  Min. MAPE Unpol BETA.: ' num2str(min(MAPE_beta)) '   Min. MAPE Pol.: ' num2str(min(MAPE_beta_pol)) '   Min. MAPE Srg.: ' num2str(min(MAPE_beta_srg))  '   Min. MAPE WLR: ' num2str(min(MAPE_wlr_beta)) '   Min. MAPE WLR Srg: ' num2str(min(MAPE_wlr_beta_srg)) ]); %'   Min. MAPE Pol. Srg.: ' num2str(min(MAPE_chl_pol_srg))
ax = gca;
ax.GridLineStyle = '-';
ax.GridColor = 'k';
ax.GridAlpha = .8;
ax.MinorGridAlpha = 1;

pos1 = get(ha(1),'position');
pos2 = get(ha(2),'position');
pos2(1) = pos2(1)-.07;
pos= [pos2(1:3) pos1(2) + pos1(4)-pos2(2)];
hAOuter=axes('position',pos,'visible','off');
hL=ylabel(hAOuter,'Remotely-Estimated','fontsize',fontsize_standard,'visible','on');

set(gca,'FontName','Times')
end


%make the adc plot 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(plot_figure_number == adc_figure)
current_figure = figure(48);
set(gcf, 'units','points','position', [ 0 100 500 425])
scatter(mean_data_all,variance_data_all,[],dark_grey)
xlabel('Mean (counts)')
ylabel('Variance (counts^2)')

coefficients = polyfit(mean_data_all,variance_data_all,1)
hold on
plot(polyval(coefficients,1:2000),'k','LineWidth',3)
hold off

gain = coefficients(1)
offset = coefficients(2)
representable_electrons = gain*4095
txt = "Y=" + num2str(gain) + "*X + " +  num2str(offset) 
txt_gain = "1/(ADC Gain)=" + num2str(gain)
ADC_gain =1/gain
txt_adc_gain = "ADC Gain=" + num2str(ADC_gain)
text(25,230,txt,'FontSize',fontsize_standard)
text(25,210,txt_gain,'FontSize',fontsize_standard)
text(25,190,txt_adc_gain,'FontSize',fontsize_standard)
legend('Camera Measured','First Order Regression')

set(gca,'FontSize',fontsize_standard)
axis([0 1900 0 240])
%grid on
ax = gca;

ax.GridLineStyle = '-';
ax.GridColor = 'k';
ax.GridAlpha = .8;
legendmarkeradjust(15)
hlegend = findobj(gcf,'Type','Legend')
set(hlegend,'Position', [.5 .17 .3 .1])
box on

set(gca,'FontName','Times')
end

%make the SNR vs. wavelength plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(plot_figure_number == snr_vs_wavelength)
    openfig('/home/flagg/Dropbox (MIT)/PikaL_Dev/CodeBase/matlab/current_data_compiler/SNR_vs_wavelength.fig')
    set(gcf, 'units','points','position', [ 0 100 500 425])

axis([390 1080 0 115])
legend_1 = legend('Theo.',['Theo./$$\sqrt{2}$$'], 'Spa. 0','Spa. 100','Spa. 1')%'Rad. Spatial SNR 30'
set(legend_1,'Interpreter', 'latex')
set(gca,'FontSize',fontsize_standard)
box on

hlegend = findobj(gcf,'Type','Legend')
set(hlegend,'Position', [.69 .74 .2 .1])
title('')
ylabel('Signal-to-noise Ratio')

grid on
ax = gca;

ax.GridLineStyle = '-';
ax.GridColor = 'k';
ax.GridAlpha = .8;

%theoretical_SNR_text_position = [.63 .75 .404 .08];
%annotation('rectangle',theoretical_SNR_text_position,'FaceColor','white')
%an = annotation('textbox',theoretical_SNR_text_position,'String','Theor. SNR','EdgeColor','none','fontsize',fontsize_standard,'FitBoxToText','off','FontWeight','normal')


set(gca,'FontName','Times')

end







%plot the Timeseries of chl and beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(plot_figure_number == chl_beta_timeseries)
figure(14)
set(gcf, 'units','points','position', [ 0 100 500 425])
% rows, columns, gap, marginh lower upper, wargin w left right)
ha = tight_subplot(2,1,[.05 .4],[.1 .023],[.19 .01])
          for ii = 1:2; axes(ha(ii));  end
          
axes(ha(1));



hold on
plot(depth_datenums_chl,depth_chl,'Color',dark_grey,'lineWidth',1.25); datetickzoom
plot(surface_datenums,surface_median_chl,'Color',black,'lineWidth',1); datetickzoom

hold on
%
ylabel({['Chlorophyll a'] [' (\mug/l)']})
set(gca,'FontSize',fontsize_standard)
legend('1 m', '5 m')
ylim([0 25])
%datetick
     set(ha(1),'XTickLabel','');
     
xlim([min(depth_datenums_chl) max(depth_datenums_chl)])
%
 yticklabels(ha(1),'auto')
 %xticklabels(ha(1),'auto')




%%%%%%%

axes(ha(2));
hold on
plot(depth_datenums_b, depth_b,'Color',dark_grey,'lineWidth',.25)
plot(surface_datenums, surface_median_backscatter,'Color',black,'lineWidth',.5)

datetick
ylabel('\beta (1/sr)')
l = legend('1 m (650 nm)','5 m (700 nm)')
set(gca,'fontsize',fontsize_standard)
set(gca,'XLim',[min(chl_datenums_depth_and_surface) max(chl_datenums_depth_and_surface)])
datetick

ylim([0 .02])
xlim([min(depth_datenums_chl) max(depth_datenums_chl)])
%xticks([min(surface_datenums):60:max(surface_datenums)])
%legendmarkeradjust(50)
 yticklabels(ha(2),'auto')



set(gca,'FontName','Times')

    
end


%plot the Timeseries of chl and beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(plot_figure_number == chl_timeseries)
figure(15)
set(gcf, 'units','points','position', [ 0 100 500 425])
% % rows, columns, gap, marginh lower upper, wargin w left right)
% ha = tight_subplot(1,1,[.05 .4],[.1 .023],[.19 .01])
%           for ii = 1; axes(ha(ii));  end
          original_data_figure = plot(chl_datenums_depth_and_surface,chl_values_depth_and_surface,'Color',light_grey,'LineWidth',1)
hold on
scatter_points_figure = scatter(unquenched_times_depth_and_surface,interpolated_chl_unquenched_depth_and_surface,1000,dark_grey,'.')

interpolated_values_figure = plot(chl_datenums_depth_and_surface,interpolated_chl_depth_and_surface_values,'Color',black,'LineWidth',2)
datetickzoom
%title('Fluorometer Measured Chlorophyll a')
ylabel('Chl a (\mug/l)')
legend_value = legend('Measured Chl a','Unquenched Samples', 'Interpolated Chl a')
%legendmarkeradjust(50)
%  hl = findobj(legend_value,'type','line');
% set(hl,'LineWidth',2);
% ht = findobj(gcf,'type','text')
%  set(ht,'FontSize',30);
% legendmarkeradjust(30)
% hlegend = findobj(gcf,'Type','Legend')



set(gca,'fontsize',fontsize_standard)
set(gca,'XLim',[min(chl_datenums_depth_and_surface) max(chl_datenums_depth_and_surface)])

%xticks([min(chl_datenums_depth_and_surface):15:max(chl_datenums_depth_and_surface)])

%xticklabels([min(chl_datenums_depth_and_surface):7:max(chl_datenums_depth_and_surface)])

%set(gca,'YLim',[min(chl_values_depth_and_surface) max(chl_values_depth_and_surface)],'XTick',[min(chl_values_depth_and_surface):5:max(chl_values_depth_and_surface)])
%datetick
ylim([0 17])
legendmarkeradjust(50)
  
set(gca,'FontName','Times')
 
end



 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %plot typical radiance above typical SNR
      if(plot_figure_number == spectral_snr_fig)
          load('/home/flagg/Dropbox (MIT)/PikaL_Dev/CodeBase/matlab/current_data_compiler/real_world_radiance_snr.mat')
       figure(9);
       set(gcf, 'units','points','position', [ 0 100 500 425])
       % rows, columns, gap, marginh lower upper, wargin w left right)
       ha = tight_subplot(2,1,[.05 .2],[.15 .03],[.19 .06])% rows, columns, gap, marginh lower upper, wargin w left right)

           for ii = 1; axes(ha(ii));  end
      axes(ha(1));
       plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),radiance_holder(3,:),'Color',black,'LineWidth',3); 
      hold on;
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),radiance_holder(1,:),'Color',darker_grey,'LineWidth',3); 
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),radiance_holder(4,:),'Color',medium_grey,'LineWidth',3); 

      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),radiance_holder(2,:),'Color',light_grey,'LineWidth',3); 
                        
      %title('Spectral Radiance')
      %xlabel('Wavelength (nm)')
      ylabel({'Radiance'  '(W/m^2/um/sr)'})
      legend({'Dark Low Chl.','Dark High Chl.','Bright Low Chl.', 'Bright High Chl.'},'Position', [.605 .807 .3 .1])%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')

      set(gca,'FontSize',20)
      axis([380 1000 0 12])
             xticklabels(ha(1),'')
box on

      axes(ha(2));
     

      hold on
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),theoretical_snr_holder(1,:)/sqrt(2),'Color',darker_grey,'LineWidth',3); 
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_snr_holder(1,:),'-.','Color',darker_grey,'LineWidth',3); 

      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),theoretical_snr_holder(2,:)/sqrt(2),'Color',light_grey,'LineWidth',3); 
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_snr_holder(2,:),'-.','Color',light_grey,'LineWidth',3); 
      
      
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),theoretical_snr_holder(3,:)/sqrt(2),'Color',black,'LineWidth',3); 
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_snr_holder(3,:),'-.','Color',black,'LineWidth',3); 
      
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),theoretical_snr_holder(4,:)/sqrt(2),'Color',medium_grey,'LineWidth',3); 
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),spatial_snr_holder(4,:),'-.','Color',medium_grey,'LineWidth',3); 
      legend({'Shot-noise limited', 'Spatially derived'},'Position', [.58 .42 .3 .1])%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')
      %title('Spectral SNR')
      xlabel('Wavelength (nm)')
      ylabel({'Signal-to-noise'  'Ratio'})
      set(gca,'FontSize',fontsize_standard)
      xticklabels(ha(2),'auto')
      yticklabels(ha(2),'auto')
      axis([380 1000 0 85])
box on
%legend({'Camera Uncorrected' ,'Camera Polarized', 'Camera SSGC',  'Radiometer Raw', 'Radiometer SSGC'},'Position', [.23 .13 .3 .1]) %'Chl a Pol. Srg', '\beta_7_0_0 Pol. Srg',

%hlegend = findobj(gcf,'Type','Legend')
%set(hlegend,'Position', [.26 .15 .3 .1])
set(gca,'FontName','Times')

      end
      
      
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if(plot_figure_number == median_snr_fig )
      %plot median SNR's
            figure(7);
                   set(gcf, 'units','points','position', [ 0 100 500 425])

      load('/home/flagg/Dropbox (MIT)/PikaL_Dev/CodeBase/matlab/current_data_compiler/real_world_median_snr.mat')
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),nanmedian(spatial_SNR_uncorrected,2),'k-.','LineWidth',3);
      hold on;
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),nanmedian(spatial_SNR_polarized,2),'-.','Color',medium_grey,'LineWidth',3);
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),nanmedian(spatial_SNR_lee,2),'-.','Color',light_grey,'LineWidth',3);
      
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),nanmedian(theoretical_snr_uncorrected,2),'k','LineWidth',3); 
      hold on;
      plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),nanmedian(theoretical_snr_polarized,2),'Color',medium_grey,'LineWidth',3);
      %plot(camera_wavelengths(camera_wavelength_row,camera_wavelength_columns),nanmedian(theoretical_snr_lee,2)); 
      


     % hold on ; plot(camera_wavelengths(1,:),median_SNR); hold off; 
     % plot(camera_wavelengths(1,:),mean_rad./STD_rad); hold off;
      legend({'Uncorr. spatial','Pol. spatial','Lee spatial','Unc. shot-noise','Pol. shot-noise'},'Position', [.578 .741 .3 .1])%'Rad. Spatial SNR 1','Rad. Spatial SNR 30','Rad. Spatial SNR 100')
      %title('Median SNR of Timeseries')
      xlabel('Wavelength (nm)')
      ylabel('Median Signal-to-noise Ratio')
      set(gca,'FontSize',20)
      hold off
      axis([400 1000 0 100 ])
      
      set(gca,'FontName','Times')

      end
