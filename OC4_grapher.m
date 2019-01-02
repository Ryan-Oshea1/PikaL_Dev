    oc4_graph = imread('oc4_plot.png');
    figure(37);
    oc4_graph = flipdim(oc4_graph,2);

     image('CData',oc4_graph,'XData',[logspace(log10(.01), log10(25))],'YData',[logspace(log10(.001), log10(100))])
     hold on  
     set(gca,'xscale','log'); set(gca,'yscale','log');
%%
    figure(37); 
    xlabel('Maximum Band Ratio')
    ylabel('Chlorophyll a (\mug/l)')
    legend(gca,'off')
    set(gca,'FontSize', 18)