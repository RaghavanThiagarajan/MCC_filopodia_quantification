% Raghavan Thiagarajan, DanStem, Copenhagen, 30th September 2018 - Version-02, 06th July 2020

% This script obtains different kinds of data and plots them in 2D and 3D.
% The data obtained are: (1) Coordinates from the track file of Tri-cellular junctions; (2) The
% ROI coordinates of cell surface contour; (3) Image dimensions and the
% corresponding (4) Time points (5) Shifted ROI coordinates and their mean intensity.

% In version 02, following updates were done: (1) In those figures where the TCJ trajectories were marked, the trajectories were color coded for
% the distance between the correpsonding TCJ and cell tip; (2) The intensity of cortex was collected in the adjoining Fiji script; this
% cortex intensity was imported and used to normalise the filopodia intensity; (3) Therefore figures with and without cortex normalisation
% were generated for filopodia intensities; Along with this, cortex intensity was also plotted; (3) Other new plots: the mean intensity of filopodia along
% each contour was averaged and this vaule was plotted against the distance (between TCJ and cell tip); (4) a new index called "filopodia
% index" was generated as a propxy for the no. of filopodia; the filopodia index is the sum of "area under the curve" above "1" after normalising with cortex;
% this is further divided by the contour length; this filopdia index was plotted against the distance (between TCJ and cell tip).

clc
clear all;
close all;

% a new directory is created to save the plots
mkdir ./plots

% Loading Image dimensions 
Imgdim = './data/Image_dimensions_&_input_parameters.csv';
Imagedim = importdata(Imgdim);
Image_dime = Imagedim.data;

% Getting pixel values to convert values into micrometers
xpixval = Image_dime(1,4); % In microns
ypixval = Image_dime(1,5); % In microns
Timevall = Image_dime(1,6); % In sec
Timeval = Timevall / 60; % Converted to minutes

% Getting Image dimensions - to be used for normalization
Img_height = Image_dime(1,3);
Image_height = Img_height * ypixval;

% Loading the tracked data from the Tri-cellular junction on the left
TCJLeft = './TCJ_left.csv'; 
TCJleft=importdata(TCJLeft); 
TCJ_left = TCJleft.data; 
TCJleftx = TCJ_left(:,3); TCJlefty = TCJ_left(:,4); 
TCJ_leftx = TCJleftx.*xpixval ; TCJ_lefty = TCJlefty.*ypixval; % x and y coordinates
TCJlefty_normalised = Image_height -  TCJ_lefty; % Normalized w.r.to the Image height for proper representation since the XY quadrants are arranged differently in ImageJ
TCJlefty_normalised_guideforeyes = TCJlefty_normalised / max(TCJlefty_normalised); % this guide for eyes is used in those plots where the intensities are plotted instead of contour;
% in this case the normalisation allows the TCJ trajectories to be present without the need for the axis to be a "distance" axis

% This gives the total number of tracked points
sz=size(TCJ_leftx); 
sze = sz(1,1);

% Loading the tracked data from the Tri-cellular junction on the right
TCJRight = './TCJ_right.csv';
TCJright=importdata(TCJRight);
TCJ_right = TCJright.data;
TCJrightx = TCJ_right(:,3); TCJrighty = TCJ_right(:,4); 
TCJ_rightx = TCJrightx.*xpixval; TCJ_righty = TCJrighty.*ypixval; % x and y coordinates
TCJrighty_normalised = Image_height -  TCJ_righty; % Normalized w.r.to the Image height for proper representation since the XY quadrants are arranged differently in ImageJ
TCJrighty_normalised_guideforeyes = TCJrighty_normalised / max(TCJrighty_normalised);% this guide for eyes is used in those plots where the intensities are plotted instead of contour;
% in this case the normalisation allows the TCJ trajectories to be present without the need for the axis to be a "distance" axis


% Here the largest "y" coordinate between the left and right TCJ is found.
% This value is used as limits in some of the plots where the inherent plot limits are smaller than the "y" coordinate of TCJ left / TCJ right; Because of this, 
% some of the TCJ points were not getting plotted. So using the maximum value of these TCJ left and TCJ right as the limits helped in plotting all the points.
if max(TCJlefty_normalised) > max(TCJrighty_normalised)
    maxTCJheight = max(TCJlefty_normalised);
else
    maxTCJheight = max(TCJrighty_normalised);
end
    
% Obtaining time points separately for TCJ left and right just for safer side
Tim = TCJ_left (:,2); Time0 = Tim.*Timeval; T=Time0.*ones(sz); 
Tim1 = TCJ_right(:,2); Time1 = Tim1.*Timeval; T1=Time1.*ones(size(TCJ_rightx));

% Defining color ranges. 
% mymap_1 is used for plotting the filopodia intensities (both on contour and separately) while the filopodia intensity is NOT normalised by Cortex; the map is designed in such
% a way that top 40 % of the intensities are coded with Magenta and the rest are coded with Cyan
mymap_1 = [0 1 1
         0 1 1
         0 1 1
         1 0 1
         1 0 1];

% mymap_2 is used for plotting the filopodia intensities (both on contour and separately) while the filopodia intensity IS NORMALISED by Cortex; the map is designed in such
% a way that all the intensity values above "1" are coded with Magenta and the rest are coded with Cyan. Therefore any intensity value that is equal
% to or more than the cortex intensity is considered a filopodia
mymap_2 = [0 1 1
           1 0 1];

% Loading all coordinates for each ROI and obtaining separately the minima of Y coordinate and the
% corresponding X coordinate for each ROI
ROI_xy = './data/ROI_Coordinates_Cell_contour/Actual_&_shifted_ROI/Actual_&_shifted_ROIcoord_with_Meanint_';

% Loading the minima values of the contour which represent the contour tip
ROI_xyminn = './data/ROI_Coordinates_Cell_contour/Actual_ROI_x_&_y_min.csv'; 
ROIxyminn = importdata(ROI_xyminn); 
ROI_xymin = ROIxyminn.data; 

% Loading the contour length values
contour_len = './data/ROI_Coordinates_Cell_contour/contour_length.csv'; 
contour_leng = importdata(contour_len); 
contour_length = contour_leng.data; 
contour_length = contour_length(:,2);
contour_length = contour_length';


%%
% Creating axis handles for those figures using multiple color gradients.
% When color coding is used for multiple axes, then separate axis handles are to be assigned for each axes with the color coding. So
% the figures are created here with the first axis handle assigned to the axis color coding for filopodia intensity with a particular color map;
% subsequent axis handle for another axis with another colormap is created later down the script.

for i = 1:6
    figure(i)
    figaxes(i) = axes;
    if i < 5
        colormap(figaxes(i), mymap_1);
    else
        colormap(figaxes(i), mymap_2);       
    end
end

%%
% In this section, the following are done: 
% (1) the distances between the cell tip and TCJleft / TCJright are computed; 
% (2) distance between the left and right TCJs are computed; 
% (3) 3D plots of contour color coded with filopodia intensity; 3D plots of
% filopodia intensity and; 3D plot of cortex intensities are generated; 
% (4) "area under the curve" i.e. filopodia index is computed from the plot where the filopodia intensity normalised by cortex intensity is plotted.

Distanceleft = zeros(1, sze); Distanceright = zeros(1, sze); Time = zeros(1,sze); DistanceTCJ = zeros(1,sze);
TopShiftedROI_meanint_avg = zeros(1,sze); area_under_the_curve = zeros(1,sze);
for i = 1:sze
    ROIxx = ROI_xymin(i,2); ROIyy = ROI_xymin(i,3);
    x1 = ROIxx * xpixval; y1 = ROIyy * ypixval;    
    
    leftxx = TCJ_left(i,3); leftyy = TCJ_left(i,4);
    x2 = leftxx * xpixval; y2 = leftyy * ypixval;    
    
    rightxx = TCJ_right(i,3); rightyy = TCJ_right(i,4);
    x3 = rightxx * xpixval; y3 = rightyy * ypixval;
    
    % Getting the distance between the tri cellular junction points and the
    % tip of each ROI (cell contour)    
    Distanceleft(i) = sqrt(((x2-x1)^2)+((y2-y1)^2)); 
    Distanceright(i) = sqrt(((x3-x1)^2)+((y3-y1)^2)); 
    
    % the most minimal value of Distanceleft and Distanceright and the most maximal value of Distanceleft and Distanceright are obtained; 
    % this is used for specifiying the lower and upper limits of the color bar that indicates the distance colorcoding in the TCJ trajectories
    if min(Distanceleft) <= min(Distanceright)
        min_distance = min(Distanceleft);    
        if max(Distanceleft) <= max(Distanceright)
            max_distance = max(Distanceright);    
        else
            max_distance = max(Distanceleft);  
        end
    else
        min_distance = min(Distanceright);
    end
    
    % Getting the distance between the tricellular junctions
    DistanceTCJ(i) = sqrt(((x3-x2)^2)+((y3-y2)^2));
       
    % Getting the time
    Time(i) = i*Timeval;
    
    % Getting the original ROI coordinates and plotting the ROI contour on the
    % already plotted Tricellular junction tracks    
    ROIXY = strcat(ROI_xy, sprintf('%01d',i), '.csv'); 
    ROI_XY=importdata(ROIXY);
    ROI_XYcrd = ROI_XY.data; 
    ROI_Xcrd = ROI_XYcrd(:,2); ROI_Ycrd = ROI_XYcrd(:,3);
    ROI_Xcoord = ROI_Xcrd * xpixval; ROI_Ycoord = ROI_Ycrd * ypixval;
    ROIYcoord_normalised = Image_height - ROI_Ycoord;
    
    % Plotting only the tip of the cell contour as a overlap on the full
    % cell contour
    Maxindex = find(ROIYcoord_normalised == max(ROIYcoord_normalised));
    szx = size(Maxindex);
    szxe = szx(1,1);
    CellTip = zeros(szxe,2);
    for j = 1:szxe
        tt = Maxindex(j,1);
        CellTip(j,1) = ROI_Xcoord(tt,1);
        CellTip(j,2) = ROIYcoord_normalised(tt,1);
    end
    % Getting the time separately for the cell tip data
    Tim2=i*Timeval; T2_1=Tim2.*ones(size(CellTip(:,1))); T2_2=Tim2.*ones(size(ROI_Xcoord));
    
    % 3D plot of the actual cell contour and the tip of cell contour
    figure(1);
    plot3(ROI_Xcoord, T2_2, ROIYcoord_normalised, '-k','LineWidth', 0.8); % 3D plot of cell surface contour
    hold on
    plot3(CellTip(:,1), T2_1, CellTip(:,2), '-m','LineWidth', 0.8); % 3D plot of only the cell contour tip
    hold on
      
    
    % Getting the top (filopodia) and down (cortex) shifted ROI coordinates and plotting the shifted ROI contour; the contour is color coded for
    % the filopodia intensity
    % Downshifted contour coordinates for cortex
    DownShiftedROI_Xcrd = ROI_XYcrd(:,2); DownShiftedROI_Ycrd = ROI_XYcrd(:,6); 
    DownShiftedROI_Xcoord = DownShiftedROI_Xcrd * xpixval; 
    DownShiftedROI_Ycoord = DownShiftedROI_Ycrd * ypixval; 
    DownShiftedROIYcoord_normalised = Image_height - DownShiftedROI_Ycoord; % This normalisation is done since in the image, the y coordinates start at the top as zero and increase while going doesn
    % Downshifted contour intensity for cortex
    DownShiftedROI_meanint = ROI_XYcrd(:,7); 
    DownShiftedROI_meanint_normalised = DownShiftedROI_meanint / max(DownShiftedROI_meanint);
    
    % Topshifted contour for filopodia
    TopShiftedROI_Xcrd = ROI_XYcrd(:,2); TopShiftedROI_Ycrd = ROI_XYcrd(:,4); 
    TopShiftedROI_Xcoord = TopShiftedROI_Xcrd * xpixval;
    TopShiftedROI_Ycoord = TopShiftedROI_Ycrd * ypixval;
    TopShiftedROIYcoord_normalised = Image_height - TopShiftedROI_Ycoord; % This normalisation is done since in the image, the y coordinates start at the top as zero and increase while going doesn
    % Topshifted contour intensity for filopodia
    TopShiftedROI_meanint = ROI_XYcrd(:,5); 
    TopShiftedROI_meanint_normalised = TopShiftedROI_meanint / max(TopShiftedROI_meanint); 
    % Topshifted contour intensity for filopodia normalised by cortex intensity
    TopShiftedROI_meanint_cortexnormalised = TopShiftedROI_meanint / nanmean(DownShiftedROI_meanint); 
    % the filopodia intensity normalised by cortex could be further normalised if needed; for now plots are done without this normalisation
    %TopShiftedROI_meanint_cortexnormalised = TopShiftedROI_meanint_cortexnormalised / max(TopShiftedROI_meanint_cortexnormalised); 
    
    % Average of the intensity values of each topshifted contour i.e. averaged filopodia intensity across every contour
    % this value is normalised by the averaged cortex mean intensity of the corresponding contour
    TopShiftedROI_meanint_avg(1,i) = nanmean(TopShiftedROI_meanint) / nanmean(DownShiftedROI_meanint);
    
    % Time values obtained for every contour   
    T3=Tim2.*ones(size(TopShiftedROI_Xcoord));
    
    % 3D plot of topshifted contour where the contour is color coded for filopodia mean intensity
    % axis handle for this figure has been created already with the color map: mymap_1
    figure(2); 
    plot3(TopShiftedROI_Xcoord, T3, TopShiftedROIYcoord_normalised); 
    patch([TopShiftedROI_Xcoord' nan], [T3' nan], [TopShiftedROIYcoord_normalised' nan], [TopShiftedROI_meanint_normalised' nan], 'EdgeColor', 'interp');
    hold on
    
    % 3D plot of filopodia mean intensity
    % axis handle for this figure has been created already with the color map: mymap_1
    figure(3); 
    plot3(TopShiftedROI_Xcoord, T3, TopShiftedROI_meanint_normalised); 
    patch([TopShiftedROI_Xcoord' nan], [T3' nan], [TopShiftedROI_meanint_normalised' nan], [TopShiftedROI_meanint_normalised' nan], 'EdgeColor', 'interp');
    hold on   
    
    % 3D plot of cortex mean intensity
    % axis handle for this figure has been created already with the color map: mymap_1
    figure(4); 
    plot3(DownShiftedROI_Xcoord, T3, DownShiftedROI_meanint_normalised); 
    patch([DownShiftedROI_Xcoord' nan], [T3' nan], [DownShiftedROI_meanint_normalised' nan], [DownShiftedROI_meanint_normalised' nan], 'EdgeColor', 'interp');
    hold on
    
    % this is to get all the filopodia intensity values above "1" (after normalisation by cortex); and this is used for specific color coding
    % in the following two plots i.e. all intensity values above "1" are color coded as magenta and the rest are coded as cyan
    filopodia_int_marking = TopShiftedROI_meanint_cortexnormalised>1;
    
    % 3D plot of topshifted contour where the contour is color coded for filopodia mean intensity normalised by cortex intensity
    % axis handle for this figure has been created already with the color map: mymap_2
    figure(5); 
    plot3(TopShiftedROI_Xcoord, T3, TopShiftedROIYcoord_normalised);
    patch([TopShiftedROI_Xcoord' nan], [T3' nan], [TopShiftedROIYcoord_normalised' nan], [filopodia_int_marking' nan], 'EdgeColor', 'interp');
    hold on
    
    % 3D plot of filopodia mean intensity normalised by cortex intensity
    % axis handle for this figure has been created already with the color map: mymap_2
    figure(6); 
    plot3(TopShiftedROI_Xcoord, T3, TopShiftedROI_meanint_cortexnormalised);
    % from this plot, the sum of "area under the curve", above the values of "1" are obtained 
    area_under_the_curve(1,i) = sum(TopShiftedROI_meanint_cortexnormalised(TopShiftedROI_meanint_cortexnormalised>1)); % area under the curve
    patch([TopShiftedROI_Xcoord' nan], [T3' nan], [TopShiftedROI_meanint_cortexnormalised' nan], [filopodia_int_marking' nan], 'EdgeColor', 'interp');
    hold on
    
end

%%
% plotting the TCJ trajectories color coded for the distance between the correspoding TCJ and the cell tip; this needs to be plotted separately
% because the corresponding axis needs to be declared as a separate handle.
% the plots are made in the following way:
        % plots are made for TCJleft and TCJright trajectories separately; in each of these plots:
        % new handle for "z axis" corresponding to the height is made with the corresponding colormap so that this colormap can be superimposed 
        % on the already existing colormap for the filopodia intensity which exists in the same / corresponding figure number.

% a new colormap is defined to have smoothly transitioning colors while coding the distance and also to not overlap with the already existing colors
mymap_3 = [0.5 0.2 0.9
               0.3 0.6 1
               0.5 1 0.5
               1 1 0
               1 0.5 0];

for i = 1:6
     
    if (i==1) || (i==2) || (i==5) % only plots that have the contours are plotted here
        
        figure(i);
        zlim([0 maxTCJheight+1]); % the maximum height is given by the maximum 'Y' coordinate of TCJleft or TCJright so that all the points of contour, 
        % and the TCJ trajectories are included within the plot; this limit is fixed only for the plots with contour.
        
        % plotting TCJleft trajectory
        figaxesTCJleft = axes; colormap(figaxesTCJleft, mymap_3); % creating new axis handle
        plot3(TCJ_leftx, T, TCJlefty_normalised, '-ob'); grid on; 
        set(gca,'color','none') % making this axis handle transparent for superimposition with the previous axis handle of the same figure
        patch([TCJ_leftx' nan], [T' nan], [TCJlefty_normalised' nan], [Distanceleft nan],'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
        
        % plotting TCJright trajectory
        figaxesTCJright = axes; colormap(figaxesTCJright, mymap_3); % creating new axis handle
        plot3(TCJ_rightx, T1, TCJrighty_normalised, '-ob'); 
        set(gca,'color','none') % making this axis handle transparent for superimposition with the previous axis handle of the same figure
        patch([TCJ_rightx' nan], [T1' nan], [TCJrighty_normalised' nan], [Distanceright nan],'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
        
        % colorbar is linked to the appropriate handle and it is positioned properly with limits to include whole of TCJleft and TCJright values
        colorbar(figaxesTCJleft,'Location','East','Direction', 'reverse','Position',[.93 .11 .036 .815]); caxis([min_distance max_distance]); 
        
        % linking the axis handles (of filopodia intensity and TCJ trajectory) so that two axis handles with different color coding
        % can be displayed in the same plot
        Link = linkprop([figaxes(i), figaxesTCJleft, figaxesTCJright],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
        setappdata(gcf, 'StoreTheLink', Link);        
               
    else  % plots show intensities without the contour       
        figure(i);
        
        % plotting TCJleft trajectory
        fig3axesTCJleft = axes; colormap(fig3axesTCJleft, mymap_3); % creating new axis handle
        plot3(TCJ_leftx, T, TCJlefty_normalised_guideforeyes, '-ob'); grid on; 
        set(gca,'color','none') % making this axis handle transparent for superimposition with the previous axis handle of the same figure
        patch([TCJ_leftx' nan], [T' nan], [TCJlefty_normalised_guideforeyes' nan], [Distanceleft nan],'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
        hold on
        
        % plotting TCJright trajectory
        fig3axesTCJright = axes; colormap(fig3axesTCJright, mymap_3); % creating new axis handle
        plot3(TCJ_rightx, T1, TCJrighty_normalised_guideforeyes, '-ob'); 
        set(gca,'color','none') % making this axis handle transparent for superimposition with the previous axis handle of the same figure
        patch([TCJ_rightx' nan], [T1' nan], [TCJrighty_normalised_guideforeyes' nan], [Distanceright nan],'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
        hold on
        
        % colorbar is linked to the appropriate handle and it is positioned properly with limits to include whole of TCJleft and TCJright values
        colorbar(fig3axesTCJleft,'Location','East','Direction', 'reverse','Position',[.93 .11 .036 .815]); caxis([min_distance max_distance]); 
        
        % linking the axis handles (of filopodia intensity and TCJ trajectory) so that two axis handles with different color coding
        % can be displayed in the same plot
        Link = linkprop([figaxes(i), fig3axesTCJleft, fig3axesTCJright],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
        setappdata(gcf, 'StoreTheLink', Link);        
    end 
    
end
%%
% directory in which the plots are to be saved
folder_save = './plots';

% setting the figure axes labels, title and then saving
figure(1);
xlabel('Lateral distance [\mum]', 'FontSize', 18);  ylabel('Time [min]', 'FontSize', 18); zlabel('Height [\mum]', 'FontSize', 18);  
title({'Shift in cell contour';'(Magenta: tip of cell contour)'}, 'FontSize', 10);
saveas(gca, fullfile(folder_save, 'Shift_in_cell_contour'), 'fig');

figure(2);
xlabel('Lateral distance [\mum]', 'FontSize', 18);  ylabel('Time [min]', 'FontSize', 18); zlabel('Height [\mum]', 'FontSize', 18);  
title({'Cell contour color coded for filopodia mean intensity'; '(Magenta: top 40% intensity)'}, 'FontSize', 10);
saveas(gca, fullfile(folder_save, 'Shift_in_cell_contour_color_coded_for_filopodia_mean_intensity'), 'fig'); 

figure(3);
xlabel('Lateral distance [\mum]', 'FontSize', 18);  ylabel('Time [min]', 'FontSize', 18); zlabel('Mean intensity [a.u]', 'FontSize', 18);  
title({'Filopodia mean intensity'; '(Magenta: top 40% intensity)'}, 'FontSize', 10);
saveas(gca, fullfile(folder_save, 'Shift_in_filopodia_meanintensity'), 'fig'); 

figure(4);
xlabel('Lateral distance [\mum]', 'FontSize', 18);  ylabel('Time [min]', 'FontSize', 18); zlabel('Mean intensity [a.u]', 'FontSize', 18);  
title({'Mean intensity of cortex'; '(Magenta: top 40% intensity)'}, 'FontSize', 10);
saveas(gca, fullfile(folder_save, 'Cortex_meanintensity'), 'fig');

figure(5);
xlabel('Lateral distance [\mum]', 'FontSize', 14);  ylabel('Time [min]', 'FontSize', 14); zlabel('Height [\mum]', 'FontSize', 14);  
title({'Cell contour color coded for filopodia mean intensity normalised with cortex'; '(Magenta: intensity values more than "1", after normalisation)'}, 'FontSize', 8);
saveas(gca, fullfile(folder_save, 'Shift_in_cell_contour_color_coded_filopodia_meanintensity_cortex_normalised'), 'fig'); 

figure(6);
xlabel('Lateral distance [\mum]', 'FontSize', 14);  ylabel('Time [min]', 'FontSize', 14); zlabel('Mean int. (cortex normalised) [a.u]', 'FontSize', 14);  
title({'Filopodia mean intensity normalised with cortex'; '(Magenta: intensity values more than "1", after normalisation)'}, 'FontSize', 10);
saveas(gca, fullfile(folder_save, 'Shift_in_filopodia_meanintensity_cortex_normalised'), 'fig'); 


% 2D plot of the distance between the left and right tri-cellular junctions and the tip of cell surface contour (or y minima and the corresponding x)
figure(7);
plot (Time, Distanceleft, '-ob'); 
hold on
plot (Time, Distanceright, '-or'); 
xlabel('Time [min]', 'FontSize', 18);  ylabel('Distance [\mum]', 'FontSize', 18);
legend({'TCJ left','TCJ right'},'FontSize',18,'TextColor','black')
title('Distance between TCJ and cell surface', 'FontSize', 18);
saveas(gca, fullfile(folder_save,'Distance_between_TCJ_and_cell_surface'), 'fig'); 

% 2D plot of the distance between the left and right tri-cellular junctions
figure(8);
plot (Time, DistanceTCJ, '-ob'); 
xlabel('Time [min]', 'FontSize', 18);  ylabel('Distance [\mum]', 'FontSize', 18);
title('Distance between TCJs', 'FontSize', 18);
saveas(gca, fullfile(folder_save,'Distance_between_TCJ'), 'fig');

% 2D plot of the average of filopodia intensity along the contour against the distance between the cell tip and TCJ
figure(9);
% the average of filopodia contour intensity is normalised if needed
% TopShiftedROI_meanint_avg_norm = TopShiftedROI_meanint_avg / max(TopShiftedROI_meanint_avg);
scatter (Distanceleft, TopShiftedROI_meanint_avg, 'bo'); 
hold on 
scatter (Distanceright, TopShiftedROI_meanint_avg, 'ro'); 
xlim([0 max_distance+1]); % max limit is fixed by the max of maximal distance of TCJleft and TCJright
set(gca, 'XDir', 'reverse');
xlabel('Distance [\mum]', 'FontSize', 18);  ylabel('Filopodia mean intensity [a.u]', 'FontSize', 18);
title('Change in filopodia mean intensity during intercalation', 'FontSize', 18);
saveas(gca, fullfile(folder_save,'Filopodia_intensity_during_intercalation'), 'fig');

%2D plot of the filopodia index (sum of the area under the curve) against the distance between the cell tip and TCJ
figure(10);
% sum of the area under the curve is normalised by the corresponding contour length in order to avoid bias based on contour length
area_under_the_curve_div_by_contour = area_under_the_curve ./ contour_length;
% sum of the area under the curve is further normalised if needed
% area_under_the_curve_div_by_contour = area_under_the_curve_div_by_contour / max(area_under_the_curve_div_by_contour);
scatter (Distanceleft, area_under_the_curve_div_by_contour, 'bo'); 
hold on 
scatter (Distanceright, area_under_the_curve_div_by_contour, 'ro'); 
xlim([0 max_distance+1]); % max limit is fixed by the max of maximal distance of TCJleft and TCJright
set(gca, 'XDir', 'reverse');
xlabel('Distance [\mum]', 'FontSize', 18);  ylabel({'Filopodia index'; '(area under the curve) [a.u]'}, 'FontSize', 18);
title({'Filopodia index'; '(total no. of filpodia estimated from area under the curve)'}, 'FontSize', 18);
saveas(gca, fullfile(folder_save,'Filopodia_index_(area_under_the_curve)'), 'fig'); 

%---------------------------------------------------------%---------------------------------------------------------------------%---------------------------------------------------------------------


