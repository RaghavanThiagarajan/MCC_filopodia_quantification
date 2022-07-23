% Raghavan Thiagarajan, DanStem, Copenhagen, 30th September 2018 - Version-03, 10th September 2020

% This script obtains different kinds of data and plots them in 2D and 3D.
% The data obtained are: (1) Coordinates from the track file of Tri-cellular junctions; (2) The
% ROI coordinates of cell surface contour; (3) Image dimensions and the
% corresponding (4) Time points (5) Shifted ROI coordinates and their mean intensity.

% In version 03, following updates were done: (1) The filopodia index that cumulative across the cell was split into left and right filopodia index
% and this was plotted against time and distance between the cell tip and left and right TCJs respectively; (2) In addition to getting the distance
% between the cell tip and the left/right TCJs, the distance between the cell tip and "top" was obtained. "Top" refers to the line connecting the
% TCJs. This was to remove the bias coming from x/y movement that lead to increase/decrease w.r.to left/right TCJs; (3) The cumulative filodpodia index was
% plotted against the difference between the distances of cell tip to left/right TCJs; (4) The ratio of left to right filopodia index and ratio of left to right 
% TCJ distances were plotted against time.

% In version 02, following updates were done: (1) In those figures where the TCJ trajectories were marked, the trajectories were color coded for
% the distance between the correpsonding TCJ and cell tip; (2) The intensity of cortex was collected in the adjoining Fiji script; this
% cortex intensity was imported and used to normalise the filopodia intensity; (3) Therefore figures with and without cortex normalisation
% were generated for filopodia intensities; Along with this, cortex intensity was also plotted; (3) Other new plots: the mean intensity of filopodia along
% each contour was averaged and this vaule was plotted against the distance (between TCJ and cell tip); (4) a new index called "filopodia
% index" was generated as a proxy for the no. of filopodia; the filopodia index is the sum of "area under the curve" above "1" after normalising with cortex;
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
Image_width_orig = Image_dime(1,2);
Image_width = Image_width_orig * xpixval;
Img_height_orig = Image_dime(1,3);
Image_height = Img_height_orig * ypixval;

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
    x1 = ROIxx * xpixval; y1 = ROIyy * ypixval;  % celltipcoord11(i,1) = x1;  celltipcoord11(i,2) = y1;
    
    leftxx = TCJ_left(i,3); leftyy = TCJ_left(i,4);
    x2 = leftxx * xpixval; y2 = leftyy * ypixval;    
    
    rightxx = TCJ_right(i,3); rightyy = TCJ_right(i,4);
    x3 = rightxx * xpixval; y3 = rightyy * ypixval;
    
    % Getting the distance between the tri cellular junction points and the
    % tip of each ROI (cell contour)    
    Distanceleft(i) = sqrt(((x2-x1)^2)+((y2-y1)^2)); 
    Distanceright(i) = sqrt(((x3-x1)^2)+((y3-y1)^2)); 
    ratio_TCJ_distance(i,1) = Distanceleft(i) / Distanceright(i); % ratio in lengths between the cell tip and the TCJs on the left and right
    
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
    Tim2=i*Timeval; T2_1=Tim2.*ones(size(CellTip(:,1))); T2_2=Tim2.*ones(size(ROI_Xcoord)); % Timeplot(i,1) = Tim2;
    
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
    
    %-----------------------------------------------------%---------------------------------------------------------%
    % start of "finding the distance between the cell tip and top"
    % here "top" indicates the line connecting the ends of the image across width, where the line passes through the two TCJ. In order to draw the line between the two TCJ, the x-coordinates of the two TCJ are taken
    % and they are linearly interpolated using "linspace". The same is done for y-coordinates. In order to extend the lines to the ends of images: for the left end of the image - x-coorindate is taken as zero and the
    % y-coordinate of left TCJ is taken; similarly for the right end of the image, the x-coordinate is taken as the last pixel of the image in "x" or the "image_width" and the y-coordinate of right TCJ is taken. Then linear
    % inrepolation is done between coordinates at the three sections: (1) at x=0 and left TCJ; (2) coordinates at left TCJ and right TCJ; (3) and right TCJ and coordinates at "x=image_width". To do this interpolation, the
    % number of points between each of these coordinates (or the spacing between any two points) is obtained in the "pix_incre_..." variable. This is computed separately for each section. This spacing will match the
    % spacing in the coordinates of the contour.
    pix_incre_start = x2/xpixval;
    pix_incre_mid = (x3/xpixval) - pix_incre_start;
    pix_incre_end = (Image_width/xpixval) - (pix_incre_start + pix_incre_mid);
    % then the line is drawn for all these sections by interpolation and based
    % on the spacing obtained above.
    xstart = linspace(1, x2, pix_incre_start);
    xmid = linspace(x2, x3, pix_incre_mid);
    xend = linspace(x3, Image_width, pix_incre_end);
    % Then these sections are combined
    xwhole = [xstart, xmid, xend];
    xwhole1 = floor(xwhole);
    % y-coordinates are obtained in the same manner
    ystart = y2 * ones(size(xstart));
    ymid = linspace(y2, y3, pix_incre_mid);
    yend = y3 * ones(size(xend));
    ywhole = [ystart, ymid, yend];
    ywhole1 = floor(ywhole);
    % at this point, we have the x and y coordinates of line connecting the image width that passes through the left and right TCJ
    % Now we take the celltip coordinates "(x1,y1)" and find the corresponding coordinates in this line
    idxsearch = find(xwhole1 == floor(x1));
    idxsearch_center_index = ceil(length(idxsearch)/2);
    idxsearch_center = idxsearch(idxsearch_center_index);
    top_x = xwhole(idxsearch_center);
    top_y = ywhole(idxsearch_center);
    % now we find the distance between the cell tip and the top
    celltip_top_distance(i,1) = sqrt((top_x - x1)^2 + (top_y - y1)^2);
    % end of "finding the distance between the cell tip and top"
    %-----------------------------------------------------%---------------------------------------------------------%
    
    %-----------------------------------------------------%---------------------------------------------------------%
    % Start of Getting the data for left - right plot
    % This is part of plotting the index of left and right filopodia index separately. For this, first we are splitting the ROI coordinates into left and right; In order to do this, the index of the midpoint 
    % (tip of the contour) of ROIcoordinate arrays i.e. (ROI_Xcrd & ROI_Ycrd) are used. Index of the elements starting from the first element to midpoint is taken as left and index of elements from after the midpoint until
    % the end of the array is taken as right. 
    celltip_center = ceil(length(Maxindex)/2);
    celltip_center_idx = Maxindex(celltip_center); 
    ROI_Xcoord_min = min(ROI_Xcoord); ROI_Xcoord_max = max(ROI_Xcoord);
    ROI_Xcoord_min_idx_list = find(ROI_Xcoord==ROI_Xcoord_min); ROI_Xcoord_max_idx_list = find(ROI_Xcoord==ROI_Xcoord_max);
    % At this point, "ROI_Xcoord_min_idx_list" and "ROI_Xcoord_max_idx_list" contain the list of indices corresponding to elements that are "min" and "max" of "ROI_Xcoord". A list is
    % stored here if there are more than one "min" and "max". Otherwise there is only one element. In case there are more than one item, then we need to find the "minimum" and "maximum" of the list depeneding on
    % how the user has drawn the ROI. -- The ROI coordinates, while obtained from the ROI in the FIJI script, are stored in the same order as the user has drawn it. For example, for a ROI starting from
    % x-coordinate 1 until x-coordinate 10, if the user draws it from left to right, then the x-coordinate list is stored as 1 to 10; if the user draws it from right to left then the x-coordinate list is stored as
    % 10 to 1.-- Then in those plots where the x-coordinates are plotted with directionality, the plot is going to look different depending on which the direction the user has drawn the ROI. This is the case while
    % plotting the left and right of filopodia index. If the user has drawn the ROI from left to right, the left filopodia index plot will contain the smaller / starting x-coordinates - and this is the proper direction i.e. the x
    % axis increases from left to right. Whereas if the user has drawn the ROI from right to left, then the left filopodia index plot will contain the larger / ending x-coordinates and the right plot will contain the smaller / 
    % starting x-coordinates - this is the opposite direction. Therefore the following "if loop" is used to make sure we always plot the smaller / starting x-coordinates in the left plot and the larger / ending x-coordinates in the right plot. 
    if max(ROI_Xcoord_min_idx_list) < min(ROI_Xcoord_max_idx_list) % normal - left to right direction
        ROI_Xcoord_min_idx = min(ROI_Xcoord_min_idx_list);
        ROI_Xcoord_max_idx = max(ROI_Xcoord_max_idx_list);    
        incre = 1;        
    else                % reverse - right to left direction
        ROI_Xcoord_min_idx = max(ROI_Xcoord_min_idx_list);
        ROI_Xcoord_max_idx = min(ROI_Xcoord_max_idx_list);
        incre = -1;        
    end      
    % Based on these indices, the filopodia mean intensities are split into left and right.
    TopShiftedROI_meanint_cortexnormalised_left = TopShiftedROI_meanint_cortexnormalised(ROI_Xcoord_min_idx:incre:celltip_center_idx, 1);
    TopShiftedROI_meanint_cortexnormalised_right = TopShiftedROI_meanint_cortexnormalised(celltip_center_idx+1:incre:ROI_Xcoord_max_idx, 1);
%     TopShiftedROI_Xcoord_left = TopShiftedROI_Xcoord(ROI_Xcoord_min_idx:incre:celltip_center_idx, 1);
%     TopShiftedROI_Xcoord_right = TopShiftedROI_Xcoord(celltip_center_idx+1:incre:ROI_Xcoord_max_idx, 1);
    T3_left = Tim2.*ones(size(TopShiftedROI_meanint_cortexnormalised_left));
    T3_right = Tim2.*ones(size(TopShiftedROI_meanint_cortexnormalised_right));
    % End of Getting the data for left - right plot
    %-----------------------------------------------------%---------------------------------------------------------%
    
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
    area_under_the_curve(1,i) = sum(TopShiftedROI_meanint_cortexnormalised(TopShiftedROI_meanint_cortexnormalised > 1)); % area under the curve
    patch([TopShiftedROI_Xcoord' nan], [T3' nan], [TopShiftedROI_meanint_cortexnormalised' nan], [filopodia_int_marking' nan], 'EdgeColor', 'interp');
    hold on
    
    % Figures(66) & (67) is only for obtaining the sum of "area under the curve" for left and right filopodia index and not for plotting
    figure(66);
    plot(T3_left, TopShiftedROI_meanint_cortexnormalised_left);
    area_under_the_curve_left(1,i) = sum(TopShiftedROI_meanint_cortexnormalised_left(TopShiftedROI_meanint_cortexnormalised_left > 1));
    hold on
    close (figure(66));
    figure(67);
    plot(T3_right, TopShiftedROI_meanint_cortexnormalised_right);
    area_under_the_curve_right(1,i) = sum(TopShiftedROI_meanint_cortexnormalised_right(TopShiftedROI_meanint_cortexnormalised_right > 1)); 
    hold on
    close (figure(67));
    
%     figure(68)
%     plot([x1,top_x],[y1, top_y]);  
%     xlim([0 Image_width]); ylim([0 Image_height]);
%     hold on
    
    
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
               
    else  % plots showing only intensities without the contour       
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
title('Distance between TCJ and cell tip', 'FontSize', 12);
saveas(gca, fullfile(folder_save,'Distance_between_TCJ_and_cell_surface'), 'fig'); 

% 2D plot of the distance between the left and right tri-cellular junctions
figure(8);
plot (Time, DistanceTCJ, '-ob'); 
xlabel('Time [min]', 'FontSize', 18);  ylabel('Distance [\mum]', 'FontSize', 18);
title('Distance between TCJs', 'FontSize', 12);
saveas(gca, fullfile(folder_save,'Distance_between_TCJ'), 'fig');

% 2D plot of the average of filopodia intensity along the contour against the distance between the cell tip and TCJ
figure(9);
% the average of filopodia contour intensity is normalised if needed
% TopShiftedROI_meanint_avg_norm = TopShiftedROI_meanint_avg / max(TopShiftedROI_meanint_avg);
plot (Distanceleft, TopShiftedROI_meanint_avg, 'bo'); 
hold on 
plot (Distanceright, TopShiftedROI_meanint_avg, 'ro'); 
xlim([0 max_distance+1]); % max limit is fixed by the max of maximal distance of TCJleft and TCJright
set(gca, 'XDir', 'reverse');
xlabel('Distance [\mum]', 'FontSize', 18);  ylabel('Filopodia mean intensity [a.u]', 'FontSize', 18);
legend({'Distance from cell tip to left TCJ','Distance from cell tip to right TCJ'},'FontSize',18,'TextColor','black')
title('Change in filopodia mean intensity during intercalation', 'FontSize', 12);
saveas(gca, fullfile(folder_save,'Filopodia_intensity_during_intercalation'), 'fig');

%2D plot of the filopodia index (sum of the area under the curve) against the distance between the cell tip and TCJ
figure(10);
% sum of the area under the curve is normalised by the corresponding contour length in order to avoid bias based on contour length
area_under_the_curve_div_by_contour = area_under_the_curve ./ contour_length;
% sum of the area under the curve is further normalised if needed
% area_under_the_curve_div_by_contour = area_under_the_curve_div_by_contour / max(area_under_the_curve_div_by_contour);
plot (Distanceleft, area_under_the_curve_div_by_contour, 'bo'); 
hold on 
plot (Distanceright, area_under_the_curve_div_by_contour, 'ro'); 
xlim([0 max_distance+1]); % max limit is fixed by the max of maximal distance of TCJleft and TCJright
set(gca, 'XDir', 'reverse');
xlabel('Distance [\mum]', 'FontSize', 18);  ylabel({'Cumulative filopodia index'; '(area under the curve) [a.u]'}, 'FontSize', 18);
legend({'Distance from cell tip to left TCJ','Distance from cell tip to right TCJ'},'FontSize',18,'TextColor','black')
title({'Filopodia index'; '(total no. of filopodia estimated from area under the curve)'}, 'FontSize', 12);
saveas(gca, fullfile(folder_save,'Filopodia_index_(area_under_the_curve)'), 'fig'); 

% 2D plot of filopodia index against the distance between the cell tip and top
figure(11);
plot(celltip_top_distance, area_under_the_curve_div_by_contour, 'ko');
set(gca, 'XDir', 'reverse');
xlabel('Distance (cell tip - line connecting TCJs) [\mum]', 'FontSize', 18);  ylabel({'Cumulative filopodia index'; '(area under the curve) [a.u]'}, 'FontSize', 18);
title('Filopodia index against the distance to top (line connecting the TCJs)', 'FontSize', 12);
saveas(gca, fullfile(folder_save,'Filopodia_index_vs_distance_to_top'), 'fig');

% 2D plot of filopodia index against time
figure(12);
plot(Time, area_under_the_curve_div_by_contour, 'ko-');
xlabel('Time [min]', 'FontSize', 18);  ylabel({'Cumulative filopodia index'; '(area under the curve) [a.u]'}, 'FontSize', 18);
title('Filopodia index against time)', 'FontSize', 12);
saveas(gca, fullfile(folder_save,'Filopodia_index_vs_Time'), 'fig');

% 2D plot of left and right filopodia index against time
figure(13);
% sum of area under the curve is normalised by the corresponding contour length in order to avoid bias based on contour length; int this case the contour length
% is half of contour length since only half of the contour is taken.
contour_length_half = contour_length/2;
area_under_the_curve_left = area_under_the_curve_left ./ contour_length_half;
area_under_the_curve_right = area_under_the_curve_right ./ contour_length_half;
plot (Time, area_under_the_curve_left, 'bo-');
hold on
plot (Time, area_under_the_curve_right, 'ro-');
xlabel('Time [min]', 'FontSize', 18);  ylabel({'Filopodia index'; '(area under the curve) [a.u]'}, 'FontSize', 18);
legend({'Left filopodia index','Right filopodia index'},'FontSize',18,'TextColor','black')
title({'Filopodia index left and right'; '(total no. of filopodia estimated from area under the curve)'}, 'FontSize', 12);
saveas(gca, fullfile(folder_save,'Filopodia_index_(area_under_the_curve)_left_&_right'), 'fig');

% 2D plot of Ratio of filopodia index and Ratio of distances between left/right TCJ against time
figure(14);
% getting the ratio of left to right filopodia indices. 
% Since some of the left and right filopodia index values are zero, some arrangement needs to be done. In general, if right (denominator) is zero or a very small value,
% then the ratio is infinity or very high, respectively. And if left (numerator) is zero then the ratio is zero. Therefore here, we make an arrangement such that if the right (denominator) is 
% zero or smaller than "1", then it is replaced with "1". In this case, instead of infinity or a very high value, the value of numerator will be the value of ratio. 
% On the other hand, if the left (numerator) is zero, we leave it as zero so that the value of ratio becomes zero.
area_under_the_curve_right(area_under_the_curve_right < 1) = 1; % changing all values below "1" to "one"
ratio_filopodia_index = area_under_the_curve_left ./ area_under_the_curve_right;
yyaxis left
plot(Time, ratio_filopodia_index, 'bo-');
yyaxis right
plot(Time, ratio_TCJ_distance, 'ro-');
xlabel('Time [min]', 'FontSize', 18);
yyaxis left; ylabel('Filopodia index ratio (left/right) [a.u]', 'FontSize', 18);
yyaxis right; ylabel('Ratio of left/right TCJ distances [a.u]) ', 'FontSize', 18);
title('Filopodia index ratio and left/right TCJ distance ratio against time', 'FontSize', 12);
saveas(gca, fullfile(folder_save,'Filopodia_index_ratio_&_TCJ_distance_ratio_vs_time'), 'fig');

% 2D plot of filopodia index against the ratio of distances between the cell tip and left/right TCJs
figure(15);
plot(celltip_top_distance, ratio_filopodia_index, 'ko');
xlabel('Distance (cell tip - line connecting TCJs) [\mum]', 'FontSize', 18);  ylabel('Filopodia index ratio (left/right) [a.u]', 'FontSize', 18);
set(gca, 'XDir', 'reverse');
title('Filopodia index ratio against the distance to top (line connecting the TCJs)', 'FontSize', 12);
saveas(gca, fullfile(folder_save,'Filopodia_index_ratio_vs_distance_to_top'), 'fig');

% saving full workspace
save('all_TCJ_contents');


%---------------------------------------------------------%---------------------------------------------------------------------%---------------------------------------------------------------------


