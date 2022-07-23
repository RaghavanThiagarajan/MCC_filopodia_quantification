% Raghavan Thiagarajan, DanStem, Copenhagen, 30th September 2018.
% This script obtains different kinds of data and plots them in 2D and 3D.
% The data obtained are: (1) Coordinates from the track file of Tri-cellular junctions; (2) The
% ROI coordinates of cell surface contour; (3) Image dimensions and the
% corresponding (4) Time points (5) Shifted ROI coordinates and their mean intensity.

clear all;
clc

% Loading Image dimensions 
Imgdim = 'Image_dimensions.txt';
Imagedim = importdata(Imgdim);
Image_dime = Imagedim.data;

% Getting pixel values to convert values into micrometers
xpixval = Image_dime(1,3); % In microns
ypixval = Image_dime(1,4); % In microns
Timevall = Image_dime(1,5); % In min
Timeval = Timevall / 60; % Converted to hour

% Getting Image dimensions - to be used for normalization
Img_height = Image_dime(1,2);
Image_height = Img_height * ypixval;

% Loading the tracked data from the Tri-cellular junction on the left
TCJLeft = 'TCJ_left.csv'; 
TCJleft=importdata(TCJLeft); 
TCJ_left = TCJleft.data; 
TCJleftx = TCJ_left(:,3); TCJlefty = TCJ_left(:,4); 
TCJ_leftx = TCJleftx.*xpixval ; TCJ_lefty = TCJlefty.*ypixval; % x and y coordinates
TCJlefty_normalised = Image_height -  TCJ_lefty; % Normalized w.r.to the Image height for proper representation since the XY quadrants are arranged differently in ImageJ
TCJlefty_normalised_guideforeyes = TCJlefty_normalised / max(TCJlefty_normalised);

sz=size(TCJ_leftx); % This gives the total number of tracked points
sze = sz(1,1);

% Loading the tracked data from the Tri-cellular junction on the right
TCJRight = 'TCJ_right.csv';
TCJright=importdata(TCJRight);
TCJ_right = TCJright.data;
TCJrightx = TCJ_right(:,3); TCJrighty = TCJ_right(:,4); 
TCJ_rightx = TCJrightx.*xpixval; TCJ_righty = TCJrighty.*ypixval; % x and y coordinates
TCJrighty_normalised = Image_height -  TCJ_righty; % Normalized w.r.to the Image height for proper representation since the XY quadrants are arranged differently in ImageJ
TCJrighty_normalised_guideforeyes = TCJrighty_normalised / max(TCJrighty_normalised);

Tim = TCJ_left (:,2); Time0 = Tim.*Timeval; T=Time0.*ones(sz); % Obtaining time points
Tim1 = TCJ_right(:,2); Time1 = Tim1.*Timeval; T1=Time1.*ones(size(TCJ_rightx));

% Defining color range
mymap = [0 1 1
         0 1 1
         0 1 1
         1 0 1
         1 0 1];

for i = 1:3
    if i==3
        figure(i);
        plot3(TCJ_leftx, T, TCJlefty_normalised_guideforeyes, '-ob'); grid on; % 3D plot of Tri-cellular junction on the left
        hold on
        plot3(TCJ_rightx, T1, TCJrighty_normalised_guideforeyes, '-ob'); % 3D plot of Tri-cellular junction on the right
        hold on
    else
        figure(i);
        plot3(TCJ_leftx, T, TCJlefty_normalised, '-ob'); grid on; % 3D plot of Tri-cellular junction on the left
        hold on
        plot3(TCJ_rightx, T1, TCJrighty_normalised, '-ob'); % 3D plot of Tri-cellular junction on the right
        hold on
    end 
end

% Loading all coordinates for each ROI and obtaining separately the minima of Y coordinate and the
% corresponding X coordinate for each ROI.
ROI_xy = './ROI_Coordinates_Cell_contour/Actual_ROI_x_&_y_values/ROI_x_&_y_';
ShiftedROI_xy = './ROI_Coordinates_Cell_contour/Shifted_ROI/Shifted_ROIcoordinates_&_Meanint_';

ROI_xyminn = './ROI_Coordinates_Cell_contour/Actual_ROI_x_&_y_min.txt'; 
ROIxyminn = importdata(ROI_xyminn); 
ROI_xymin = ROIxyminn.data; 

Distanceleft = zeros(1, sze); Distanceright = zeros(1, sze); Time = zeros(1,sze); DistanceTCJ = zeros(1,sze);
for i = 1:sze
    ROIxx = ROI_xymin(i,1); ROIyy = ROI_xymin(i,2);
    x1 = ROIxx * xpixval; y1 = ROIyy * ypixval;    
    
    leftxx = TCJ_left(i,3); leftyy = TCJ_left(i,4);
    x2 = leftxx * xpixval; y2 = leftyy * ypixval;    
    
    rightxx = TCJ_right(i,3); rightyy = TCJ_right(i,4);
    x3 = rightxx * xpixval; y3 = rightyy * ypixval;
    
    % Getting the distance between the tri cellular junction points and the
    % tip of each ROI (cell contour)    
    Distanceleft(i) = sqrt(((x2-x1)^2)+((y2-y1)^2)); 
    Distanceright(i) = sqrt(((x3-x1)^2)+((y3-y1)^2)); 
    
    % Getting the distance between the tricellular junctions
    DistanceTCJ(i) = sqrt(((x3-x2)^2)+((y3-y2)^2));
    
    % Getting the time
    Time(i) = i*Timeval;
    
    % Getting the ROI coordinates and plotting the ROI contour on the
    % already plotted Tricellular junction tracks    
    ROIXY = strcat(ROI_xy, sprintf('%01d',i), '.txt'); 
    ROI_XY=importdata(ROIXY);
    ROI_XYcrd = ROI_XY.data; 
    ROI_Xcrd = ROI_XYcrd(:,1); ROI_Ycrd = ROI_XYcrd(:,2);
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
    
    Tim2=i*Timeval; T2_1=Tim2.*ones(size(CellTip(:,1))); T2_2=Tim2.*ones(size(ROI_Xcoord));
        
    figure(1);
    plot3(ROI_Xcoord, T2_2, ROIYcoord_normalised, '-k','LineWidth', 0.8); % 3D plot of cell surface contour
    %view(gca, [10, 45, 45]);
    hold on
    plot3(CellTip(:,1), T2_1, CellTip(:,2), '-m','LineWidth', 0.8); % 3D plot of only the cell contour tip
    hold on
    
    
    % Getting the shifted ROI coordinates and plotting the shifted ROI contour on the
    % already plotted Tricellular junction tracks where the shifted ROI
    % contour is color coded for its mean intensity on the image
    ShiftedROIXY = strcat(ShiftedROI_xy, sprintf('%01d',i), '.txt'); 
    ShiftedROI_XY=importdata(ShiftedROIXY);
    ShiftedROI_XYcrd = ShiftedROI_XY.data; 
    ShiftedROI_Xcrd = ShiftedROI_XYcrd(:,1); ShiftedROI_Ycrd = ShiftedROI_XYcrd(:,2); 
    ShiftedROI_meanint = ShiftedROI_XYcrd(:,3); ShiftedROI_meanint = ShiftedROI_meanint / (max(ShiftedROI_meanint)); 
    ShiftedROI_Xcoord = ShiftedROI_Xcrd * xpixval; ShiftedROI_Ycoord = ShiftedROI_Ycrd * ypixval;
    ShiftedROIYcoord_normalised = Image_height - ShiftedROI_Ycoord;
    
    T3=Tim2.*ones(size(ShiftedROI_Xcoord));
    
    figure(2);
    plot3(ShiftedROI_Xcoord, T3, ShiftedROIYcoord_normalised); % 3D plot of cell surface contour where the contour is color coded for mean intensity
    surface([ShiftedROI_Xcoord(:), ShiftedROI_Xcoord(:)], [T3(:), T3(:)], [ShiftedROIYcoord_normalised(:), ShiftedROIYcoord_normalised(:)], [ShiftedROI_meanint(:), ShiftedROI_meanint(:)], 'EdgeColor','interp', 'FaceColor','none','LineWidth', 0.8);
    colormap(mymap);
    hold on
    
    figure(3);
    plot3(ShiftedROI_Xcoord, T3, ShiftedROI_meanint); % 3D plot of filopodia mean intensity
    surface([ShiftedROI_Xcoord(:), ShiftedROI_Xcoord(:)], [T3(:), T3(:)], [ShiftedROI_meanint(:), ShiftedROI_meanint(:)], [ShiftedROI_meanint(:), ShiftedROI_meanint(:)], 'EdgeColor','interp', 'FaceColor','none','LineWidth', 0.8);
    colormap(mymap);
    %view(gca, [10, 45, 45]);
    hold on
end
figure(1);
xlabel('Lateral distance [\mum]');  ylabel('Time [h]'); zlabel('Height [\mum]');  
saveas(gca, 'Biased_shift_in_cell_contour.fig'); 

figure(2);
xlabel('Lateral distance [\mum]');  ylabel('Time [h]'); zlabel('Height [\mum]');  
saveas(gca, 'Biased_shift_in_cell_contour_&_filopodia_meanintensity.fig'); 

figure(3);
xlabel('Lateral distance [\mum]');  ylabel('Time [h]'); zlabel('Mean intensity [a.u]');  
saveas(gca, 'Biased_shift_in_filopodia_meanintensity.fig'); 


% 2D plot of the distance between the left and right tri-cellular junctions and the tip of cell surface contour (or y minima and the corresponding x)
figure(4);
plot (Time, Distanceleft, '-ob'); 
hold on
plot (Time, Distanceright, '-or'); 
xlabel('Time [h]');  ylabel('Distance [\mum]'); legend
saveas(gca, 'Distance between TCJ and cell surface.fig'); 
legend({'TCJ left','TCJ right'},'FontSize',12,'TextColor','black')

% 2D plot of the distance between the left and right tri-cellular junctions
figure(5);
plot (Time, DistanceTCJ, '-ob'); 
xlabel('Time [h]');  ylabel('Distance [\mum]');
saveas(gca, 'Distance between TCJ.fig');



