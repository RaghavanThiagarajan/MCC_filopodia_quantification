// Raghavan Thiagarajan, DanStem, Copenhagen, 30th September 2018.
// This code extracts the coordinates from the saved ROI's (of cell-surface contour).
// Also it finds the minima of y coordinate and the corresponding x coordinate for each ROI, which will be used for plotting distance 
// with other points obtained elsewhere from tracking (of Tri-cellular junctions). And it shifts the ROI of the cell contour a little above
// to measure the mean intensity of the filopodia.

// Directory is defined for saving the files
Path = getInfo("image.directory"); 
Dir_Roi_Background_intensity=Path+"/ROI_Coordinates_Background/"; File.makeDirectory(Dir_Roi_Background_intensity); 
Dir_Roi_Cell_contour=Path+"/ROI_Coordinates_Cell_contour/"; File.makeDirectory(Dir_Roi_Cell_contour);
Shifted_ROI=Dir_Roi_Cell_contour+"/Shifted_ROI/"; File.makeDirectory(Shifted_ROI); 
ROI_x_y=Dir_Roi_Cell_contour+"/Actual_ROI_x_&_y_values/"; File.makeDirectory(ROI_x_y); 

// Image dimensions are obtained
getDisplayedArea(x, y, width, height);
getPixelSize(unit, pixelWidth, pixelHeight);
Time_interval = Stack.getFrameInterval();
run("Set Measurements...", "area mean perimeter integrated redirect=None decimal=3");
print("Image"+"\t"+"Image"+"\t"+"Pixel"+"\t"+"Pixel"+"\t"+"Time"+"\n"+"width"+"\t"+"height"+"\t"+"width"+"\t"+"height"+"\t"+"Interval"+"\n"+ width+"\t"+ height+"\t"+pixelWidth+"\t"+pixelHeight+"\t"+Time_interval);
selectWindow("Log"); saveAs("Text", Path +"Image_dimensions" + ".txt"); run("Close"); 

// Setting flexible array lengths
setOption("ExpandableArrays", true); xpoints_min = newArray; ypoints_min = newArray; 

// Selecting an area and measuring the backgorund mean intensity
if (isOpen("ROI Manager")) {
     selectWindow("ROI Manager");
     run("Close");
}
waitForUser("Select an area outside the ROI for calculating \nmean backgorund intensity and then click \"OK\"");
roiManager("add"); 
Background_meanintensity=newArray;
for (i=1; i<=nSlices; i++) {
	roiManager("select", 0);
	setSlice(i);
	run("Measure");
	l=getResult("Mean", i-1);
    Background_meanintensity[i-1] = l;
    if (i==nSlices) {
    	Array.show(Background_meanintensity);
    	}
}
selectWindow("ROI Manager"); roiManager("save", Dir_Roi_Background_intensity +"ROI_Backgroundintensity" + ".zip"); run("Close");
selectWindow("Results"); saveAs("results", Dir_Roi_Background_intensity +"Measure_background_results" + ".txt"); close("Results");

// ROI manager is used to fetch and process the ROI's
if (isOpen("ROI Manager")) {
     selectWindow("ROI Manager");
     run("Close");
}
roiManager("open", Path +"ROI-contour" + ".zip");
ROI_count = roiManager("count");
for (j=1; j<=ROI_count; j++) {
		roiManager("select", j-1);
		setSlice(j);
		yshift = newArray; Shifted_meanintensity = newArray;

		// Shifting the ROI of cell contour and getting the intensity of filopodia.
		Roi.getCoordinates(xpoints, ypoints); 
		for (jj=1; jj<=xpoints.length; jj++) {
			yshift[jj-1] = ypoints[jj-1] - 10; // Distance (in pixels) by which the contour is moved
			if (jj==xpoints.length) {
				Array.show(xpoints); Array.show(yshift);
			}
		}
		makeSelection("freeline", xpoints, yshift);
		run("Line Width...", "line=8"); // Increasing the thickness of the contour to cover the filopodia and subsequently measuring the mean intensity
		run("Plot Profile");
		Plot.getValues(x, y);
		run("Close");
		selectWindow("Background_meanintensity");
		Bkgrnd = Background_meanintensity[j-1];
		for (b=1; b<=x.length; b++) {
			Shifted_meanintensity[b-1] = y[b-1] - Bkgrnd; // Background is subtracted from the mean intensity of the contour
			if (b==x.length) {
				Table.showArrays("Shifted_ROI_Coordinates_&_Mean_intensity_"+j, xpoints, yshift, Shifted_meanintensity);			
			}
		}
		selectWindow("Shifted_ROI_Coordinates_&_Mean_intensity_"+j); Table.save(Shifted_ROI+"Shifted_ROIcoordinates_&_Meanint_"+j+".txt"); run("Close");
		selectWindow("xpoints"); run("Close"); selectWindow("yshift"); run("Close");
		
		// Minima of the Y coordinate is obtained since in this image, this Y coordinate will define the 
		// topmost point of cell surface which will be used for measuring distance with Tri-cellular junction points.
		Roi.getCoordinates(xpoints, ypoints); 
		ypoints_sort = Array.sort(ypoints); 
		ypoints_minn = ypoints_sort[0]; 
		Roi.getCoordinates(xpoints, ypoints); 

		// This "for loop" ensures that if there is a set of Y minima (i.e more than one Y minima), then the median of this set is taken as the Y minima.
		for (r=0; r<=ypoints.length-1; r++) { 
			if (ypoints[r]==ypoints_minn) {
				setResult("Index_of_min_values", nResults, r); 
			}
		}
		updateResults;
		ypoints_minindex=getResult("Index_of_min_values", ((nResults/2)-1)); 
		selectWindow("Results"); /*saveAs("results", Path +"ROI_results" + ".txt");*/ close("Results");

		//Corresponding X coordinate of the above obtained Y (minima) coordinate is fetched
		ypoints_min[j-1] = ypoints[ypoints_minindex];
		xpoints_min[j-1] = xpoints[ypoints_minindex];
		Table.showArrays("ROI_x_&_y_min", xpoints_min, ypoints_min);
		Table.showArrays("ROI_x_&_y_", xpoints, ypoints);			
		selectWindow("ROI_x_&_y_"); Table.save(ROI_x_y+"ROI_x_&_y_" +j+ ".txt"); run("Close");
}
selectWindow("ROI_x_&_y_min"); Table.save(Dir_Roi_Cell_contour+"Actual_ROI_x_&_y_min"+".txt"); run("Close");
selectWindow("Background_meanintensity"); saveAs("text", Dir_Roi_Background_intensity +"Background_meanintensity" + ".txt"); run("Close");
selectWindow("ROI Manager"); run("Close"); 
