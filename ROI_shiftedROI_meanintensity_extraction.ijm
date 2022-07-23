/* Raghavan Thiagarajan, DanStem, Copenhagen, 30th September 2018 - Version-02, 06th July 2020

This code extracts the coordinates from the saved ROI's (of cell-surface contour).
Also it finds the minima of y coordinate and the corresponding x coordinate for each ROI, which will be used for plotting distance 
with other points obtained elsewhere from tracking (of Tri-cellular junctions). And it shifts the ROI of the cell contour a little above
to measure the mean intensity of the filopodia. */

/* In this version-02, following major modifications were made: 
 * In addition to shifting the contour up for filopodia measurement, the contour is now also shifted down to measure the cortex intensity; This intensity
 * will be used for normalising the filopodia intensity. The original contour length is measured. The ROI's of the newly generated contours for filopodia 
 * and cortex are saved to allow future checking.
 * 
 * Following mistakes from version-01 were corrected:
 * Median used to find the contour tip was changed from (nResults/2)-1 to (nResults/2)+1
 * The contour chosen for mid point finding was corrected to be the oringinal contour
 * Roi.getcoordinates is replaced by Roi.getcontainedpoints
 * line thickness is carefully adjusted to the appropriate value
 * 
 * In addition to this, following changes were made to improve the efficiency of the code:
 * Array display is prevented
 * Data is now genereated as tables and concatenated instead of generating them as arrays
 * data is now stored as .csv files instead of .txt files 
 */

//----------------------------------------------------------------------------------------------------------------------------------------------//
// Fetching and saving inputs

id = getImageID();

// Getting values for contour shifting and line thickness separately for filopodia and cortex
Dialog.create("Pixel value to be shifted");
Dialog.addNumber("No. of pixels to be shifted for filopodia [in pixels]:", 0000);
Dialog.addNumber("No. of pixels to be shifted for cortex [in pixels]:", 0000);
Dialog.addNumber("Line thickness for filopoda [in pixels]:", 0000);
Dialog.addNumber("Line thickness for cortex [in pixels]:", 0000);
Dialog.show;
shift_for_filopodia = Dialog.getNumber();
shift_for_cortex = Dialog.getNumber();
line_thickness_filopodia = Dialog.getNumber();
line_thickness_cortex = Dialog.getNumber();

// Directory is defined for saving the files
Path = getInfo("image.directory"); 
Main_dir = Path+"/data/"; File.makeDirectory(Main_dir);
Dir_Roi_Background_intensity=Main_dir+"/ROI_Coordinates_Background/"; File.makeDirectory(Dir_Roi_Background_intensity); 
Dir_Roi_Cell_contour=Main_dir+"/ROI_Coordinates_Cell_contour/"; File.makeDirectory(Dir_Roi_Cell_contour);
Shifted_ROI=Dir_Roi_Cell_contour+"/Actual_&_shifted_ROI/"; File.makeDirectory(Shifted_ROI); 

// Image dimensions are obtained and saved along with input parameters
getDisplayedArea(x, y, width, height);
getPixelSize(unit, pixelWidth, pixelHeight);
Time_interval = Stack.getFrameInterval();
run("Set Measurements...", "area mean perimeter integrated redirect=None decimal=3");
setOption("ShowRowIndexes", true); // setting the table index column - this will be the first column in all tables that indicates the table index and is useful while inspecting the table
Table.set("Image_width", 0, width); Table.set("Image_height", 0, height); 
Table.set("Pixel_width", 0, pixelWidth); Table.set("Pixel_height", 0, pixelHeight); 
Table.set("Time_interval", 0, Time_interval); Table.set("shift_for_filopodia", 0, shift_for_filopodia); 
Table.set("shift_for_cortex", 0, shift_for_cortex); Table.set("line_thickness_filopodia", 0, line_thickness_filopodia); 
Table.set("line_thickness_cortex", 0, line_thickness_cortex); 
Table.save(Main_dir + "Image_dimensions_&_input_parameters.csv"); 
setOption("ShowRowIndexes", false);

// Setting flexible array lengths
setOption("ExpandableArrays", true); 
xpoints_min = newArray; ypoints_min = newArray; contour_length = newArray; 

//----------------------------------------------------------------------------------------------------------------------------------------------//
// Background intensity measurements

if (isOpen("ROI Manager")) {
     selectWindow("ROI Manager");
     run("Close");
}

// Selecting an area and measuring the backgorund mean intensity
waitForUser("Select an area outside the ROI for calculating \nmean backgorund intensity and then click \"OK\"");
roiManager("add"); 
Background_meanintensity=newArray;
setBatchMode(true);
for (i=1; i<=nSlices; i++) {
	selectImage(id);
	roiManager("select", 0);
	setSlice(i);
	run("Measure");
	l = getResult("Mean", i-1); 
	Background_meanintensity[i-1] = l;    	
}

selectWindow("ROI Manager"); roiManager("save", Dir_Roi_Background_intensity +"ROI_Backgroundintensity" + ".zip"); run("Close");
selectWindow("Results"); //saveAs("results", Dir_Roi_Background_intensity +"Measure_background_results" + ".csv"); 
close("Results");
setBatchMode(false);

//----------------------------------------------------------------------------------------------------------------------------------------------//
// Roi shifting, intensity measurements of shifted roi's and finding the midpoint of original contour for 
// distance measurement (from tri-cellular junctions)

if (isOpen("ROI Manager")) {
     selectWindow("ROI Manager");
     run("Close");
}

// Shifting roi towards top (for filopodia) and down (for cortex) and measuring the respective mean intensities 
roiManager("open", Path +"ROI-contour" + ".zip"); // opening the roi of the contour that was manually drawn
ROI_count = roiManager("count");
setBatchMode(true);
for (j=1; j<=ROI_count; j++) {
	
	// selecting the corresponding roi inorder to measure the contour length of the original contour
	selectImage(id);
	roiManager("select", 0); setSlice(j);  
	run("Measure"); // this contour length is used to normalise the "area under the curve" in the matlab script
	contour_length [j-1] = getResult("Length", 0); // the index is set to zero since this results window is closed immediately
	close("Results");

	// selecting the corresponding roi and making sure the line width is "1" to use the roi.getcontainedpoints function
	selectImage(id);
	roiManager("select", 0);
	setSlice(j); run("Properties... ", "  width=1"); 

	//--------------------------------------------------------------------//

	// Shifting the ROI of cell contour up and down and getting the intensity of filopodia and cortex respectively
	shift_ROI_and_get_meanintensity(); // this function is called

	//--------------------------------------------------------------------//

	// Minima of the Y coordinate is obtained since in this image, this Y coordinate will define the 
	// topmost point of cell surface which will be used for measuring distance with Tri-cellular junction points
	// getting the pixel coordinates of the original roi
	selectImage(id);
	roiManager("select", 0); // selecting the original contour that was manually drawn.
	setSlice(j); run("Line Width...", "line=1"); run("Properties... ", "  width=1"); 
	Roi.getContainedPoints(xpoints, ypoints);
	ypoints_sort = Array.sort(ypoints); 
	ypoints_minn = ypoints_sort[0]; 
	Roi.getContainedPoints(xpoints, ypoints); run("Select None"); 
	
	// This "for loop" ensures that if there is a set of Y minima (i.e more than one Y minima), then the median of this set is taken as the Y minima.
	for (r=0; r<=ypoints.length-1; r++) { 
		if (ypoints[r]==ypoints_minn) {
			setResult("Index_of_min_values", nResults, r); 
		}
	}
	updateResults;
	 
	 // getting the median pixel of the original contour
	if ((nResults == 1) || (nResults == 2)) {
		ypoints_minindex=getResult("Index_of_min_values", 0); 
	}
	else {
		ypoints_minindex=getResult("Index_of_min_values", ((nResults/2)+1)); 		
	}	
	selectWindow("Results"); /*saveAs("results", Path +"ROI_results"+j + ".txt"); */ close("Results");
	
	//Corresponding X coordinate of the above obtained Y (minima) coordinate is fetched
	ypoints_min[j-1] = ypoints[ypoints_minindex];
	xpoints_min[j-1] = xpoints[ypoints_minindex];
	
	roiManager("select", 0); // selecting the current roi
	roiManager("delete"); // deleting the current roi inorder to automatically move to the next roi
}

// saving all the shifted roi's - these roi's will be sorted in the steps below
selectWindow("ROI Manager"); roiManager("save", Main_dir +"ROI_shifted_contour" + ".zip"); run("Close"); 
// saving the x and y coordinates of the median pixel of the original contour i.e.  minimal point (top most point / tip of the contour) from actual roi coordinates 
setOption("ShowRowIndexes", true);
Table.setColumn("xpoints_min", xpoints_min); Table.setColumn("ypoints_min", ypoints_min); 
Table.save(Dir_Roi_Cell_contour+"Actual_ROI_x_&_y_min"+".csv"); run("Close"); 
// saving the contour length of the original contour - this contour length is used to normalise the "area under the curve" in the matlab script
Table.setColumn("contour_length", contour_length); 
Table.save(Dir_Roi_Cell_contour+"contour_length"+".csv"); run("Close"); 
// saving the backgrouund intensities and the corresponding roi used for subtraction 
Table.setColumn("Background_meanintensity", Background_meanintensity);
Table.save(Dir_Roi_Background_intensity + "Background_meanintensity"+".csv"); run("Close");
setOption("ShowRowIndexes", false);
setBatchMode(false);

//----------------------------------------------------------------------------------------------------------------------------------------------//
// Opening the saved "ROI_shifted_contour" and splitting the roi for filopodia and roi for cortex.

if (isOpen("ROI Manager")) {
     selectWindow("ROI Manager");
     run("Close");
}

// Total number of roi's in "ROI_shifted_contour".
ROI_count *= 2; // The number of roi's in the "ROI_shifted_contour" will be twice the number of roi's in the original "ROI-contour"

// Based on the total number of roi's in the "ROI_shifted_contour", two arrays "odd" and "even" are created to split the roi's into 
// roi for filopodia and roi for cortex.
even = newArray; odd = newArray;
for (roi = 0; roi <= (ROI_count / 2)-1; roi++) {
	even[roi] = roi * 2;
	odd[roi] = (roi * 2) + 1;
}

// Based on the "odd" / "even" arrays, the "ROI_shifted_contour" is split into "roi for cortex" and "roi for filopodia" and
// saved as separate roi files.
for (i = 1; i <= 2; i++) {
	roiManager("open", Main_dir +"ROI_shifted_contour" + ".zip");
	if (i == 1) {
		roiset = even;		
		suffix = "cortex";
	}
	else {
		roiset = odd;
		suffix = "filopodia";
	}
	roiManager("select", roiset);
	roiManager("delete");
	roiManager("save", Main_dir +"ROI_shifted_contour_"+ suffix + ".zip"); run("Close");
}
run("Select None"); 
File.delete(Main_dir + "ROI_shifted_contour" + ".zip");	
selectWindow("Log"); run("Close");
close("Results");

//---------------------------------END------------------------------------------------------END-----------------------------------------END-------------------------------------------//





//----------------------------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------------//
// function for shifting the ROI (in both directions) and fetching the meanintensity of the new ROI contour.

function shift_ROI_and_get_meanintensity() {
	y_top_shift = newArray; top_shifted_meanintensity = newArray; // top_shift corresponds to the shift in contour for filopodia
	y_down_shift = newArray; down_shifted_meanintensity = newArray; // down_shift corresponds to the shift in contour for cortex
	Roi.getContainedPoints(xpoints, ypoints); // refetching the pixels of the original contour
	// making sure the same roi selection of original contour is not confused with the roi's to be fetched in future during roi shifting
	run("Select None"); 
	
	// creating arrays with the shift values - by which the original contour has to be shifted for getting filopodia and cortex intensities
	for (jj=1; jj<=xpoints.length; jj++) {
		y_top_shift[jj-1] = ypoints[jj-1] - shift_for_filopodia; // Distance (in pixels) by which the contour is moved top - for filopodia measurement
		y_down_shift[jj-1] = ypoints[jj-1] + shift_for_cortex; // Distance (in pixels) by which the contour is moved down - for cortex measurement
	}

	// shifting the roi twice (for filopodia anc cortex) and measuring the intensities
	setBatchMode(true);
	// shift value is 2 since the roi is shifted twice
	for (shift_no = 1; shift_no <= 2; shift_no++) { 
		// first shift is for filopodia
		if (shift_no == 1) {
			yshift = y_top_shift;
			line_thickness = line_thickness_filopodia;
		}
		// second shift is for cortex
		else {
			yshift = y_down_shift;
			line_thickness = line_thickness_cortex;
		}
		// creating selection based on the shifted values
		selectImage(id);
		makeSelection("freeline", xpoints, yshift);
		// getting meanintensity from shifted ROI
		// Increasing the thickness of the contour to cover the filopodia and cortex and subsequently measuring the mean intensity
		run("Line Width...", "line="+line_thickness); 
		run("Properties... ", "  width="+line_thickness); // making sure the line thickness is assigned
		run("Plot Profile"); 
		Plot.getValues(x, y); 
		run("Close");
		// adding the newly created shifted contour to the roi manager - for it to be saved later
		roiManager("add");
		selectImage(id); 
		// making sure the line width is assigned back to "1" in order to use the roi.getcontainedpoints function later
		run("Line Width...", "line=1"); run("Properties... ", "  width=1"); run("Select None"); 

		// subtracting the background mean intensity from the intensity measured for the shifted contours
		Bkgrnd = Background_meanintensity[j-1];
		for (b=1; b<=x.length; b++) {
			if (shift_no == 1) {
				top_shifted_meanintensity[b-1] = y[b-1] - Bkgrnd; // Background is subtracted from the mean intensity of the contour
			}
			else {
				down_shifted_meanintensity[b-1] = y[b-1] - Bkgrnd; // Background is subtracted from the mean intensity of the contour
			}
		}
	}

	setOption("ShowRowIndexes", true);
	// saving the shifted values of the original contour for both filopodia and cortex
	Table.setColumn("xpoints", xpoints); Table.setColumn("ypoints", ypoints); // columns for actual roi coordinates
	// columns for top shifted roi coordinates and mean intensity for filopodia
	Table.setColumn("y_top_shift", y_top_shift); Table.setColumn("top_shifted_meanintensity", top_shifted_meanintensity); 
	// columns for down shifted roi coordinates and mean intensity for cortex
	Table.setColumn("y_down_shift", y_down_shift); 	Table.setColumn("down_shifted_meanintensity", down_shifted_meanintensity);
	// saving all these columns as a table
	Table.save(Shifted_ROI+"Actual_&_shifted_ROIcoord_with_Meanint_"+j+".csv"); run("Close"); 
	setOption("ShowRowIndexes", true);
	
	
	setBatchMode(false);
}
