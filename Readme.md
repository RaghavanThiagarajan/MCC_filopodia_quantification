# **Filopodia analysis pipeline**

This analysis pipeline is part of the paper _“Multiciliated cells use filopodia to probe tissue mechanics during epithelial integration in vivo”_ by Ventura & Amiri et al., XXX., XXX.

The pipeline is used to quantify the filopodia interaction with tricellular junctions (TCJ) during the radial intercalation of multiciliated cells (MCC) in _Xenopus_ embryonic epithelia. Below, the steps for practical execution of the pipeline are outliined. For the related context, see the ***Materials and Methods*** section of the paper mentioned above.

***

### **The pipeline is executed in 3 steps:**
1. Manual annotation of the leading cell edge and manual tracking of TCJs.
2. Extracting the region of interests (ROIs) and intensities of cortex and filopodia.
3. Plotting.


### **Softwares required:**
* [Fiji](https://imagej.net/software/fiji/downloads)
* MATLAB (preferrably 2017b) (https://se.mathworks.com/products/get-matlab.html?s_tid=gn_getml)
***

## **Analysis steps:**
### **1. Manual annotation of the leading cell edge and manual tracking of TCJs**
In this step, we get the ROIs of the cell boundary and the track coordinates of the TCJs.
* [Reslice](https://imagej.net/imaging/z-functions#stack-reslice) the 3D movies of intercalating cells in X-Z direction and manually annotate the cell boundary or the leading edge of the cell using "freehand tool" in Fiji. All annotations should be registered as ROI (by pressing "t" after every annotation). Then save all the annotations of the full time stack as ***Roi-contour.zip*** using the [ROI manager](https://imagej.nih.gov/ij/docs/guide/146-30.html#sub:ROI-Manager...) in Fiji. 

* The two tricellular junctions (TCJs) close to the Multiciliated cells (MCC) (see Fig 1b and 1d in the paper) should be tracked using the [manual tracking](https://imagej.nih.gov/ij/plugins/track/track.html) plugin in Fiji. Save the resulting track files as ***TCJ_left.csv*** and ***TCJ_right.csv***, respectively. 

* At this point, download the ***sample_data*** folder into a new folder ***test_filopodia_analysis*** in your local drive. Inside the ***sample_data*** folder, you will find: 
    1. a resliced movie (***C1.tif***)
    2. a ROI list consisting the annotations of the leading cell edge (***Roi-contour.zip***)
    3. the track files of TCJs (***TCJ_left.csv*** & ***TCJ_right.csv***)

### **2. Extracting the ROIs and intensities of cortex and filopodia**
In this step, we will use the ROI of leading cell edge to obtain the ROIs of cell cortex and filopodia, and their intensities.

* Download the Fiji script `ROI_shiftedROI_meanintensity_extraction.ijm`” to the ***test_filopodia_analysis*** folder.

* Drag and drop the script file into Fiji.

* Open your resliced data or ***C1.tif***

* Check if the data (***C1.tif***) is calibrated for pixel size and time interval. If not, select the data, and using “Properties” option under “Image” in Fiji menu, input all the parameters. Save and close the file.

* Open the saved file in Fiji. 

* Run the script.

* A dialog will ask for the following values in pixels:
    * No. of pixels to be shifted for filopodia
    * No. of pixels to be shifted for cortex
    * Line thickness for filopodia
    * Line thickness for cortex.

    The first two values define the location of filopodia and cortex w.r.to the leading cell edge. The next two values define their spatial spread in the image i.e. how thick a line should be to cover the entire width of the cortex or filopodia. These values should be manually estimated using Fiji. 
* Input these values and click “OK”.

* A dialog will ask for area selection to estimate the background intensity. Select the "rectangle or oval tool" in Fiji and select / draw a small area in the image where the actual signal (filopodia / cortex / cell boundary) is not present.
* Click “OK”.

When the Fiji script is finished running, it will generate a folder called ***data***, inside which you will find: 
* ***Image_dimensions_&_input_parameters.csv*** file which stores the image property values and the shift & thickness values that you had entered.

* Two ROI files, one each for cortex and filopodia. Check these ROIs (drag and drop into Fiji) to see if they correspond to the area of filopodia and cortex of the cell. If not, Run the script again with modified values of “shift” and “thickness” for filopodia and cortex. 
* Two folders, one each for Background and cell contour that contain the x-y coordinates and intensity values. If you navigate through the following folder: ***ROI_Coordinates_Cell_contour / Actual_&_shifted_ROI***, you will find the data sheets for individual timepoints. On opening one of the .csv files, you will find x-y coordinates (xpoints & ypoints) for the original contour that was manually drawn, and the y coordinates and background subtracted intensity values of shifted contours for filopodia (top_shift) and cortex (down-shift).

### **3. Plotting**
In this step, we will use the data generated from the previous step to make the plots.
* Download the `TCJ_Cell_contour_plots.m` file into the ***test_filopodia_analysis*** folder. 
* Open the script in MATLAB.
* Run the script.

When the script is finished running, you will find a ***all_TCJ_contents.mat*** file and a folder called ***plots***. The .mat file is a collection of all the MATLAB workspace variables. The ***plots*** folder contains various plots of filopodia, TCJs and cortex, in .fig (MATLAB figure) format. Before plotting, the filopodia intensities were normalised by the cortex intensity. The plot ***Shift_in_cell_contour_color_coded_filopodia_meanintensity_cortex_normalised.fig*** is similar to what you will see in Fig 1g, Supplementary figure 1d and Supplementary movie 2 in the paper. The name and title of the plots are self-explanatory. You can save the .fig files into different formats like .tif, .jpg etc.
***