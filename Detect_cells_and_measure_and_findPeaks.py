#@ File(label="Input directory", style="directory") inputDir
#@ File(label="Output directory", style="directory") outputDir
#@ String (label = "File extension", value=".tif") ext
#@ String (label = "File name contains", value = "crop") containString

# ImageJ Jython script to measure intensity and inclusions within cell area in a 2-channel Z-stack

# TODO:
#1. Create a function to Find Maxima
#2. Decide what data to collect from maxima (count, intensity of the maxima at center)
#3. Organize results in a table and save that in a SEPARATE table


# ---- IMPORTS

from ij import IJ, ImagePlus, ImageStack
from ij.plugin import ChannelSplitter
from ij.io import DirectoryChooser
import os
from ij.plugin import ImageCalculator
from ij.gui import PointRoi
from ij.plugin.frame import RoiManager
from net.imglib2.img.display.imagej import ImageJFunctions as IL
from net.imglib2.algorithm.dog import DogDetection
from net.imglib2.view import Views
from jarray import zeros
from ij.measure import ResultsTable

# --- FUNCTION DEFINITIONS
# --- This batch file processor includes 2 basic functions:
# 		1) "run": traverses the input directory
#		2) "process": executes the actual job of the script,
#			including calling additional functions if desired

def run(): # this is the function that walks through the directory
	""" Walk through the directory and identify files matching our criteria """

	# Setup
	IJ.run("Close All", "") # close all open windows
	IJ.run("Clear Results", "") # clear the Results window
	IJ.log("\\Clear") # clear the ImageJ log
	
	IJ.log("Beginning batch processing")
	
	# Convert user-selected directories into Python paths
	srcDir = inputDir.getAbsolutePath()
	tarDir = outputDir.getAbsolutePath()
	
	#Tables of measurements and peaks
	peaksTable = ResultsTable()
	
	# Loop through directories and sub-directories
	for root, directories, filenames in os.walk(srcDir):
		filenames.sort();
		for filename in filenames:
			IJ.log("Checking file " + filename)
			# Check for file extension
			# The string ext is defined in the script parameters
			if not filename.endswith(ext):
				continue
			# Check for file name pattern
			# The string containString is defined in the script parameters
			if containString not in filename:
				continue
			# Process the file
			process(srcDir, tarDir, root, filename, peaksTable)
			
			#Go to next row in peaksTable
			peaksTable.incrementCounter()
	
	#Show peaksTable
	peaksTable.show("Peaks")		
	
	# Save the accumulated results and peak tables in one file
	#ResultsTable.save(os.path.join(tarDir, "Results.csv"))
	peaksTable.save(os.path.join(tarDir, "Peaks.csv"))
	IJ.selectWindow("Results")
	IJ.saveAs("Results", os.path.join(tarDir, "Results.csv"))
	
	IJ.log("Done")
	
	#UN-COMMENT LINES IF YOU WANT THE TABLES TO CLOSE:
	IJ.run("Close")
	IJ.selectWindow("Peaks")
	#IJ.run("Close")
	
	#UN-COMMENT LINES IF YOU WANT ROI MANAGER TO CLOSE:
	IJ.selectWindow("ROI Manager")
	IJ.run("Close")
	
	#UN-COMMENT LINES IF YOU WANT THE LOG TO CLOSE:
	IJ.selectWindow("Log")
	IJ.run("Close")

def process(srcDir, tarDir, root, filename, peaksTable):
	""" Run processing functions on a file and save output """
	
	IJ.log("Processing " + filename)
	
	# Get the image path
	imagePath = os.path.join(srcDir, filename)
		
	# Open the image as an ImagePlus
	imp = IJ.openImage(imagePath)
	
	# Get the original image name
	origName = imp.getTitle()
	
	#get the right channel from split
	imp = IJ.openImage(imagePath)
	channels = ChannelSplitter.split(imp)
	imp = channels[0]
	imp.setTitle(origName)
	
	# Generate a masked image where areas outside the cell volume are 0
	maskedImg = threshold_stack(imp, "Huang")
	
	# Measure the stack
	measure_stack(maskedImg)
	
	#Find peaks in masked image
	findPeaks(imp, peaksTable, srcDir, tarDir)
	
	# Save output: tarDir defined in script parameters
	maskedFilename = filename + "_Masked.tif" 
	
	#resultsFilename = filename + "_Results.csv"
	IJ.log("Saving output for " + filename)
	IJ.saveAs(maskedImg, "Tiff", os.path.join(tarDir, maskedFilename)) 
	#IJ.saveAs("Results", os.path.join(tarDir, resultsFilename))
	imp.close()
	return


def threshold_stack(imp, method):
	""" Apply the indicated threshold method to an ImagePlus, clean up, and return the masked image """
	
	# Duplicate the image
	origName = imp.getTitle()
	mask = imp.duplicate()
	mask.setTitle("Thresholded")
	
	# Threshold the stack
	IJ.run(mask, "Auto Threshold", "method="+method+" white stack use_stack_histogram")
	
	# Clean up the cell border using binary operations: close, remove outliers, fill
	IJ.run(mask, "Close-", "stack")
	IJ.run(mask, "Remove Outliers...", "radius=2 threshold=50 which=Bright stack")
	IJ.run(mask, "Fill Holes", "stack")
	
	# Divide the mask by 255 so that values are either 0 or 1
	IJ.run(mask, "Divide...", "value=255 stack")
	
	# mask the original image by multiplying it by the mask
	ic = ImageCalculator()
	maskedImg = ic.run("Multiply create stack", imp, mask)		
	maskedImg.setTitle(origName + "_Masked")
	
	return maskedImg


def measure_stack(imp):

	""" Measure intensity within cell volume defined by a mask """
	
	# set up measurements to collect
	IJ.run("Set Measurements...", "area mean integrated limit display decimal=2")
	
	# set the threshold so pixels=0 are not included
	
	# TODO: Figure out how to determine the max value for the image 
	# and use that as the max threshold instead of 65535
	# maxValue = 
	IJ.setRawThreshold(imp, 1, 65535)
	
	# measure the stack
	IJ.run(imp, "Measure Stack...", "");
	
	return

def findPeaks(imp, peaksTable, srcDir, tarDir):
		# Access its pixel data from an ImgLib2 data structure: a RandomAccessibleInterval  
	ipRAI = IL.wrapReal(imp)
	
	# View as an infinite image, mirrored at the edges which is ideal for Gaussians  
	ipExtended = Views.extendMirrorSingle(ipRAI)
	
	
	# sigmaSmaller ==> Size of the smaller dots (in calibrated units)
	# sigmaLarger ==> Size of the bigger dots (in calibrated units)
	# minPeakValue ==> Intensity above which to look for dots
	# calibration = [1.0 for i in range(ip1_1.numDimensions())]
	cal = imp.getCalibration()
	calibration = [cal.pixelWidth] # must be a double array 
	extremaType = DogDetection.ExtremaType.MINIMA # bright spots on black background
	sigmaSmaller = 1
	sigmaLarger = 1.475
	minPeakValue = 8.1
	normalizedMinPeakValue = False
	
	detection = DogDetection(ipExtended, ipRAI, calibration, sigmaSmaller, sigmaLarger,
	  				   extremaType, minPeakValue, normalizedMinPeakValue)
	
	
	peaks = detection.getPeaks()
	#IJ.log(str(len(peaks)) + " peaks were found.")
	
	# A temporary array of integers, one per dimension the image has
	peakCoords = zeros(ipRAI.numDimensions(), 'i')
	
	#roiPeaks = new PointRoi()   
	

	
	# Add points to ROI manager
	rm = RoiManager.getInstance()
	if not rm:
	    rm = RoiManager()
	rm.reset()
	
	for peak in peaks:
		# Read peak coordinates into an array of integers
		peak.localize(peakCoords) # in pixels
		proi = PointRoi(peak.getDoublePosition(0), peak.getDoublePosition(1))
		proi.setPosition(int(peak.getDoublePosition(2)))
		rm.addRoi(proi)
		
	peaksTable.addValue("File",imp.getTitle())
	peaksTable.addValue("Peak Number", str(len(peaks)))
	
	
	# Show all ROIs on the image
	rm.runCommand(imp, "Show All")
	
	#COMMENT THIS LINE IN OR OUT DEPENDING ON IF YOU WANT TO SEE THE ROIS ON THE IMAGES
	imp.close()
	
	
run()

