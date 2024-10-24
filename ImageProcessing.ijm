// Constants
/*micronPerPixel=0.0000023; // Used to normalize optic flow*/
// User input
#@ File (label = "Image directory", style = "directory") input
#@ Boolean (label = "Overwrite", value = false) bOverwrite
#@ Boolean (label = "Show Images", value = false) bShowImages
#@ String (choices={".tif", ".ome.tif", ".avi"}, style="radioButtonHorizontal") inputType

run("ROI Manager...");
run("Clear Results");
run("Set Measurements...", "bounding integrated area_fraction limit redirect=None decimal=6");
setBatchMode(!bShowImages);
processFolder(input);
setBatchMode(false);
print("Finished image processing");
run("Clear Results");
close("ROI Manager");
close("Results");

function processFolder(input)
{	// Scan folders/subfolders/files to find files with correct suffix
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++)
	{	if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + substring(list[i],0,lengthOf(list[i])-1));
		if(endsWith(list[i], inputType) && !endsWith(list[i], "_processed.tif") && !endsWith(list[i], "_graph.tif"))
			processFile(input, list[i]);
	}
}

function processFile(input, file)
{	// Process the images
	fileNm=substring(file, 0, lengthOf(file)-lengthOf(inputType));
	if(File.exists(input + File.separator + fileNm + "_processed.tif") && !bOverwrite)
	{	print("Skipped " + input + File.separator + file);
		wait(500);
		return;
	}
	
	print("\\Clear");
	print("Processing: " + input + File.separator + file);
	/*img=*/open(input + File.separator + file);

	// Get the micron/pixel sizes
	/*getPixelSize(unit,pw,ph,pd);
	if(unit=="µm") pw=pw/micronPerPixel;*/
	
	resetMinAndMax();
	
	// Apply optical flow analysis
	//run("Gaussian Window MSE", "sigma="(4*pw) + " maximal_distance="(7*pw));
	run("Gaussian Window MSE", "sigma=4 maximal_distance=7");
	
	// Split channels
	selectImage(file + " flow vectors");
	run("Split Channels");
	
	// Convert to 8-bit image
	selectImage("C1-" + file + " flow vectors");
    run("8-bit");

	// Analyse the raw pixel intensity and save to file
	run("Set Scale...", "distance=1 known=1 pixel=1 unit=µm");
	setThreshold(1, 255);
    run("Select All");
	roiManager("Add");
	roiManager("Select", 0);
	roiManager("Multi Measure");
	Table.update;
	/*if(Table.columnExists("BX1"))			*/Table.deleteColumn("BX1");
	/*if(Table.columnExists("BY1"))			*/Table.deleteColumn("BY1");
	/*if(Table.columnExists("IntDen1"))		*/Table.deleteColumn("IntDen1");
	/*if(Table.columnExists("MinThr1"))		*/Table.deleteColumn("MinThr1");
	/*if(Table.columnExists("MaxThr1"))		*/Table.deleteColumn("MaxThr1");
	/*if(Table.columnExists("Width1"))		*/Table.renameColumn("Width1","Width");
	/*if(Table.columnExists("Height1"))		*/Table.renameColumn("Height1","Height");
	/*if(Table.columnExists("%Area1"))		*/Table.renameColumn("%Area1","%Area");
	/*if(Table.columnExists("RawIntDen1"))	*/Table.renameColumn("RawIntDen1","RawIntDen");
	Table.update;
	
	// Set NaN to 0.
	ii = Table.size;
	while (ii>0)
	{	ii--;
   		if ( isNaN(Table.get("RawIntDen", ii)) ) Table.set("RawIntDen", ii, 0);
	}
	
	saveAs("Results", input + File.separator + fileNm + "_processed.csv");
	roiManager("Delete");
	
    // Save under filename file_processed
    print("Saving to: " + input + File.separator + fileNm + "_processed.tif");
    saveAs("TIFF", input + File.separator + fileNm + "_processed");
    
    run("Clear Results");
    run("Close All");
    
    wait(1000);
}

function calculateVariables()
{	// Calculate the variables.
	ResultsTable rt = Analyzer.getResultsTable();
	ResultsTable ct = new ResultsTable();
	Analyzer.setResultsTable(ct);
	
	//ct.addValue("Beat Duration",beatDur);
	
	Analyzer.setResultsTable(rt);
}
