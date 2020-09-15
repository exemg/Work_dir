// The script iterates over the active windows. For each one of them:
// The focal plane is identified, and the features are localized through maxima search
// for each of these features, a stack-of-interest (SOI) is selected and the
// MetroloJ plugin is called upon

//selects the destination directory
savepath = getDirectory("Select a saving directory...");				
// for the maxima search
tolcoeff = 0.4;
wvl = 515;
NA = 1.4;
ph = 1;
sc = 1;
// for the MetroloJ call
micstr = "microscope=Confocal";
wvstr = " wavelength="+toString(wvl);
nastr = " na="+toString(NA);
pinstr = " pinhole="+toString(ph);
txtstr = " text1=[Sample infos:\n] text2=Comments:\n";
scstr = "scale="+toString(sc);

IMLIST = getList("image.titles");
for (j = 0; j < IMLIST.length; j++) {
	id = IMLIST[j];
	tifix = indexOf(id, ".tif");
	//name = "Zprof_pks_raw";
	name = substring(id,0,tifix);
	L = 60;
	//tol0 = 50;
	selectWindow(id);
	getDimensions(width, height, channels, slices, frames);
	getVoxelSize(pxW, pxH, pxZ, unit);
	
	tol = 0;
	imstd = 0;
	midn = 0;
	// Now gets the stdev of each slice and returns the focal plane index
	for (i = 1; i < slices+1; i++) {
		setSlice(i);
		run("Measure");
		std0 = getResult("StdDev",0);
		Ma0 = getResult("Max", 0);
		Mi0 = getResult("Min", 0);
		D0 = Ma0-Mi0;
		if (std0 > imstd) {
			imstd = std0;
			tol = tolcoeff*(D0-(5*std0));
			midn = i;
		}
		close("Results");
	}
	//midn = 6; //floor((nSlices+1)/2);
	setSlice(midn);

	run("Duplicate...", "title=Imagecopy");
	selectWindow("Imagecopy");
	run("Gaussian Blur...", "sigma=2");
	// WE FIND THE MAXIMA COORDINATES AND SAVE THE TABLE VALUES IN TWO ARRAYS
	maxlist = "prominence=" + toString(tol) + " output=List";
	run("Find Maxima...", maxlist);
	Npoints = Table.size; //returns the number of point ROIs
	X = Table.getColumn("X");
	Y = Table.getColumn("Y");
	close("Results");

	tablename = "Coord_xy_"+id;
	Table.create(tablename);		//creates an empty table for x/y position values
	Table.setColumn("X", X);		//fills the table i-th x values
	Table.setColumn("Y", Y);		//fills the table i-th y values
	Table.save(savepath + tablename + ".csv");
	close(tablename);
	// Table.create("Coord_z_positions");		//creates an empty table for z position values
	for (i = 0; i < Npoints; i++) {	//for each ROI...
		istr = toString(i+1);
		if (lengthOf(istr)==1){
			istr = "0"+istr;
			}
			
		XX = X[i];
		YY = Y[i];
		c = 0;
	
		// small section to remove edge maxima
		if(XX>L){
			c = c+1;
		}
		if(YY>L){
			c = c+1;
		}
		if(XX<(width-L)){
			c = c+1;
		}
		if(YY<(height-L)){
			c = c+1;
		}
		
		if (c==4){
			//selectWindow(id);
			
			// small snippet to centre the ROI to the peak center of mass
			Xoff = XX-L;
			Yoff = YY-L;
			makeRectangle(Xoff, Yoff, 2*L, 2*L);
			run("Measure");
			XX = floor(getResult("XM", 0)/pxW);
			YY = floor(getResult("YM", 0)/pxH);
			close("Results");
			close("Imagecopy");

			selectWindow(id);
			makePoint(XX, YY);
			Xoff = XX-L;
			Yoff = YY-L;
			makeRectangle(Xoff, Yoff, 2*L, 2*L);

			run("Duplicate...", "title=imtemp duplicate range=1-"+toString(slices));
			selectWindow("imtemp");
			//----------------------------------------------------
			//Runs the MetroloJ plugin
			svstr = " save save=["+savepath+"PSFrep_"+id+"_"+istr+".pdf]";
			Metrolostring = micstr + wvstr + nastr + pinstr + txtstr + scstr + svstr;
			run("Generate PSF report", Metrolostring);
			close("imtemp");
			
			//----------------------------------------------------
			selectWindow(id);
			setSlice(midn);
		}
	}
	
	close(id);
}