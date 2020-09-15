// TO-DO
// 1. Include tolerance calibration X--DONE--X
// 2. Introduce some testing routine to check everything is going all right
// 3. Give the possibility to manually specify the line ROI (with a live graph)
// 4. write also the position of the points, not only the widths

savepath = getDirectory("Select a saving directory...");	//select the destination directory
Table.create("Comprehensive_Widths_values");			//creates an empty table for z values
tolcoeff = 0.4;

IMLIST = getList("image.titles");
for (j = 0; j < IMLIST.length; j++) {
	id = IMLIST[j];
	//name = "Zprof_pks_raw";
	name = substring(id,0,9) + "sp_";
	L = 25;
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
	run("Gaussian Blur...", "sigma=1");
	// WE FIND THE MAXIMA COORDINATES AND SAVE THE TABLE VALUES IN TWO ARRAYS
	maxlist = "prominence=" + toString(tol) + " output=List";
	run("Find Maxima...", maxlist);
	Npoints = Table.size; //returns the number of point ROIs
	X = Table.getColumn("X");
	Y = Table.getColumn("Y");
	close("Results");
	
	Table.create("Coord_xyplot_values");		//creates an empty table for x/y values
	Table.create("Coord_zplot_values");			//creates an empty table for z values
	Table.create("Widths_values");			//creates an empty table for z values
	FitEq = "y = a + (b-a)*exp(-(x-c)*(x-c)/(2*d*d))";  // defines the fitting equation
	for (i = 0; i < Npoints; i++) {	//for each ROI...
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
			//----------------------------------------------------
			//Finds and saves the Z profile
			run("Plot Z-axis Profile");	//gets the z-axis profile
			rename("plot");
			Plot.getValues(zv,Ivz);			//gets the z/Int values
			close("plot");
			//IMI = Array.findMaxima(Ivz,0.5*tol);
			//zM = zv[IMI];
			//setSlice(IMI[0]+1);
			Array.getStatistics(Ivz, zmin, zmax, zavg, zstd);
			zIGs = newArray(zmin,zmax,zv[floor(zv.length/2)],zv[floor(zv.length/6)]);
			
			//Finds and saves the X/Y profiles
			selectWindow(id);
			makeLine(XX-L,YY,XX+L,YY); 
			run("Plot Profile");
			rename("plot");
			Plot.getValues(xv,Ivx);			//gets the x/Int values
			close("plot");
			Array.getStatistics(Ivx, xmin, xmax, xavg, xstd);
			xIGs = newArray(xmin,xmax,xv[floor(xv.length/2)],xv[floor(xv.length/6)]);
			
			selectWindow(id);
			makeLine(XX,YY-L,XX,YY+L); 
			run("Plot Profile");
			rename("plot");
			Plot.getValues(yv,Ivy);			//gets the y/Int values
			close("plot");
			Array.getStatistics(Ivy, ymin, ymax, yavg, ystd);
			yIGs = newArray(ymin,ymax,yv[floor(yv.length/2)],yv[floor(yv.length/6)]);
			
			//----------------------------------------------------
			// Fills the tables (x/y and z) with the plot profiles
			selectWindow("Coord_zplot_values");
			if (i==0) {
				Table.setColumn("Z", zv);
			}
			Table.setColumn("IZ" + toString(i+1), Ivz);		//fills the table i-th z values
			
			selectWindow("Coord_xyplot_values");
			if (i==0) {
				Table.setColumn("X", xv);
			}
			Table.setColumn("IX" + toString(i+1), Ivx);		//fills the table i-th x values
			Table.setColumn("IY" + toString(i+1), Ivy);		//fills the table i-th y values
			
			//----------------------------------------------------
			// Fits the profiles
			Fid = 3;
			
			Fit.doFit(FitEq,zv,Ivz,zIGs);
			sigz = (Fit.p(Fid))/pxZ;		// selects the "d" parameter, which is the sigma
			FWHMz = 2.35482*sigz;	// calculates the FWHM of the profile
			
			Fit.doFit(FitEq,xv,Ivx,xIGs);
			sigx = (Fit.p(Fid))/pxW;		// selects the "d" parameter, which is the sigma
			FWHMx = 2.35482*sigx;	// calculates the FWHM of the profile
			
			Fit.doFit(FitEq,yv,Ivy,yIGs);
			sigy = (Fit.p(Fid))/pxH;		// selects the "d" parameter, which is the sigma
			FWHMy = 2.35482*sigy;	// calculates the FWHM of the profile
			
			//----------------------------------------------------
			// Saves the width parameter (plus the FWHM)
			selectWindow("Widths_values");
			Table.set("Point", 	i, i+1);	// Creates the index
			Table.set("X", 		i, XX);		// i-th X position
			Table.set("Y", 		i, YY);		// i-th Y position
			Table.set("SX", 	i, sigx);	// i-th sigma x value
			Table.set("FWHMX", 	i, FWHMx);	// i-th FWHM x value
			Table.set("SY", 	i, sigy);	// i-th sigma x value
			Table.set("FWHMY", 	i, FWHMy);	// i-th FWHM x value
			Table.set("SZ", 	i, sigz);	// i-th sigma x value
			Table.set("FWHMZ", 	i, FWHMz);	// i-th FWHM x value
			
			selectWindow(id);
			setSlice(midn);
		}
	}
	
	selectWindow("Coord_xyplot_values");
	Table.save(savepath + "xyplots/" + name + "xyplot.csv");
	close("Coord_xyplot_values");
	
	selectWindow("Coord_zplot_values");
	Table.save(savepath + "zplots/" + name + "zplot.csv");
	close("Coord_zplot_values");
	
	selectWindow("Widths_values");
	Array.getStatistics(Table.getColumn("SX"), 		minSX, 	maxSX, 		avgSX, 	stdSX);
	Array.getStatistics(Table.getColumn("FWHMX"), minFWHMX, maxFWHMX, avgFWHMX, stdFWHMX);
	Array.getStatistics(Table.getColumn("SY"), 		minSY, 	maxSY, 		avgSY, 	stdSY);
	Array.getStatistics(Table.getColumn("FWHMY"), minFWHMY, maxFWHMY, avgFWHMY, stdFWHMY);
	Array.getStatistics(Table.getColumn("SZ"), 		minSZ, 	maxSZ, 		avgSZ, 	stdSZ);
	Array.getStatistics(Table.getColumn("FWHMZ"), minFWHMZ, maxFWHMZ, avgFWHMZ, stdFWHMZ);
	Table.set("Point", 	Npoints, "");	// Additional empty row
	Table.set("Point", 	Npoints+1, "Average");	// Additional "average" row
	Table.set("Point", 	Npoints+2, "Std Dev");	// Additional "Std Dev" row
	Table.set("SX", 	Npoints+1, avgSX);		// Additional avg value
	Table.set("SX", 	Npoints+2, stdSX);		// Additional std value
	Table.set("FWHMX", 	Npoints+1, avgFWHMX);	
	Table.set("FWHMX", 	Npoints+2, stdFWHMX);	
	Table.set("SY", 	Npoints+1, avgSY);		
	Table.set("SY", 	Npoints+2, stdSY);		
	Table.set("FWHMY", 	Npoints+1, avgFWHMY);	
	Table.set("FWHMY", 	Npoints+2, stdFWHMY);	
	Table.set("SZ", 	Npoints+1, avgSZ);		
	Table.set("SZ", 	Npoints+2, stdSZ);		
	Table.set("FWHMZ", 	Npoints+1, avgFWHMZ);	
	Table.set("FWHMZ", 	Npoints+2, stdFWHMZ);	
	Table.save(savepath + "widths/" + name + "widths.csv");
	close("Widths_values");
	
	selectWindow("Comprehensive_Widths_values");
	Table.set("Sample", j, substring(id,0,8));	
	Table.set("Points", j, Npoints);
	Table.set("avgSX", 	j, avgSX);		// Additional avg value
	Table.set("stdSX", 	j, stdSX);		// Additional std value
	Table.set("avgFWHMX", 	j, avgFWHMX);	
	Table.set("stdFWHMX", 	j, stdFWHMX);	
	Table.set("avgSY", 	j, avgSY);		
	Table.set("stdSY", 	j, stdSY);		
	Table.set("avgFWHMY", 	j, avgFWHMY);	
	Table.set("stdFWHMY", 	j, stdFWHMY);	
	Table.set("avgSZ", 	j, avgSZ);		
	Table.set("stdSZ", 	j, stdSZ);		
	Table.set("avgFWHMZ", 	j, avgFWHMZ);	
	Table.set("stdFWHMZ", 	j, stdFWHMZ);
	
	close(id);
}

Table.save(savepath + name + "comprehensive_widths.csv");
close("Comprehensive_Widths_values");