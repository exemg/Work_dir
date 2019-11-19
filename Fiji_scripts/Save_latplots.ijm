//loadmainpath = getDirectory("Choose a loading directory..."); //select the folder containing the files or the folder files
//maindir_list = getFileList(loadmainpath);
savepath = getDirectory("Select a saving directory...");	//select the destination directory
id = getImageID();
//name = "Zprof_pks_raw";
name = "Latplot_camPC0_net";
L = 7;
tol = 1000;

getDimensions(width, height, channels, slices, frames);
//print("Frames: " + toString(nSlices));
midn = floor((nSlices+1)/2);
setSlice(midn);

// WE FIND THE MAXIMA COORDINATES AND SAVE THE TABLE VALUES IN TWO ARRAYS
//run("Find Maxima...", "prominence=8000 output=List");
maxlist = "prominence=" + toString(tol) + " output=List";
run("Find Maxima...", maxlist);
Npoints = Table.size; //returns the number of point ROIs
X = Table.getColumn("X");
Y = Table.getColumn("Y");
close("Results");

// CONTROL SNIPPET TO CHECK THAT THE COORDINATES ARE 
// ACTUALLY DIFFERENT BETWEEN EACH ELEMENT
//print("Let's print X");
//for(idx=1;idx<5;idx++){
//	print(X[idx]);
//}
//print("Let's print Y");
//for(idx=1;idx<5;idx++){
//	print(Y[idx]);
//}

Table.create("Coord_plot_values");			//creates an empty table for x/y values
for (i = 0; i < Npoints; i++) {	//for each ROI...
	XX = X[i];
	YY = Y[i];
	c = 0;
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
		selectImage(id);
		makePoint(XX, YY);
		run("Plot Z-axis Profile");	//gets the z-axis profile
		rename("plot");
		Plot.getValues(zv,Iv);			//gets the z/Int values
		close("plot");
		IMI = Array.findMaxima(Iv,tol);
		//zM = zv[IMI];
		setSlice(IMI[0]+1);
		
		makeLine(XX-L,YY,XX+L,YY); 
		run("Plot Profile");
		rename("plot");
		Plot.getValues(xv,Ivx);			//gets the z/Int values
		close("plot");
		
		selectImage(id);
		makeLine(XX,YY-L,XX,YY+L); 
		run("Plot Profile");
		rename("plot");
		Plot.getValues(yv,Ivy);			//gets the z/Int values
		close("plot");
		
		selectWindow("Coord_plot_values");
		if (i==0) {
			Table.setColumn("X", xv);
		}
		Table.setColumn("IX" + toString(i+1), Ivx);		//fills the table i-th y values
		Table.setColumn("IY" + toString(i+1), Ivy);
		//resets the plot table
		//Table.reset("Values");		//resets the table values
		//close();					//closes the plot window
		selectImage(id);
		setSlice(midn);
	}
}
//close("Values");
Table.save(savepath + name + ".csv");
//close("Results");