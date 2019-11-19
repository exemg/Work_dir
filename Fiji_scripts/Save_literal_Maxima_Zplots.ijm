//loadmainpath = getDirectory("Choose a loading directory..."); //select the folder containing the files or the folder files
//maindir_list = getFileList(loadmainpath);
savepath = getDirectory("Select a saving directory...");	//select the destination directory
id = getImageID();
//name = "Zprof_pks_raw";
name = "Zprof_camPC0_pks_net";
tol = 1000;

getDimensions(width, height, channels, slices, frames);
//print("Frames: " + toString(frames));
midn = floor((nSlices+1)/2);
//print(midn);
setSlice(midn);

// WE FIND THE MAXIMA COORDINATES AND SAVE THE TABLE VALUES IN TWO ARRAYS
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
	selectImage(id);
	makePoint(X[i], Y[i]);
	run("Plot Z-axis Profile");	//gets the z-axis profile
	rename("plot");
	Plot.getValues(xv,yv);			//gets the x/y values
	close("plot");
	selectWindow("Coord_plot_values");
	if (i==0) {						//if we start the loop, fill the table x values
		Table.setColumn("X", xv);
		}
	Table.setColumn("Y" + toString(i+1), yv);		//fills the table i-th y values
	//resets the plot table
	//Table.reset("Values");		//resets the table values
	//close();					//closes the plot window
}
//close("Values");
Table.save(savepath + name + ".csv");
//close("Results");