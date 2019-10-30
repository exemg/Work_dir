path = getDirectory("Select a directory");	//select the destination directory
id = getImageID();
run("Find Maxima...", "prominence=1000 output=List");
Npoints = Table.size; //returns the number of point ROIs
X = Table.getColumn("X");
Y = Table.getColumn("Y");
close("Results");
name = "Z_profiles_table";
D1 = 26; //ROI diameter 1
D2 = 20; //ROI diameter 2
D3 = 14; //ROI diameter 3
D4 = 6; //ROI diameter 4
Table.create("ROI_values");			//creates an empty table for x/y values
for (i = 0; i < Npoints; i++) {	//for each ROI...
	selectImage(id);
	makeRectangle(X[i]-D1/2, Y[i]-D1/2, D1, D1);
	run("Measure");
	XM = getResult("XM",0); 
	YM = getResult("YM",0);
	close("Results");
	makeRectangle((XM-D2/2), (YM-D2/2), D2, D2);
	run("Measure");
	XM = getResult("XM",0); 
	YM = getResult("YM",0);
	close("Results");
	makeRectangle((XM-D3/2), (YM-D3/2), D3, D3);
	run("Measure");
	XM = getResult("XM",0); 
	YM = getResult("YM",0);
	close("Results");
	makeOval((XM-D4/2), (YM-D4/2), D4, D4);
	run("Plot Z-axis Profile");	//gets the z-axis profile
	rename("plot");
	Plot.getValues(x,y);			//gets the x/y values
	close("plot");
	selectWindow("ROI_values");
	if (i==0) {						//if we start the loop, fill the table x values
		Table.setColumn("X", x);
		}
	Table.setColumn("Y" + toString(i+1), y);		//fills the table i-th y values
	//resets the plot table
	//Table.reset("Values");		//resets the table values
	//close();					//closes the plot window
}
//close("Values");
Table.save(path + name + "_ROIs.csv");
//close("Results");