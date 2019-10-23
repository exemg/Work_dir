path = getDirectory("Select a directory");	//select the destination directory
Nroi = roiManager("count");		//returns the number of ROIs
name = "Denoised_BL_net_beattrack";
//Table.create("Values");			//creates an empty table for x/y values
for (i = 0; i < Nroi; i++) {	//for each ROI...
	roiManager("select", i);	//selects it
	Roiname = Roi.getName;		//gets the name of the ROI
	run("*Plot Z-axis Profile");	//gets the z-axis profile
	Plot.getValues(x, y);			//gets the x/y values
	//if (i==0) {						//if we start the loop, fill the table x values
	//	Table.setColumn("x", x);
	//	
	//}
	Table.setColumn("y " + Roiname, y);		//fills the table i-th y values
	//saves the plot table
	//Table.reset("Values");		//resets the table values
	close();					//closes the plot window
}
//close("Values");
Table.save(path + name + "_ROIs.csv"); 
close("Results");