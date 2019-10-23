name = "Denoised_strm_BL_net";
run("*Plot Z-axis Profile");		//gets the z-axis profile
Plot.getValues(x, y);				//gets the x/y values
close();
Roi.remove;							//removes the ROI
for (i = 1; i < nSlices+1; i++) {		//for each frame...
	setSlice(i);
	vl = toString(y[i-1]);
	run("Subtract...", "value="+vl); 
	//subtract the average value from   
	//the bkg ROI of the same frame 
}
rename(name);