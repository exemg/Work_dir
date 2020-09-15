PSFname = getTitle();
title = 'test_clean_image';
type = '16-bit black';
width = 1260;
height = 1260;
depth = 13;
z0 = floor((depth-1)/2);


I0 = 100;	// background intensity
N = 50;	// noise SRM
NP = 50;	// number of points

// array with the (net) intensities of the beads
Iarr = newArray(50, 100, 150, 200, 350, 500, 700, 1000, 1500, 2000, 3000, 5000, 7000, 10000);


// array with the resampling factors
resarr = newArray(1,2,3,4,5,7,10);


//select the destination directory
savepath = getDirectory("Select a saving directory...");	


//generate empty image
newImage(title, type, width, height, depth);   
   

//draw PSFs at random positions (all intensities = 1)
setSlice(z0);
for (i = 0; i < NP; i++) {
	x = floor(width*(0.1+0.8*random()));
   	y = floor(height*(0.1+0.8*random()));
   	setPixel(x, y, 1);
}


//convolve with PSF
title2 = title+"_deconv";
run("Convolve 3D", "image="+title+" psf="+PSFname+" extension=[Zero Pad (usually best)] create output="+title2);
selectWindow(title2);
close(title);

// "for" loop: intensities 
for(j=0; j<Iarr.length; j++){
	
	IPSF = Iarr[j];			// select intensity
		
	if (IPSF<10000){		// create intensity suffix
		prex = "0";
   	} 
   	if (IPSF<1000){
		prex = prex+"0";
   	}
   	if (IPSF<100){
		prex = prex+"0";
   	}
   	if (IPSF>9999){
		prex = "";
   	}
    	
	
	// create scaled copy (with IPSF intensity)
	selectWindow(title2);
	outname = "Noisetest_"+prex+toString(IPSF);
	run("Duplicate...", "title="+outname+" duplicate");
	selectWindow(outname);
	run("Multiply...", "value="+toString(IPSF)+" stack");
	//saveAs('.tif', savepath+outname+"_clear");
	
	
   	//add background
   	run("Add...", "value="+toString(I0)+" stack");
   	  	
	
	// "for" loop: rescaling
	for(k=0; k<resarr.length; k++){
		
		// selects the desired rescaling
		res = resarr[k];
		
		if (res<10){
			sux = "0";
		}
	   	if (res>10){
			sux = "";
   		}
		
		
		// creates rescaled duplicate
		resoutname = outname+"_"+sux+toString(res);
		selectWindow(outname);
		run("Duplicate...", "title="+resoutname+" duplicate");
		selectWindow(resoutname);
		run("Bin...", "x=res y=res z=1 bin=Average");
		
		
		//add Poisson noise
   		run("Insert Poisson noise", "offset=0.000 choose=Poisson");
		
		
		//apply read noise (Gaussian)
   		run("Add Specified Noise...", "stack standard="+toString(N));
		
		
		// saves and closes image in .tif and .jpg formats		
   		saveAs('.tif', savepath+resoutname);
		jpgoutname = "Example_"+prex+toString(IPSF)+"_"+sux+toString(res);
		run("Duplicate...", "title="+jpgoutname);
		saveAs('.jpg', savepath+jpgoutname);
   		close(resoutname+".tif");
		close(jpgoutname+".jpg");
	}
	close(outname);
}

close(title2);
