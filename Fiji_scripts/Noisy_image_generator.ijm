PSFname = getTitle();
title = 'test_clean_image.tif';
outname = "test_noisy_image.tif";
type = '16-bit black';
width = 1024;
height = 1024;
depth = 13;
z0 = floor((depth-1)/2);

I0 = 100;
N = 25;

NP = 100;

//generate image
newImage(title, type, width, height, depth);

//draw PSFs at random intensities 
for (i = 0; i < NP; i++) {
	x = floor(width*random);
	y = floor(height*random);
	v = floor(10000*random);
	setPixel(x, y, v);
}

//convolve with PSF
run("Convolve 3D", "image="+title+" psf="+PSFname+" extension=[Zero Pad (usually best)] normalize create output="+outname);
selectWindow(outname);

//add background
run("Add...", "value="+toString(I0)+" stack");

//add Poisson noise
run("Insert Poisson noise", "offset=0.000 choose=Poisson");

//apply read noise (Gaussian)
run("Add Specified Noise...", "stack standard="+toString(N));
