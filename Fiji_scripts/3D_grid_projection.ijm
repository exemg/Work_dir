Imagename = getTitle();
dirname = getDirectory("Image");

run("Gaussian Blur...", "sigma=1 stack");
run("Z Project...", "projection=[Min Intensity]");
selectWindow(Imagename);
run("Z Project...", "projection=[Max Intensity]");
Imname1 = "MAX_" + Imagename;
Imname2 = "MIN_" + Imagename;
imageCalculator("Subtract create", Imname1, Imname2);

selectWindow("Result of " + Imname1);
rename("MAX-MIN_" + Imagename)

selectWindow(Imname1);
close();
selectWindow(Imname2);
close();
selectWindow(Imagename);
close();