//Selects the image inside the folder
dir = getDirectory("C:\Users\mg588\Michele\Work\Experiments\20190618_Cairntests");
Filelist0 = getFileList(dir);
print(dir);
suffix= ".tif";
imagename1 = "";
for (i = 0; i < Filelist0.length; i++) {
        if(endsWith(Filelist0[i], suffix)){
            open(Filelist0[i]);
            print("File found");
            imagename1 = imagename1 + Filelist0[i];
        }
        else {
        	print("File not found!");
        }
    }

selectWindow(imagename1);

//Crops the image 
makeRectangle(204, 108, 902, 902);
run("Crop");

//Generates a copy, removes the outliers and subtract
//it to the original to get a background free image
run("Duplicate...", "duplicate");
imagename2 = getTitle();
run("Remove Outliers...", "radius=10 threshold=5 which=Bright stack");
selectWindow(imagename1);
run("Gaussian Blur...", "sigma=1 stack");
imageCalculator("Subtract create stack", imagename1, imagename2);
selectWindow("Result of " + imagename1);
imagename3 = "net_" + imagename1;
rename(imagename3);
Stack.getDimensions(width, height, channels, slices, frames);
Stack.setSlice(25);

//saves the image
saveAs("Tiff", dir+imagename3);

//little bit of cleaning
selectWindow(imagename3);
close("\\Others");

//Finds a few maxima and takes the coordinates of one of them 
run("Find Maxima...", "prominence=50 exclude output=List");
X0 = getResult("X", 0);
Y0 = getResult("Y", 0);
selectWindow("Results");
close();

//uses those coordinates to find the z-Axis Profile, then saves 
//the values and fits them in the middle range (where we think 
//the middle peak is)
makePoint(X0, Y0, "small yellow hybrid");
run("Plot Z-axis Profile");
Plot.showValues();
Plot.getValues(xpoints, ypoints);
xpointsM = Array.slice(xpoints,floor(0.25*slices),floor(0.75*slices)); 
ypointsM = Array.slice(ypoints,floor(0.25*slices),floor(0.75*slices));
Fit.doFit("Gaussian (no offset)", xpointsM, ypointsM);
Fit.logResults;

//run("Plot Profile");
//run("Z Project...", "projection=[Max Intensity]");
//selectWindow("Zstk_SEP00mm_F00_0_1_MMStack_Default.ome.tif");
//run("Z Project...", "projection=[Min Intensity]");

//selectWindow("Zstk_SEP00mm_F00_0_1_MMStack_Default.ome.tif");
//run("Remove Outliers...", "radius=10 threshold=5 which=Bright stack");

//run("Find Maxima...", "prominence=13 exclude output=List");
//String.copyResults();
//selectWindow("Results");
//String.copyResults();
//run("Distribution...");
//run("Input/Output...");

