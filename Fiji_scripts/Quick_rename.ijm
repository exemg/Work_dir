IMLIST = getList("image.titles");
for (i = 0; i < IMLIST.length; i++) {
	name0 = IMLIST[i];
	name1 = replace(name0, ".", "_");
	subname1 = substring(name1,0,indexOf(name1, "lif - "));
	subname2 = substring(name1,indexOf(name1, "lif - ")+6,lengthOf(name1));
	pathname = "C:/Users/mg588/Michele/Work/PROJECT - PSF check/PSFchk_QCtif/";
	filename = subname1 + subname2 + ".tif";
	finalname = pathname + filename;
	//print(filename);
	saveAs("Tiff", finalname);
	close();
}