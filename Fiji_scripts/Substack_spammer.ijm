Nframes = 20;
Nstk = 401;
offset = 10
Nsubstk = floor((Nstk-Nframes)/offset); 
namebase = getTitle();
LIO = lastIndexOf(namebase, ".tif");
namebase = substring(namebase, 0, LIO);
for(i=1;i<=Nsubstk;i++){
	n1 = 1+offset*(i-1); //The starting point moves by "offset" points
	n2 = n1+Nframes-1;	 //The stack is big a "Nframes" tall file
	sn1 = padstr(toString(n1-1),3);
	sn2 = padstr(toString(n2-1),3);
	promptstring = "  slices=" + toString(n1) + "-" + toString(n2);
	run("Make Substack...", promptstring);
	newname = namebase + "_RF" + sn1 + "_" + sn2 + ".tif";
	rename(newname);
	newpath = "F:\\Aurox Clarity\\Paperable\\0819_calibration\\Hi_sec\\" + newname;
	save(newpath); 
	close(newname);
}
print("Finished!");

function padstr(string,FinalLength) {
	L = lengthOf(string);
	D = FinalLength - L;
	if (D < 0) {
		exit("error: string input must be shorter than the FinalLength parameter")
	} else {
		if (D > 0) {
			prefix = "";
			for (i = 0; i < D; i++) {
				prefix = prefix+"0";
			}
			string = prefix + string;
		}
		return string;
	}
}