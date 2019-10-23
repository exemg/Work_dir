fname = File.name;
selectWindow(fname);
Nframes = 11;
for(i=0;i<Nframes;i++){
	arg = pow((i-(Nframes-1)/2)/3,2);
	factor = pow(exp(1), -arg);
	print(factor);
	selectWindow(fname);
	run("Duplicate...", "title="+fname+"_"+toString(i));
	run("Multiply...","value="+toString(factor));
}
close(fname);
run("Images to Stack","name="+fname+"_pseudostack title=[] use");