// N is the number of frames in the stack
// Nref1/2 are the reference frames from where 
//    the feature avg intensity was extracted

N = 71; //*For Worm_fixed\WEfixed (all 4 stacks)*
//Nref1 = 26; //*For Worm_fixed\WEfixed (filter 1)*
//Nref2 = 51; //*For Worm_fixed\WEfixed (filter 1)*
//Nref1 = 6; //*For Worm_fixed\WEfixed (filter 3)*
//Nref2 = 60; //*For Worm_fixed\WEfixed (filter 3)*
Nref1 = 40; //*For Worm_fixed\WEfixed (filter 2)*
Nref2 = 49; //*For Worm_fixed\WEfixed (filter 2)*
//Nref1 = 27; //*For Worm_fixed\WEfixed (filter 4)*
//Nref2 = 50; //*For Worm_fixed\WEfixed (filter 4)*

//---- d0 = avgI(Nref2)/avgI(Nref1) -- d = avgI(n-1)/avgI(n) --------
//d0 = 20742.003/2753.299; //*For Worm_fixed\WEfixed (filter 3)*
//d0 = 1217.255/945.675; //*For Worm_fixed\WEfixed (filter 3/1?)*
//d0 = 597.101/135.458; //*For Worm_fixed\WEfixed (filter 4)*
d0 = 768.078/410.285; //*For Worm_fixed\WEfixed (filter 2)*
d = exp(log(d0)/(Nref2-Nref1)); //It uses log because I didn't find the n-root operator
//print(d0);
//print(d);
//stackname = "WEfixed_stack_filter1_medfilt2_enhanced.tif";
//stackname = "WEfixed_stack_filter2_medfilt2_enhanced.tif";
//stackname = "WEfixed_stack_filter3_medfilt2_enhanced.tif";
//stackname = "WEfixed_stack_filter4_medfilt2_enhanced.tif";

//selectWindow(stackname);
for(i=1;i<N+1;i++){
	setSlice(i);
	factor = pow(d, Nref2-i); //factor=1 for i = Nref2
	//factor = pow(d, Nref1-i); //factor=1 for i = Nref1
	//print(factor);
	run("Multiply...","value="+toString(factor));
}