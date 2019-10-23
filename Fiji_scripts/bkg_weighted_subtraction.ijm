//N = 130; //*For Worm_live/stack/WE7*
//N = 23; //*For Worm_live/stack/WE7*
//N = 76; //*For Worm_fixed/Spherical_larva_3_stk (all 3 stacks)*
N = 70; //*For Worm_fixed\WEfixed2 (filter 2)

//---- d = avgIstart/avgIend --------
//d = ( (943.518/757.142)-1 ) / 89; //*For GFP*
//d = ( (943.518/757.166)-1 ) / 89; //*For GFP* (denoised)
//d = ( (776.674/656.198)-1 ) / 89; //*For blue*
//d = ( (35.032/26.901)-1 ) / 89; //*For DAPI*
//d = ( (776.674/656.198)-1 ) / 89; //*For red*
//d0 = 313.911/408.006; //*For Worm_live/stack/WE7*
//d0 = 577.007/647.048; //*For Worm_live/stack/WL*
//d0 = 510.993/576.240; //*For Worm_fixed/Spherical_larva_3_stk (DAPI)*
//d0 = 126.815/140.510; //*For Worm_fixed/Spherical_larva_3_stk (green)*
//d0 = 120.092/121.082; //*For Worm_fixed/Spherical_larva_3_stk (red)*
d0 = 202.479/190.663; //*For Worm_fixed\WEfixed2 (filter 2)
d = (1-d0)/(N-1);
//stackname = "WE6fixed_stk_RF122_275_dz02c00um_exp1000ms_lowsec_red_medfilt2.tif";
//stackname = "WE6fixed_stack_DAPI_medfilt2_net.tif";
//stackname = "WE6fixed_stack_green_medfilt2_net.tif";
//stackname = "WE6fixed_stack_red_medfilt2_net.tif";
stackname = "WEfixed2_stk_RF175_315_dz02c00um_exp0500ms_lowsec_2_goodstack_filter2_medfilt2.tif";

selectWindow(stackname);
for(i=0;i<N;i++){
	setSlice(i+1);
	factor = d0 + d*i;
	selectWindow("filter2_bkg");
	run("Duplicate...", "title=bkg1");
	run("Multiply...","value="+toString(factor));
	imageCalculator("Subtract", stackname, "bkg1");
	close("bkg1");
	selectWindow(stackname);
}