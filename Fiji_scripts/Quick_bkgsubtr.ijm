suffix = "mnfilt40net";

title = getTitle();
run("Duplicate...", "duplicate");

title2 = getTitle();
run("Mean...", "radius=40 stack");

imageCalculator("Subtract create stack", title,title2);
close(title);
close(title2);

title3 = getTitle();
LST3 = lengthOf(title3);
title3 = substring(title3, 10, LST3-20)+ suffix + ".tif";
//replace(title3, "MMStack_Pos0.ome", "mnfilt40net");
rename(title3);
wait(50);
run("Properties...");