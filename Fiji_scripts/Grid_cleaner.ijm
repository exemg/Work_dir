Name1 = getTitle();
run("Median...", "radius=1 stack");
run("Duplicate...", "duplicate");
Name2 = getTitle();
run("Remove Outliers...", "radius=15 threshold=6 which=Bright stack");
imageCalculator("Subtract create stack", Name1,Name2);
close(Name1);
close(Name2);
run("Remove Outliers...", "radius=1 threshold=0 which=Bright stack");