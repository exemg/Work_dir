a=1434.54639+186.8785;
b=-0.67769;
c= 0.014547;
d=-0.0010037;
e= 8.34032e-6;
f=-2.99323e-8;
g= 5.11164e-11;
h=-3.40674e-14;

for(i=1; i<=nSlices; i++){
	v = a+b*i+c*pow(i,2)+d*pow(i,3)+e*pow(i,4)+f*pow(i,5)+g*pow(i,6)+h*pow(i,7);
	run("Subtract...", "value="+v); 
}
