// Recursively lists the files in a user-specified directory and change their names to match the name of the folder.
// Open a file on the list by double clicking on it.

  dir = getDirectory("Choose a Directory ");
  count = 1;
  listFiles(dir); 
  newfilename = "";
  print(newfilename);
  function listFiles(dir) {
     list = getFileList(dir);
     for (i=0; i<list.length; i++) {
        if (endsWith(list[i], "/")){
           // needs to take the name of the directory and remove the "/"
           newfilename = list[i];
           print(newfilename);
           LNF = lengthOf(newfilename);
           newfilename = substring(newfilename,0,LNF-1);
           print(newfilename);
           listFiles(""+dir+list[i]);
        }
        else{
       		if (endsWith(list[i], ".tif")){
            	name = newfilename + ".tif";
            	status = File.rename( dir + list[i], dir + name );
   		     	if ( status == false ) { exit( "Could not rename" ); }
   		     	print((count++) + ": " + dir + name);
		     }
		     else if (endsWith(list[i], ".txt")){
		     	name = newfilename + ".txt";
		     	status = File.rename( dir + list[i], dir + name );
   		        if ( status == false ) { exit( "Could not rename" ); }
   		        print((count++) + ": " + dir + name);
		     }
        }
     }
  }