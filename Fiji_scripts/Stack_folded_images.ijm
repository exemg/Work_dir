// Recursively lists the files in a user-specified directory.
// Open a file on the list by double clicking on it.
// Possible causes of misbehaviour:
// 1. The starting folder doesn't only have subfolders, but image files
// 2. The "Stack" stack doesn't reopen on the blank image when called 

  dir = getDirectory("Choose a Directory ");
  count = 1;
  j = 0;
  stackFiles(dir,j);
  
  
  function stackFiles(dir,j) {
     list = getFileList(dir);  
	 // Be careful with this list! It doesn't refer to the list 
	 // of all the images required for the stack, but to the list
	 // of the files inside the current folder, meaning that "list"
	 // will usually contain a configuration file, an image and a text file,
	 // unless it will contain a group of folders
     
     for (i=0; i<list.length; i++) {
     	// if we have a sub-folder: add 1 to j and invoke the function again
        if (endsWith(list[i], "/")){  
        	j++ 
           stackFiles(""+dir+list[i], j);
        }
        // if we have a regular file: if it is a .tif image, 
        // 			1) open
        // 			2) create Stack
        // 	3 or more) add to stack
        else{
       		if (endsWith(list[i], ".tif")){
            	open(dir+list[i]);
            	selectWindow(list[i]);
            	if (j==2){
            		run("Images to Stack","name=Stack");
            		// no need to close anything, the images 
            		// are already deleted when added to Stack
            	}
            	else{
            		if(j>2){
            			selectWindow("Stack");
            			run("Add Slice");
						selectWindow(dir+list[i]);
						run("Copy");
						selectWindow("Stack"); // Stack should reopen on the blank image
						run("Paste");
						close(dir + list[i]);
            		}
            		
            	}
            	print((count++) + ": " + dir + name);
            	}
            else {print("No images in directory " + dir)}
        }    
     }
  }