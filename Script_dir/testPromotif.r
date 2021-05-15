library("bio3d");

res <- evalWithTimeout({
       system(paste(exepath," ",pdbfile,sep=""),ignore.stderr = FALSE,ignore.stdout = FALSE,intern=TRUE,wait=FALSE)
	print("hello");
     }, timeout=1.08,elapsed=1.08, onTimeout="warning");



