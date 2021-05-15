# Automute.r




automute=function(pdbpref,mutvec,automute_mode=0)
{
	file1=paste(automutdir,"dinputs.txt",sep="",collapse="");
	file2=paste(automutdir,"doutputs.txt",sep="",collapse="");
	file3=paste(automutdir,"test_vector.arff",sep="",collapse="");
	#print(file1);
	if(file.exists(file1)){unlink(file1)}
	if(file.exists(file2)){unlink(file2)}
	if(file.exists(file3)){unlink(file3)}
	outpref=paste(pdbpref,"_result.txt",sep="",collapse="");
	mutfile=paste(automutdir,pdbpref,"_mutant.txt",sep="",collpase="")
	writeMutant(mutvec,mutfile=mutfile);
	print(paste("perl ",automutpath," -t 1 -f ",mutfile," -m ",automute_mode," -o ",outpref,sep="",collapse=""));
	result=system(paste("perl ",automutpath," -t 1 -f ",mutfile," -m ",automute_mode," -o ",outpref,sep="",collapse=""),ignore.stderr=FALSE,intern=TRUE)
	#result=system(paste("perl ",automutpath," -t 1 -f ",mutfile," -m ",automute_mode, " >",sep="",collapse=""),ignore.stdout=FALSE,ignore.stderr=FALSE,intern=TRUE)
	
	#print(result);
	#print("---");
	#autodata=readLines(outpref);
	#print(autodata);
	#print("+++");
	
	autodata=read.table(outpref,header=TRUE,stringsAsFactors=FALSE);
	print(autodata);
	if(nrow(autodata)==0){return(c("-","-"))}
	templabel="I";
	if(autodata[1,4]=="Decreased"){templabel="D"};
	unlink(outpref);
	unlink(mutfile);
	return(c(autodata[1,5],templabel));
}



writeMutant=function(mutvec="",ph=7.0,temp=25.0,mutfile="")
{
	write.table(cbind(mutvec,temp,ph),file=mutfile,sep=" ",eol="",quote=FALSE,col.names=FALSE,row.names=FALSE);	
}





