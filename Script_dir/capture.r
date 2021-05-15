capture=function(prefix="",exepath="",outfile="",filepath="",mode=1)
{
	# prefix is tempfile which will be deleted, outfile is the individual cation-pi detailed file
	# filepath will be "1abc.pdb" or "absolutepath of file to be uploaded"
	file1=paste(prefix,".txt",sep="",collapse="");
	
	# mode 1 means pdb mode else upload mode
	if(mode==1)
	{
	# in this case the filepath should be "1abc.pdb" so that using substring(x,1,4) we get the pdb id for which capture to run
	system(paste("perl ",exepath," ",substr(filepath,1,4)," ",file1," ", outfile," pdb",sep="",collapse=""));

	}else{
	# in this case the filepath should be "absolute path of pdb file" so that it will be uploaded for which capture to run
	system(paste("perl ",exepath," ",filepath," ",file1," ",outfile," XXX",sep="",collapse=""));
	}
	raw.lines <- readLines(outfile);
	unlink(c(file1));
	
	type<- which(substring(raw.lines, 1, 4)=="Type")
	sigf=as.numeric(t(matrix(unlist(strsplit(raw.lines[(type+1):(type+6)],"\t")),nrow=3))[,3])
	# If no cation-pi found then sigf=NANANA etc so if sum is zero then we can return a vector of 0 of length 6
	if(sum(sigf,na.rm=TRUE)==0){sigf=rep(0,times=6);}

	ind=grep("There is a total",raw.lines)
	x=regexpr("There is a total of (\\d+).*",raw.lines[ind],perl=TRUE)
	tot=as.numeric(substr(raw.lines[ind],attributes(x)$capture.start,(attributes(x)$capture.start+attributes(x)$capture.length-1)))
	# If no cation-pi found then tot=numeric(0)
	if(length(tot)==0){tot=0;}

	ind=grep('E(es) <= -2.0',raw.lines,fixed=TRUE)
	x=regexpr("\\* There are (\\d+).*",raw.lines[ind],perl=TRUE)
	less=as.numeric(substr(raw.lines[ind],attributes(x)$capture.start,(attributes(x)$capture.start+attributes(x)$capture.length-1)))
	if(length(less)==0){less=0;}

	ind=grep('between -2.0',raw.lines,fixed=TRUE)
	x=regexpr("\\* There are (\\d+).*",raw.lines[ind],perl=TRUE)
	betw=as.numeric(substr(raw.lines[ind],attributes(x)$capture.start,(attributes(x)$capture.start+attributes(x)$capture.length-1)))
	if(length(betw)==0){betw=0;}

	ind=grep('There is an average of',raw.lines,fixed=TRUE)
	x=regexpr("There is an average of (\\d+).*",raw.lines[ind],perl=TRUE)
	ave=as.numeric(substr(raw.lines[ind],attributes(x)$capture.start,(attributes(x)$capture.start+attributes(x)$capture.length-1)))
	if(length(ave)==0){ave=0;}
	
	return(c(tot,less,betw,ave,sigf));

}
capture_parse=function(prefix="",exepath="",outfile="",filepath="",mode=1)
{
	catpimat="";
	# prefix is tempfile which will be deleted, outfile is the individual cation-pi detailed file
	# filepath will be "1abc.pdb" or "absolutepath of file to be uploaded"
	file1=paste(prefix,".txt",sep="",collapse="");
	
	# mode 1 means pdb mode else upload mode
	if(mode==1)
	{
	# in this case the filepath should be "1abc.pdb" so that using substring(x,1,4) we get the pdb id for which capture to run
	system(paste("perl ",exepath," ",substr(filepath,1,4)," ",file1," ", outfile," pdb",sep="",collapse=""));

	}else{
	# in this case the filepath should be "absolute path of pdb file" so that it will be uploaded for which capture to run
	system(paste("perl ",exepath," ",filepath," ",file1," ",outfile," XXX",sep="",collapse=""));
	}
	raw.lines <- readLines(outfile);
	unlink(c(file1,outfile));
	type<- which(substring(raw.lines, 1, 4)=="Type")
	catpimat=raw.lines[grep("^[A-Za-z]{3}\t\\d+",raw.lines[1:type],perl=TRUE)];
	
	if(length(catpimat)>0) # If more than one catpi is found
	{
	catpimat=matrix(unlist(strsplit(catpimat,"\t")),nrow=length(catpimat),byrow=TRUE);
	colnames(catpimat)=c("Resid1","Resno1","Chain1","Resid2","Resno2","Chain2","Ees","-","Evdw");
	}
	return(catpimat);
}

getCatpi=function(chres="",catpi)
{
	chain=chres[1];
	resno=chres[2];
	if(is.na(chain)){chain="~"}
	query=paste(chain,resno,sep="",collapse="");
	catpi_res1=paste(catpi[,3],catpi[,2],sep="");
	catpi_res2=paste(catpi[,6],catpi[,5],sep="");
	catpi.ind=which(catpi_res1==query | catpi_res2==query);
	if(length(catpi.ind)==0){return("SORRY");} # If not found return SORRY else an matrix
	catpimat=catpi[catpi.ind,];
	catpimat=matrix(catpimat,ncol=ncol(catpi));
	colnames(catpimat)=colnames(catpi);
	return(catpimat)
}

callCatpi=function(pdb)
{
			# This is the workflow that is used in Pattern_analysis code
			catpi="";
			status.catpi=0;
			if(catpimode==1)
			{
			#PDB mode
			#status.catpi=1;
			#catpi.filename=paste(captureprefix,".out",sep="",collapse="");
			#catpi=capture(prefix=captureprefix,exepath=capturepath,outfile=catpi.filename,filepath=f,mode=1);
			#if(is.matrix(catpi)){status.catpi=1;}
			}else{
			# Upload features
			
			}
return(list(catpi=catpi,status.catpi=status.catpi));			
}

checkCatpi=function(patmat,catpi)
{
	# This is the workflow that is used in Pattern_analysis code
	catpimat=do.call(rbind,lapply(apply(patmat,1,getCatpi,catpi),matrix,ncol=ncol(catpi),byrow=FALSE));
	sorry.ind=which(apply(catpimat,1,function(x){sum(x=="SORRY")>0}));
	if(length(sorry.ind)>0){catpimat=catpimat[-sorry.ind,];}
	catpimat=matrix(catpimat,ncol=ncol(catpi));
	return(catpimat);					
}

