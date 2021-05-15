#mupro.r

fetch_seq=function(pdb,chain)
{
	return(paste(pdbseq(pdb,inds=atom.select(pdb,chain=chain,elety="CA",verbose=FALSE)),collapse="")	);
}

fetch_wildresInd=function(pdb,chain,resno)
{
	inds=atom.select(pdb,chain=chain,elety="CA",verbose=FALSE);
	return(which(inds$atom==atom.select(pdb,chain=chain,elety="CA",resno=resno,verbose=FALSE)$atom));
}

mupro=function(pdb,chain,resno,wildres,mutres,runfile)
{
	temp.seq=fetch_seq(pdb,chain);
	cat("Trial",file=runfile,sep="",eol="\n");
	cat(temp.seq,file=runfile,sep="",eol="\n",append=TRUE);
	resno=fetch_wildresInd(pdb,chain,resno);
	cat(resno,file=runfile,sep="",eol="\n",append=TRUE);
	cat(wildres,file=runfile,sep="",eol="\n",append=TRUE);
	cat(mutres,file=runfile,sep="",eol="\n",append=TRUE);
	#print(readLines(runfile));
	command=paste(mupropath," ./",runfile,sep="",collapse="");
	#print(getwd());
	#print(command);
	result=system(command,ignore.stderr=FALSE,ignore.stdout=FALSE,intern=TRUE)
	#print(result);
	if(is.list(attributes(result))){return(c("-","-"));}
	labelpos=regexec("The mutation (.*) the stability of the protein.",result[1])[[1]][2];
	mupro.label=substr(result[1],labelpos,labelpos+7);
	
	labelpos=regexec("Confidence score: ([-\\.0-9]+)",result[2])[[1]][2];
	mupro.value=substr(result[2],labelpos,labelpos+5);
	
	return(c(mupro.label,mupro.value));
}
