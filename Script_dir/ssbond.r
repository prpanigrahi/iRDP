convertPDB=function(pdb,outpdb,outfile)
{
	resno=as.numeric(pdb$atom[,"resno"])
	ind=resno[1]
	temp=0;
	val=1;
	for(i in resno)
	{
		if(i == ind)
		{
		temp=cbind(temp,val);
		}
		else
		{
		val=val+1;
		ind=i;
		temp=cbind(temp,val);
		}
	}

	temp=temp[,-1]
map=pdb$atom[,"chain"]
map=rbind(map,pdb$atom[,"resno"])
map=rbind(map,temp)
map=t(map);
write.pdb(pdb,resno=map[,3],file=outpdb)
cat("",file=outfile,sep="",eol="");
cat("Chain\tOriginal\tRenumbered\n",file=outfile,sep="",eol="",append=TRUE);
write.table(unique(map),file=outfile,sep="\t",eol="\n",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
}
runSSBOND=function(infile,exepath="",chi3dev=0,maxener=10,maxenerdif=5.0)
{
	cat(c(infile,"mutant.pdb","NO",chi3dev,maxener,maxenerdif),file="SSBOND_OPTION.txt",sep="\n",eol="")
	result=system(	paste(exepath," < SSBOND_OPTION.txt"),ignore.stderr = FALSE, intern=TRUE)
	result=paste(result,sep="",collapse="\n");
	cat(result,file="SSBOND.out",sep="",eol="")
	unlink("SSBOND_OPTION.txt");
	unlink(infile);
	SSBlist=parseSSBOND("SSBOND.out");
	return(SSBlist);
}



convertSSBOND=function(SSBlist="",mapfile)
{
	map=read.table(mapfile,skip=1,stringsAsFactors=FALSE)
	unlink(mapfile);

	for(i in 1:length(SSBlist[,1]))
	{
  
		map.index=which(map[,3]==SSBlist[i,"Resno_1"])
		SSBlist[i,"Chain_1"]=map[map.index,1]
		SSBlist[i,"Resno_1"]=map[map.index,2]
		
		map.index=which(map[,3]==SSBlist[i,"Resno_2"])
		SSBlist[i,"Chain_2"]=map[map.index,1]
		SSBlist[i,"Resno_2"]=map[map.index,2]
	}
	return(SSBlist);
}
parseSSBOND=function(infile)
{
	SSB=readLines(infile);
	unlink("SSBOND.out");
	SSB.st=grep("\\s+\\d+\\s+CONFORMATION BETWEEN :",SSB,perl=TRUE)
  if(length(SSB.st)==0)
  {
    return("SORRY");
  }
  
	SSBlist=matrix(NA,ncol=18)
	colnames(SSBlist,do.NULL=FALSE);
	colnames(SSBlist)=c("DS_no","Chain_1","Resno_1","Resid_1","Chain_2","Resno_2","Resid_2","Conf_no","SGdist","Chi1_1","Chi2_1","Chi3","Chi1_2","Chi2_2","Chinrg","Taunrg","Disnrg","Totnrg");
	resno1=as.numeric(substring(SSB[SSB.st],29 ,32))
	resid1=substring(SSB[SSB.st],44 ,46);
	resno2=as.numeric(substring(SSB[SSB.st],36 ,39));
	resid2=substring(SSB[SSB.st],48 ,50)

	for(i in 1:length(SSB.st))
	{
		st=SSB.st[i]	
		end=SSB.st[i+1]
		if(is.na(end))	conf=SSB[(st+2)] else conf=SSB[(st+2):(end-3)]
		tempmat=matrix(NA,nrow=length(conf),ncol=18)
		
		tempmat[,1]=i; # DS pair number
		tempmat[,2]="X"; # chain of res1
		tempmat[,3]= resno1[i] #resno1		
		tempmat[,4]= resid1[i] #resid1
		tempmat[,5]="X"	;# chain of res2
		tempmat[,6]= resno2[i] #resno2
		tempmat[,7]= resid2[i] #resid2
		tempmat[,8]=as.numeric(substring(conf,5 ,6)) #conf no
		tempmat[,9]=as.numeric(substring(conf,9 ,13)) # sgdist
		tempmat[,10]=as.numeric(substring(conf,15 ,20)) #chi11
		tempmat[,11]=as.numeric(substring(conf,22 ,27)) #chi21
		tempmat[,12]=as.numeric(substring(conf,29 ,34)) #chi3
		tempmat[,13]=as.numeric(substring(conf,36 ,41)) #chi12
		tempmat[,14]=as.numeric(substring(conf,43 ,48)) #chi22
		tempmat[,15]=as.numeric(substring(conf,51 ,55)) #chinrg
		tempmat[,16]=as.numeric(substring(conf,58 ,62)) #taunrg
		tempmat[,17]=as.numeric(substring(conf,65 ,69)) #disnrg
		tempmat[,18]=as.numeric(substring(conf,72 ,76)) #totnrg
	
		SSBlist=rbind(SSBlist,tempmat);	
	}
SSBlist=SSBlist[-1,]
return(SSBlist);
}



