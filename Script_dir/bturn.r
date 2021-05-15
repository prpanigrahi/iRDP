detectBturn=function(BT,chain="",pos="",len)
{
	pos=as.numeric(pos);
	start=as.numeric(BT[,"start"])
	chain=BT[,"chain"]
	temp="";
	if(len<=4)
	{
		if(len==2)
		{
		index=which(chain==chain & (start==pos | start==(pos-1) | start ==(pos-2)))	
		temp=BT[index,];
		}else{
			if(len==3)
			{
			index=which(chain==chain & (start==pos | start==(pos-1) ))			
			temp=BT[index,];
			}else{
			index=which(chain==chain & start==pos)			
			temp=BT[index,];			
			}
		}	
	return(temp);
	}
	else
	{
		stop("Pattern cant be longer than 4 residue");
	}
}
parsePROMOTIF_BT=function(infile)
{
### POSSIBLE bugs: delete first 2 line is universally applied or not?

	raw.lines <- readLines(infile)
	raw.lines <- raw.lines[-c(1:2)]  # Delete first 2 line
	chain <- substring(raw.lines, 1, 1)
	BT=chain;	
	start <- substring(raw.lines, 2, 6)
	end = substring(raw.lines, 23, 27)
start=gsub(" ","",start)
end=gsub(" ","",end)

	seq=	paste(substring(raw.lines, 7,7 ),substring(raw.lines, 14,14),substring(raw.lines, 21,21 ),substring(raw.lines, 28,28 ),sep="")
	type=substring(raw.lines,30 ,33 )
	R2_phi=substring(raw.lines,38 ,43 )
	R2_psi=substring(raw.lines,44 ,49 )
	R3_phi=substring(raw.lines,50 ,55 )
	R3_psi=substring(raw.lines,56 ,61 )
	R1R4_CAdist=substring(raw.lines,64 ,66 )
	Hbond=substring(raw.lines,68 ,68 )
	R2_chi1=substring(raw.lines,69 ,74 )
	R3_chi1=substring(raw.lines,75 ,80 )

	BT=cbind(BT,start=start)
	BT=cbind(BT,end=end)
	BT=cbind(BT,seq=seq)
	BT=cbind(BT,type=type)
	BT=cbind(BT,R2_phi=R2_phi)
	BT=cbind(BT,R2_psi=R2_psi)
	BT=cbind(BT,R2_chi1=R2_chi1)

	BT=cbind(BT,R3_phi=R3_phi)
	BT=cbind(BT,R3_psi=R3_psi)
	BT=cbind(BT,R3_chi1=R3_chi1)
	BT=cbind(BT,R1R4_CAdist=R1R4_CAdist)
	BT=cbind(BT,Hbond=Hbond)
	colnames(BT)[1]="chain";
	return(BT);
}
runPROMOTIF=function(pdb="",exepath="",prefix="")
{
	pdbfile=paste(prefix,".pdb",sep="",collapse="");
	bturnfile=paste(prefix,".bturns",sep="",collapse="");
	BT="";
	write.pdb(pdb,file=pdbfile);
	tryCatch({error=system(paste(exepath," ",pdbfile,sep=""),ignore.stderr = FALSE,ignore.stdout = TRUE,intern=TRUE)},error=function(e){return("SORRY")});
	files=list.files(".",pattern=paste("^",prefix,sep="",collapse=""))
	# In 1CGD.pdb promotif does not outputs .bturn file so we cant get BT
	if(sum(files==bturnfile)>0)
	{
	unlink(files[-which(files==bturnfile)]);
	unlink(c("phipsi.mat","promotif2.prm"))
	BT=parsePROMOTIF_BT(bturnfile);
	unlink(bturnfile);
	}
	return(BT);
}

