helixCompo=function(DSSP="")
{
	chain=unique(DSSP$ch)
	Hfreq="";
	for(ch in chain)
	{
	index=which(DSSP$cha==ch & DSSP$ss=="H")
	Hres=matrix(DSSP$aa[index],ncol=1); # Helix residues
	Hresfreq=entropy(Hres)$freq # 22 freq
	aaLabel=rownames(Hresfreq)[1:20] #20 aa label
	Hresfreq=round(Hresfreq[1:20],2) # Helix 20 residue frequencies
	Hfreq=cbind(Hfreq,Hresfreq);
	}
	Hfreq=Hfreq[,-1];
	rownames(Hfreq)=aaLabel;
	colnames(Hfreq)=chain;
return(t(Hfreq));
}
outSecstr=function(secstr,outsscomp)
{
	cat("# Secondary structure composition\n\nTotRes = Total number of residue\nss = Type of secondary structure\nH = alpha helix\nB = residue in isolated beta-bridge\nE = extended strand, participates in beta ladder\nG = 3-helix (3/10 helix)\nT = hydrogen bonded turn\nS = bend\nX = Any non standard amino acid residue\n\nType\tPercentage\n-------------------------",file=outsscomp,sep="\t",eol="\n");
	write.table(cbind(names(secstr),secstr),file=outsscomp,sep="\t",eol="\n",append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE);
}
outSSaadistr=function(sscomp,outsscomp)
{
	cat("\nAmino acid distribution in each secondary structure\n---------------------------------------------------------\n",file=outsscomp,sep="",eol="",append=TRUE);
	write.table(sscomp,file=outsscomp,sep="\t",eol="\n",append=TRUE,col.names=FALSE,quote=FALSE);

}
plotSecstr=function(secstr,outimg)
{
#print(outimg);
png(file=outimg); 
mx=max(secstr[2:8])+10
barplot(secstr[2:8],xlab="Secondary structure type",ylim=c(0,mx),ylab="Percentage",main="DSSP secondary structure composition",sub=paste("Number of residue: ",secstr[1],sep="",collapse=""));
}
plotSSaadistr=function(ssaacomp,ssaacompimg,type)
{
png(file=ssaacompimg);
mx=max(as.numeric(ssaacomp[3:22]));
mx=ceiling(mx/5)*5
barplot(as.numeric(ssaacomp[3:22]),names.arg=names(ssaacomp[3:22]),ylim=c(0,mx),xlim=c(0,22),xlab="Amino acid residues",ylab="Percentage",main=paste(type," composition",sep="",collapse=""),sub=paste("Total residue in ",type,": ",ssaacomp[2],sep="",collapse=""));
}

secstr_aadistr=function(DSSP="",sstype)
{
	if(length(DSSP$ss)!=length(DSSP$aa)){stop("length(DSSP$ss)!=length(DSSP$aa), check str(DSSP)\n");}
	ssindex=which(DSSP$ss %in% sstype) # indexes of DSSP$ss which is of type sstype
	aa=DSSP$aa[ssindex]; # aa corresponding to that sec. str type
	ss=DSSP$ss[ssindex]
	#aassfreq=round(((tapply(aa,factor(aa),length))*100)/length(aa),2) # For a given sstype, calculates freq by tapply, a vector
	
	aassfreq=tapply(aa,factor(aa),length) # For a given sstype, calculates freq by tapply, a vector
	
	aassfreq=c(ss=paste(sstype,"",sep="",collapse=""),totres=length(aa),aassfreq[c("V","I","L","M","F","W","Y","S","T","N","Q","H","K","R","D","E","A","G","P","C","X")]);
	aassfreq[is.na(aassfreq)]="0";
	names(aassfreq)=c("ss","totres","V","I","L","M","F","W","Y","S","T","N","Q","H","K","R","D","E","A","G","P","C","X")
	return(aassfreq);
}

secstrCompo=function(DSSP="",pdb="")
{
	# if seqres available then use that
	len=1;
	if(length(pdb$seqres)>0)
	{
	len=length(pdb$seqres);
	}else{
	len=length(DSSP$ss); #no of res
	}

sscomp=round((tapply(DSSP$ss,factor(DSSP$ss),length)*100)/len,2)
# keep the order always BCEGHST, if some structure not found, it returns NA so whereever NA comes make it 0
compo=c(len=len,sscomp["B"],sscomp["C"],sscomp["E"],sscomp["G"],sscomp["H"],sscomp["I"],sscomp["S"],sscomp["T"])
compo[is.na(compo)]=0;
names(compo)=c("TotRes","B","C","E","G","H","I","S","T");
return(compo);
}

secstrCompo_chainwise=function(DSSP="")
{
	chain=unique(DSSP$ch)
	compo="";
	for(ch in chain)
	{
	index=which(DSSP$ch==ch)
	len=length(index); #no of res
	#sscomp=round(tapply(DSSP$ss[index],factor(DSSP$ss[index]),length)/len,2)
	sscomp=tapply(DSSP$ss[index],factor(DSSP$ss[index]),length)/len
	
	
	# keep the order always BCEGHST, if some structure not found, it returns NA so whereever NA comes make it 0
	temp=c(ch,sscomp["B"],sscomp["C"],sscomp["E"],sscomp["G"],sscomp["H"],sscomp["S"],sscomp["T"])
	temp[is.na(temp)]="0"
	compo=rbind(compo,temp);
	}
	compo=compo[-1,]
	compo=matrix(compo,ncol=8)
	colnames(compo,do.NULL=FALSE)
	colnames(compo)=c("chain","B","C","E","G","H","S","T");
	return(compo);
}

print_HTMLsecstr_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Secondary structure composition</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("Filename","Total number of residues",paste("Number of residues with secondary structure type: ",c("B","C","E","G","H","I","S","T"),
	c("(Isolated beta-bridge)","(Random coil)","(Extended strand)","(3/10 helix)","(Alpha helix)","(5 helix/pi helix)","(Bend)","(Hydrogen bonded turn)"),sep=""));
	write2Html.tableHeader(c("Filename","Length","B","C","E","G","H","I","S","T"),filepath,append=TRUE,tipvec);
}

print_HTMLsecstrdist_sum_start=function(filepath,type)
{
	write2Html("<html>\n",filepath);
	temptag=paste("<b><u>",type," composition</u><br><br></b>\n",sep="",collapse="");
	write2Html(temptag,filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c(
	"Filename",paste("Secondary structure of type: ",type," (DSSP notation)",sep="",collapse=""),paste("Total number of residues present in ",type,sep="",collapse=""),
	paste("Total number of residue of type:(",c("V","I","L","M","F","W","Y","S","T","N","Q","H","K","R","D","E","A","G","P","C"),") present in ",type,sep=""),
	paste("Total number of non-standard aa present in ",type,sep="",collapse=""));
	write2Html.tableHeader(c("Filename","SS","T.Res","V","I","L","M","F","W","Y","S","T","N","Q","H","K","R","D","E","A","G","P","C","X"),filepath,append=TRUE,tipvec);
}

