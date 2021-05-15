AAcomp=function(pdb,from="",filename="")
{
# This function is latest one. # It takes only one pdb object and returns freq as vector
# Output is Length,20 aa percentage and other combined
# You need to loop over all pdb files
#print("Running AA compo");
	
	print(str(pdb))
	print("*");
  print(length(pdb$seqres))
	if(length(pdb$seqres)>0)
	{
	temp=aa321(pdb$seqres) #seqres to oneletter
	print(temp);
	}
	else
	{
	#temp=aa321(pdb$atom[which(pdb$calpha),"resid"]) #Take atom tags residue information
	  print(atom.select(pdb,elety="CA",verbose=FALSE)$atom);
    temp=aa321(pdb$atom[atom.select(pdb,elety="CA",verbose=FALSE)$atom,"resid"]);
  
	}
	
	####### For percentage
		#dim(temp)=c(length(temp),1) #convert vector to array
		#tempEnt=entropy(temp) #input array to entropy
		#freq=tempEnt$freq #store entropy_freq to freq
		#freq=freq[-21,]; #remove gap
		#freq=c(len=length(temp),freq*100); #make percentage
		
	####### For number
		
		len=length(temp)
		tempfreq=tapply(temp,factor(temp),length);
		freq=c(Total=len,tempfreq["V"],tempfreq["I"],tempfreq["L"],tempfreq["M"],tempfreq["F"],tempfreq["W"],tempfreq["Y"],tempfreq["S"],tempfreq["T"],tempfreq["N"],tempfreq["Q"],tempfreq["H"],tempfreq["K"],tempfreq["R"],tempfreq["D"],tempfreq["E"],tempfreq["A"],tempfreq["G"],tempfreq["P"],tempfreq["C"],tempfreq["X"])
		names(freq)=c("len","V","I","L","M","F","W","Y","S","T","N","Q","H","K","R","D","E","A","G","P","C","X");
		freq[is.na(freq)]=0;
		freq=c(freq,aro=(freq["F"]+freq["W"]+freq["Y"]),NQST=(freq["N"]+freq["Q"]+freq["S"]+freq["T"]),RKH=(freq["R"]+freq["K"]+freq["H"]),DE=(freq["D"]+freq["E"]),RtoK=(freq["R"]/freq["K"]))
		freq=round(freq,2)
		if(is.infinite(freq["RtoK"])){freq["RtoK"]="-"}
		#freq=c(filename,freq);
		names(freq)=c("len","V","I","L","M","F","W","Y","S","T","N","Q","H","K","R","D","E","A","G","P","C","X","aro","NQST","RKH","DE","RtoK")
		return(freq);
}
outAAcomp=function(aacomp,outaacomp)
{
	cat("# Amino acid composition\n\nAmino acid\tPercentage\n-------------------------\n",file=outaacomp,sep="\t",eol="\n");
	write.table(cbind(names(aacomp),aacomp),file=outaacomp,sep="\t",eol="\n",append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE);
}

plotAAcomp=function(aacomp,outimg)
{
	png(file=outimg);
	mx=max(aacomp[3:22])+2;
	barplot(aacomp[3:22],xlim=c(0,22),ylim=c(0,mx),xlab="Amino acid residues",ylab="Percentage",main="Amino acid composition",sub=paste("Total residue: ",aacomp[2],sep="",collapse=""));
	dev.off();
}

print_HTMLaacomp_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Amino acid composition</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("Filename","The chains considered for aa composition analysis","Total number of residues",paste("Total number of residue of type:",c("V","I","L","M","F","W","Y","S","T","N","Q","H","K","R","D","E","A","G","P","C"),sep=""),"Total number of non-standard aa",
	"Total number of aromatic aa [FYW]","Uncharged polar residues [N,Q,S,T] content","R+K+H content","DE content","R to K ratio");
	write2Html.tableHeader(c("Filename","Chain","Length","V","I","L","M","F","W","Y","S","T","N","Q","H","K","R","D","E","A","G","P","C","X","Aro","NQST","RKH","DE","R/K"),filepath,append=TRUE,tipvec);
}
