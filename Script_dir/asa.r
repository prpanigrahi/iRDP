asa_summary=function(nacc="")
{
	temppdb=nacc$atomacc;
	asaC=sum(as.numeric(temppdb$atom[atom.select(temppdb,verbose=FALSE,elety=as.character(unique(temppdb$atom[grep("^C",temppdb$atom[,2]),2])))$atom,11]))
	asaN=sum(as.numeric(temppdb$atom[atom.select(temppdb,verbose=FALSE,elety=as.character(unique(temppdb$atom[grep("^N",temppdb$atom[,2]),2])))$atom,11]))
	asaO=sum(as.numeric(temppdb$atom[atom.select(temppdb,verbose=FALSE,elety=as.character(unique(temppdb$atom[grep("^O",temppdb$atom[,2]),2])))$atom,11]))
	asaS=sum(as.numeric(temppdb$atom[atom.select(temppdb,verbose=FALSE,elety=as.character(unique(temppdb$atom[grep("^S",temppdb$atom[,2]),2])))$atom,11]))
	tempNpP=round(as.numeric(nacc$tot[5])/as.numeric(nacc$tot[6]),3);
	if(is.infinite(tempNpP)){tempNpP="-"}
	asasum=c(as.numeric(nacc$tot[-1]),tempNpP,asaC,asaN,asaO,asaS);
	return(asasum);
}

print_HTMLasa_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Solvent Accessible surface area (ASA) summary</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("PDB Filename",
	paste("Absolute ASA of ",c("all atoms","all side chain atoms","all main-chain atoms","all non-polar side chain atoms","all polar side chain atoms"),sep=""),
	"Ratio of absolute ASA of non-polar to polar side chain atoms",paste("Absolute ASA of all ",c("C","N","O","S")," atoms",sep=""),"Relative solvent ASA file by NACCESS"
	);
	write2Html.tableHeader(c("Filename","All-atoms","Total-Side","Main-Chain","Non-polar(NP)","Polar(P)","NP/P","ASA_C","ASA_N","ASA_O","ASA_S","RSA_File"),filepath,append=TRUE,tipvec);
}


asa_perres=function(nacc="",acclow=20)
{
	nacctemp=matrix(nacc$asa[,c(3,5)],ncol=2);
	nacctemp[,1]=aa321(nacctemp[,1]);
	nacctemp[,2]=apply(matrix(nacctemp[,2],ncol=1),1,function(x){if(as.numeric(x)<=acclow){return("B")}else{return("E")};})

	naccres=t(apply(matrix(c("V","I","L","M","F","W","Y","S","T","N","Q","H","K","R","D","E","A","G","P","C"),ncol=1),1,function(x){
	tempvec=nacctemp[which(nacctemp[,1]==x),2];
	tempres=tapply(tempvec,factor(tempvec),length);
	tempres=c(tempres["B"],tempres["E"])
	tempres[which(is.na(tempres))]=0;
	names(tempres)=c("B","E")
	return(c(AA=x,TotRes=length(tempvec),tempres));
	}));

return(matrix(naccres,ncol=4));	
}


print_HTMLasa=function(asa_perres,filepath,append=FALSE,html=TRUE)
{
	if(html){write2Html("<html>\n",filepath,append);}else{write2Html("",filepath,append);}
	write2Html("<br><br><b><u> Per-residue ASA profile</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("The amino acid residue","Total number of aa residue","Number of residues buried","Number of residues exposed");
	write2Html.tableHeader(c("Amino acid","Total","Buried","Exposed"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(asa_perres,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);

	if(html){write2Html("</html>\n",filepath,TRUE);}else{write2Html("",filepath,TRUE);}
}


