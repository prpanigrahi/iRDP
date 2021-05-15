Deamidation_addDSSP=function(DSSP="",deam="")
{
#Args are dssp object, deamidation object
#> deam[1,]
# chain  resno    pat NGdist 
#   "B"  "175"   "NS" "4.27" 

# If there is no chain then deam[,"chain"] are NA but DSSP use " " (space). Hence NA are converted to space.
deam[which(is.na(deam[,"chain"])),"chain"]=" ";
deamChRes=paste(deam[,"chain"],deam[,"resno"],sep="") #combine chain and resno of deamidation object
DSSPChRes=paste(DSSP$cha,DSSP$res,sep="") #combine chain and resno of DSSP object
#acc=0;
SS="";
	for(i in deamChRes)
	{
		index=which(DSSPChRes==i);
		if(length(index)==0){SS=cbind(SS,"-");}else{
		SS=cbind(SS,DSSP$ss[index]);
		}
#		acc=cbind(acc,DSSP$acc[index]);		
		
	}
#	acc=acc[,-1]
	SS=SS[,-1]
	deam=cbind(deam,SS)
#	
	return(deam);
}

Deamidation_addNaccess=function(nacc="",deam="")
{
deam[which(is.na(deam[,"chain"])),"chain"]=" ";
deamChRes=paste(deam[,"chain"],deam[,"resno"],sep="") #combine chain and resno of deamidation object
chres=paste(nacc$asa[,1],nacc$asa[,2],sep=""); #ChainResnoResid format
acc="-";
	for(i in deamChRes)
	{
		index=which(chres==i);
		if(length(index)==0){acc=cbind(acc,"-");}else{
		acc=cbind(acc,nacc$asa[index,5]);
		}
	}
	acc=acc[,-1];
	deam=cbind(deam,acc);
	return(deam);
}

Deamidation_PatFind=function(pdb,chain="INCORRECT_CHAIN",pat="DOESNOTEXIST")
{
	
	pdb$atom[which(is.na(pdb$atom[,"insert"])),"insert"]="";
	#The comments are explained with respect to "NG"
	if(is.na(chain)){chain=""}
	temp=atom.select(pdb,chain=chain,elety="CA",verbose=FALSE)$atom; #pdb$atom CA row/atom index for a given chain
	seq=paste(aa321(pdb$atom[temp,"resid"]),"",sep="",collapse="") # create aa sequence of a given chain ch
	PATPOS=as.vector(gregexpr(pat, seq, fixed=T)[[1]]) #find position of pattern in seq
	# Suppose seq="abcdbc" pat="bc" it returns a vector [2 5]. if pat="xy" returns -1 So no match returns -1  
	# length(temp) and nchar(seq) will be same, so temp[PATPOS] will give pdb$atom index of CA atom of Asn residue assuming pat is "NG"
	# if no pattern found, PATPOS will be -1 hence temp[PATPOS] or temp[-1] will include all index of temp other than 1st index.
	# To prevent this we have put condition that PATPOS[1] if != -1 then only u analyze that pattern.

	if(PATPOS[1]!=-1)
	{
	# Analysis of a pattern is problem when missing residue is there
	# Considering pat="NG" if missing residue e.g. N=93 and G=100 and 94-99 sequence is missing then seq = "aa1...N93G100...aaLAST"
	# N93G100 will be wrongly interpreted as "NG". we can overcome by simply checking whether resno of N and G difference is 1 or not

	Nresno=pdb$atom[temp[PATPOS],"resno"]	#PDB resno of all N, as vector
	Gresno=pdb$atom[temp[PATPOS+1],"resno"] #PDB resno of all G, as vector
	Nresno.ins=pdb$atom[temp[PATPOS],"insert"];
	Gresno.ins=pdb$atom[temp[PATPOS+1],"insert"];
	
	TakeIndex=which(as.numeric(Nresno)-as.numeric(Gresno) == -1); #consider only these indexes where N-G residue number difference is 1

	# Sometimes there is one pattern so it skips PATPOS[1]!=-1 condition then we calculate Nresno,Gresno.
	# in 1H5W.pdb NS exists in chain A, Nresno=229. Gresno=245. Takeindex becomes integer(0) so indirectly no pattern. so we skip further.
		if(length(TakeIndex) > 0)
		{
			Nresno=Nresno[TakeIndex] 
			Gresno=Gresno[TakeIndex]
			Nresno.ins=Nresno.ins[TakeIndex]
			Gresno.ins=Gresno.ins[TakeIndex]
			
			
	
			# Now Nresno and Gresno are the true pattern's resno of N and G i.e first and second residue of pat
			# To calculate the distance we need to consider main chain Nitrogen atom's XYZ of Glycine residue i.e Gresno
			
			Gxyz.tempind=atom.select(pdb,chain=chain,resno=Gresno,elety="N",verbose=FALSE)$atom;

			tempvec=paste(pdb$atom[Gxyz.tempind,"resno"],pdb$atom[Gxyz.tempind,"insert"],sep="");
			
tosearch=matrix(paste(Gresno,Gresno.ins,sep=""),ncol=1);

tempvec.list=apply(tosearch,1,function(x){tempvec.ind=which(tempvec==x); if(length(tempvec.ind)>0){return(tempvec.ind);}else{return("-")} })
Gxyz.tempindNew=rep("-",times=nrow(tosearch));
Gxyz.tempindNew[which(tempvec.list!="-")]=Gxyz.tempind[as.numeric(tempvec.list[which(tempvec.list!="-")])]

			
			GxyzNatom=matrix(as.numeric(pdb$atom[as.numeric(Gxyz.tempindNew[which(Gxyz.tempindNew!="-")]),c("x","y","z")]),ncol=3) #xyz matrix of n*3 where n= no of glycine/ser/ala

			# and CG atoms XYZ coordinate in case of Asn else CD if in case of Gln
			if(substr(pat,1,1)=="N"){ele="CG"}else{ele="CD"}
			
			Nxyz.tempind=atom.select(pdb,chain=chain,resno=Nresno,elety=ele,verbose=FALSE)$atom
			tempvec=paste(pdb$atom[Nxyz.tempind,"resno"],pdb$atom[Nxyz.tempind,"insert"],sep="");
			#Nxyz.tempind=Nxyz.tempind[apply(matrix(paste(Nresno,Nresno.ins,sep=""),ncol=1),1,function(x){tempvec.ind=which(tempvec==x); if(length(tempvec.ind)>0){return(tempvec.ind);}  })]
	tosearch=matrix(paste(Nresno,Nresno.ins,sep=""),ncol=1);
	tempvec.list=apply(tosearch,1,function(x){tempvec.ind=which(tempvec==x); if(length(tempvec.ind)>0){return(tempvec.ind);}else{return("-")} })

Nxyz.tempindNew=rep("-",times=nrow(tosearch));
Nxyz.tempindNew[which(tempvec.list!="-")]=Nxyz.tempind[as.numeric(tempvec.list[which(tempvec.list!="-")])]

			NxyzCGatom=matrix(as.numeric(pdb$atom[as.numeric(Nxyz.tempindNew[which(Nxyz.tempindNew!="-")]),c("x","y","z")]),ncol=3) #xyz matrix of n*3 where n= no of asn/gln
			# Sometimes either N atom of Gly missing or sidechain CG atom of 
	
				if(!(sum(Nxyz.tempindNew=="-")>0 || sum(Gxyz.tempindNew=="-")>0))
				{
					# If there is only one pattern we get then dist.xyz of bio3d package does not work and cbind() also doesnot work
					if(length(Nresno)==1)
					{
					if(nrow(GxyzNatom)==0 || nrow(NxyzCGatom)==0){NGdist="-";}else{NGdist=round(dist(rbind(NxyzCGatom,GxyzNatom)),3)}
					# even if residue is there but if sidechain is missing then
					# Try for 1G3K.pdb andpat="QG"
					deamtemp=matrix(pdb$atom[as.numeric(Nxyz.tempindNew),c("chain","resno","insert")],ncol=3)
					deamtemp[,2]=paste(deamtemp[,2],deamtemp[,3],sep="")
					deamtemp=deamtemp[,c(1:2)]

					deamtemp=c(deamtemp,pat=pat);
					deamtemp=c(deamtemp,dist=NGdist);
					}
					else
					{
					if(nrow(GxyzNatom)==0 || nrow(NxyzCGatom)==0){NGdist="-";}else{NGdist=round(dist.xyz(NxyzCGatom,GxyzNatom,all.pairs=FALSE),3)}
					deamtemp=pdb$atom[as.numeric(Nxyz.tempindNew),c("chain","resno","insert")]
					deamtemp[,2]=paste(deamtemp[,2],deamtemp[,3],sep="")
					deamtemp=deamtemp[,c(1:2)]
										
					deamtemp=cbind(deamtemp,pat);
					deamtemp=cbind(deamtemp,NGdist);
					}
					return(deamtemp);
				}
				else
				{
					for(i in 1:length(Nresno))
					{
					GxyzNatomtemp="";
					if(Gxyz.tempindNew[i]!="-") {GxyzNatomtemp=matrix(as.numeric(pdb$atom[as.numeric(Gxyz.tempindNew[i]),c("x","y","z")]),ncol=3) #xyz matrix of n*3 where n= no of glycine/ser/ala
}else{
GxyzNatomtemp=matrix(numeric(0), 0,0);
}
					NxyzCGatomtemp="";
					if(Nxyz.tempindNew[i]!="-") {NxyzCGatomtemp=matrix(as.numeric(pdb$atom[as.numeric(Nxyz.tempindNew[i]),c("x","y","z")]),ncol=3) }else{
NxyzCGatomtemp=matrix(numeric(0), 0,0);
}
						if(dim(NxyzCGatomtemp)[1] == dim(GxyzNatomtemp)[1])
						{
						NGdist=round(dist(rbind(NxyzCGatomtemp,GxyzNatomtemp)),3)
						deamtemp=matrix(pdb$atom[Nxyz.tempind[i],c("chain","resno","insert")],ncol=3)
						deamtemp[,2]=paste(deamtemp[,2],deamtemp[,3],sep="")
						deamtemp=deamtemp[,c(1:2)]

						deamtemp=c(deamtemp,pat=pat);
						deamtemp=c(deamtemp,dist=NGdist);				
						}
						else
						{
						NGdist="-";
						deamtemp=matrix(pdb$atom[Nxyz.tempind[i],c("chain","resno","insert")],ncol=3)
						deamtemp[,2]=paste(deamtemp[,2],deamtemp[,3],sep="")
						deamtemp=deamtemp[,c(1:2)]
						deamtemp=c(deamtemp,pat=pat);
						deamtemp=c(deamtemp,dist=NGdist);
						}
					}
					return(deamtemp);
				}	
	
		}else{return(NA);} # Even if pattern found but N, G are not sequential, after removing those cases if nothing left then return NA}
	}
	else
	{
	return(NA);
	}
}
Deamidation_protein=function(pdb="")
{
	deam=c(0,0,0,0); #Final structure stores here as a matrix
	chain=unique(pdb$atom[,"chain"]) #chain wise we need to analyze NG
	for(ch in chain)
	{
	print(ch);
	deamtemp=Deamidation_PatFind(pdb,chain=ch,pat="NG"); if(!sum(is.na(deamtemp))>0){deam=rbind(deam,deamtemp);}
	deamtemp=Deamidation_PatFind(pdb,chain=ch,pat="NA"); if(!sum(is.na(deamtemp))>0){deam=rbind(deam,deamtemp);}
	deamtemp=Deamidation_PatFind(pdb,chain=ch,pat="NS"); if(!sum(is.na(deamtemp))>0){deam=rbind(deam,deamtemp);}
	deamtemp=Deamidation_PatFind(pdb,chain=ch,pat="QG"); if(!sum(is.na(deamtemp))>0){deam=rbind(deam,deamtemp);}
	deamtemp=Deamidation_PatFind(pdb,chain=ch,pat="QA"); if(!sum(is.na(deamtemp))>0){deam=rbind(deam,deamtemp);}
	deamtemp=Deamidation_PatFind(pdb,chain=ch,pat="QS"); if(!sum(is.na(deamtemp))>0){deam=rbind(deam,deamtemp);}
	}
	# If no pattern found then deam will be a vector 0,0,0,0. If atleast a single pattern found then it ill be a matrix
	if(!is.matrix(deam))
	{
	return("SORRY");
	}
	deam=deam[-1,];
	deam=matrix(deam,ncol=4);
	colnames(deam)=c("chain","resno","pat","dist");
	return(deam); # returns a matrix
	# tapply(deam[,3],factor(deam[,1]),length) chain wise count of deamidating residue
}

Deamidation_summary=function(deamid)
{
	return(c(Total=nrow(deamid),NGly=sum(deamid[,3]=="NG"),NAla=sum(deamid[,3]=="NA"),NSer=sum(deamid[,3]=="NS"),QGly=sum(deamid[,3]=="QG"),QAla=sum(deamid[,3]=="QA"),QSer=sum(deamid[,3]=="QS"),DistLt4=sum(deamid[,4]<4.0)));
}

print_HTMLhelixDipole=function(hdipo,filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Helix dipole stabilization profile</u><br><br></b>\n",filepath,TRUE);
#	write2Html("Chain1,Resno1,Resid1: Corresponds to Acidic residue details<br>\n",filepath,TRUE);
#	write2Html("Chain2,Resno2,Resid2: Corresponds to Basic residue details<br>\n",filepath,TRUE);
#	write2Html("Resid1_Acc and Resid2_Acc corresponds to relative ASA of acidic and basic residue respectively<br>\n",filepath,TRUE);
#	write2Html("Resid1_SS and Resid2_SS corresponds to secondary structure of the acidic and basic residues respectively<br><br><br>\n",filepath,TRUE);

	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html.tableHeader(c("H.Chain","H.Start","H.End","H.type","H.length","NcapAcid","CcapBase"),filepath,append=TRUE);
	write2Html.tableBody(hdipo,filepath,TRUE,"<td>");
	write2Html("\n</table>\n</html>",filepath,TRUE);
}

print_HTMLdeamid=function(deam,filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Details of Thermolabile bonds</u><br><br></b>\n",filepath,TRUE);
#	write2Html("Chain1,Resno1,Resid1: Corresponds to Acidic residue details<br>\n",filepath,TRUE);
#	write2Html("Chain2,Resno2,Resid2: Corresponds to Basic residue details<br>\n",filepath,TRUE);
#	write2Html("Resid1_Acc and Resid2_Acc corresponds to relative ASA of acidic and basic residue respectively<br>\n",filepath,TRUE);
#	write2Html("Resid1_SS and Resid2_SS corresponds to secondary structure of the acidic and basic residues respectively<br><br><br>\n",filepath,TRUE);

	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("Chain","Start residue","Thermolabile bond","Distance","Secondary structure [DSSP notation]","Relative Solvent accessibility value");
	write2Html.tableHeader(c("Chain","Position","Res.Pair","Distance","SS","Acc"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(deam,filepath,TRUE,"<td>");
	write2Html("\n</table>\n</html>",filepath,TRUE);
}

print_HTMLdeam_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Thermolabile bond profile </u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("PDB Filename","Total number of thermolabile bond",
	paste("Number of thermolabile bond of type: ",c("NG","NA","NS","QG","QA","QS"),sep=""),
	"Number of thermolabile bonds having distance between the n +1 peptide nitrogen and the carbon (CG in case of Asn and CD in case of Gln) containing the side chain amide group less than 4 A ̊"
	);
	write2Html.tableHeader(c("Filename","Total","NG","NA","NS","QG","QA","QS","<4"),filepath,append=TRUE,tipvec);
}

