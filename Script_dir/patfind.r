patFind_old=function(pdb,ch,pat)
{
	# If no chain information is there, then ch=NA, for that different way of atom.select
	if(is.na(ch) | ch=="~" | ch==" "){atomsel=atom.select(pdb,elety="CA",verbose=FALSE);}else{
	atomsel=atom.select(pdb,chain=ch,elety="CA",verbose=FALSE);
	}

	# no CA atoms then no pattern so u return NA
	if(length(atomsel$atom)==0){return(NA);} 
	
	# seq stores the Main sequence derived from ATOM tag
	seq=paste(seq.pdb(pdb,inds=atomsel),collapse="");

	patind=unlist(gregexpr(pat,seq,perl=TRUE));
	
	# If no pattern found then it returns -1, if found it returns index of first residue of pattern
	if(sum(patind)<1){return(NA);}
	# patind is the position in the seq where pattern is found
	# As we are getting the sequence from atom tags, if any missing residues are there the seq obtained is wrong 
	# Consider pat="NG". if missing residue e.g. N=93 and G=100 and 94-99 sequence is missing then seq = "aa1...N93G100...aaLAST"
	# N93G100 will be wrongly interpreted as "NG". we can overcome by simply checking whether resno of N and G difference is 1 or not

	patlen=nchar(pat); #length of pattern
	patstresno=as.numeric(pdb$atom[atomsel$atom[patind],"resno"])	#PDB resno of all N, as vector
	# If pattern length is 1, then directly return the residue information
	if(patlen==1){return(patstresno);}
	
	patendresno=as.numeric(pdb$atom[atomsel$atom[patind+patlen-1],"resno"])	#PDB resno of all N, as vector
	# If pattern length >1, then we need to check whether st-end=patlen-1 ?
	# inds stores those indexes where st-end is not equal to patlen-1
	# If any such inds found then return the trimmed array else return the whole array
	inds=which((patstresno-patendresno)!= -(patlen-1))
	if(length(inds)>0){return(patstresno[-inds])}else{return(patstresno);}
}

patFind=function(pdb,ch,pat)
{

	# If no chain information is there, then ch=NA, for that different way of atom.select
	if(is.na(ch) | ch=="~" | ch==" "){atomsel=atom.select(pdb,elety="CA",verbose=FALSE);}else{
	atomsel=atom.select(pdb,chain=ch,elety="CA",verbose=FALSE);
	}
	
	# no CA atoms then no pattern so u return NA
	if(length(atomsel$atom)==0){return(NA);} 
	
	# seq stores the Main sequence derived from ATOM tag
	seq=paste(pdbseq(pdb,inds=atomsel),collapse="");
	patlist=gregexpr(pat,seq,perl=TRUE);
	#print(patlist);
	
	patind=unlist(patlist);
	patlen=attributes(patlist[[1]])$match.length
	
	#print(patind);
	#print(patlen);
	
	
	# If no pattern found then it returns -1, if found it returns index of first residue of pattern
	if(sum(patind)<1){return(NA);}
	# patind is the position in the seq where pattern is found
	# As we are getting the sequence from atom tags, if any missing residues are there the seq obtained is wrong 
	# Consider pat="NG". if missing residue e.g. N=93 and G=100 and 94-99 sequence is missing then seq = "aa1...N93G100...aaLAST"
	# N93G100 will be wrongly interpreted as "NG". we can overcome by simply checking whether resno of N and G difference is 1 or not
	pdb$atom[which(is.na(pdb$atom[,"insert"])),"insert"]="";
	
	#print("+++");
	
	tempList=list();
	list.ind=1;

# Loop over each pattern		
for(i in 1:length(patlen))
{
	print(c(patlen[i],patind[i]));
	tempind=atomsel$atom[patind[i]:(patind[i]+patlen[i]-1)];
	#print(tempind);
	noofind=length(tempind);
	tempmat=matrix(pdb$atom[tempind,c("resno","insert","resid")],ncol=3);
	#print(tempmat);
	if((as.numeric(tempmat[noofind,1])-as.numeric(tempmat[1,1])) > (patlen[i]-1))
	{
	print("missing");
	}else{
		if((as.numeric(tempmat[noofind,1])-as.numeric(tempmat[1,1])) < (patlen[i]-1))
		{
		#insertion
		tempList[[list.ind]]=tempmat;
		list.ind=list.ind+1;
		print("insertion")
		}else{
		mat1=as.numeric(tempmat[,1]);
		mat2=mat1[2:length(mat1)]
			if(sum(mat2-mat1[1:length(mat1)-1]>1)==0)
			{
			tempList[[list.ind]]=tempmat;
			list.ind=list.ind+1;
			}
	}	}

}



	#patstresno=as.numeric(pdb$atom[atomsel$atom[patind],"resno"])	#PDB resno of all N, as vector
	#patstins=pdb$atom[atomsel$atom[patind],"insert"]	#PDB resno of all N, as vector

	# If pattern length is 1, then directly return the residue information
	#print(patstresno);
	#patendresno=as.numeric(pdb$atom[atomsel$atom[patind+patlen-1],"resno"])	#PDB resno of all N, as vector
	#patendins=pdb$atom[atomsel$atom[patind+patlen-1],"insert"]	#PDB resno of all N, as vector

	#print(patendresno);
	# If pattern length >1, then we need to check whether st-end=patlen-1 ?
	# inds stores those indexes where st-end is not equal to patlen-1
	# If any such inds found then return the trimmed array else return the whole array
	

	#inds=which((patlen-1)-(patendresno-patstresno)!=0)
	#if(length(inds)>0)
	#{
	#	if(length(inds)==length(patstresno))
	#	{
	#	return(NA);
	#	}else{
	#	return(list(patst=patstresno[-inds],patlen=patlen[-inds]));
	#	}; 
	#}else{
	#return(list(patst=patstresno,patlen=patlen));
	#}
if(length(tempList)>0){return(tempList);}else{return(NA);}
	
}



		
