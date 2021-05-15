# modeller.r




modeller=function(pdb,exepath="",prefix,mutstr)
{
  #tryCatch({
    
	mutlist=unlist(strsplit(mutstr,",")) # Split with comma the input mutation string
	mutlist=gsub(";","",mutlist); #remove any trainling ; from each mutation
	mutlist=mutlist[grep("^[A-Z][A-Z0-9]\\d+[A-Z]$",mutlist)]; #Only Strict format will be accepted and Format: [A-Z] [A-Z] at least one digit and [A-Z]
  if(length(mutlist)==0){cat(paste(c(" (Some error while handeling mutation)","\n\n"),sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);return(FALSE);}
  mut_patmat=cbind(substr(mutlist,2,2),substr(mutlist,3,nchar(mutlist)-1),substr(mutlist,nchar(mutlist),nchar(mutlist)),substr(mutlist,1,1)); # [chain,resno,mutres,wildres]

	print(mut_patmat);
	mutfilename=paste(prefix,".pdb",sep="",collapse="");
  print(mutfilename);
	write.pdb(pdb,file=mutfilename); # write wild type pdb file with prefix name, it will be 1abc_modeller_mut1.pdb for 1st mutant(single/multi)
	mut_patmat[which(mut_patmat[,1]=="_"),1]="";
	
	flag=0;
	# Loop over each mutation of mutant list.
	for(mut in 1:nrow(mut_patmat))
	{
		
		chres.ind=atom.select(pdb,chain=mut_patmat[mut,1],resno=mut_patmat[mut,2],verbose=FALSE)$atom; # Obtain the pdb atom index using chain and resno information
		# if either chain is not present in pdb (e.g. user input as chain X but in pdb chain A is there)
		# or if chain information is correct but that residue is missing then chres.in will be zero
		# In that case we skip to next mutation
		
		if(length(chres.ind)==0)
		{
		cat(paste(c(" Either the input chain or residue number is missing in the PDB file, so skipping...","\n\n"),sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
		flag=1;
		next;
		}
		
		actres=aa321(unique(pdb$atom[chres.ind,"resid"])); #actual residue name

		if(actres!=mut_patmat[mut,4])
		{
			cat(paste(c("\n\n",mut_patmat[mut,4],mut_patmat[mut,1],mut_patmat[mut,2],mut_patmat[mut,3],":wild type residue is mismatching. True residue is:",actres," Skipping\n"),sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
			#mut_patmat[mut,4]=actres; #replace the wrong wildres to actual residue
			flag=1;
			next;
		}
		
		todomut=paste("python ",exepath," ",prefix," ",mut_patmat[mut,2]," ",aa123(mut_patmat[mut,3])," ",mut_patmat[mut,1]);
		print(todomut);
		system(todomut,ignore.stderr=FALSE,ignore.stdout=TRUE);
		modelpref=paste(prefix,aa123(mut_patmat[mut,3]),mut_patmat[mut,2],mut_patmat[mut,1],sep="",collapse="");
		print(modelpref);
		unlink(paste(prefix,".pdb",sep="",collapse=""));
		file.rename(from=paste(modelpref,".pdb",sep="",collapse=""),to=paste(prefix,".pdb",sep="",collapse=""));
	}
  #},error=function(e){cat(paste(c(" (Some error while handeling mutation)","\n\n"),sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE); unlink(mutfilename);return(FALSE);},finally=print("ERROR"));
	
	if(flag==1)
	{
	unlink(mutfilename);
	return(FALSE);
	}
	
	if(!file.exists(paste(prefix,".pdb",sep="",collapse="")))
	{
	cat(paste(c(" (Invalid mutation)","\n\n"),sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
	return(FALSE);
	}else{
	cat(paste(c(" (Valid mutation)","\n\n"),sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
	return(TRUE);
	}
}

print_HTMLionic_mutation_start=function(filepath)
{
	write2Html("<h3>Mutation summary</h3>\n",filepath,TRUE);
	write2Html("The result table below describes the possible loss/gain of interactions due to mutation.
The loss or gain of interaction is quantitatively summarized as <b>Interaction profile</b>. 
The values in the interaction profile represents the number of various interactions the wild/mutant residue is involved.<br> 
By comparing the wild and mutant interaction profile one can get quick overview of possible loss/gain of interactions due to mutation.
<br>Click on the link given in <b>Local Interactions</b> column to get a detailed interaction analysis. <br>\nClick on the link given in <b>Download</b> column to download the mutant structures.<br>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr> <th rowspan=2 title=\" Detailed Local interactions of wild/mutant residue\"> Local Interactions </th> <th rowspan=2 title=\" Click to download the mutant PDB file\"> Download </th> <th rowspan=2 title=\" Type of residue [Wild/Mutant]\"> Type </th> <th colspan=4 title=\" The Wild or mutant residue details\"> Residue details </th> <th colspan=11  title=\" The interaction profile of the wild/mutant residue\"> Interaction profile </th> </tr>",filepath,append=TRUE);
	tipvec=c("Chain","Residue number","Residue name","Conservation Score",
	  "Number of ion-pairs formed by wild/mutant residue",
	  "Number of ion-pair networks in which wild/mutant residue is involved in",
	  "Number of aromatic interactions formed by wild/mutant residue",
	  "Number of aromatic-aromatic interaction networks in which wild/mutant residue is involved in",
	  "Number of Aro-S interaction formed by wild/mutant residue",
	  "Number of Aro-S interaction networks in which wild/mutant residue is involved in",
	  "Number of hydrogen bonds formed by wild/mutant residues",
	  "Number of disulfide bonds formed by wild/mutant residues",
	  "Number of Cat-pi interaction formed by wild/mutant residues",
	  "Number of Cat-pi interaction networks in which wild/mutant residue is involved in",
	  "Number of Hydrophobic formed by wild/mutant residues");
	 write2Html.tableHeader(c("Chain","Res.No","Res.ID","CScore","IP","IP.Net","AP","AP.Net","AS","AS.Net","HB","Disul","Cat-pi","Cat-pi.Net","Hphob"),filepath,append=TRUE,tipvec);
}



