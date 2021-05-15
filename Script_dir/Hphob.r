# Hphob.r

hphob=function(pdb="",cutoff=5.0)
{
	hphobmat="";
	# All atom index of hydrophobic residues
	hphobres.ind=atom.select(pdb,resid=c("ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO","TYR"),verbose=FALSE)$atom
	# At least 2 atoms hould be there to calculate distance.
	if(! length(hphobres.ind)>1){return(hphobmat);}
	hphob.sideatoms=setdiff(unique(pdb$atom[hphobres.ind,"elety"]),c("N","CA","C","O")); # Sidechain atomnames of these resdiues
	hphobside.ind=atom.select(pdb,resid=c("ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO","TYR"),elety=hphob.sideatoms,verbose=FALSE)$atom;
	hphob.detail=pdb$atom[hphobside.ind,c("chain","resno","resid","elety","x","y","z","insert")];
	hphob.detail[which(is.na(hphob.detail[,8])),8]="";
	hphob.detail[,2]=paste(hphob.detail[,2],hphob.detail[,8],sep="")

	#calcualte  distance
	hphobdist=dist.xyz(matrix(as.numeric(hphob.detail[,c("x","y","z")]),ncol=3))
	
	dist.Cutoff=matrix(which(hphobdist <cutoff,arr.ind=TRUE),ncol=2);
	# Whatever cutoff u use as distance between same atom is 0, so we will always get a matrix.
	# Select only upper triangle if after selcting nothing pass the filter then nrow(dist.Cutoff)=0
	dist.Cutoff=matrix(dist.Cutoff[which(dist.Cutoff[,1]<dist.Cutoff[,2]),],ncol=2)
	if(nrow(dist.Cutoff)==0){return(hphobmat);}
	
	hphobTemp=matrix(hphob.detail[dist.Cutoff[,1],c("chain","resno","resid")],ncol=3);
	hphobTemp=cbind(hphobTemp,matrix(hphob.detail[dist.Cutoff[,2],c("chain","resno","resid")],ncol=3));
	
	uniq=unique(paste(hphobTemp[,1],hphobTemp[,2],hphobTemp[,3],hphobTemp[,4],hphobTemp[,5],hphobTemp[,6],sep="_"))
	hphobmat=matrix(unlist(strsplit(uniq,"_")),ncol=6,byrow=TRUE);
	# Bcoz we are doing all atom distance, for a given residue if there is 6 sidcechain atoms and if they fall within 5A then we get
	# False positive i.e same residue interacting with same.
	# We neeed to remove those interactions where residue1 is interacting with itself;
	
	resch1=paste(hphobmat[,1],hphobmat[,2],hphobmat[,3],sep="");
	resch2=paste(hphobmat[,4],hphobmat[,5],hphobmat[,6],sep="");
	hphobmat=hphobmat[- which(resch1==resch2),];
	hphobmat=matrix(hphobmat,ncol=6)
	
	colnames(hphobmat)=c("Chain1","Resno1","Resid1","Chain2","Resno2","Resid2");
	return(hphobmat);	
	
}

hphob_addDSSP=function(DSSP="",hp="")
{
#     Chain1 Resno1 Resid1 Chain2 Resno2 Resid2
# [1,] "A"    "5"    "TYR"  "A"    "5"    "TYR" 

	hp[which(hp[,1]=="~"),1]=" ";  # Replace ~ of chain1 with " "
	hp[which(hp[,4]=="~"),4]=" "; # Replace ~ chain2 with " "
	
	hp1ChRes=paste(hp[,1],hp[,2],sep="")
hp2ChRes=paste(hp[,4],hp[,5],sep="")
DSSPChRes=paste(DSSP$cha,DSSP$res,sep="") #combine chain and resno of DSSP object
hp1SS="";hp2SS="";
	for(i in 1:length(hp1ChRes))
	{
		index=which(DSSPChRes==hp1ChRes[i]);
		if(length(index)==0){hp1SS=cbind(hp1SS,"-");}else{
		hp1SS=cbind(hp1SS,DSSP$ss[index]);
		}
		index=which(DSSPChRes==hp2ChRes[i]);
		if(length(index)==0){hp2SS=cbind(hp2SS,"-");}else{
		hp2SS=cbind(hp2SS,DSSP$ss[index]);
		}
		
	}
	hp1SS=hp1SS[-1];
	hp2SS=hp2SS[-1];
	hp=cbind(hp,resid1_SS=hp1SS);
	hp=cbind(hp,resid2_SS=hp2SS);

hp[which(hp[,1]==" "),1]="~";  # Replace ~ of chain1 with " "
hp[which(hp[,4]==" "),4]="~"; # Replace ~ chain2 with " "
return(hp);
}

hphob_addNaccess=function(nacc="",hp="")
{
#chain1     resno1     resid1     chain2     resno2     resid2   CentDist      Dihed
# "A"       "29"      "TYR"        "A"       "26"      "PHE"    "6.429"   "35.859"
	chres=paste(nacc$asa[,1],nacc$asa[,2],nacc$asa[,3],sep=""); #ChainResnoResid format
	hp[which(hp[,1]=="~"),1]=" ";
	hp[which(hp[,4]=="~"),4]=" ";
	hp1ChRes=paste(hp[,1],hp[,2],hp[,3],sep="") #Chain1, Resno1, Resid1
	hp2ChRes=paste(hp[,4],hp[,5],hp[,6],sep="") #Chain2, Resno2, Resid2
	hp1Acc=0;hp2Acc=0;
	for(i in 1:length(hp1ChRes))
	{
		index=which(chres==hp1ChRes[i]);
		if(length(index)==0){hp1Acc=cbind(hp1Acc,"-");}else{
		hp1Acc=cbind(hp1Acc,nacc$asa[index,5]); #asa_r_rel ASA residue relative
		}
		
		index=which(chres==hp2ChRes[i]);
		if(length(index)==0){hp2Acc=cbind(hp2Acc,"-");}else{
		hp2Acc=cbind(hp2Acc,nacc$asa[index,5]);
		}
	}
	hp1Acc=hp1Acc[-1];
	hp2Acc=hp2Acc[-1];
	hp=cbind(hp,resid1_Acc=hp1Acc);
	hp=cbind(hp,resid2_Acc=hp2Acc);
hp[which(hp[,1]==" "),1]="~";  # Replace " " of chain1 with "~"
hp[which(hp[,4]==" "),4]="~"; # Replace " "  chain2 with "~"
return(hp);	
}


hphob_summary=function(hp,acclow=20)
{
	Total=nrow(hp);
	
	#pIntra=round((sum(hp[,1]==hp[,4])/nrow(hp))*100,2);
	#pInter=round((sum(hp[,1]!=hp[,4])/nrow(hp))*100,2);
	
	pIntra=sum(hp[,1]==hp[,4]);
	pInter=sum(hp[,1]!=hp[,4]);
	
	
	resid1_acc=as.numeric(hp[,"resid1_Acc"])
	resid2_acc=as.numeric(hp[,"resid2_Acc"])
	avg_Acc=(resid1_acc+resid2_acc)/2;
	pB=round((sum(avg_Acc<=acclow)*100)/length(avg_Acc),2);
	pE=round((sum(avg_Acc>acclow)*100)/length(avg_Acc),2);
	hpsum=c(Total=Total,pIntra=pIntra,pInter=pInter,pBur=pB,pExpo=pE);
	
	# hppair Network
	#dhppair=hphob_network(hp); #decomposed ionpair graph
	#network=paste(sapply(dhppair, vcount),sapply(dhppair,ecount),sep="_");
	#hpsum=c(hpsum,pIso=round((sum(network=="2_1")*100)/length(hp[,1]),2));
	
	#netfreq=tapply(network,factor(network),length)		
	#hpsum=c(hpsum,totNet=sum(sapply(dhppair, vcount)>2));
	#hpsum=c(hpsum,NetDet=paste(names(netfreq),netfreq,sep=":",collapse=", "));
	
	return(hpsum);
}

hphob_network=function(hp)
{
	hp1=paste(hp[,1],hp[,2],hp[,3],sep=""); #index 1,2,3 correspond to chain1, resno1,resid1
	hp2=paste(hp[,4],hp[,5],hp[,6],sep=""); #index 4,5,6 correspond to chain2, resno2,resid2
	if(length(hp1)!=length(hp2)){stop("No of hp1 and hp2 residues are not equal,check function hphp_network: length(hp1)!=length(hp2)");}	
	# One residue can form interaction more than one residue
	# Hence we get the list of all residues involved in hphp interaction by unique()
	# If n such residues are there then size of adjancy matrix wil;l be n*n
	total=unique(c(hp1,hp2)); 
	adjmat=matrix(0,nrow=length(total),ncol=length(total));
	# Loop over each hppair
	for(i in 1:length(hp1))
	{
	row=which(total==hp1[i]);col=which(total==hp2[i]);
	adjmat[row,col]=adjmat[row,col]+1;
	adjmat[col,row]=adjmat[col,row]+1;
	}
	rownames(adjmat)=total;
	#colnames(adjmat)=total;
	
	library("igraph");
	ghppair=graph.adjacency(adjmat,mode="upper",add.rownames=NULL); #adjancy matyrix to graph, upper triangle only considered
	dhppair=decompose.graph(ghppair,mode="weak",max.comps=NA) # decomposition of graph into sub graph
	return(dhppair);
}
print_hpobnetwork=function(hp,dhppair,outhphob)
{
	size=sapply(dhppair, vcount)
	noint=sapply(dhppair, ecount);
	ind=which(size>0);
	# If network present then only proceed
	if(length(ind)>0)
	{
	size3=size[ind];
	noint3=noint[ind];
	# order them
	order=order(size3)
	size3=size3[order]
	noint3=noint3[order]
	ind=ind[order]

	hpres=c(paste(hp[,1],hp[,2],hp[,3],sep=""),paste(hp[,4],hp[,5],hp[,6],sep=""))
	hpacc=c(hp[,9],hp[,10])
	cat("No\tN.Res\tN.Int\tBuried\tPartial\tExposed\tNetwork_Details\n-----------------------------------",file=outhphob,sep="",eol="\n",append=TRUE);

	for(tempi in 1:length(ind))
	{
	edgelist=get.edgelist(dhppair[[ind[tempi]]]);
	hpres_net=unique(c(edgelist[,1],edgelist[,2]));
		acc=0;
		for(thp in 1:length(hpres_net))
		{
			thp.ind=which(hpres==hpres_net[thp])[1]
			acc=c(acc,hpacc[thp.ind]);
		}
		acc=acc[-1];
		acc=as.numeric(acc); # We need to do bcoz acc>30 wont treat differently if acc is string
	cat(c(tempi,size3[tempi],noint3[tempi],sum(acc<=30),sum(acc>30 & acc<70),sum(acc>=70),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")),file=outhphob,sep="\t",eol="\n",append=TRUE)	
	}
  }
}


getHphob=function(chres="",hp="")
{
	chain=chres[1];
	resno=chres[2];
	if(is.na(chain)){chain="~"}
	query=paste(chain,resno,sep="",collapse="");
	hp_res1=paste(hp[,1],hp[,2],sep="");
	hp_res2=paste(hp[,4],hp[,5],sep="");
	hp.ind=which(hp_res1==query | hp_res2==query);
	if(length(hp.ind)==0){return("SORRY");} # If not found return SORRY else an matrix
	hpmat=hp[hp.ind,];
	hpmat=matrix(hpmat,ncol=10);
	colnames(hpmat)=colnames(hp);
	return(hpmat)
}

callHphob=function(pdb="",DSSP="",nacc="",status.DSSP=0,status.nacc=0,hpcut)
{
			# This is the workflow that is used in Pattern_analysis code
			hp=hphob(pdb,cutoff=hpcut);
			status.hp=0;
			if(is.matrix(hp) && nrow(hp)>0)
			{
			status.hp=1;
			hp=Nochain(hp,1,4); #1st col is chain1 and 2nd col is 
			if(status.DSSP){hp=hphob_addDSSP(DSSP,hp);}else{hp=cbind(hp,resid1_SS="~");hp=cbind(hp,resid2_SS="~");}
			if(status.nacc){hp=hphob_addNaccess(nacc,hp);}else{hp=cbind(hp,resid1_Acc="~");hp=cbind(hp,resid2_Acc="~");}
			}
return(list(hp=hp,status.hp=status.hp));			
}
checkHphob=function(patmat,hp)
{
	# This is the workflow that is used in Pattern_analysis code
	hpmat=do.call(rbind,lapply(apply(patmat,1,getHphob,hp),matrix,ncol=10,byrow=FALSE));
	sorry.ind=which(apply(hpmat,1,function(x){sum(x=="SORRY")>0}));
	if(length(sorry.ind)>0){hpmat=hpmat[-sorry.ind,];}
	hpmat=matrix(hpmat,ncol=10);
	return(hpmat);					
}

print_HTMLhphob=function(hp,filepath,append=FALSE,html=TRUE)
{
	if(html){write2Html("<html>\n",filepath,append);}else{write2Html("",filepath,append);}
	write2Html("<br><br><b><u>Hydrophobic interactions</u><br><br></b>\n",filepath,TRUE);
#	write2Html("Chain1,Resno1,Resid1: Corresponds to Acidic residue details<br>\n",filepath,TRUE);
#	write2Html("Chain2,Resno2,Resid2: Corresponds to Basic residue details<br>\n",filepath,TRUE);
#	write2Html("Resid1_Acc and Resid2_Acc corresponds to relative ASA of acidic and basic residue respectively<br>\n",filepath,TRUE);
#	write2Html("Resid1_SS and Resid2_SS corresponds to secondary structure of the acidic and basic residues respectively<br><br><br>\n",filepath,TRUE);

	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr ><th colspan=3 title=\"First hydrophobic residue partner\"> Residue1 </th> <th colspan=3 title=\"Second hydrophobic residue partner\"> Residue2 </th> <th colspan=2 title=\"Secondary structure of residue (DSSP notation)\"> SS </th> <th colspan=2 title=\"Relative solvent accessibility of the residue\"> ASA </th> </tr>",filepath,append=TRUE);
	tipvec=c("Chain","Residue number","Three letter code of residue","Chain","Residue number","Three letter code of residue","","","","");
	write2Html.tableHeader(c("Chain","Res.No","Res.ID","Chain","Res.No","Res.ID","Res1","Res2","Res1","Res2"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(hp,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	if(html){write2Html("</html>\n",filepath,TRUE);}else{write2Html("",filepath,TRUE);}
}

print_HTMLhphob_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Hydrophobic interaction (HPI) summary</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("PDB Filename","Total number of HPI","Number of intra-chain HPI","Number of inter-chain HPI",
	"Percentage of buried HPI","Percentage of exposed HPI"
	);
	write2Html.tableHeader(c("Filename","Total","IntraCh","InterCh","B","E"),filepath,append=TRUE,tipvec);
}

