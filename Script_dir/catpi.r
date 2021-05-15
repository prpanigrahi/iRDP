# catpi
catpi=function(pdb="",catpicut=6.0)
{
	catpimat="SORRY";
	# All atom index LYS/ARG residues
	cation.ind=atom.select(pdb,resid=c("LYS","ARG"),verbose=FALSE)$atom;
	if(! length(cation.ind)>1){return(catpimat);} # If no RK then no catpi

	# All atom index FWY residues
	aro.ind=atom.select(pdb,resid=c("PHE","TYR","TRP"),verbose=FALSE)$atom;
	if(! length(aro.ind)>1){return(catpimat);}
	
	# cation side chain details
	cation.sideatoms=setdiff(unique(pdb$atom[cation.ind,"elety"]),c("N","CA","C","O")); # Sidechain atomnames of these resdiues
	cationside.ind=atom.select(pdb,resid=c("LYS","ARG"),elety=cation.sideatoms,verbose=FALSE)$atom;
	
	# aro side chain details
	aro.sideatoms=setdiff(unique(pdb$atom[aro.ind,"elety"]),c("N","CA","C","O")); # Sidechain atomnames of these resdiues
	aroside.ind=atom.select(pdb,resid=c("PHE","TYR","TRP"),elety=aro.sideatoms,verbose=FALSE)$atom;
	
	# If n cation atoms  and m aro atoms
	# str(acidXYZmat) will be n*3 3:xyz and basicXYZmat=m*3
	cationXYZmat=matrix(as.numeric(pdb$atom[cationside.ind,c("x","y","z")]),ncol=3)
	aroXYZmat=matrix(as.numeric(pdb$atom[aroside.ind,c("x","y","z")]),ncol=3)
	# We calculate distance between each row of cationXYZmat with each row of aroXYZmat
	# So size of distmat will be n*m str(distmat)= length(cationside.in) * length(aroside.int)
	distmat=dist.xyz(cationXYZmat,aroXYZmat)
	arrind=which(distmat<catpicut,arr.ind=TRUE);
	
	# If no interaction is within cutoff then arrind will be empty
	if(nrow(arrind)==0){return(catpimat);}
	
	# Select only those index where rowno < colno i.e. upper triangle
	arrind=matrix(arrind[which(arrind[,1]<arrind[,2]),],ncol=2)
	if(nrow(arrind)==0){return(catpimat);}
	
	# Fetch all residue information of cationic atom using arrind[,1]
	catpimat=matrix(pdb$atom[cationside.ind[arrind[,1]],c("chain","resno","resid","insert")],ncol=4)
	temp=matrix(pdb$atom[aroside.ind[arrind[,2]],c("chain","resno","resid","insert")],ncol=4);
	catpimat=cbind(catpimat,temp);
	
	catpimat[which(is.na(catpimat[,4])),4]="";
	catpimat[,2]=paste(catpimat[,2],catpimat[,4],sep="");

	catpimat[which(is.na(catpimat[,8])),8]="";
	catpimat[,6]=paste(catpimat[,6],catpimat[,8],sep="");
	
	catpimat=catpimat[,c(1:3,5:7)]
return(catpimat);
}

catpi_dupRemove=function(catpimat="")
{
	uniq=unique(paste(catpimat[,1],catpimat[,2],catpimat[,3],catpimat[,4],catpimat[,5],catpimat[,6],sep="_"))
	catpimat=matrix(unlist(strsplit(uniq,"_")),ncol=6,byrow=TRUE);
	colnames(catpimat)=c("chain1","resno1","resid1","chain2","resno2","resid2");
	return(catpimat);
}

catpi_Nochain=function(catpimat,chain1index,chain2index)
{
# chain1index =1, chain2index=5
catpimat[is.na(catpimat[,chain1index]),chain1index]="~";
catpimat[is.na(catpimat[,chain2index]),chain2index]="~";
return(catpimat);
}


catpi_addDSSP=function(DSSP="",catpimat="")
{
catpimat[which(catpimat[,"chain1"]=="~"),"chain1"]=" ";
catpimat[which(catpimat[,"chain2"]=="~"),"chain2"]=" ";
# acidChres corresponds to cationic residue where as baseChres: aromatic
acidChRes=paste(catpimat[,"chain1"],catpimat[,"resno1"],sep="")
baseChRes=paste(catpimat[,"chain2"],catpimat[,"resno2"],sep="")
DSSPChRes=paste(DSSP$cha,DSSP$res,sep="") #combine chain and resno of DSSP object
acidSS="";baseSS="";
	for(i in 1:length(acidChRes))
	{
		index=which(DSSPChRes==acidChRes[i]);
		if(length(index)==0){acidSS=cbind(acidSS,"-");}else{
		acidSS=cbind(acidSS,DSSP$ss[index]);
		}
		index=which(DSSPChRes==baseChRes[i]);
		if(length(index)==0){baseSS=cbind(baseSS,"-");}else{
		baseSS=cbind(baseSS,DSSP$ss[index]);
		}
	}
	acidSS=acidSS[-1];
	baseSS=baseSS[-1];
	#Convert back from " " to "~" chain
	catpimat[which(catpimat[,"chain1"]==" "),"chain1"]="~";
	catpimat[which(catpimat[,"chain2"]==" "),"chain2"]="~";
	catpimat=cbind(catpimat,resid1_SS=acidSS);
	catpimat=cbind(catpimat,resid2_SS=baseSS);
return(catpimat);	
}

catpi_addNaccess=function(nacc="",catpimat="")
{
	chres=paste(nacc$asa[,1],nacc$asa[,2],nacc$asa[,3],sep=""); #ChainResnoResid format
	# we have called the function ionic_Nochain which convert NA to ~
	# But dssp and naccess doesnot detect ~, first we need to convert to " ", then after adding ACC, reconvert back to ~
	catpimat[which(catpimat[,"chain1"]=="~"),"chain1"]=" ";
	catpimat[which(catpimat[,"chain2"]=="~"),"chain2"]=" ";
	acidChRes=paste(catpimat[,"chain1"],catpimat[,"resno1"],catpimat[,"resid1"],sep="")
	baseChRes=paste(catpimat[,"chain2"],catpimat[,"resno2"],catpimat[,"resid2"],sep="")
	acidAcc=0;baseAcc=0;
	for(i in 1:length(acidChRes))
	{
		index=which(chres==acidChRes[i]);
		if(length(index)==0){acidAcc=cbind(acidAcc,"-");}else{
		acidAcc=cbind(acidAcc,nacc$asa[index,5]); #asa_r_rel ASA residue relative
		}
		index=which(chres==baseChRes[i]);
		if(length(index)==0){baseAcc=cbind(baseAcc,"-");}else{
		baseAcc=cbind(baseAcc,nacc$asa[index,5]);
		}
	}
	acidAcc=acidAcc[-1];
	baseAcc=baseAcc[-1];
	#Convert back from " " to "~" chain
	catpimat[which(catpimat[,"chain1"]==" "),"chain1"]="~";
	catpimat[which(catpimat[,"chain2"]==" "),"chain2"]="~";
	catpimat=cbind(catpimat,resid1_Acc=acidAcc);
	catpimat=cbind(catpimat,resid2_Acc=baseAcc);
return(catpimat);	
}

catpi_network=function(catpimat)
{
	#acidic means cation, basic means aromatic
	acidic=paste(catpimat[,"chain1"],catpimat[,"resno1"],catpimat[,"resid1"],sep="");
	basic=paste(catpimat[,"chain2"],catpimat[,"resno2"],catpimat[,"resid2"],sep="");
	if(length(acidic)!=length(basic)){stop("No of acidic and basic residues are not equal,check function ionic_network: length(acidic)!=length(basic)");}	
	total=unique(c(acidic,basic));
	adjmat=matrix(0,nrow=length(total),ncol=length(total));
	# Loop over each catpimat
	for(i in 1:length(acidic))
	{
	adjmat[which(total==acidic[i]),which(total==basic[i])]=1;
	adjmat[which(total==basic[i]),which(total==acidic[i])]=1;
	}
	rownames(adjmat)=total;
	#colnames(adjmat)=total;
	
	library("igraph");
	gcatpimat=graph.adjacency(adjmat,mode="upper",add.rownames=NULL); #adjancy matyrix to graph, upper triangle only considered
	dcatpimat=decompose.graph(gcatpimat,mode="weak",max.comps=NA) # decomposition of graph into sub graph
	return(dcatpimat);
}

get_catpinetmat=function(catpimat,dcatpimat,acclow=20)
{
	size=sapply(dcatpimat, vcount);
	noint=sapply(dcatpimat, ecount);
	ind=which(size>2);
	catpinetmat="";
	
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
	
	acidic=paste(catpimat[,1],catpimat[,2],catpimat[,3],sep="");
	basic=paste(catpimat[,4],catpimat[,5],catpimat[,6],sep="");
	
	for(tempi in 1:length(ind))
	{
	edgelist=get.edgelist(dcatpimat[[ind[tempi]]]);
		acid_net=unique(edgelist[,1]);
		basic_net=unique(edgelist[,2]);
		acc=0;
		for(tacid in 1:length(acid_net))
		{
			tacid.ind=which(acidic==acid_net[tacid])[1]
			acc=c(acc,catpimat[tacid.ind,7]);
		}
		for(tbase in 1:length(basic_net))
		{
			tbase.ind=which(basic==basic_net[tbase])[1]
			acc=c(acc,catpimat[tbase.ind,8]);
		}
		acc=acc[-1];
		acc=as.numeric(acc); # We need to do bcoz acc>30 wont treat differently if acc is string
	#catpinetmat=rbind(catpinetmat,c(size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow & acc<acchigh),sum(acc>=acchigh),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")))
	catpinetmat=rbind(catpinetmat,c(size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")))		
	}
	catpinetmat=catpinetmat[-1,];
	catpinetmat=matrix(catpinetmat,ncol=5);
	colnames(catpinetmat)=c("TotRes","TotInt","BurRes","ExpoRes","Network");
   }
   return(catpinetmat);
}

print_catpinetwork=function(catpimat,dcatpimat,outionic,acclow=20)
{
	size=sapply(dcatpimat, vcount);
	noint=sapply(dcatpimat, ecount);
	ind=which(size>2);
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
	
	acidic=paste(catpimat[,1],catpimat[,2],catpimat[,3],sep="");
	basic=paste(catpimat[,4],catpimat[,5],catpimat[,6],sep="");
	cat("No\tN.Res\tN.Int\tBuried\tExposed\tNetwork_Details\n-----------------------------------",file=outionic,sep="",eol="\n",append=TRUE);
	
	for(tempi in 1:length(ind))
	{
	edgelist=get.edgelist(dcatpimat[[ind[tempi]]]);
		acid_net=unique(edgelist[,1]);
		basic_net=unique(edgelist[,2]);
		acc=0;
		for(tacid in 1:length(acid_net))
		{
			tacid.ind=which(acidic==acid_net[tacid])[1]
			acc=c(acc,catpimat[tacid.ind,7]);
		}
		for(tbase in 1:length(basic_net))
		{
			tbase.ind=which(basic==basic_net[tbase])[1]
			acc=c(acc,catpimat[tbase.ind,8]);
		}
		acc=acc[-1];
		acc=as.numeric(acc); # We need to do bcoz acc>30 wont treat differently if acc is string
		
	#cat(c(tempi,size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow & acc<acchigh),sum(acc>=acchigh),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")),file=outionic,sep="\t",eol="\n",append=TRUE)
	cat(c(tempi,size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")),file=outionic,sep="\t",eol="\n",append=TRUE)
	}
   }
}

getcatpi=function(chres="",catpimat)
{
	chain=chres[1];
	resno=chres[2];
	
	if(is.na(chain)){chain="~"}
	query=paste(chain,resno,sep="",collapse="");
	ip_res1=paste(catpimat[,1],catpimat[,2],sep="");
	ip_res2=paste(catpimat[,4],catpimat[,5],sep="");
	ip.ind=which(ip_res1==query | ip_res2==query);
	if(length(ip.ind)==0){return("SORRY");} # If not found return SORRY else an matrix
	ipmat=catpimat[ip.ind,];
	ipmat=matrix(ipmat,ncol=10);
	colnames(ipmat)=colnames(catpimat);
	return(ipmat)
}

get_catpinetwork=function(chres="",catpimat,dcatpimat,acclow=20)
{
	chain=chres[1];
	resno=chres[2];
	query=paste(chain,resno,sep="",collapse="");
	size=sapply(dcatpimat, vcount);
	noint=sapply(dcatpimat, ecount);
	ind=which(size>2);
	Net="SORRY";
	
	# If network (size is 3 or more) present then only proceed
	if(length(ind)>0)
	{
	size3=size[ind];
	noint3=noint[ind];
	# order them
	order=order(size3)
	size3=size3[order]
	noint3=noint3[order]
	ind=ind[order]
	
	acidic=paste(catpimat[,1],catpimat[,2],catpimat[,3],sep="");
	basic=paste(catpimat[,4],catpimat[,5],catpimat[,6],sep="");
	
	for(tempi in 1:length(ind))
	{
	edgelist=get.edgelist(dcatpimat[[ind[tempi]]]); 
	#	[,1]      [,2]     
	#[1,] "~341GLU" "~344LYS"
	#[2,] "~341GLU" "~411ARG"
	edgeChres=gsub("...$","",as.vector(edgelist));  #"~341" "~341" "~410" "~344" "~411" "~411"  

		if(length(grep(paste(query,"$",sep="",collapse=""),edgeChres,perl=TRUE))>0)
		{
			acid_net=unique(edgelist[,1]);
			basic_net=unique(edgelist[,2]);
			acc=0;
			for(tacid in 1:length(acid_net))
			{
				tacid.ind=which(acidic==acid_net[tacid])[1]
				acc=c(acc,catpimat[tacid.ind,7]);
			}
			for(tbase in 1:length(basic_net))
			{
				tbase.ind=which(basic==basic_net[tbase])[1]
				acc=c(acc,catpimat[tbase.ind,8]);
			}
			acc=acc[-1];
			acc=as.numeric(acc); # We need to do bcoz acc>30 wont treat differently if acc is string
			Net=rbind(Net,c(tempi,size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")));
		}
	}
	# If in any single network query is found then Net will be a matrix and first row is "" "" ""
	if(is.matrix(Net)){Net=Net[-1,];Net=matrix(Net,ncol=6);}
	return(Net);
	}else{
	return(Net); # Returns "";
	}
}

callCatpi=function(pdb,DSSP="",nacc="",status.DSSP=0,status.nacc=0,catpicut)
{
	# This is the workflow that is used in Pattern_analysis code
	# catpimat
	catpimat=catpi(pdb,catpicut=catpicut);
	dcatpimat="";
	status.catpi=0;
	if(is.matrix(catpimat))
	{
	status.catpi=1;
	catpimat=Nochain(catpimat,1,4);
	catpimat=catpi_dupRemove(catpimat);
	if(status.nacc){catpimat=catpi_addNaccess(nacc,catpimat);}else{catpimat=cbind(catpimat,resid1_Acc="-");catpimat=cbind(catpimat,resid2_Acc="-");}
	if(status.DSSP){catpimat=catpi_addDSSP(DSSP,catpimat);}else{catpimat=cbind(catpimat,resid1_SS="-");catpimat=cbind(catpimat,resid2_SS="-");}
	dcatpimat=catpi_network(catpimat);
	}
return(list(catpimat=catpimat,status.catpi=status.catpi,dcatpimat=dcatpimat));	
}

checkCatpi=function(patmat,catpimat)
{
	# This is the workflow that is used in Pattern_analysis code
	# refer to hbond.r checkHB function for details
	ionmat=do.call(rbind,lapply(apply(patmat,1,getcatpi,catpimat),matrix,ncol=10,byrow=FALSE));
	sorry.ind=which(apply(ionmat,1,function(x){sum(x=="SORRY")>0}));
	if(length(sorry.ind)>0){ionmat=ionmat[-sorry.ind,];}
	ionmat=matrix(ionmat,ncol=10);
	# For 3 residue pattern, if all 3 residue doesnot form any inoic interaction, we get 3*10 matrix with all value sorry
	# And after ionmat[-,] i.e. removing those SORRY row we get an matrix but nrow(ionmat)=0
	return(ionmat);				
}

checkcatpinet=function(patmat,catpimat,dcatpimat,acclow=20)
{
	# This is the workflow that is used in Pattern_analysis code
	#check for network
	ionnetmat=do.call(rbind,lapply(apply(patmat,1,get_catpinetwork,catpimat,dcatpimat,acclow),matrix,ncol=6,byrow=FALSE));
	sorry.ind=which(apply(ionnetmat,1,function(x){sum(x=="SORRY")>0}));
	if(length(sorry.ind)>0){ionnetmat=ionnetmat[-sorry.ind,];}
	ionnetmat=matrix(ionnetmat,ncol=6)
	ionnetmat=matrix(ionnetmat[,2:6],ncol=5)
	if(nrow(ionnetmat)>0)
	{
	ionnetmat=matrix(unlist(strsplit(unique(apply(ionnetmat,1,function(x){paste(x,sep="",collapse=":")})),":")),ncol=5,byrow=TRUE) #remove duplicate
	}
return(ionnetmat);
}

catpiSummary=function(catpimat,acclow=20)
{
	ionsum=c(total=length(catpimat[,1]));
	
	#ionsum=c(ionsum,pIntra=round((sum(catpimat[,1]==catpimat[,4])*100)/length(catpimat[,1]),2));
	#ionsum=c(ionsum,pInter=round((sum(catpimat[,1]!=catpimat[,4])*100)/length(catpimat[,1]),2));	

	ionsum=c(ionsum,pIntra=sum(catpimat[,1]==catpimat[,4]));
	ionsum=c(ionsum,pInter=sum(catpimat[,1]!=catpimat[,4]));	

	# Percentage of catpimat involving D/E/R/K/H
	resid1=aa321(catpimat[,3]) #resid1
	resid2=aa321(catpimat[,6]) #resid2
	
	# Percentage of ionpair of type RF,RY,RW,KF,KY,KW
	iontype=paste(resid1,resid2,sep="")
	
	#ionsum=c(ionsum,pKF=round((sum(iontype=="KF")*100)/length(iontype),2));
	#ionsum=c(ionsum,pKY=round((sum(iontype=="KY")*100)/length(iontype),2));
	#ionsum=c(ionsum,pKW=round((sum(iontype=="KW")*100)/length(iontype),2));
	#ionsum=c(ionsum,pRF=round((sum(iontype=="RF")*100)/length(iontype),2));
	#ionsum=c(ionsum,pRY=round((sum(iontype=="RY")*100)/length(iontype),2));
	#ionsum=c(ionsum,pRW=round((sum(iontype=="RW")*100)/length(iontype),2));
	
	ionsum=c(ionsum,pKF=sum(iontype=="KF"));
	ionsum=c(ionsum,pKY=sum(iontype=="KY"));
	ionsum=c(ionsum,pKW=sum(iontype=="KW"));
	ionsum=c(ionsum,pRF=sum(iontype=="RF"));
	ionsum=c(ionsum,pRY=sum(iontype=="RY"));
	ionsum=c(ionsum,pRW=sum(iontype=="RW"));
	

	
	resid1_acc=as.numeric(catpimat[,"resid1_Acc"])
	resid2_acc=as.numeric(catpimat[,"resid2_Acc"])
	avg_Acc=(resid1_acc+resid2_acc)/2;
	ionsum=c(ionsum,pB=round((sum(avg_Acc<=acclow)*100)/length(avg_Acc),2));
	ionsum=c(ionsum,pE=round((sum(avg_Acc>acclow)*100)/length(avg_Acc),2));
		
	# catpimat Network
	dcatpimat=catpi_network(catpimat); #decomposed catpimat graph
	network=paste(sapply(dcatpimat, vcount),sapply(dcatpimat,ecount),sep="_");
	ionsum=c(ionsum,pIso=round((sum(network=="2_1")*100)/length(catpimat[,1]),2));
	netfreq=tapply(network,factor(network),length)		
	ionsum=c(ionsum,totNet=sum(sapply(dcatpimat, vcount)>2));
	netfreq=netfreq[names(netfreq)!="2_1"];
	NetDet="-";
	if(length(netfreq)>0){NetDet=paste(names(netfreq),netfreq,sep=":",collapse=", ")}
	ionsum=c(ionsum,NetDet=NetDet);
	return(ionsum);
}

print_HTMLcatpi=function(catpimat,catpinetmat,filepath,append=FALSE,html=TRUE)
{
	if(html){write2Html("<html>\n",filepath,append);}else{write2Html("",filepath,append);}
	write2Html("<br><br><br><br><b><u>Cation-pi interactions (CPI)</u><br><br></b>\n",filepath,TRUE);
#	write2Html("Chain1,Resno1,Resid1: Corresponds to Acidic residue details<br>\n",filepath,TRUE);
#	write2Html("Chain2,Resno2,Resid2: Corresponds to Basic residue details<br>\n",filepath,TRUE);
#	write2Html("Resid1_Acc and Resid2_Acc corresponds to relative ASA of acidic and basic residue respectively<br>\n",filepath,TRUE);
#	write2Html("Resid1_SS and Resid2_SS corresponds to secondary structure of the acidic and basic residues respectively<br><br><br>\n",filepath,TRUE);

	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr ><th colspan=3 title=\"Cationic residue [ KR ]\"> Cationic residue </th> <th colspan=3 title=\"Aromatic residue [ FYW ]\"> Aromatic residue </th> <th colspan=2 title=\"Relative solvent accessibility of the residue\"> ASA </th> <th colspan=2 title=\"Secondary structure of residue (DSSP notation)\" > SS </th> </tr>",filepath,append=TRUE);
	tipvec=c("Chain","Residue number","Three letter code of residue","Chain","Residue number","Three letter code of residue","","","","");
	write2Html.tableHeader(c("Chain","Res.No","Res.ID","Chain","Res.No","Res.ID","Cation","Aromatic","Cation","Aromatic"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(catpimat,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	
	if(is.matrix(catpinetmat) && nrow(catpinetmat) > 0)
	{
	write2Html("<br><br><br><b><u>Cation-pi interaction (CPI)  network</u></b><br><br>\n",filepath,TRUE);
	#write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; width:1000px\">\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed;\">\n",filepath,TRUE);
	tipvec=c("Total number of residues present in the network","Total number of interactions present in the network","Number of residues of the network that are buried",
	"Number of residues of the network that are exposed","Format: Partner1 (Chain-Resno-Resid) -- Partner2 (Chain-Resno-Resid)"
	);
	
	write2Html.tableHeader(c("N.Res","N.Int","Buried","Exposed","Network_Details"),filepath,TRUE,tipvec);
	#write2Html.tableBody(ionnetmat,filepath,TRUE,"<td style=\"word-wrap: break-word\">");
	catpinetmat[,5]=apply(matrix(catpinetmat[,5],ncol=1),1,function(x){gsub(",","<br>",x)})
	write2Html.tableBody(catpinetmat,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	if(html){write2Html("</html>\n",filepath,TRUE);}else{write2Html("",filepath,TRUE);}
	}else{
	#write2Html("<br><br><br> <b>No Cation-pi interaction network is observed</b>\n",filepath,TRUE);
	if(html){write2Html("</html>\n",filepath,TRUE);}else{write2Html("",filepath,TRUE);}
	}
}

print_HTMLcatpi_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Cation-pi interaction (CPI) summary</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("Filename","Total number of CPIs","Number of intra-chain CPIs","Number of inter-chain CPIs",
	paste("Number of CPIs of type [Cation-pi]: ",c("KF","KY","KW","RF","RY","RW"),sep=""),
	"Percentage of buried CPIs","Percentage of exposed CPIs","Percentage of CPIs that are not part of any network (Isolated)",
	"Total number of CPI-networks","Network details [SIZE_STRENGTH:NUMBER OF SUCH NETWORK] (SIZE:Number of residue in the network, STRENGTH:Number of interaction in the network)"
	);
	write2Html.tableHeader(c("Filename","Total","IntraCh","InterCh","KF","KY","KW","RF","RY","RW","B","E","Iso","T.Net","Net_Details"),filepath,append=TRUE,tipvec);
}

