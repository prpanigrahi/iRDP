aroaro=function(pdb="",low=4.5,high=7.0,angle=30.0,usedihed=1)
{
	phe=aroCentroid(pdb,restype="PHE");
	tyr=aroCentroid(pdb,restype="TYR");
	trp=aroCentroid(pdb,restype="TRP");
	FYW=matrix(0,ncol=9);
	# If no phe/y/w found aroCentroid function returns "SORRY" . If one F/Y/W present, it returns matrix, if more than one then also returns matrix
	# If returns a matrix then only u do rbind or bind the centroid coordinates
	if(is.matrix(phe)){FYW=rbind(FYW,phe);}
	if(is.matrix(tyr)){FYW=rbind(FYW,tyr);}
	if(is.matrix(trp)){FYW=rbind(FYW,trp);}
	
	FYW=FYW[-1,]
	FYW=matrix(FYW,ncol=9);
	colnames(FYW)=c("chain","resno","resid","Cx","Cy","Cz","EqnA","EqnB","EqnC");
	
	# If no phe/y/w then no interaction so return
	if(nrow(FYW)==0){return("SORRY");}
	
	#calcualte centroid distance
	FYWdist=dist.xyz(matrix(as.numeric(FYW[,c("Cx","Cy","Cz")]),ncol=3))
	dist.Cutoff=matrix(which(FYWdist>=low & FYWdist <=high,arr.ind=TRUE),ncol=2);
	# FYWdist is a matrix of size n*3 n:no of aromatic res, 3:xyz
	# dist.xyz(FYWdist) gives n*n matrix with distances between centroids
	# We should only consider upper triangle
	# Select only those index where rowno < colno 2,29
	dist.Cutoff=matrix(dist.Cutoff[which(dist.Cutoff[,1]<dist.Cutoff[,2]),],ncol=2)
	# If no residues are falling within cutoff then dist.Cutoff gives no row and column. so we proceed 
	# further for dihedrals or returning if any residue are withing cutoff. Else no interaction is seen.
	if(dim(dist.Cutoff)[1] > 0)
	{
	# dist.Cutoff is n*2 array with n=no of aromatic pairs that satisfies the condition
	# 2: row and column index (Row indicates one aro and col indicates index of second aro). So dist.Cutoff[,1] will give all row index of FYWdist matrix which is same as FYW matrix
	# Using the row index (dist.Cutoff[,1]) we fetch "chain","resno","resid" of first partner
	# Using the row index (dist.Cutoff[,2]) we fetch "chain","resno","resid" of 2nd partner
		aroaroTemp=matrix(FYW[dist.Cutoff[,1],c("chain","resno","resid")],ncol=3);
		aroaroTemp=cbind(aroaroTemp,matrix(FYW[dist.Cutoff[,2],c("chain","resno","resid")],ncol=3));
		aroaroTemp=cbind(aroaroTemp,round(as.numeric(FYWdist[dist.Cutoff]),2));
		#calculate dihedrals
		dihed=0;
		for(i in 1:length(dist.Cutoff[,1]))
		{
			tempEQN1=FYW[dist.Cutoff[i,1],c("EqnA","EqnB","EqnC")];
			tempEQN2=FYW[dist.Cutoff[i,2],c("EqnA","EqnB","EqnC")];
			dihed=rbind(dihed,ringANGL(tempEQN1,tempEQN2));	
		}
		dihed=dihed[-1,]
		acute=which(dihed>90.0);
		dihed[acute]=180-dihed[acute]; # convert all obtuse angle to acute
		dihed=round(dihed,2);
		aroaroTemp=cbind(aroaroTemp,dihed); 
		# use dihed cutoff or not? Default is to use
		if(usedihed){aroaroTemp=aroaroTemp[which(as.numeric(aroaroTemp[,8])>= angle),]}
		
		# if after putting dihedral cutoff, no aroaro pair remains then we return "sorry"
		tempname=c("chain1","resno1","resid1","chain2","resno2","resid2","CentDist","Dihed");
		if(is.matrix(aroaroTemp)& length(aroaroTemp)> 0 ){colnames(aroaroTemp)=tempname;rownames(aroaroTemp)=1:length(aroaroTemp[,1]);return(aroaroTemp);}
		else
		{
		#aroaroTemp=matrix(aroaroTemp,nrow=1);colnames(aroaroTemp)=tempname;
		return("SORRY");
		}
		
	}else{
	return("SORRY");
	}
}
aroaro_addDSSP=function(DSSP="",aro="")
{
#colnames of aro: chain1   resno1   resid1   chain2   resno2   resid2 CentDist    Dihed
# "A"     "29"    "TYR"      "A"     "26"    "PHE"  "6.429" "35.859"

aro[which(aro[,1]=="~"),1]=" ";  # Replace ~ of chain1 with " "
aro[which(aro[,4]=="~"),4]=" "; # Replace ~ chain2 with " "
	
aro1ChRes=paste(aro[,1],aro[,2],sep="")
aro2ChRes=paste(aro[,4],aro[,5],sep="")
DSSPChRes=paste(DSSP$cha,DSSP$res,sep="") #combine chain and resno of DSSP object
aro1SS="";aro2SS="";
	for(i in 1:length(aro1ChRes))
	{
		index=which(DSSPChRes==aro1ChRes[i]);
		if(length(index)==0){aro1SS=cbind(aro1SS,"-");}else{
		aro1SS=cbind(aro1SS,DSSP$ss[index]);
		}
		index=which(DSSPChRes==aro2ChRes[i]);
		if(length(index)==0){aro2SS=cbind(aro2SS,"-");}else{
		aro2SS=cbind(aro2SS,DSSP$ss[index]);
		}
	}
	aro1SS=aro1SS[-1];
	aro2SS=aro2SS[-1];
	aro=cbind(aro,resid1_SS=aro1SS);
	aro=cbind(aro,resid2_SS=aro2SS);

aro[which(aro[,1]==" "),1]="~";  # Replace ~ of chain1 with " "
aro[which(aro[,4]==" "),4]="~"; # Replace ~ chain2 with " "

return(aro);	
}
aroaro_addNaccess=function(nacc="",aro="")
{
#chain1     resno1     resid1     chain2     resno2     resid2   CentDist      Dihed
# "A"       "29"      "TYR"        "A"       "26"      "PHE"    "6.429"   "35.859"
	chres=paste(nacc$asa[,1],nacc$asa[,2],nacc$asa[,3],sep=""); #ChainResnoResid format
	aro[which(aro[,1]=="~"),1]=" ";
	aro[which(aro[,4]=="~"),4]=" ";
	aro1ChRes=paste(aro[,1],aro[,2],aro[,3],sep="") #Chain1, Resno1, Resid1
	aro2ChRes=paste(aro[,4],aro[,5],aro[,6],sep="") #Chain2, Resno2, Resid2
	aro1Acc=0;aro2Acc=0;
	for(i in 1:length(aro1ChRes))
	{
		index=which(chres==aro1ChRes[i]);
		if(length(index)==0){aro1Acc=cbind(aro1Acc,"-");}else{
		aro1Acc=cbind(aro1Acc,nacc$asa[index,5]); #asa_r_rel ASA residue relative
		}
		index=which(chres==aro2ChRes[i]);
		if(length(index)==0){aro2Acc=cbind(aro2Acc,"-");}else{
		aro2Acc=cbind(aro2Acc,nacc$asa[index,5]);
		}
	}
	aro1Acc=aro1Acc[-1];
	aro2Acc=aro2Acc[-1];
	aro=cbind(aro,resid1_Acc=aro1Acc);
	aro=cbind(aro,resid2_Acc=aro2Acc);
aro[which(aro[,1]==" "),1]="~";  # Replace " " of chain1 with "~"
aro[which(aro[,4]==" "),4]="~"; # Replace " "  chain2 with "~"
return(aro);	
}
aroaro_network=function(aro)
{
	aro1=paste(aro[,1],aro[,2],aro[,3],sep=""); #index 1,2,3 correspond to chain1, resno1,resid1
	aro2=paste(aro[,4],aro[,5],aro[,6],sep=""); #index 4,5,6 correspond to chain2, resno2,resid2
	if(length(aro1)!=length(aro2)){stop("No of aro1 and aro2 residues are not equal,check function aroaro_network: length(aro1)!=length(aro2)");}	
	# One residue can form interaction more than one residue
	# Hence we get the list of all residues involved in aroaro interaction by unique()
	# If n such residues are there then size of adjancy matrix wil;l be n*n
	total=unique(c(aro1,aro2)); 
	adjmat=matrix(0,nrow=length(total),ncol=length(total));
	# Loop over each aropair
	for(i in 1:length(aro1))
	{
	row=which(total==aro1[i]);col=which(total==aro2[i]);
	adjmat[row,col]=adjmat[row,col]+1;
	adjmat[col,row]=adjmat[col,row]+1;
	}
	rownames(adjmat)=total;
	#colnames(adjmat)=total;
	
	library("igraph");
	garopair=graph.adjacency(adjmat,mode="upper",add.rownames=NULL); #adjancy matyrix to graph, upper triangle only considered
	daropair=decompose.graph(garopair,mode="weak",max.comps=NA) # decomposition of graph into sub graph
	return(daropair);
}
aroaroSummary=function(aro="",acclow=20)
{
#colnames(aro):"chain1","resno1","resid1","chain2","resno2","resid2","CentDist","Dihed","resid1_SS","resid2_SS",resid1_Acc" "resid2_Acc"
	aroarosum=c(total=length(aro[,1]));
	
	#aroarosum=c(aroarosum,pintra=round((sum(aro[,1]==aro[,4])*100)/length(aro[,1]),2));
	#aroarosum=c(aroarosum,pInter=round((sum(aro[,1]!=aro[,4])*100)/length(aro[,1]),2));	

	aroarosum=c(aroarosum,pintra=sum(aro[,1]==aro[,4]));
	aroarosum=c(aroarosum,pInter=sum(aro[,1]!=aro[,4]));	


	#Total aromatic residues involved. For a pair FY, we say 2 aromatic residue involved.
	#For n interaction, 2n residues involved. Out of 2n % of F/Y/W are calculated.
	tot=c(aa321(aro[,3]),aa321(aro[,6])); # index 3:resid1 6:resid2
	#aroarosum=c(aroarosum,paste(c("F","Y","W"),c(round((sum(tot=="F")*100)/length(tot),2),round((sum(tot=="Y")*100)/length(tot),2),round((sum(tot=="W")*100)/length(tot),2)),sep=":",collapse=", "));
	
	#aroarosum=c(aroarosum,pF=round((sum(tot=="F")*100)/length(tot),2));
	#aroarosum=c(aroarosum,pY=round((sum(tot=="Y")*100)/length(tot),2));
	#aroarosum=c(aroarosum,pW=round((sum(tot=="W")*100)/length(tot),2));

	aroarosum=c(aroarosum,pF=sum(tot=="F"));
	aroarosum=c(aroarosum,pY=sum(tot=="Y"));
	aroarosum=c(aroarosum,pW=sum(tot=="W"));

	# Total aroaro pair of certain type FF,FW,FY...
	totpair=paste(aa321(aro[,3]),aa321(aro[,6]),sep="");
	
	#aroarosum=c(aroarosum,pFF=round((sum(totpair=="FF")*100)/length(totpair),2));
	#aroarosum=c(aroarosum,pFY=round((sum(totpair=="FY")*100)/length(totpair),2));
	#aroarosum=c(aroarosum,pFW=round((sum(totpair=="FW")*100)/length(totpair),2));
	#aroarosum=c(aroarosum,pYY=round((sum(totpair=="YY")*100)/length(totpair),2));
	#aroarosum=c(aroarosum,pYW=round((sum(totpair=="YW")*100)/length(totpair),2));
	#aroarosum=c(aroarosum,pWW=round((sum(totpair=="WW")*100)/length(totpair),2));
	
	aroarosum=c(aroarosum,pFF=sum(totpair=="FF"));
	aroarosum=c(aroarosum,pFY=sum(totpair=="FY"));
	aroarosum=c(aroarosum,pFW=sum(totpair=="FW"));
	aroarosum=c(aroarosum,pYY=sum(totpair=="YY"));
	aroarosum=c(aroarosum,pYW=sum(totpair=="YW"));
	aroarosum=c(aroarosum,pWW=sum(totpair=="WW"));
	
	
	
	
	
	#totpairfreq=round((tapply(totpair,factor(totpair),length)*100)/length(totpair),2);
	#aroarosum=c(aroarosum,paste(names(totpairfreq),totpairfreq,sep=":",collapse=", "));
	resid1_acc=as.numeric(aro[,"resid1_Acc"])
	resid2_acc=as.numeric(aro[,"resid2_Acc"])
	avg_Acc=(resid1_acc+resid2_acc)/2;
	aroarosum=c(aroarosum,pB=round((sum(avg_Acc<=acclow)*100)/length(avg_Acc),2));
	aroarosum=c(aroarosum,pE=round((sum(avg_Acc>acclow)*100)/length(avg_Acc),2));
	
	
	# Percentage of aroaropair of type BB,BE,EE. BE is same as EB where one is exposed and another buried
	#bb=round((sum(as.numeric(aro[,"resid1_Acc"]) < 30.0 & as.numeric(aro[,"resid2_Acc"]) < 30.0)*100)/length(aro[,1]),2)
	#be=round((sum(as.numeric(aro[,"resid1_Acc"]) < 30.0 & as.numeric(aro[,"resid2_Acc"]) > 30.0)*100)/length(aro[,1]),2)
	#eb=round((sum(as.numeric(aro[,"resid1_Acc"]) > 30.0 & as.numeric(aro[,"resid2_Acc"]) < 30.0)*100)/length(aro[,1]),2)
	#ee=round((sum(as.numeric(aro[,"resid1_Acc"]) > 30.0 & as.numeric(aro[,"resid2_Acc"]) > 30.0)*100)/length(aro[,1]),2)
	#aroarosum=c(aroarosum,paste("BB:",bb,", BE:",be+eb,", EE:",ee,sep=""));

	# Aropair Network
	daropair=aroaro_network(aro); #decomposed ionpair graph
	network=paste(sapply(daropair, vcount),sapply(daropair,ecount),sep="_");
	aroarosum=c(aroarosum,pIso=round((sum(network=="2_1")*100)/length(aro[,1]),2));
	
	netfreq=tapply(network,factor(network),length)		
	aroarosum=c(aroarosum,totNet=sum(sapply(daropair, vcount)>2));
	netfreq=netfreq[names(netfreq)!="2_1"];
	NetDet="-";
	if(length(netfreq)>0){NetDet=paste(names(netfreq),netfreq,sep=":",collapse=", ")}
	aroarosum=c(aroarosum,NetDet=NetDet);

return(aroarosum);
}

get_aroaronetmat=function(aro,daropair,acclow=20)
{
# This function returns aromatic network matrix
	size=sapply(daropair, vcount)
	noint=sapply(daropair, ecount);
	ind=which(size>2);
	aronetmat="";
	
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

	arores=c(paste(aro[,1],aro[,2],aro[,3],sep=""),paste(aro[,4],aro[,5],aro[,6],sep=""))
	aroacc=c(aro[,9],aro[,10])
	
	for(tempi in 1:length(ind))
	{
	edgelist=get.edgelist(daropair[[ind[tempi]]]);
	arores_net=unique(c(edgelist[,1],edgelist[,2]));
		acc=0;
		for(taro in 1:length(arores_net))
		{
			taro.ind=which(arores==arores_net[taro])[1]
			acc=c(acc,aroacc[taro.ind]);
		}
		acc=acc[-1];
		acc=as.numeric(acc); # We need to do bcoz acc>30 wont treat differently if acc is string
	#aronetmat=rbind(aronetmat,c(size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow & acc<acchigh),sum(acc>=acchigh),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")));
	aronetmat=rbind(aronetmat,c(size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")));
	}
	aronetmat=aronetmat[-1,];
	aronetmat=matrix(aronetmat,ncol=5);
	colnames(aronetmat)=c("TotRes","TotInt","BurRes","ExpoRes","Network");
  }
  return(aronetmat);
}


print_aroaronetwork=function(aro,daropair,outaroaro,acclow=20)
{
	size=sapply(daropair, vcount)
	noint=sapply(daropair, ecount);
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

	arores=c(paste(aro[,1],aro[,2],aro[,3],sep=""),paste(aro[,4],aro[,5],aro[,6],sep=""))
	aroacc=c(aro[,9],aro[,10])
	cat("No\tN.Res\tN.Int\tBuried\tExposed\tNetwork_details\n-----------------------------------",file=outaroaro,sep="",eol="\n",append=TRUE);

	for(tempi in 1:length(ind))
	{
	edgelist=get.edgelist(daropair[[ind[tempi]]]);
	arores_net=unique(c(edgelist[,1],edgelist[,2]));
		acc=0;
		for(taro in 1:length(arores_net))
		{
			taro.ind=which(arores==arores_net[taro])[1]
			acc=c(acc,aroacc[taro.ind]);
		}
		acc=acc[-1];
		acc=as.numeric(acc); # We need to do bcoz acc>30 wont treat differently if acc is string
#	cat(c(tempi,size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow & acc<acchigh),sum(acc>=acchigh),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")),file=outaroaro,sep="\t",eol="\n",append=TRUE)	
	cat(c(tempi,size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")),file=outaroaro,sep="\t",eol="\n",append=TRUE)	
	}
  }
}

print_aroaroSummary=function(aroarosum,outaroaro)
{
	cat("\n\nSummary\n----------\n",file=outaroaro,sep="",eol="",append=TRUE);
	aroarosumnames=c("Total no. of. aropairs: ","Percentage of intrachain aropairs: ",
	"Percentage of interchain aropairs: ","Percentage of aropairs involving F,Y,W: ",
	"Percentage of aropair of types: ","Percentage of aropair types BB,BE,EE (B:Buried,E:Exposed): ",
	"Percentage of Isolated aroaro pairs","No. of aropair network involving more than 2 residues: ","Aropair network distribution: ");
	newaroarosum=aroarosum[1:3];
	newaroarosum=c(newaroarosum,paste(c("F","Y","W"),aroarosum[4:6],sep=":",collapse=", "));
	newaroarosum=c(newaroarosum,paste(c("FF","FY","FW","YY","YW","WW"),aroarosum[7:12],sep=":",collapse=", "));
	newaroarosum=c(newaroarosum,paste(c("B","E"),aroarosum[13:14],sep=":",collapse=", "));
	newaroarosum=c(newaroarosum,aroarosum[15:17]);
	write.table(cbind(aroarosumnames,newaroarosum),file=outaroaro,sep=" ",eol="\n",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
}

getAroAro=function(chres="",aroaro)
{
	chain=chres[1];
	resno=chres[2];
	if(is.na(chain)){chain="~"}
	query=paste(chain,resno,sep="",collapse="");
	aroaro_res1=paste(aroaro[,1],aroaro[,2],sep="");
	aroaro_res2=paste(aroaro[,4],aroaro[,5],sep="");
	aroaro.ind=which(aroaro_res1==query | aroaro_res2==query);
	if(length(aroaro.ind)==0){return("SORRY");} # If not found return SORRY else an matrix
	aroaromat=aroaro[aroaro.ind,];
	aroaromat=matrix(aroaromat,ncol=12);
	colnames(aroaromat)=colnames(aroaro);
	return(aroaromat)
}

get_aroaronetwork=function(chres="",aro,daropair,acclow=20)
{
	chain=chres[1];
	resno=chres[2];
	query=paste(chain,resno,sep="",collapse="");
	size=sapply(daropair, vcount)
	noint=sapply(daropair, ecount);
	ind=which(size>2);
	Net="SORRY";
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
	arores=c(paste(aro[,1],aro[,2],aro[,3],sep=""),paste(aro[,4],aro[,5],aro[,6],sep=""))
	aroacc=c(aro[,9],aro[,10])
	for(tempi in 1:length(ind))
	{
	edgelist=get.edgelist(daropair[[ind[tempi]]]);
#	     [,1]      [,2]     
#[1,] "A221PHE" "A199TRP"
#[2,] "A170TYR" "A199TRP"
	edgeChres=gsub("...$","",as.vector(edgelist));  #"~341" "~341" "~410" "~344" "~411" "~411"  

		if(length(grep(paste(query,"$",sep="",collapse=""),edgeChres,perl=TRUE))>0)
		{
			arores_net=unique(c(edgelist[,1],edgelist[,2]));
			acc=0;
			for(taro in 1:length(arores_net))
			{
				taro.ind=which(arores==arores_net[taro])[1]
				acc=c(acc,aroacc[taro.ind]);
			}
			acc=acc[-1];
			acc=as.numeric(acc); # We need to do bcoz acc>30 wont treat differently if acc is string
		Net=rbind(Net,c(tempi,size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")));
		}
	
	}
	if(is.matrix(Net)){Net=Net[-1,];Net=matrix(Net,ncol=6);}
	return(Net);
	}else{
	return(Net); # Returns "";
	}
}

callAroAro=function(pdb,DSSP="",nacc="",status.DSSP=0,status.nacc=0,aroarolow,aroarohigh,aroarodihed,usedihed)
{
			# This is the workflow that is used in Pattern_analysis code
			# Aropair
			aro=aroaro(pdb,aroarolow,aroarohigh,aroarodihed,usedihed);
			daropair="";
			status.aroaro=0;
			if(is.matrix(aro))
			{
			status.aroaro=1;
			aro=Nochain(aro,1,4); #1st col is chain1 and 2nd col is 
			if(status.nacc){aro=aroaro_addNaccess(nacc,aro);}else{aro=cbind(aro,resid1_Acc="~");aro=cbind(aro,resid2_Acc="~");}
			if(status.DSSP){aro=aroaro_addDSSP(DSSP,aro);}else{aro=cbind(aro,resid1_SS="~");aro=cbind(aro,resid2_SS="~");}
			daropair=aroaro_network(aro);
			}
return(list(aro=aro,status.aroaro=status.aroaro,daropair=daropair));			
}

checkAroAro=function(patmat,aro)
{
	# This is the workflow that is used in Pattern_analysis code
	aroaromat=do.call(rbind,lapply(apply(patmat,1,getAroAro,aro),matrix,ncol=12,byrow=FALSE));
	sorry.ind=which(apply(aroaromat,1,function(x){sum(x=="SORRY")>0}));
	if(length(sorry.ind)>0){aroaromat=aroaromat[-sorry.ind,];}
	aroaromat=matrix(aroaromat,ncol=12);
	return(aroaromat);
}

checkAroAronet=function(patmat,aro,daropair,acclow=20)
{
	# This is the workflow that is used in Pattern_analysis code
	# Check for network
	aroaronetmat=do.call(rbind,lapply(apply(patmat,1,get_aroaronetwork,aro,daropair,acclow),matrix,ncol=6,byrow=FALSE));
	sorry.ind=which(apply(aroaronetmat,1,function(x){sum(x=="SORRY")>0}));
	if(length(sorry.ind)>0){aroaronetmat=aroaronetmat[-sorry.ind,];}
	aroaronetmat=matrix(aroaronetmat,ncol=6)
	aroaronetmat=matrix(aroaronetmat[,2:6],ncol=5)
	if(nrow(aroaronetmat)>0){
	aroaronetmat=matrix(unlist(strsplit(unique(apply(aroaronetmat,1,function(x){paste(x,sep="",collapse=":")})),":")),ncol=5,byrow=TRUE)
	}
return(aroaronetmat);							
}


print_HTMLaroaro=function(aro,aronetmat,filepath,append=FALSE,html=TRUE)
{
	if(html){write2Html("<html>\n",filepath,append);}else{write2Html("",filepath,append);}
	write2Html("<br><br><b><u>Aromatic-Aromatic interactions (AAI)</u><br><br></b>\n",filepath,TRUE);
#	write2Html("Chain1,Resno1,Resid1: Corresponds to Acidic residue details<br>\n",filepath,TRUE);
#	write2Html("Chain2,Resno2,Resid2: Corresponds to Basic residue details<br>\n",filepath,TRUE);
#	write2Html("Resid1_Acc and Resid2_Acc corresponds to relative ASA of acidic and basic residue respectively<br>\n",filepath,TRUE);
#	write2Html("Resid1_SS and Resid2_SS corresponds to secondary structure of the acidic and basic residues respectively<br><br><br>\n",filepath,TRUE);

	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr><th colspan=3 title=\"The first aromatic partner\"> Aromatic residue 1 </th> <th colspan=3 title=\"The second aromatic partner\"> Aromatic residue 2 </th> <th title=\"Distance between the centroid of the two aromatic ring\"> Centroid Distance </th><th title=\"The angle between two aromatic ring\"> Dihedral </th> <th colspan=2 title=\"Relative solvent accessibility of the residue\"> ASA </th> <th colspan=2 title=\"Secondary structure of residue (DSSP notation)\"> SS </th> </tr>",filepath,append=TRUE);
	tipvec=c("Chain","Residue number","Three letter code of residue","Chain","Residue number","Three letter code of residue","","","","","","");
	write2Html.tableHeader(c("Chain","Res.No","Res.ID","Chain","Res.No","Res.ID","(Ang)","(Deg)","Aro 1","Aro 2","Aro 1","Aro 2"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(aro,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	
	
	if(is.matrix(aronetmat) && nrow(aronetmat) > 0)
	{
	write2Html("<br><br><br><b><u>Aromatic-Aromatic interaction (AAI) network</u></b><br><br>\n",filepath,TRUE);
	#write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; width:1000px\">\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed;\">\n",filepath,TRUE);
	tipvec=c("Total number of residues present in the network","Total number of interactions present in the network","Number of residues of the network that are buried",
	"Number of residues of the network that are exposed","Format: Partner1 (Chain-Resno-Resid) -- Partner2 (Chain-Resno-Resid)"
	);
	write2Html.tableHeader(c("N.Res","N.Int","Buried","Exposed","Network_details"),filepath,TRUE,tipvec);
	#write2Html.tableBody(ionnetmat,filepath,TRUE,"<td style=\"word-wrap: break-word\">");
	aronetmat[,5]=apply(matrix(aronetmat[,5],ncol=1),1,function(x){gsub(",","<br>",x)})
	write2Html.tableBody(aronetmat,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	if(html){write2Html("</html>\n",filepath,TRUE);}else{write2Html("",filepath,TRUE);}
	}else{
	#write2Html("<br><br><br> <b>No Aromatic network is observed</b>\n",filepath,TRUE);
	if(html){write2Html("</html>\n",filepath,TRUE);}else{write2Html("",filepath,TRUE);}
	}
}

print_HTMLaroaro_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Aromatic-Aromatic interaction (AAI) summary</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("Filename","Total number of AAIs","Number of intra-chain AAIs","Number of inter-chain AAIs",
	paste("Number of AAIs involving residue: ",c("F","Y","W"),sep=""),
	paste("Number of AAIs of type [Aromatic residue 1-Aromatic residue 2]: ",c("FF","FY","FW","YY","YW","WW"),sep=""),
	"Percentage of buried AAIs","Percentage of exposed AAIs","Percentage of AAIs that are not part of any network (Isolated)",
	"Total number of AAI-networks","Network details [SIZE_STRENGTH:NUMBER OF SUCH NETWORK] (SIZE:Number of residue in the network, STRENGTH:Number of interaction in the network)"
	);
	write2Html.tableHeader(c("Filename","Total","IntraCh","InterCh","F","Y","W","FF","FY","FW","YY","YW","WW","B","E","Iso","T.Net","Net_Details"),filepath,append=TRUE,tipvec);
}
	
