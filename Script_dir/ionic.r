ionic=function(pdb,cutoff=3.2)
{
	#selects index of pdb$atom any oxygen atoms of Asp or Glu of given chain
	acidic=atom.select(pdb,resid=c("ASP","GLU"),elety=c("OD1","OD2","OE1","OE2"),verbose=FALSE)$atom
	#selects index of pdb$atom any nitrogen atoms of Asp or Glu of given chain
	basic=atom.select(pdb,resid=c("LYS","ARG","HIS"),elety=c("NZ","NE","NH1","NH2","ND1","NE2"),verbose=FALSE)$atom
	# Skip if either no DE or no RK exists, if any of acidic or basic group is missing, no point of ionic 
	if(length(acidic) & length(basic))
	{
		# If n acidic residues (length(acidic)=n) and m basic residues(length(basic)=m)
		# str(acidXYZmat) will be n*3 3:xyz and basicXYZmat=m*3
		acidXYZmat=matrix(as.numeric(pdb$atom[acidic,c("x","y","z")]),ncol=3)
		basicXYZmat=matrix(as.numeric(pdb$atom[basic,c("x","y","z")]),ncol=3)
		# We calculate distance between each row of acidXYZmat with each row of basicXYZmat
		# So size of distmat will be n*m str(distmat)= length(acidic) * length(basic)
		distmat=dist.xyz(acidXYZmat,basicXYZmat)
		arrind=which(distmat<=cutoff,arr.ind=TRUE);
		# If no interaction is within cutoff then arrind will be empty
		if(nrow(arrind)==0){return("SORRY");}
		# Fetch all residue information of acidic pair using arrind[,1]
		ionicInt=matrix(pdb$atom[acidic[arrind[,1]],c("chain","resno","resid","elety","insert")],ncol=5)
		temp=matrix(pdb$atom[basic[arrind[,2]],c("chain","resno","resid","elety","insert")],ncol=5);
		ionicInt=cbind(ionicInt,temp);
		ionicInt=cbind(ionicInt,round(distmat[arrind],2))
		colnames(ionicInt)=c("chain1","resno1","resid1","elety1","insert1","chain2","resno2","resid2","elety2","insert2","dist")
		return(ionicInt);
	}
	else{return("SORRY");}
	
}
ionic_dupRemove=function(ionpair="")
{

	ionpair[which(is.na(ionpair[,5])),5]="~";
	ionpair[which(is.na(ionpair[,10])),10]="~";
	uniq=unique(paste(ionpair[,1],ionpair[,2],ionpair[,3],ionpair[,5],ionpair[,6],ionpair[,7],ionpair[,8],ionpair[,10],sep="_"))
	ionpair=matrix(unlist(strsplit(uniq,"_")),ncol=8,byrow=TRUE);
	ionpair[which(ionpair[,4]=="~"),4]="";
	ionpair[which(ionpair[,8]=="~"),8]="";
	ionpair[,2]=paste(ionpair[,2],ionpair[,4],sep="")
	ionpair[,6]=paste(ionpair[,6],ionpair[,8],sep="")
	ionpair=matrix(ionpair[,c(1:3,5:7)],ncol=6)
	colnames(ionpair)=c("chain1","resno1","resid1","chain2","resno2","resid2");
	return(ionpair);
}
ionic_network=function(ionpair)
{
	acidic=paste(ionpair[,"chain1"],ionpair[,"resno1"],ionpair[,"resid1"],sep="");
	basic=paste(ionpair[,"chain2"],ionpair[,"resno2"],ionpair[,"resid2"],sep="");
	if(length(acidic)!=length(basic)){stop("No of acidic and basic residues are not equal,check function ionic_network: length(acidic)!=length(basic)");}	
	total=unique(c(acidic,basic));
	adjmat=matrix(0,nrow=length(total),ncol=length(total));
	# Loop over each ionpair
	for(i in 1:length(acidic))
	{
	adjmat[which(total==acidic[i]),which(total==basic[i])]=1;
	adjmat[which(total==basic[i]),which(total==acidic[i])]=1;
	}
	rownames(adjmat)=total;
	#colnames(adjmat)=total;
	
	library("igraph");
	gionpair=graph.adjacency(adjmat,mode="upper",add.rownames=NULL); #adjancy matyrix to graph, upper triangle only considered
	dionpair=decompose.graph(gionpair,mode="weak",max.comps=NA) # decomposition of graph into sub graph
	return(dionpair);
}

ionic_Nochain=function(ionpair,chain1index,chain2index)
{
# chain1index =1, chain2index=5
ionpair[is.na(ionpair[,chain1index]),chain1index]="~";
ionpair[is.na(ionpair[,chain2index]),chain2index]="~";
return(ionpair);
}

ionicSummary=function(ionpair,acclow=20)
{
	ionsum=c(total=length(ionpair[,1]));
	
	#ionsum=c(ionsum,pIntra=round((sum(ionpair[,1]==ionpair[,4])*100)/length(ionpair[,1]),2));
	#ionsum=c(ionsum,pInter=round((sum(ionpair[,1]!=ionpair[,4])*100)/length(ionpair[,1]),2));	

	ionsum=c(ionsum,pIntra=sum(ionpair[,1]==ionpair[,4]));
	ionsum=c(ionsum,pInter=sum(ionpair[,1]!=ionpair[,4]));	

	# Percentage of ionpair involving D/E/R/K/H
	resid1=aa321(ionpair[,3]) #resid1
	resid2=aa321(ionpair[,6]) #resid2
	
	########### PERCENTAGE ###############
	#ionsum=c(ionsum,pD=round((sum(resid1=="D")*100)/length(resid1),2));
	#ionsum=c(ionsum,pE=round((sum(resid1=="E")*100)/length(resid1),2));
	#ionsum=c(ionsum,pR=round((sum(resid2=="R")*100)/length(resid2),2));
	#ionsum=c(ionsum,pK=round((sum(resid2=="K")*100)/length(resid2),2));
	#ionsum=c(ionsum,pH=round((sum(resid2=="H")*100)/length(resid2),2));
	
	ionsum=c(ionsum,pD=sum(resid1=="D"));
	ionsum=c(ionsum,pE=sum(resid1=="E"));
	ionsum=c(ionsum,pR=sum(resid2=="R"));
	ionsum=c(ionsum,pK=sum(resid2=="K"));
	ionsum=c(ionsum,pH=sum(resid2=="H"));
	
	# Percentage of ionpair of type DR,DK,DH,ER,EK,EH
	iontype=paste(resid1,resid2,sep="")

	#ionsum=c(ionsum,pDK=round((sum(iontype=="DK")*100)/length(iontype),2));
	#ionsum=c(ionsum,pDR=round((sum(iontype=="DR")*100)/length(iontype),2));
	#ionsum=c(ionsum,pDH=round((sum(iontype=="DH")*100)/length(iontype),2));
	#ionsum=c(ionsum,pEK=round((sum(iontype=="EK")*100)/length(iontype),2));
	#ionsum=c(ionsum,pER=round((sum(iontype=="ER")*100)/length(iontype),2));
	#ionsum=c(ionsum,pEH=round((sum(iontype=="EH")*100)/length(iontype),2));
	
	ionsum=c(ionsum,pDK=sum(iontype=="DK"));
	ionsum=c(ionsum,pDR=sum(iontype=="DR"));
	ionsum=c(ionsum,pDH=sum(iontype=="DH"));
	ionsum=c(ionsum,pEK=sum(iontype=="EK"));
	ionsum=c(ionsum,pER=sum(iontype=="ER"));
	ionsum=c(ionsum,pEH=sum(iontype=="EH"));
	
	
	
	
	
	
	#iontype=paste(resid1,resid2,sep="")
	#piontype=round((tapply(iontype,factor(iontype),length)*100)/length(iontype),2);
	#ionsum=c(ionsum,paste(names(piontype),piontype,sep=":",collapse=", "));
	
	# Percentage of ionpair of type BB,BE,EB,EE for -ve and +ve pair
	#bb=round((sum(as.numeric(ionpair[,"resid1_Acc"]) < 30.0 & as.numeric(ionpair[,"resid2_Acc"]) < 30.0)*100)/length(ionpair[,1]),2)
	#be=round((sum(as.numeric(ionpair[,"resid1_Acc"]) < 30.0 & as.numeric(ionpair[,"resid2_Acc"]) > 30.0)*100)/length(ionpair[,1]),2)
	#eb=round((sum(as.numeric(ionpair[,"resid1_Acc"]) > 30.0 & as.numeric(ionpair[,"resid2_Acc"]) < 30.0)*100)/length(ionpair[,1]),2)
	#ee=round((sum(as.numeric(ionpair[,"resid1_Acc"]) > 30.0 & as.numeric(ionpair[,"resid2_Acc"]) > 30.0)*100)/length(ionpair[,1]),2)
	#ionsum=c(ionsum,paste("BB:",bb,", BE:",be,", EB:",eb,", EE:",ee,sep=""));
	
	resid1_acc=as.numeric(ionpair[,"resid1_Acc"])
	resid2_acc=as.numeric(ionpair[,"resid2_Acc"])
	avg_Acc=(resid1_acc+resid2_acc)/2;
	
	ionsum=c(ionsum,pB=round((sum(avg_Acc<=acclow)*100)/length(avg_Acc),2));
	ionsum=c(ionsum,pE=round((sum(avg_Acc>acclow)*100)/length(avg_Acc),2));
		
	# Ionpair Network
	dionpair=ionic_network(ionpair); #decomposed ionpair graph
	network=paste(sapply(dionpair, vcount),sapply(dionpair,ecount),sep="_");
	ionsum=c(ionsum,pIso=round((sum(network=="2_1")*100)/length(ionpair[,1]),2));
	netfreq=tapply(network,factor(network),length)		
	ionsum=c(ionsum,totNet=sum(sapply(dionpair, vcount)>2));
	netfreq=netfreq[names(netfreq)!="2_1"];
	NetDet="-";
	if(length(netfreq)>0){NetDet=paste(names(netfreq),netfreq,sep=":",collapse=", ")}
	ionsum=c(ionsum,NetDet=NetDet);
	return(ionsum);
}

ionpair_addDSSP=function(DSSP="",ionpair="")
{
ionpair[which(ionpair[,"chain1"]=="~"),"chain1"]=" ";
ionpair[which(ionpair[,"chain2"]=="~"),"chain2"]=" ";
acidChRes=paste(ionpair[,"chain1"],ionpair[,"resno1"],sep="")
baseChRes=paste(ionpair[,"chain2"],ionpair[,"resno2"],sep="")
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
	ionpair[which(ionpair[,"chain1"]==" "),"chain1"]="~";
	ionpair[which(ionpair[,"chain2"]==" "),"chain2"]="~";
	ionpair=cbind(ionpair,resid1_SS=acidSS);
	ionpair=cbind(ionpair,resid2_SS=baseSS);
return(ionpair);	
}


ionpair_addNaccess=function(nacc="",ionpair="")
{
	chres=paste(nacc$asa[,1],nacc$asa[,2],nacc$asa[,3],sep=""); #ChainResnoResid format
	# we have called the function ionic_Nochain which convert NA to ~
	# But dssp and naccess doesnot detect ~, first we need to convert to " ", then after adding ACC, reconvert back to ~
	ionpair[which(ionpair[,"chain1"]=="~"),"chain1"]=" ";
	ionpair[which(ionpair[,"chain2"]=="~"),"chain2"]=" ";
	acidChRes=paste(ionpair[,"chain1"],ionpair[,"resno1"],ionpair[,"resid1"],sep="")
	baseChRes=paste(ionpair[,"chain2"],ionpair[,"resno2"],ionpair[,"resid2"],sep="")
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
	ionpair[which(ionpair[,"chain1"]==" "),"chain1"]="~";
	ionpair[which(ionpair[,"chain2"]==" "),"chain2"]="~";
	ionpair=cbind(ionpair,resid1_Acc=acidAcc);
	ionpair=cbind(ionpair,resid2_Acc=baseAcc);
return(ionpair);	
}

print_ionicSummary=function(ionsum="",outionic)
{
cat("\n\nSummary\n----------\n",file=outionic,sep="",eol="",append=TRUE);

	ionsumname=c("Total number of ionpairs: ","Percentage of intrachain ionpairs: ",
	"Percentage of interchain ionpairs: ","Percentage of ionpairs involving the residues D,E,R,K,H: ",
	"Percentage of ionpair of type DR,DK,DH,ER,EK,EH: ","Percentage of Buried(B) or Exposed(E) ionpairs: ",
	"Percentage of Isolated ionpairs","No. of ionpair network involving more than 2 residues: ","Ionpair network distribution (Format: 'No. of. residue _ No. of. interaction : Number of such network'): ");
	newionsum=ionsum[1:3];
	newionsum=c(newionsum,paste(c("D","E","R","K","H"),ionsum[4:8],sep=":",collapse=", "));
	newionsum=c(newionsum,paste(c("DK","DR","DH","EK","ER","EH"),ionsum[9:14],sep=":",collapse=", "));
	newionsum=c(newionsum,paste(c("B","E"),ionsum[15:16],sep=":",collapse=", "));
	newionsum=c(newionsum,ionsum[17:19]);
	write.table(cbind(ionsumname,newionsum),file=outionic,sep=" ",eol="\n",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
	
}

get_ionnetmat=function(ionpair,dionpair,acclow=20)
{
# This function returns ionic network matrix
	size=sapply(dionpair, vcount);
	noint=sapply(dionpair, ecount);
	ind=which(size>2);
	ionnetmat="";
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
	
	acidic=paste(ionpair[,1],ionpair[,2],ionpair[,3],sep="");
	basic=paste(ionpair[,4],ionpair[,5],ionpair[,6],sep="");
	
	for(tempi in 1:length(ind))
	{
	edgelist=get.edgelist(dionpair[[ind[tempi]]]);
		acid_net=unique(edgelist[,1]);
		basic_net=unique(edgelist[,2]);
		acc=0;
		for(tacid in 1:length(acid_net))
		{
			tacid.ind=which(acidic==acid_net[tacid])[1]
			acc=c(acc,ionpair[tacid.ind,7]);
		}
		for(tbase in 1:length(basic_net))
		{
			tbase.ind=which(basic==basic_net[tbase])[1]
			acc=c(acc,ionpair[tbase.ind,8]);
		}
		acc=acc[-1];
		acc=as.numeric(acc); # We need to do bcoz acc>30 wont treat differently if acc is string
	#ionnetmat=rbind(ionnetmat,c(size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow & acc<acchigh),sum(acc>=acchigh),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")));
	ionnetmat=rbind(ionnetmat,c(size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")));		
	}
	ionnetmat=ionnetmat[-1,];
	ionnetmat=matrix(ionnetmat,ncol=5);
	colnames(ionnetmat)=c("TotRes","TotInt","BurRes","ExpoRes","Network");
   }
 return(ionnetmat);  
}

print_ionnetwork=function(ionpair,dionpair,outionic,acclow=20)
{
# This function print ionic network matrix
	size=sapply(dionpair, vcount);
	noint=sapply(dionpair, ecount);
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
	
	acidic=paste(ionpair[,1],ionpair[,2],ionpair[,3],sep="");
	basic=paste(ionpair[,4],ionpair[,5],ionpair[,6],sep="");
	cat("No\tN.Res\tN.Int\tBuried\tExposed\tNetwork_Details\n-----------------------------------",file=outionic,sep="",eol="\n",append=TRUE);
	
	for(tempi in 1:length(ind))
	{
	edgelist=get.edgelist(dionpair[[ind[tempi]]]);
		acid_net=unique(edgelist[,1]);
		basic_net=unique(edgelist[,2]);
		acc=0;
		for(tacid in 1:length(acid_net))
		{
			tacid.ind=which(acidic==acid_net[tacid])[1]
			acc=c(acc,ionpair[tacid.ind,7]);
		}
		for(tbase in 1:length(basic_net))
		{
			tbase.ind=which(basic==basic_net[tbase])[1]
			acc=c(acc,ionpair[tbase.ind,8]);
		}
		acc=acc[-1];
		acc=as.numeric(acc); # We need to do bcoz acc>30 wont treat differently if acc is string
		
	cat(c(tempi,size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")),file=outionic,sep="\t",eol="\n",append=TRUE)
	}
   }
}

getIonic=function(chres="",ionpair)
{
# You pass a chain and resno, it will return local ionpair matrix
	chain=chres[1];
	resno=chres[2];
	
	if(is.na(chain)){chain="~"}
	query=paste(chain,resno,sep="",collapse="");
	ip_res1=paste(ionpair[,1],ionpair[,2],sep="");
	ip_res2=paste(ionpair[,4],ionpair[,5],sep="");
	ip.ind=which(ip_res1==query | ip_res2==query);
	if(length(ip.ind)==0){return("SORRY");} # If not found return SORRY else an matrix
	ipmat=ionpair[ip.ind,];
	ipmat=matrix(ipmat,ncol=10);
	colnames(ipmat)=colnames(ionpair);
	return(ipmat)
}

get_ionnetwork=function(chres="",ionpair,dionpair,acclow=20)
{
# You pass a chain and resno, it will return local ionnet matrix
	chain=chres[1];
	resno=chres[2];
	query=paste(chain,resno,sep="",collapse="");
	size=sapply(dionpair, vcount);
	noint=sapply(dionpair, ecount);
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
	
	acidic=paste(ionpair[,1],ionpair[,2],ionpair[,3],sep="");
	basic=paste(ionpair[,4],ionpair[,5],ionpair[,6],sep="");
	
	for(tempi in 1:length(ind))
	{
	edgelist=get.edgelist(dionpair[[ind[tempi]]]); 
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
				acc=c(acc,ionpair[tacid.ind,7]);
			}
			for(tbase in 1:length(basic_net))
			{
				tbase.ind=which(basic==basic_net[tbase])[1]
				acc=c(acc,ionpair[tbase.ind,8]);
			}
			acc=acc[-1];
			acc=as.numeric(acc); # We need to do bcoz acc>30 wont treat differently if acc is string
			#Net=rbind(Net,c(tempi,size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow & acc<acchigh),sum(acc>=acchigh),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")));
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

callIonic=function(pdb,DSSP="",nacc="",status.DSSP=0,status.nacc=0,ioncut=6.0)
{
	# This is the workflow that is used in Pattern_analysis code
	# Ionpair
	ionpair=ionic(pdb,cutoff=ioncut);
	dionpair="";
	status.ionic=0;
	if(is.matrix(ionpair))
	{
	status.ionic=1;
	ionpair=Nochain(ionpair,1,6);
	ionpair=ionic_dupRemove(ionpair);
	if(status.nacc){ionpair=ionpair_addNaccess(nacc,ionpair);}else{ionpair=cbind(ionpair,resid1_Acc="-");ionpair=cbind(ionpair,resid2_Acc="-");}
	if(status.DSSP){ionpair=ionpair_addDSSP(DSSP,ionpair);}else{ionpair=cbind(ionpair,resid1_SS="-");ionpair=cbind(ionpair,resid2_SS="-");}
	dionpair=ionic_network(ionpair);
	}
return(list(ionpair=ionpair,status.ionic=status.ionic,dionpair=dionpair));	
}

checkIonic=function(patmat,ionpair)
{
	# This is the workflow that is used in Pattern_analysis code
	# refer to hbond.r checkHB function for details
	ionmat=do.call(rbind,lapply(apply(patmat,1,getIonic,ionpair),matrix,ncol=10,byrow=FALSE));
	sorry.ind=which(apply(ionmat,1,function(x){sum(x=="SORRY")>0}));
	if(length(sorry.ind)>0){ionmat=ionmat[-sorry.ind,];}
	ionmat=matrix(ionmat,ncol=10);
	# For 3 residue pattern, if all 3 residue doesnot form any inoic interaction, we get 3*10 matrix with all value sorry
	# And after ionmat[-,] i.e. removing those SORRY row we get an matrix but nrow(ionmat)=0
	return(ionmat);				
}
checkIonnet=function(patmat,ionpair,dionpair,acclow=20)
{
	# This is the workflow that is used in Pattern_analysis code
	#check for network

	ionnetmat=do.call(rbind,lapply(apply(patmat,1,get_ionnetwork,ionpair,dionpair,acclow),matrix,ncol=6,byrow=FALSE));
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


print_HTMLionic=function(ionpair,ionnetmat,filepath,append=FALSE,html=TRUE)
{
	if(html){write2Html("<html>\n",filepath,append);}else{write2Html("",filepath,append);}
	write2Html("<br><br><b><u>Ionpair interactions</u><br><br></b>\n",filepath,TRUE);
	#write2Html("Chain1,Resno1,Resid1: Corresponds to Acidic residue details<br>\n",filepath,TRUE);
	#write2Html("Chain2,Resno2,Resid2: Corresponds to Basic residue details<br>\n",filepath,TRUE);
	#write2Html("Resid1_Acc and Resid2_Acc corresponds to relative ASA of acidic and basic residue respectively<br>\n",filepath,TRUE);
	#write2Html("Resid1_SS and Resid2_SS corresponds to secondary structure of the acidic and basic residues respectively<br><br><br>\n",filepath,TRUE);

	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr><th colspan=3 title=\" Residue type: Acidic [ Asp | Glu ]\" > Acidic </th> <th colspan=3 title=\" Residue type:Basic [ Lys | Arg | His ]\" > Basic </th> <th colspan=2 title=\"Relative solvent accessibility of the residue\"> ASA </th> <th colspan=2 title=\"Secondary structure of residue (DSSP notation)\"> SS </th> </tr>",filepath,append=TRUE);
	tipvec=c("Chain","Residue number","Three letter code of residue","Chain","Residue number","Three letter code of residue","","","","");
	write2Html.tableHeader(c("Chain","Res.No","Res.ID","Chain","Res.No","Res.ID","Acidic","Basic","Acidic","Basic"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(ionpair,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	
	
	if(is.matrix(ionnetmat) && nrow(ionnetmat) > 0)
	{
	write2Html("<br><br><br><b><u>Ionic network</u></b><br><br>\n",filepath,TRUE);
	#write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; width:1000px\">\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed;\">\n",filepath,TRUE);
	tipvec=c("Total number of residues present in the network","Total number of interactions present in the network","Number of residues of the network that are buried",
	"Number of residues of the network that are exposed","Format: Partner1 (Chain-Resno-Resid) -- Partner2 (Chain-Resno-Resid)"
	);
	write2Html.tableHeader(c("N.Res","N.Int","Buried","Exposed","Network_Details"),filepath,TRUE,tipvec);
	#write2Html.tableBody(ionnetmat,filepath,TRUE,"<td style=\"word-wrap: break-word\">");
	ionnetmat[,5]=apply(matrix(ionnetmat[,5],ncol=1),1,function(x){gsub(",","<br>",x)})
	write2Html.tableBody(ionnetmat,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	if(html){write2Html("</html>\n",filepath,TRUE);}else{write2Html("",filepath,TRUE);}
	}else{
	#write2Html("<br><br><br>No ionic network is observed\n",filepath,TRUE);
	if(html){write2Html("</html>\n",filepath,TRUE);}else{write2Html("",filepath,TRUE);}
	}
}

print_HTMLionic_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Ionpair summary</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("Filename","Total number of ion-pairs","Number of intra-chain IPs","Number of inter-chain IPs",
	paste("Number of IPs involving residue: ",c("D","E","R","K","H"),sep=""),
	paste("Number of IPs of type [Acidic-Basic]: ",c("DK","DR","DH","EK","ER","EH"),sep=""),
	"Percentage of buried IPs","Percentage of exposed IPs","Percentage of IPs that are not part of any network (Isolated)",
	"Total number of IP-networks","Network details [SIZE_STRENGTH:NUMBER OF SUCH NETWORK] (SIZE:Number of residue in the network, STRENGTH:Number of interaction in the network)"
	);
	write2Html.tableHeader(c("Filename","Total","IntraCh","InterCh","D","E","R","K","H","DK","DR","DH","EK","ER","EH","B","E","Iso","T.Net","Net_Details"),filepath,append=TRUE,tipvec);
}
print_HTMLionic_sum_end=function(filepath)
{
	write2Html("\n</table>\n</html>",filepath,TRUE);
}
	
