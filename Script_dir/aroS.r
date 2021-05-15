aroS=function(pdb="",cutoff=5.3)
{
	aroSint="";
	phe=aroCentroid(pdb,restype="PHE");
	tyr=aroCentroid(pdb,restype="TYR");
	trp=aroCentroid(pdb,restype="TRP");
	
FYW=matrix(0,ncol=9);
	
	# If no phe/y/w found we get NA. so if no NA then only u do rbind
	if(is.matrix(phe)){FYW=rbind(FYW,phe);}
	if(is.matrix(tyr)){FYW=rbind(FYW,tyr);}
	if(is.matrix(trp)){FYW=rbind(FYW,trp);}
	FYW=FYW[-1,]
	FYW=matrix(FYW,ncol=9);
	colnames(FYW)=c("chain","resno","resid","Cx","Cy","Cz","EqnA","EqnB","EqnC");
	# If no phe/y/w then no interaction so return
	if(nrow(FYW)==0){return("SORRY");}

	sul.detail=matrix(pdb$atom[atom.select(pdb,resid=c("CYS","MET"),elety=c("SG","SD"),verbose=FALSE)$atom,c("chain","resno","resid","insert","x","y","z")],ncol=7)
	sul.detail[which(is.na(sul.detail[,4])),4]="";
	sul.detail[,2]=paste(sul.detail[,2],sul.detail[,4],sep="")

	sul.detail=matrix(sul.detail[,c(1:3,5:7)],ncol=6)

	# For no chain info NA to be changed to ~
	sul.detail[which(is.na(sul.detail[,1])),1]="~";
	# if no aromatic residues or no S, then no aro-S interaction.
	if(dim(sul.detail)[1]==0) {return("SORRY");}
	
	aroSdist=dist.xyz(matrix(as.numeric(FYW[,c("Cx","Cy","Cz")]),ncol=3),matrix(as.numeric(sul.detail[,4:6]),ncol=3))
	aroSdist.Cutoff=which(aroSdist<cutoff,arr.ind=TRUE);
	# If no interaction is observed
	if(dim(aroSdist.Cutoff)[1]==0)
	{
	return("SORRY");
	}
	aroSint=matrix(FYW[aroSdist.Cutoff[,1],1:3],ncol=3)
	aroSint=cbind(aroSint,matrix(sul.detail[aroSdist.Cutoff[,2],1:3],ncol=3))
	aroSint=cbind(aroSint,round(aroSdist[aroSdist.Cutoff],2))
	
	colnames(aroSint)=c("chain1","resno1","resid1","chain2","resno2","resid2","Dist")
	#rownames(aroSint)=1:length(aroSint[,1]);
	
	return(aroSint);
}
aros_addDSSP=function(DSSP="",aros="")
{
	aros[which(aros[,1]=="~"),1]=" ";
	aros[which(aros[,4]=="~"),4]=" ";
	aroChRes=paste(aros[,1],aros[,2],sep="")
	SChRes=paste(aros[,4],aros[,5],sep="")
	DSSPChRes=paste(DSSP$cha,DSSP$res,sep="") #combine chain and resno of DSSP object
	aroSS="";SSS="";
	for(i in 1:length(aroChRes))
	{
		index=which(DSSPChRes==aroChRes[i]);
		if(length(index)==0){aroSS=cbind(aroSS,"-");}else{
		aroSS=cbind(aroSS,DSSP$ss[index]);
		}
		index=which(DSSPChRes==SChRes[i]);
		if(length(index)==0){SSS=cbind(SSS,"-");}else{
		SSS=cbind(SSS,DSSP$ss[index]);
		}
	}
	aroSS=aroSS[-1];
	SSS=SSS[-1];
	aros=cbind(aros,resid1_SS=aroSS);
	aros=cbind(aros,resid2_SS=SSS);
	aros[which(aros[,1]==" "),1]="~";  # Replace " " of chain1 with "~"
	aros[which(aros[,4]==" "),4]="~"; # Replace " "  chain2 with "~"
	return(aros);
}
aros_addNaccess=function(nacc="",aros="")
{
	chres=paste(nacc$asa[,1],nacc$asa[,2],nacc$asa[,3],sep=""); #ChainResnoResid format
	aros[which(aros[,1]=="~"),1]=" ";
	aros[which(aros[,4]=="~"),4]=" ";
	aroChRes=paste(aros[,1],aros[,2],aros[,3],sep="") #Chain1, Resno1, Resid1
	SChRes=paste(aros[,4],aros[,5],aros[,6],sep="") #Chain2, Resno2, Resid2
	aroAcc=0;SAcc=0;
	for(i in 1:length(aroChRes))
	{
		index=which(chres==aroChRes[i]);
		if(length(index)==0){aroAcc=cbind(aroAcc,"-");}else{
		aroAcc=cbind(aroAcc,nacc$asa[index,5]); #asa_r_rel ASA residue relative
		}
		index=which(chres==SChRes[i]);
		if(length(index)==0){SAcc=cbind(SAcc,"-");}else{
		SAcc=cbind(SAcc,nacc$asa[index,5]);
		}
	}
	aroAcc=aroAcc[-1];
	SAcc=SAcc[-1];
	aros=cbind(aros,resid1_Acc=aroAcc);
	aros=cbind(aros,resid2_Acc=SAcc);
	aros[which(aros[,1]==" "),1]="~";  # Replace " " of chain1 with "~"
	aros[which(aros[,4]==" "),4]="~"; # Replace " "  chain2 with "~"
	return(aros);
}

aros_network=function(aros)
{
	aro=paste(aros[,1],aros[,2],aros[,3],sep=""); #index 1,2,3 correspond to chain1, resno1,resid1
	S=paste(aros[,4],aros[,5],aros[,6],sep=""); #index 4,5,6 correspond to chain2, resno2,resid2
	total=unique(c(aro,S)); 
	adjmat=matrix(0,nrow=length(total),ncol=length(total));
	# Loop over each aropair
	for(i in 1:length(aro))
	{
	row=which(total==aro[i]);
	col=which(total==S[i]);
	adjmat[row,col]=1;
	adjmat[col,row]=1;
	}
	rownames(adjmat)=total;
	#colnames(adjmat)=total;
	
	library("igraph");
	garospair=graph.adjacency(adjmat,mode="upper",add.rownames=NULL); #adjancy matyrix to graph, upper triangle only considered
	darospair=decompose.graph(garospair,mode="weak",max.comps=NA) # decomposition of graph into sub graph
	return(darospair);
}
arosSummary=function(aros,acclow=20)
{
	len=length(aros[,1]); #No of aroS interaction
	arossum=c(total=len);
	#arossum=c(arossum,pIntra=round((sum(aros[,1]==aros[,4])*100)/len,2)); # % Intra chain
	#arossum=c(arossum,pInter=round((sum(aros[,1]!=aros[,4])*100)/len,2)); # % Inter chain
	
	arossum=c(arossum,pIntra=sum(aros[,1]==aros[,4])); # % Intra chain
	arossum=c(arossum,pInter=sum(aros[,1]!=aros[,4])); # % Inter chain
	
	# Percentage involvinf F/Y/W
	#arossum=c(arossum,pF=round((sum(aa321(aros[,3])=="F")*100)/len,2));
	#arossum=c(arossum,pY=round((sum(aa321(aros[,3])=="Y")*100)/len,2));
	#arossum=c(arossum,pW=round((sum(aa321(aros[,3])=="W")*100)/len,2));
	
	arossum=c(arossum,pF=sum(aa321(aros[,3])=="F"));
	arossum=c(arossum,pY=sum(aa321(aros[,3])=="Y"));
	arossum=c(arossum,pW=sum(aa321(aros[,3])=="W"));
	
	
	#arossum=c(arossum,pC=round((sum(aa321(aros[,6])=="C")*100)/len,2));
	#arossum=c(arossum,pM=round((sum(aa321(aros[,6])=="M")*100)/len,2));
	
	arossum=c(arossum,pC=sum(aa321(aros[,6])=="C"));
	arossum=c(arossum,pM=sum(aa321(aros[,6])=="M"));
	
	resid1=aa321(aros[,3]) #resid1
	resid2=aa321(aros[,6]) #resid2
	arostype=paste(resid1,resid2,sep="")
	
	arossum=c(arossum,pFC=sum(arostype=="FC"));
	arossum=c(arossum,pYC=sum(arostype=="YC"));
	arossum=c(arossum,pWC=sum(arostype=="WC"));
	arossum=c(arossum,pFM=sum(arostype=="FM"));
	arossum=c(arossum,pYM=sum(arostype=="YM"));
	arossum=c(arossum,pWM=sum(arostype=="WM"));
	
	
	#arossum=c(arossum,paste(c("F","Y","W"),c(round((sum(aa321(aros[,3])=="F")*100)/len,2),round((sum(aa321(aros[,3])=="Y")*100)/len,2),round((sum(aa321(aros[,3])=="W")*100)/len,2)),sep=":",collapse=", "));
	
	resid1_acc=as.numeric(aros[,"resid1_Acc"])
	resid2_acc=as.numeric(aros[,"resid2_Acc"])
	avg_Acc=(resid1_acc+resid2_acc)/2;
	arossum=c(arossum,pB=round((sum(avg_Acc<=acclow)*100)/length(avg_Acc),2));
	arossum=c(arossum,pE=round((sum(avg_Acc>acclow)*100)/length(avg_Acc),2));
	
	
	darospair=aros_network(aros);
	network=paste(sapply(darospair, vcount),sapply(darospair,ecount),sep="_");
	arossum=c(arossum,pIso=round((sum(network=="2_1")*100)/len,2));
	
	netfreq=tapply(network,factor(network),length)		
	arossum=c(arossum,totNet=sum(sapply(darospair, vcount)>2));
	netfreq=netfreq[names(netfreq)!="2_1"];
	NetDet="-";
	if(length(netfreq)>0){NetDet=paste(names(netfreq),netfreq,sep=":",collapse=", ")}
	
	arossum=c(arossum,NetDet=NetDet);
	return(arossum);
}
print_arosnetwork=function(aros,darospair,outaros,acclow=20)
{
	size=sapply(darospair, vcount)
	noint=sapply(darospair, ecount);
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
	ind=ind[order];

	aro=paste(aros[,1],aros[,2],aros[,3],sep="");
	S=paste(aros[,4],aros[,5],aros[,6],sep="");

	
	cat("No\tN.Res\tN.Int\tBuried\tExposed\tNetwork_Details\n-----------------------------------",file=outaros,sep="",eol="\n",append=TRUE);
	
	for(tempi in 1:length(ind))
	{
	edgelist=get.edgelist(darospair[[ind[tempi]]]);
		aro_net=unique(edgelist[,1]);
		S_net=unique(edgelist[,2]);
		acc=0;
		for(taro in 1:length(aro_net))
		{
			taro.ind=which(aro==aro_net[taro])[1]
			acc=c(acc,aros[taro.ind,8]);
		}
		for(tS in 1:length(S_net))
		{
			tS.ind=which(S==S_net[tS])[1]
			acc=c(acc,aros[tS.ind,9]);
		}
		acc=acc[-1];
		acc=as.numeric(acc); # We need to do bcoz acc>30 wont treat differently if acc is string
		
	cat(c(tempi,size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")),file=outaros,sep="\t",eol="\n",append=TRUE)
	}
   }
}

get_arosnetmat=function(aros,darospair,acclow=20)
{
	size=sapply(darospair, vcount)
	noint=sapply(darospair, ecount);
	ind=which(size>2);
	arosnetmat="";
	
	# If network present then only proceed
	if(length(ind)>0)
	{
	size3=size[ind];
	noint3=noint[ind];
	# order them
	order=order(size3)
	size3=size3[order]
	noint3=noint3[order]
	ind=ind[order];

	aro=paste(aros[,1],aros[,2],aros[,3],sep="");
	S=paste(aros[,4],aros[,5],aros[,6],sep="");

	
	
	
	for(tempi in 1:length(ind))
	{
	edgelist=get.edgelist(darospair[[ind[tempi]]]);
		aro_net=unique(edgelist[,1]);
		S_net=unique(edgelist[,2]);
		acc=0;
		for(taro in 1:length(aro_net))
		{
			taro.ind=which(aro==aro_net[taro])[1]
			acc=c(acc,aros[taro.ind,8]);
		}
		for(tS in 1:length(S_net))
		{
			tS.ind=which(S==S_net[tS])[1]
			acc=c(acc,aros[tS.ind,9]);
		}
		acc=acc[-1];
		acc=as.numeric(acc); # We need to do bcoz acc>30 wont treat differently if acc is string
	arosnetmat=rbind(arosnetmat,c(size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")));	
	}
	arosnetmat=arosnetmat[-1,];
	arosnetmat=matrix(arosnetmat,ncol=5);
	colnames(arosnetmat)=c("TotRes","TotInt","BurRes","ExpoRes","Network");
   }
   return(arosnetmat);
}

print_arosSummary=function(arossum,outaros)
{
	cat("\n\nSummary\n----------\n",file=outaros,sep="",eol="",append=TRUE);
	arossumnames=c("Total no. of Aro-S interaction: ","Percentage of intrachain Aro-S interaction: ","Percentage of interchain Aro-S interaction: ",
	"Percentage of Aro-S interaction involving residues F/Y/W: ","Percentage of Isolated Aro-S pairs","No. of Aro-S network involving more than 2 residues: ",
	"Aro-S network distribution: ");
	newarossum=arossum[1:3];
	newarossum=c(newarossum,paste(c("F","Y","W"),arossum[4:6],sep=":",collapse=", "));
	newarossum=c(newarossum,arossum[7:9]);
	write.table(cbind(arossumnames,newarossum),file=outaros,sep=" ",eol="\n",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE);
}

getAroS=function(chres="",aros)
{
	chain=chres[1];
	resno=chres[2];
	if(is.na(chain)){chain="~"}
	query=paste(chain,resno,sep="",collapse="");
	aros_res1=paste(aros[,1],aros[,2],sep=""); #A112
	aros_res2=paste(aros[,4],aros[,5],sep=""); #A112
	aros.ind=which(aros_res1==query | aros_res2==query);
	if(length(aros.ind)==0){return("SORRY");} # If not found return SORRY else an matrix
	arosmat=aros[aros.ind,];
	arosmat=matrix(arosmat,ncol=11);
	colnames(arosmat)=colnames(aros);
	return(arosmat)
}


get_aroSnetwork=function(chres="",aros,darospair,acclow=20)
{
	#print(aros);
	chain=chres[1];
	resno=chres[2];
	query=paste(chain,resno,sep="",collapse="");
	#print(c("query",query));
	size=sapply(darospair, vcount)
	noint=sapply(darospair, ecount);
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
	arosres=c(paste(aros[,1],aros[,2],aros[,3],sep=""),paste(aros[,4],aros[,5],aros[,6],sep=""))
	#print(c("arosres",arosres));
	arosacc=c(aros[,8],aros[,9])
	#print(c("arosacc",arosacc));
	for(tempi in 1:length(ind))
	{
	edgelist=get.edgelist(darospair[[ind[tempi]]]);
#	     [,1]      [,2]     
#[1,] "A221PHE" "A199TRP"
#[2,] "A170TYR" "A199TRP"
	edgeChres=gsub("...$","",as.vector(edgelist));  #"~341" "~341" "~410" "~344" "~411" "~411"  
	
		if(length(grep(paste(query,"$",sep="",collapse=""),edgeChres,perl=TRUE))>0)
		{
	#		print(c("edgeChres",edgeChres));
			arosres_net=unique(c(edgelist[,1],edgelist[,2]));
	#		print(c("arosres_net",arosres_net));
			acc=0;
			for(taros in 1:length(arosres_net))
			{
				taros.ind=which(arosres==arosres_net[taros])[1]
				acc=c(acc,arosacc[taros.ind]);
			}
			acc=acc[-1];
			acc=as.numeric(acc); # We need to do bcoz acc>30 wont treat differently if acc is string
	#		print(c("----------",acc));
		Net=rbind(Net,c(tempi,size3[tempi],noint3[tempi],sum(acc<=acclow),sum(acc>acclow),paste(edgelist[,1],edgelist[,2],sep="-",collapse=", ")));
		}
	
	}
	if(is.matrix(Net)){Net=Net[-1,];Net=matrix(Net,ncol=6);}
	return(Net);
	}else{
	return(Net); # Returns "";
	}
}

callAroS=function(pdb="",DSSP="",nacc="",status.DSSP=0,status.nacc=0,aroscut)
{
			aros=aroS(pdb,cutoff=aroscut);
			darospair="";
			status.aros=0;
			if(is.matrix(aros))
			{
				status.aros=1;
				aros=Nochain(aros,1,4);
				if(status.nacc){aros=aros_addNaccess(nacc,aros);}else{aros=cbind(aros,resid1_Acc="~");aros=cbind(aros,resid2_Acc="~");}
				if(status.DSSP){aros=aros_addDSSP(DSSP,aros);}else{aros=cbind(aros,resid1_SS="~");aros=cbind(aros,resid2_SS="~");}
				darospair=aros_network(aros);
			}
return(list(aros=aros,status.aros=status.aros,darospair=darospair));			
}

checkAroS=function(patmat,aros)
{
	# This is the workflow that is used in Pattern_analysis code
	arosmat=do.call(rbind,lapply(apply(patmat,1,getAroS,aros),matrix,ncol=11,byrow=FALSE));
	sorry.ind=which(apply(arosmat,1,function(x){sum(x=="SORRY")>0}));
	if(length(sorry.ind)>0){arosmat=arosmat[-sorry.ind,];}
	arosmat=matrix(arosmat,ncol=11);
	return(arosmat);
}
checkAroSnet=function(patmat,aros,darospair,acclow=20)
{
	# This is the workflow that is used in Pattern_analysis code
	# Check for network
	arosnetmat=do.call(rbind,lapply(apply(patmat,1,get_aroSnetwork,aros,darospair,acclow),matrix,ncol=6,byrow=FALSE));
	sorry.ind=which(apply(arosnetmat,1,function(x){sum(x=="SORRY")>0}));
	if(length(sorry.ind)>0){arosnetmat=arosnetmat[-sorry.ind,];}
	arosnetmat=matrix(arosnetmat,ncol=6)
	arosnetmat=matrix(arosnetmat[,2:6],ncol=5)
	if(nrow(arosnetmat)>0){
	arosnetmat=matrix(unlist(strsplit(unique(apply(arosnetmat,1,function(x){paste(x,sep="",collapse=":")})),":")),ncol=5,byrow=TRUE)
	}
return(arosnetmat);							
}

print_HTMLaroS=function(aros,arosnetmat,filepath,append=FALSE,html=TRUE)
{
	if(html){write2Html("<html>\n",filepath,append);}else{write2Html("",filepath,append);}
	write2Html("<br><br><b><u>Aromatic-Sulphur interactions (ASI)</u><br><br></b>\n",filepath,TRUE);
#	write2Html("Chain1,Resno1,Resid1: Corresponds to Acidic residue details<br>\n",filepath,TRUE);
#	write2Html("Chain2,Resno2,Resid2: Corresponds to Basic residue details<br>\n",filepath,TRUE);
#	write2Html("Resid1_Acc and Resid2_Acc corresponds to relative ASA of acidic and basic residue respectively<br>\n",filepath,TRUE);
#	write2Html("Resid1_SS and Resid2_SS corresponds to secondary structure of the acidic and basic residues respectively<br><br><br>\n",filepath,TRUE);

	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr><th colspan=3 title=\"The aromatic residue [ FWY ]\"> Aromatic residue </th> <th colspan=3 title=\"Sulphur containing residues [ CM ]\"> Sulphur residue </th> <th title=\"Distance between centroid of aromatic ring and Sulphur atom\"> Dist</th> <th colspan=2 title=\"Relative solvent accessibility of the residue\"> ASA </th> <th colspan=2 title=\"Secondary structure of residue (DSSP notation)\"> SS </th> </tr>",filepath,append=TRUE);
	tipvec=c("Chain","Residue number","Three letter code of residue","Chain","Residue number","Three letter code of residue","","","","","");
	write2Html.tableHeader(c("Chain","Res.No","Res.ID","Chain","Res.No","Res.ID","(Ang)","Aro","S","Aro","S"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(aros,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	
	#print(is.matrix(arosnetmat) & nrow(arosnetmat)>0);
	if(is.matrix(arosnetmat) && nrow(arosnetmat)>0)
	{
	write2Html("<br><br><b><u>Aromatic-sulphur interaction (ASI) network</u></b><br><br>\n",filepath,TRUE);
	#write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; width:1000px\">\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed;\">\n",filepath,TRUE);
	tipvec=c("Total number of residues present in the network","Total number of interactions present in the network","Number of residues of the network that are buried",
	"Number of residues of the network that are exposed","Format: Partner1 (Chain-Resno-Resid) -- Partner2 (Chain-Resno-Resid)"
	);
	
	write2Html.tableHeader(c("N.Res","N.Int","Buried","Exposed","Network_details"),filepath,TRUE,tipvec);	
	#write2Html.tableBody(ionnetmat,filepath,TRUE,"<td style=\"word-wrap: break-word\">");
	arosnetmat[,5]=apply(matrix(arosnetmat[,5],ncol=1),1,function(x){gsub(",","<br>",x)})
	write2Html.tableBody(arosnetmat,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	if(html){write2Html("</html>\n",filepath,TRUE);}else{write2Html("",filepath,TRUE);}
	}else{
	#write2Html("<br><br><br> <b>No Aromatic-sulphor network is observed</b>\n",filepath,TRUE);
	#print("///");
	if(html){write2Html("</html>\n",filepath,TRUE);}else{write2Html("",filepath,TRUE);}
	}
}

print_HTMLaros_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Aromatic Sulphur interaction (ASI) summary</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("PDB Filename","Total number of ASIs","Number of intra-chain ASIs","Number of inter-chain ASIs",
	paste("Number of ASIs involving residue: ",c("F","Y","W","C","M"),sep=""),
	paste("Number of ASIs of type [Aromatic-Sulphur]: ",c("FC","YC","WC","FM","YM","WM"),sep=""),
	"Percentage of buried ASIs","Percentage of exposed ASIs","Percentage of ASIs that are not part of any network (Isolated)",
	"Total number of ASAI-networks","Network details [SIZE_STRENGTH:NUMBER OF SUCH NETWORK] (SIZE:Number of residue in the network, STRENGTH:Number of interaction in the network)"
	);
	write2Html.tableHeader(c("Filename","Total","IntraCh","InterCh","F","Y","W","C","M","FC","YC","WC","FM","YM","WM","B","E","Iso","T.Net","Net_Details"),filepath,append=TRUE,tipvec);
}

