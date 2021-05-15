disulfide=function(pdb="",low=2.2)
{
cys=atom.select(pdb,resid="CYS",elety="SG",verbose=FALSE)$atom
# If less than 2 cys then u return sorry
if(length(cys)<2){return("SORRY");}
cysXYZmat=matrix(as.numeric(pdb$atom[cys,c("x","y","z")]),ncol=3);
cysdist=dist.xyz(cysXYZmat);
dist.Cutoff=matrix(which(cysdist<=low & cysdist >0,arr.ind=TRUE),ncol=2);
dist.Cutoff=matrix(dist.Cutoff[which(dist.Cutoff[,1]<dist.Cutoff[,2]),],ncol=2)
# If no cys pair satisfies disulfide bond then return sorry
if(nrow(dist.Cutoff)==0){return("SORRY");}

cystemp=matrix(pdb$atom[cys[dist.Cutoff[,1]],c("chain","resno","insert")],ncol=3);
cystemp=cbind(cystemp,matrix(pdb$atom[cys[dist.Cutoff[,2]],c("chain","resno","insert")],ncol=3));

cystemp[which(is.na(cystemp[,3])),3]="";
cystemp[,2]=paste(cystemp[,2],cystemp[,3],sep="");
cystemp[which(is.na(cystemp[,6])),6]="";
cystemp[,5]=paste(cystemp[,5],cystemp[,6],sep="");
cystemp=matrix(cystemp[,c(1,2,4,5)],ncol=4)
colnames(cystemp)=c("chain1","resno1","chain2","resno2");
return(cystemp);
}

disulfide_addDSSP=function(DSSP="",disul="")
{
#change ~ chain to space, we have called Nochain() function before
disul[which(disul[,"chain1"]=="~"),"chain1"]=" ";
disul[which(disul[,"chain2"]=="~"),"chain2"]=" ";
cys1ChRes=paste(disul[,"chain1"],disul[,"resno1"],sep="")
cys2ChRes=paste(disul[,"chain2"],disul[,"resno2"],sep="")

DSSPChRes=paste(DSSP$cha,DSSP$res,sep="") #combine chain and resno of DSSP object
cys1SS="";cys2SS="";
	for(i in 1:length(cys1ChRes))
	{
		index=which(DSSPChRes==cys1ChRes[i]);
		if(length(index)==0){cys1SS=cbind(cys1SS,"-");}else{
		cys1SS=cbind(cys1SS,DSSP$ss[index]);
		}
		index=which(DSSPChRes==cys2ChRes[i]);
		if(length(index)==0){cys2SS=cbind(cys2SS,"-");}else{
		cys2SS=cbind(cys2SS,DSSP$ss[index]);
		}
	}
	cys1SS=cys1SS[-1];
	cys2SS=cys2SS[-1];
disul[which(disul[,"chain1"]==" "),"chain1"]="~";
disul[which(disul[,"chain2"]==" "),"chain2"]="~";

disul=cbind(disul,resid1_SS=cys1SS);
disul=cbind(disul,resid2_SS=cys2SS);
return(disul);
}

disulfide_addNaccess=function(nacc="",disul="")
{
chres=paste(nacc$asa[,1],nacc$asa[,2],sep=""); #ChainResnoResid format
#change ~ chain to space, we have called Nochain() function before
disul[which(disul[,"chain1"]=="~"),"chain1"]=" ";
disul[which(disul[,"chain2"]=="~"),"chain2"]=" ";

cys1ChRes=paste(disul[,"chain1"],disul[,"resno1"],sep="")
cys2ChRes=paste(disul[,"chain2"],disul[,"resno2"],sep="")
cys1Acc=0;cys2Acc=0;
	for(i in 1:length(cys1ChRes))
	{
		index=which(chres==cys1ChRes[i]);
		if(length(index)==0){cys1Acc=cbind(cys1Acc,"-");}else{
		cys1Acc=cbind(cys1Acc,nacc$asa[index,5]); #asa_r_rel ASA residue relative
		}
		index=which(chres==cys2ChRes[i]);
		if(length(index)==0){cys2Acc=cbind(cys2Acc,"-");}else{
		cys2Acc=cbind(cys2Acc,nacc$asa[index,5]);
		}
	}
	cys1Acc=cys1Acc[-1];
	cys2Acc=cys2Acc[-1];

	disul[which(disul[,"chain1"]==" "),"chain1"]="~";
	disul[which(disul[,"chain2"]==" "),"chain2"]="~";

disul=cbind(disul,resid1_Acc=cys1Acc);
disul=cbind(disul,resid2_Acc=cys2Acc);
return(disul);		
}

disulfide_summary=function(dis,acclow=20)
{
	tot=nrow(dis);
	intra=sum(dis[,1]==dis[,3])
	inter=sum(dis[,1]!=dis[,3])
	bur="";exp="";
	# If nacc is failed then dis[,5] and dis[,6] will be "-" so u cant calculate Bur and Exp
	if(sum(dis[,5:6]!="-")>0)
	{
	bur=sum(((as.numeric(dis[,5])+as.numeric(dis[,6]))/2)<=acclow); #cys1,cys2 avg acc if <=20 then buried
	exp=tot-bur;
	}else{
	bur="-";exp="-";
	}
	# loop lemgth can only be calculated for intadisulfide bonds only
	intra.ind=which(dis[,1]==dis[,3])
	loopstat="";
	if(length(intra.ind)>0)
	{
	intra.looplen=(as.numeric(dis[intra.ind,4])-as.numeric(dis[intra.ind,2]))-1; #If C1: 10, C2: 12, then loop length is 1. Exclude the Cys positions
	# range are (0,10] (10,20] (20,30] (30,40] (40,50] (50,1000]
	loopstat=apply(matrix(c(0,10,10,20,20,30,30,40,40,50,50,100000),ncol=2,byrow=TRUE),1,function(x){sum(intra.looplen>x[1] & intra.looplen<=x[2])});
	}else{
	# If no intra disulfide bonds then no question of loop length. so
	loopstat=rep(0,times=6);
	}
	names(loopstat)=c("len0_10","len10_20","len20_30","len30_40","len40_50","len>50");
	pp="";nn="";pn=""; #P: periodic secondary structure like H/G/I/E else T/C/S/B are non-periodic
	# Total number of disulfide connecting 2 periodic sec str / NN/NP
	 
	# If DSSP is true then only pp,pn,nn
	if(sum(dis[,7:8]!="-")>0)
	{
	period1=dis[,7];
	period1=gsub("[HGIE]","P",period1,perl=TRUE)
	period1=gsub("[TCSB]","N",period1,perl=TRUE)
	period2=dis[,8];
	period2=gsub("[HGIE]","P",period2,perl=TRUE)
	period2=gsub("[TCSB]","N",period2,perl=TRUE)
	paste(period1,period2,collapse="")
	pp=sum(paste(period1,period2,sep="")=="PP");
	nn=sum(paste(period1,period2,sep="")=="NN");
	pn=sum(period1!=period2);
	}else{
	pp="-";nn="-";pn="-";
	}
	return(c(Total=tot,intraNo=intra,interNo=inter,totBur=bur,totExp=exp,loopstat,PP=pp,PN=pn,NN=nn));
}

getDisulfide=function(chres="",dis="")
{
	chain=chres[1];
	resno=chres[2];
	if(is.na(chain)){chain="~"}
	query=paste(chain,resno,sep="",collapse="");
	dis_res1=paste(dis[,1],dis[,2],sep="");
	dis_res2=paste(dis[,3],dis[,4],sep="");
	dis.ind=which(dis_res1==query | dis_res2==query);
	if(length(dis.ind)==0){return("SORRY");} # If not found return SORRY else an matrix
	dismat=dis[dis.ind,];
	dismat=matrix(dismat,ncol=8);
	colnames(dismat)=colnames(dis);
	return(dismat)
}

callDis=function(pdb,DSSP="",nacc="",status.DSSP=0,status.nacc=0)
{
			# This is the workflow that is used in Pattern_analysis code
			dis=disulfide(pdb);
			status.dis=0;
			if(is.matrix(dis))
			{
				status.dis=1;
				dis=Nochain(dis,1,3);
				if(status.nacc){dis=disulfide_addNaccess(nacc,dis);}else{dis=cbind(dis,resid1_Acc="-");dis=cbind(dis,resid2_Acc="-");}
				if(status.DSSP){dis=disulfide_addDSSP(DSSP,dis);}else{	dis=cbind(dis,resid1_SS="-");dis=cbind(dis,resid2_SS="-");}
			}
return(list(dis=dis,status.dis=status.dis));				
}

checkDis=function(patmat,dis)
{
	# This is the workflow that is used in Pattern_analysis code
	dismat=do.call(rbind,lapply(apply(patmat,1,getDisulfide,dis),matrix,ncol=8,byrow=FALSE));
	sorry.ind=which(apply(dismat,1,function(x){sum(x=="SORRY")>0}));
	if(length(sorry.ind)>0){dismat=dismat[-sorry.ind,];}
	dismat=matrix(dismat,ncol=8);
	return(dismat);
					
}

print_HTMLdis=function(dis,filepath,append=FALSE,html=TRUE)
{
	if(html){write2Html("<html>\n",filepath,append);}else{write2Html("",filepath,append);}
	write2Html("<br><br><b><u>Disulfide bonds</u><br><br></b>\n",filepath,TRUE);
#	write2Html("Chain1,Resno1,Resid1: Corresponds to Acidic residue details<br>\n",filepath,TRUE);
#	write2Html("Chain2,Resno2,Resid2: Corresponds to Basic residue details<br>\n",filepath,TRUE);
#	write2Html("Resid1_Acc and Resid2_Acc corresponds to relative ASA of acidic and basic residue respectively<br>\n",filepath,TRUE);
#	write2Html("Resid1_SS and Resid2_SS corresponds to secondary structure of the acidic and basic residues respectively<br><br><br>\n",filepath,TRUE);

	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr><th colspan=2 title=\"First Cys partner\"> Cys1 </th> <th colspan=2 title=\"Second Cys partner\"> Cys2 </th> <th colspan=2 title=\"Relative solvent accessibility of the residue\"> ASA </th> <th colspan=2 title=\"Secondary structure of residue (DSSP notation)\"> Secstr </th> </tr>",filepath,append=TRUE);
	tipvec=c("Chain","Residue number","Chain","Residue number","","","","");
	write2Html.tableHeader(c("Chain","Res1_no","Chain","Res2_no","Cys1","Cys2","Cys1","Cys2"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(dis,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	if(html){write2Html("</html>\n",filepath,TRUE);}else{write2Html("",filepath,TRUE);}
}

print_HTMLdis_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Disulfide bond summary</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("Filename","Total number of disulfide bond present","Number of intra chain disulfide bond",
	"Number of inter chain disulfide bond","Number of disulfide bond buried","Number of disulfide bond exposed",
	paste("Number of disulfide bond with loop size (sequence gap between two cys) between (",c("0-10","10-20","20-30","30-40","40-50"),"]",sep=""),
	"Number of disulfide bond with loop size (sequence gap between two cys)>50",
	"Number of disulfide bond connecting two periodic secondary structure [DSSP notation:HGIE]",
	"Number of disulfide bond connecting a periodic [DSSP notation:HGIE] and a non-periodic [DSSP notation:TCSB] secondary structure",
	"Number of disulfide bond connecting two non-periodic secondary structure [DSSP notation:TCSB]"
	);
	write2Html.tableHeader(c("Filename","Total","IntraCh","InterCh","B","E","0_10","10_20","20_30","30_40","40_50",">50","PP","PN","NN"),filepath,append=TRUE,tipvec);
}

