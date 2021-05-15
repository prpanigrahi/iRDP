proline=function(pdb="",DSSP="",nacc="",acclow=20,prefix="")
{
	pro.ind=atom.select(pdb,resid="PRO",elety="CA",verbose=FALSE)$atom;
	totP=length(pro.ind); #total prolines
	p.prop=rep("-",times=6);
	names(p.prop)=c("hel","str","turn","coil","bur","expo");
	totbtp=0; #total beta-tuen proline
	totBT2P=0; # total beta-turn 2nd proline
	p.ncap=""; # Ncap proline helix details
	p.ncapno=0; # Total no of ncap prolie residues
	
	BT=""; # BT details
	BT2P=""; # BT2P details
	
	if(is.list(DSSP))
	{
		dsspchres=paste(DSSP$cha,DSSP$res,sep="");
   	p.ind=which(DSSP$aa=="P");
		p.distr=tapply(DSSP$ss[p.ind],factor(DSSP$ss[p.ind]),length)

		p.prop["hel"]=sum(p.distr["H"],p.distr["G"],p.distr["I"],na.rm=T) # Total prolines in helices
		p.prop["str"]=sum(p.distr["B"],p.distr["E"],na.rm=T) # Total prolines in strands
		p.prop["turn"]=p.distr["T"] # Total prolines in turns
		p.prop["coil"]=sum(p.distr["C"],p.distr["S"],na.rm=T) # Total prolines in bands and coils
		
		p.ssdet= rbind(
		cbind("Helix",paste(dsspchres[p.ind][which(DSSP$ss[p.ind] %in% c("H","G","I"))],collapse=",")),
		cbind("Strand",paste(dsspchres[p.ind][which(DSSP$ss[p.ind] %in% c("B","E"))],collapse=",")),
		cbind("Turn",paste(dsspchres[p.ind][which(DSSP$ss[p.ind] %in% c("T"))],collapse=",")),
		cbind("Coil",paste(dsspchres[p.ind][which(DSSP$ss[p.ind] %in% c("S","C"))],collapse=","))
		);
		
		# If any helix is present then only proceeed
		# Even if no helix is present dssp$helix is a list, so we cant apply function is.list(DSSP4helix)
		if(length(DSSP$helix$start)>0){
			ncap.ind=apply(matrix(paste(DSSP$helix$chain,DSSP$helix$start,sep=""),ncol=1),1,function(x){which(paste(DSSP$cha,DSSP$res,sep="")==x)}) #DSSP index of all Ncap helix position
			p.ncapind=ncap.ind[which(DSSP$aa[ncap.ind]=="P")]  #DSSP index of only Ncap helix position where N-cap is proline
			if(length(p.ncapind)>0){ # Proceed if any Ncap proline is present then only give details
			p.ncap=matrix(cbind(chain=DSSP$cha[p.ncapind],resno=DSSP$res[p.ncapind],resid=DSSP$aa[p.ncapind],ss=DSSP$ss[p.ncapind],acc=DSSP$acc[p.ncapind]),ncol=5); #a matrix
			p.ncapno=nrow(p.ncap);
			}
		}
	
	}
	if(is.list(nacc))
	{	
		naccChres=paste(nacc$asa[,1],nacc$asa[,2],sep="");
		insert.res=pdb$atom[pro.ind,"resno"];
		insert.ins=pdb$atom[pro.ind,"insert"];
		insert.ind=which(!is.na(insert.ins))
		insert.res[insert.ind]=paste(insert.res[insert.ind],insert.ins[insert.ind],sep="")
		nacc.ind=apply(matrix(paste(pdb$atom[pro.ind,"chain"],insert.res,sep=""),ncol=1),1,function(x){which(naccChres==x)})
		p.prop["bur"]=sum(nacc$asa[nacc.ind,5]<=acclow) # Total buried prolines
		p.prop["expo"]=sum(nacc$asa[nacc.ind,5]>acclow) # Total exposed prolines
	}
	p.prop[is.na(p.prop)]=0;
	
	BT=runPROMOTIF(pdb,promotifpath,prefix);
	if(is.matrix(BT))
	{
	totbtp=length(grep("P",BT[,4])); # Total no BT where proline is present irrespective of any position
	BT2P.ind=which(substring(BT[,4],2,2)=="P"); # Row index in BT where 2nd position is PRO
	totBT2P=length(BT2P.ind); # Total no of 2nd-beta turn proline
	if(totBT2P>0)
	{
	BT2P=matrix(BT[BT2P.ind,],ncol=13);
	}
		if(totBT2P>0)
		{
			BT2P.ss="";
			if(is.list(DSSP))
			{
			# generate Bt2P detilas
			BT2P.chres=paste(BT[BT2P.ind,1],BT[BT2P.ind,2],sep="") # "A183" "B52"  "B86"  "B204" "B495" "B512"
			DSSP.chres=paste(DSSP$cha,DSSP$res,sep=""); # "A183" "B52"  "B86"
				for(i in 1:length(BT2P.chres))
				{
				#print(i);	
				DSSP.BT1P.ind=which(DSSP.chres==BT2P.chres[i]) # In DSSP where the beta-turn first position is located, its index.
				if(length(DSSP.BT1P.ind)>0){BT2P.ss=c(BT2P.ss,paste(DSSP$ss[DSSP.BT1P.ind:(DSSP.BT1P.ind+3)],sep="",collapse=""))}else{BT2P.ss=c(BT2P.ss,"----");}
				}
			BT2P.ss=BT2P.ss[-1];
			BT2P=cbind(BT2P,secstr=BT2P.ss);
			}else{
			BT2P=cbind(BT2P,secstr="----");
			}

			BT2P.acc=""; # stores buried or exposed string for each beta-turn
			if(is.list(nacc))
			{
			BT2P.chres=paste(BT[BT2P.ind,1],BT[BT2P.ind,2],sep="") # "A183" "B52"  "B86"  "B204" "B495" "B512"
			nacc.chres=paste(nacc$asa[,1],nacc$asa[,2],sep="");
				for(i in 1:length(BT2P.chres))
				{	
				nacc.BT1P.ind=which(nacc.chres==BT2P.chres[i]) # In DSSP where the beta-turn first position is located, its index.
				if(length(nacc.BT1P.ind)>0){BT2P.acc=c(BT2P.acc,paste(apply(matrix(nacc$asa[nacc.BT1P.ind:(nacc.BT1P.ind+3),5],ncol=1),1,function(x){if(x <= acclow){return("B")}else{return("E")}; }),sep="",collapse=""))}else{BT2P.acc=c(BT2P.acc,"----");}
				}
			BT2P.acc=BT2P.acc[-1];
			BT2P=cbind(BT2P,BurExp=BT2P.acc);
			}else{
			BT2P=cbind(BT2P,BurExp="----");
			}
		}
	}
return(list(prop=c(tot=totP,p.prop,TotBtP=totbtp,TotBt2P=totBT2P,TotNcapP=p.ncapno),Bt2P=BT2P,NcapP=p.ncap,ssdet=p.ssdet));	
}

proline_old=function(pdb,DSSP="",SDM,FOLDX,IMUT)
{
	totP=length(atom.select(pdb,resid="PRO",elety="CA",verbose=FALSE)$atom); #total prolines
	p.prop="";
	if(is.list(DSSP))
	{
		p.ind=which(DSSP$aa=="P");
		p.distr=tapply(DSSP$ss[p.ind],factor(DSSP$ss[p.ind]),length)

		p.hel=sum(p.distr["H"],p.distr["G"],p.distr["I"],na.rm=T) # Total prolines in helices
		p.str=sum(p.distr["B"],p.distr["E"],na.rm=T) # Total prolines in strands
		p.turn=p.distr["T"] # Total prolines in turns
		p.coil=sum(p.distr["C"],p.distr["S"],na.rm=T) # Total prolines in bands and coils
	
		p.bur=sum(DSSP$acc[p.ind]<=20) # Total buried prolines
		p.expo=sum(DSSP$acc[p.ind]>20) # Total exposed prolines
		p.prop=c(p.hel,p.str,p.turn,p.coil,p.bur,p.expo);
	}else{
	p.prop=rep("",times=6);
	}
	
	BT=runPROMOTIF(pdb,promotifpath)
	# Total no BT where proline is present irrespective of any position
	if(!is.matrix(BT))
	{
	
	totbtp=length(grep("P",BT[,4]));
	
	BT2P.ind=which(substring(BT[,4],2,2)=="P"); # Row index in BT where 2nd position is PRO
	totBT2P=length(BT2P.ind); # Total no of 2nd-beta turn proline
	
	BT2P.chres=paste(BT[BT2P.ind,1],as.numeric(BT[BT2P.ind,2]),sep="") # "A183" "B52"  "B86"  "B204" "B495" "B512"
	DSSP.chres=paste(DSSP$cha,DSSP$res,sep=""); # "A183" "B52"  "B86"

	BT2P.acc=""; # stores buried or exposed string for each beta-turn
	BT2P.sdm="";
	BT2P.fold="";
	BT2P.imut="";
	BT2P.ss="";
	#DSSP index of all Ncap helix position
	ncap.ind=apply(matrix(paste(DSSP$helix$chain,DSSP$helix$start,sep=""),ncol=1),1,function(x){which(paste(DSSP$cha,DSSP$res,sep="")==x)})
	#DSSP index of only Ncap helix position
	p.ncapind=ncap.ind[which(DSSP$aa[ncap.ind]=="P")]
	p.ncap=cbind(chain=DSSP$cha[p.ncapind],resno=DSSP$res[p.ncapind],resid=DSSP$aa[p.ncapind],ss=DSSP$ss[p.ncapind],acc=DSSP$acc[p.ncapind]);
	
		
	if(SDM)
	{
		write.pdb(pdb, file = sdmpdbfile);
		
	}
	if(FOLDX)
	{
		write.pdb(pdb, file = foldxtemppdb);
	}
	
	for(i in 1:length(BT2P.chres))
	{
		DSSP.BT1P.ind=which(DSSP.chres==BT2P.chres[i]) # In DSSP where the beta-turn first position is located, its index.
		BT2P.acc=c(BT2P.acc,paste(apply(as.matrix(DSSP$acc[DSSP.BT1P.ind:(DSSP.BT1P.ind+3)]<=20),1,function(x){if(x){return("B")}else{return("E")}}),collapse=""));
		BT2P.ss=c(BT2P.ss,paste(DSSP$ss[DSSP.BT1P.ind:(DSSP.BT1P.ind+3)],sep="",collapse=""))
		if(SDM)
		{
		print(paste("perl SDM.pl ", sdmpdbfile," ",DSSP$cha[DSSP.BT1P.ind+1]," G ", DSSP$res[DSSP.BT1P.ind+1],sep=""))
		BT2P.sdm=rbind(BT2P.sdm,checkSDM(sdmpdbfile,DSSP$cha[DSSP.BT1P.ind+1],"G",DSSP$res[DSSP.BT1P.ind+1]));
		}
		if(FOLDX)
		{
		mut=paste(DSSP$aa[DSSP.BT1P.ind+1],DSSP$cha[DSSP.BT1P.ind+1],DSSP$res[DSSP.BT1P.ind+1],"G;",collapse="",sep="");
		print(mut);
		stab=runFoldx(foldxtemppdb,mut,foldxpath,foldxrotamer);
		BT2P.fold=rbind(BT2P.fold,stab);
		}
		if(IMUT)
		{
			print("running iMut mut");
			BT2P.imut=rbind(BT2P.imut,iMut(pdb,DSSP$cha[DSSP.BT1P.ind+1],DSSP$res[DSSP.BT1P.ind+1],"G"));
		}
	}
	
	BT2P.acc=BT2P.acc[-1];
	BT2P.ss=BT2P.ss[-1];
	BT2P=cbind(BT[BT2P.ind,],BurExp=BT2P.acc,secstr=BT2P.ss);
	
	if(SDM){BT2P.sdm=BT2P.sdm[-1,];BT2P=cbind(BT2P,BT2P.sdm);}
	if(FOLDX){BT2P.fold=BT2P.fold[-1,];BT2P=cbind(BT2P,BT2P.fold);}
	if(IMUT){BT2P.imut=BT2P.imut[-1,];BT2P=cbind(BT2P,BT2P.imut);}
	return(list(prop=c(tot=totP,hel=p.hel,str=p.str,turn=p.turn,coil=p.coil,bur=p.bur,expo=p.expo,totbtp=totbtp,bt2p=nrow(BT2P),ncapp=p.ncap),ncapP=p.ncap,BT2P=BT2P));
}
}

#library("bio3d");
#source("dssp.r");
#source("PATH.r");
#source("bturn.r");
#source("SDM.r");
#source("foldx.r");
#source("iMut.r");
#pdb=read.pdb("1gk9.pdb");
#DSSP=dssp_new(pdb,exepath=dssppath ,infile="temp.pdb",outfile="temp.dssp");

#proline(pdb,DSSP,0,0,0);

print_HTMLproline=function(pro,filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u> Proline distribution</u><br><br></b>\n\n",filepath,TRUE);
#	write2Html("Chain1,Resno1,Resid1: Corresponds to Acidic residue details<br>\n",filepath,TRUE);
#	write2Html("Chain2,Resno2,Resid2: Corresponds to Basic residue details<br>\n",filepath,TRUE);
#	write2Html("Resid1_Acc and Resid2_Acc corresponds to relative ASA of acidic and basic residue respectively<br>\n",filepath,TRUE);
#	write2Html("Resid1_SS and Resid2_SS corresponds to secondary structure of the acidic and basic residues respectively<br><br><br>\n",filepath,TRUE);
	
	write2Html("\n<br><br> <b><u>Proline residue distributon in each type of secondary structure<u><b><br></table>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("Secondary structure type","The details of residue [Format: Chain Resno]");
	write2Html.tableHeader(c("Secondary structure","Residue details"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(pro$ssdet,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);	
	
	if(is.matrix(pro$Bt2P))
	{
	write2Html("\n<br><br> <b><u>Proline at 2nd position of beta-turns<u><b><br></table>\n",filepath,TRUE);	
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("Chain","Start residue number of beta-turn","End residue number of beta-turn","Beta-turn sequence","Beta-turn type",
	"phi of residue R2 (Promotif notation)","psi of residue R2 (Promotif notation)","chi1 of residue R2 (Promotif notation)",
	"phi of residue R3 (Promotif notation)","psi of residue R3 (Promotif notation)","chi1 of residue R3 (Promotif notation)",
	"CA distance between R1 and R4 residue","Does hydrogen bond exists? Y:Yes, N:No","Secondary structure of beta-turn residues [DSSP notation]",
	"Solvent accessibility of beta-turn residues [B:Buried,E:Exposed]"
	);
	write2Html.tableHeader(c("Chain","Start","End","Sequence","BT_type","R2_phi","R2_psi","R2_chi1","R3_phi","R3_psi","R3_chi1","R1R4_CAdist","Hbond","SS","B/E"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(pro$Bt2P,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);	
	}
	write2Html("\n<br><br>\n",filepath,TRUE);
	
	if(is.matrix(pro$NcapP))
	{
	write2Html("\n<br><br> <b><u>Proline at N-cap of helices <u><b><br></table>\n",filepath,TRUE);	
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("Chain","Start of helix","Residue at start/N-cap position of helix","Helix type [DSSP notation H/G/I]","Solvent accessibility of N-cap position");
	write2Html.tableHeader(c("Chain","Start","Res.ID","Type","Acc"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(pro$NcapP,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);	
	}
	write2Html("\n</html>\n",filepath,TRUE);

}

print_HTMLproline_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Proline profile</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("PDB filename","Total number of prolines",
	paste("Number of prolines in ",c("Helix","Strand","Turn","Coil"), " (DSSP notation: ",c("HGI","BE","T","SC"),")",sep=""),
	"Number of prolines buried","Number of prolines exposed",
	"Total number of prolines in beta-turns","Number of prolines at 2nd-position of beta-turns","Number of prolines at N-cap of helices (DSSP notation: HGI)"
	);
	write2Html.tableHeader(c("Filename","T.Proline","Helix","Strand","Turn","Coil","B","E","T.BtP","T.Bt2P","T.NcapP"),filepath,append=TRUE,tipvec);
}


