helixDipole=function(DSSP="")
{

hinfo=""; #stores final information as matrix
hstat=c(0,0,0,0,0,0,0,0,0,0); #N-2 to N+2 and C-2 to C+2, total 10
hstatLabel=c("N-2","N-1","N","N+1","N+2","C-2","C-1","C","C+1","C+2"); 
		
# If at least one helix is there in the protein then only go for checking helix dipole
if(length(DSSP$helix$start)>0)
{

# Suppose DSSP$helix$start=101,122,123 and DSSP$helix$chain=A,B,C then Hst_chres=A101 B122 C123
# str(DSSP) will return DSSP$ch: all chain info, DSSP$res: all resno, we join them in DSSP_chres A1 A2 A3....

	Hst_chres=paste(DSSP$helix$chain,DSSP$helix$start,sep="") 
	DSSP_chres=paste(DSSP$ch,DSSP$res,sep="") 
	ext=2; #external span of helix i.e N-ext to N and C to C+ext
	int=0; #internal span to helix i.e N to N+int and C to C-int
	
	
	# loop over each helix
	# To store no of helix where Ndipole and Cdipole stab present.
	
	
	for(hno in 1:length(DSSP$helix$start))
	{
	#print(hno)
	hlen=as.numeric(DSSP$helix$length[hno]); #helix length, [0-5]: int=0,  [6,7]: int=1, [8,inf]:int=2
	
	if(hlen<6)
	{
	int=0;
	NcapsLabel=c("N-2","N-1","N"); #Ncap labels
	CcapsLabel=c("C","C+1","C+2"); #Ccap labels
	
	}else{
		if(hlen>5 & hlen<8 )
		{
			int=1;
			NcapsLabel=c("N-2","N-1","N","N+1");
			CcapsLabel=c("C-1","C","C+1","C+2");
			}else{
				int=2;
				NcapsLabel=c("N-2","N-1","N","N+1","N+2");
				CcapsLabel=c("C-2","C-1","C","C+1","C+2");
			}
	}
	
		DSSP.Hstindex=which(DSSP_chres==Hst_chres[hno]); # lets consider helix1 (hno=1), st=105, chain=A, then Hst_chres[hno]="A105". which DSSP.chres is "A105". the index we store as same as DSSP.Hstindex
		DSSP.Hendindex=DSSP.Hstindex+hlen-1; # DSSP.Hstindex +helix length -1 is end index
	#	print(DSSP$aa[(DSSP.Hstindex-ext):(DSSP.Hendindex+ext)]);
		Ncaps=DSSP$aa[(DSSP.Hstindex-ext):(DSSP.Hstindex+int)] # by knowing DSSP.Hstindex, we can easily fetch Ncaps and Ccaps
		Ccaps=DSSP$aa[(DSSP.Hendindex-int):(DSSP.Hendindex+ext)]
	
		NcapsDE=which(Ncaps=="D" | Ncaps == "E")
		CcapsRK=which(Ccaps=="R" | Ccaps == "K")
		flag=FALSE; # if for a helix dipole stabilization is there, then only print.
		if(length(NcapsDE)>0 | length(CcapsRK)>0){flag=TRUE;}
		
		if(flag)
		{
			hinfotemp=c(DSSP$helix$chain[hno],DSSP$helix$start[hno],DSSP$helix$end[hno],DSSP$helix$type[hno],hlen);
			if(length(NcapsDE)>0)
			{	
				hinfotemp=c(hinfotemp,paste(NcapsLabel[NcapsDE],Ncaps[NcapsDE],sep=":",collapse=","));
				#count of positions of dipole stabilization
				for(ncaptemp in NcapsLabel[NcapsDE]){ncaptemp.index=which(hstatLabel==ncaptemp);hstat[ncaptemp.index]=hstat[ncaptemp.index]+1}
				#print(NcapsLabel[NcapsDE]);
			}else{hinfotemp=c(hinfotemp,"-");}
		
			if(length(CcapsRK)>0)
			{
				hinfotemp=c(hinfotemp,paste(CcapsLabel[CcapsRK],Ccaps[CcapsRK],sep=":",collapse=","))
				for(ccaptemp in CcapsLabel[CcapsRK]){ccaptemp.index=which(hstatLabel==ccaptemp);hstat[ccaptemp.index]=hstat[ccaptemp.index]+1}
				#print(CcapsLabel[CcapsRK]);
			}else{hinfotemp=c(hinfotemp,"-");}
		hinfo=rbind(hinfo,hinfotemp);
		}
	}
	# If no helixdipole then hinfo="";
	if(is.matrix(hinfo))
	{
	hinfo=hinfo[-1,];
	hinfo=matrix(hinfo,ncol=7);
	colnames(hinfo)=c("chain","start","end","type","len","NcapAcid","CcapBase")
	#rownames(hinfo)=1:length(hinfo[,1])
	N_cap_StabHelices=sum(hinfo[,7]=="-" & hinfo[,6]!="-")
	C_cap_StabHelices=sum(hinfo[,7]!="-" & hinfo[,6]=="-")	
	NC_cap_StabHelices=sum(hinfo[,7]!="-" & hinfo[,6]!="-")
	tot_StabHelices=N_cap_StabHelices+C_cap_StabHelices+NC_cap_StabHelices
	return(list(hdipo=hinfo,tot_helices=length(DSSP$helix$start),tot_StabHelices=tot_StabHelices,N_cap_StabHelices=N_cap_StabHelices,C_cap_StabHelices=C_cap_StabHelices,NC_cap_StabHelices=NC_cap_StabHelices,hstat=rbind(hstatLabel,hstat)));
	}else{
	# Dssp success but no stabilized helices then we print 0 instead of "-" to indicate that actually it is not present rather than dssp failure
	return(list(hdipo=hinfo,tot_helices=length(DSSP$helix$start),tot_StabHelices=0,N_cap_StabHelices=0,C_cap_StabHelices=0,NC_cap_StabHelices=0,hstat=rbind(hstatLabel,hstat)));
	}
	
	
}else{
# If dssp fails then instead of 0 we print "-"
return(list(hdipo=hinfo,tot_helices="-",tot_StabHelices="-",N_cap_StabHelices="-",C_cap_StabHelices="-",NC_cap_StabHelices="-",hstat=rbind(hstatLabel,rep("-",times=10))));
}

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
	tipvec=c("Chain","Start of helix","End residue of helix","Type of helix [DSSP notation H/G/I]","Helix length","Details of dipole stabilization at N-terminal of helix [Format: Cap Position-aa Residue]",
	"Details of dipole stabilization at C-terminal of helix [Format: Cap Position-aa Residue]"
	);
	write2Html.tableHeader(c("Chain","H.start","H.end","H.type","Length","NcapAcid","CcapBase"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(hdipo,filepath,TRUE,"<td>");
	write2Html("\n</table>\n</html>",filepath,TRUE);
}

print_HTMLhdipo_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Helix dipole stabilization profile</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("PDB filename","Total number of Helix","Total number of helices whose dipole are stabilized",
	"Number of helices whose dipoles are stabilized only at N-terminal side","Number of helices whose dipoles are stabilized only at C-terminal side",
	"Number of helices whose dipoles are stabilized both at N- and C-terminal side",
	paste("Number of helices in which the dipole stabilization is seen at ",c("N-2","N-1","N","N+1","N+2","C-2","C-1","C","C+1","C+2")," position",sep="")
	);
	write2Html.tableHeader(c("Filename","T.Helix","St.Helix","N.St.Helix","C.St.Helix","NC.St.Helix","St_N-2","St_N-1","St_N","St_N+1","St_N+2","St_C-2","St_C-1","St_C","St_C+1","St_C+2"),filepath,append=TRUE,tipvec);
}


