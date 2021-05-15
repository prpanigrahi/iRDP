#intprofile.r

#########################################################################################################################################################
calc_interaction=function(pdb,DSSP="",nacc="",status.DSSP=0,status.nacc=0,calc_ionic=0,calc_aroaro=0,calc_aros=0,calc_catpi=0,calc_hphob=0,calc_hbond=0,calc_disul=0,ioncut=6.0,aroarolow,aroarohigh,aroarodihed,usedihed,aroscut,catpicut,hpcut)
{
# This function is the function which will generate list of all interaction.
# We need to pass a pdb object it will generate a list of 11 elements
# 3: ionList$ionpair,ionList$status.ionic,ionList$dionpair
# 3: aroaroList$aro,aroaroList$daropair,aroaroList$status.aroaro,
# 3: arosList$aros,arosList$darospair,arosList$status.aros,
# 2: hbList$hb,hbList$status.hb,
# 2: disList$dis,disList$status.dis,
# 3: catpiList$catpimat,catpiList$dcatpimat,catpiList$status.catpi,
# 2: hpList$hp,hpList$status.hp
# Total list of 18 element.
# ionList$ionpair can be an matrix if ionpair found else "SORRY" so we can apply is.matrix() function this is what the ionic() function returns
# ionList$status.ionic either 0/1
# ionList$dionpair: a dionpair object
ionList="";aroaroList="";arosList="";hbList="";disList="";catpiList="";hpList="";
			if(calc_ionic){ionList=callIonic(pdb,DSSP,nacc,status.DSSP,status.nacc,ioncut);}else{ionList=list(ionpair="",status.ionic=0,dionpair="");} # Calculate Ionpair
			if(calc_aroaro){aroaroList=callAroAro(pdb,DSSP,nacc,status.DSSP,status.nacc,aroarolow,aroarohigh,aroarodihed,usedihed);}else{aroaroList=list(aro="",status.aroaro=0,daropair="");} #  Calculate Aropair
			if(calc_aros){arosList=callAroS(pdb,DSSP,nacc,status.DSSP,status.nacc,aroscut);}else{arosList=list(aros="",status.aros=0,darospair="");} # Calculate  AroS
			if(calc_hbond){hbList=callHbond(pdb);}else{hbList=list(hb="",status.hb=0)} # Calculate  Hbond
			
			if(calc_disul){disList=callDis(pdb,DSSP,nacc,status.DSSP,status.nacc);}else{disList=list(dis="",status.dis=0);} # Calculate  Disulfide
			if(calc_catpi){catpiList=callCatpi(pdb,DSSP,nacc,status.DSSP,status.nacc,catpicut);}else{catpiList=list(catpimat="",status.catpi=0,dcatpimat="");} # Calculate  Catpi
			if(calc_hphob){hpList=callHphob(pdb,DSSP,nacc,status.DSSP,status.nacc,hpcut);}else{hpList=list(hp="",status.hp=0);} #  Calculate Hydrophobic
			
			return(list( 
			ionpair=ionList$ionpair,status.ionic=ionList$status.ionic,dionpair=ionList$dionpair,
			aro=aroaroList$aro,daropair=aroaroList$daropair,status.aroaro=aroaroList$status.aroaro,
			aros=arosList$aros,darospair=arosList$darospair,status.aros=arosList$status.aros,
			hb=hbList$hb,status.hb=hbList$status.hb,dis=disList$dis,status.dis=disList$status.dis,
			catpimat=catpiList$catpimat,dcatpimat=catpiList$dcatpimat,status.catpi=catpiList$status.catpi,hp=hpList$hp,status.hp=hpList$status.hp));
}


calc_intProfile=function(patmat,intList,acclow=20)
{
# This is a function to estimate local interaction profile of a residue
# You pass residue information as patmat and it will return you a list (intprofile) of local interaction
# This list can be used to generate profile (count of interactions)
# intprofile:a list of 11 element.
# ip,ipnet,ap,apnet,as,asnet,hb,dis,catpi,catpinet,hp
# The values can be either 0 if that residue doesnot invovle in that type-interaction or in that type-network
# Else it returns a matrix (local ionpair matrix, local ion network matrix)	
	
#Input: patmat, it is a 1*3 matrix where col1:chain,col2:resno,col3:resid
#we append dummy residue to patmat by rbind so that all functions work well, a trick
	
	patmat=rbind(patmat,"-"); #anyways - - - wont be detected as ionpair so no worry
	dmat=0;
	intprofile=list(ip=dmat,ipnet=dmat,ap=dmat,apnet=dmat,as=dmat,asnet=dmat,hb=dmat,dis=dmat,catpi=dmat,catpinet=dmat,hp=dmat); # IP	IPnet	AP	APnet	AS	ASnet	HB	Dis	Catpi	Hphob
	
	if(intList$status.ionic & length(grep("[DERKH]",patmat[,3],ignore.case = TRUE))>0) #3rd column contains resid information
	{
		ionmat=checkIonic(patmat,intList$ionpair);
		if(nrow(ionmat)!=0)
		{
		ionmat=matrix(unlist(strsplit(unique(apply(ionmat,1,function(x){paste(x,sep="",collapse=":")})),":")),ncol=10,byrow=TRUE) #remove duplicate
		}
		if(nrow(ionmat)>0){ionnetmat=checkIonnet(patmat,intList$ionpair,intList$dionpair,acclow); 
		intprofile$ip=ionmat;if(nrow(ionnetmat)==0){ionnetmat=0}; intprofile$ipnet=ionnetmat;
		}
	}
	
	if(intList$status.aroaro & length(grep("[FWY]",patmat[,3],ignore.case = TRUE))>0)
	{
		aroaromat=checkAroAro(patmat,intList$aro);
		if(nrow(aroaromat)!=0)
		{
		aroaromat=matrix(unlist(strsplit(unique(apply(aroaromat,1,function(x){paste(x,sep="",collapse=":")})),":")),ncol=12,byrow=TRUE) #remove duplicate
		}
		if(nrow(aroaromat)>0){
		aroaronetmat=checkAroAronet(patmat,intList$aro,intList$daropair,acclow);
		intprofile$ap=aroaromat;if(nrow(aroaronetmat)==0){aroaronetmat=0}; intprofile$apnet=aroaronetmat;					
		}
	}
	
	if(intList$status.aros & length(grep("[FWYMC]",patmat[,3],ignore.case = TRUE))>0)
	{
		arosmat=checkAroS(patmat,intList$aros);
		if(nrow(arosmat)!=0)
		{
		arosmat=matrix(unlist(strsplit(unique(apply(arosmat,1,function(x){paste(x,sep="",collapse=":")})),":")),ncol=11,byrow=TRUE) #remove duplicate
		}
		if(nrow(arosmat)>0){
		arosnetmat=checkAroSnet(patmat,intList$aros,intList$darospair,acclow);
		intprofile$as=arosmat;if(nrow(arosnetmat)==0){arosnetmat=0}; intprofile$asnet=arosnetmat;				
		}
	}
	if(intList$status.hb)
	{
		hbmat=checkHB(patmat,intList$hb);
		if(nrow(hbmat)!=0)
		{
		hbmat=matrix(unlist(strsplit(unique(apply(hbmat,1,function(x){paste(x,sep="",collapse=":")})),":")),ncol=14,byrow=TRUE) #remove duplicate
		}
		if(nrow(hbmat)>0){
		intprofile$hb=hbmat;
		}
	}
	if(intList$status.dis & length(grep("C",patmat[,3],ignore.case = TRUE))>0)
	{
		dismat=checkDis(patmat,intList$dis);
		if(nrow(dismat)!=0)
		{
		dismat=matrix(unlist(strsplit(unique(apply(dismat,1,function(x){paste(x,sep="",collapse=":")})),":")),ncol=8,byrow=TRUE) #remove duplicate
		}
		
		if(nrow(dismat)>0){
		intprofile$dis=dismat;
		}
	}
	if(intList$status.catpi & length(grep("[FWYRK]",patmat[,3],ignore.case = TRUE))>0)
	{
		catpimat=checkCatpi(patmat,intList$catpimat);
		if(nrow(catpimat)!=0)
		{
		catpimat=matrix(unlist(strsplit(unique(apply(catpimat,1,function(x){paste(x,sep="",collapse=":")})),":")),ncol=10,byrow=TRUE) #remove duplicate
		}
		if(nrow(catpimat)>0){
		intprofile$catpi=catpimat;
		catpinet=checkcatpinet(patmat,intList$catpimat,intList$dcatpimat,acclow);
		if(nrow(catpinet)==0){catpinet=0}; intprofile$catpinet=catpinet;
		}
	}
	if(intList$status.hp & length(grep("[AVLIMFWPY]",patmat[,3],ignore.case = TRUE)))
	{
		hpmat=checkHphob(patmat,intList$hp);
		if(nrow(hpmat)!=0)
		{
		hpmat=matrix(unlist(strsplit(unique(apply(hpmat,1,function(x){paste(x,sep="",collapse=":")})),":")),ncol=10,byrow=TRUE) #remove duplicate
		}
		if(nrow(hpmat)>0){
		intprofile$hp=hpmat
		}
	}
return(intprofile); # a list
}

write_intprofile=function(intprof,sumpath.int)
{
# This function was written to write the intprofile
# You give the filename, it will write the local interaction profile

	#intprofile is a list of 11 element
	#  	List of 11 # $ ip,$ ipnet...
colhead=c("\nIONIC\n-----------\nChain1\tResno1\tResid1\tChain2\tResno2\tResid2\tResid1_Acc\tResid2_Acc\tResid1_SS\tResid2_SS\n",
"\nIONIC NETWORK\n----------------\n#res	#int	#Buried	#Exposed	Network detail\n----------------------------------------------------------------------------------\n",
"\nAromatic_Aromatic_Interaction\n--------------------------------\nChain1\tResno1\tResid1\tChain2\tResno2\tResid2\tCentDist\tDihed\tResid1_Acc\tResid2_Acc\tResid1_SS\tResid2_SS\n",
"\nAromatic_Aromatic_Network\n------------------------\n#res	#int	#Buried	#Exposed	Network detail\n----------------------------------------------------------------------------------\n",
"\nAromatic_sulphor Interaction\n--------------------------------\nChain1\tResno1\tResid1\tChain2\tResno2\tResid2\tDist\tResid1_Acc\tResid2_Acc\tResid1_SS\tResid2_SS\n",
"\nAromatic_sulphor Interaction Network\n----------------------------------\n#res	#int	#Buried	#Exposed	Network detail\n----------------------------------------------------------------------------------\n",
"\nHydrogen Bonds\n-----------------\nDch\tDresno\tDresid\tDatom\tAch\tAresno\tAresid\tAatom\tHtype\tDAdist\tHAdist\tDHAang\tHAAAang\tDAAAang\n",
"\nDISULFIDE\n------------------\nChain1\tResno1\tChain2\tResno2\tResid1_Acc\tResid2_Acc\tResid1_SS\tResid2_SS\n",
"\nCat_pi interaction\n------------------\nResid1\tResno1\tChain1\tResid2\tResno2\tChain2\tEes\t-\tEvdw\n",
"\nCation-pi interaction Network\n------------------------\n#res	#int	#Buried	#Exposed	Network detail\n----------------------------------------------------------------------------------\n",
"\nHydrophobic interaction\n------------------\nChain1\tResno1\tResid1\tChain2\tResno2\tResid2\tResid1_SS\tResid2_SS\tResid1_Acc\tResid2_Acc\n"
);	
	for(i in 1:length(intprof))
	{
		if(is.matrix(intprof[[i]]))
		{
			cat(colhead[i],file=sumpath.int,append=TRUE);
			write.table(intprof[[i]],sep="\t",eol="\n",file=sumpath.int,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE);
		}
	}
}


write_intprofile_HTML=function(intprof,filepath)
{
# This function was written to write the intprofile
# You give the filename, it will write the local interaction profile

	#intprofile is a list of 11 element
	#  	List of 11 # $ ip,$ ipnet...
	
	#print(intprof);
	
	if(is.matrix(intprof$ip))
	{
		write2Html("<br><b><u>Ionpair</b></u>",filepath,TRUE);
		write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
		write2Html("<tr><th colspan=3 title=\"Acidic residue\"> Acidic </th> <th colspan=3> Basic </th> <th colspan=2> ASA </th> <th colspan=2> SS </th> </tr>",filepath,append=TRUE);
		write2Html.tableHeader(c("Chain","Res.No","Res.ID","Chain","Res.No","Res.ID","Acidic","Basic","Acidic","Basic"),filepath,append=TRUE);
		write2Html.tableBody(intprof$ip,filepath,TRUE,"<td>");
		write2Html("\n</table>\n",filepath,TRUE);
	}
	if(is.matrix(intprof$ipnet))
	{
		write2Html("<br><b><u>Ionic network</b></u>",filepath,TRUE);
		write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed;\">\n",filepath,TRUE);
		write2Html.tableHeader(c("N.Res","N.Int","Buried","Exposed","Network_Details"),filepath,TRUE);
		#write2Html.tableBody(ionnetmat,filepath,TRUE,"<td style=\"word-wrap: break-word\">");
		intprof$ipnet[,5]=apply(matrix(intprof$ipnet[,5],ncol=1),1,function(x){gsub(",","<br>",x)})
		write2Html.tableBody(intprof$ipnet,filepath,TRUE,"<td>");
		write2Html("\n</table>\n",filepath,TRUE);
	}
	if(is.matrix(intprof$ap))
	{
	write2Html("<br><b><u>Aromatic interactions</b></u>",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr><th colspan=3 title=\"First\"> Aromatic residue 1 </th> <th colspan=3> Aromatic residue 2 </th> <th> Centroid Distance</th><th> Dihedral </th> <th colspan=2> ASA </th> <th colspan=2> SS </th> </tr>",filepath,append=TRUE);
	write2Html.tableHeader(c("Chain","Res.No","Res.ID","Chain","Res.No","Res.ID","(Ang)","(Deg)","Res1","Res2","Res1","Res2"),filepath,append=TRUE);
	write2Html.tableBody(intprof$ap,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	}
	if(is.matrix(intprof$apnet))
	{
	write2Html("<br><b><u>Aromatic-aromatic interaction network</b></u>",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed;\">\n",filepath,TRUE);
	write2Html.tableHeader(c("N.Res","N.Int","Buried","Exposed","Network_Details"),filepath,TRUE);
	#write2Html.tableBody(ionnetmat,filepath,TRUE,"<td style=\"word-wrap: break-word\">");
	intprof$apnet[,5]=apply(matrix(intprof$apnet[,5],ncol=1),1,function(x){gsub(",","<br>",x)})
	write2Html.tableBody(intprof$apnet,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	}
	if(is.matrix(intprof$as))
	{
	write2Html("<br><b><u>Aromatic-Sulphur interactions</b></u>",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr><th colspan=3 title=\"First\"> Aromatic residue </th> <th colspan=3> Sulphur residue </th> <th> Distance </th> <th colspan=2> ASA </th> <th colspan=2> SS </th> </tr>",filepath,append=TRUE);
	write2Html.tableHeader(c("Chain","Res.No","Res.ID","Chain","Res.No","Res.ID","(Ang)","Aro","S","Aro","S"),filepath,append=TRUE);
	write2Html.tableBody(intprof$as,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	}
	if(is.matrix(intprof$asnet))
	{
	write2Html("<br><b><u>Aromatic-sulphor interaction network</b></u>",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed;\">\n",filepath,TRUE);
	write2Html.tableHeader(c("N.Res","N.Int","Buried","Exposed","Network_Details"),filepath,TRUE);
	#write2Html.tableBody(ionnetmat,filepath,TRUE,"<td style=\"word-wrap: break-word\">");
	intprof$asnet[,5]=apply(matrix(intprof$asnet[,5],ncol=1),1,function(x){gsub(",","<br>",x)})
	write2Html.tableBody(intprof$asnet,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	}
	if(is.matrix(intprof$hb))
	{
	write2Html("<br><b><u>Hydrogen bonding interaction</b></u>",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr style=\"background-color:rgb(235,235,235)\"><th colspan=4 title=\"First\"> Donor </th> <th colspan=4> Acceptor </th> <th rowspan=2> Type </th><th rowspan=2> Distance_DA </th> <th rowspan=2> Distance_HA </th> <th rowspan=2> DHA_Angle </th> <th rowspan=2> HAAA_Angle </th> <th rowspan=2> DAAA_Angle </th> </tr>",filepath,append=TRUE);
	write2Html.tableHeader(c("Chain","Res.No","Res.ID","Atom","Chain","Res.No","Res.ID","Atom"),filepath,append=TRUE);
	write2Html.tableBody(intprof$hb,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	}
	if(is.matrix(intprof$dis))
	{
	write2Html("<br><b><u>Disulfide bonds</b></u>",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr><th colspan=2 title=\"First\"> Cys1 </th> <th colspan=2> Cys2 </th> <th colspan=2> ASA </th> <th colspan=2> SS </th> </tr>",filepath,append=TRUE);
	write2Html.tableHeader(c("Chain","Res.No","Chain","Res.No","Cys1","Cys2","Cys1","Cys2"),filepath,append=TRUE);
	write2Html.tableBody(intprof$dis,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	}
	if(is.matrix(intprof$catpi))
	{
	write2Html("<br><b><u>Cation-pi interactions</b></u>",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr ><th colspan=3 title=\"First\"> Cationic residue </th> <th colspan=3> Aromatic residue </th> <th colspan=2> ASA </th> <th colspan=2> SS </th> </tr>",filepath,append=TRUE);
	write2Html.tableHeader(c("Chain","Res.No","Res.ID","Chain","Res.No","Res.ID","Cation","Aromatic","Cation","Aromatic"),filepath,append=TRUE);
	write2Html.tableBody(intprof$catpi,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	}
	if(is.matrix(intprof$catpinet))
	{
	write2Html("<br><b><u>Cation-pi interaction network</b></u>",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed;\">\n",filepath,TRUE);
	write2Html.tableHeader(c("N.Res","N.Int","Buried","Exposed","Network_Details"),filepath,TRUE);
	#write2Html.tableBody(ionnetmat,filepath,TRUE,"<td style=\"word-wrap: break-word\">");
	intprof$catpinet[,5]=apply(matrix(intprof$catpinet[,5],ncol=1),1,function(x){gsub(",","<br>",x)})
	write2Html.tableBody(intprof$catpinet,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	}
	if(is.matrix(intprof$hp))
	{
	write2Html("<br><b><u>Hydrophobic interactions</b></u>",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr ><th colspan=3 title=\"First\"> Residue 1 </th> <th colspan=3> Residue 2 </th> <th colspan=2> SS </th> <th colspan=2> ASA </th> </tr>",filepath,append=TRUE);
	write2Html.tableHeader(c("Chain","Res.No","Res.ID","Chain","Res.No","Res.ID","Res1","Res2","Res1","Res2"),filepath,append=TRUE);
	write2Html.tableBody(intprof$hp,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	}
}

