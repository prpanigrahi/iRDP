HBclassify=function(hb="")
{
	if(ncol(hb)!=14){stop("No of column in hb is not 14, check parseHB()");}
	# We need to convert this into matrix becoz if one MS entry found, it treats it as vector and not matrix
	ms=matrix(hb[which(hb[,"Htype"]=="MS"),],ncol=14);
	sm=matrix(hb[which(hb[,"Htype"]=="SM"),],ncol=14);
	ss=matrix(hb[which(hb[,"Htype"]=="SS"),],ncol=14);
	# matrix() leads to loss of colnames. Hence need to rename
	colnames(ms)=colnames(hb);colnames(sm)=colnames(hb);colnames(ss)=colnames(hb);
	nn=matrix("",ncol=14);cn=matrix("",ncol=14);
	# If no MS types of Hbond exists then no need of searching there
	if(dim(ms)[1]>0)
	{
	# Donor: Mainchain atom of any residue and Acceptor: Sidechain atom of Neutal residues
	nn=matrix(ms[grep("ASP|GLU|ARG|LYS|HIS",ms[,"Aresid"],invert=TRUE),],ncol=14)
	# Donor: Mainchain atom of any residue and Acceptor: Sidechain atom of Charged residues
	cn=matrix(ms[grep("ASP|GLU|ARG|LYS|HIS",ms[,"Aresid"]),],ncol=14)
	}
	# If no SM types of Hbond exists then no need of searching there
	if(dim(sm)[1]>0)
	{
	# Donor: Sidechain atom of neutral residue and Acceptor: Mainchain atom of any residues
	nn=rbind(nn,matrix(sm[grep("ASP|GLU|ARG|LYS|HIS",sm[,"Dresid"],invert=TRUE),],ncol=14));
	# Donor: Sidechain atom of Charged residue and Acceptor: Mainchain atom of any residues
	cn=rbind(cn,matrix(sm[grep("ASP|GLU|ARG|LYS|HIS",sm[,"Dresid"]),],ncol=14));
	}
	if(dim(ss)[1]>0)
	{
	comb=paste(c("ARG","LYS","HIS","ASP","GLU"),rep(c("ASP","GLU","HIS","ASP","GLU"),each=5),collapse="|")
	# SS contains CC/CN=NC/NN combination of residues.
	# First remove CC combinations by filtering SS having not comb, so we get ss_cn_nn
	# ss_cn_nn contains entry where donor and acceptor can be CN or NC or NN
	ss_cn_nn=matrix(ss[grep(comb,paste(ss[,"Dresid"],ss[,"Aresid"]),invert=TRUE),],ncol=14)
	colnames(ss_cn_nn)=colnames(hb);
	
	# As CC removed, by filtering for charged , we Identify CN or NC (both are same), remaining will be NN
	cn=rbind(cn,matrix(ss_cn_nn[c(grep("ASP|GLU|ARG|LYS|HIS",ss_cn_nn[,"Dresid"]),grep("ASP|GLU|ARG|LYS|HIS",ss_cn_nn[,"Aresid"])),],ncol=14));
	nn=rbind(nn,matrix(ss_cn_nn[-c(grep("ASP|GLU|ARG|LYS|HIS",ss_cn_nn[,"Dresid"]),grep("ASP|GLU|ARG|LYS|HIS",ss_cn_nn[,"Aresid"])),],ncol=14));
	}
	nn=matrix(nn[-1,],ncol=14); cn=matrix(cn[-1,],ncol=14);
	
	colnames(cn)=colnames(hb);
	colnames(nn)=colnames(hb);
	
	#pcn=round((length(cn[,1])*100)/length(hb[,1]),3) #percentage of cnhb
	#pnn=round((length(nn[,1])*100)/length(hb[,1]),3)
	
	pcn=length(cn[,1]) #percentage of cnhb
	pnn=length(nn[,1])
	
	return(list(cnhbond=cn,nnhbond=nn,pcnhb=pcn,pnnhb=pnn));
}
HBSummary=function(hb="",nncnhb="")
{
	hbsum=c(total=length(hb[,1]));
	
	#hbsum=c(hbsum,pIntra=round((sum(hb[,"Dch"]==hb[,"Ach"])*100)/length(hb[,1]),3));
	#hbsum=c(hbsum,pInter=round((sum(hb[,"Dch"]!=hb[,"Ach"])*100)/length(hb[,1]),3));
	
	hbsum=c(hbsum,pIntra=sum(hb[,"Dch"]==hb[,"Ach"]));
	hbsum=c(hbsum,pInter=sum(hb[,"Dch"]!=hb[,"Ach"]));
	
	len=length(hb[,"Htype"])
	
	#hbsum=c(hbsum,pMM=round((sum(hb[,"Htype"]=="MM")*100)/len,2));
	#hbsum=c(hbsum,pMS=round((sum(hb[,"Htype"]=="MS")*100)/len,2));
	#hbsum=c(hbsum,pSM=round((sum(hb[,"Htype"]=="SM")*100)/len,2));
	#hbsum=c(hbsum,pSS=round((sum(hb[,"Htype"]=="SS")*100)/len,2));

	hbsum=c(hbsum,pMM=sum(hb[,"Htype"]=="MM"));
	hbsum=c(hbsum,pMS=sum(hb[,"Htype"]=="MS"));
	hbsum=c(hbsum,pSM=sum(hb[,"Htype"]=="SM"));
	hbsum=c(hbsum,pSS=sum(hb[,"Htype"]=="SS"));
	
	#hbtype=round((tapply(hb[,"Htype"],factor(hb[,"Htype"]),length)*100/length(hb[,1])),3);
	#hbsum=c(hbsum,paste(names(hbtype),hbtype,sep=":",collapse=", "));
	hbsum=c(hbsum,pCNHB=nncnhb$pcnhb)
	hbsum=c(hbsum,pNNHB=nncnhb$pnnhb)
	return(hbsum);
}
parseHB=function(infile)
{
print("parsehb");

HB=readLines(infile);
HB <- HB[-(1:(which(substring(HB,1 ,6) == "n    s")))]
m="";
#print("parsehb1");
if(length(HB)>0)
{
	#search for donor
	Dch=substring(HB,1 ,1)
	Dresno=substring(HB,2 ,6)
	Dresno=gsub("^0+","",Dresno) # remove starting zeros
	Dresno=gsub("-$","",Dresno) # insertion code, its - when no insertion code, remove
	Dresid=substring(HB,7,9)
	Datom=substring(HB,11,13)
	Datom=sub(' +$', '', Datom)

	#search for acceptor
	Ach=substring(HB,15 ,15)
	Aresno=substring(HB,16 ,20)
	Aresno=gsub("^0+","",Aresno) # remove starting zeros
	Aresno=gsub("-$","",Aresno) # insertion code, its - when no insertion code, remove
	Aresid=substring(HB,21,23)
	Aatom=substring(HB,25,27)
	Aatom=sub(' +$', '', Aatom)

	Htype=substring(HB,34,35)

	DAdist=as.numeric(substring(HB,29,32))
	HAdist=as.numeric(substring(HB,53,57))

	DHAang=as.numeric(substring(HB,46,51))
	HAAAang=as.numeric(substring(HB,58,63))
	DAAAang=as.numeric(substring(HB,64,69))

	m=matrix(c(Dch,Dresno,Dresid,Datom,Ach,Aresno,Aresid,Aatom,Htype,DAdist,HAdist,DHAang,HAAAang,DAAAang),nrow=length(HB),ncol=14)
	# HOW TO REMOVE OTHER HETATM RESIDUES INTERACTION ?grep("ALA|CYS|GLU|PHE|ALA",hb[,"Dresid"],perl=TRUE)
	index=which(m[,3]=="HOH" | m[,7]=="HOH")
	if(length(index)>0){m=m[-index,]} ## To remove any interaction with water
	colnames(m)=c("Dch","Dresno","Dresid","Datom","Ach","Aresno","Aresid","Aatom","Htype","DAdist","HAdist","DHAang","HAAAang","DAAAang");
	rownames(m)=1:length(m[,1])
}
return(m);
}
print_HBSummary=function(hbsum,outhb)
{
	hbsumnames=c("Total no. of. hbond: ","Percentage of intrachain hbond: ","Percentage of interchain hbond: ",
	"Percentage of hbond type(Donor-Acceptor as MM,MS,SM,SS where M=Mainchain,S=Sidechain): ","Percentage of Charged-Neutral hydrogen bond: ",
	"Percentage of neutral-neutral hydrogen bond: ");
cat("\n\nSummary\n----------\n",file=outhb,sep="",eol="",append=TRUE);
newhbsum=hbsum[1:3];
newhbsum=c(newhbsum,paste(c("MM","MS","SM","SS"),c(hbsum[4:7]),sep=":",collapse=", "));
newhbsum=c(newhbsum,c(hbsum[8:9]));
write.table(cbind(hbsumnames,newhbsum),file=outhb,sep=" ",eol="\n",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE);
}

runHB=function(pdb,exepath="",prefix="")
{
	print("Running hbplus");
	#print(exepath);
	hbtemppdb=paste(prefix,".pdb",sep="",collapse="");
	hbtemphb2=paste(prefix,".hb2",sep="",collapse="");
	write.pdb(pdb,file=hbtemppdb);
	#print(paste(exepath,hbtemppdb,sep=" "));
	#print(getwd());
	
	tryCatch({system(paste(exepath,hbtemppdb, sep = " "),ignore.stderr = FALSE,ignore.stdout=TRUE)},error=function(e){return("SORRY")})
	
	#system("/data2/hbplus/exe/hbplus.exe /var/www/workdir/100/100_hbplus1.pdb",ignore.stderr = FALSE,ignore.stdout=FALSE)
	#system("ls");
	
	#unlink(hbtemppdb);
	if(!file.exists(hbtemphb2)){return("SORRY")}
	HB=parseHB(hbtemphb2);
	#unlink(hbtemphb2);
	return(HB);
}



getHB=function(chres="",hb="")
{
	chain=chres[1];
	resno=chres[2];
	if(is.na(chain)){chain="-"} # HBplus takes - as nochain info
	query=paste(chain,resno,sep="",collapse="");
	hb_res1=paste(hb[,1],hb[,2],sep="");
	hb_res2=paste(hb[,5],hb[,6],sep="");
	
	hb.ind=which(hb_res1==query | hb_res2==query);
	if(length(hb.ind)==0){return("SORRY");} # If not found return SORRY else an matrix
	hbmat=hb[hb.ind,];
	hbmat=matrix(hbmat,ncol=14);
	colnames(hbmat)=colnames(hb);
	return(hbmat)
}

callHbond=function(pdb)
{
			# This is the workflow that is used in Pattern_analysis code
			# Hbond
			hb=runHB(pdb,exepath=hbpluspath,prefix=hbplusprefix);
			status.hb=0;
			if(is.matrix(hb))
			{
			status.hb=1;
			}
return(list(hb=hb,status.hb=status.hb));				
}

checkHB=function(patmat,hb)
{
	# This is the workflow that is used in Pattern_analysis code
	# apply(patmat,1,getHB,hb) returns a list
	# lapply(list,matrix,args) applys to a list object and we apply matrix func so we convert the list to matrix
	# To work well patmat must have a empty row so that apply() will always return a list
#case1: both 176 and 177 form hbond. apply function retunrs a list and it works fine
#> patmat:
#[1,] "A"  "176" "D" 
#[2,] "A"  "177" "R" 
#[3,] "-"  "-"   "-" 

# case2: patmat1=patmat[c(1,3),] where one residue forms hbond and other is dummy residue
# case3: where 10000 doesnot form hbond and we have dummy we get empty hbmat, nrow(hbmat)=0;
# [1,] "A"  "1000000" "D" 
# [2,] "-"  "-"       "-"

	hbmat=do.call(rbind,lapply(apply(patmat,1,getHB,hb),matrix,ncol=14,byrow=FALSE));
	sorry.ind=which(apply(hbmat,1,function(x){sum(x=="SORRY")>0}));
	if(length(sorry.ind)>0){hbmat=hbmat[-sorry.ind,];}
	hbmat=matrix(hbmat,ncol=14);
	return(hbmat);					
}

print_HTMLhb=function(hb,filepath,append=FALSE,html=TRUE,tag="")
{
	if(html){write2Html("<html>\n",filepath,append);}else{write2Html("",filepath,append);}
	write2Html(tag,filepath,TRUE);
#	write2Html("Chain1,Resno1,Resid1: Corresponds to Acidic residue details<br>\n",filepath,TRUE);
#	write2Html("Chain2,Resno2,Resid2: Corresponds to Basic residue details<br>\n",filepath,TRUE);
#	write2Html("Resid1_Acc and Resid2_Acc corresponds to relative ASA of acidic and basic residue respectively<br>\n",filepath,TRUE);
#	write2Html("Resid1_SS and Resid2_SS corresponds to secondary structure of the acidic and basic residues respectively<br><br><br>\n",filepath,TRUE);

	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	write2Html("<tr style=\"background-color:rgb(235,235,235)\"><th colspan=4 title=\"Hydrogen bond donor residue\"> Donor </th> <th colspan=4 title=\"Hydrogen bond acceptor residue\"> Acceptor </th> <th rowspan=2 title=\"Type of hydrogen bond [Format:Donor atomtype acceptor atomtype] (M:Mainchain atom,S:Sidechain atom)\"> Bond_type </th><th rowspan=2 title=\"Distance between Donor(D),Acceptor(A) atom\" > Distance_DA </th> <th rowspan=2 title=\"Distance between Hydrogen(H),Acceptor(A) atom\"> Distance_HA </th> <th rowspan=2 title=\"Angle formed by Donor(D),Hydrogen(H),Acceptor(A) atom\"> DHA_Angle </th> <th rowspan=2 title=\"Angle formed by Hydrogen(H),Acceptor(A),Acceptor Antecednets(AA) atom\"> HAAA_Angle </th> <th rowspan=2 title=\"Angle formed by Donor(D),Acceptor(A),Acceptor Antecednets(AA) atom\"> DAAA_Angle </th> </tr>",filepath,append=TRUE);
	tipvec=c("Chain","Residue number","Three letter code of residue","Atom","Chain","Residue number","Three letter code of residue","Atom");
	write2Html.tableHeader(c("Chain","Res.No","Res.ID","Atom","Chain","Res.No","Res.ID","Atom"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(hb,filepath,TRUE,"<td>");
	write2Html("\n</table>\n",filepath,TRUE);
	if(html){write2Html("</html>\n",filepath,TRUE);}else{write2Html("",filepath,TRUE);}
}

print_HTMLhbond_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Hydrogen bond (HB) summary</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("PDB Filename","Total number of HB","Number of intra-chain HB","Number of inter-chain HB",
	paste("Number of HB of type : ",c("MM","MS","SM","SS")," (Donor-Acceptor atomtypes; M:Main chain, S:Side chain atom)",sep=""),
	"Number of Charged-Neutral HB","Number of Neutral-Neutral HB"
	);

	write2Html.tableHeader(c("Filename","Total","IntraCh","InterCh","MM","MS","SM","SS","CNHB","NNHB"),filepath,append=TRUE,tipvec);
}

