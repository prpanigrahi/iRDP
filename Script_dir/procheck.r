get_procheck=function(pdb,prefix="",prochResolution=2.0)
{
  print("Running get Procheck")
	prochtemppdb=paste(prefix,".pdb",sep="",collapse="");
	prochout=paste(prefix,".out",sep="",collapse="");
	write.pdb(pdb,file=prochtemppdb);
	tryCatch({error=system(paste(procheckpath,prochtemppdb,prochResolution),ignore.stderr = FALSE,ignore.stdout = TRUE,intern=TRUE)},error=function(e){print("+++");return("SORRY");});
  #	if(is.list(attributes(error))){return("SORRY");}
	files=list.files(".",pattern=prefix)
print("*****")
	unlink(c(files[which(files!=prochout)],"anglen.log","bplot.log","clean.log","fort.27","nb.log","pplot.log","procheck.prm","secstr.log","tplot.log"));
raw.lines=""
#print("#####")
if(file.exists(prochout)){raw.lines <- readLines(prochout);}else{return("SORRY");}
#	print("$$$$$")
	unlink(prochout);
	raw.lines=raw.lines[1:which(substring(raw.lines,1,33)=="                          M A I N")]
	raw.ind=which(regexpr("^\\s*\\d+.\\s+[A-Za-z]{3}\\s*",raw.lines,perl=TRUE)!=-1)
	serial=substring(raw.lines[raw.ind],1,4);
	chain=substring(raw.lines[raw.ind],5,5);
	resid=aa321(substring(raw.lines[raw.ind],7,9));
	resno=sub(" +","",substring(raw.lines[raw.ind],10,14)); #insertion code too is considered here
	resno=gsub(" ","",resno)	
	ss=substring(raw.lines[raw.ind],17,17);
	ss[ss==" "]="-"
	rama=sub(" +","",substring(raw.lines[raw.ind],23,24));
	return(cbind(serial,chain,resno,resid,ss,rama));
}

procheck=function(pdb,DSSP,prefix="",prochResolution=2.0)
{
prochtemppdb=paste(prefix,".pdb",sep="",collapse="");
prochout=paste(prefix,".out",sep="",collapse="");
	
	write.pdb(pdb,file=prochtemppdb);
	pdb$atom[which(is.na(pdb$atom[,"insert"])),"insert"]="";
	tryCatch({error=system(paste(procheckpath,prochtemppdb,prochResolution),ignore.stderr = FALSE,ignore.stdout = TRUE,intern=TRUE)},error=function(){return("SORRY");})
	if(is.list(attributes(error))){return("SORRY");}
	files=list.files(".",pattern=prefix)
	unlink(c(files[which(files!=prochout)],"anglen.log","bplot.log","clean.log","fort.27","nb.log","pplot.log","procheck.prm","secstr.log","tplot.log"));
	raw.lines=""	
	
if(file.exists(prochout)){raw.lines <- readLines(prochout);}else{return("SORRY");}
	unlink(prochout);

	raw.lines=raw.lines[1:which(substring(raw.lines,1,33)=="                          M A I N")]
	raw.ind=which(regexpr("^\\s*\\d+.\\s+[A-Za-z]{3}\\s*",raw.lines,perl=TRUE)!=-1)
	
	serial=substring(raw.lines[raw.ind],1,4);
	chain=substring(raw.lines[raw.ind],5,5);
	resid=aa321(substring(raw.lines[raw.ind],7,9));
	resno=sub(" +","",substring(raw.lines[raw.ind],10,14)); #insertion code too is considered here
	resno=gsub(" ","",resno)	
	ss=substring(raw.lines[raw.ind],17,17);
	ss[ss==" "]="-"
	rama=sub(" +","",substring(raw.lines[raw.ind],23,24));
	left.ind=which(rama=="l" | rama=="L" | rama=="~l")
	# If no left residues then left.ind=integer(0)

	if(length(left.ind)>0)
	{
		if(is.list(DSSP))
		{
		DSSPres=paste(DSSP$cha,DSSP$res,DSSP$aa,sep="")
		left.phipsi=t(apply(cbind(chain[left.ind],resno[left.ind],resid[left.ind]),1,function(x){DSSP.ind=which(DSSPres==paste(x,sep="",collapse=""));if(length(DSSP.ind)>0){return(c(DSSP$phi[DSSP.ind],DSSP$psi[DSSP.ind]))}else{return(c("-","-"));} }))
		colnames(left.phipsi)=c("phi","psi");
		}else{
		# For some pdb if DSSP fails then you cant estimate the phi,psi
		left.phipsi=matrix(rep("-",times=2*length(left.ind)),ncol=2)
		colnames(left.phipsi)=c("phi","psi");
		}
  
########
left.dist=apply(cbind(chain[left.ind],resno[left.ind]),1,function(x)
{
print(x);
tempreg=regexpr("[-0-9]+",x[2]); #get the resno
tempreg.resno=substr(x[2],tempreg[1],attributes(tempreg)$match.length);
tempreg=regexpr("[A-Za-z]",x[2],perl=TRUE);
tempreg.ins="";
if(tempreg[1]!=-1){tempreg.ins=substr(x[2],tempreg[1],tempreg[1])}
CBO=atom.select(pdb,verbose=FALSE,chain=x[1],resno=as.numeric(tempreg.resno),elety=c("CB","O"))$atom;
tempvec=paste(pdb$atom[CBO,"resno"],pdb$atom[CBO,"insert"],sep="")

tosearch=x[2];
CBO=CBO[which(tempvec==x[2])];

if(length(CBO)==2)
{
return(dist(matrix(as.numeric(pdb$atom[CBO,c("x","y","z")]),ncol=3)));
}else{return("-");}
});

left.dist[which(left.dist!="-")]=round(as.numeric(left.dist[which(left.dist!="-")]),2); # those values which are not "-" round it
########
	
		left=cbind(chain=chain[left.ind],resno=resno[left.ind],resid=resid[left.ind],ss=ss[left.ind],rama=rama[left.ind],left.phipsi,dist=left.dist);
		return(left);
	}else{return("SORRY");}
}

procheck_summary=function(left)
{
	resfreq=tapply(left[,3],factor(left[,3]),length);
	return(c(Total=nrow(left),Total_L=sum(left[,5]=="L"),Total_l=sum(left[,5]=="l"),Total_la=sum(left[,5]=="~l"),ResFreq=paste(names(resfreq),resfreq,sep=":",collapse=",")));
}

#if(IMUT)
#{
#	left=cbind(left,t(apply(left,1,function(x){return(iMut(pdb,x[1],x[2],"G"));})));
#}

print_HTMLprocheck=function(left,filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Details of conformationally strained residues</u><br><br></b>\n",filepath,TRUE);
#	write2Html("Chain1,Resno1,Resid1: Corresponds to Acidic residue details<br>\n",filepath,TRUE);
#	write2Html("Chain2,Resno2,Resid2: Corresponds to Basic residue details<br>\n",filepath,TRUE);
#	write2Html("Resid1_Acc and Resid2_Acc corresponds to relative ASA of acidic and basic residue respectively<br>\n",filepath,TRUE);
#	write2Html("Resid1_SS and Resid2_SS corresponds to secondary structure of the acidic and basic residues respectively<br><br><br>\n",filepath,TRUE);

	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("Chain","Residue number","Residue","Secondary structure [DSSP notation]","Ramachandran plot region [Procheck notation]",
	"phi angle","psi angle","Distance between the beta-carbon and the carbonyl oxygen within the residue"
	);
	write2Html.tableHeader(c("Chain","Res.No","Res.ID","SS","R.Plot_region","Phi","Psi","Distance"),filepath,append=TRUE,tipvec);
	write2Html.tableBody(left,filepath,TRUE,"<td>");
	write2Html("\n</table>\n</html>",filepath,TRUE);
}


print_HTMLleft_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Conformationally strained (CS) residue profile</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("PDB filename","Total number of CS residues",
	paste("Number of CS residues occupying the ",c("L","l","~l")," region of Ramachandran plot [Procheck notation]",sep=""),
	"Format: Residue:count"
	);
	write2Html.tableHeader(c("Filename","Total","L_region","l_region","la_region","Frequency"),filepath,append=TRUE,tipvec);
}

