########## FOLDX #########

foldXcreateMutList=function(foldxMutfile,s)
{
	cat(s,file=foldxMutfile,sep="\n");
}

foldXcreateRunFile=function(pdbname,foldxMutfile,foldxRunfile,prefix,Tmp=298,ph=7,ion=0.05)
{
	cat("<TITLE>FOLDX_runscript;",file=foldxRunfile,sep="",eol="\n");
	cat("<JOBSTART>#;",file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat(paste("<PDBS>",paste(pdbname,";",sep=""),sep=""),file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat("<BATCH>#;\n<COMMANDS>FOLDX_commandfile;",file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat(paste(c("<BuildModel>",prefix,",",foldxMutfile,";"),sep="",collapse=""),file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat("<END>#;\n<OPTIONS>FOLDX_optionfile;",file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat(paste(c("<Temperature>",Tmp,";"),sep="",collapse=""),file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat("<R>#;",file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat(paste(c("<pH>",ph,";"),sep="",collapse=""),file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat(paste(c("<IonStrength>",ion,";"),sep="",collapse=""),file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat("<OutPDB>true;",file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat("<water>-IGNORE;\n<metal>-CRYSTAL;\n<VdWDesign>2;\n<OutPDB>false;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;",file=foldxRunfile,sep="",eol="\n",append=TRUE);
}

foldXcreateRunFileStability=function(pdbname,foldxRunfile,prefix,Tmp=298,ph=7,ion=0.05)
{
	cat("<TITLE>FOLDX_runscript;",file=foldxRunfile,sep="",eol="\n");
	cat("<JOBSTART>#;",file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat(paste("<PDBS>",pdbname,";",sep="",collapse=""),file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat("<BATCH>#;\n<COMMANDS>FOLDX_commandfile;",file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat(paste(c("<Stability>",prefix,";"),sep="",collapse=""),file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat("<END>#;\n<OPTIONS>FOLDX_optionfile;",file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat(paste(c("<Temperature>",Tmp,";"),sep="",collapse=""),file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat("<R>#;",file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat(paste(c("<pH>",ph,";"),sep="",collapse=""),file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat(paste(c("<IonStrength>",ion,";"),sep="",collapse=""),file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat("<OutPDB>true;",file=foldxRunfile,sep="",eol="\n",append=TRUE);
	cat("<water>-IGNORE;\n<metal>-CRYSTAL;\n<VdWDesign>2;\n<OutPDB>false;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;",file=foldxRunfile,sep="",eol="\n",append=TRUE);
}


# exepath n rotabase is imp from user.
#older version prior to foldx drastically modified their version

runFoldx_oldversion=function(mut,exepath="",foldxprefix)
{
foldXcreateMutList(foldxMutfile,mut);
result=system(paste(exepath, " -runfile ", foldxRunfile,sep="",collapse=""),ignore.stdout = TRUE,ignore.stderr = FALSE,intern=TRUE)
return(parseFoldx(foldxprefix));
}

# new code after folx updated their code like command line

runFoldx=function(mut,exepath="",foldxprefix,pdbname)
{
  foldXcreateMutList(foldxMutfile,mut);
  result=system(paste(exepath, " -c BuildModel --mutant-file  ", foldxMutfile," --pdb ",pdbname,sep="",collapse=""),ignore.stdout = TRUE,ignore.stderr = FALSE,intern=TRUE)
  return(parseFoldx("ProtEng.fxout"));
}


# Older version 
runFoldxStability_oldversion=function(pdb,exepath="",foldxprefix)
{
result=system(paste(exepath, " -runfile ", foldxRunfile,sep="",collapse=""),ignore.stdout = TRUE,ignore.stderr = FALSE,intern=TRUE)
return(parseFoldxStability(foldxprefix));
}

# new version after folx changed the code
runFoldxStability=function(pdb,exepath="",foldxprefix,pdbname)
{
  result=system(paste(exepath, " -c Stability  --pdb ", pdbname,sep="",collapse=""),ignore.stdout = TRUE,ignore.stderr = FALSE,intern=TRUE)
  temppdbname=sub(".pdb","",pdbname)
  return(parseFoldxStability(paste(temppdbname,"_0_ST.fxout",sep="",collapse="")));
}

# New version after foldx modified drastically

parseFoldxStability=function(foldxprefix)
{
  print(foldxprefix);
  if(file.exists(foldxprefix))
  {
    foldxdata=readLines(foldxprefix);
    print(foldxdata);
    stabdata=unlist(strsplit(foldxdata[1],"\\t"));
    
    if(is.na(stabdata[1]))
    {
      return("SORRY");
    }
    return(stabdata[2:length(stabdata)]);
  }else{
    return("SORRY");
  }
}



parseFoldxStability_old=function(foldxprefix)
{
  print(foldxprefix);
	if(file.exists(foldxprefix))
	{
	foldxdata=readLines(foldxprefix);
	print(foldxdata);
	stabdata=unlist(strsplit(foldxdata[(which(substr(foldxdata,2,6)=="total"))+1],"\\t"));
	
	if(is.na(stabdata[1]))
	{
	return("SORRY");
	}
	
	return(stabdata[2:length(stabdata)]);
	}else{
	return("SORRY");
	}
}


parseFoldx=function(prefix)
{
	diffile=paste("Dif_",prefix,sep="");
	print(getwd());
	print(list.files(getwd()));
	print(file.exists(diffile));
	# Need to check whether all output begions with line 9 so that 8 lines to be skipped.
	if(!file.exists(diffile))
	{
	print("No diff file");
	return(c("-","-")); # If foldx error then no diff file will be produced
	} 
	foldxdata=read.table(diffile,skip=9,sep="\t",stringsAsFactor=FALSE);
	unlink(paste(c("*",prefix,"*"),sep="",collapse=""));
	return(foldxdata[,2]); #2nd row, 2nd column is fold change in stability
}

protEng_runFoldx=function(protengpref,mut,foldxpath)
{
	stab=round(as.numeric(runFoldx(mut,foldxpath)),2)
	mutmat=cbind(mutmat,FX=stab);
	stabLab=apply(matrix(as.numeric(stab),ncol=1),1,function(x){if(x>0){return("I")}else{return("D")}; })
	mutmat=cbind(mutmat,FX_type=stabLab);
	return(mutmat);
}

print_HTMLgibbs_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Gibbs energy of folding using FOLDX force field</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("Filename",
	"Predicted overall stability of your protein",
	"Contribution of backbone Hbonds","Contribution of sidechain-sidechain and sidechain-backbone Hbonds","Contribution of the VanderWaals",
	"Electrostatic interactions","Penalization for burying polar groups","Contribution of hydrophobic groups",
	"Energy penalization due to VanderWaals’ clashes (interresidue)","Entropy cost of fixing the side chain","Entropy cost of fixing the main chain",
	"","","Cost of having a cis peptide bond","VanderWaals torsional clashes (intraresidue)","Backbone-backbone VanderWaals.These are not considered in the total energy",
	"Electrostatic contribution of the helix dipole","Contribution of water bridges","Contributiuon of disulphide bonds","Electrostatic interaction between molecules in the precomplex",
	"Interactions with bound metals","","",""
	);
	write2Html.tableHeader(c("Filename","TE","BB_HB","SC_HB","VDW","Elec","Sol_P","Sol_H","VDW_C","Ent_SC","Ent_MC","Ent_SL","Ent_ML","Cis","Tor_C","BB_C","HD","WB","DB","Elec_kon","PCB","IE","Ent_Cx","N.Res"),filepath,append=TRUE,tipvec);
}


