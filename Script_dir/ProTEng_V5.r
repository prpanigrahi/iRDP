# This version is includes PSI-BLAST features

# Module4.r
tryCatch(
{
Sys.setenv(http_proxy=""); # Env variable set
####################### SET All global directories ###########################
script_dir="/var/www/Script_dir/";
work_dir="/var/www/workdir/"; # It has to be shorter coz procheck fails if a longer path comes, donno why.
pdb_dir="/data/ent";
output_dir="/data/out/";
hlinkpath="/out/";
consdir="/data/consurf/";
cslinkpath="/CS/";
upload_dir="/var/www/upload";
############################ Load all Libraries and sources ###########################
loadsource=function(script_dir,file)
{
  source(paste(script_dir,file,sep="",collapse=""));
}

library("bio3d");
library("parallel",quietly=TRUE,verbose=FALSE);
loadsource(script_dir,"bturn.r");
loadsource(script_dir,"PATH.r");
loadsource(script_dir,"foldx.r");
loadsource(script_dir,"automute.r");
loadsource(script_dir,"iMut.r");
loadsource(script_dir,"mupro.r");
loadsource(script_dir,"general.r");
loadsource(script_dir,"writehtml.r");
loadsource(script_dir,"ssbond.r");
loadsource(script_dir,"dssp.r");
loadsource(script_dir,"procheck.r");

######################## Fetch arguments ##############################
args=commandArgs(trailingOnly = TRUE);
optionfile=args[1]; # Arg1 is option file

########################################################### Get prefix ########################################################
options=read.table(optionfile,stringsAsFactors=FALSE) # Read the option file
prefix=options[which(options$V1=="prefix"),2] # First row of option file is the session id to be used as prefix
print(c("Prefix entered",prefix));

########################################################### Check mode of operation ########################################################
modeofop=""; #mode of operation [bypdb, byupload]

if(options[which(options$V1=="pdbmode"),2]!="-")
{
  modeofop="bypdb";
}else{
  modeofop="byupload";
}

########################################################### Obtain filenames as per mode of operation ###################################################
filenames="";
pdbid="";

if(modeofop=="bypdb")
{
  pdbids=options[which(options$V1=="pdbmode"),2];
  pdbid=tolower(pdbids);
  filenames=tolower(unique(paste("pdb",pdbids,".ent",sep="")));
  print(filenames); #print(length(filenames));
}else{
  uploaddir=paste(upload_dir,"/",prefix,sep="",collapse="");
  filenames=list.files(uploaddir,pattern=".ent$|.pdb$"); #accept pdb and ent files only  
  pdb_dir=uploaddir;
}

print(pdb_dir);
print(modeofop);
print(filenames);
#count=70477; #filenames=list.files(pdb_dir,pattern=".ent$");

########################################################### Settings for which programs to run ###################################################
usetool=options[which(options$V1=="tool"),2]; #which tools to use
conserv=as.numeric(options[which(options$V1=="conserv"),2]);
foldph=options[which(options$V1=="foldph"),2]; # Foldx variables
foldtemp=options[which(options$V1=="foldtemp"),2];
foldion=options[which(options$V1=="foldion"),2];
automute_mode=options[which(options$V1=="AM"),2];
imph=options[which(options$V1=="imph"),2]; ## I-mutant variables
imtemp=options[which(options$V1=="imtemp"),2]; ## I-mutant variables
option=options[which(options$V1=="option"),2]; # which rules to use
sstop=options[which(options$V1=="sstop"),2]; # no of ssbond to predict
if(option==1){usetool=0;}
########################################### If mode is custom mutation then parse mutation file ############################################
mutationfile=NA;
if(option==5)
{
      mutationfile=args[2]; # Arg1 is option file
      mutlist="";
      if(!is.na(mutationfile))
      {
      mutlist=toupper(readLines(mutationfile));
      }
#      print(mutlist);

}

########################################## Set the working directory ########################################
work_dir=paste(work_dir,prefix,sep="",collapse="");
createDir(work_dir);
setwd(work_dir);
orig_wd=getwd();

##################### Create Output directory and log, status files ###########################
subpref1="HTML";
subpref2="TEXT";
dirpath=paste(output_dir,prefix,sep="",collapse="");
htmldir=paste(dirpath,"/",subpref1,sep="",collapse="");
txtdir=paste(dirpath,"/",subpref2,sep="",collapse="");
htmloutfile=paste(htmldir,"/result.html",sep="",collapse="");
txtoutfile=paste(txtdir,"/result.txt",sep="",collapse="");
ssbondMut=paste(hlinkpath,prefix,"/TEXT/result.txt",sep="",collapse="");

#Status file
statusfile=paste(htmldir,"/status.txt",sep="",collapse="");
logfile=paste(htmldir,"/log.txt",sep="",collapse="");

createDir(dirpath);
createDir(htmldir);
createDir(txtdir);

################ Read the pdb file ############################
pdbfile=paste(pdb_dir,"/",filenames[1],sep="");

if(!file.exists(pdbfile))
{
  print("checking file in PDB");
  
  pdblines="";
  tryCatch({pdblines=readLines(paste("http://www.rcsb.org/pdb/files/",pdbid,".pdb",sep="",collapse=""))},error = function(e){print("pdb cant download");
cat(paste("\n!!! INVALID PDB !!! ","\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);stop();
  }, finally = print("Hello"))
  
  if(length(grep("ATOM",pdblines))>0)
  {
    print("Downloaded file from PDB");
    cat(pdblines,sep="\n",eol="\n",file=filenames[1]);
    pdbfile=filenames[1];
  }else{
    print(c("PDB file not found",pdbfile));
    cat(paste("\n!!! PDB FILE NOT FOUND !!! ",filenames[1],"\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
    stop(); 
  }
}
cat("Job started....\n\n",file=logfile,sep="",eol="");
# After download or obtain of uploaded pdb file, if no ATOM record then?

if(length(grep("^ATOM",readLines(pdbfile)))==0)
{
  cat(paste("\n*** NO ATOM RECORDS *** ",filenames[1],"\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
  stop(); 
}

pdb=read.pdb(pdbfile,maxlines=100000,verbose=TRUE); # read the wild pdb file filenames[1]
print("+++");
uniqch=uniqueChain(pdb);

if(!sum(is.na(uniqch))>0)
{
	inds=atom.select(pdb,chain=uniqch,verbose=FALSE); 
	if(length(inds$atom)<=1){next;} 
	pdb=trim.pdb(pdb,inds=inds);
}
				

################################### Declare All Temporary variables #########################
protengpref="ProtEng.pdb"; # the prefix for temp pdb files on whioch all mutation to be worked out
protengdssp="ProtEng.dssp"; #used by i-mutant
status.dssp=0;
mutmat=""; # the mutation matrix
totmut=500; #maximum number of mutations to process
qtag="chain"; #query tag

#FoldX variables
foldxpdbout="ProtEng"; # Foldx mutant pdb prefix generated by foldx
foldxMutfile="individual_list_FoldXMutfile";
foldxRunfile="FoldxRunfile";
foldxprefix="foldxtemp";

#Mupro variables
muproRunfile="muprorun.txt";

## SSBOND variables ###
sspdb="sstemp.pdb";
ssout="ssmap.txt";
ssbondmutfile="mutant.pdb"; #DO NOT CHANGE THIS ELSE YOU HAVE TO CHANGE IN SSBOND.R FILE


############## Write Wild pdb file ####################
write.pdb(pdb,file=protengpref); # write the processed pdb file on which the mutations to be carrid out

############## If FoldX then create Runfile and Copy Rotamer file to workdir ####################
if(usetool==0)
{  
  #foldXcreateRunFile(protengpref,foldxMutfile,foldxRunfile,foldxprefix,foldtemp,foldph,foldion);
  system(paste("cp ",foldxrotamer," .")); #copy the foldx rotamer library
}

############## If I-mutant then run dssp, set protengdssp and status.dssp ####################
if(usetool==2)
{
error=system(paste(dssppath, " -i ", protengpref, " -o ", protengdssp, sep = ""),ignore.stderr = FALSE,intern=TRUE)
if(! is.list(attributes(error))){status.dssp=1;}
}

############## If Automute then generate random pdbid, copy wild pdb file to automute directory ####################
pdbpref=paste(paste(sample(letters,4,replace=TRUE),collapse=""),sep="");
pdbfile=paste(automutdir,"pdb",pdbpref,".ent",sep="",collapse="");
system(paste("cp ",protengpref," ",pdbfile)); #copy the foldx rotamer library

########  nodetect function ##################
nodetect=function(txt)
{
write2Html("<html>\n<body>",htmloutfile);		
write2Html(txt,htmloutfile,append=TRUE);
write2Html("</body></html>",htmloutfile,append=TRUE);
}

#################### Write HTML header Start ###################
write2Html("<html>\n",htmloutfile);


##################### Generate mutmat ##########################
if(option==1)
{
      cat("Running SSBOND....\n",file=logfile,sep="",eol="",append=TRUE);
      print("Running SSBOND for disulfide insertion...");
      convertPDB(pdb,sspdb,ssout);
      SSBlist=runSSBOND(sspdb,exepath=ssbondpath);
      #print(SSBlist);
      if(is.matrix(SSBlist))
      {
          SS=convertSSBOND(SSBlist,ssout);
          SS_top=SS[which(SS[,"Conf_no"]==1),]
          SS_print=SS_top[order(as.numeric(SS_top[,"Totnrg"])),]
          SS_print=SS_print[,-c(1,8)];
          
          cat("\nFinished Running SSBOND....\n",file=logfile,sep="",eol="",append=TRUE);
          cat(paste("\n\nTotal number of residue pair that are predicted to be ideal for disulfide bond insertion: ",nrow(SS_print),"\n\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
          
          cat("Respair for Disulfide insertion\n------------------------------\n\nNo",file=txtoutfile,sep="",eol="\t");
          write.table(SS_print,file=txtoutfile,sep="\t",eol="\n",append=TRUE,quote=FALSE)
          
          write2Html(paste(c("<b><a href=\"", ssbondMut, "\" target=\"_blank\" > Download all predicted respairs </a><br><br>\n"),sep="",collapse=""),htmloutfile,TRUE);
          
          write2Html("<b><u>Respair for Disulfide insertion</u><br></b>\n",htmloutfile,TRUE);
          write2Html("Table below gives the possible residue pairs predicted by SSBOND program which are ideal for disulfide bond insertion.\n
<br>Disulfide bonds were introduced and its effect on overall protein stability was predicted using FoldX as I (Increasing Stability), D (Decreasing stability).<br>
Residue conservation scores are also estimated in the scale of 1-9 (9: Highly conserved)<br>For more details <a href=\"http://irdp.ncl.res.in/help_istability.html#results\" target=\"_blank\" >click here</a><br>",htmloutfile,TRUE);
          write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",htmloutfile,TRUE);
          write2Html("<tr><td colspan=16 title=\"Results of SSBOND\" >SSBOND results</td><td colspan=2 title=\"Stability prediction due to disulfide bond insertion\">FoldX Stability Prediction</td> <td colspan=2 title=\"Conservation score of the two mutation sites\" >Conservation scores </td> </tr>\n",htmloutfile,TRUE);
          write2Html.tableHeader(c("Chain","Res1_No","Res1_ID","Chain","Res2_No","Res2_ID","SG_Distance","Chi1_1","Chi2_1","Chi3","Chi1_2","Chi2_2","Chi.Energy","Tau.Energy","Dis.Energy","T.Energy","Score","Stability Change","CScore_1","CScore_2"),htmloutfile,append=TRUE);
          
          ssfoldpred=as.numeric(sstop); #How manu mutations to predict stability by foldx
          if(ssfoldpred>nrow(SS_print)){ssfoldpred=nrow(SS_print)}
          if(ssfoldpred>totmut){ssfoldpred=totmut}
          ssfoldfinal=matrix(SS_print[1:ssfoldpred,],ncol=ncol(SS_print));# only tht top matrix
          ss1=paste(aa321(SS_print[1:ssfoldpred,3]),SS_print[1:ssfoldpred,1],SS_print[1:ssfoldpred,2],"C",sep="");
          ss2=paste(aa321(SS_print[1:ssfoldpred,6]),SS_print[1:ssfoldpred,4],SS_print[1:ssfoldpred,5],"C",sep="");
          ssfold=paste(ss1,ss2,sep=",");
          ssfold=paste(ssfold,";",sep=""); #foldx format
          print(ssfold);
    
                        cat(paste("\n Predicting Stability for only top ",ssfoldpred," respair by FOLDX\n\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
                        
                        for(ssrow in 1:length(ssfold))
                        {
                            #### Consurf
                            wildcons1="-";
                            wildcons2="-";
                            consfile1=paste(consdir,pdbid,"_",ssfoldfinal[ssrow,1],".txt",sep="",collapse=""); # for cys1
                            consfile2=paste(consdir,pdbid,"_",ssfoldfinal[ssrow,4],".txt",sep="",collapse=""); # for cys2
                        
                            if(file.exists(consfile1))
                            {
                              conscmd=paste("cat ",consfile1," | grep ",ssfoldfinal[ssrow,3],ssfoldfinal[ssrow,2],":",ssfoldfinal[ssrow,1],sep="",collapse="");
                              wildcons1=NA;
                              tryCatch({wildcons1=system(conscmd,ignore.stdout=FALSE,ignore.stderr=FALSE,intern=TRUE);},error=function(e){});
                            }
                            if(file.exists(consfile2))
                            {
                              conscmd=paste("cat ",consfile2," | grep ",ssfoldfinal[ssrow,6],ssfoldfinal[ssrow,5],":",ssfoldfinal[ssrow,4],sep="",collapse="");
                              wildcons2=NA;
                              tryCatch({wildcons2=system(conscmd,ignore.stdout=FALSE,ignore.stderr=FALSE,intern=TRUE);},error=function(e){});
                            }
                            if(!is.na(wildcons1)){wildcons1=substr(wildcons1,32,32)}
                            if(!is.na(wildcons2)){wildcons2=substr(wildcons2,32,32)}
                            
                            ### Consurf  
                            cat(paste("\n Predicting Stability for pair: ",ssrow," (",ssfold[ssrow],")\n\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
                            foldxres=runFoldx(ssfold[ssrow],foldxpath,foldxprefix,protengpref);
                            ssfoldres=c("","");
                            #print(foldxres);
                            if(sum(foldxres=="-")!=2)
                            {
                              stab=round(as.numeric(foldxres),2);
                              stablab="";
                              if(stab>0){stablab="D";}else{stablab="I"};
                              ssfoldres=c(stab,stablab);
                            }else{
                              ssfoldres=c("-","-");
                            }
                            file1=paste(foldxpdbout,"_1.pdb",sep="",collapse="");
                            file2=paste("WT_",foldxpdbout,"_1.pdb",sep="",collapse="");
                            if(file.exists(file1)){unlink(file1);   }
                            if(file.exists(file2)){unlink(file2);   }
                            write2Html.tableBody(matrix(c(ssfoldfinal[ssrow,],ssfoldres,wildcons1,wildcons2),ncol=(ncol(ssfoldfinal)+4)),htmloutfile,TRUE,"<td>");
                        }
                        
                        write2Html("</html>\n",htmloutfile,TRUE);
                        #	if(file.exists(ssbondmutfile))
                        #	{
                        #	file.copy(ssbondmutfile,htmldir);
                        #	}
      }else{
        cat("\n No disulfide bond detected\n ",file=logfile,sep="",eol="",append=TRUE);
      }
}else{
	if(option==2)
	{
	# Beta-turn 2nd position proline insertion
		BT=runPROMOTIF(pdb,promotifpath,prefix);
		
		if(is.matrix(BT))
		{
			BT2P.ind=which(substring(BT[,4],2,2)!="P");
			if(length(BT2P.ind)>0)
			{
			BT1=matrix(BT[BT2P.ind,],ncol=13);
			mutmat=cbind(BT1[,1],(as.numeric(BT1[,2])+1),substr(BT1[,4],2,2),"P");
			print(mutmat)
			cat("Chain\tRes.No\tWild_Residue\tMut_Residue\tScore\tStability\tCScore\n",file=txtoutfile,sep="",eol="");
			write2Html("<b><u>Insertion of proline residues at 2nd position of beta-turn</u><br><br></b>\n",htmloutfile,TRUE);
			write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",htmloutfile,TRUE);
			write2Html.tableHeader(c("MutantPDB","Chain","Res.No","Wild_Residue","Mut_Residue","Score","Stability","CScore"),htmloutfile,append=TRUE,tipvec=c("Download the Mutant PDB","Chain","Residue  number","Wild residue","Mutant residue","Stability prediction scores","Stability change [I: Increase,D: Decrease]","Conservation scores"));
			}else{
			nodetect("No beta-turns were detected in this protein");
			}
		}else{
		nodetect("No beta-turns were detected in this protein");
		}
		
	}else{
		if(option==3)
		{
		dsspprefix=paste(prefix,"_dssp",sep="",collapse="");
		DSSP=dssp_new(pdb,exepath=dssppath,prefix=dsspprefix); #Runs DSSP
			if(is.list(DSSP))
			{
				if(length(DSSP$helix$length)>0)
				{
				ncap.ind=apply(matrix(paste(DSSP$helix$chain,DSSP$helix$start,sep=""),ncol=1),1,function(x){which(paste(DSSP$cha,DSSP$res,sep="")==x)}) #DSSP index of all Ncap helix position
				p.ncapind=ncap.ind[which(DSSP$aa[ncap.ind]!="P")]  #DSSP index of only Ncap helix position where N-cap is not proline
					if(length(p.ncapind)>0)
					{ 
					mutmat=matrix(cbind(chain=DSSP$cha[p.ncapind],resno=DSSP$res[p.ncapind],resid=DSSP$aa[p.ncapind],"P"),ncol=4); #a matrix	
					cat("Chain\tRes.No\tWild_Residue\tMut_Residue\tScore\tStability\tCScore\n",file=txtoutfile,sep="",eol="");
					write2Html("<b><u>Insertion of proline residues at N-cap of helices</u><br><br></b>\n",htmloutfile,TRUE);
					write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",htmloutfile,TRUE);
					write2Html.tableHeader(c("MutantPDB","Chain","Res.No","Wild_Residue","Mut_Residue","Score","Stability","CScore"),htmloutfile,append=TRUE,tipvec=c("Download the Mutant PDB","Chain","Residue  number","Wild residue","Mutant residue","Stability prediction scores","Stability change [I: Increase,D: Decrease]","Conservation scores"));
					}else{
					nodetect("All helices N-cap positions are with Proline residues only");
					}
						
				}else{
				nodetect("No helix is found on this protein");
				}
			}else{
			nodetect("DSSP could not run on this protein");
			}
		}else{
			if(option==4)
			{
			# conf strain
			procheckprefix=paste(prefix,"_procheck",sep="",collapse="");
			dsspprefix=paste(prefix,"_dssp",sep="",collapse="");
			DSSP=dssp_new(pdb,exepath=dssppath,prefix=dsspprefix); #Runs DSSP
			left=procheck(pdb,DSSP,procheckprefix); # returns a matrix or SORRY
      if(is.matrix(left))
				{
        if(nrow(left)==1){mutmat=cbind(matrix(left[,c(1,2,3)],ncol=3),"G");}else{mutmat=cbind(left[,c(1,2,3)],"G");}
        cat("Chain\tRes.No\tWild_Residue\tMut_Residue\tScore\tStability\tCScore\n",file=txtoutfile,sep="",eol="");
				write2Html("<b><u>Conformationally strained residues mutated to Glycine</u><br><br></b>\n",htmloutfile,TRUE);
				write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",htmloutfile,TRUE);
				write2Html.tableHeader(c("MutantPDB","Chain","Res.No","Wild_Residue","Mut_Residue","Score","Stability","CScore"),htmloutfile,append=TRUE,tipvec=c("Download the Mutant PDB","Chain","Residue  number","Wild residue","Mutant residue","Stability prediction scores","Stability change [I: Increase,D: Decrease]","Conservation scores"));
				}else{
				nodetect("No conformationally strained residues found");
				}
			}
      else{
			print("Doing custom mutation");
        # Do custom mutation
        print(mutationfile);
				if(!is.na(mutationfile))
				{
				mutlist=gsub(";","",mutlist);
				mutlist=mutlist[grep("^[A-Z][A-Z0-9]\\d+[A-Z]$",mutlist)]; #Only Single mutant will be accepted and Format: [A-Z] [A-Z] at least one digit and [A-Z]
        print(mutlist);
				mutmat=matrix(cbind(substr(mutlist,2,2),substr(mutlist,3,nchar(mutlist)-1),substr(mutlist,1,1),substr(mutlist,nchar(mutlist),nchar(mutlist))),ncol=4);
				print(mutmat);
				cat("Chain\tRes.No\tWild_Residue\tMut_Residue\tScore\tStability\tCScore\n",file=txtoutfile,sep="",eol="");
				write2Html("<b><u>The effect of mutations you have suggested are given below below</u><br><br></b>\n",htmloutfile,TRUE);
				write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",htmloutfile,TRUE);
				write2Html.tableHeader(c("MutantPDB","Chain","Res.No","Wild_Residue","Mut_Residue","Score","Stability","CScore"),htmloutfile,append=TRUE,tipvec=c("Download the Mutant PDB","Chain","Residue  number","Wild residue","Mutant residue","Stability prediction scores","Stability change [I: Increase,D: Decrease]","Conservation scores"));
				}else{
				nodetect("Error: please mail author");
				}
			}
		}

	}
}

############### Conservation analysis #######################


res.foldx=""; #result of foldx
print(mutmat);
print(is.matrix(mutmat))

if(is.matrix(mutmat)) 
{
          ############ Those mutations where insertion codes are there, it should be removed (promotif, procheck detects insertion code).###################
          insert.ind=grep("[A-Za-z]",mutmat[,2]);
          insert.ind=c(insert.ind,which(is.na(mutmat[,2]))); #Sometimes we have seen that resno is NA, in that case exclude that too
          
          
          if(length(insert.ind)>0)
          {
           mutmat=mutmat[-insert.ind,]; #Remove those rows where insertion code is there at resno
           mutmat=matrix(mutmat,ncol=4);
          }
          ################### Final generation of mutmat if it exceeds 50 (totmut) ########
          # Except ssbond for rest we will use this code
          cat("Following are the mutations to be done for the parameter you have selected\n\n",file=logfile,sep="",eol=""); print("*")
          cat(paste(c("Total number mutations:",nrow(mutmat)),sep="",collapse=""),file=logfile,sep="",eol="\n\n",append=TRUE);print("**")
          cat("Chain\tRes.No\tWild_Residue\tMut_Residue\n----------------------------\n",file=logfile,sep="",eol="",append=TRUE);print("***")
          write.table(mutmat,sep="\t",eol="\n",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file=logfile);print("****")
          print(nrow(mutmat));
          if(nrow(mutmat)>totmut)
          {
          	cat(paste("Total number mutations exceeds ",totmut,"\n Only Top ",totmut," mutations are considered\n",sep="",collapse=""),file=logfile,sep="",eol="\n\n",append=TRUE);
          	mutmat=mutmat[1:totmut,];
          }
if(nrow(mutmat)>0)
{
          if(conserv==1)
          {
            cat("##### Now Doing Conservation Analysis ####\n\n",file=logfile,sep="",eol="",append=TRUE);
            mutchains=unique(mutmat[,1]); # The chains in which mutations is to be generated
            if(modeofop=="byupload")
            {
              chmat=uniqueChain_map(pdb); #obtain unique chain matrix where col1 enlists all chain where col2 the unique chain that it corresponds to
              chmat.ind=which(chmat[,1] %in% mutchains);
              mutchains.uniq=unique(chmat[chmat.ind,2]);
              
              for(mutch in mutchains.uniq)
              {
                cat(paste("\nConservation analysis for chain ",mutch,sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
                print(c("Conservation analysis for chain ",mutch));
                mutch.ind=atom.select(pdb,chain=mutch,elety="CA",verbose=FALSE)$atom; #the CA indexes
                mutch.seq=paste(aa321(pdb$atom[mutch.ind,"resid"]),sep="",collapse=""); # The aa sequence of given chain obtained from CA tags
                
                #print(c("Mutchseq",mutch.seq));
                mutch.map=cbind(1:length(mutch.ind),paste(pdb$atom[mutch.ind,"chain"],pdb$atom[mutch.ind,"resno"],pdb$atom[mutch.ind,"resid"],sep="")); #To be written in mapfile col1: serial number of residue, col2:[Chain.Resno.Resid]
                seqfile=paste(qtag,mutch,".txt",sep="",collapse=""); #will be used as query fasta sequence
                seqpref=paste(qtag,mutch,"_cons",sep="",collapse=""); # The prefix [seq_chain.txt]
                cat(paste(">",qtag,mutch,"\n",mutch.seq,"\n",sep="",collapse=""),file=seqfile,sep="",eol=""); #write sequence to seq_A.txt file, will be used for blast
                
                seqcmd=paste("perl ",conspsiblast," ",seqfile," ",seqpref," ",logfile,sep="",collapse="");
                print(seqcmd);
                system(seqcmd,ignore.stdout=FALSE,ignore.stderr=FALSE);
                
                # Suppose seqpref is chainA_cons.Blast runs successfully, we will get chainA_cons_conserv.txt output file created by psiblast.pl program
                # This file has 5 column. Col1: seqno, Col2: sequence, Col3: probability, Col4: color grade, Col5: variability
                # We take seq from col2 of this file and store as pssmseq. The original query seq obtained from CA atom tags and used as query.txt to psiblast is stored in mutch.seq.
                # If mutch.seq and pssmseq are of same length then output is perfect, we can generate mapfile
                # Else no map file is generated but still we have the chainA_cons_conserv.txt file
                # In those case where map files avaibale, then we give this as output and CScore can be given
                # In those case due to sequence ambiguity between mutch.seq and pssmseq, no map files are generated, only chainA_cons_conserv.txt file can be given and CScore will be blank
                
                pssmfile=paste(seqpref,"_conserv.txt",sep="",collapse="");
                
                if(file.exists(pssmfile))
                {
                  pssmmat=read.table(pssmfile,stringsAsFactors=FALSE,skip=6);  
                  pssmseq=paste(pssmmat[,2],sep="",collapse=""); #sequence of pssm file in col2
                  pssm.match=regexpr(pssmseq,mutch.seq);
                  
                  if(pssm.match==1 && nchar(pssmseq)==nchar(mutch.seq))
                  {
                    mapfile=paste(qtag,mutch,"_map.txt",sep="",collapse=""); #map file
                    cat("Col1:Serian number\nCol2:Residue details as [Chain Resno Resid]\nCol3:Residue\nCol4:Weighted observed percentages of occurance of the residue (0-1). 1 means occuring all hits, highly conserved\nCol5: Value of Col4 scaled between 0-9 where 9 is fully conserved residue\nCol6: Other kind of observed residues in the position\n",file=mapfile,sep="",eol="");
                    write.table(cbind(mutch.map,pssmmat[,2:5]),file=mapfile,sep="\t",eol="\n",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE); #write map file
                  }
                }
              }
              
              pssmfilenames= list.files(".","_conserv.txt"); # Original 5column output of psiblast.pl program
              mapfilenames= list.files(".","_map.txt"); # If no sequence ambiguity between mutch.seq and pssmseq, map files are generated
              psiblastfilenames= list.files(".","_cons.html"); # If psiblast runs succcessfully, we get html outputs
              pssmoutfilenames=list.files(".","_cons_pssm.txt"); # if psiblast tuns successfully we get pssm outputs
              
              if(length(pssmfilenames)>0)
              {
                write2Html("<h3>Conservation analysis</h3>\n1. Conservation scores:  ",htmloutfile,TRUE);
                #  write2Html(paste("<a href='",hlinkpath,prefix,"/",subpref1,"/",entfilenames,"' target=\"_blank\" >",entfilenames,"</a>",sep="",collapse="&nbsp;&nbsp;"),htmloutfile,TRUE);
                
                # If map files available then give map files else give Original 5column output of psiblast.pl program
                apply(matrix(mutchains,ncol=1),1,function(x){
                  if(file.exists(paste(qtag,chmat[which(chmat[,1]==x),2],"_map.txt",sep="",collapse="")))
                  {
                    file.copy(from=paste(qtag,chmat[which(chmat[,1]==x),2],"_map.txt",sep="",collapse=""),to=htmldir);  
                    write2Html(paste("<a href='",hlinkpath,prefix,"/",subpref1,"/",qtag,chmat[which(chmat[,1]==x),2],"_map.txt","' target=\"_blank\" > ",qtag,x," conservation","</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",sep="",collapse=""),htmloutfile,TRUE);  
                  }else{
                    file.copy(from=paste(qtag,chmat[which(chmat[,1]==x),2],"_conserv.txt",sep="",collapse=""),to=htmldir);  
                    write2Html(paste("<a href='",hlinkpath,prefix,"/",subpref1,"/",qtag,chmat[which(chmat[,1]==x),2],"_conserv.txt","' target=\"_blank\" > ",qtag,x," conservation","</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",sep="",collapse=""),htmloutfile,TRUE);    
                  }
                });
                
                write2Html("<br>2. Blast Outputs:  ",htmloutfile,TRUE);
                #  write2Html(paste("<a href='",hlinkpath,prefix,"/",subpref1,"/",blastfiltfilenames,"' target=\"_blank\" >",blastfiltfilenames,"</a>",sep="",collapse="&nbsp;&nbsp;"),htmloutfile,TRUE);
                
                apply(matrix(mutchains,ncol=1),1,function(x){
                  file.copy(from=paste(qtag,chmat[which(chmat[,1]==x),2],"_cons.html",sep="",collapse=""),to=htmldir);  
                  write2Html(paste("<a href='",hlinkpath,prefix,"/",subpref1,"/",qtag,chmat[which(chmat[,1]==x),2],"_cons.html","' target=\"_blank\" >",qtag,x," blast hits","</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",sep="",collapse=""),htmloutfile,TRUE);
                });
                write2Html("<br>3. PSSM Matrix: ",htmloutfile,TRUE);
                #  write2Html(paste("<a href='",hlinkpath,prefix,"/",subpref1,"/",msafilenames,"' target=\"_blank\" >",msafilenames,"</a>",sep="",collapse="&nbsp;&nbsp;"),htmloutfile,TRUE);
                
                apply(matrix(mutchains,ncol=1),1,function(x){
                  file.copy(from=paste(qtag,chmat[which(chmat[,1]==x),2],"_cons_pssm.txt",sep="",collapse=""),to=htmldir);  
                  write2Html(paste("<a href='",hlinkpath,prefix,"/",subpref1,"/",qtag,chmat[which(chmat[,1]==x),2],"_cons_pssm.txt","' target=\"_blank\" >",qtag,x," PSSM","</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",sep="",collapse=""),htmloutfile,TRUE);
                });
                
                write2Html("<br><br>\n",htmloutfile,TRUE);
              }
              
                            
            }else{
              
              write2Html("<h3>Consurf conservation analysis</h3>\n",htmloutfile,TRUE);
              mutchains.temp=paste(pdbid,"_",mutchains,".txt",sep="");
              
              apply(matrix(mutchains.temp,ncol=1),1,function(x)
              {
                if(file.exists(paste(consdir,x,sep="",collapse="")))
                {
                  write2Html(paste("<a href='",cslinkpath,x,"' target=\"_blank\" >",x,"</a>&nbsp;&nbsp;&nbsp;&nbsp;",sep="",collapse=""),htmloutfile,TRUE);
                }
              }
              );
              write2Html("<br><br>\n",htmloutfile,TRUE);
            }
            cat("##### Finished conservation analysis\n\n",file=logfile,sep="",eol="",append=TRUE);
          }
         
          ###################### Loop over Mutation to predict stability #################
          for(mutrow in 1:nrow(mutmat))
          {
          	temp.result="";
          	foldxnewpdb="";
          	cat(paste(c("\n###############\n",mutrow,": Processing Mutation:",mutmat[mutrow,3],mutmat[mutrow,1],mutmat[mutrow,2],mutmat[mutrow,4],"\n################\n"),sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
            print(usetool);
          	if(usetool==0)
          	{
          #  tryCatch(
          #  {
          	# FOLDX
            print("Running foldx");
          	mut=paste(mutmat[mutrow,3],mutmat[mutrow,1],mutmat[mutrow,2],mutmat[mutrow,4],";",sep="",collapse="\n");
          	print(c(mutrow,mut));
          	cat("    Running FoldX....\n",file=logfile,sep="",eol="",append=TRUE);
          	foldxres=runFoldx(mut,foldxpath,foldxprefix,protengpref);
            print("999");
            print(foldxres);
            muttempname=paste(mutmat[mutrow,3],mutmat[mutrow,1],mutmat[mutrow,2],mutmat[mutrow,4],sep="",collapse="\n"); #will be used in naming of mutant pdb file generated by foldx
          	
          	#print(foldxres);
          		if(sum(foldxres=="-")!=2)
          		{
          			foldxoldpdb=paste(foldxpdbout,"_1.pdb",sep="",collapse="");
          			foldxnewpdb=paste(foldxpdbout,"_",muttempname,".pdb",sep="",collapse="");
          			if(file.exists(foldxoldpdb))
          			{
          			file.rename(from=foldxoldpdb,to=foldxnewpdb);
          			file.copy(from=foldxnewpdb,to=htmldir);
                file.remove(foldxnewpdb);
          			}
          		stab=round(as.numeric(foldxres),2);
          		stablab="";
          		if(stab>0){stablab="D";}else{stablab="I"};
          		temp.result=c(mutmat[mutrow,],stab,stablab);
          		print(temp.result);
          		}else{
          		temp.result=c(mutmat[mutrow,],"-","-");
          		}
          #	},function(e){print(e);cat("    FoldX failed....Please recheck the inputs or mail us\n",file=logfile,sep="",eol="",append=TRUE); next;},finally=print("Error exception"));
          }
          	
          	if(usetool==1)
          	{
          #  tryCatch(
           # {
          	# AUTOMUTE
          	print("Running automute");
          	cat("    Running Automute....\n",file=logfile,sep="",eol="",append=TRUE);
          	setwd(automutdir);
          	print(automutdir);
          	mutvec=paste(mutmat[mutrow,3],mutmat[mutrow,2],mutmat[mutrow,4],sep="");
          	mutvec=cbind(paste(pdbpref,mutmat[mutrow,1],sep=""),mutvec);
          	print(mutvec);
          	temp.result=c(mutmat[mutrow,],automute(pdbpref,mutvec,automute_mode=automute_mode));	
          	print(temp.result);
          	setwd(orig_wd);
          	print(getwd());	
            #},error=function(e){  cat("    Automute failed....Please recheck the inputs or mail us\n",file=logfile,sep="",eol="",append=TRUE); next;},finally=print(""));
            
            }
          	
          	if(usetool==2)
          	{
          	# I-mutant
          #    tryCatch(
          #{
          	print("Running i-mutant");
          	cat("    Running I-mutant....\n",file=logfile,sep="",eol="",append=TRUE);
          	print(getwd());
          	
          	imuttemp=iMut(protengpref,protengdssp,mutmat[mutrow,1],mutmat[mutrow,2],mutmat[mutrow,4],imph,imtemp);
          
          	print(str(imuttemp[2]))
          	imuttemp1=imuttemp;
          	imuttemp1[1]=imuttemp[2];
          	
          	if(imuttemp[2]<0){imuttemp1[2]="D"}
          	if(imuttemp[2]>0){imuttemp1[2]="I"}
          	if(imuttemp[2]==0){imuttemp1[2]="N"}
          
          	temp.result=c(mutmat[mutrow,],imuttemp1);
          #},error=function(e){  cat("    I-mutant failed....Please recheck the inputs or mail us\n",file=logfile,sep="",eol="",append=TRUE); next;},finally=print(""));
          }
          	
          	if(usetool==3)
          	{
          	# Mupro
            #  tryCatch(
          #{
          	print("Running mupro");
          	cat("    Running Mupro....\n",file=logfile,sep="",eol="",append=TRUE);
          	mupro.result=mupro(pdb,mutmat[mutrow,1],mutmat[mutrow,2],mutmat[mutrow,3],mutmat[mutrow,4],muproRunfile);
          	print(mupro.result);
          	if(mupro.result[1]=="DECREASE"){mupro.result[1]="D"}
          	if(mupro.result[1]=="INCREASE"){mupro.result[1]="I"}
          	mupro.result1=c(mupro.result[2],mupro.result[1]);
          	temp.result=c(mutmat[mutrow,],mupro.result1);
          	print(temp.result);
          #	},error=function(e){  cat("    Mupro failed....Please recheck the inputs or mail us\n",file=logfile,sep="",eol="",append=TRUE); next;},finally=print(""));
          }
          
          wildcons1="-";
          
          if(conserv & modeofop=="bypdb")
          {
                      consfile1=paste(consdir,pdbid,"_",mutmat[mutrow,1],".txt",sep="",collapse=""); # for cys1
                      
                      if(file.exists(consfile1))
                      {
                        conscmd=paste("cat ",consfile1," | grep ",aa123(mutmat[mutrow,3]),mutmat[mutrow,2],":",mutmat[mutrow,1],sep="",collapse="");
                        print(conscmd);
                      #  tryCatch({
                          wildcons1=system(conscmd,ignore.stdout=FALSE,ignore.stderr=FALSE,intern=TRUE);
                      #  },error=function(e){next;},finally=print("Error CS analysis mupro"));
                      }
                      print(nchar(wildcons1));
                      if(!is.list(attributes(wildcons1)))
                      {
                        wildcons1=substr(wildcons1,32,32)
                      }
          }
          
          if(conserv & modeofop=="byupload")
          {
            consfile=paste(qtag,chmat[which(chmat[,1]==mutmat[mutrow,1]),2],"_map.txt",sep="",collapse="");
            print("INSIDE FETCH");
            print(consfile);
            print(file.exists(consfile));
            if(file.exists(consfile))
            {
              conscmd=paste("grep ",chmat[which(chmat[,1]==mutmat[mutrow,1]),2],mutmat[mutrow,2],aa123(mutmat[mutrow,3])," ",consfile," | awk '{print $5}'",sep="",collapse="");
              print(conscmd);
              consres="";
              #tryCatch({
              consres=system(conscmd,ignore.stdout=FALSE,ignore.stderr=FALSE,intern=TRUE);
              #},error=function(e){});
              print(consres);
              
              if(nchar(consres)>0)
              {
                wildcons1=consres;
                wlink=wildcons1;
                #wlink.name=paste(hlinkpath,prefix,"/",subpref1,"/",consfile,sep="",collapse="");
                #wlink=paste("<a href=\"",wlink.name,"\" target=\"_blank\" >",wildcons,"</a>");
              }
            }
          }
          
          	res.foldx=rbind(res.foldx,temp.result);
          	print(temp.result);
          	cat(c(temp.result,wildcons1),file=txtoutfile,sep="\t",eol="\n",append=TRUE);
          	
          	hlink.name=paste(hlinkpath,prefix,"/HTML/",foldxnewpdb,sep="",collapse="");
          	print(hlink.name);
          	hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >","mutantpdb","</a>");
          	print(hlink);
          
          	write2Html.tableBody(matrix(c(hlink,temp.result,wildcons1),nrow=1),htmloutfile,TRUE,"<td>");
          }

}
}

# Remove temp files from automute directory.
if(usetool==1)
{
    file1=paste(automutdir,"dinputs.txt",sep="",collapse="");
    file2=paste(automutdir,"doutputs.txt",sep="",collapse="");
    file3=paste(automutdir,"test_vector.arff",sep="",collapse="");
    #print(file1);
    if(file.exists(file1)){unlink(file1)}
    if(file.exists(file2)){unlink(file2)}
    if(file.exists(file3)){unlink(file3)}
    if(file.exists(pdbfile)){unlink(pdbfile)}
}

write2Html("\n</table>\n",htmloutfile,append=TRUE);
write2Html("\n</html>\n",htmloutfile,append=TRUE);

#cat("Your job is finished\nThank you for submitting\n",file=logfile,sep="",eol="");
cat("100",file=statusfile,sep="",eol="");
print("THANKS");

#unlink(pdbfile);

},
error=function(e)
{
  print("--ERROR--");
  cat("-111",file=statusfile,sep="",eol="");
}

)



