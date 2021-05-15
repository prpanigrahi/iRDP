# Test_catpi

prefix="catpi_test";
####################### SET All global directories
script_dir="/home/priyabrata/Work/ProtTS/SASMAT_ason_29sept2013_before4th5thPres/Script_dir/";
# It has to be shorter coz procheck fails if a longer path comes, donno why.
work_dir="/home/priyabrata/Work/ProtTS/workdir/";
#pdb_dir="/data1/pdb";
pdb_dir="/data1/testing_dataset";
output_dir="/home/priyabrata/Work/ProtTS/SASMAT_ason_29sept2013_before4th5thPres/CompStabAna/out/";


# Load all Libraries and sources
loadsource=function(script_dir,file)
{
	source(paste(script_dir,file,sep="",collapse=""));
}

library("bio3d");
library("igraph");
loadsource(script_dir,"PATH.r");
loadsource(script_dir,"general.r");
loadsource(script_dir,"dssp.r");
loadsource(script_dir,"naccess.r");
loadsource(script_dir,"catpi.r");
######################## Set the working directory
work_dir=paste(work_dir,prefix,sep="",collapse="");
createDir(work_dir);
setwd(work_dir);
calc_uniqch=1
# To use capture
calc_catpi=1;
use_capture=0;
catpimode=1; # CATion-pi mode=1: pdb 2:upload
catpicut=6.0;

filenames=list.files(pdb_dir,pattern=".pdb$");
#filenames=filenames[50:length(filenames)]

# All output directory paths
# Prefix: For a given run (i.e. for a set of pdb files given in pdb directory) the output root folder name

dirpath=paste(output_dir,prefix,sep="",collapse="");
dirpath.catpi=paste(output_dir,prefix,"/","Catpi",sep="",collapse="");
sumpath.catpi=paste(dirpath.catpi,"/","Catpi_summary.txt",sep="",collapse="");
detpath.catpi=paste(dirpath.catpi,"/","Catpi_detail.txt",sep="",collapse="");
cat("",file=detpath.catpi);

# create root directory
createDir(dirpath);

if(calc_catpi)
{
	createDir(dirpath.catpi);
	if(use_capture)
	{ #use capture to calcualte catpi interaction
	cat("Filename\tTotal\tEesLow\tEesHigh\tAvg\tRF\tRY\tRW\tKF\tKY\tKW\n",file=sumpath.catpi);
	}else{ #use catpi() to calculate catpi interaction.
	cat("Filename\tTotal\tpIntra\tpInter\tpKF\tpKY\tpKW\tpRF\tpRY\tpRW\tpB\tpE\tpIso\tTotNet\tNetDet\n",file=sumpath.catpi);
	}
}


for(f in filenames)
	{
	pdb="";nacc="";DSSP="";
	count=which(filenames==f);

	# All temp pdbs names
	dsspprefix=paste(prefix,"_dssp",count,sep="",collapse="");
	naccprefix=paste(prefix,"_nacc",count,sep="",collapse="");
	hbplusprefix=paste(prefix,"_hbplus",count,sep="",collapse="");
	promotifprefix=paste(prefix,"_promotif",count,sep="",collapse="");
	procheckprefix=paste(prefix,"_procheck",count,sep="",collapse="");
	captureprefix=paste(prefix,"_capture",count,sep="",collapse="");

	print(paste("Processing ProtTS for file ",f,sep="",collapse=""));
	pdbfile=paste(pdb_dir,"/",f,sep="");
	print(pdbfile);
	if(!file.exists(pdbfile)){print(c("PDB file not found",pdbfile));next; }
	pdb=read.pdb(pdbfile,maxlines=100000,verbose=FALSE); # read the pdb file
		if(calc_uniqch)
		{ #Do uniquechain calculation 
		uniqch=uniqueChain(pdb);
		if(!sum(is.na(uniqch))>0){pdb=trim.pdb(pdb,inds=atom.select(pdb,chain=uniqch,verbose=FALSE));}
		}
	nacc=naccess(pdb,naccpath,prefix=naccprefix); # Runs naccess
	DSSP=dssp_new(pdb,exepath=dssppath,prefix=dsspprefix); #Runs DSSP
	status.DSSP=0;if(is.list(DSSP)){status.DSSP=1;}
	status.nacc=0;if(is.list(nacc)){status.nacc=1;}
	
	if(calc_catpi)
	{
		print("Running catpi");
		if(use_capture)
		{
		}else{
		catpiList=callCatpi(pdb,DSSP,nacc,status.DSSP,status.nacc);
			if(is.matrix(catpiList$catpimat))
			{
			catpisum=catpiSummary(catpiList$catpimat);
			cat(paste(c(f,catpisum),collapse="\t"),file=sumpath.catpi,eol="\n",append=TRUE);
			cat(f,file=detpath.catpi,append=TRUE); # Create individual ionpair table file
			cat("\n-----\n\nChain1\tResno1\tResid1\tChain2\tResno2\tResid2\tResid1_Acc\tResid2_Acc\tResid1_SS\tResid2_SS\n",file=detpath.catpi,append=TRUE); # Create individual ionpair table file
			write.table(catpiList$catpimat,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=detpath.catpi);
			# Note: If no network then print_ionnetwork fn does not print network details :)
			cat("\n\nNetwork Summary\n-------------------------------------\n",file=detpath.catpi,sep="",eol="",append=TRUE)
			print_catpinetwork(catpiList$catpimat,catpiList$dcatpimat,detpath.catpi);
			}else{
			cat(paste(c(f,rep("-",times=14)),collapse="\t"),file=sumpath.catpi,eol="\n",append=TRUE);
			}		
		}
	}

}

