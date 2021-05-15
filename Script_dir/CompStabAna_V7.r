# A new version that incorporates file upload features

tryCatch(
{
now <- Sys.time()
Sys.setenv(http_proxy=""); # Env variable set

########################################################### SET All global directories #################################################################
script_dir="/var/www/Script_dir/";
work_dir="/var/www/workdir/"; # It has to be shorter coz procheck fails if a longer path comes, donno why.
pdb_dir="/data/ent";
output_dir="/data/out/";
hlinkpath="/out/";
upload_dir="/var/www/upload";
totfile=100; #maximum number of files to process
subpref1="HTML";
subpref2="TEXT";

########################################################### Load all Libraries and sources ########################################################
loadsource=function(script_dir,file)
{
	source(paste(script_dir,file,sep="",collapse=""));
}

library("bio3d",quietly=TRUE,verbose=FALSE);
library("parallel",quietly=TRUE,verbose=FALSE);
library("igraph",quietly=TRUE,verbose=FALSE);
loadsource(script_dir,"PATH.r");
loadsource(script_dir,"dssp.r");
loadsource(script_dir,"naccess.r");
loadsource(script_dir,"aacomp.r");
loadsource(script_dir,"secstr.r");
loadsource(script_dir,"ionic.r");
loadsource(script_dir,"general.r");
loadsource(script_dir,"aroaro.r");
loadsource(script_dir,"aroS.r");
loadsource(script_dir,"deamidate.r");
loadsource(script_dir,"hdipole.r");
loadsource(script_dir,"bturn.r");print("*");
loadsource(script_dir,"procheck.r");
loadsource(script_dir,"disulfide.r");
loadsource(script_dir,"hbond.r");
loadsource(script_dir,"asa.r");
loadsource(script_dir,"proline.r");
loadsource(script_dir,"Hphob.r");
loadsource(script_dir,"catpi.r");
loadsource(script_dir,"intprofile.r");
loadsource(script_dir,"writehtml.r");
loadsource(script_dir,"findgeo.r");
loadsource(script_dir,"foldx.r");

########################################################### Fetch arguments & parse option file ########################################################
args=commandArgs(trailingOnly = TRUE);
optionfile=args[1]; # Arg1 is option file
options=read.table(optionfile,stringsAsFactors=FALSE) # Read the option file
prefix=options[which(options$V1=="prefix"),2] # First row of option file is the session id to be used as prefix
print(c("Prefix entered",prefix));

modeofop=""; #mode of operation [bypdb, byupload]
if(options[which(options$V1=="pdbmode"),2]=="args")
{
modeofop="bypdb";
}else{
modeofop="byupload";
}

filenames="";

if(modeofop=="bypdb")
{
separate=",";
pdbids=paste(readLines(args[2]),sep="",collapse="");
pdbids=gsub("\n","",pdbids); #chomp any \n
pdbids=gsub(" ","",pdbids); #chomp any space
filenames=tolower(unique(paste("pdb",unlist(strsplit(pdbids,separate)),".ent",sep="")));
#print(filenames); #print(length(filenames));

	if(length(filenames)>totfile)
	{
	  filenames=filenames[1:totfile]; #only 100 files will be accepted
	}
}else{
	uploaddir=paste(upload_dir,"/",prefix,sep="",collapse="");
	filenames=list.files(uploaddir,pattern=".ent$|.pdb$"); #accept pdb and ent files only	
	pdb_dir=uploaddir;
}

print(modeofop);
print(filenames);

######################## Set the working directory
work_dir=paste(work_dir,prefix,sep="",collapse="");
createDir(work_dir);
setwd(work_dir);



# Setting sfor which programs to run
calc_aacomp=as.numeric(options[which(options$V1=="aa"),2]);
calc_secstr=as.numeric(options[which(options$V1=="ss"),2]);
calc_secstrdist=as.numeric(options[which(options$V1=="aass"),2]);
calc_disul=as.numeric(options[which(options$V1=="dis"),2]);
calc_ionic=as.numeric(options[which(options$V1=="ip"),2]);

calc_hbond=as.numeric(options[which(options$V1=="hb"),2]);
calc_aroaro=as.numeric(options[which(options$V1=="ap"),2]);
calc_aros=as.numeric(options[which(options$V1=="aros"),2]);
calc_catpi=as.numeric(options[which(options$V1=="catpi"),2]);
calc_hphob=as.numeric(options[which(options$V1=="hp"),2]);
calc_asa=as.numeric(options[which(options$V1=="asa"),2]);
calc_hdipo=as.numeric(options[which(options$V1=="hdipo"),2]);
calc_proline=as.numeric(options[which(options$V1=="pro"),2]);
calc_left=as.numeric(options[which(options$V1=="left"),2]);
calc_deam=as.numeric(options[which(options$V1=="deam"),2]);
calc_deam=as.numeric(options[which(options$V1=="deam"),2]);
calc_metal=as.numeric(options[which(options$V1=="metal"),2]);
calc_gibbs=as.numeric(options[which(options$V1=="stab"),2]);
calc_uniqch=as.numeric(options[which(options$V1=="uniq"),2]);

# Program parameters
ioncut=as.numeric(options[which(options$V1=="ioncut"),2]);
aroarolow=as.numeric(options[which(options$V1=="aroarolow"),2]);
aroarohigh=as.numeric(options[which(options$V1=="aroarohigh"),2]);
usedihed=as.numeric(options[which(options$V1=="usedihed"),2]);
aroarodihed=as.numeric(options[which(options$V1=="aroarodihed"),2]);
aroscut=as.numeric(options[which(options$V1=="aroscut"),2]);
hpcut=as.numeric(options[which(options$V1=="hpcut"),2]);
acclow=as.numeric(options[which(options$V1=="acclow"),2]);
catpicut=as.numeric(options[which(options$V1=="catpicut"),2]);

#count=79678; #filenames=list.files(pdb_dir,pattern=".ent$"); #filenames=filenames[count:length(filenames)]

# All output directory paths
# Prefix: For a given run (i.e. for a set of pdb files given in pdb directory) the output root folder name


makedirPath=function(output_dir,prefix,subpref,suffix)
{
print("call to makedirPath");
dirpath=paste(output_dir,prefix,sep="",collapse="");
subdirpath=paste(output_dir,prefix,"/",subpref,sep="",collapse="");
dirpath.aacomp=paste(output_dir,prefix,"/",subpref,"/","AminoComp",sep="",collapse="");
dirpath.secstr=paste(output_dir,prefix,"/",subpref,"/","SecStrCompo",sep="",collapse="");
dirpath.secstrdist=paste(output_dir,prefix,"/",subpref,"/","SecStrDistr",sep="",collapse="");
dirpath.disul=paste(output_dir,prefix,"/",subpref,"/","Disulfide",sep="",collapse="");
dirpath.ionic=paste(output_dir,prefix,"/",subpref,"/","Ionpair",sep="",collapse="");
dirpath.aroaro=paste(output_dir,prefix,"/",subpref,"/","AroAro",sep="",collapse="");
dirpath.aros=paste(output_dir,prefix,"/",subpref,"/","AroS",sep="",collapse="");
dirpath.hbond=paste(output_dir,prefix,"/",subpref,"/","Hbond",sep="",collapse="");
dirpath.catpi=paste(output_dir,prefix,"/",subpref,"/","Catpi",sep="",collapse="");
dirpath.asa=paste(output_dir,prefix,"/",subpref,"/","ASA",sep="",collapse="");
dirpath.proline=paste(output_dir,prefix,"/",subpref,"/","Proline",sep="",collapse="");
dirpath.hdipo=paste(output_dir,prefix,"/",subpref,"/","HelixDipole",sep="",collapse="");
dirpath.deam=paste(output_dir,prefix,"/",subpref,"/","Deamidation",sep="",collapse="");
dirpath.left=paste(output_dir,prefix,"/",subpref,"/","ConfStrain",sep="",collapse="");
dirpath.hphob=paste(output_dir,prefix,"/",subpref,"/","Hydrophobic",sep="",collapse="");
dirpath.summ=paste(output_dir,prefix,"/",subpref,"/","Summary",sep="",collapse="");
dirpath.metal=paste(output_dir,prefix,"/",subpref,"/","Metal",sep="",collapse="");
dirpath.gibbs=paste(output_dir,prefix,"/",subpref,"/","Gibbs",sep="",collapse="");

templist=list(
dirpath=dirpath,
subdirpath=subdirpath,
dirpath.aacomp=dirpath.aacomp,
dirpath.secstr=dirpath.secstr,
dirpath.secstrdist=dirpath.secstrdist,
dirpath.disul=dirpath.disul,
dirpath.ionic=dirpath.ionic,
dirpath.aroaro=dirpath.aroaro,
dirpath.aros=dirpath.aros,
dirpath.hbond=dirpath.hbond,
dirpath.catpi=dirpath.catpi,
dirpath.asa=dirpath.asa,
dirpath.proline=dirpath.proline,
dirpath.hdipo=dirpath.hdipo,
dirpath.deam=dirpath.deam,
dirpath.left=dirpath.left,
dirpath.hphob=dirpath.hphob,
dirpath.summ=dirpath.summ,
dirpath.metal=dirpath.metal,
dirpath.gibbs=dirpath.gibbs,
sumpath.aacomp=paste(dirpath.aacomp,"/","aacomp_summary.",suffix,sep="",collapse=""),
sumpath.secstr=paste(dirpath.secstr,"/","secstr_summary.",suffix,sep="",collapse=""),
sumpath.hcomp=paste(dirpath.secstrdist,"/","hcomp_summary.",suffix,sep="",collapse=""),
sumpath.scomp=paste(dirpath.secstrdist,"/","scomp_summary.",suffix,sep="",collapse=""),
sumpath.tcomp=paste(dirpath.secstrdist,"/","tcomp_summary.",suffix,sep="",collapse=""),
sumpath.ccomp=paste(dirpath.secstrdist,"/","ccomp_summary.",suffix,sep="",collapse=""),
sumpath.disul=paste(dirpath.disul,"/","Disulfide_summary.",suffix,sep="",collapse=""),
sumpath.ionic=paste(dirpath.ionic,"/","Ionpair_summary.",suffix,sep="",collapse=""),
sumpath.aroaro=paste(dirpath.aroaro,"/","Aromatic_summary.",suffix,sep="",collapse=""),
sumpath.hbond=paste(dirpath.hbond,"/","Hbond_summary.",suffix,sep="",collapse=""),
sumpath.aros=paste(dirpath.aros,"/","AroS_summary.",suffix,sep="",collapse=""),
sumpath.catpi=paste(dirpath.catpi,"/","Catpi_summary.",suffix,sep="",collapse=""),
sumpath.asa=paste(dirpath.asa,"/","ASA_summary.",suffix,sep="",collapse=""),
sumpath.proline=paste(dirpath.proline,"/","Proline_summary.",suffix,sep="",collapse=""),
sumpath.hdipo=paste(dirpath.hdipo,"/","Hdipole_summary.",suffix,sep="",collapse=""),
sumpath.deam=paste(dirpath.deam,"/","Deamidation_summary.",suffix,sep="",collapse=""),
sumpath.left=paste(dirpath.left,"/","ConfStrain_summary.",suffix,sep="",collapse=""),
sumpath.hphob=paste(dirpath.hphob,"/","Hydrophobic_summary.",suffix,sep="",collapse=""),
sumpath.metal=paste(dirpath.metal,"/","Metal_summary.",suffix,sep="",collapse=""),
sumpath.gibbs=paste(dirpath.gibbs,"/","Gibbs_summary.",suffix,sep="",collapse="")
);
return(templist);
}

htmldir=makedirPath(output_dir,prefix,subpref1,"html"); # a list
txtdir=makedirPath(output_dir,prefix,subpref2,"txt"); # a list

#print("***"); #print(txtdir$subdirpath);

# hlink path
hlinktempdir=makedirPath(hlinkpath,prefix,subpref1,"html");
hlink.ionic=paste(hlinktempdir$dirpath.ionic,"/Details/",sep="",collapse="");
hlink.aroaro=paste(hlinktempdir$dirpath.aroaro,"/Details/",sep="",collapse="");
hlink.aros=paste(hlinktempdir$dirpath.aros,"/Details/",sep="",collapse="");
hlink.disul=paste(hlinktempdir$dirpath.disul,"/Details/",sep="",collapse="");
hlink.hbond=paste(hlinktempdir$dirpath.hbond,"/Details/",sep="",collapse="");
hlink.catpi=paste(hlinktempdir$dirpath.catpi,"/Details/",sep="",collapse="");
hlink.hphob=paste(hlinktempdir$dirpath.hphob,"/Details/",sep="",collapse="");
hlink.proline=paste(hlinktempdir$dirpath.proline,"/Details/",sep="",collapse="");
hlink.hdipo=paste(hlinktempdir$dirpath.hdipo,"/Details/",sep="",collapse="");
hlink.deam=paste(hlinktempdir$dirpath.deam,"/Details/",sep="",collapse="");
hlink.left=paste(hlinktempdir$dirpath.left,"/Details/",sep="",collapse="");
hlink.asa=paste(hlinktempdir$dirpath.asa,"/Details/",sep="",collapse="");
hlink.metal=paste(hlinktempdir$dirpath.metal,"/Details/",sep="",collapse="");
hlink.secstr=paste(hlinktempdir$dirpath.secstr,"/Details/",sep="",collapse="");

# create root directory
createDir(htmldir$dirpath);
createDir(htmldir$subdirpath);
createDir(txtdir$subdirpath);

#Status file
readmefile=paste(txtdir$subdirpath,"/README.txt",sep="",collapse="");
statusfile=paste(htmldir$subdirpath,"/status.txt",sep="",collapse="");
logfile=paste(htmldir$subdirpath,"/log.txt",sep="",collapse="");
cat(paste(c("Total number of files to be processed:",length(filenames),"\n"),sep="",collapse=""),file=logfile,sep="",eol="");

# Summary HTML and TEXT files
sumhtmlfile=paste(htmldir$subdirpath,"/Summary.html",sep="",collapse="");
sumtxtfile=paste(txtdir$subdirpath,"/Summary.txt",sep="",collapse="");
#Summary vector
sumvec=c(Filename="-",Chain="-",Length="-",Aro="-",NQST="-",RKH="-",DE="-","RtoK"="-",TotalIP="-",TotalAAI="-",TotalASI="-",TotalCPI="-",TotalDB="-",TotalHB="-",TotalHP="-",TotalPro="-",TotalBt2P="-",TotalNcapP="-","NP/P"="-",TotalTL="-",St.Helix="-",TotalCS="-",TotalMetal="-",TotalGibbsEnergy="-"); 

print_HTMLglobal_sum_start(sumhtmlfile);
cat("Summary of Parameters",sep="\t",eol="\n\n",file=sumtxtfile); 
cat(names(sumvec),sep="\t",eol="\n",file=sumtxtfile,append=T); 
 
cat("The result folder contains many subfolder depending upon the choice of the input features while submission.
Below are the details of each folder

AminoComp  : Results of Amino acid composition
SecStrCompo	: Results of Secondary structure composition. It contains Details folder containing raw dssp output files
SecStrDistr	: Results of Helix/Strand/Turn/Coil composition. It contains 4 files
Disulfide	: Results of disulfide bond analysis.
Ionpair		: Results of Ionnic interactions.
AroAro		: Results of Aromatic-aromatic interaction analysis.
AroS		: Results of Aromatic-sulphur interaction analysis.
Catpi		: Results of Cation-pi interaction analysis.
Hbond		: Results of Hydrogen bond analysis.
Hydrophobic	: Results of Hydrophobic interaction analysis.
ASA		: Results of Accessible surface Area analysis.
Proline		: Results of Proline residue profile analysis.
Deamidation	: Results of Thermolabile bond analysis.
ConfStrain	: Results of Conformationally strained residue profile analysis.
HelixDipole	: Results of Helix Dipole Stabilization profile analysis.
Metal		: Results of Metal binding analysis.
Gibbs		: Results of Gibbs free energy of folding analysis.

In each folder, you may find a sub-folder (Details) containing details of each analysis\n",file=readmefile,sep="",eol="");



if(calc_aacomp)
{
	createDir(htmldir$dirpath.aacomp);
	print_HTMLaacomp_sum_start(htmldir$sumpath.aacomp);
	createDir(txtdir$dirpath.aacomp);
	cat("Amino acid composition\n----------------------\nFilename\tChain\tLength\tV\tI\tL\tM\tF\tW\tY\tS\tT\tN\tQ\tH\tK\tR\tD\tE\tA\tG\tP\tC\tX\tAro[FWY]\tNQST\tRKH\tDE\tR/K\n",file=txtdir$sumpath.aacomp); # Create aacomp_summary.txt file
	
}
if(calc_secstr)
{
	createDir(htmldir$dirpath.secstr);
	print_HTMLsecstr_sum_start(htmldir$sumpath.secstr);
	createDir(txtdir$dirpath.secstr);
	cat("Secondary Structure Composition [Estimated by DSSP]\n---------------------------------------------------\nFilename\tLength\tB\tC\tE\tG\tH\tI\tS\tT\n",file=txtdir$sumpath.secstr); # Create secstr summary file
	
}
if(calc_secstrdist)
{
	createDir(htmldir$dirpath.secstrdist);
	print_HTMLsecstrdist_sum_start(htmldir$sumpath.hcomp,"Helix");
	print_HTMLsecstrdist_sum_start(htmldir$sumpath.scomp,"Strand");
	print_HTMLsecstrdist_sum_start(htmldir$sumpath.tcomp,"Turn");
	print_HTMLsecstrdist_sum_start(htmldir$sumpath.ccomp,"Coil");
	createDir(txtdir$dirpath.secstrdist);
	cat("Helix [DSSP assignments:H+G+I] Composition\n---------------------------------------\nFilename\tSS\tT.Res\tV\tI\tL\tM\tF\tW\tY\tS\tT\tN\tQ\tH\tK\tR\tD\tE\tA\tG\tP\tC\tX\n",file=txtdir$sumpath.hcomp); # Create hcomp summary.txt file
	cat("Strand [DSSP assignments:B+E] Composition\n---------------------------------------\nFilename\tSS\tT.Res\tV\tI\tL\tM\tF\tW\tY\tS\tT\tN\tQ\tH\tK\tR\tD\tE\tA\tG\tP\tC\tX\n",file=txtdir$sumpath.scomp); # Create scomp summary.txt file
	cat("Turn [DSSP assignments:T] Composition\n---------------------------------------\nFilename\tSS\tT.Res\tV\tI\tL\tM\tF\tW\tY\tS\tT\tN\tQ\tH\tK\tR\tD\tE\tA\tG\tP\tC\tX\n",file=txtdir$sumpath.tcomp); # Create tcomp summary.txt file
	cat("Coil [DSSP assignments:S+C] Composition\n---------------------------------------\nFilename\tSS\tT.Res\tV\tI\tL\tM\tF\tW\tY\tS\tT\tN\tQ\tH\tK\tR\tD\tE\tA\tG\tP\tC\tX\n",file=txtdir$sumpath.ccomp); # Create ccomp summary.txt file
	
}
if(calc_disul)
{
	createDir(htmldir$dirpath.disul);
	print_HTMLdis_sum_start(htmldir$sumpath.disul);
	createDir(txtdir$dirpath.disul);
	cat("Disulfide bond (DB) summary\n-------------------------------------------\nFilename\tTotal\tIntraCh\tInterCh\tB\tE\t0_10\t10_20\t20_30\t30_40\t40_50\t>50\tPP\tPN\tNN\n",file=txtdir$sumpath.disul);
}
if(calc_ionic)
{
	createDir(htmldir$dirpath.ionic);
	print_HTMLionic_sum_start(htmldir$sumpath.ionic);
	createDir(txtdir$dirpath.ionic);
	cat("Ionic Interaction (IP) summary\n-----------------------------------------\nFilename\tTotal\tIntraCh\tInterCh\tD\tE\tR\tK\tH\tDK\tDR\tDH\tEK\tER\tEH\tB\tE\tIso\tT.Net\tNet_Details\n",file=txtdir$sumpath.ionic);
}
if(calc_aroaro)
{
	createDir(htmldir$dirpath.aroaro);
	print_HTMLaroaro_sum_start(htmldir$sumpath.aroaro);
	createDir(txtdir$dirpath.aroaro);
	cat("Aromatic-Aromatic interaction (AAI) summary\n---------------------------------------------------\nFilename\tTotal\tIntraCh\tInterCh\tF\tY\tW\tFF\tFY\tFW\tYY\tYW\tWW\tB\tE\tIso\tT.Net\tNet_Details\n",file=txtdir$sumpath.aroaro);
}

if(calc_hbond)
{
	createDir(htmldir$dirpath.hbond);
	print_HTMLhbond_sum_start(htmldir$sumpath.hbond);
	createDir(txtdir$dirpath.hbond);
	cat("Hydrogen bond (HB) summary\n---------------------------------------\nFilename\tTotal\tIntraCh\tInterCh\tMM\tMS\tSM\tSS\tCNHB\tNNHB\n",file=txtdir$sumpath.hbond);
	
}

if(calc_aros)
{
	createDir(htmldir$dirpath.aros);
	print_HTMLaros_sum_start(htmldir$sumpath.aros);
	createDir(txtdir$dirpath.aros);
	cat("Aromatic Sulphur interaction (ASI) summary\n----------------------------------------------\nFilename\tTotal\tIntraCh\tInterCh\tF\tY\tW\tC\tM\tFC\tYC\tWC\tFM\tYM\tWM\tB\tE\tIso\tT.Net\tNet_Details\n",file=txtdir$sumpath.aros);
	
}

if(calc_asa)
{
	createDir(htmldir$dirpath.asa);
	print_HTMLasa_sum_start(htmldir$sumpath.asa);
	createDir(txtdir$dirpath.asa);
	cat("Solvent Accessible surface area (ASA) summary\n----------------------------------------------------\nFilename\tAll-atoms\tTotal-Side\tMain-Chain\tNon-polar(NP)\tPolar(P)\tNP/P\tASA_C\tASA_N\tASA_O\tASA_S\n",file=txtdir$sumpath.asa);
}

if(calc_catpi)
{
	createDir(htmldir$dirpath.catpi);
	createDir(txtdir$dirpath.catpi);
	cat("Cation-pi interaction (CPI) summary\n-------------------------------------------\nFilename\tTotal\tIntraCh\tInterCh\tKF\tKY\tKW\tRF\tRY\tRW\tB\tE\tIso\tT.Net\tNet_Details\n",file=txtdir$sumpath.catpi);
	print_HTMLcatpi_sum_start(htmldir$sumpath.catpi);
}

if(calc_proline)
{
	createDir(htmldir$dirpath.proline);
	print_HTMLproline_sum_start(htmldir$sumpath.proline);
	createDir(txtdir$dirpath.proline);
	cat("Proline residue profile\n--------------------------------\nFilename\tT.Proline\tHelix\tStrand\tTurn\tCoil\tB\tE\tT.BtP\tT.Bt2P\tT.NcapP\n",file=txtdir$sumpath.proline);
	
}

if(calc_hdipo)
{
	createDir(htmldir$dirpath.hdipo);
	print_HTMLhdipo_sum_start(htmldir$sumpath.hdipo);
	createDir(txtdir$dirpath.hdipo);
	cat("Helix dipole stabilization profile\n-----------------------------------------\nFilename\tT.Helix\tSt.Helix\tN.St.Helix\tC.St.Helix\tNC.St.Helix\tSt_N-2\tSt_N-1\tSt_N\tSt_N+1\tSt_N+2\tSt_C-2\tSt_C-1\tSt_C\tSt_C+1\tSt_C+2\n",file=txtdir$sumpath.hdipo);
	
}

if(calc_deam)
{
	createDir(htmldir$dirpath.deam);
	print_HTMLdeam_sum_start(htmldir$sumpath.deam);
	createDir(txtdir$dirpath.deam);
	cat("Thermolabile bond profile\n---------------------------------------------\nFilename\tTotal\tNG\tNA\tNS\tQG\tQA\tQS\t<4\n",file=txtdir$sumpath.deam);
}

if(calc_left)
{
	createDir(htmldir$dirpath.left);
	print_HTMLleft_sum_start(htmldir$sumpath.left);
	createDir(txtdir$dirpath.left);
	cat("Conformationally strained (CS) residue profile\n-------------------------------------------\nFilename\tTotal\tL_region\tl_region\tla_region\tFrequency\n",file=txtdir$sumpath.left);
	
}

if(calc_hphob)
{
	createDir(htmldir$dirpath.hphob);
	print_HTMLhphob_sum_start(htmldir$sumpath.hphob);
	createDir(txtdir$dirpath.hphob);
	cat("Hydrophobic interaction (HP) summary\n--------------------------------------------\nFilename\tTotal\tIntraCh\tInterCh\tB\tE\n",file=txtdir$sumpath.hphob);
}

if(calc_metal)
{
	createDir(htmldir$dirpath.metal);
	print_HTMLmetal_sum_start(htmldir$sumpath.metal);
	createDir(txtdir$dirpath.metal);
	cat("Metal binding summary\n--------------------------------\nFilename\tTotal\tDetails\n",file=txtdir$sumpath.metal);
}

if(calc_gibbs)
{
	createDir(htmldir$dirpath.gibbs);
	print_HTMLgibbs_sum_start(htmldir$sumpath.gibbs);
	createDir(txtdir$dirpath.gibbs);
	cat("Gibbs energy of folding using FOLDX force field\n-------------------------------------------------------------\nFilename\tTE\tBB_HB\tSC_HB\tVDW\tElec\tSol_P\tSol_H\tVDW_C\tEnt_SC\tEnt_MC\tEnt_SL\tEnt_ML\tCis\tTor_C\tBB_C\tHD\tWB\tDB\tElec_kon\tPCB\tIE\tEnt_Cx\tN.Res\n",file=txtdir$sumpath.gibbs); # Create aacomp_summary.txt file
system(paste("cp ",foldxrotamer," .")); #copy the foldx rotamer library to current workdir
}
	
### Program begins here
filecount=0;
totalfile=length(filenames);
cat(paste(floor(filecount*100/totalfile), sep="",collapse=""),file=statusfile,sep="",eol="");

	for(f in filenames)
	{
		pdb="";
		nacc="";
		DSSP="";
		count=which(filenames==f);
		if(modeofop=="bypdb")
		{
		filepref=substr(f,4,7); # for ent format
		}else{
		ftemp=f;
		filepref=gsub(".ent|.pdb","",ftemp);
		}
		
		# All temp pdbs names
		dsspprefix=paste(filepref,"_dssp",count,sep="",collapse="");
		naccprefix=paste(filepref,"_nacc",count,sep="",collapse="");
		hbplusprefix=paste(filepref,"_hbplus",count,sep="",collapse="");
		promotifprefix=paste(filepref,"_promotif",count,sep="",collapse="");
		procheckprefix=paste(filepref,"_procheck",count,sep="",collapse="");
		captureprefix=paste(filepref,"_capture",count,sep="",collapse="");
		
		# Read the pdb file
		print(paste("Processing ProtTS for file ",f,sep="",collapse=""));
		cat(paste("\n",count,": Processing file ",f,"\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
		pdbfile=paste(pdb_dir,"/",f,sep="");
		print(pdbfile);
		
		if(!file.exists(pdbfile))
		{
		print("checking file in PDB");
		
		pdblines="";
		tryCatch({pdblines=readLines(paste("http://www.rcsb.org/pdb/files/",filepref,".pdb",sep="",collapse=""))},error = function(e){print("pdb cant download");
		cat(paste("\n!!! INVALID PDB !!! ",filepref,"\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
		}, finally = print("Hello"))

			if(length(grep("ATOM",pdblines))>0)
			{
			print("Downloaded file from PDB");
			cat(pdblines,sep="\n",eol="\n",file=f);
			pdbfile=f;
			}else{
			print(c("PDB file not found",pdbfile));
			cat(paste("\n!!! PDB FILE NOT FOUND !!! ",f,"\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
			next; 
			}
		}
		
    if(length(grep("^ATOM",readLines(pdbfile)))==0)
		{
		  cat(paste("\n!!! NO ATOM RECORDS !!! ",filenames[1],"\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
		  next; 
		}
    
		pdb=read.pdb(pdbfile,maxlines=1000000,verbose=FALSE); # read the pdb file
		
		#paste(pdb$seqres,collapse="");
		
		#Do uniquechain calculation if required
		uniqch=uniqueChain(pdb);
			
		print(uniqch);
		if(calc_uniqch)
		{
		cat(paste("Doing unique chain calculation ","\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE); 
		if(!sum(is.na(uniqch))>0)
		{
			inds=atom.select(pdb,chain=uniqch,verbose=FALSE); 
			if(length(inds$atom)<=1)
			{next;} 
			trimseq=NULL;
			#print(length(pdb$seqres))
			if(length(pdb$seqres)>0)
			{
			trimseq=pdb$seqres[names(pdb$seqres) %in% uniqch];
			}
			pdb=trim.pdb(pdb,inds=inds);
			#print(trimseq);
			if(length(trimseq)>0){pdb=list(atom=pdb$atom,xyz=pdb$xyz,seqres=trimseq)}
			#print(str(pdb));
			
		}
		}

		if(length(uniqch)==0){print(c("No protein chain present",pdbfile)); next;}
		if(nrow(pdb$atom)==0){print(c("No atoms records present",pdbfile)); next;}

		# perform DSSP and naccess calculations
		status.DSSP=0;status.nacc=0;
		if(calc_secstr | calc_secstrdist | calc_ionic | calc_aroaro | calc_aros | calc_disul | calc_catpi | calc_hphob | calc_asa | calc_proline | calc_hdipo | calc_deam | calc_left)
		{
		nacc=naccess(pdb,naccpath,prefix=naccprefix); # Runs naccess
		print("***");
		DSSP=dssp_new(pdb,exepath=dssppath,prefix=dsspprefix); #Runs DSSP
		if(is.list(DSSP)){status.DSSP=1;}
		if(is.list(nacc)){status.nacc=1;}
		}
		
		# Running amino acid composition
		if(calc_aacomp)
		{
		print("Running aacomp");
		aacomp=AAcomp(pdb,"",f);
		tempuniqch=paste(unique(pdb$atom[,"chain"]),sep="",collapse=",")
		cat(paste(c(filepref,tempuniqch,aacomp),sep="",collapse="\t"),file=txtdir$sumpath.aacomp,eol="\n",append=TRUE);
		write2Html.tableBody(matrix(c(filepref,tempuniqch,aacomp),nrow=1),htmldir$sumpath.aacomp,TRUE,"<td>");
		cat("   Finished aa composition analysis...\n",file=logfile,sep="",eol="",append=TRUE);
    sumvec[1:8]=c(filepref,tempuniqch,aacomp[c(1,23:27)])
		}
		# Secstr
		if(calc_secstr)
		{
		print("Running secstr");
			hlink.name=paste(hlink.secstr,dsspprefix,".dssp",sep="");
			hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");

			if(status.DSSP)
			{
			secstr=secstrCompo(DSSP,pdb);
		
				# create DSSP detail directory
				secstr.dirname=paste(htmldir$dirpath.secstr,"/Details",sep="",collapse="");
				secstr.dirnametxt=paste(txtdir$dirpath.secstr,"/Details",sep="",collapse="");
				createDir(secstr.dirname);
				createDir(secstr.dirnametxt);
					
				if(file.exists(paste(dsspprefix,".dssp",sep="")))
				{
				file.copy(from=paste(dsspprefix,".dssp",sep=""),to=secstr.dirname);
				file.copy(from=paste(dsspprefix,".dssp",sep=""),to=secstr.dirnametxt);
				}
				write2Html.tableBody(matrix(c(hlink,secstr),nrow=1),htmldir$sumpath.secstr,TRUE,"<td>");
				cat(paste(c(filepref,secstr),collapse="\t"),file=txtdir$sumpath.secstr,eol="\n",append=TRUE);
			}else{
			cat(paste(c(filepref,rep("-",times=9)),collapse="\t"),file=txtdir$sumpath.secstr,eol="\n",append=TRUE); # For those pdb DSSP fails, no secstr calculation
			write2Html.tableBody(matrix(c(filepref,rep("-",times=9)),nrow=1),htmldir$sumpath.secstr,TRUE,"<td>");
			}
		cat("   Finished ss composition analysis...\n",file=logfile,sep="",eol="",append=TRUE);	
		}
		
		# Secstr distribution
		if(calc_secstrdist)
		{
			print("Running secstr dist");
			if(status.DSSP)
			{
			htemp=secstr_aadistr(DSSP,c("H","G","I"));
			stemp=secstr_aadistr(DSSP,c("B","E"));
			ttemp=secstr_aadistr(DSSP,c("T"));
			ctemp=secstr_aadistr(DSSP,c("S","C"));
			
			cat(paste(c(filepref,htemp),collapse="\t"),file=txtdir$sumpath.hcomp,eol="\n",append=TRUE);
			cat(paste(c(filepref,stemp),collapse="\t"),file=txtdir$sumpath.scomp,eol="\n",append=TRUE);
			cat(paste(c(filepref,ttemp),collapse="\t"),file=txtdir$sumpath.tcomp,eol="\n",append=TRUE);
			cat(paste(c(filepref,ctemp),collapse="\t"),file=txtdir$sumpath.ccomp,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(filepref,htemp),nrow=1),htmldir$sumpath.hcomp,TRUE,"<td>");
			write2Html.tableBody(matrix(c(filepref,stemp),nrow=1),htmldir$sumpath.scomp,TRUE,"<td>");
			write2Html.tableBody(matrix(c(filepref,ttemp),nrow=1),htmldir$sumpath.tcomp,TRUE,"<td>");
			write2Html.tableBody(matrix(c(filepref,ctemp),nrow=1),htmldir$sumpath.ccomp,TRUE,"<td>");
			}else{
			cat(paste(c(filepref,rep("-",times=23)),collapse="\t"),file=txtdir$sumpath.hcomp,eol="\n",append=TRUE);
			cat(paste(c(filepref,rep("-",times=23)),collapse="\t"),file=txtdir$sumpath.scomp,eol="\n",append=TRUE);
			cat(paste(c(filepref,rep("-",times=23)),collapse="\t"),file=txtdir$sumpath.tcomp,eol="\n",append=TRUE);
			cat(paste(c(filepref,rep("-",times=23)),collapse="\t"),file=txtdir$sumpath.ccomp,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(filepref,rep("-",times=23)),nrow=1),htmldir$sumpath.hcomp,TRUE,"<td>");
			write2Html.tableBody(matrix(c(filepref,rep("-",times=23)),nrow=1),htmldir$sumpath.scomp,TRUE,"<td>");
			write2Html.tableBody(matrix(c(filepref,rep("-",times=23)),nrow=1),htmldir$sumpath.tcomp,TRUE,"<td>");
			write2Html.tableBody(matrix(c(filepref,rep("-",times=23)),nrow=1),htmldir$sumpath.ccomp,TRUE,"<td>");
			}
		cat("   Finished ss distribution analysis...\n",file=logfile,sep="",eol="",append=TRUE);	
		}
		
		# Obtain global interaction list object
		intList=calc_interaction(pdb,DSSP,nacc,status.DSSP,status.nacc,calc_ionic=calc_ionic,calc_aroaro=calc_aroaro,calc_aros=calc_aros,calc_catpi=calc_catpi,calc_hphob=calc_hphob,calc_hbond=calc_hbond,calc_disul=calc_disul,ioncut,aroarolow,aroarohigh,aroarodihed,usedihed,aroscut,catpicut,hpcut);
		
		if(calc_disul)
		{
			print("Disulfide");
			hlink.name=paste(hlink.disul,filepref,"_disulfide.html",sep="");
			hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
	
			if(intList$status.dis)
			{
			dissum=disulfide_summary(intList$dis,acclow);
			cat(paste(c(filepref,dissum),collapse="\t"),file=txtdir$sumpath.disul,eol="\n",append=TRUE); # write disulfide summary
			write2Html.tableBody(matrix(c(hlink,dissum),nrow=1),htmldir$sumpath.disul,TRUE,"<td>");
			sumvec[13]=dissum[1];	
      # write disulfide detail
				dis.dirname=paste(htmldir$dirpath.disul,"/Details",sep="",collapse="");
				createDir(dis.dirname);
				dis.filename=paste(dis.dirname,"/",filepref,"_disulfide.html",sep="",collapse="");
				print_HTMLdis(intList$dis,dis.filename);
				
				dis.dirname=paste(txtdir$dirpath.disul,"/Details",sep="",collapse="");
				createDir(dis.dirname);
				dis.filename=paste(dis.dirname,"/",filepref,"_disulfide.txt",sep="",collapse="");
				cat("Disulfide bonds\n-------------------------------------\nChain\tRes1_No\tChain\tRes2_No\tAcc_Res1\tAcc_Res2\tSS_Res1\tSS_Res2\n",file=dis.filename); # Create individual disulfide bond detail file
				write.table(intList$dis,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=dis.filename);
				
			}else{
			cat(paste(c(filepref,rep("-",times=14)),collapse="\t"),file=txtdir$sumpath.disul,eol="\n",append=TRUE);  # No disulfide bond present
			write2Html.tableBody(matrix(c(filepref,rep("-",times=14)),nrow=1),htmldir$sumpath.disul,TRUE,"<td>");
			}
		cat("   Finished decting disulfide bond ...\n",file=logfile,sep="",eol="",append=TRUE);	
		}
		
		if(calc_ionic)
		{
			print("Running ionic");
			hlink.name=paste(hlink.ionic,filepref,"_ionpair.html",sep="");
			hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
				
			if(intList$status.ionic)
			{
				ionsum=ionicSummary(intList$ionpair,acclow);
				cat(paste(c(filepref,ionsum),collapse="\t"),file=txtdir$sumpath.ionic,eol="\n",append=TRUE); # write ip summary
				write2Html.tableBody(matrix(c(hlink,ionsum),nrow=1),htmldir$sumpath.ionic,TRUE,"<td>");
				sumvec[9]=ionsum[1];
				# ip detail file				
				ionic.dirname=paste(htmldir$dirpath.ionic,"/Details",sep="",collapse="");
				createDir(ionic.dirname);
				ionic.filename=paste(ionic.dirname,"/",filepref,"_ionpair.html",sep="",collapse="");
				ionnetmat=get_ionnetmat(intList$ionpair,intList$dionpair,acclow);
				print_HTMLionic(intList$ionpair,ionnetmat,ionic.filename)
				
				ionic.dirname=paste(txtdir$dirpath.ionic,"/Details",sep="",collapse="");
				createDir(ionic.dirname);
				ionic.filename=paste(ionic.dirname,"/",filepref,"_ionpair.txt",sep="",collapse="");
				# write ip detail
				cat("Ionic Interactions\n------------------------\nChain\tRes1_No\tRes1_ID\tChain\tRes2_No\tRes2_ID\tAcc_Res1\tAcc_Res2\tSS_Res1\tSS_Res2\n",file=ionic.filename); # Create individual ionpair table file
				write.table(intList$ionpair,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=ionic.filename);
				# write ip net details
				cat("\n\nNetwork\n-------------------------------------\n",file=ionic.filename,sep="",eol="",append=TRUE)
				print_ionnetwork(intList$ionpair,intList$dionpair,ionic.filename,acclow); #Note: If no network then print_ionnetwork fn does not print network details :)
			}else{
			cat(paste(c(filepref,rep("-",times=19)),collapse="\t"),file=txtdir$sumpath.ionic,eol="\n",append=TRUE); # No ionpairs; In case of no acidic/basic residues or not nearer.
			write2Html.tableBody(matrix(c(filepref,rep("-",times=19)),nrow=1),htmldir$sumpath.ionic,TRUE,"<td>");
			}
		cat("   Finished detecting ionic interactions...\n",file=logfile,sep="",eol="",append=TRUE);		
		}
		
		if(calc_aroaro)
		{
				hlink.name=paste(hlink.aroaro,filepref,"_aroaropair.html",sep="");
				hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
				
			if(intList$status.aroaro)
			{
				print("Running aroaro");
				aroarosum=aroaroSummary(intList$aro,acclow);
				cat(paste(c(filepref,aroarosum),collapse="\t"),file=txtdir$sumpath.aroaro,eol="\n",append=TRUE); # write aroaro summary
				
				write2Html.tableBody(matrix(c(hlink,aroarosum),nrow=1),htmldir$sumpath.aroaro,TRUE,"<td>");
				sumvec[10]=aroarosum[1];
				aroaro.dirname=paste(htmldir$dirpath.aroaro,"/Details",sep="",collapse="");
				createDir(aroaro.dirname);
				aroaro.filename=paste(aroaro.dirname,"/",filepref,"_aroaropair.html",sep="",collapse="");
				aronetmat=get_aroaronetmat(intList$aro,intList$daropair,acclow);
				print_HTMLaroaro(intList$aro,aronetmat,aroaro.filename);
				
				aroaro.dirname=paste(txtdir$dirpath.aroaro,"/Details",sep="",collapse="");
				createDir(aroaro.dirname);
				aroaro.filename=paste(aroaro.dirname,"/",filepref,"_aroaropair.txt",sep="",collapse="");
				
				cat("Aromatic aromatic interactions\n----------------------------------------------\nChain\tRes1_No\tRes1_ID\tChain\tRes2_No\tRes2_ID\tCentroid_Distance(Ang)\tDihedral\tAcc_Res1\tAcc_Res2\tSS_Res1\tSS_Res2\n",file=aroaro.filename); # Create individual ionpair table file
				write.table(intList$aro,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=aroaro.filename);
				
				# Note: If no network then print_ionnetwork fn does not print network details :)
				cat("\n\nNetwork\n-------------------------------------\n",file=aroaro.filename,sep="",eol="",append=TRUE)
				print_aroaronetwork(intList$aro,intList$daropair,aroaro.filename,acclow);
			
			}else{
			# No aroaro pair then
			cat(paste(c(filepref,rep("-",times=17)),collapse="\t"),file=txtdir$sumpath.aroaro,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(filepref,rep("-",times=17)),nrow=1),htmldir$sumpath.aroaro,TRUE,"<td>");
			}
		cat("   Finished detecting aro-aro interactions...\n",file=logfile,sep="",eol="",append=TRUE);	
		}
		if(calc_hbond)
		{
		print("Running hb");
		hlink.name=paste(hlink.hbond,filepref,"_hbond.html",sep="");
		hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
			if(intList$status.hb)
			{
			nncnhb=HBclassify(intList$hb);
			hbsum=HBSummary(intList$hb,nncnhb);
			cat(paste(c(filepref,hbsum),collapse="\t"),file=txtdir$sumpath.hbond,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(hlink,hbsum),nrow=1),htmldir$sumpath.hbond,TRUE,"<td>");
			sumvec[14]=hbsum[1];
			hbond.dirname=paste(htmldir$dirpath.hbond,"/Details",sep="",collapse="");
			createDir(hbond.dirname);
			hbond.filename=paste(hbond.dirname,"/",filepref,"_hbond.html",sep="",collapse="");
			if(nrow(intList$hb)){print_HTMLhb(intList$hb,hbond.filename,append=FALSE,html="FALSE","<br><br><b><u>Hydrogen bonding interaction</u><br><br></b>\n");print("*");}
			if(nrow(nncnhb$cnhbond)!=0){print_HTMLhb(nncnhb$cnhbond,hbond.filename,append=TRUE,html="FALSE","<br><br><b><u>Charged neutral hydrogen bonds </u><br><br></b>\n");print("**");}
			if(nrow(nncnhb$nnhbond)!=0){print_HTMLhb(nncnhb$nnhbond,hbond.filename,append="TRUE",html=TRUE,"<br><br><b><u>Neutral-neutral Hydrogen bonds</u><br><br></b>\n");print("***");}
			
			
			hbond.dirname=paste(txtdir$dirpath.hbond,"/Details",sep="",collapse="");
			createDir(hbond.dirname);
			hbond.filename=paste(hbond.dirname,"/",filepref,"_hbond.txt",sep="",collapse="");
			cat("Hydrogen bonds\n-------------------------------------\nD.Chain\tD.Resno\tD.Resid\tD.atom\tA.Chain\tA.Resno\tA.Resid\tA.atom\tBond_type\tDistance_DA\tDistance_HA\tDHA_Angle\tHAAA_Angle\tDAAA_Angle\n",file=hbond.filename); # Create individual ionpair table file
			write.table(intList$hb,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=hbond.filename);
			
			
			if(nrow(nncnhb$cnhbond)!=0)
			{
			cat("\n\nCharged-neutal hydrogen bonds\n\nD.Chain\tD.Resno\tD.Resid\tD.atom\tA.Chain\tA.Resno\tA.Resid\tA.atom\tBond_type\tDistance_DA\tDistance_HA\tDHA_Angle\tHAAA_Angle\tDAAA_Angle\n",file=hbond.filename,append=TRUE); # Create individual ionpair table file
			write.table(nncnhb$cnhbond,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=hbond.filename);
			}
			
			if(nrow(nncnhb$nnhbond)!=0)
			{
			cat("\n\nNeutral-neutal hydrogen bonds\n\nD.Chain\tD.Resno\tD.Resid\tD.atom\tA.Chain\tA.Resno\tA.Resid\tA.atom\tBond_type\tDistance_DA\tDistance_HA\tDHA_Angle\tHAAA_Angle\tDAAA_Angle\n",file=hbond.filename,append=TRUE); # Create individual ionpair table file
			write.table(nncnhb$nnhbond,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=hbond.filename);
			}
			
			}else{# No hbond cases
			cat(paste(c(filepref,rep("-",times=14)),collapse="\t"),file=txtdir$sumpath.hbond,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(filepref,rep("-",times=14)),nrow=1),htmldir$sumpath.hbond,TRUE,"<td>");
			}
		cat("   Finished detecting hbond interactions...\n",file=logfile,sep="",eol="",append=TRUE);	
		}
		
		if(calc_aros)
		{
			print("Running aros");
			hlink.name=paste(hlink.aros,filepref,"_aroS.html",sep="");
			hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
			
			#print(intList$status.aros);
			if(intList$status.aros)
			{
				
				arossum=arosSummary(intList$aros,acclow);
				cat(paste(c(filepref,arossum),collapse="\t"),file=txtdir$sumpath.aros,eol="\n",append=TRUE);
				
				write2Html.tableBody(matrix(c(hlink,arossum),nrow=1),htmldir$sumpath.aros,TRUE,"<td>");
				sumvec[11]=arossum[1];
				aros.dirname=paste(htmldir$dirpath.aros,"/Details",sep="",collapse="");
				createDir(aros.dirname);
				aros.filename=paste(aros.dirname,"/",filepref,"_aroS.html",sep="",collapse="");
				#print(aros.filename);
				arosnetmat=get_arosnetmat(intList$aros,intList$darospair,acclow);
				print_HTMLaroS(intList$aros,arosnetmat,aros.filename);
				
				aros.dirname=paste(txtdir$dirpath.aros,"/Details",sep="",collapse="");
				createDir(aros.dirname);
				aros.filename=paste(aros.dirname,"/",filepref,"_aroS.txt",sep="",collapse="");
				
				cat("Aromatic-Sulphur interactions\n------------------------------------------\nChain\tRes1_No\tRes1_ID\tChain\tRes2_No\tRes2_ID\tDistance(Ang)\tAcc_Res1\tAcc_Res2\tSS_Res1\tSS_Res2\n",file=aros.filename); # Create individual ionpair table file
				write.table(intList$aros,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=aros.filename);
			
				cat("\n\nNetwork\n-------------------------------------\n",file=aros.filename,sep="",eol="",append=TRUE)
				print_arosnetwork(intList$aros,intList$darospair,aros.filename,acclow);
			}else{
			# No AroS then
			cat(paste(c(filepref,rep("-",times=19)),collapse="\t"),file=txtdir$sumpath.aros,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(filepref,rep("-",times=19)),nrow=1),htmldir$sumpath.aros,TRUE,"<td>");
			}
		cat("   Finished detecting Aro-S interactions...\n",file=logfile,sep="",eol="",append=TRUE);	
		}
		
		if(calc_asa)
		{
		print("Running ASA");
			if(status.nacc)
			{
			
			hlink.name=paste(hlink.asa,filepref,"_ASAreswise.html",sep="");
			hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
			
			rsa.name=paste(hlink.asa,naccprefix,".rsa",sep="");
			rsalink=paste("<a href=\"",rsa.name,"\" target=\"_blank\" >","RSA","</a>");
			
			asasum=asa_summary(nacc);
			naccres=asa_perres(nacc);
			
			asa.dirname=paste(htmldir$dirpath.asa,"/Details",sep="",collapse="");
			createDir(asa.dirname);
			asa.dirname1=paste(txtdir$dirpath.asa,"/Details",sep="",collapse="");
			createDir(asa.dirname1);
			
				if(file.exists(paste(naccprefix,".rsa",sep="")))
				{
				file.copy(from=paste(naccprefix,".rsa",sep=""),to=asa.dirname);
				file.copy(from=paste(naccprefix,".rsa",sep=""),to=asa.dirname1);
				}
			
			asa.filename=paste(asa.dirname,"/",filepref,"_ASAreswise.html",sep="",collapse="");
			print_HTMLasa(naccres,asa.filename);
			
			asa.filename=paste(asa.dirname1,"/",filepref,"_ASAreswise.txt",sep="",collapse="");
			cat("Per-residue ASA profile\n-----------------------------------\nAA\tTotal\tB\tE\n",file=asa.filename); # Create individual ionpair table file
			write.table(naccres,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=asa.filename);
			
			cat(paste(c(filepref,asasum),collapse="\t"),file=txtdir$sumpath.asa,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(hlink,asasum,rsalink),nrow=1),htmldir$sumpath.asa,TRUE,"<td>");
			sumvec[19]=asasum[6];
			
			}else{
			cat(paste(c(filepref,rep("-",times=10)),collapse="\t"),file=txtdir$sumpath.asa,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(filepref,rep("-",times=11)),nrow=1),htmldir$sumpath.asa,TRUE,"<td>"); #why 11 bcoz one extra column for rsa file link
			}
		cat("   Finished ASA calculations...\n",file=logfile,sep="",eol="",append=TRUE);	
		}
		
		
		if(calc_catpi)
		{
				hlink.name=paste(hlink.catpi,filepref,"_catPi.html",sep="");
				hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
			if(intList$status.catpi)
			{
				catpisum=catpiSummary(intList$catpimat,acclow);
				cat(paste(c(filepref,catpisum),collapse="\t"),file=txtdir$sumpath.catpi,eol="\n",append=TRUE);
				write2Html.tableBody(matrix(c(hlink,catpisum),nrow=1),htmldir$sumpath.catpi,TRUE,"<td>");
				sumvec[12]=catpisum[1];
				catpi.dirname=paste(htmldir$dirpath.catpi,"/Details",sep="",collapse="");
				createDir(catpi.dirname);
				catpi.filename=paste(catpi.dirname,"/",filepref,"_catPi.html",sep="",collapse="");
				catpinetmat=get_catpinetmat(intList$catpimat,intList$dcatpimat,acclow);
				print_HTMLcatpi(intList$catpimat,catpinetmat,catpi.filename);
				
				catpi.dirname=paste(txtdir$dirpath.catpi,"/Details",sep="",collapse="");
				createDir(catpi.dirname);
				catpi.filename=paste(catpi.dirname,"/",filepref,"_catPi.txt",sep="",collapse="");
				cat("Cation-pi interactions\n---------------------------------\nChain\tRes1_No\tRes1_ID\tChain\tRes2_No\tRes2_ID\tAcc_Res1\tAcc_Res2\tSS_Res1\tSS_Res2\n",file=catpi.filename); # Create individual ionpair table file
				write.table(intList$catpimat,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=catpi.filename);
	
				# Note: If no network then print_ionnetwork fn does not print network details :)
				cat("\n\nNetwork\n-------------------------------------\n",file=catpi.filename,sep="",eol="",append=TRUE)
				print_catpinetwork(intList$catpimat,intList$dcatpimat,catpi.filename,acclow);
				
				
			}else{			# No aroaro pair then
			cat(paste(c(filepref,rep("-",times=14)),sep="",collapse="\t"),file=txtdir$sumpath.catpi,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(filepref,rep("-",times=14)),nrow=1),htmldir$sumpath.catpi,TRUE,"<td>");
			}
		cat("   Finished detecting Cation-pi interactions...\n",file=logfile,sep="",eol="",append=TRUE);	
		}
		
		
		
		if(calc_hphob)
		{
				hlink.name=paste(hlink.hphob,filepref,"_Hydrophobic.html",sep="");
				hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
			if(intList$status.hp)
			{
				hpsum=hphob_summary(intList$hp,acclow);
				cat(paste(c(filepref,hpsum),collapse="\t"),file=txtdir$sumpath.hphob,eol="\n",append=TRUE);
				write2Html.tableBody(matrix(c(hlink,hpsum),nrow=1),htmldir$sumpath.hphob,TRUE,"<td>");
				sumvec[15]=hpsum[1];
				hphob.dirname=paste(htmldir$dirpath.hphob,"/Details",sep="",collapse="");
				createDir(hphob.dirname);
				hphob.filename=paste(hphob.dirname,"/",filepref,"_Hydrophobic.html",sep="",collapse="");
				print_HTMLhphob(intList$hp,hphob.filename);
				
				hphob.dirname=paste(txtdir$dirpath.hphob,"/Details",sep="",collapse="");
				createDir(hphob.dirname);
				hphob.filename=paste(hphob.dirname,"/",filepref,"_Hydrophobic.txt",sep="",collapse="");
				cat("Details of Hydrophobic interaction\n--------------------------------\n\nChain\tRes1_No\tRes1_ID\tChain\tRes2_No\tRes2_ID\tDistance\tSS_Res1\tSS_Res2\tAcc_Res1\tAcc_Res2\n",file=hphob.filename); # Create individual ionpair table file
				write.table(intList$hp,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=hphob.filename);
				
			}else{
			# No aroaro pair then
			cat(paste(c(filepref,rep("-",times=5)),collapse="\t"),file=txtdir$sumpath.hphob,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(filepref,rep("-",times=5)),nrow=1),htmldir$sumpath.hphob,TRUE,"<td>"); 
			}
		cat("   Finished detecting Hphob interactions...\n",file=logfile,sep="",eol="",append=TRUE);	
		}
		
		if(calc_proline)
		{
			
			hlink.name=paste(hlink.proline,filepref,"_proline.html",sep="");
			hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
			#print("Running proline profile");
			pro=proline(pdb,DSSP,nacc,acclow=acclow,prefix=promotifprefix);
			cat(paste(c(filepref,pro$prop),collapse="\t"),file=txtdir$sumpath.proline,eol="\n",append=TRUE);
			
			
			if(pro$prop[1] > 0)
			{
			  sumvec[16]=pro$prop[1];
			  sumvec[17]=pro$prop[9];
			  sumvec[18]=pro$prop[10];
			  
			# If either BT2P or Ncap proline present then create individual file
			write2Html.tableBody(matrix(c(hlink,pro$prop),nrow=1),htmldir$sumpath.proline,TRUE,"<td>");
			proline.dirname=paste(htmldir$dirpath.proline,"/Details",sep="",collapse="");
			createDir(proline.dirname);
			proline.filename=paste(proline.dirname,"/",filepref,"_proline.html",sep="",collapse="");
			print_HTMLproline(pro,proline.filename);
			
			proline.dirname=paste(txtdir$dirpath.proline,"/Details",sep="",collapse="");
			createDir(proline.dirname);
			proline.filename=paste(proline.dirname,"/",filepref,"_proline.txt",sep="",collapse="");
			
			cat("Proline residue distributon in each type of secondary structure\n----------------------------------------------------\nSecondary structure type\tResidue details[Format: Chain Resno]\n",file=proline.filename); # Create individual ionpair table file
			write.table(pro$ssdet,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=proline.filename);
				
				if(pro$prop[9]>0){
				cat("\nDetails of Beta Turn 2nd position proline residue\n----------------------------------------------------\nChain\tStart\tEnd\tSequence\tBT_Type\tR2_phi\tR2_psi\tR2_chi1\tR3_phi\tR3_psi\tR3_chi1\tR1R4_CAdist\tHbond\tSS\tB/E\n",file=proline.filename,append=TRUE); # Create individual ionpair table file
				write.table(pro$Bt2P,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=proline.filename);
				}
				
				if(pro$prop[10]>0){
				cat("\nDetails of N-cap proline residue\n--------------------------------\nChain\tStart\tResid\tType\tAcc\n",file=proline.filename,append=TRUE); # Create individual ionpair table file
				write.table(pro$NcapP,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=proline.filename);
				}
			}else{
			write2Html.tableBody(matrix(c(filepref,pro$prop),nrow=1),htmldir$sumpath.proline,TRUE,"<td>");
			}
		cat("   Finished proline residue profile ...\n",file=logfile,sep="",eol="",append=TRUE);	
		}
		if(calc_hdipo)
		{
			#print("Running dipole profile");
			hlink.name=paste(hlink.hdipo,filepref,"_Hdipole.html",sep="");
			hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
			if(status.DSSP)
			{
			hdip=helixDipole(DSSP);
			hdipsum=c(hdip$tot_helices,hdip$tot_StabHelices,hdip$N_cap_StabHelices,hdip$C_cap_StabHelices,hdip$NC_cap_StabHelices,hdip$hstat[2,]);
			cat(paste(c(filepref,hdipsum),collapse="\t"),file=txtdir$sumpath.hdipo,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(hlink,hdipsum),nrow=1),htmldir$sumpath.hdipo,TRUE,"<td>");
			sumvec[21]=hdipsum[2];
			if(is.matrix(hdip$hdipo))
			{
			hdipo.dirname=paste(htmldir$dirpath.hdipo,"/Details",sep="",collapse="");
			createDir(hdipo.dirname);
			hdipo.filename=paste(hdipo.dirname,"/",filepref,"_Hdipole.html",sep="",collapse="");
			print_HTMLhelixDipole(hdip$hdipo,hdipo.filename);
			
			hdipo.dirname=paste(txtdir$dirpath.hdipo,"/Details",sep="",collapse="");
			createDir(hdipo.dirname);
			hdipo.filename=paste(hdipo.dirname,"/",filepref,"_Hdipole.txt",sep="",collapse="");
			cat("Details of Helix dipole stabilized helices\n--------------------------------\nChain\tH.start\tH.end\tH.type\tLength\tNcapAcid\tCcapBase\n",file=hdipo.filename); # Create individual ionpair table file
			write.table(hdip$hdipo,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=hdipo.filename);
			
			}
			}else{
			cat(paste(c(filepref,rep("-",times=15)),collapse="\t"),file=txtdir$sumpath.hdipo,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(filepref,rep("-",times=15)),nrow=1),htmldir$sumpath.hdipo,TRUE,"<td>");
			}
		cat("   Finished helix dipole stabilization profile...\n",file=logfile,sep="",eol="",append=TRUE);	
		}
		
		if(calc_deam)
		{
		#print("Running deamidation profile");
		deamid=Deamidation_protein(pdb); # we get either a matrix else SORRY
		deamidsum="";
		hlink.name=paste(hlink.deam,filepref,"_deamidation.html",sep="");
		hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
				
			if(is.matrix(deamid))
			{
			if(status.DSSP){deamid=Deamidation_addDSSP(DSSP,deamid);}else{deamid=cbind(deamid,SS="-");}
			if(status.nacc){deamid=Deamidation_addNaccess(nacc,deamid);}else{deamid=cbind(deamid,acc="-");}
			deamidsum=Deamidation_summary(deamid);
			cat(paste(c(filepref,deamidsum),collapse="\t"),file=txtdir$sumpath.deam,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(hlink,deamidsum),nrow=1),htmldir$sumpath.deam,TRUE,"<td>");
			sumvec[20]=deamidsum;
			deam.dirname=paste(htmldir$dirpath.deam,"/Details",sep="",collapse="");
			createDir(deam.dirname);
			deam.filename=paste(deam.dirname,"/",filepref,"_deamidation.html",sep="",collapse="");
			print_HTMLdeamid(deamid,deam.filename);
			
			deam.dirname=paste(txtdir$dirpath.deam,"/Details",sep="",collapse="");
			createDir(deam.dirname);
			deam.filename=paste(deam.dirname,"/",filepref,"_deamidation.txt",sep="",collapse="");
			cat("Details of Thermolabile bonds\n--------------------------------\nChain\tPosition\tRes.Pair\tDistance\tSS\tAcc\n",file=deam.filename); # Create individual ionpair table file
			write.table(deamid,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=deam.filename);
			
			}else{
			cat(paste(c(filepref,rep("-",times=8)),collapse="\t"),file=txtdir$sumpath.deam,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(filepref,rep("-",times=8)),nrow=1),htmldir$sumpath.deam,TRUE,"<td>");
			}
		cat("   Finished thermolabile bond profile...\n",file=logfile,sep="",eol="",append=TRUE);	
		}
		
		if(calc_left)
		{
			#print("Detecting Conformationally strained residue profile");
			left=procheck(pdb,DSSP,procheckprefix); # returns a matrix or SORRY
			
			hlink.name=paste(hlink.left,filepref,"_ConfStrain.html",sep="");
			hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
			
			leftsum="";
			if(is.matrix(left))
			{
			leftsum=procheck_summary(left);
			cat(paste(c(filepref,leftsum),collapse="\t"),file=txtdir$sumpath.left,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(hlink,leftsum),nrow=1),htmldir$sumpath.left,TRUE,"<td>");
			sumvec[22]=leftsum[1];
			left.dirname=paste(htmldir$dirpath.left,"/Details",sep="",collapse="");
			createDir(left.dirname);
			left.filename=paste(left.dirname,"/",filepref,"_ConfStrain.html",sep="",collapse="");
			print_HTMLprocheck(left,left.filename);
			
			left.dirname=paste(txtdir$dirpath.left,"/Details",sep="",collapse="");
			createDir(left.dirname);
			left.filename=paste(left.dirname,"/",filepref,"_ConfStrain.txt",sep="",collapse="");
			cat("Details of conformationally strained residues\n--------------------------------\nChain\tRes.No\tRes.ID\tSS\tR.Plot_region\tPhi\tPsi\tDistance\n",file=left.filename); # Create individual ionpair table file
			write.table(left,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=left.filename);
			
			}else{
			cat(paste(c(filepref,rep("-",times=5)),collapse="\t"),file=txtdir$sumpath.left,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(filepref,rep("-",times=5)),nrow=1),htmldir$sumpath.left,TRUE,"<td>");
			}
		cat("   Finished conf. strained residue profile...\n",file=logfile,sep="",eol="",append=TRUE);	
		}
		
		if(calc_metal)
		{
			print("Running metal");
			metaldir=paste(work_dir,"/",prefix,"_metal",count,sep="",collapse="");	# /var/www/workdir/prefix/prefix_metal1, /var/www/workdir/prefix/prefix_metal2...
			createDir(metaldir);
			setwd(metaldir); #change to individual metal directory
			#print(getwd());
      print(pdbfile);
			#file.copy(from=paste(work_dir,"/",pdbfile,sep="",collapse=""),to=".");
			file.copy(from=pdbfile,to=".");
			print(system("ls"));
			fgeoout=paste(metaldir,"/fgeoout",count,".txt",sep="",collapse=""); #fgeo output increments as per each pdb count variable
			print(fgeoout);
      metalcmd=paste("perl ",fgeoperl," ",fgeopython," ",fgeopath," ",f," ",fgeoout,sep="",collapse="");
      print(metalcmd);
			system(metalcmd,ignore.stdout=FALSE,ignore.stderr=FALSE);
			print(file.exists(fgeoout));
			#print(readLines(fgeoout));
			if(file.exists(fgeoout))
			{
			fgeodata=read.table(fgeoout,sep="\t",na.strings="~",stringsAsFactors=FALSE);
			if(nrow(fgeodata)==0){fgeodata=matrix(apply(matrix(readLines(fgeoout),ncol=1),1,function(x){unlist(strsplit(x,"\\t"))}),ncol=8,byrow=TRUE)}
			#print(fgeodata);

			metals=tapply(fgeodata[,2],factor(fgeodata[,2]),length);
			metalsum=c(nrow(fgeodata),paste(names(metals),metals,sep=":",collapse=", "));

			fgeofinal=fgeodata[,c(2:5,7:8)];
			fgeocode=fgeodata[,6];
			
			#write summary
			hlink.name=paste(hlink.metal,filepref,"_metal.html",sep="");
			hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
			cat(paste(c(filepref,metalsum),collapse="\t"),file=txtdir$sumpath.metal,eol="\n",append=TRUE); # write metal summary
			write2Html.tableBody(matrix(c(hlink,metalsum),nrow=1),htmldir$sumpath.metal,TRUE,"<td>");
			sumvec[23]=metalsum[1];
			#detailed file
			
			metal.dirname=paste(htmldir$dirpath.metal,"/Details",sep="",collapse="");
			createDir(metal.dirname);
			metal.filename=paste(metal.dirname,"/",filepref,"_metal.html",sep="",collapse="");
			print_HTMLmetal(fgeofinal,fgeocode,metal.filename);
			
			metal.dirname=paste(txtdir$dirpath.metal,"/Details",sep="",collapse="");
			createDir(metal.dirname);
			metal.filename=paste(metal.dirname,"/",filepref,"_metal.txt",sep="",collapse="");
			
			cat("Metal geometry and binding site details\n--------------------------------\nMetal\tChain\tRes.No\tAtom.No\tGeometry\tSite_details\n",file=metal.filename); # Create individual ionpair table file
			write.table(fgeofinal,sep="\t",eol="\n",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,file=metal.filename);
					
			# change to original work dir	
			}else{
			# No metal binding case
			cat(paste(c(filepref,rep("-",times=2)),collapse="\t"),file=txtdir$sumpath.metal,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(filepref,rep("-",times=2)),nrow=1),htmldir$sumpath.metal,TRUE,"<td>");
			}
			setwd(work_dir); 
		cat("   Finished metal binding analysis...\n",file=logfile,sep="",eol="",append=TRUE);	

		}
		
		if(calc_gibbs)
		{
			print("Estimating Gibbs free energy of folding");

			foldxpdb=paste(filepref,"_foldx",count,".pdb",sep="",collapse=""); #the foldx prefix which is unique to each pdb file
			write.pdb(pdb,file=foldxpdb); #write pdb file
			
			foldxRunfile=paste(filepref,"_foldx",count,"_FoldxRunfile",sep="",collapse=""); #using the foldx prefix name the runfile so for each pdb runfile name will be different
			foldxprefix=paste(filepref,"_foldx",count,"_foldxtemp.txt",sep="",collapse=""); #foldx stability output file name, each pdb will have different name
			foldXcreateRunFileStability(foldxpdb,foldxRunfile,foldxprefix); # This will generate the run file

			stabdata=runFoldxStability(pdb,foldxpath,foldxprefix,foldxpdb); # foldxprefix argument is used for parsing the output
			print(stabdata);
			
			if(stabdata[1]!="SORRY")
			{
			stabdata=round(as.numeric(stabdata),2);
      cat(paste(c(filepref,stabdata),sep="",collapse="\t"),file=txtdir$sumpath.gibbs,eol="\n",append=TRUE);
			write2Html.tableBody(matrix(c(filepref,stabdata),nrow=1),htmldir$sumpath.gibbs,TRUE,"<td>");
			sumvec[24]=stabdata[1];
      }else{
			cat(paste(c(filepref,rep("-",times=23)),collapse="\t"),file=txtdir$sumpath.gibbs,eol="\n",append=TRUE); # For those pdb DSSP fails, no secstr calculation
			write2Html.tableBody(matrix(c(filepref,rep("-",times=23)),nrow=1),htmldir$sumpath.gibbs,TRUE,"<td>");
			
			}
		cat("   Finished Gibbs free energy estimation...\n",file=logfile,sep="",eol="",append=TRUE);	
		}

		write2Html.tableBody(matrix(sumvec,nrow=1),sumhtmlfile,TRUE,"<td>");
    write.table(matrix(sumvec,nrow=1),file=sumtxtfile,sep="\t",eol="\n",quote=F,row.names=F,append=T,col.names=F);

#print("****");
filecount=filecount+1;
cat(paste(floor((filecount*100)/totalfile), sep="",collapse=""),file=statusfile,sep="",eol="");

}


# dont confuse with the name of print_HTMLionic_sum_end function..it does same job for all and it is kept under ionic.r
if(calc_ionic){print_HTMLionic_sum_end(htmldir$sumpath.ionic);	}
if(calc_aroaro){print_HTMLionic_sum_end(htmldir$sumpath.aroaro);}
if(calc_aacomp){print_HTMLionic_sum_end(htmldir$sumpath.aacomp);}
if(calc_secstr){print_HTMLionic_sum_end(htmldir$sumpath.secstr);}

if(calc_secstrdist)
{
print_HTMLionic_sum_end(htmldir$sumpath.hcomp);
print_HTMLionic_sum_end(htmldir$sumpath.scomp);
print_HTMLionic_sum_end(htmldir$sumpath.tcomp);
print_HTMLionic_sum_end(htmldir$sumpath.ccomp);
}
if(calc_disul){print_HTMLionic_sum_end(htmldir$sumpath.disul);}
if(calc_hbond){print_HTMLionic_sum_end(htmldir$sumpath.hbond);}
if(calc_aros){print_HTMLionic_sum_end(htmldir$sumpath.aros);}
if(calc_catpi){print_HTMLionic_sum_end(htmldir$sumpath.catpi);}
if(calc_hphob){print_HTMLionic_sum_end(htmldir$sumpath.hphob);}

if(calc_asa){print_HTMLionic_sum_end(htmldir$sumpath.asa);}
if(calc_proline){print_HTMLionic_sum_end(htmldir$sumpath.proline);}
if(calc_hdipo){print_HTMLionic_sum_end(htmldir$sumpath.hdipo);}
if(calc_deam){print_HTMLionic_sum_end(htmldir$sumpath.deam);}
if(calc_left){print_HTMLionic_sum_end(htmldir$sumpath.left);}
if(calc_metal){print_HTMLionic_sum_end(htmldir$sumpath.metal);}
if(calc_gibbs){print_HTMLionic_sum_end(htmldir$sumpath.gibbs);}

print_HTMLionic_sum_end(sumhtmlfile);


# remove the working directory
setwd(htmldir$dirpath);
#unlink(work_dir, recursive = TRUE);
system(paste("zip -r ","./result.zip ", "./",subpref2,sep="",collapse=""),ignore.stderr=FALSE);
cat("100",file=statusfile,sep="",eol="");
},
error=function(e)
{
  print("--ERROR--");
  cat("-111",file=statusfile,sep="",eol="");
}

)



