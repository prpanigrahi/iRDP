############################################
# This is the version for doing html and txt output
############################################

print("USAGE: /usr/bin/Rscript Pattern_analysis_V7.r optionfile");
args=commandArgs(trailingOnly = TRUE);
#print(getwd());
optionfile=args[1]; # Arg1 is option file
print(optionfile);
options=read.table(optionfile,stringsAsFactors=FALSE) # Read the option file
#print(length(options));
prefix=options[which(options$V1=="prefix"),2] # First row of option file is the session id to be used as prefix
print(c("Prefix entered",prefix));
#prefix="temp";

# PDB fetch 
separate=",";

#pdbids=options[which(options$V1=="pdbids"),2]; # second row of option file is pdbids comma separated
#pdbids=gsub("\n","",pdbids); #chomp any \n
#filenames=unique(paste("pdb",unlist(strsplit(pdbids,separate)),".ent",sep=""));
#filenames=tolower(unique(paste("pdb",unlist(strsplit(pdbids,separate)),".ent",sep="")));

pdbids=paste(readLines(args[2]),sep="",collapse="");
pdbids=gsub("\n","",pdbids); #chomp any \n
pdbids=gsub(" ","",pdbids); #chomp any space
filenames=tolower(unique(paste("pdb",unlist(strsplit(pdbids,separate)),".ent",sep="")));

print(filenames);
print(length(filenames));

#if(length(filenames)>100)
#{
#  filenames=filenames[1:100]; #only 100 files will be accepted
#}


pat=toupper(options[which(options$V1=="pattern"),2]);

####################### All directories
script_dir="/var/www/Script_dir/";
work_dir="/var/www/workdir/"; # It has to be shorter coz procheck fails if a longer path comes, donno why.
#pdb_dir="/home/redesign/structures/ent";
pdb_dir="/var/www/ent";
#output_dir="/var/www/html/out/";
output_dir="/data/out/";
hlinkpath="/out/";

# Load all Libraries and sources
loadsource=function(script_dir,file)
{
	source(paste(script_dir,file,sep="",collapse=""));
}

library("bio3d");
library("parallel");
library("igraph");
loadsource(script_dir,"PATH.r");
loadsource(script_dir,"dssp.r");
loadsource(script_dir,"naccess.r");
loadsource(script_dir,"ionic.r");
loadsource(script_dir,"general.r");
loadsource(script_dir,"aroaro.r");
loadsource(script_dir,"aroS.r");
loadsource(script_dir,"catpi.r");
loadsource(script_dir,"procheck.r");
loadsource(script_dir,"disulfide.r");
loadsource(script_dir,"hbond.r");
loadsource(script_dir,"Hphob.r");
loadsource(script_dir,"patfind.r");
loadsource(script_dir,"torpat.r");
loadsource(script_dir,"intprofile.r");
loadsource(script_dir,"writehtml.r");

# Setting sfor which programs to run
catpimode=1;
use_capture=0;
calc_uniqch=as.numeric(options[which(options$V1=="uniq"),2]);
catpicut=as.numeric(options[which(options$V1=="catpicut"),2]);
ioncut=as.numeric(options[which(options$V1=="ioncut"),2]);
aroarolow=as.numeric(options[which(options$V1=="aroarolow"),2]);
aroarohigh=as.numeric(options[which(options$V1=="aroarohigh"),2]);
usedihed=as.numeric(options[which(options$V1=="usedihed"),2]);
aroarodihed=as.numeric(options[which(options$V1=="aroarodihed"),2]);
aroscut=as.numeric(options[which(options$V1=="aroscut"),2]);
hpcut=as.numeric(options[which(options$V1=="hpcut"),2]);
acclow=as.numeric(options[which(options$V1=="acclow"),2]);
print(c("acclow",acclow));
######################## Set the working directory
work_dir=paste(work_dir,prefix,sep="",collapse="");
createDir(work_dir);
setwd(work_dir);

#count=70477
#filenames=list.files(pdb_dir,pattern=".ent$");
#filenames=filenames[3];
#filenames=filenames[which(filenames=="1A5T.pdb")]
#filenames=filenames[count:length(filenames)]

subpref1="HTML";
subpref2="TEXT";

makedirPath=function(output_dir,prefix,subpref,suffix)
{
dirpath=paste(output_dir,prefix,sep="",collapse="");
subdirpath=paste(dirpath,"/",subpref,sep="",collapse="");
subsubdirpath=paste(subdirpath,"/Details",sep="",collapse="");
sumpath=paste(subdirpath,"/Pattern_summary.",suffix,sep="",collapse="");
sumpath.int="";
return(list(dirpath=dirpath,subdirpath=subdirpath,subsubdirpath=subsubdirpath,sumpath=sumpath,sumpath.int=sumpath.int));
}

htmldir=makedirPath(output_dir,prefix,subpref1,"html");
txtdir=makedirPath(output_dir,prefix,subpref2,"txt");

# create root directory
createDir(htmldir$dirpath); #/100
createDir(htmldir$subdirpath); # /100/HTML
createDir(htmldir$subsubdirpath); # /100/HTML/Details
createDir(txtdir$subdirpath); # /100/TXT
createDir(txtdir$subsubdirpath); # /100/TXT/Details

cat("Filename\tChain\tStart\tPattern\tAcc\tSS\tR.Plot_region\tIP\tIP_Net\tAP\tAP_Net\tAS\tAS_Net\tHB\tDisul\tCat-pi\tCat-pi_Net\tHphob",file=txtdir$sumpath);
print_HTMLionic_torpat_start(htmldir$sumpath);

#Status file
statusfile=paste(htmldir$subdirpath,"/status.txt",sep="",collapse="");
logfile=paste(htmldir$subdirpath,"/log.txt",sep="",collapse="");
cat(paste(c("Total number of files to be processed:",length(filenames),"\n"),sep="",collapse=""),file=logfile,sep="",eol="");


#########################################################################################################################################################
# PROGRAM STARTS FROM HERE
#########################################################################################################################################################
	for(f in filenames)
	{
	pdb="";nacc="";DSSP="";proch="";
	filepref=substr(f,4,7);
	#filepref=substr(f,1,4);
	
	# All temp pdbs names
	count=which(filenames==f);
	#print(c("Count: ",count));
	#count=count+1;
	dsspprefix=paste(prefix,"_dssp",count,sep="",collapse="");
	naccprefix=paste(prefix,"_nacc",count,sep="",collapse="");
	hbplusprefix=paste(prefix,"_hbplus",count,sep="",collapse="");
	promotifprefix=paste(prefix,"_promotif",count,sep="",collapse="");
	procheckprefix=paste(prefix,"_procheck",count,sep="",collapse="");
	captureprefix=paste(prefix,"_capture",count,sep="",collapse="");

	print(paste("Processing ProtTS for file ",f,sep="",collapse=""));
	cat(paste(c("Processing file: ",f,"\n"),sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
	
	pdbfile=paste(pdb_dir,"/",f,sep="");
	
	
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
			cat("Downloaded file from PDB",file=logfile,sep="",eol="",append=TRUE);
			cat(pdblines,sep="\n",eol="\n",file=f);
			pdbfile=f;
			}else{
			print(c("PDB file not found",pdbfile));
			cat(paste("\n!!! PDB FILE NOT FOUND !!! ",f,"\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE);
			next; 
		}
	}
		
	pdb=read.pdb(pdbfile,maxlines=1000000,verbose=FALSE); # read the pdb file
	now <- Sys.time()
		print(now);
		
	#Do uniquechain calculation	
		uniqch=uniqueChain(pdb);
		if(calc_uniqch)
		{
		cat("Doing unique chain calculations\n\n",file=logfile,sep="",eol="",append=TRUE);  
		if(!(sum(is.na(uniqch))>0 | sum(uniqch==" ")>0)){pdb=trim.pdb(pdb,inds=atom.select(pdb,chain=uniqch,verbose=FALSE));}
		}
		print(uniqch);
print("*****");
	# Discard the PDB if no pattern found
		flag=0;
		for(ch in uniqch)
		{
		patlist=patFind(pdb,ch,pat);
		print(patlist);
			if(! is.na(patlist[1]))
			{
				flag=1;
			}	
		}
		if(flag==0)
		{
		
		cat("No Pattern Found, Skipping\n\n",file=logfile,sep="",eol="",append=TRUE);  
		print("No pattern skiiping by next"); next;
		
		}	
	
		print(unique(pdb$atom[,"chain"]))	
		nacc=naccess(pdb,naccpath,prefix=naccprefix); # Runs naccess
	
	DSSP=dssp_new(pdb,exepath=dssppath,prefix=dsspprefix); #Runs DSSP
	proch=get_procheck(pdb,procheckprefix);
	#print(c("Procheck",proch));
	status.DSSP=0;if(is.list(DSSP)){status.DSSP=1;}
	status.nacc=0;if(is.list(nacc)){status.nacc=1;}
	status.procheck=0;if(is.matrix(proch)){status.procheck=1;}
	#print("starting tors calculation");
	tors=newtorsion.pdb(pdb); # finds torsion
	tors$tbl[!is.na(tors$tbl)]=round(tors$tbl[!is.na(tors$tbl)],2);
	#print("starting intList calculation");
	intList=calc_interaction(pdb,DSSP,nacc,status.DSSP,status.nacc,calc_ionic=0,calc_aroaro=0,calc_aros=0,calc_catpi=0,calc_hphob=0,calc_hbond=1,calc_disul=0,ioncut,aroarolow,aroarohigh,aroarodihed,usedihed,aroscut,catpicut,hpcut); # calculate interaction list


			
		#Loop over each chain to find the pattern
		# We go chainwise.
		# For each chain, call to patFind() will return pdb resno of start residue of pattern.
		# If no pattern found it returns NA, if more than one found then it returns a vector
		# So we loop over each pattern position and extract its details
		# First we need to fetch the torsion information.
		# Inside the Loop, for each pattern position in a given chain we get "chain","pos","pattern"
		# Then we call torpat which returns a vector of all torsion detail, 8 per residue, phi,psi,omega, chi1 to chi5
	for(ch in uniqch)
		{
			print(ch);
			#print(pat);
			patlist=patFind(pdb,ch,pat) # It returns actual residue no information in pdb file, if not matching returns NA
			
			# If no NA or a list  then proceed 	# is.list(patlist) & length(patlist$patst)>0
			if(! is.na(patlist[1]))
			{
				cat(paste("Number of motifs found in chain ",ch,": ",length(patlist),"\n\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE); 
				
				for(tempi in 1:length(patlist)) # Loop over each pattern
				{
					
					cat(paste("Processing motif:",tempi,"\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE); 
					pos=paste(patlist[[tempi]][1,1],patlist[[tempi]][1,2],sep="")
					patseq=paste(aa321(patlist[[tempi]][,3]),sep="",collapse="")
					tempdetail=c(filename=filepref,chain=ch,pos,patseq=patseq);
					cat(paste("Chain:",ch," Motif:",patseq,"\n\n",sep="",collapse=""),file=logfile,sep="",eol="",append=TRUE); 
					
					patlen=nrow(patlist[[tempi]]);
					
          # interaction and dihedral files per pattern
					txtdir$sumpath.int=paste(txtdir$subsubdirpath,"/","details.txt",sep="",collapse="");
					htmldir$sumpath.int=paste(htmldir$subsubdirpath,"/","details.html",sep="",collapse="");
					
					
					
          
					hlink.name=paste(hlinkpath,prefix,"/",subpref1,"/Details/","details.html",sep="",collapse="");
					hlink=paste("<a href=\"",hlink.name,"\" target=\"_blank\" >",filepref,"</a>");
					print(hlink);
					pat.acc=paste(rep("-",times=patlen),sep="",collapse="");
					pat.ss=pat.acc;
					pat.rama=pat.acc;
					tosearch=paste(ch,patlist[[tempi]][,1],patlist[[tempi]][,2],sep="")
	
					if(status.nacc){pat.acc=paste(apply(matrix(tosearch,ncol=1),1,patAddNaccess,nacc=nacc),sep="",collapse="")}

					if(status.DSSP){pat.ss=paste(apply(matrix(tosearch,ncol=1),1,patAddDSSP,DSSP=DSSP),sep="",collapse="")}
					if(status.procheck){pat.rama=paste(apply(matrix(tosearch,ncol=1),1,patAddRama,proch=proch),sep="",collapse=",")}

					tempdetail=c(tempdetail,pat.acc,pat.ss,pat.rama);
					
					cat("\n",file=txtdir$sumpath,sep="",eol="",append=TRUE);
					cat(paste(tempdetail,sep="",collapse="\t"),file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
					
					tempdetail.html=tempdetail;tempdetail.html[1]=hlink;
					write2Html("<tr>\n",htmldir$sumpath,append=TRUE);
					
          write2Html.tableBodyOnlytd(matrix(tempdetail.html,nrow=1),htmldir$sumpath,TRUE,"<td>");
					
          
					torsmat=t(apply(matrix(paste(patlist[[tempi]][,1],patlist[[tempi]][,2],ch,sep="."),ncol=1),1,torPat,pdb,tors));
					
					patmat=	cbind(ch,paste(patlist[[tempi]][,1],patlist[[tempi]][,2],sep=""),aa321(patlist[[tempi]][,3]));				
					pattors=cbind(patmat,torsmat);
					
					#cat("Dihedrals\n--------------\n\nChain\tRes.No\tRes.ID\tPhi\tPsi\tChi1\tChi2\tChi3\tChi4\tChi5\n",file=txtdir$sumpath.int);								
					#write.table(pattors,file=txtdir$sumpath.int,sep="\t",eol="\n",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE);
					
					
					#print_HTMLdihed(htmldir$sumpath.int,append=FALSE,html=TRUE);
					#write2Html.tableBody(pattors,htmldir$sumpath.int,TRUE,"<td>");
					#write2Html("</table>",htmldir$sumpath.int,TRUE);
					
          patmat=patmat[c(1,nrow(patmat)),];
					patmat=rbind(patmat,"-"); # just a trick
          
					print(patmat);					
					# Check for Interaction
					
					# Ionic # If any residue of pattern is D/E/R/K/H then only check for ionic
					if(intList$status.ionic & length(grep("[DERKH]",tempdetail[4],ignore.case = TRUE))>0)
					{ 
						ionmat=checkIonic(patmat,intList$ionpair);
						
            # Any intrachain hydrogen bonds will be resulted twice, so 
						uniq=unique(paste(ionmat[,1],ionmat[,2],ionmat[,3],ionmat[,4],ionmat[,5],ionmat[,6],ionmat[,7],ionmat[,8],ionmat[,9],ionmat[,10],sep="_"))
						ionmat=matrix(unlist(strsplit(uniq,"_")),ncol=10,byrow=TRUE);
            
						cat(nrow(ionmat),file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
						#cat("I",file=sumpath,sep="",eol="\t",append=TRUE);
						write2Html.tableBodyOnlytd(matrix(nrow(ionmat),nrow=1),htmldir$sumpath,TRUE,"<td>");
						
						if(nrow(ionmat)==0)
						{
						print("No ionic interaction");
						cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE); # no inopair network if no ionic interaction
						write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
						}else{
						# write ionmat
						cat("\nIONIC\n-----------\nChain\tRes1.No\tRes1.ID\tChain\tRes2.No\tRes2.ID\tRes1_Acc\tRes2_Acc\tRes1_SS\tRes2_SS\n",file=txtdir$sumpath.int,append=TRUE); # Create individual ionpair table file
						write.table(ionmat,sep="\t",eol="\n",file=txtdir$sumpath.int,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE);
						ionnetmat=checkIonnet(patmat,intList$ionpair,intList$dionpair,acclow); #check netowrk
						print("***");
						print(ionnetmat);
						print_HTMLionic(ionmat,ionnetmat,htmldir$sumpath.int,append=TRUE,FALSE);	
							if(nrow(ionnetmat)==0)
							{
							print("No ionic network"); cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
							write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
							}else{
							cat(nrow(ionnetmat),file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
							write2Html.tableBodyOnlytd(matrix(nrow(ionnetmat),nrow=1),htmldir$sumpath,TRUE,"<td>");
							cat("\nIONIC NETWORK\n----------------\nN.Res\tN.Int\tBuried\tExposed\tNetwork_Details\n----------------------------------------------------------------------------------\n",file=txtdir$sumpath.int,append=TRUE);
							write.table(ionnetmat,sep="\t",eol="\n",file=txtdir$sumpath.int,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE);
							}
						}
					
					}else{
					cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE); # for ip
					cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE); # for network
					write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
					write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
					print("NO IONPAIR DETECTED");
					}
					
					## Aro-Aro # If any residue of pattern is FWY then only check for aroaro
					if(intList$status.aroaro & length(grep("[FWY]",tempdetail[4],ignore.case = TRUE))>0)
					{ 
						aroaromat=checkAroAro(patmat,intList$aro);
						cat(nrow(aroaromat),file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
						#cat("A",file=sumpath,sep="",eol="\t",append=TRUE);
						write2Html.tableBodyOnlytd(matrix(nrow(aroaromat),nrow=1),htmldir$sumpath,TRUE,"<td>");
						
						if(nrow(aroaromat)==0)
						{
						print("No AroAro interaction"); 
						cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE); 
						write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
						}else{
						# To do ionmat
						cat("\nAromatic_Aromatic_Interaction\n--------------------------------\nChain\tRes1.No\tRes1.ID\tChain\tRes2.No\tRes2.ID\tC.Distance\tDihedral\tAcc_Res1\tAcc_Res2\tSS_Res1\tSS_Res2\n",file=txtdir$sumpath.int,append=TRUE); # Create individual ionpair table file
						write.table(aroaromat,sep="\t",eol="\n",file=txtdir$sumpath.int,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE);
						
							aroaronetmat=checkAroAronet(patmat,intList$aro,intList$daropair,acclow);
							print_HTMLaroaro(aroaromat,aroaronetmat,htmldir$sumpath.int,TRUE,FALSE);
							
							if(nrow(aroaronetmat)==0)
							{
							print("No aroaro network");cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
							write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
							}else{
							cat(nrow(aroaronetmat),file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
							write2Html.tableBodyOnlytd(matrix(nrow(aroaronetmat),nrow=1),htmldir$sumpath,TRUE,"<td>");
							cat("\nAromatic_Aromatic_Network\n------------------------\nN.Res\tN.Int\tBuried\tExposed\tNetwork_Detail\n----------------------------------------------------------------------------------\n",file=txtdir$sumpath.int,append=TRUE);
							write.table(aroaronetmat,sep="\t",eol="\n",file=txtdir$sumpath.int,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE);
							}
						}
					}else{
					cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
					cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
					write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
					write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
					print("NO AROPAIR DETECTED");
					}
					
					## Aro-S # If any residue of pattern is FWY then only check for aroaro
					if(intList$status.aros & length(grep("[FWYS]",tempdetail[4],ignore.case = TRUE))>0)
					{ 
					arosmat=checkAroS(patmat,intList$aros);
					cat(nrow(arosmat),file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
					#cat("S",file=sumpath,sep="",eol="\t",append=TRUE);
					write2Html.tableBodyOnlytd(matrix(nrow(arosmat),nrow=1),htmldir$sumpath,TRUE,"<td>");
						
						if(nrow(arosmat)==0){print("No AroS interaction"); cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE);  
						write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
						}else{
						# To do ionmat
						cat("\nAromatic_sulphor Interaction\n--------------------------------\nChain\tRes1.No\tRes1.ID\tChain\tRes2.No\tRes2.Id\tDistance\tAcc_Res1\tAcc_Res2\tSS_Res1\tSS_Res2\n",file=txtdir$sumpath.int,append=TRUE); # Create individual ionpair table file
						write.table(arosmat,sep="\t",eol="\n",file=txtdir$sumpath.int,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE);
						
						arosnetmat=checkAroSnet(patmat,intList$aros,intList$darospair,acclow);
						print_HTMLaroS(arosmat,arosnetmat,htmldir$sumpath.int,TRUE,FALSE);
						
							if(nrow(arosnetmat)==0){print("No aros network"); cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE); 
							write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
							}else{
							cat(nrow(arosnetmat),file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
							write2Html.tableBodyOnlytd(matrix(nrow(arosnetmat),nrow=1),htmldir$sumpath,TRUE,"<td>");
							cat("\nAromatic_sulphor Interaction Network\n----------------------------------\nN.Res	N.Int	Buried	Partial	Exposed	Network_Detail\n----------------------------------------------------------------------------------\n",file=txtdir$sumpath.int,append=TRUE);
							write.table(arosnetmat,sep="\t",eol="\n",file=txtdir$sumpath.int,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE);
							}	
						}
					}else{
					cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
					cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
					write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
					write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
					print("NO AROS DETECTED");
					}
					
				  if(intList$status.hb)
					{
					  hbmat=checkHB(patmat,intList$hb);
					  
					  hbmat=matrix(unlist(strsplit(unique(apply(hbmat,1,function(x){paste(x,sep="",collapse=":")})),":")),ncol=14,byrow=TRUE) #remove duplicate
            print(patmat);
            print(hbmat);
					  ind1=which(hbmat[,2] %in% patmat[1:2,2]);
            print(ind1);
					  ind2=which(hbmat[,6] %in% patmat[1:2,2]);
					  print(ind2);
					  ind3=intersect(ind1,ind2);
					  print(ind3);
					  
					  if(length(ind3)==0){hbmat=matrix(numeric(0), 0,0)}else{hbmat=matrix(hbmat[ind3,],ncol=ncol(hbmat))}
					  cat(nrow(hbmat),file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
					  write2Html.tableBodyOnlytd(matrix(nrow(hbmat),nrow=1),htmldir$sumpath,TRUE,"<td>");
					  #cat("H",file=sumpath,sep="",eol="\t",append=TRUE);
					  if(nrow(hbmat)==0){print("No HB interaction");}else{
					    write2Html.tableBody(matrix(tempdetail,nrow=1),htmldir$sumpath.int,TRUE,"<td>");
					    cat("\n",file=txtdir$sumpath.int,sep="",eol="\n",append=TRUE);
              cat(paste(tempdetail,sep="",collapse="\t"),file=txtdir$sumpath.int,sep="",eol="\n",append=TRUE);
              # To do ionmat
					    cat("\nHydrogen Bonds\n-----------------\nD.Chain\tD.Resno\tD.Resid\tD.atom\tA.Chain\tA.Resno\tA.Resid\tA.atom\tBond_type\tDistance_DA\tDistance_HA\tDHA_Angle\tHAAA_Angle\tDAAA_Angle\n",file=txtdir$sumpath.int,append=TRUE); # Create individual ionpair table file
					    write.table(hbmat,sep="\t",eol="\n",file=txtdir$sumpath.int,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE);
					    print_HTMLhb(hbmat,htmldir$sumpath.int,append=TRUE,html=FALSE,"<br><br><b><u>Hydrogen bonding interaction</u><br><br></b>\n");
					  }
					}else{
					  cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
					  write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
					  print("NO HB DETECTED");
					}
					
					# Disulfide
					
					if(intList$status.dis & length(grep("C",tempdetail[4],ignore.case = TRUE))>0)
					{
					dismat=checkDis(patmat,intList$dis);
					cat(nrow(dismat),file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
					#cat("D",file=sumpath,sep="",eol="\t",append=TRUE);
					write2Html.tableBodyOnlytd(matrix(nrow(dismat),nrow=1),htmldir$sumpath,TRUE,"<td>");
						if(nrow(dismat)==0){print("No dis interaction");}else{
						# To do ionmat
						cat("\nDISULFIDE\n------------------\nChain\tRes1_No\tChain\tRes2_No\tAcc_Res1\tAcc_Res2\tSS_Res1\tSS_Res2\n",file=txtdir$sumpath.int,append=TRUE); # Create individual ionpair table file
						write.table(dismat,sep="\t",eol="\n",file=txtdir$sumpath.int,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE);
						print_HTMLdis(dismat,htmldir$sumpath.int,TRUE,FALSE);
						}
					}else{
					cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
					write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
					print("NO DIS DETECTED");
					}
					
					# Catpi
					if(intList$status.catpi & length(grep("[FWYRK]",tempdetail[4],ignore.case = TRUE))>0)
					{
					catpimat_local=checkCatpi(patmat,intList$catpimat);
						cat(nrow(catpimat_local),file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
						#cat("P",file=sumpath,sep="",eol="\t",append=TRUE);
						write2Html.tableBodyOnlytd(matrix(nrow(catpimat_local),nrow=1),htmldir$sumpath,TRUE,"<td>");
						
						if(nrow(catpimat_local)==0){print("No Catpi interaction"); cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE); 
						write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
						}else{
						cat("\nCat_pi interaction\n------------------\nChain\tRes1_No\tRes1_ID\tChain\tRes2_No\tRes2_ID\tDistance\tAcc_Res1\tAcc_Res2\tSS_Res1\tSS_Res2\n",file=txtdir$sumpath.int,append=TRUE); # Create individual ionpair table file
						write.table(catpimat_local,sep="\t",eol="\n",file=txtdir$sumpath.int,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE);
						
						catpinet=checkcatpinet(patmat,intList$catpimat,intList$dcatpimat,acclow);
						print_HTMLcatpi(catpimat_local,catpinet,htmldir$sumpath.int,TRUE,FALSE);
							if(nrow(catpinet)==0){print("No catpi network"); cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE); 
							write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
							}else{
							cat(nrow(catpinet),file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
							write2Html.tableBodyOnlytd(matrix(nrow(catpinet),nrow=1),htmldir$sumpath,TRUE,"<td>");
							cat("\nCatpi Interaction Network\n----------------------------------\nN.Res	N.Int	Buried	Partial	Exposed	Network_Detail\n----------------------------------------------------------------------------------\n",file=txtdir$sumpath.int,append=TRUE);
							write.table(catpinet,sep="\t",eol="\n",file=txtdir$sumpath.int,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE);
							}
						}
					}else{
					cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE); #no catpi
					cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE); #no catpi-net
					write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
					write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
					print("NO CATPI DETECTED");
					}
					
					# Hydrophobic
					if(intList$status.hp & length(grep("[AVLIMFWPY]",tempdetail[4],ignore.case = TRUE)))
					{
					hpmat=checkHphob(patmat,intList$hp);
					cat(nrow(hpmat),file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
					#cat("B",file=sumpath,sep="",eol="\t",append=TRUE);	
					write2Html.tableBodyOnlytd(matrix(nrow(hpmat),nrow=1),htmldir$sumpath,TRUE,"<td>");
						if(nrow(hpmat)==0){print("No Hydrophobic interaction");}else{
						cat("\nHydrophobic interaction\n------------------\nChain\tRes1_No\tRes1_ID\tChain\tRes2_No\tRes2_ID\tDistance\tAcc_Res1\tAcc_Res2\tSS_Res1\tSS_Res2\n",file=txtdir$sumpath.int,append=TRUE); # Create individual ionpair table file
						write.table(hpmat,sep="\t",eol="\n",file=txtdir$sumpath.int,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE);
						print_HTMLhphob(hpmat,htmldir$sumpath.int,TRUE,FALSE);
						}
					}else{
					cat("0",file=txtdir$sumpath,sep="",eol="\t",append=TRUE);
					write2Html.tableBodyOnlytd(matrix("0",nrow=1),htmldir$sumpath,TRUE,"<td>");
					print("NO Hphob DETECTED");
					}
					
					write2Html("\n</tr>\n",htmldir$sumpath,append=TRUE);
					write2Html("</html>",htmldir$sumpath.int,append=TRUE);
				}
			}
		}
	}
print("hii");
# remove the working directory
setwd(htmldir$dirpath);
#unlink(work_dir, recursive = TRUE);
system(paste("zip -r ","./result.zip ", "./",subpref2,sep="",collapse=""),ignore.stderr=FALSE);

cat("100",file=statusfile,sep="",eol="");


