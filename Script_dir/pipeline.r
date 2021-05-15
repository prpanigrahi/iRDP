####################### All directories
script_dir="/var/www/SASMAT/";
work_dir="/var/www/upload";
pdb_dir="/var/www/pdb";
output_dir="/var/www/html/";

######################## Arguments
args=commandArgs(trailingOnly = TRUE);
filename=tolower(paste(args[1],".pdb",sep="",collapse=""));
print(filename);
pat=args[2];

######################## All loadings
library("bio3d");
source(paste(script_dir,"patfind.r",sep=""));
source(paste(script_dir,"torpat.r",sep=""));
source(paste(script_dir,"general.r",sep=""));
source(paste(script_dir,"PATH.r",sep=""));
source(paste(script_dir,"dssp.r",sep=""));
source(paste(script_dir,"naccess.r",sep=""));

######################## Set working directory
setwd(work_dir);
#print(getwd());

######################## List of all pdb files
files=list.files(pdb_dir,pattern=".pdb$");
files_low=tolower(files);
#print(files);

######################## All output files
finalout="";
outfile=paste(output_dir,"Result.txt",sep="",collapse="");
outimg=paste(output_dir,"temp.png",sep="",collapse="");
ramafile="Rama.txt";
#print(outfile);

file.ind=which(files_low==filename);
#print(length(file.ind)>0);

if(length(file.ind)>0)
{
# Empty the work_dir then copy the pdb from pdb dir
delfiles=list.files(work_dir);
delfiles=delfiles[-c(which(delfiles=="Ramachandran.java"),which(delfiles=="Ramachandran.class"),which(delfiles=="Ramachandran$1.class"))]
if(length(delfiles)>0){file.remove(delfiles);}

file.copy(from=paste(pdb_dir,"/",files[file.ind],sep="",collapse=""),to=work_dir);
#print(list.files(pdb_dir,pattern=".pdb$"));
pdb=read.pdb(files[file.ind],verbose=FALSE);
 
#print(str(pdb));
chains=uniqueChain(pdb); # List of unique chains in which pattern is to be found, if no chain then NA will be there
print(chains);

tors=torsion.pdb(pdb);
#print(str(tors));
DSSP=dssp_new(pdb,exepath=dssppath ,infile="temp.pdb",outfile="temp.dssp");
#print(str(DSSP));
nacc=naccess(pdb,naccpath,infile="temp.pdb",rsafile="temp.rsa");
# We go chainwise.
# For each chain, call to patFind() will return pdb resno of start residue of pattern.
# If no pattern found it returns NA, if more than one found then it returns a vector
# So we loop over each pattern position and extract its details
# First we need to fetch the torsion information.
# Inside the Loop, for each pattern position in a given chain we get "chain","pos","pattern"
# Then we call torpat which returns a vector of all torsion detail, 8 per residue, phi,psi,omega, chi1 to chi5

ramaPlot=function(finalout="")
{
phi=seq(from=5,to=ncol(finalout),by=9)
psi=seq(from=6,to=ncol(finalout),by=9)
rama=cbind(as.vector(finalout[,phi]),as.vector(finalout[,psi]));
#print(rama);
# Store the phi psi value in script_dir, execute Ramachandran java program, move the png to output_dir
write.table(rama,file=ramafile,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t");
system(paste("java Ramachandran -f ", ramafile,sep="",collapse=""),ignore.stderr = FALSE);
file.rename("ramachandran.png",paste(output_dir,"ramachandran.png",sep="",collapse=""));
#png(file=outimg);
#plot(as.vector(finalout[,phi]),as.vector(finalout[,psi]),xlim=c(-180,180),ylim=c(-180,180),pch=19,xlab="Phi",ylab="Psi",main="Ramachandran plot for pattern residues");
#lines(c(-200,200),c(0,0));
#lines(c(0,0),c(-200,200));
#dev.off();
}

for(ch in chains)
{
	#print(pat);
	patlist=patFind(pdb,ch,pat)# It returns actual residue no information in pdb file, if not matching returns NA
	patresno=patlist$patst;
	patlength=patlist$patlen; #pattern length
	# Sort the pattern according to length
	patlength[order(patlength)]
	
	#print(patresno);
	# If no NA then proceed for torsion finding
	if(!is.na(patresno[1]))
	{
		
		for(tempi in 1:length(patresno))
		{
			pos=patresno[tempi];
			patlen=patlength[tempi]
			
			tempdetail=c(filename,ch,pos);
			if(is.na(ch) | ch=="~" | ch==" "){atomsel=atom.select(pdb,elety="CA",resno=pos:(pos+patlen-1),verbose=FALSE);}else{
			atomsel=atom.select(pdb,chain=ch,resno=pos:(pos+patlen-1),elety="CA",verbose=FALSE);
			}
			tempdetail=c(tempdetail,paste(aa321(pdb$atom[atomsel$atom,"resid"]),collapse="")); # pattern
			for(st in pos:(pos+patlen-1))
			{
			tempdetail=c(tempdetail,round(torPat(pdb,tors,ch,st),2));
			if(is.list(nacc)){tempdetail=c(tempdetail,patAddNaccess(pdb,nacc,ch,st));}else{tempdetail=c(tempdetail,"NA");}
			if(is.list(DSSP)){tempdetail=c(tempdetail,patAddDSSP(pdb,DSSP,ch,st));}else{tempdetail=c(tempdetail,"NA");}
			}
			#names(tempdetail)=
			finalout=rbind(finalout,tempdetail);
			
			#cat(tempdetail,file=outfile,sep="\t",eol="\n",append=TRUE);
		}
	}

}


}else{print("PDB FILE NOT IN OUR DATABASE, PLEASE USE UPLOAD FEATURE");}

#print(is.matrix(finalout));
if(is.matrix(finalout))
{
#print(finalout);
cat(c("file","chain","pos","pat",rep(c("phi","psi","chi1","chi2","chi3","chi4","chi5","ACC","SS"),times=patlen)),file=outfile,sep="\t",eol="\n",append=TRUE);
nc=ncol(finalout);
finalout=finalout[-1,];
finalout=matrix(finalout,ncol=nc);
write.table(finalout,file=outfile,sep="\t",eol="\n",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE);
ramaPlot(finalout);
}else{
}

#seq(from=5,to=ncol(finalo))

#phi=
#options(device="png")
#png(file=outimg)
#plot(finalout[,5],finalout[,6]);
#title("hiiiii");
#plot(finalout[,6],finalout[,6]);
#dev.off()
