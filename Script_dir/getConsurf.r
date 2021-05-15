library("bio3d");
setwd("/var/www/ent/");
filenames=list.files(pattern=".ent$");

for(f in filenames[1:2])
{
  pdb="";
  count=which(filenames==f);
  pdbfile=f;
  filepref=substr(f,4,7); # for ent format
  pdb=read.pdb(pdbfile,maxlines=1000000,verbose=FALSE); # read the pdb file
  chains=unique(names(pdb$seqres));
  #cat(paste(c(filepref,chains),sep="",collapse="\t"),file="list_for_getConsurf.txt",eol="\n",append=TRUE);
  for(ch in chains)
  {
  print(c(count,filepref,ch));
  system();
  }

}

