iMut=function(iMuttemppdb,iMuttempdssp,chain,resno,new,pH=7.0,temp=25)
{
	 print(c("Mutation running belongs to ",chain,resno,new));
	 command=paste(iMutPythonPath,"-O",iMutPath,"-pdb",iMuttemppdb,iMuttempdssp,chain,resno,pH,temp,sep=" ");
	 print(command);
	 result=system(command,ignore.stderr = FALSE,ignore.stdout=FALSE,intern=TRUE);
	 print(result);
	 st=(grep("  Position",result)+1)
	 temp.result=result[st:(st+18)];
	 final.result=temp.result[which(substring(temp.result,20,20)==new)]
# I-mutant predicts either "sign of DDG" or "DDG value" by using 2 different svm program
# Therefore sometimes they give different result when RI is very low
# So we should use either Sign+RI or DDG value
# Reliability index and DDG but rel index is used for sign 
# If we use DDG value then if DDG<0, decrease stability else increase
# So we have decided to use DDG value and from that value we decide I/D of stability
	 temp.imut=c(substring(final.result,36,36),substring(final.result,40,44));
	 temp.imut=as.numeric(gsub(" ","",temp.imut));
	 
	 return(temp.imut); 
}

#print(iMut(pdb,"B",1,"A"));

