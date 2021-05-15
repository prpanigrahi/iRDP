naccess=function(pdb,exepath="",prefix="")
{
print("Running NACCESS")    
print(unique(pdb$atom[,"chain"]))
infile=paste(prefix,".pdb",sep="",collapse="");
    rsafile=paste(prefix,".rsa",sep="",collapse="");
    asafile=paste(prefix,".asa",sep="",collapse="");
    logfile=paste(prefix,".log",sep="",collapse="");
    write.pdb(pdb, file = infile);
    tryCatch({error=system(paste(exepath, " ",infile, sep = ""),ignore.stderr = FALSE,ignore.stdout =FALSE,intern=TRUE)},error=function(e){errorflag=1;return("SORRY");});
    #In case of 1VTZ.pdb nacc fails with error "STOP SOLVA_ERROR: max cubes exceeded statement executed"
    # But it doesnot captured in error variable instead we get nacc list object, however nacc$asa doesnot have any information
    # 2 way of error handeling: 
    # 1) check error variable which is implemented below
    # 2) nacc$asa which is same as finalasa variable, which is implemented at the last in return statement

    #print(attributes(error));
    if(is.list(attributes(error))){return("SORRY");}  

    if(!(file.exists(rsafile) & file.exists(asafile)))
    {
     # print("sorry")
    	return("SORRY");
    }
    
    raw.lines <- readLines(rsafile);
    atom_acc=read.pdb(asafile);
    unlink(c(infile,logfile,asafile))
    id=substring(raw.lines,1,3);
    TOT=unlist(strsplit(raw.lines[which(id=="TOT")]," +",perl=TRUE));
    RES=raw.lines[which(id=="RES")];

    resid=substring(RES, 5, 7)
    ch=substring(RES, 9, 9)
    resno=substring(RES, 10, 14) #resno.insert
    resno=gsub(" ","",resno) # replace space with ""
  
    # Residue ASA
    asa_r_abs=as.numeric(substring(RES, 16, 22))
    asa_r_rel=as.numeric(substring(RES, 23, 28))
    # SideChain ASA
    asa_s_abs=as.numeric(substring(RES, 29, 35))
    asa_s_rel=as.numeric(substring(RES, 36, 41))
    # MainChain ASA
    asa_m_abs=as.numeric(substring(RES, 42, 48))
    asa_m_rel=as.numeric(substring(RES, 49, 54))
    # NP ASA
    asa_np_abs=as.numeric(substring(RES, 55, 61))
    asa_np_rel=as.numeric(substring(RES, 62, 67))
    # NP ASA
    asa_p_abs=as.numeric(substring(RES, 68, 74))
    asa_p_rel=as.numeric(substring(RES, 75, 80))

finalasa= ch;
finalasa=cbind(finalasa,resno) ;
finalasa=cbind(finalasa,resid) ;
finalasa=cbind(finalasa,asa_r_abs) ;
finalasa=cbind(finalasa,asa_r_rel) ;
finalasa=cbind(finalasa,asa_s_abs) ;
finalasa=cbind(finalasa,asa_s_rel) ;
finalasa=cbind(finalasa,asa_m_abs) ;
finalasa=cbind(finalasa,asa_m_rel) ;
finalasa=cbind(finalasa,asa_np_abs) ;
finalasa=cbind(finalasa,asa_np_rel) ;
finalasa=cbind(finalasa,asa_p_abs) ;
finalasa=cbind(finalasa,asa_p_rel) ;
colnames(finalasa)[1]="Chain";
print("Finished NACCESS")
if(length(finalasa)==0){return("SORRY")}else{return(list(asa=finalasa,tot=TOT,atomacc=atom_acc));}
}

