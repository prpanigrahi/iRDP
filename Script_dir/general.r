aroCentroid=function(pdb="",restype="")
{
	if(restype=="PHE" | restype=="TYR" | restype=="TRP")
	{
	phe.resIndex=atom.select(pdb,resid=restype,elety="CA",verbose=FALSE)$atom #SELECT pdb$atom index corresponding to CA atom of PHE/Y/W
	
		if(length(phe.resIndex)==0){
		#sometimes no Phe/Tyr/Trp, then no chance of centroid. Hence return with "SORRY"
		return("SORRY");
		}
	# Get the pdb.resIndex, call the getRing.index function by apply. It takes 3 Arg, chain,resno,restype and returns either 6 ring atom index else returns rep("-",6) in case of anything missing	
	# So we check if 6 gaps are there that means for that residue, mainchain there but sidechain not there
	
	phe.ringIndex=apply(cbind(matrix(pdb$atom[phe.resIndex,c("chain","resno")],ncol=2),restype),1,getRing.index,pdb);
	phe.ringIndexGap=apply(phe.ringIndex,2,function(x){sum(x!="-")>0});
	
	phe.resIndex=as.numeric(phe.resIndex[phe.ringIndexGap]);
	phe.ringIndex=as.numeric(as.vector(phe.ringIndex[,phe.ringIndexGap]));
	
	if(0)
	{
	# This was the older version which worked well. But suddernly in 2VS7.pdb file, In D chain although residue Tyr 13 have sidechain but it doesnot have CA atom
	# Therefore  in phe.resIndex atomselect option, D13Tyr could not come but as it does have side chain in phe.ringIndex atomselect it had come as indexes
	# We had implemented the condition where residue having CA but if no side chain then discard
	# But surpriosingly we found that residue with sidechain but no mainchain atom [Actually mainchain does present but it doesnot parse, jut see pdb file u will understand]	
	# So a better way it was solved, described above
	
	
		if(restype=="TRP")
		{
		phe.ringIndex=atom.select(pdb,resid=restype,elety=c("CD2","CE2","CE3","CZ2","CZ3","CH2"),verbose=FALSE)$atom #select pdb$atom index of ring atoms of TRP
		}else{
		phe.ringIndex=atom.select(pdb,resid=restype,elety=c("CG","CD1","CD2","CE1","CE2","CZ"),verbose=FALSE)$atom #select pdb$atom index of ring atoms of F/Y
		}
		
		# Sometimes although residue is there but side chain ring atoms are missing.
		# Suppose 3 Phe residue, we should get 18 ring atoms, if not then definitely some atoms are missing.
		# With any one atom missing centroid calcualtion may go wrong. Hence we must identify that residue and discard.
		# condition below will check if length(phe.resIndex)*6 {no of phe res*6 atoms == no of ring atoms}
		# If condition fails, we loop over each residue
		if((length(phe.resIndex)*6) != length(phe.ringIndex)){
			i=1;
			while(i <= length(phe.resIndex))
			{
			#print(i);
				# we have phe.resIndex, which stores pdb$atom index for all CA atom of restype PHE/Y/W.
				# using that pdb$atom index, we fetch this residue's chain and resno and using which we
				# try to fetch all ring atoms, if length is not 6 then we are sure that side chain ring atoms are missing.
				# So we delete that index
				tempcha=pdb$atom[phe.resIndex[i],"chain"]	
				tempres=pdb$atom[phe.resIndex[i],"resno"]
				print(length(atom.select(pdb,chain=tempcha,resno=tempres,elety=c("CG","CD1","CD2","CE1","CE2","CZ"),verbose=FALSE)$atom));
				if(length(atom.select(pdb,chain=tempcha,resno=tempres,elety=c("CG","CD1","CD2","CE1","CE2","CZ"),verbose=FALSE)$atom)!=6){
					phe.resIndex=phe.resIndex[-i]; #remove ith index
					print(i);
				}else{
				i=i+1;
				}
			}
		}
	}

	if(length(phe.resIndex)==0){
		# after removing some index or F/W/Y residue whose side chain missing, if no Phe/Tyr/Trp, then no chance of centroid. Hence return with "SORRY"
		return("SORRY");
		}
		
		phe.detail=matrix(pdb$atom[phe.resIndex,c("chain","resno","insert","resid")],ncol=4) # using the above phe.resIndex get details of chain and resno
		phe.detail[which(is.na(phe.detail[,3])),3]="";
		phe.detail[,2]=paste(phe.detail[,2],phe.detail[,3],sep="")
		phe.detail=matrix(phe.detail[,c(1,2,4)],ncol=3);
		

		# If any chain is not there then we get NA, so make it "-"
		# phe.detail[which(is.na(phe.detail[,1])),1]="~";
		
		# using phe.ringIndex, fetch n*3 (n: no of atom, 3:x,y,z) matrix, convert to as.numeric which create a 1D vector, convert back to matrix
		phe.ringXYZ=matrix(as.numeric(pdb$atom[phe.ringIndex,c("x","y","z")]),ncol=3) 
		# Suppose there are 10 Phe residue, phe.resIndex size will be 10, phe.details will be 10*3, phe.ringIndex will be 60 (10 residue and 6 ring atoms)
		# phe.ringXYZ will be 60*3 matrix [60 ring atoms, 3 xyz coordinates] i.e 180 elements
		# we split the 60*3 matrix into a 3D 3*6*10 phe.ringXYZmat where 3: x,y,z 6: 6 ring atoms, 10: 10 phenyl residue
		# phe.ringXYZmat[,,1] is a 3*6 matrix corresponding to xyz coordinates of 6 ring atoms of phe residue1
		phe.ringXYZmat=array(t(phe.ringXYZ),dim=c(3,6,length(phe.ringXYZ)/18))

		# Loop over each phe residue i.e phe.ringXYZmat[,,i], input the 3*6 matrix to ringCentroid() function and get a 3 element vector, rbind to centroidMat	
		# Same 3*6 matrix we pass to ringEQN function, we get another 3 element 1D vector.
		centroidMat=0;
		EQNmat=0;
		for(i in 1:dim(phe.ringXYZmat)[3])
		{
		centroidMat=rbind(centroidMat,ringCentroid(phe.ringXYZmat[,,i]))
		EQNmat=rbind(EQNmat,ringEQN(phe.ringXYZmat[,,i]));
		}
		centroidMat=centroidMat[-1,]
		EQNmat=EQNmat[-1,];
		# We have to convert to matrix bcoz if one residue is there then one ring, so one centroid x,y,z and it will be 1D vector rather than a matrix, hence after
		# removal of -1, we get a vector of 3 element. so we must convert to matrix else rbind gives error.
		centroidMat=matrix(centroidMat,ncol=3)
		EQNmat=matrix(EQNmat,ncol=3)
	
		#centroidMat=round(centroidMat,3) #round to 3 decimal place
		#EQNmat=round(EQNmat,3)
	
	
		phe.detail=cbind(phe.detail,centroidMat)
		phe.detail=cbind(phe.detail,EQNmat)
	
		colnames(phe.detail)=c("chain","resno","resid","Cx","Cy","Cz","EqnA","EqnB","EqnC")
		return(phe.detail);
	}else{
	stop("Your input should be PHE or TYR or TRP, case sensitive. Check aroaro() function");
	}
	
}
createDir=function(dirpath)
{
	if (!file.exists(dirpath)){dir.create(dirpath);}
}
getRing.index=function(x,pdb="")
{
	chain=x[1];resno=x[2];restype=x[3];
	atomlabels="";
	ring.atoms="";
	if(restype=="TRP"){atomlabels=c("CD2","CE2","CE3","CZ2","CZ3","CH2")}else{atomlabels=c("CG","CD1","CD2","CE1","CE2","CZ")}

	if(!is.na(chain))
	{
	ring.atoms=atom.select(pdb,chain=chain,resno=resno,elety=atomlabels,verbose=FALSE)$atom
	}else{
	ring.atoms=atom.select(pdb,resno=resno,elety=atomlabels,verbose=FALSE)$atom
	}

	if(length(ring.atoms)==6)
	{
	return(ring.atoms);
	}else{
	return(rep("-",times=6));
	}
}

Nochain=function(aro,chain1index,chain2index)
{
#for aroaro chain1index=1, chain2index=4
aro[is.na(aro[,chain1index]),chain1index]="~";
aro[is.na(aro[,chain2index]),chain2index]="~";
return(aro);
}
ringANGL=function(data1,data2)
{
	#data1 and data2 are 1*4 vectors where 4 corresponds to A,B, C, D coefficients of plane
	#data1 is first plane coeff and data2 is second plane coeff.
	return(57.2957795*acos(  sum(as.numeric(data1)*as.numeric(data2))  / (sqrt(sum(as.numeric(data1)**2)) * sqrt(sum(as.numeric(data2)**2)))));
}
ringCentroid=function(data)
{
	#data must be 3*n dimension where 3 xyz and n for n atoms
	# x1+x2+x3+...+xn / n so data[1,] is x-coordinats of all 6 ring atoms.sum(data[1,]) is their sum. length(data[1,]) is 6 ring atoms
	# so x1+x2+...+x6/6 is x-coordinates of centroid.
	# It returns a 1D vector of x,y,z component of centroid
	return(c(sum(data[1,])/length(data[1,]),sum(data[2,])/length(data[2,]),sum(data[3,])/length(data[3,])))
}
ringEQN=function(data)
{
	# data[1,1]=x1	data[1,2]=x2	data[1,3]=x3
	# data[2,1]=y1	data[2,2]=y2	data[2,3]=y3
	# data[3,1]=z1	data[3,2]=z2	data[3,3]=z3
	
	A= (data[2,1]*(data[3,2]-data[3,3])) + (data[2,2]*(data[3,3]-data[3,1])) + (data[2,3]*(data[3,1]-data[3,2]))
	B= (data[3,1]*(data[1,2]-data[1,3])) + (data[3,2]*(data[1,3]-data[1,1])) + (data[3,3]*(data[1,1]-data[1,2]))
	C= (data[1,1]*(data[2,2]-data[2,3])) + (data[1,2]*(data[2,3]-data[2,1])) + (data[1,3]*(data[2,1]-data[2,2]))
#	D= -((data[1,1]* ((data[2,2]*data[3,3])-(data[2,3]*data[3,2]))) + (data[1,2]* ((data[2,3]*data[3,1])-(data[2,1]*data[3,3]))) + (data[1,3]* ((data[2,1]*data[3,2])-(data[2,2]*data[3,1]))))
	return(c(A,B,C))
}

screenDNA=function(pdb)
{
# If returns the name of DNA chains
	if(length(pdb$seqres)>0)
	{
	return(unique(names(which(nchar(pdb$seqres)!=3))));
	}else{
	return(unique(pdb$atom[which(nchar(pdb$atom[,"resid"])!=3),"chain"]));
	}
}

screenHetChain=function(pdb)
{
#screen any hetero atom chains
# example: 176d in which chain A is a heteratom chain
return(unique(pdb$het[,"chain"]));
}
TrimDNA=function(pdb)
{
# If any DNA chains are there, it will trim those chains and return a newpdb object
# It first lists all chains, identifies dna chains then lists protein chains only
# Then it checks if no chain info i.e. NA, then it returns entire pdb object as it is
# Else trim and return trimmed pdb
	chains=unique(pdb$atom[,"chain"]);
	dnach=screenDNA(pdb);
	chains=setdiff(chains,dnach);
	# After trimming DNA, if remaining doesnot have any chain info then return pdb object
	if(length(chains)==0){return("SORRY")}
	# if(sum(is.na(chains))>0){return(pdb);}
	
	newpdb=trim.pdb(pdb,inds=atom.select(pdb,chain=chains,verbose=FALSE));
}
uniqueChain=function(pdb)
{
#print("uniqch");
	# Returns unique chains
	chains="";
	seqres=0;
	if(length(pdb$seqres)>0)
	{
	seqres=1;
	chains=unique(names(pdb$seqres));
	}else{
	chains=unique(pdb$atom[,"chain"]) #Lists all chains
	}

	dnach=screenDNA(pdb); # Lists all DNA chains
	chains=setdiff(chains,dnach); # Filters DNA chains from all chains
	#hetch=screenHetChain(pdb);
	#chains=setdiff(chains,hetch); # Filters hetero atomchains from all chains

	# If after trimming only one chain lefts then no need to check for unique chains
	# we can return that one chain information.
	if(length(chains)<2){return(chains);}
	
	# After filtering of dna chains if no of chains are more than 1 then we go for uchain identification
	uchains=chains[1] # array of uchains
	
	for(i in 2:length(chains))
	{
			flag=FALSE;
			for(j in 1:length(uchains))
			{
				#print(c(i,j));	#print(uchains);
				if(!seqres)
				{
				# Take the coordinates from atom
				#print("Atoms");
				seqres1=paste(aa321(pdb$atom[atom.select(pdb,chain=chains[i],elety="CA",verbose=FALSE)$atom,"resid"]),sep="",collapse="");
				seqres2=paste(aa321(pdb$atom[atom.select(pdb,chain=uchains[j],elety="CA",verbose=FALSE)$atom,"resid"]),sep="",collapse="");	
				}else{
				#Take the coordinate from seqres
				#print("Seqres");
				 seqres1=paste(aa321(pdb$seqres[which(names(pdb$seqres)==chains[i])]),sep="",collapse="")
				 seqres2=paste(aa321(pdb$seqres[which(names(pdb$seqres)==uchains[j])]),sep="",collapse="")
				}
				
				if(seqres1==seqres2)
				{
				#print(c(chains[i],chains[j],"diff"));
					flag=TRUE;
					break;
				}
				
				
			}
			if(flag==FALSE)
			{
				uchains=c(uchains,chains[i])		
			}
	}
		return(uchains);	
}

uniqueChain_map=function(pdb)
{
  #print("uniqch");
  # Returns unique chains
  chains="";
  seqres=0;
  if(length(pdb$seqres)>0)
  {
    seqres=1;
    chains=unique(names(pdb$seqres));
  }else{
    chains=unique(pdb$atom[,"chain"]) #Lists all chains
  }
  
  dnach=screenDNA(pdb); # Lists all DNA chains
  chains=setdiff(chains,dnach); # Filters DNA chains from all chains
  #hetch=screenHetChain(pdb);
  #chains=setdiff(chains,hetch); # Filters hetero atomchains from all chains
  
  # If after trimming only one chain lefts then no need to check for unique chains
  # we can return that one chain information.
  if(length(chains)<2){return(matrix(c(chains,chains),ncol=2));}
  
  # After filtering of dna chains if no of chains are more than 1 then we go for uchain identification
  
  uchains=chains[1] # array of uchains
  chmat=c(chains[1],chains[1]); #chain matrix
  for(i in 2:length(chains))
  {
    flag=FALSE;
    for(j in 1:length(uchains))
    {
      #print(c(i,j));	#print(uchains);
      if(!seqres)
      {
        # Take the coordinates from atom
        #print("Atoms");
        seqres1=paste(aa321(pdb$atom[atom.select(pdb,chain=chains[i],elety="CA",verbose=FALSE)$atom,"resid"]),sep="",collapse="");
        seqres2=paste(aa321(pdb$atom[atom.select(pdb,chain=uchains[j],elety="CA",verbose=FALSE)$atom,"resid"]),sep="",collapse="");	
      }else{
        #Take the coordinate from seqres
        #print("Seqres");
        seqres1=paste(aa321(pdb$seqres[which(names(pdb$seqres)==chains[i])]),sep="",collapse="")
        seqres2=paste(aa321(pdb$seqres[which(names(pdb$seqres)==uchains[j])]),sep="",collapse="")
      }
      
      if(seqres1==seqres2)
      {
        #print(c(chains[i],chains[j],"diff"));
        chmat=rbind(chmat,c(chains[i],chains[j]));
        flag=TRUE;
        break;
      }
      
      
    }
    if(flag==FALSE)
    {
      uchains=c(uchains,chains[i])	
      chmat=rbind(chmat,c(chains[i],chains[i]));
    }
  }
  return(chmat);	
}

uniqueChain_old=function(pdb)
{
# Returns unique chains
	chains=unique(pdb$atom[,"chain"]) #Lists all chains
	dnach=screenDNA(pdb); # Lists all DNA chains
	chains=setdiff(chains,dnach); # Filters DNA chains from all chains
	# If after trimming only one chain lefts then no need to check for unique chains
	# we can return that one chain information.
	if(length(chains)<2){return(chains);}
	
	# After filtering of dna chains if no of chains are more than 1 then we go for uchain identification
	uchains=chains[1] # array of uchains
	for(i in 2:length(chains))
	{
			flag=FALSE;
			for(j in 1:length(uchains))
			{
				#print(c(i,j));	#print(uchains);

				# Take the coordinates from atom
				seqres1=paste(aa321(pdb$atom[atom.select(pdb,chain=chains[i],elety="CA",verbose=FALSE)$atom,"resid"]),sep="",collapse="");
				seqres2=paste(aa321(pdb$atom[atom.select(pdb,chain=uchains[j],elety="CA",verbose=FALSE)$atom,"resid"]),sep="",collapse="");	
				#Take the coordinate from seqres
				# seqres1=paste(aa321(pdb$seqres[which(names(pdb$seqres)==chains[i])]),sep="",collapse="")
				# seqres2=paste(aa321(pdb$seqres[which(names(pdb$seqres)==chains[j])]),sep="",collapse="")
				
				if(seqres1==seqres2)
				{
				#print(c(chains[i],chains[j],"diff"));
					flag=TRUE;
					break;
				}
				
				
			}
			if(flag==FALSE)
			{
				uchains=c(uchains,chains[i])		
			}
	}
		return(uchains);	
}

      pdb.bounds <- function(nums) {
        nums <- as.numeric(nums)
        bounds <- nums[1]
        diff.i <- 1
        j <- 1
        nums.start <- nums[1]
        store.inds <- NULL
        for (i in 2:length(nums)) {
            if (nums[i] != nums[j]) {
                if ((nums[i] - diff.i) != nums.start) {
                  bounds <- c(bounds, nums[i - 1], nums[i])
                  nums.start <- nums[i]
                  diff.i <- 1
                }
                else {
                  diff.i <- diff.i + 1
                }
                store.inds <- c(store.inds, i)
            }
            j <- j + 1
        }
        bounds <- c(bounds, nums[length(nums)])
        bounds <- matrix(bounds, ncol = 2, byrow = TRUE, dimnames = list(NULL, 
            c("start", "end")))
        bounds <- cbind(bounds, length = (bounds[, 2] - bounds[, 
            1]) + 1)
        return(list(bounds = bounds, r.ind = store.inds))
    }    
pdbsum=function(pdb,filepath,doprint=TRUE,append=FALSE)
{
            sum.segid <- unique(pdb$atom[, "segid"])
            sum.chain <- matrix(unique(pdb$atom[, "chain"]),nrow=1);
            sum.rnum <- pdb.bounds(pdb$atom[, "resno"])
            sum.resno <- sum.rnum$bounds;
            
            sum.resid <- table(pdb$atom[sum.rnum$r.ind, "resid"])
            sum.eleno <- pdb.bounds(pdb$atom[, "eleno"])$bounds
            sum.elety <- table(pdb$atom[, "elety"])
            if(doprint)
            {
            cat(" * Structure Summary *", sep = "\n", file=filepath,append=append)
            cat("---- segid ----", sep = "\n", file=filepath,append=TRUE)
            cat(sum.segid,sep="\t",eol="\n",file=filepath,append=TRUE)
            cat("\n---- chain ----", sep = "\n", file=filepath,append=TRUE)
            write.table(sum.chain,sep="\t",eol="\n",quote=FALSE, col.names=FALSE,row.names=FALSE, file=filepath,append=TRUE)
            cat("\n---- resno ----", sep = "\n", file=filepath,append=TRUE)
            write.table(sum.resno,sep="\t",eol="\n",row.names=FALSE,quote=FALSE, file=filepath,append=TRUE)
            cat("\n---- resid ----", sep = "\n", file=filepath,append=TRUE)
            write.table(sum.resid,sep="\t",eol="\n",quote=FALSE,row.names=FALSE, file=filepath,append=TRUE)
            cat("\n---- eleno ----", sep = "\n", file=filepath,append=TRUE)
            write.table(sum.eleno,sep="\t",eol="\n",quote=FALSE, row.names=FALSE,file=filepath,append=TRUE)
            cat("\n---- elety ----", sep = "\n", file=filepath,append=TRUE)
            write.table(sum.elety,sep="\t",eol="\n",quote=FALSE,row.names=FALSE, file=filepath,append=TRUE)
	    }
            return(list(sum.segid=sum.segid,sum.chain=sum.chain,sum.resno=sum.resno,sum.resid=sum.resid,sum.eleno=sum.eleno,sum.elety=sum.elety));
}

print_HTMLglobal_sum_start=function(filepath)
{
  write2Html("<html>\n",filepath);
  write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
  tipvec=c("Filename","The chains considered for the analysis","Total number of residues",paste("Total number of residue of type: ",c("Aromatic (FWY)","Uncharged Polar (NQST)","Positively charged (RKH)","Negatively charged (DE)"),sep=""),
           "Arg/Lys ratio",paste("Total number interactions: ",c("Ionpair","Aromatic Aromatic interaction","Aromatic sulphur interaction","Cation pi interation","Disulfide bond","Hydrogen bond","Hydrophobic interactions"),sep=""),           
           "Total number of Proline residues","Total number of proline at 2nd position of beta turns","Total number of prolines at N-cap of helix",
           "Ratio of Nonpolar to Polar accessible surface areas","Total number of thermolabile bonds","Total number of dipole stabilised helices","Total number of conformationally strained residues","Total number of metal","Total energy, an approximation of stability");
  
  write2Html.tableHeader(c("Filename","Chain","Length","Aro","NQST","RKH","DE","R/K","TotalIP","TotalAAI","TotalASI","TotalCPI","TotalDB","TotalHB","TotalHP","TotalPro","TotalBt2P","TotalNcapP","NP/P","TotalTL","St.Helix","TotalCS","TotalMetal","TotalGibbsEnergy"),filepath,append=TRUE,tipvec);

}
