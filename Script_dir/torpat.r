torPat=function(tempstr,pdb,tors)
{
#tempstr: resno.insertion.chain
	tempinds=which(colnames(tors$coords[,1,])==tempstr);
	return(tors$tbl[tempinds,]);
}

patAddNaccess=function(tempchres,nacc,acclow=20)
{
	#print(pos);
	chres=paste(nacc$asa[,1],nacc$asa[,2],sep=""); #ChainResnoResid format
	index=which(chres==tempchres);
	if(length(index)>0)
	{
		if(nacc$asa[index,5]<=acclow)
		{
		return("B");
		}else{return("E");}
	}else{
	return("-");
	}
}

patAddDSSP=function(tempchres,DSSP)
{
DSSPChRes=paste(DSSP$cha,DSSP$res,sep="") #combine chain and resno of DSSP object
index=which(DSSPChRes==tempchres);
if(length(index)==0){return("-");}
return(DSSP$ss[index]);
}

patAddRama=function(tempchres,proch)
{
prochChRes=paste(proch[,2],proch[,3],sep="") #combine chain and resno of DSSP object
index=which(prochChRes==tempchres);
if(length(index)==0){return("-");}
return(proch[index,6]);
}


posfreq=function(finalout,col=4,freqfile="temp.seq")
{
tempseq=paste(paste(paste("\n>pat",1:length(finalout[,4]),"\n",sep="")),finalout[,4],sep="",collapse="");
cat(tempseq,file=freqfile);
tempali=read.fasta(freqfile);
unlink(freqfile);
tempent=entropy(tempali$ali);
return(tempent$freq);
}

print_HTMLionic_torpat_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<h3>Motif summary</h3>\n",filepath,TRUE);

	tipvec=c(
"Number of ion-pairs formed by motif residues",
"Number of ion-pair networks in which motif residues are involved in",
"Number of aromatic interactions formed by motif residues",
"Number of aromatic-aromatic interaction networks in which motif residues are involved in",
"Number of Aro-S interaction formed by motif residues",
"Number of Aro-S interaction networks in which motif residues are involved in",
"Number of hydrogen bonds formed by motif residues",
"Number of disulfide bonds formed by motif residues",
"Number of Cat-pi interaction formed by motif residues",
"Number of Cat-pi interaction networks in which motif residues are involved in",
"Number of Hydrophobic formed by motif residues");

write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
write2Html("<tr> <th rowspan=2 title=\"PDB file name\" > Filename </th> <th rowspan=2 title=\"Chain in which the motif is present\"> Chain </th>  <th rowspan=2 title=\"Start residue number of motif\"> Start </th> <th rowspan=2 title=\"The amino acid sequence of motif\"> Motif </th> <th rowspan=2 title=\"The solvent accessibility of motif residues B:Buried, E:Exposed\"> Acc </th><th rowspan=2 title=\"The secondary structure of motif residues (DSSP notation)\"> SS </th><th rowspan=2 title=\"Ramachandran plot region of motif residues (Procheck notation)\"> R.plot_region </th> <th colspan=11 title=\"Interaction profile of motif\"> Interaction profile </th> </tr>",filepath,append=TRUE);
write2Html.tableHeader(c("IP","IP.Net","AP","AP.Net","AS","AS.Net","HB","Disul","Cat-pi","Cat-pi.Net","Hphob"),filepath,append=TRUE,tipvec);
}

print_HTMLdihed=function(filepath,append=FALSE,html=TRUE)
{
	if(html){write2Html("<html>\n",filepath,append);}else{write2Html("",filepath,append);}
	write2Html("<b><u>Main chain and Side chain dihedral angles</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c(
	"Chain","Residue number","Three letter code of residue","The main chain dihedral angle: Phi","The main chain dihedral angle: Psi",
	"The side chain dihedral angle: Chi1","The side chain dihedral angle: Chi2","The side chain dihedral angle: Chi3","The side chain dihedral angle: Chi4","The side chain dihedral angle: Chi5"
	);
	
	write2Html.tableHeader(c("Chain","Res.No","Res.ID","Phi","Psi","Chi1","Chi2","Chi3","Chi4","Chi5"),filepath,append=TRUE,tipvec);
}

newtorsion.pdb=function (pdb) 
{
pdb$atom[which(is.na(pdb$atom[,"insert"])),"insert"]="";

    colpaste <- function(x, col.names = colnames(x)) {
        apply(x, 1, function(row) paste(row[col.names], collapse = "."))
    }
    getinds <- function(atoms, ref = atom.names) {
        sort(atom2xyz(charmatch(atoms, ref)))
    }
    repadd <- function(num, nrep = nres, toadd = nxyz) {
        c(num, rep(num, (nrep - 1)) + rep(cumsum(rep(toadd, (nrep - 
            1))), each = length(num)))
    }
    atom.data <- colpaste(pdb$atom, c("elety", "resno", "insert","chain"))
    atom.list <- matrix(unlist(strsplit(atom.data, "\\.")), ncol = 4, 
        byrow = TRUE)
    res.data <- colpaste(pdb$atom, c("resno", "insert","chain"))
    res.list <- unique(res.data)
    atom.names <- c("N", "CA", "C", "O", "CB", "*G", "*G1", "*G2", 
        "*D", "*D1", "*D2", "*E", "*E1", "*E2", "*Z", "NH1", 
        "NH2")
    atom.greek <- c("N", "CA", "C", "O", "CB", "G", "G1", "G2", 
        "D", "D1", "D2", "E", "E1", "E2", "Z", "*", "*")
    coords <- NULL
    blank <- matrix(NA, nrow = 3, ncol = length(atom.names))
    for (i in 1:length(res.list)) {
        res.blank <- blank
        res.ind <- which(res.list[i] == res.data)
        atoms.noh <- atom.list[res.ind, 1]
        atoms.noh[grep("H", atoms.noh)] = "H"
        blank.ind <- charmatch(atoms.noh, atom.names, nomatch = 0) + 
            charmatch(substr(atoms.noh, 2, 4), atom.greek, nomatch = 0)
        res.blank[, blank.ind[blank.ind != 0]] <- matrix(pdb$xyz[atom2xyz(res.ind[blank.ind != 
            0])], nrow = 3)
        coords <- cbind(coords, res.blank)
    }
    natm <- length(atom.names)
    nxyz <- 3 * natm
    nres <- length(coords)/(nxyz)
    dim(coords) <- c(3, natm, nres)
    dimnames(coords) = list(xyz = c("x", "y", "z"), atm = atom.names, 
        res = res.list)
    co <- c(coords)
    chi1 <- torsion.xyz(co[repadd(getinds(c("N", "CA", "CB", 
        "*G")))])
    chi11 <- torsion.xyz(co[repadd(getinds(c("N", "CA", "CB", 
        "*G1")))])
    chi2 <- torsion.xyz(co[repadd(getinds(c("CA", "CB", "*G", 
        "*D")))])
    chi21 <- torsion.xyz(co[repadd(getinds(c("CA", "CB", "*G", 
        "*D1")))])
    chi2.ILE <- torsion.xyz(co[repadd(getinds(c("CA", "CB", "*G1", 
        "*D1")))])
    chi3 <- torsion.xyz(co[repadd(getinds(c("CB", "*G", "*D", 
        "*E")))])
    chi31 <- torsion.xyz(co[repadd(getinds(c("CB", "*G", "*D", 
        "*E1")))])
    chi4 <- torsion.xyz(co[repadd(getinds(c("*G", "*D", "*E", 
        "*Z")))])
    chi51 <- torsion.xyz(co[repadd(getinds(c("*D", "*E", "*Z", 
        "NH1")))])
    omega <- torsion.xyz(co[repadd(c(4:9, 52:57))])
    alpha <- c(NA, torsion.xyz(co[repadd(c(4:6, 55:57, 106:108, 
        157:159))]))
    phi <- c(NA, torsion.xyz(co[repadd(c(7:9, 52:60))]))
    psi <- torsion.xyz(co[repadd(c(1:9, 52:54))])
    tor.collapse <- function(a1, a11) {
        a <- a1
        got.a11 <- !(is.na(a11))
        a[got.a11] <- a11[got.a11]
        return(a)
    }
    chi1.F <- tor.collapse(chi1, chi11)
    chi2.F <- tor.collapse(chi2, chi21)
    chi2.F <- tor.collapse(chi2.F, chi2.ILE)
    chi3.F <- tor.collapse(chi3, chi31)
    tbl = cbind(phi[-(nres + 1)], psi, chi1.F, chi2.F, chi3.F, 
        chi4, chi51)
    colnames(tbl) <- c("phi", "psi", "chi1", "chi2", "chi3", 
        "chi4", "chi5")
    out <- list(psi = psi, phi = phi[-(nres + 1)], omega = omega, 
        chi1 = chi1.F, chi2 = chi2.F, chi3 = chi3.F, chi4 = chi4, 
        chi5 = chi51, alpha = alpha[-(nres + 1)], coords = coords, 
        tbl = tbl)
}

