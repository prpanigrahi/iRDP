# pdb summary

pdbsum=function(pdb)
{
	chains=unique(pdb$atom[,"chain"]); # Total chains
	dnach=screenDNA(pdb);  # Total DNA chains
	uniqch=uniqueChain(pdb); # Total unique chains

	# Total CA residues with coordinates
	apply(matrix(uniqch))
	
	
	
}	
