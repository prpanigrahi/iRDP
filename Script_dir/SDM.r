checkSDM=function(infile="",chain="",mutation="",resno="")
{
	# infile has to be the exact pdb file name
	tempsdm=system(paste("perl SDM.pl ", infile," ",chain," ",mutation," ",resno,sep=""),intern=TRUE);
	return(c(tempsdm,checkSDMclass(tempsdm)));
}


checkSDMclass=function(tempsdm)
{
	if(tempsdm>=-0.5 & tempsdm <= 0.5)
	{
	return("N");
	}else if(tempsdm>=-1.0 & tempsdm < -0.5){
		return("SD");
		}else if(tempsdm>0.5 & tempsdm<=1.0){
			return("SS");
			}else if(tempsdm>=-2.0 & tempsdm< -1){
				return("D");	
				}else if(tempsdm>1.0 & tempsdm<=2){
					return("S")
					}else{
					 return("HD");
						}
					
}			
			
		
	

