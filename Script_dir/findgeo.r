print_HTMLmetal=function(fgeofinal,fgeocode,filepath,append=FALSE,html=TRUE)
{
	if(html){write2Html("<html>\n",filepath,append);}else{write2Html("",filepath,append);}
	write2Html("<br><br><b><u>Metal geometry and binding site details</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("Metal name","Chain","Residue number","Atom number","Metal geometry [Findgeo]","Coordinating residues","Geometry figure");
	write2Html.tableHeader(c("Metal","Chain","Res.No","Atom.No","Geometry","Site_Details","Figure"),filepath,append=TRUE,tipvec);
	fgeofinal[,6]=apply(matrix(fgeofinal[,6],ncol=1),1,function(x){gsub(",","<br>",x)}); #Replace , with <br>
	

	fgeotemp=apply(matrix(fgeocode,ncol=1),1,function(x){paste("<img src='","/icons/findgeo/",x,".png' alt='some_text'  height='96' width='132'>",sep="",collapse="")});
	
	for(i in 1:nrow(fgeofinal))
	{
		write2Html("<tr>",filepath,TRUE);
		write2Html.tableBodyOnlytd(matrix(fgeofinal[i,],ncol=1),filepath,TRUE,"<td>");
		write2Html("<td>",filepath,TRUE);
		write2Html(fgeotemp[i],filepath,TRUE);
		write2Html("</td>",filepath,TRUE);
		write2Html("</tr>",filepath,TRUE);
	}
	write2Html("\n</table>\n",filepath,TRUE);
	if(html){write2Html("</html>\n",filepath,TRUE);}else{write2Html("",filepath,TRUE);}
}


print_HTMLmetal_sum_start=function(filepath)
{
	write2Html("<html>\n",filepath);
	write2Html("<b><u>Metal binding summary</u><br><br></b>\n",filepath,TRUE);
	write2Html("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \">\n",filepath,TRUE);
	tipvec=c("PDB Filename","Total number of metails","Format: Metal name- Count");
	write2Html.tableHeader(c("Filename","Total","Details"),filepath,append=TRUE,tipvec);
}

