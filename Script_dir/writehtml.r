
			
write2Html=function(x,filepath,append=FALSE)
{
 # a function to write any html tags to a file
	cat(x,sep="",eol="",file=filepath,append=append)
}

write2Html.tableBody=function(tempmat,filepath,append=TRUE,td)
{
# a function dedicated to write table body as html tag
# it writes as <tr> all td stuffs </tr>
#  so a complete row i.e tr and td will be printed...
# what is u want only td to print, then use write2Html.tableBodyOnlytd function
	apply(
		matrix(1:nrow(tempmat),ncol=1),
		1,
		function(ind)
		{
		temprow=paste("<tr>\n",paste("\t",td,tempmat[ind,],"</td>",collapse=" "),"\n</tr>\n",collapse="");	
		write2Html(temprow,filepath,append=append);
		} 
	);

}

write2Html.tableBodyOnlytd=function(tempmat,filepath,append=TRUE,td)
{
# a function dedicated to write table body as html tag
# it writes as only td part
	apply(
		matrix(1:nrow(tempmat),ncol=1),
		1,
		function(ind)
		{
		temprow=paste("\t",td,tempmat[ind,],"</td>",collapse=" ");	
		write2Html(temprow,filepath,append=append);
		} 
	);

}

write2Html.tableHeader=function(tempvec,filepath,append=TRUE,tipvec="")
{
# a function dedicated to write table header as html tag
# you pass a vector of same no of column as table body

# If tooltip is to be provided for the table headers
# you have to pass tooltip=TRUE and the tipvec vector (Vector of tooltips)
# NOTE: tipvec and tempvec should be of equal length so that 1 to 1 correspondance of tooltip will be generated

write2Html(paste("<tr>",paste("<th style=\"background-color:rgb(235,235,235)\" title=\"", tipvec,"\" >",tempvec, "</th>",sep="",collapse=" "),"</tr>",sep="",collapse=""),filepath,append=append);

}


