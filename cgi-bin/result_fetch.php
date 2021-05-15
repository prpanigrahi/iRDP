<?php
$rand=$_GET["ID"];
$subpref1="HTML";
$status = 0;	
$chkStatus = "refresh";
$statusfile="../html/out/".$rand."/".$subpref1."/status.txt";
sscanf(file_get_contents($statusfile),"%d",$status);
if($status == 100 or $status == -111) $chkStatus = "";
?>


<!DOCTYPE html>
<html lang="en">
<head>
<title>iRDP | iCAPS</title>
<meta charset="utf-8">
<meta http-equiv="<?php echo $chkStatus ?>" content=" 10;">
<link rel="stylesheet" href="../css/reset.css" type="text/css" media="screen">
<link rel="stylesheet" href="../css/style.css" type="text/css" media="screen">
<link rel="stylesheet" href="../css/grid.css" type="text/css" media="screen">
<script src="../js/jquery-1.7.1.min.js"></script>
<script src="../js/cufon-yui.js"></script>
<script src="../js/cufon-replace.js"></script>
<script src="../js/Vegur_500.font.js"></script>
<script src="../js/Ropa_Sans_400.font.js"></script>
<script src="../js/FF-cash.js"></script>
<script src="../js/script.js"></script>
<!--[if lt IE 9]>
<script src="js/html5.js"></script>
<link rel="stylesheet" href="css/ie.css" type="text/css" media="screen">
<![endif]-->
</head>
<body id="page3">
<!--==============================header=================================-->
<header>
  <div class="border-bot">
    <div class="main">
      <h1><a href="../index.html">iRDP</a></h1>
      <nav>
        <ul class="menu">
          <li><a href="../index.html">Home</a></li>
          <li><a class="active" href="../CompStabAna.html">iCAPS</a></li>
          <li><a href="../Pattern_analysis.html">iMotifs</a></li>
          <!-- <li><a class="active" href="Mutation_analysis.html">iMutants</a></li> -->
          <li><a href="../Mutation_analysis.html">iMutants</a></li>
          <li><a href="../Proteng.html">iStability</a></li>
          <li><a href="../iATMs.html">iATMs</a></li>
          <li><a href="../help.html">Help</a></li>
          <li><a href="../contact.html">Contact Us</a></li>
        </ul>
      </nav>
      <div class="clear"></div>
    </div>
  </div>
</header>
<!--==============================content================================-->
<section id="content">
  <div class="main">
    <div class="container_12">
      <div >
        <article >	
        


	<?php

	function write2Html_tableBody($col1,$link)
		{
			# $rand,$dir,$sum (args)
		#echo "<tr> <td>$col1 </td> <td> <a href='/out/$rand/$dir/$sum' target=\"_blank\" > Result </a> </td> </tr>";
		echo("<tr> <td><b> $col1 </b> </td> <td> <a href='".$link."'> Result </a> </td> </tr>");
		}
	
		function printSummary($name)
		{
		$temp= $_SERVER['DOCUMENT_ROOT']."$name";
		$html=file_get_contents($temp);
		echo("$html");
		}
	
		function gototop()
		{
		echo('<a href="#top">TOP   </a> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp');
		}
		
		function help($name)
		{
		$temp= "http://irdp.ncl.res.in/"."$name";
		echo("<a href=\"".$temp."\" target=\"_blank\">HELP</a> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp");
		}
	
	
	$handle = fopen('../option/'.$rand."_option.txt", "r");
	$cntLn=0;
	$cntA=0;
	$selectLine=array(4,5,6,7,8,10,15,17,19,20,22,23,24,25,26,27,28);
	$arrayContent=array(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	if ($handle) {
	    while (($line = fgets($handle)) !== false) {
		if(($selectLine[$cntLn]-1) == $cntA)
		{
			sscanf($line,"%s\t%d",$tmp,$arrayContent[$cntLn]);	
			$cntLn++;
		}
		$cntA++;
	    }
	} else {
	    echo "Error opening options file.";
	    exit("Error opening options file.");
	}


	#$htmlFile='/out/'.$rand."/".$rand."_result.html";
		# Results starts here	
	

	echo('<a id=top> </a>');
	#echo('hii');
	#echo('bye');


	echo "<br>";
	#echo $chkStatus;	
	
	
	function printLog($name)
	{
		echo "<br>";
		echo "<img src='../ajax-loader_circle_ebebeb.gif' alt='some_text'><br><br>";
		echo "<b>Your job is running<br>This page will be automatically updated every 10 seconds. You can also reload it manually.<br>"; 
		echo "If you wish to view these results at a later time without recalculating them, please bookmark this page. <br>The results will be kept on the server for 7 days<br><br>";
		$temp= "http://irdp.ncl.res.in/$name";
		$html=file_get_contents($temp);
		print("<textarea cols='50' rows='20' name='textdesc'>$html</textarea><br><br>");
	}

	if($status != 100 and $status != -111)
	{
	echo "<h3>Status</h3><h6><b>".$status."% of files are processed.</b></h6>";
	printLog("/out/$rand/$subpref1/log.txt");
	}else{
	
	if($status != -111)
	{
		print("<br><b>Your job is finished. Thank you for submitting </b><br>");
		echo "<a href="."/out/$rand/$subpref1/log.txt"." target=\"_blank\" >Check Log file</a><br><br><br>";
	}else{
		print("<br><b>We are extremely sorry !! We encountered some error </b><br> Please <a href=\"/contact.html\"> contact </a> the author. <br> Thank you for submitting<br>");
		echo "<a href="."/out/$rand/$subpref1/log.txt"." target=\"_blank\" >Check Log file</a><br><br><br>";
	}

	}

	#echo "$statusfile";


		echo("<h3>Results available for</h3><br>");
		echo("<table border=1 cellspacing=0 cellpadding=5 style=\"table-layout:fixed; \" >");
	#echo ""<tr style=\"background-color:pink\"><th colspan=4 title=\"First\"> Donor </th> <th colspan=4> Acceptor </th> <th rowspan=2> Type </th><th rowspan=2> DAdist </th> <th rowspan=2> HAdist </th> <th rowspan=2> DHAang </th> <th rowspan=2> HAAAang </th> <th rowspan=2> DAAAang </th> </tr>";
		echo("<tr style=\"background-color:rgb(235,235,235)\" > <th> Parameter </th> <th> Result </th> </tr>");
	
	write2Html_tableBody("Summary","#summary");
	if($arrayContent[0]==1){write2Html_tableBody("Aminoacid composition","#aatag");}
	if($arrayContent[1]==1){write2Html_tableBody("Secondary structure composition","#sstag");}
	if($arrayContent[2]==1){write2Html_tableBody("Amino acid distribution in  Helices","#hcomptag");}
	if($arrayContent[2]==1){write2Html_tableBody("Amino acid distribution in  Strands","#scomptag");}
	if($arrayContent[2]==1){write2Html_tableBody("Amino acid distribution in  Turns","#tcomptag");}
	if($arrayContent[2]==1){write2Html_tableBody("Amino acid distribution in  Coils","#ccomptag");}
	if($arrayContent[4]==1){write2Html_tableBody("Ion pairs","#ionsum");}
	if($arrayContent[5]==1){write2Html_tableBody("Aromatic-Aromatic interaction","#arosum");}
	if($arrayContent[6]==1){write2Html_tableBody("Aromatic-Sulphur interaction","#arossum");}

	if($arrayContent[7]==1){write2Html_tableBody("Cation-pi interaction","#catpitag");}
	if($arrayContent[8]==1){write2Html_tableBody("Hydrogen bonds","#hbtag");}
	if($arrayContent[9]==1){write2Html_tableBody("Hydrophobic interaction","#hptag");}
	if($arrayContent[3]==1){write2Html_tableBody("Disulfide bonds","#distag");}
	if($arrayContent[13]==1){write2Html_tableBody("Proline profile","#protag");}
	if($arrayContent[10]==1){write2Html_tableBody("ASA summary","#asatag");}

	if($arrayContent[12]==1){write2Html_tableBody("Thermolabile bond profile","#deamtag");}
	if($arrayContent[11]==1){write2Html_tableBody("Helix dipole stabilization profile","#hdipotag");}
	if($arrayContent[14]==1){write2Html_tableBody("Conformationally strained residue profile","#lefttag");}
	if($arrayContent[15]==1){write2Html_tableBody("Metal binding sites","#metal");}
	if($arrayContent[16]==1){write2Html_tableBody("Gibbs energy of folding","#gibbs");}
	
		
	
		echo("</table>");
	
		echo("<br><br>");
		if($status == 100) echo("<h3> Download Results </h3> Results are available to <a href=\"/out/$rand/result.zip\"> download from here. </a>");
		echo("<br><br> <h3> Results </h3>");
		#echo("</html>");
		#echo "<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>";
	
		# Load aacomp
		if($arrayContent[0]==1)
		{
		echo('<a id="aatag"> </a> <br>');
		gototop();
		help("help_icaps.html#aacomp");
		printSummary("/out/$rand/$subpref1/AminoComp/aacomp_summary.html");
		if($status!=100){echo "</table></html>";}
		}

		if($arrayContent[1]==1)
		{
		echo('<a id="sstag"> </a> <br>');
		gototop();
		help("help_icaps.html#sscomp");
		printSummary("/out/$rand/$subpref1/SecStrCompo/secstr_summary.html");
		if($status!=100){echo "</table></html>";}
		}
	
		if($arrayContent[2]==1)
		{
				echo('<a id="hcomptag"> </a> <br>');
				gototop();
				help("help_icaps.html#hstccomp");
				printSummary("/out/$rand/$subpref1/SecStrDistr/hcomp_summary.html");
				if($status!=100){echo "</table></html>";}
				echo('<a id="scomptag"> </a> <br>');
				gototop();
				help("help_icaps.html#hstccomp");
				printSummary("/out/$rand/$subpref1/SecStrDistr/scomp_summary.html");
				if($status!=100){echo "</table></html>";}
				echo('<a id="tcomptag"> </a> <br>');
				gototop();
				help("help_icaps.html#hstccomp");
				printSummary("/out/$rand/$subpref1/SecStrDistr/tcomp_summary.html");
				if($status!=100){echo "</table></html>";}
				echo('<a id="ccomptag"> </a> <br>');
				gototop();
				help("help_icaps.html#hstccomp");
				printSummary("/out/$rand/$subpref1/SecStrDistr/ccomp_summary.html");
				if($status!=100){echo "</table></html>";}
		}
	
		# Load ionic summary
		if($arrayContent[4]==1)
		{
		echo('<a id="ionsum"> </a> <br>');
		gototop();
		help("help_icaps.html#ip");
		printSummary("/out/$rand/$subpref1/Ionpair/Ionpair_summary.html");
		if($status!=100){echo "</table></html>";}
		}
		if($arrayContent[5]==1)
		{
		# Load aromatic summary
		echo('<a id="arosum"> </a> <br>');
		gototop();
		help("help_icaps.html#aai");
		printSummary("/out/$rand/$subpref1/AroAro/Aromatic_summary.html");
		if($status!=100){echo "</table></html>";}
		}
	
		if($arrayContent[6]==1)
		{
		# Load aromatic summary
		echo('<a id="arossum"> </a> <br>');
		gototop();
		help("help_icaps.html#asi");
		printSummary("/out/$rand/$subpref1/AroS/AroS_summary.html");
		if($status!=100){echo "</table></html>";}
		}
	
		if($arrayContent[7]==1)
		{
		# Load aromatic summary
		echo('<a id="catpitag"> </a> <br>');
		gototop();
		help("help_icaps.html#cpi");
		printSummary("/out/$rand/$subpref1/Catpi/Catpi_summary.html");
		if($status!=100){echo "</table></html>";}
		}
	
		if($arrayContent[3]==1)
		{
		# Load aromatic summary
		echo('<a id="distag"> </a> <br>');
		gototop();
		help("help_icaps.html#dis");
		printSummary("/out/$rand/$subpref1/Disulfide/Disulfide_summary.html");
		if($status!=100){echo "</table></html>";}
		}
	
		if($arrayContent[8]==1)
		{
		# Load aromatic summary
		echo('<a id="hbtag"> </a> <br>');
		gototop();
		help("help_icaps.html#hb");
		printSummary("/out/$rand/$subpref1/Hbond/Hbond_summary.html");
		if($status!=100){echo "</table></html>";}
		}
	
		if($arrayContent[9]==1)
		{
		# Load aromatic summary
		echo('<a id="hptag"> </a> <br>');
		gototop();
		help("help_icaps.html#hp");
		printSummary("/out/$rand/$subpref1/Hydrophobic/Hydrophobic_summary.html");
		if($status!=100){echo "</table></html>";}
		}
	
		if($arrayContent[13]==1)
		{
		# Load aromatic summary
		echo('<a id="protag"> </a> <br>');
		gototop();
		help("help_icaps.html#pp");
		printSummary("/out/$rand/$subpref1/Proline/Proline_summary.html");
		if($status!=100){echo "</table></html>";}
		}
		if($arrayContent[10]==1)
		{
		# Load aromatic summary
		echo('<a id="asatag"> </a> <br>');
		gototop();
		help("help_icaps.html#asa");
		printSummary("/out/$rand/$subpref1/ASA/ASA_summary.html");
		if($status!=100){echo "</table></html>";}
		}
	
		if($arrayContent[12]==1)
		{
		# Load aromatic summary
		echo('<a id="deamtag"> </a> <br>');
		gototop();
		help("help_icaps.html#tbp");
		printSummary("/out/$rand/$subpref1/Deamidation/Deamidation_summary.html");
		if($status!=100){echo "</table></html>";}
		}
	
		if($arrayContent[11]==1)
		{
		# Load aromatic summary
		echo('<a id="hdipotag"> </a> <br>');
		gototop();
		help("help_icaps.html#hdsp");
		printSummary("/out/$rand/$subpref1/HelixDipole/Hdipole_summary.html");
		if($status!=100){echo "</table></html>";}
		}

		if($arrayContent[14]==1)
		{
		# Load aromatic summary
		echo('<a id="lefttag"> </a> <br>');
		gototop();
		help("help_icaps.html#csrp");
		printSummary("/out/$rand/$subpref1/ConfStrain/ConfStrain_summary.html");
		if($status!=100){echo "</table></html>";}
		}
	
		if($arrayContent[15]==1)
		{
		# Load metal summary
		echo('<a id="metal"> </a> <br>');
		gototop();
		help("help_icaps.html#mba");
		printSummary("/out/$rand/$subpref1/Metal/Metal_summary.html");
		if($status!=100){echo "</table></html>";}
		}
		if($arrayContent[16]==1)
		{
		# Load metal summary
		echo('<a id="gibbs"> </a> <br>');
		gototop();
		help("help_icaps.html#gibbs");
		printSummary("/out/$rand/$subpref1/Gibbs/Gibbs_summary.html");
		if($status!=100){echo "</table></html>";}
		}
		
		echo('<a id="summary"> </a> <br>');
		gototop(); echo('<b><u>Summary</b></u>');
		printSummary("/out/$rand/$subpref1/Summary.html");
		if($status!=100){echo "</table></html>";}
		
	?>

	</article>
        
      </div>
    </div>
  </div>
</section>
<!--==============================footer=================================-->
<footer>
  <div class="main">
    <div class="container_12">
      <div class="wrapper">
        <div class="grid_5">
          <div align="center" class="spacer-1"> <a href="../index.html"><img src="../images/csir.png" alt=""></a> </div>
        </div>
        <div class="grid_5">
          <div class="indent-top2">
            <br><p class="prev-indent-bot">&copy; 2014 All rights reserved by CSIR-NCL, Pune <a target="_blank" href="http://www.ncl-india.org/">ncl-india.org</a></p>
            Phone: +91 20 2590 2711 Email: <a href="mailto://irdpncl@gmail.com">irdpncl@gmail.com</a> </div>
        </div>
        
    </div>
  </div>
</footer>
<script>Cufon.now();</script>
</body>
</html>

