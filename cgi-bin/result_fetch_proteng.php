<?php
$rand=$_GET["ID"];
$subpref1="HTML";
$status = 0;	
$chkStatus = "refresh";
$statusfile="http://irdp.ncl.res.in/out/".$rand."/".$subpref1."/status.txt";
sscanf(file_get_contents($statusfile),"%d",$status);
if($status == 100 or $status == -111) $chkStatus = "";
?>

<!DOCTYPE html>
<html lang="en">
<head>
<title>iRDP | iStability</title>
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
          <li><a href="../CompStabAna.html">iCAPS</a></li>
          <li><a href="../Pattern_analysis.html">iMotifs</a></li>
          <!-- <li><a class="active" href="Mutation_analysis.html">iMutants</a></li> -->
          <li><a href="../Mutation_analysis.html">iMutants</a></li>
          <li><a class="active" href="../Proteng.html">iStability</a></li>
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

function printSummary($name)
{
	$temp= $_SERVER['DOCUMENT_ROOT']."$name";
	$html=file_get_contents($temp);
	echo("$html");
}

function printLog($name)
{
	echo "<img src='../ajax-loader_circle_ebebeb.gif' alt='some_text'><br><br>";
	echo "<b>Your job is running<br>This page will be automatically updated every 10 seconds. You can also reload it manually.<br>"; 
	echo "If you wish to view these results at a later time without recalculating them, please bookmark this page. <br>The results will be kept on the server for 7 days<br><br>";
	$temp= $_SERVER['DOCUMENT_ROOT']."$name";
	$html=file_get_contents($temp);
	print("<textarea cols='50' rows='20' name='textdesc'>$html</textarea><br><br>");
}

function help($name)
{
		$temp= "http://irdp.ncl.res.in/"."$name";
		echo("<a href=\"".$temp."\" target=\"_blank\">HELP</a> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp");
}

if($status != 100 and $status != -111)
	{
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
	

printSummary("/out/$rand/$subpref1/result.html");	
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

