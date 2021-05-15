<!DOCTYPE html>
<html lang="en">
<head>
<title>iRDP | iStability</title>
<meta charset="utf-8">
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
<body id="page17">
<!--==============================header=================================-->
<header>
  <div class="border-bot">
    <div class="main">
      <h1><a href="../index.html">ReDesign</a></h1>
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
      <div class="wrapper">
        <article class="grid_12">
	
	
	<?php
	
	function checknumeric($value,$error1,$error2)
	{
		if(!(is_numeric($value))){
		exit($error1);
		}
		else if(floatval($value)<=0)
		{
		exit($error2);
		}
	}
	$pdbmode="-";
	$uploadmode="-";
	$rand=$_POST["randID"]; #generate random number
	$subpref1="HTML";
	$pdbids=$_POST["pdbid"]; #pdbids comma separated
	$optionFile='../option/'.$rand."_option.txt";
	$mutFile='../option/'.$rand."_mut.txt";
	
	if($_POST["pdbids"]=="pdbids")
	{
	#user has selected to go by pdb ids rather than file upload
		if($pdbids=="")
		{
		exit("pdb id cant be blank");
		}else{
		if(strlen($pdbids)!=4){exit("Invalid pdb id");}
		$pdbmode=$pdbids;
		$uploadmode="-";
		}
	}
	else{
	#user has selected to go by file uplaod rather than pdb ids
	$pdbmode="-";
	$uploadmode=$rand;
	}
	
	$conserv=0;
	if($_POST["conserv"]=="conserv"){$conserv=1;}
	
	$opt=0;
	$mut="";
	$am=$_POST["AM"];
	$imph=$_POST["IMpH"];
	$imtemp=$_POST["IMtemp"];
	
	$foldph=$_POST["FoldpH"];
	$foldtemp=$_POST["Foldtemp"];
	$foldion=$_POST["Foldion"];
	$tool=$_POST["tool"];
	$sstop=$_POST["sstop"];
	if($_POST["rule"]=="ssbond"){$opt=1;}
	if($_POST["rule"]=="bt2p"){$opt=2;}
	if($_POST["rule"]=="Ncapp"){$opt=3;}
	if($_POST["rule"]=="Confstr"){$opt=4;}
	
	if($_POST["rule"]=="Custom")
	{
		$opt=5;
		$mut=$_POST["mutation"];
		if(strlen($mut)==0)
		{
		exit("Mutation field cant be blank");
		}else{
		$fm = fopen($mutFile, 'w') or die('fopen mut failed');
		fwrite($fm,$_POST["mutation"]);
		fclose($fm);
		}
	}else{
	$mutFile="";
	}

	checknumeric($_POST["sstop"],"############################## <br><br> Number of SSBOND to predict for stability must be numeric","############################## <br><br> Number of SSBOND to predict for stability must be greater than 0");
	
	$fp = fopen($optionFile, 'w') or die('fopen failed');
	fwrite($fp, "prefix\t$rand\npdbmode\t$pdbmode\nuploadmode\t$uploadmode\noption\t$opt\nAM\t$am\nimph\t$imph\nimtemp\t$imtemp\nfoldph\t$foldph\nfoldtemp\t$foldtemp\nfoldion\t$foldion\ntool\t$tool\nsstop\t$sstop\nconserv\t".$conserv."\n"
	);
	
	fclose($fp);
	
	$url='/usr/bin/Rscript /var/www/Script_dir/ProTEng_V5.r ' . $optionFile .' '.$mutFile. ' >/dev/null &';
	#$url='/usr/bin/Rscript /var/www/Script_dir/ProTEng_V5.r ' . $optionFile .' '.$mutFile;
	exec($url,$x);
	
	echo "<h3>Congratulations!</h3>";
	echo '<h6 class="prev-indent-bot" align="justify">Your job is submmited to the server.</h6>';
	echo '<p class="indent-bot">Please <a target="_top" href="http://irdp.ncl.res.in/cgi-bin/result_fetch_proteng.php?ID='  . $rand . '"> <b>click here</b></a> to check the result status.</p><br>';

#			for($i=0;$i<count($x);$i++)
#	 		{
#			  echo $x[$i];
#			  echo "<br>";
#			}
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

