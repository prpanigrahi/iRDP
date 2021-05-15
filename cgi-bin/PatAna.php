<!DOCTYPE html>
<html lang="en">
<head>
<title>iRDP | iMotifs</title>
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
<body id="page15">
<!--==============================header=================================-->
<header>
  <div class="border-bot">
    <div class="main">
      <h1><a href="../index.html">ReDesign</a></h1>
      <nav>
        <ul class="menu">
          <li><a href="../index.html">Home</a></li>
          <li><a href="../CompStabAna.html">iCAPS</a></li>
          <li><a class="active" href="../Pattern_analysis.html">iMotifs</a></li>
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
      <div class="wrapper">
        <article class="grid_12">


	<?php
	#echo error_reporting(E_ALL);
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
	$pdbidFile='../option/'.$rand."_pdbid.txt";

	
	$dihed=0;$uniq=0;
	if($_POST["uniq"]=="uniq"){$uniq=1;}
	if($_POST["dihed"]=="dihed"){$dihed=1;}
	
	if($_POST["pdbids"]=="pdbids")
	{
	#user has selected to go by pdb ids rather than file upload
		if($pdbids=="")
		{
		exit("pdb id cant be blank");
		}else{
		$fpdb = fopen($pdbidFile, 'w') or die('fopen failed for pdbid file'); #write the pdb ids into a file
		fwrite($fpdb,"$pdbids");
		fclose($fpdb);
		$pdbmode="args";
		$uploadmode="-";
		}
	
	}
	else{
	#user has selected to go by file uplaod rather than pdb ids
	$pdbmode="-";
	$pdbidFile="";
	$uploadmode=$rand;
	}
	
	
	if($_POST["pattern"]==""){exit("pattern cant be blank");}
	
	checknumeric($_POST["ioncut"],"############################## <br><br> Ion pair cutoff value must be numeric","############################## <br><br> Ion pair cutoff value must be greater than 0");
	checknumeric($_POST["aroarolow"],"############################## <br><br> Aromatic aromatic interaction lower cutoff value must be numeric","############################## <br><br> Aromatic aromatic interaction lower cutoff value must be greater than 0");		
	checknumeric($_POST["aroarohigh"],"############################## <br><br> Aromatic aromatic interaction higher cutoff value must be numeric","############################## <br><br> Aromatic aromatic interaction higher cutoff value must be greater than 0");		
	checknumeric($_POST["aroarodihed"],"############################## <br><br> Aromatic aromatic interaction dihedral cutoff value must be numeric","############################## <br><br> Aromatic aromatic interaction dihedral cutoff value must be greater than 0");		
	checknumeric($_POST["aroscut"],"############################## <br><br> Aromatic-Sulphor interaction distance cutoff value must be numeric","############################## <br><br> Aromatic-Sulphor interaction distance cutoff value must be greater than 0");
	checknumeric($_POST["catpicut"],"############################## <br><br> Cation-pi interaction distance cutoff value must be numeric","############################## <br><br> Cation-pi interaction distance cutoff value must be greater than 0");
	checknumeric($_POST["hpcut"],"############################## <br><br> Hydrophobic interaction distance cutoff value must be numeric","############################## <br><br> Hydrophobic interaction distance cutoff value must be greater than 0");								
	checknumeric($_POST["acclow"],"<br>Relative Solvent Accessible Area cutoff must be numeric<br>","<br>Relative Solvent Accessible Area cutoff must be positive<br>");
	
	if(floatval($_POST["aroarodihed"])>90.0)
	{
	exit("Dihedral cut-off must be acute angle <=90.0");
	}
	if(floatval($_POST["aroarohigh"])-floatval($_POST["aroarolow"])<=0)
	{
	exit("Aromatic aromatic centroid distance range invalid");
	}
	
	#write option file
	$fp = fopen($optionFile, 'w') or die('fopen failed');
	fwrite($fp, "prefix\t$rand\npdbmode\t$pdbmode\nuploadmode\t$uploadmode\npattern\t".$_POST["pattern"]."\nioncut\t".$_POST["ioncut"]."\naroarolow\t".$_POST["aroarolow"]."\naroarohigh\t".$_POST["aroarohigh"]."\nusedihed\t$dihed\naroarodihed\t"
	.$_POST["aroarodihed"]."\naroscut\t".$_POST["aroscut"]."\ncatpicut\t".$_POST["catpicut"]."\nhpcut\t".$_POST["hpcut"]."\nuniq\t$uniq\nacclow\t".$_POST["acclow"]."\n");
	fclose($fp);
	
	#$url='/usr/bin/Rscript /var/www/Script_dir/DXR_analysis_V8.r ' . $optionFile.' '.$pdbidFile. ' >/dev/null &';
	#$url='/usr/bin/Rscript /var/www/Script_dir/DXR_analysis_V8.r ' . $optionFile.' '.$pdbidFile;
	
	$url='/usr/bin/Rscript /var/www/Script_dir/Pattern_analysis_V11.r ' . $optionFile.' '.$pdbidFile. ' >/dev/null &';
#	$url='/usr/bin/Rscript /var/www/Script_dir/Pattern_analysis_V11.r ' . $optionFile.' '.$pdbidFile;
	exec($url,$x);
	
	echo "<h3>Congratulations!</h3>";
	echo '<h6 class="prev-indent-bot" align="justify">Your job is submmited to the server.</h6>';
	echo '<p class="indent-bot">Please <a target="_top" href="http://irdp.ncl.res.in/cgi-bin/result_fetch_PatAna.php?ID=' . $rand . '"> <b>click here</b></a> to check the result status.</p><br>';
	#echo "Your job is submmited on server.<br> Please click OK to check the result status.";
	#echo('<a target="_top" href="http://172.28.0.171/cgi-bin/result_fetch_PatAna.php?ID=' . $rand . '"> <b>Check Results</b> </a>');
	#echo('');

			for($i=0;$i<count($x);$i++)
			{
			  echo $x[$i];
			  echo "<br>";
			}
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

