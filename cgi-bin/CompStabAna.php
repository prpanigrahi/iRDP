<!DOCTYPE html>
<html lang="en">
<head>
<title>iRDP | iCAPS</title>
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
<body id="page14">
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
      <div class="wrapper">
        <article class="grid_12">

	<?php
	#echo error_reporting(E_ALL);
####### Note: $_POST["pdbids"] is for radio button (Its value can be either pdbids[if mode is pdbid mode] or upload[mode is via file upload] )
#######       $_POST["pdbid"] is for text area where user iputs pdbids comma separated

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
	$aa=0;$ss=0;$aass=0;$dis=0;$ip=0;$ap=0;$dihed=0;$aros=0;$catpi=0;$hp=0;$hb=0;$asa=0;$hdipo=0;$deam=0;$pro=0;$left=0;$metal=0;$stab=0;$uniq=0;
	$pdbmode="-";
	$uploadmode="-";
	$subpref1="HTML";
	$rand=$_POST["randID"]; #generate random number
	#$rand="1".$rand;
	#$rand="50000";
	$pdbids=$_POST["pdbid"]; #pdbids comma separated
	$optionFile='../option/'.$rand."_option.txt";
	$pdbidFile='../option/'.$rand."_pdbid.txt";

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
	
	if($_POST["aa"]=="aa"){$aa=1;}
	if($_POST["ss"]=="ss"){$ss=1;}
	if($_POST["aass"]=="aass"){$aass=1;}
	if($_POST["dis"]=="dis"){$dis=1;}
	if($_POST["ip"]=="ip"){$ip=1;}
	if($_POST["ap"]=="ap"){$ap=1;}
	if($_POST["dihed"]=="dihed"){$dihed=1;}
	if($_POST["aros"]=="aros"){$aros=1;}
	if($_POST["catpi"]=="catpi"){$catpi=1;}
	if($_POST["hp"]=="hp"){$hp=1;}
	if($_POST["hb"]=="hb"){$hb=1;}
	if($_POST["asa"]=="asa"){$asa=1;}
	if($_POST["hdipo"]=="hdipo"){$hdipo=1;}
	if($_POST["deam"]=="deam"){$deam=1;}
	if($_POST["pro"]=="pro"){$pro=1;}
	if($_POST["left"]=="left"){$left=1;}
	if($_POST["metal"]=="metal"){$metal=1;}
	if($_POST["stab"]=="stab"){$stab=1;}
	if($_POST["uniq"]=="uniq"){$uniq=1;}
	
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
	
	if(($aa+$ss+$aass+$dis+$ip+$ap+$aros+$catpi+$hb+$hp+$asa+$hdipo+$deam+$pro+$left+$metal+$stab)==0)
	{
	exit("Choose at least one feature to analyze");
	}
	
	#write option file
	$fp = fopen($optionFile, 'w') or die('fopen failed');
	fwrite($fp, "prefix\t$rand\npdbmode\t$pdbmode\nuploadmode\t$uploadmode\naa\t$aa\nss\t$ss\naass\t$aass\ndis\t$dis\nip\t$ip\nioncut\t".$_POST["ioncut"]."\nap\t$ap\naroarolow\t".$_POST["aroarolow"]."\naroarohigh\t".$_POST["aroarohigh"]."\nusedihed\t$dihed\naroarodihed\t"
	.$_POST["aroarodihed"]."\naros\t$aros\naroscut\t".$_POST["aroscut"]."\ncatpi\t$catpi\ncatpicut\t".$_POST["catpicut"]."\nhb\t$hb\nhp\t$hp\nhpcut\t".$_POST["hpcut"]."\nasa\t$asa\nhdipo\t$hdipo\ndeam\t$deam\npro\t$pro\nleft\t$left\nmetal\t$metal\nstab\t$stab\nuniq\t$uniq\nacclow\t".$_POST["acclow"]."\n"
	);
	fclose($fp);
	
	#call CompStabana.r script
	$url='/usr/bin/Rscript /var/www/Script_dir/CompStabAna_V7.r ' . $optionFile.' '.$pdbidFile. ' >/dev/null &';
	#$url='/usr/bin/Rscript /var/www/Script_dir/CompStabAna_V7.r ' . $optionFile.' '.$pdbidFile;
	exec($url,$x);
	
	echo "<h3>Congratulations!</h3>";
	echo '<h6 class="prev-indent-bot" align="justify">Your job is submmited to the server.</h6>';
	echo '<p class="indent-bot">Please <a target="_top" href="http://irdp.ncl.res.in/cgi-bin/result_fetch.php?ID=' . $rand . '"> <b>click here</b></a> to check the result status.</p><br>';
	
	
	
#	for($i=0;$i<count($x);$i++)
#	{
#	  echo $x[$i];
#	  echo "<br>";
#	}
		
	#delete option file
	# unlink($optionFile);	
	
#	echo "hii";
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

