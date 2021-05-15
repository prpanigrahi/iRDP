# Run psi blast and parse pssm file and ontain result file
# perl consurf_blast.pl sequenceFile outputprefix logfile
use strict;
use warnings;

sub psiblast
{
print LOG "\n\nRunning PSI-Blast against Uniref90 database\n";
print LOG "\nPlease wait..This may take between 5-15 minutes to complete depending on your sequence length\nYour job is Running...\n\n";

	print `time /programs/ncbi-blast-2.2.29+/bin/psiblast -db $_[0] -query $_[1] -out $_[2] -out_ascii_pssm $_[3] -num_iterations  $_[4] -evalue $_[5] -num_threads $_[6] -html`;
}


sub parsePSSM
{
	my %poshash=qw (A 22  R 23  N  24 D  25 C  26 Q 27  E 28  G 29  H 30  I 31  L 32  K 33  M 34  F  35 P  36 S  37 T 38  W 39  Y 40  V 41);
	my %seqhash=qw (22 A 23  R 24  N  25 D  26 C  27 Q 28  E 29  G 30  H 31  I 32  L 33  K 34  M 35  F  36 P  37 S  38 T 39  W 40  Y 41  V);
	open(PM,$_[0]) or die("cant open $_[0]\n");
	open(OUT,">$_[1]") or die("cant open $_[1]\n");
	print OUT "Col1:Serian number\nCol2:Residue details as [Chain Resno Resid]\nCol3:Residue\n";
  print OUT "Col4:Weighted observed percentages of occurance of the residue (0-1). 1 means occuring all hits, highly conserved\n";
  print OUT "Col5: Value of Col4 scaled between 1-9 where 9 is fully conserved residue\n";
  print OUT "Col6: Other kind of observed residues in the position\n";

	map{
		my $line=$_;
		if($line=~/^\s*\d+ [A-Za-z]/)
		{
			#print "$line";
			my @temp=split / +/,$line;
			shift(@temp);
			#print scalar @temp,"\n";
			my $seq="";
			map{if($temp[$_]!=0){$seq.=$seqhash{$_};}}(22..41);
			if(length($seq)==0){$seq="-";}
			#print OUT "$temp[0]\t$temp[1]\t",($temp[$poshash{$temp[1]}]/100),"\t",sprintf("%.0f",($temp[$poshash{$temp[1]}]/100)*9),"\t$seq\n";
    
    print OUT "$temp[0]\t$temp[1]\t",($temp[$poshash{$temp[1]}]/100),"\t",sprintf("%.0f",(($temp[$poshash{$temp[1]}]*0.08)+1)),"\t$seq\n";
    #The even more general case (when you have a value c between a and b and you want a value x between y and z), x is calculated as follows:
    #x := (c - a) * (z - y) / (b - a) + y

		}
	}<PM>;	
	close(PM);
	close(OUT);
}

# perl psiblast.pl query.txt prefix
my $db="/data/BlastDB/uniref90";
#$db="/data/BlastDB/swissprot";

my $query=$ARGV[0];
my $iter=2;
my $evalue=1;
my $out=$ARGV[1].".html";
my $pssmout=$ARGV[1]."_pssm.txt";
my $pssmfinal=$ARGV[1]."_conserv.txt";
my $nthread=4;
my $logfile="$ARGV[2]";

open(LOG,">>$logfile") or die("cant open logfile");

psiblast($db,$query,$out,$pssmout,$iter,$evalue,$nthread,$logfile);
print LOG "\n\nFinished PSI-BLAST\n   Now parsing PSSM matrix and generating position specific conservation scores\n\n";
parsePSSM($pssmout,$pssmfinal);
print LOG "\n\nFinished generating position specific conservation scores\n\n";
