# perl consurf_blast.pl sequenceFile outputprefix logfile
use strict;
use warnings;

my $qcov_cut=80; #qcov should be atleast 80
my $piden_low=30; #piden should be >30
my $piden_high=100; # piden should be <95
my $tothomo=150; # Total homologous sequence to consider
my $evalcut="0.0001";
my $cdhitcut="0.9";
my $cdhitN="5";
my $dbchoose=30; #cutoff for switching from sp to uniref to nr

sub blastp
{
# -db nr
# query query file
# -out outfile
# evalue 0.0001
#-max_target_seqs 150
	print `time /programs/ncbi-blast-2.2.29+/bin/blastp -db $_[0] -query $_[1] -out $_[2] -evalue $_[3] -num_threads 4 -outfmt "7 qseqid sseqid sgi pident length mismatch gapopen qstart qend sstart send qcovs"`;
}

sub parseBlast
{
open(BL,$_[0]) or die("cant open $_[0]");
open(WL,">$_[1]") or die("cant open $_[1]");


my @blastres=<BL>;
print LOG "*Number of Blast hit:",scalar(@blastres),"\n";
print LOG "*Applying filter: %identity between $piden_low to $piden_high and query coverage > $qcov_cut %\n\n";
my @id;

map
{
	
	my $line=$_;
#	print $line;
	chomp($line);
	if($line!~/^\#/)
	{
	my @temp=split /\t/,"$line";
#	print "$temp[1]\n";
		if($_[2] eq "uniref90")
		{
			#$temp[2]=$temp[1];
			$temp[1]=~/gnl\|unk\|(.*)/;
			$temp[2]=$1;
		}
		if($temp[11]>$qcov_cut && $temp[3]>$piden_low && $temp[3] < $piden_high)
		{
		print WL "$line\n";
		push @id,$temp[2];	
		}
	
	}else{
    if($line!~/\# Database/){print WL "$line\n";}
    if($line!~/hits found/){print WL "$line\n";}

	}
}@blastres;

#if(@id>$tothomo){@id=@id[1..$tothomo]}
#print @id;
return([@id]);
}


sub fetchFasta
{
`/programs/ncbi-blast-2.2.29+/bin/blastdbcmd -db $_[2] -dbtype prot -entry $_[0] > $_[1]`;
}

sub cdhit
{
	`/programs/cd-hit-v4.6.1-2012-08-27/cd-hit -i $_[0] -o $_[1] -c $_[2] -n $_[3] -T 4`;
}

sub mafft
{
	`mafft  --quiet --thread 4 $_[0] > $_[1]`;
}

my $prefix=$ARGV[1];
my $blastout="$prefix".".txt";
my $blastfiltfile="$prefix"."_filt.txt";
my $logfile=$ARGV[2];

=head
First we search against swissprot database and obtain the blast tabular result
We filter based on qcovs and iden cutoff
[Not yet implemented]: if less than 50 query, we can then retry using uniref90 database
we use all ids to fetch their fasta sequence
We then cluster all hits at 90% identity
We then do MSA
=cut

=head
blastp($ARGV[2],$ARGV[0],$blastout,"0.0001");
my @id=@{parseBlast($blastout,$blastfiltfile)};
my $seqid=join ",",@id;
print $seqid;
fetchFasta($seqid,$blastseq);
cdhit($blastseq,$blastseq90);
mafft($blastseq90,$msa);
=cut

open(LOG,">>$logfile") or die("cant open logfile");
# First Blast against Swissprot database
my $db="/data/BlastDB/swissprot";
print LOG "\n\nRunning Blast against SwissProt database\n";
print LOG "\nPlease wait..This may take appx 2 minutes to complete depending on sequence length\nYour job is Running...\n\n";

blastp($db,$ARGV[0],$blastout,$evalcut); #perform blast
my @id=@{parseBlast($blastout,$blastfiltfile,"swissprot")}; #parse blast output and apply filter to get @id and filter file;
print LOG "*Number of filtered hits found:",scalar(@id),"\n";

if(@id<$dbchoose)
{
	print LOG "*Number of filtered hits less than cutoff:",$dbchoose,"\n        Therefore\n";
  @id=();
	#`rm $blastout`;
	#`rm $blastfiltfile`;
	#$prefix=$prefix."_ur";
  $blastout="$prefix".".txt";$blastfiltfile="$prefix"."_filt.txt";

	print("Removed old blast output files\n");
	print LOG "Running Blast against Uniref90 database\n";
  print LOG "\nPlease wait..This may take appx 2-4 minutes to complete depending on sequence length.\n Your job is running....\n\n";

	$db="/data/BlastDB/uniref90";
	blastp($db,$ARGV[0],$blastout,$evalcut); #perform blast
	@id=@{parseBlast($blastout,$blastfiltfile,"uniref90")}; #parse blast output and apply filter to get @id and filter file;
	print LOG "*Number of filtered hits found:",scalar(@id),"\n";
	
	if(@id<$dbchoose)
	{
	@id=();
  print LOG "*Number of filtered hits less than cutoff:",$dbchoose,"\n        Therefore\n";
	#`rm $blastout`;
	#`rm $blastfiltfile`;
	#$prefix=$prefix."_nr";
  $blastout="$prefix".".txt";$blastfiltfile="$prefix"."_filt.txt";
	print("Removed old blast output files\n");
	print LOG "Running Blast against nr database\n";
  print LOG "\nPlease wait..This may take appx 5 minutes to complete depending on sequence length. \nYour job is running.....\n\n";

	$db="/data/BlastDB/nr";
	blastp($db,$ARGV[0],$blastout,$evalcut); #perform blast
	@id=@{parseBlast($blastout,$blastfiltfile,"nr")}; #parse blast output and apply filter to get @id and filter file;
	print LOG "*Number of filtered hits found:",scalar(@id),"\n";
	}
}

#print "\n\n\n@id\n\n";
my $seqid=join ",",@id;
my $blastseq="$prefix"."_seq.fasta";
my $blastseq90="$prefix"."_seq90.fasta";

# At least 2 sequences should be there so that we can fetch fasta sequences
# Else no cdhit and no msa

if(@id>1)
{
          fetchFasta($seqid,$blastseq,$db);
          
          print LOG "Fetched sequences of all hits\n";
          print LOG "Running CD-HIT at cutoff:",$cdhitcut,"\n";
          cdhit($blastseq,$blastseq90,$cdhitcut,$cdhitN);
          my $cdhitcount=`grep \">\" $blastseq90 \| wc -l`;
          print LOG "Number of hits after CD-HIT:",$cdhitcount,"\n";
          
          while($cdhitcount>$tothomo && $cdhitcut > 0.7)
          {
          	
          	$cdhitcut=$cdhitcut-0.1;
          	if($cdhitcut == 0.6){print("hiii");last;}
          	
          	#if($cdhitcut==0.6){$cdhitN=4;}elsif($cdhitcut==0.5){$cdhitN=3;}elsif($cdhitcut==0.4){$cdhitN=2;}
          	print LOG "Running CD-HIT at cutoff:",$cdhitcut,"\n";
          	if($cdhitcut>=0.7)
          	{
          	cdhit($blastseq,$blastseq90,$cdhitcut,$cdhitN);
          	$cdhitcount=`grep \">\" $blastseq90 \| wc -l`;
          	print LOG "Number of hits after CD-HIT:",$cdhitcount,"\n";
          	}
          }
          
          print LOG "\n\nFinal Number of hits after CD-HIT:",$cdhitcount,"\n";
          
          my $blastfinal="$prefix"."_final.fasta";
          
          if($cdhitcount>$tothomo)
          {
          	print LOG "Taking only $tothomo sequences only\n";
          	open(FA,$blastseq90)or die("cant open $blastseq90");
          	open(FB,">$blastfinal")or die("cant open $blastfinal");
          	
          	my $gicount=0;
          
          	foreach my $line (<FA>)
          	{
          		if($line=~/^>/){$gicount++;}
          		if($gicount==$tothomo+1)
          		{
          		last;
          		}
          		print FB "$line";
          	}
          	close(FA);
          }else{
            `cp $blastseq90 $blastfinal`;
          }
          
          
          `cat $ARGV[0] >> $blastfinal`;

        # If after cdhit we get 1 sequence we dont do msa
        if($cdhitcount>1)
        {
          my $msa="$prefix".".fastaaln";
          print LOG "Running MSA\n";
          mafft($blastfinal,$msa);
        }
}

close(LOG);












