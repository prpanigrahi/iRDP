open(F,"$ARGV[0]");
open(F1,">$ARGV[1]");

my $file=join "",<F>;


############## RF ################
$file=~/Number of ARG\/PHE interacting pairs:\s+(\d+)/;
my $RF=$1;
$file=~/Number of energetically significant ARG\/PHE cation-pi interactions:\s+(\d+)/;
my $RFsig=$1;
$file=~/The following are the.*ARG\/PHE.*\n.*\nCation.*\n.*\n.*\n((.*\n){$RFsig})/;
my $RFdet=$1;


############## RY ################
$file=~/Number of ARG\/TYR interacting pairs:\s+(\d+)/;
my $RY=$1;
$file=~/Number of energetically significant ARG\/TYR cation-pi interactions:\s+(\d+)/;
my $RYsig=$1;
$file=~/The following are the.*ARG\/TYR.*\n.*\nCation.*\n.*\n.*\n((.*\n){$RYsig})/;
my $RYdet=$1;


############## RW ################
$file=~/Number of ARG\/TRP interacting pairs:\s+(\d+)/;
my $RW=$1;
$file=~/Number of energetically significant ARG\/TRP cation-pi interactions:\s+(\d+)/;
my $RWsig=$1;
$file=~/The following are the.*ARG\/TRP.*\n.*\nCation.*\n.*\n.*\n((.*\n){$RWsig})/;
my $RWdet=$1;

############## KF ################
$file=~/Number of LYS\/PHE interacting pairs:\s+(\d+)/;
my $KF=$1;
$file=~/Number of energetically significant LYS\/PHE cation-pi interactions:\s+(\d+)/;
my $KFsig=$1;
$file=~/The following are the.*LYS\/PHE.*\n.*\nCation.*\n.*\n.*\n((.*\n){$KFsig})/;
my $KFdet=$1;


############## KY ################
$file=~/Number of LYS\/TYR interacting pairs:\s+(\d+)/;
my $KY=$1;
$file=~/Number of energetically significant LYS\/TYR cation-pi interactions:\s+(\d+)/;
my $KYsig=$1;
$file=~/The following are the.*LYS\/TYR.*\n.*\nCation.*\n.*\n.*\n((.*\n){$KYsig})/;
my $KYdet=$1;

############## KW ################
$file=~/Number of LYS\/TRP interacting pairs:\s+(\d+)/;
my $KW=$1;
$file=~/Number of energetically significant LYS\/TRP cation-pi interactions:\s+(\d+)/;
my $KWsig=$1;
$file=~/The following are the.*LYS\/TRP.*\n.*\nCation.*\n.*\n.*\n((.*\n){$KWsig})/;
my $KWdet=$1;

print F1 "Cation-pi interaction details\n---------------------------------------------------\n";
print F1 "Resid\tResno\tChain\tResid\tResno\tChain\tEes\tEvdw\n";
print F1 "$RFdet";
print F1 "$RYdet";
print F1 "$RWdet";
print F1 "$KFdet";
print F1 "$KYdet";
print F1 "$KWdet";

print F1 "\n\n\nType\tTotal\tSignificant\n";
print F1  "RF\t$RF\t$RFsig\n";
print F1 "RY\t$RY\t$RYsig\n";
print F1  "RW\t$RW\t$RWsig\n";
print F1  "KF\t$KF\t$KFsig\n";
print F1  "KY\t$KY\t$KYsig\n";
print F1  "KW\t$KW\t$KWsig\n\n\n";


$file=~/There is a total of\s+(\d+).*\n/;
print F1 "$&";
$file=~/\* There are\s+(\d+)\s+.*<=.*\n/;
print F1 "$&";
$file=~/\* There are\s+(\d+)\s+.*between.*\n.*\n/;
print F1 "$&";
$file=~/There is an average of (\d+).*\n/;
print F1 "$&";

close(F);
