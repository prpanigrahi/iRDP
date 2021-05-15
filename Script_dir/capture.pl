BEGIN 
{
push @INC, '/opt/ActivePerl-5.16/lib'
}
print "@INC";
use WWW::Mechanize;

my $url = "http://capture.caltech.edu/";
my $mech = WWW::Mechanize->new();
$mech->proxy(['http','ftp'],'');

my @id=("1GK9");

open(F,">captureout.txt");

foreach my $id(@id)
{
print $id;

### First form #####
$mech->get($url);
$mech->form_number(1); #select 1st form
$mech->set_fields("data",$id);
$mech->submit();
print F $mech->content();
}


### second form #####
=head
$mech->form_number(2); #select 2nd form
$mech->set_fields("upfile","/home/priyabrata/Work/ProtTS/Test30Jan2013/PdPGA.pdb");
$mech->submit();
print $mech->content();
=cut

