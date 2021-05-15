BEGIN 
{
push @INC, '/home/priyabrata/Software/ActivePerl-5.14/lib'
}

#print "@INC\n";
use WWW::Mechanize;

#=head
my $url = "http://mordred.bioc.cam.ac.uk/sdm/sdm.php";
my $mech = WWW::Mechanize->new();
$mech->proxy(['http','ftp'],'');

$mech->get($url);

# '/home/priyabrata/Work/ProtTS/SASMAT_ason_29sept2013_before4th5thPres/1gk9.pdb'
# 'A'
# 'G'
# 184

$mech->submit_form(
               form_number => 3,
               fields      => {
                   wt_pdbfile    => $ARGV[0],
                   #pdb_code => '1pnk',
                   pdb_chain    => $ARGV[1],
                   mt_res  => $ARGV[2],
                   mt_pos => $ARGV[3],
               }
           );
my $res=$mech->text();        
$res=~/Pseudo DELTA DELTA G =\s+([\+\-0-9\.]+)/;
print "$1";

     
=head
#$mech->form_number(2); #select 1st form
$mech->form_with_fields("mt_res");
$mech->set_fields("wt_pdbfile","/home/priyabrata/Work/ProtTS/SASMAT_ason_29sept2013_before4th5thPres/1gk9.pdb");
$mech->set_fields("pdb_chain","A");
$mech->set_fields("mt_res","G");
$mech->set_fields("mt_pos","184");

#print($mech->field("pdb_code"));

$mech->submit_form(form_number=> 2);
print $mech->content();

#print $mech->forms;




