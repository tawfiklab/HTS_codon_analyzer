package DNACode;

use strict;
use warnings;

use base 'Exporter';
our @EXPORT_OK = qw(CodonNumber Translation getCodon getRate);


sub CodonNumber
{
 my $ID= $_[0];

 my $numCodons = {
                  "A" => 4,
                  "R" => 6,
                  "N" => 2,
                  "D" => 2,
                  "C" => 2,
                  "E" => 2,
                  "Q" => 2,
                  "G" => 4,
                  "H" => 2,
                  "I" => 3,
                  "L" => 6,
                  "K" => 2,
                  "M" => 1,
                  "F" => 2,
                  "P" => 4,
                  "S" => 6,
                  "T" => 4,
                  "W" => 1,
                  "Y" => 2,
                  "V" => 4,
                  "Stop" => 3 };

 return($numCodons->{$ID});
}

sub Translation
{
 my $Codon=$_[0];
 my $CodonTable =
 {
  'GCG' => 'A',
  'GCC' => 'A',
  'GCA' => 'A',
  'GCT' => 'A',
  'CGC' => 'R',
  'CGT' => 'R',
  'CGG' => 'R',
  'CGA' => 'R',
  'AGA' => 'R',
  'AGG' => 'R',
  'AAC' => 'N',
  'AAT' => 'N',
  'GAT' => 'D',
  'GAC' => 'D',
  'TGC' => 'C',
  'TGT' => 'C',
  'GAG' => 'E',
  'GAA' => 'E',
  'CAG' => 'Q',
  'CAA' => 'Q',
  'GGC' => 'G',
  'GGT' => 'G',
  'GGG' => 'G',
  'GGA' => 'G',
  'CAT' => 'H',
  'CAC' => 'H',
  'ATT' => 'I',
  'ATC' => 'I',
  'ATA' => 'I',
  'CTG' => 'L',
  'TTG' => 'L',
  'TTA' => 'L',
  'CTT' => 'L',
  'CTC' => 'L',
  'CTA' => 'L',
  'AAA' => 'K',
  'AAG' => 'K',
  'ATG' => 'M',
  'TTT' => 'F',
  'TTC' => 'F',
  'CCG' => 'P',
  'CCA' => 'P',
  'CCT' => 'P',
  'CCC' => 'P',
  'AGC' => 'S',
  'AGT' => 'S',
  'TCG' => 'S',
  'TCT' => 'S',
  'TCC' => 'S',
  'TCA' => 'S',
  'TAA' => 'Stop',
  'TGA' => 'Stop',
  'TAG' => 'Stop',
  'ACC' => 'T',
  'ACG' => 'T',
  'ACT' => 'T',
  'ACA' => 'T',
  'TGG' => 'W',
  'TAT' => 'Y',
  'TAC' => 'Y',
  'GTG' => 'V',
  'GTC' => 'V',
  'GTA' => 'V',
  'GTT' => 'V'
 };
 if ($CodonTable->{"$Codon"})
 {
  return ($CodonTable->{"$Codon"});
 }
 else {return 1000;}
}

sub getCodon
{
my $DNACodons;
my $ID= $_[0];
my $c= $_[1];

$DNACodons->{"A"}->{"CODON"}->[0]='GCG';
$DNACodons->{"A"}->{"CODON"}->[1]='GCC';
$DNACodons->{"A"}->{"CODON"}->[2]='GCA';
$DNACodons->{"A"}->{"CODON"}->[3]='GCT';
$DNACodons->{"R"}->{"CODON"}->[0]='CGC';
$DNACodons->{"R"}->{"CODON"}->[1]='CGT';
$DNACodons->{"R"}->{"CODON"}->[2]='CGG';
$DNACodons->{"R"}->{"CODON"}->[3]='CGA';
$DNACodons->{"R"}->{"CODON"}->[4]='AGA';
$DNACodons->{"R"}->{"CODON"}->[5]='AGG';
$DNACodons->{"N"}->{"CODON"}->[0]='AAC';
$DNACodons->{"N"}->{"CODON"}->[1]='AAT';
$DNACodons->{"D"}->{"CODON"}->[0]='GAT';
$DNACodons->{"D"}->{"CODON"}->[1]='GAC';
$DNACodons->{"C"}->{"CODON"}->[0]='TGC';
$DNACodons->{"C"}->{"CODON"}->[1]='TGT';
$DNACodons->{"E"}->{"CODON"}->[0]='GAA';
$DNACodons->{"E"}->{"CODON"}->[1]='GAG';
$DNACodons->{"Q"}->{"CODON"}->[0]='CAG';
$DNACodons->{"Q"}->{"CODON"}->[1]='CAA';
$DNACodons->{"G"}->{"CODON"}->[0]='GGC';
$DNACodons->{"G"}->{"CODON"}->[1]='GGT';
$DNACodons->{"G"}->{"CODON"}->[2]='GGG';
$DNACodons->{"G"}->{"CODON"}->[3]='GGA';
$DNACodons->{"H"}->{"CODON"}->[0]='CAT';
$DNACodons->{"H"}->{"CODON"}->[1]='CAC';
$DNACodons->{"I"}->{"CODON"}->[0]='ATT';
$DNACodons->{"I"}->{"CODON"}->[1]='ATC';
$DNACodons->{"I"}->{"CODON"}->[2]='ATA';
$DNACodons->{"L"}->{"CODON"}->[0]='CTG';
$DNACodons->{"L"}->{"CODON"}->[1]='TTG';
$DNACodons->{"L"}->{"CODON"}->[2]='TTA';
$DNACodons->{"L"}->{"CODON"}->[3]='CTT';
$DNACodons->{"L"}->{"CODON"}->[4]='CTC';
$DNACodons->{"L"}->{"CODON"}->[5]='CTA';
$DNACodons->{"K"}->{"CODON"}->[0]='AAA';
$DNACodons->{"K"}->{"CODON"}->[1]='AAG';
$DNACodons->{"M"}->{"CODON"}->[0]='ATG';
$DNACodons->{"F"}->{"CODON"}->[0]='TTT';
$DNACodons->{"F"}->{"CODON"}->[1]='TTC';
$DNACodons->{"P"}->{"CODON"}->[0]='CCG';
$DNACodons->{"P"}->{"CODON"}->[1]='CCA';
$DNACodons->{"P"}->{"CODON"}->[2]='CCT';
$DNACodons->{"P"}->{"CODON"}->[3]='CCC';
$DNACodons->{"S"}->{"CODON"}->[0]='AGC';
$DNACodons->{"S"}->{"CODON"}->[1]='AGT';
$DNACodons->{"S"}->{"CODON"}->[2]='TCG';
$DNACodons->{"S"}->{"CODON"}->[3]='TCT';
$DNACodons->{"S"}->{"CODON"}->[4]='TCC';
$DNACodons->{"S"}->{"CODON"}->[5]='TCA';
$DNACodons->{"Stop"}->{"CODON"}->[0]='TAA';
$DNACodons->{"Stop"}->{"CODON"}->[1]='TGA';
$DNACodons->{"Stop"}->{"CODON"}->[2]='TAG';
$DNACodons->{"T"}->{"CODON"}->[0]='ACC';
$DNACodons->{"T"}->{"CODON"}->[1]='ACG';
$DNACodons->{"T"}->{"CODON"}->[2]='ACT';
$DNACodons->{"T"}->{"CODON"}->[3]='ACA';
$DNACodons->{"W"}->{"CODON"}->[0]='TGG';
$DNACodons->{"Y"}->{"CODON"}->[0]='TAT';
$DNACodons->{"Y"}->{"CODON"}->[1]='TAC';
$DNACodons->{"V"}->{"CODON"}->[0]='GTG';
$DNACodons->{"V"}->{"CODON"}->[1]='GTT';
$DNACodons->{"V"}->{"CODON"}->[2]='GTC';
$DNACodons->{"V"}->{"CODON"}->[3]='GTA';

return ($DNACodons->{"$ID"}->{"CODON"}->[$c]);
}

sub getRate
{
my $ID= $_[0];
my $c= $_[1];
my $DNACodons;

$DNACodons->{"A"}->{"RATE"}->[0]=0.36;
$DNACodons->{"A"}->{"RATE"}->[1]=0.27;
$DNACodons->{"A"}->{"RATE"}->[2]=0.21;
$DNACodons->{"A"}->{"RATE"}->[3]=0.16;
$DNACodons->{"R"}->{"RATE"}->[0]=0.40;
$DNACodons->{"R"}->{"RATE"}->[1]=0.38;
$DNACodons->{"R"}->{"RATE"}->[2]=0.10;
$DNACodons->{"R"}->{"RATE"}->[3]=0.06;
$DNACodons->{"R"}->{"RATE"}->[4]=0.04;
$DNACodons->{"R"}->{"RATE"}->[5]=0.02;
$DNACodons->{"N"}->{"RATE"}->[0]=0.55;
$DNACodons->{"N"}->{"RATE"}->[1]=0.45;
$DNACodons->{"D"}->{"RATE"}->[0]=0.63;
$DNACodons->{"D"}->{"RATE"}->[1]=0.37;
$DNACodons->{"C"}->{"RATE"}->[0]=0.55;
$DNACodons->{"C"}->{"RATE"}->[1]=0.45;
$DNACodons->{"E"}->{"RATE"}->[0]=0.69;
$DNACodons->{"E"}->{"RATE"}->[1]=0.31;
$DNACodons->{"Q"}->{"RATE"}->[0]=0.65;
$DNACodons->{"Q"}->{"RATE"}->[1]=0.35;
$DNACodons->{"G"}->{"RATE"}->[0]=0.40;
$DNACodons->{"G"}->{"RATE"}->[1]=0.34;
$DNACodons->{"G"}->{"RATE"}->[2]=0.15;
$DNACodons->{"G"}->{"RATE"}->[3]=0.11;
$DNACodons->{"H"}->{"RATE"}->[0]=0.57;
$DNACodons->{"H"}->{"RATE"}->[1]=0.43;
$DNACodons->{"I"}->{"RATE"}->[0]=0.51;
$DNACodons->{"I"}->{"RATE"}->[1]=0.42;
$DNACodons->{"I"}->{"RATE"}->[2]=0.07;
$DNACodons->{"L"}->{"RATE"}->[0]=0.50;
$DNACodons->{"L"}->{"RATE"}->[1]=0.13;
$DNACodons->{"L"}->{"RATE"}->[2]=0.13;
$DNACodons->{"L"}->{"RATE"}->[3]=0.10;
$DNACodons->{"L"}->{"RATE"}->[4]=0.10;
$DNACodons->{"L"}->{"RATE"}->[5]=0.04;
$DNACodons->{"K"}->{"RATE"}->[0]=0.77;
$DNACodons->{"K"}->{"RATE"}->[1]=0.23;
$DNACodons->{"M"}->{"RATE"}->[0]=1.00;
$DNACodons->{"F"}->{"RATE"}->[0]=0.57;
$DNACodons->{"F"}->{"RATE"}->[1]=0.43;
$DNACodons->{"P"}->{"RATE"}->[0]=0.52;
$DNACodons->{"P"}->{"RATE"}->[1]=0.19;
$DNACodons->{"P"}->{"RATE"}->[2]=0.16;
$DNACodons->{"P"}->{"RATE"}->[3]=0.12;
$DNACodons->{"S"}->{"RATE"}->[0]=0.28;
$DNACodons->{"S"}->{"RATE"}->[1]=0.15;
$DNACodons->{"S"}->{"RATE"}->[2]=0.15;
$DNACodons->{"S"}->{"RATE"}->[3]=0.15;
$DNACodons->{"S"}->{"RATE"}->[4]=0.15;
$DNACodons->{"S"}->{"RATE"}->[5]=0.12;
$DNACodons->{"Stop"}->{"RATE"}->[0]=0.64;
$DNACodons->{"Stop"}->{"RATE"}->[1]=0.29;
$DNACodons->{"Stop"}->{"RATE"}->[2]=0.07;
$DNACodons->{"T"}->{"RATE"}->[0]=0.44;
$DNACodons->{"T"}->{"RATE"}->[1]=0.27;
$DNACodons->{"T"}->{"RATE"}->[2]=0.17;
$DNACodons->{"T"}->{"RATE"}->[3]=0.13;
$DNACodons->{"W"}->{"RATE"}->[0]=1.00;
$DNACodons->{"Y"}->{"RATE"}->[0]=0.57;
$DNACodons->{"Y"}->{"RATE"}->[1]=0.43;
$DNACodons->{"V"}->{"RATE"}->[0]=0.37;
$DNACodons->{"V"}->{"RATE"}->[1]=0.26;
$DNACodons->{"V"}->{"RATE"}->[2]=0.22;
$DNACodons->{"V"}->{"RATE"}->[3]=0.15;

return($DNACodons->{"$ID"}->{"RATE"}->[$c]);

}
       
