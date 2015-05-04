#!/usr/bin/perl

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Program HTS_filtered_codon_mutations.pl
#
# Reads, filters and analyzes high-throughput sequencing (HTS) data (aligned by http://www.novocraft.com/products/novoalign/)
#
# Input: alignment of the Illumina sequencing reads, WT sequence (fasta), gene_position_start, gene_position_end
#
# Output: *.allmut (all mutations per codons per position), *.freq (frequencies of mutations (number of occurences/number of reads))
# *.ns1 (single numcleotide changes), *.ns2 (double nucleotide changes), *.ns3 (triple nucleotide changes),
# *.syn (synonymous changes), *.indels (insertions and deletions)
# *.codoposinreads (position within the reads (3-37 out of 40 bp long reads) where the mutation occurs)
#
# Created by AGNES Toth-Petroczy, 31.05.2011
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Example:
# ./HTS_filtered_codon_mutations.pl novo.align G0.txt 31 1151

use strict;
use warnings;
my $scriptDirectory="./";
require "$scriptDirectory/DNACode.pm";

die "usage: perl filename WT_seq genepos1 genepos2\n" unless (scalar (@ARGV) ==4);


my $timings=time();
print "HTS_filtered_codon_mutations.pl program started :)\n";

my $file=$ARGV[0];
my $wtseq=$ARGV[1];
my $gp1=$ARGV[2];
my $gp2=$ARGV[3];
my $line;
my $AllSequences=[];
my $nseq=0;

#Read WT seq
my @WTseq;
open (WT, $wtseq) or die "Could not open file $wtseq\n";
while (my $line=<WT>)
{
 if ($line !~ m/>/)
 {
  chomp($line);
  my @junk=split(//,$line);
  push (@WTseq,@junk);
 }
}

# Assign WT codons
my @WTcodon;
my $counter=0;
for (my $i=$gp1; $i <=$gp2; $i=$i+3)
{
 $WTcodon[$counter]=$WTseq[$i-1].$WTseq[$i].$WTseq[$i+1];
 $counter++;
}
close WT;

#Read and analyze sequencing data (novo.align)#############################################################################################################

my %SubstitutionMatrix;
my %CodonPerPos;
my @coverage;
my %PositionInReads;
open (INPUT_FILE_HANDLE, $file) or die "Could not open file $file\n";
open (INDELS, '>', $file.'.indels');

print "Analyzing $file alignment\n";

for(my $i=1; $line = <INPUT_FILE_HANDLE>; $i++)
{
 if ($line=~ m/^@/)
 { 
  my $Sequence;
  my @info;
  $nseq++;
  chomp($line);
  @info = split(/\t/,$line);
   my $columns = @info;
   die if ($columns > 14);
   $Sequence->{"ID"} = $info[0];
   $Sequence->{"SEQUENCE"} = $info[2];
   $Sequence->{"QUALITY"} = $info[3];
   $Sequence->{"INFO"} = $info[4];
   next if ($Sequence->{"INFO"} =~ m/NM|QC/);
   $Sequence->{"REFSEQ_ID"} = $info[7];
   $Sequence->{"OFFSET"} = $info[8];
   $Sequence->{"DIRECTION"} = $info[9];
   $Sequence->{"MUTATIONS"}=$info[13];

   #Remove reads where quiality < 20 in any of the bases
   my @quality=parseQualityLine($Sequence->{"QUALITY"});
   my @qchar= split(/\s+/,$Sequence->{"QUALITY"});

  if (scalar(@info) > 6)
  {

   #Checking indels - skipping reads with indels for SNP anal
   if (($Sequence->{"MUTATIONS"}) and ($Sequence->{"MUTATIONS"} =~ m/\-|\+/))
    {
     print INDELS "Indel $Sequence->{OFFSET} $Sequence->{MUTATIONS}\n";
     next;
    }
   #Ignore read if it contains N (non-identified base)
   next if ($Sequence->{"SEQUENCE"} =~ m/N/);
   #Chech where does the first full codon start
   # ignore incomplete codons
   my ($firstbasepos, $firstreadpos, $firstcodon);
   if (($Sequence->{"OFFSET"} >= $gp1) and ($Sequence->{"OFFSET"} <=$gp2))    
   {
    if ((int($Sequence->{"OFFSET"} - $gp1)/3)!~ m/\D/)
    {
     $firstbasepos=$Sequence->{"OFFSET"};
     $firstreadpos=0;
    }
    elsif ((int($Sequence->{"OFFSET"} - $gp1 -1 )/3 )!~ m/\D/)
    {
     $firstbasepos=$Sequence->{"OFFSET"} + 2;
     $firstreadpos=2;
    }
    elsif ((int($Sequence->{"OFFSET"} - $gp1 - 2)/3 )!~ m/\D/)
    {
     $firstbasepos=$Sequence->{"OFFSET"} + 1;
     $firstreadpos=1;
    }
    $firstcodon=($firstbasepos - $gp1)/3;
   }
   else { next;}
   #print "$firstcodon\n";

   #Identify sequenced codons
    my $cc=$firstcodon;
    my @mutcodon;
    my @readseq= split(//, $Sequence->{"SEQUENCE"});
    #if forward strand
    if ($Sequence->{"DIRECTION"} eq 'F')
    {
     for (my $f=$firstreadpos; $f < 38; $f=$f+3)
     {
      if ((($quality[$f] < 20) or ($quality[$f+1] < 20)) or ($quality[$f+2] < 20)) {$cc++; next;} #skip codon

      $mutcodon[$cc]=$readseq[$f].$readseq[$f+1].$readseq[$f+2];
      #if it contains not identified base (N), drop it
      if ($mutcodon[$cc] =~ m/N/) {next;}
      #print "$cc, $firstreadpos, $mutcodon[$cc], $WTcodon[$cc]\n"; 
      $coverage[$cc]++;
      $PositionInReads{$cc}{$mutcodon[$cc]}{$f}++;      

      if ($WTcodon[$cc])
      {
       if ($mutcodon[$cc] ne $WTcodon[$cc])
       {
        $CodonPerPos{"$cc"}{$mutcodon[$cc]}++;
        $SubstitutionMatrix{$WTcodon[$cc].'->'.$mutcodon[$cc]}++;
       }
       $cc++;
      }
      else {last;}
     }
    }
    else
    #if reverse strand
    {
     my $firstreadreverse;
     if ($firstreadpos == 0) {$firstreadreverse=39;} #1
     elsif ($firstreadpos== 1) {$firstreadreverse=38;} #0
     elsif ($firstreadpos== 2) {$firstreadreverse=37;} #2

     for (my $f=$firstreadreverse; $f > 1; $f=$f-3)
     {
      if ((($quality[$f] < 20) or ($quality[$f-1] < 20)) or ($quality[$f-2] < 20)) {$cc++; next;}

      $mutcodon[$cc]=$readseq[$f].$readseq[$f-1].$readseq[$f-2];
      #if it contains not identified base (N), drop it
      if ($mutcodon[$cc] =~ m/N/) {next;}
      $mutcodon[$cc]=ReverseStrand($mutcodon[$cc]);
      $coverage[$cc]++;

      $PositionInReads{$cc}{$mutcodon[$cc]}{$f}++;
 
     #print "$cc, $WTcodon[$cc], $mutcodon[$cc]\n";
      if ($WTcodon[$cc])
      {
       if ($mutcodon[$cc] ne $WTcodon[$cc])
       {
        $CodonPerPos{"$cc"}{$mutcodon[$cc]}++;
        $SubstitutionMatrix{$WTcodon[$cc].'->'.$mutcodon[$cc]}++;
       }
       $cc++;
      }
      else {last;}
     }

    }


  } #long enough line
 }
}
close(INPUT_FILE_HANDLE);


#Print output ##########################################################################################################################################################
my @codons = ("TAA", "TAG","TGA","GCA","GCC","GCG","GCT","TGC","TGT","GAC","GAT","GAA","GAG","TTC","TTT","GGA","GGC","GGG","GGT","CAC","CAT","ATA","ATC","ATT","AAA","AAG","CTA","CTC","CTG","CTT","TTA","TTG","ATG","AAC","AAT","CCA","CCC","CCG","CCT","CAA","CAG","AGA","AGG","CGA","CGC","CGG","CGT","AGC","AGT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","GTA","GTC","GTG","GTT","TGG","TAC","TAT");

open (OUTPOS, '>', $wtseq.'.allmut');
open (OUT1, '>', $wtseq.'.ns1');
open (OUT2, '>', $wtseq.'.ns2');
open (OUT3, '>', $wtseq.'.ns3');
open (FREQ, '>', $wtseq.'.freq');
open (SYN, '>', $wtseq.'.syn');
open (POSINREADS, '>', $wtseq.'.codonposinreads');

#for all codon positions
for (my $i=0; $i < ($gp2-$gp1+1)/3; $i++)
{
 my $totalfreq=0;
 my $synfreq=0;
 my $wtaa=DNACode::Translation($WTcodon[$i]);
 printf OUTPOS ("%3d %3s %1s %d", $i-20, $WTcodon[$i], $wtaa, $coverage[$i]);
 printf OUT1 ("%3d %3s %1s", $i-20, $WTcodon[$i], $wtaa);
 printf OUT2 ("%3d %3s %1s", $i-20, $WTcodon[$i], $wtaa);
 printf OUT3 ("%3d %3s %1s", $i-20, $WTcodon[$i], $wtaa);

 foreach my $key (sort {$CodonPerPos{"$i"}{$b} <=> $CodonPerPos{"$i"}{$a}} keys %{$CodonPerPos{"$i"}})
 {
  my $aa=DNACode::Translation($key);
  my $freq=($CodonPerPos{"$i"}{$key}/$coverage[$i])*100;
  #is it single, double or triple nucleotide mutation?
  my $nsub=0;
  my @wtjunk=split(//, $WTcodon[$i]);
  my @mutjunk=split(//, $key);
  for (my $j=0; $j < 3; $j++)
  {
   if ($wtjunk[$j] ne $mutjunk[$j])
   {
    $nsub++;
   }
  }
  #Single nucleotide mutation
  if ($nsub == 1)
  {
   printf OUT1 (" %6.4f %3s %1s",$freq, $key,$aa);
  }
  #Double nucleotide mutation
  if ($nsub == 2)
  {
   printf OUT2 (" %6.4f %3s %1s",$freq, $key,$aa);
  }
  #Triple nucelotide mutation
  if ($nsub == 3)
  {
   printf OUT3 (" %6.4f %3s %1s",$freq, $key,$aa);
  }
  #synonymous changes
  if ($aa eq $wtaa)
  {
   printf SYN (" %6.4f %3s %1s",$freq, $key,$aa);
   $synfreq=$synfreq+$freq;
  }
  else
  #nonsynomynous changes
  {
   if ($freq > 0.005)
   {
#    printf OUTPOS (" %6.4f %3s %1s",$freq, $key,$aa);
    $totalfreq=$totalfreq+$freq;
   }
  }
 printf OUTPOS (" %6.4f %3s %1s",$freq, $key,$aa);
#printf OUTPOS (" %6.4f %3s %1s",$CodonPerPos{"$i"}{$key}, $key,$aa);
 }
 print OUTPOS "\n";
 print OUT1 "\n";
 print OUT2 "\n";
 print OUT3 "\n";
 print SYN "\n";
 printf FREQ ("%3d %3s %1s %6.4f %6.4f\n", $i-20, $WTcodon[$i], $wtaa, $totalfreq, $synfreq); 

 #Print positions in reads reverse (from 39 to 2) and forward (from 0 to 37)
  for (my $p=0; $p <= 39; $p++)
  {
    printf POSINREADS ("%3d %3s %d %d",$i-20, $WTcodon[$i], $coverage[$i], $p);
    foreach my $c (@codons)
    {
        if ($PositionInReads{$i}{$c}{$p})
        {
        print POSINREADS " $PositionInReads{$i}{$c}{$p}";
        }
        else
        {
          print POSINREADS " 0";
        }
    }
 
  print POSINREADS "\n";
  }
} #i codons
close OUTPOS;
close FREQ;
close OUT1;
close OUT2;
close OUT3;
close POSINREADS;

#Print subsitution matrix
open (OUT, ">", $file."_SubMat.dat");
for my $key (keys %SubstitutionMatrix)
{
 my $wt=substr($key,0,3);
 my $mut=substr($key,5,3);
 my $wtaa=DNACode::Translation($wt);
 my $mutaa=DNACode::Translation($mut);
 printf OUT ("%s %s%s%s %10d\n", $key, $wtaa, "->", $mutaa, $SubstitutionMatrix{$key});
}
close OUT;

#Subroutines ###################################################################################################################

#Illumina sequencing quality information
sub parseQualityLine
{
        my ($line) = @_;
        my @qchar = split("", $line);
        my @quals      = (0) x ( scalar @qchar);
        for (my $i=0; $i <= $#qchar; ++$i)
        {
                $quals[$i] = unpack("C",$qchar[$i]) - 33;
        }
        return @quals;
}


sub ReverseStrand
{
 my ($string)=@_;
 my @junk=split(//,$string);
 for (my $i=0; $i < scalar(@junk); $i++)
 {
  if ($junk[$i] eq 'A') {$junk[$i] = 'T';}
  elsif ($junk[$i] eq 'T') {$junk[$i] = 'A';}
  elsif ($junk[$i] eq 'G') {$junk[$i] = 'C';}
  elsif ($junk[$i] eq 'C') {$junk[$i] = 'G';}
 }
 my $reverse=join('',@junk);
 return ($reverse);
}
