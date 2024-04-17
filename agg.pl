#! /usr/bin/perl

use strict;
use warnings;


my %j;
my $co = 0;
foreach my $sa(glob "*stastic.fasta")
{
    print $sa."\n";
            my $filename = $sa;
        open(my $fh, "<", $filename) or die "Cannot open file: $!";
        open IO,">".$sa.".sort.fasta";
        my %hash;
        while (<$fh>) {
            chomp;
            my ($key, $value) = split("\t");
            $hash{$key} = $value;
            
            if($j{$key})
            {
                $j{$key} = $j{$key} + $value;
            }
            else
            {
                $j{$key} = $value;
            }
            $co = $co + $value;
        }
        close($fh);

        # 按照第二列的值进行排序并输出
        my $cou = 0;
        foreach my $key (sort { $hash{$b} <=> $hash{$a} } keys %hash) {
            #print IO "$key\t$hash{$key}\n";

            next if(length($key) < 19);

            print IO ">count$hash{$key}.$cou.\n$key\n";
            $cou++;
            #last if($cou > 10);
            last if($hash{$key} < 50);
        }
        close IO;



}

open IO,">"."all.merge.sort.fasta";
my $cou = 1;
foreach my $key (sort { $j{$b} <=> $j{$a} } keys %j) {
    #print IO "$key\t$hash{$key}\n";
    my $hl = sprintf('%.05f',($j{$key} / $co) * 100);
    print IO ">count.$j{$key}\_$hl\_$cou\n$key\n";
    $cou++;
    #last if($cou > 100);
    last if($j{$key} < 100);
}
close IO;


# foreach my $sa(glob "*.fq")
# {
#     print $sa."\n";
#     my $filename = $sa;
#     # my $ou= "temp.".$filename;
#     # my $finou = $filename.".1000.fasta";
#     # system("seqtk sample -s100 $filename 1000 > $ou");
#     # system("seqtk seq -a $ou > $finou");
#     # system("rm $ou");
#     ####
#     #system("seqtk seq -a $filename > ./blast/$filename");
#     system("blastn -query ./blast/$filename -out ./blast/$filename.blast -db ./blast/selc -outfmt 6");

# }