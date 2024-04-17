#usr/bin/perl
use strict;
#use PerlIO::gzip;
#TGGGAGCTCTCCGGATCCAACCCGGGGTACCGAATTCCTCGAG
my  %gg;
my @i;
#foreach my $a (glob "*.fq")
#foreach my $a (glob "*45*.fq.gz")
foreach my $a (glob "*.fq.gz")
{
		#next if($a !~ m/35/);
		#print $a."\n";
	
	#open F,"<:gzip(gzip)", $a or die ("can not open $_\n");

	open(F, "gunzip -c $a |");
	#open F,"$a";

	my ($fff) = (split/\./,$a)[0];
	print $fff."\n";
	open S,">".$fff.".fasta";
	open SQ,">".$fff.".revese.fasta";
	open SAS,">".$fff.".forwarever.fasta";
	open K,">".$fff."stastic.fasta";
    %gg = ();
	while(<F>)
	 {

	 	my $se =<F>;
	 	my $thir =<F>;
	 	my $fo =<F>;
	 	#my ($aa) = $se =~ m/aCTCGAGGAATTCGGTACCCCGGGT(\w+)TGGATCCGGAGAGCTCCCAacgcgt/i;
	 	my ($aa) = $se =~ m/CCCCGGGT(\w+)TGGATCCGGA/i;

	 	my ($aa2) = $se =~ m/TCCGGATCCA(\w+)ACCCGGGG/i;
	 	#my ($qa) = $se =~ m/acgcgtTGGGAGCTCTCCGGATCCA(\w+)ACCCGGGGTACCGAATTCCTCGAGt/i;
            s/[\n\r]//;

	 	 #if($aa=~ m/[ATGC]/ && length($aa) <= 50 && length($aa) > 5)
	 	 #if($aa=~ m/[ATGC]/ && length($aa) == 45)

	 	if($aa=~ m/[ATGC]/)
	 	  {
	 	  	 print S ">".$_."\n".$aa."\n";
	 	  	 print SAS ">".$_."\n".$aa."\n";
	 	  	$gg{$aa}++;
	 	  }
	 	if($aa2=~ m/[ATGC]/)
	 	{
	 		my $reverse_complement_sequence = reverse_complement($aa2);
	 		
	 		print SQ ">".$_."\n".$reverse_complement_sequence."\n";
	 		print SAS ">".$_."\n".$reverse_complement_sequence."\n";
	 		$gg{$reverse_complement_sequence}++;
	 	}

	 }

		foreach my $key (sort {$gg{$b}<=>$gg{$a}} keys %gg) {
		    print K "$key\t$gg{$key}\n";
		}
	close F;
	close S;
	close SQ;
	close SAS;
	 close K;

	 # open KA,">".$fff."stastic.fasta.txt";
		# my $filename = $fff."stastic.fasta";
		# open(my $fh, "<", $filename) or die "Cannot open file: $!";
		# my %hash;
		# while (<$fh>) {
		#     chomp;
		#     my ($key, $value) = split("\t");
		#     $hash{$key} = $value;
		# }
		# close($fh);

		# # 按照第二列的值进行排序并输出
		# foreach my $key (sort { $hash{$b} <=> $hash{$a} } keys %hash) {
		#     print KA "$key\t$hash{$key}\n";
		# }
		# close KA;


}

sub reverse_complement {
    my ($sequence) = @_;

    # 使用正则表达式将序列反转
    my $reverse_sequence = reverse($sequence);

    # 使用替换操作将互补碱基替换
    $reverse_sequence =~ tr/ACGTacgt/TGCAtgca/;

    return $reverse_sequence;
}

