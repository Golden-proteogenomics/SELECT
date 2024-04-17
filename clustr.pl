use Algorithm::KMeans;
use strict;
use warnings;

# DNA sequences data


my @data;

open IN,"SELEX-13-LFJ14896_L3_1stastic.fasta";
while (<IN>) {
    s/[\n\r]//;
    my @r = split/\t/,$_;
    if($r[-1] > 100 & length($r[0]) == 45)
    {
        push @data,$r[0];
    }
    # body...
}
close IN;

my $max_k = 500; # Try different values of k
my @inertia_values;

for my $k (1..$max_k) {
    my $kmeans = Algorithm::KMeans->new(
        data => \@data,
        k    => $k,
    );
    $kmeans->kmeans;
    push @inertia_values, $kmeans->inertia;
}

# Find the "elbow" point using the elbow method
my $elbow_point = find_elbow_point(\@inertia_values);

print "The suggested k value is $elbow_point\n";

sub find_elbow_point {
    my ($inertia_values) = @_;
    my @diffs;
    for my $i (1..$#{$inertia_values}) {
        push @diffs, $inertia_values->[$i-1] - $inertia_values->[$i];
    }
    for my $i (1..$#diffs) {
        if ($diffs[$i] < 0) {
            return $i + 1;
        }
    }
    return 1;
}
