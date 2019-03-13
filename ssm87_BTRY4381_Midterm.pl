#!/usr/bin/perl
use warnings; use strict;

# Samuel Moijueh
# BTRY 4381 Bioinformatics Programming
# April 5, 2013
# MIDTERM

open(F1, "gene_expression.txt") or die $!;
open(PROBLEM1, "+>ssm87_BTRY4381_Midterm_problem1.txt") or die$!;
open(PROBLEM2, "+>ssm87_BTRY4381_Midterm_problem2.txt") or die$!;

# Hash of Arrays where the keys are the genes and the value is an array that
# of the 114 gene expression values at different conditions
my %gene_centers = ();
my %all_genes = ();
my $go = 1;
my $k = 10;

############## K-means Clustering STEP 1 ############################
# parse the data file; create two hash of arrays:
# one containing all the genes and their gene expression profile
# the second containing the initial centers
# keys represent the gene and the value is an array of its normalized
# gene expression values

while (<F1>){
    $_=~ s/\r//;
    chomp;
    my @cols = split(/\t/, $_);
    push( @{$all_genes{$cols[0]}}, @cols[1..$#cols]);
    
    # HoA: assign the first 10 genes in the file as initial centers
    if ($go <= 10) {
	push(@{$gene_centers{$cols[0]}}, @cols[1..$#cols]);
	$go++;
    }
}

close F1;

print "\t\tBeginning K-means clustering algorithm. Please wait...\n";

############## K-means Clustering STEP 2 ############################
# loop through all the genes and assign each one to the cluster with the
# closest center
# key is the center
# value is an array of genes that cluster to that point
# loop through the keys from lowest to highest

my $start = time();
my %clusters;  # a HoA that will contain the gene centroids and the genes that converged to that point
foreach my $key ( sort {$a <=> $b} (keys %all_genes) ) {
    push( @{ $clusters{pearson_coefficient($all_genes{$key}, \%gene_centers)}}, $key );
}

my $end = time();
my $elapsed = $end - $start;

print "\tStep 2 Benchmark: $elapsed seconds\n";

%gene_centers = ();  # empty this hash
$go = 1; # reinitilize

############## K-means Clustering STEP 3 ############################
# Compute the new centers of the data points
foreach my $center ( keys %clusters ) {
    # returns a 114-dimensional vector of the gene expression values
    # for the new gene center
    # pass in an array_ref to the genes in $center
    push(@{$gene_centers{"C$go"}}, compute_centers($clusters{$center}));
    $go++;  # iteratively generate new key for %gene_centers
}

# ASSERT: $gene_centers is a HoA that contains the gene centers and their gene expression values

my $iteration = 1;
print "Now on iteration: $iteration of k = $k means clustering\n";
############## K-means Clustering STEP 4 ############################
# Recursive, Iterative Approach
# Do Steps 2 and 3 iteratively until no object changes cluster
# in other words, until the array values in %cluster remain the same

my $temp_ref = [];

# record the size of each cluster
foreach my $key ( keys %clusters ) {
    push ( @{ $temp_ref },  scalar( @{$clusters{$key}} ) );
}

$iteration++;
# check if the number of genes per cluster has changed
while (1){
    unless (repeat($k, $temp_ref)){
	last;
    } else {
	$iteration++;
	$temp_ref = []; # reinitialize array
	foreach my $key ( keys %clusters ) {
	    my $size = scalar( @{$clusters{$key}} );
	    push(@{ $temp_ref }, $size );    
	}
    }
}

print "DONE! Data points have convergenced in $iteration Iterations!\n";
printf PROBLEM1 "%-11s %-11s %-11s\n",
    "Gene Center |", "Number of Points in Cluster |", "Within Cluster Sum of Squares";

my $WCSS = 0;

# for each cluster in the data set
foreach my $key ( sort { @{$clusters{$b}} <=> @{$clusters{$a}}} keys %clusters ) {
    my $size = scalar( @{$clusters{$key}} );
    $WCSS = within_sos(\@{$clusters{$key}}, $key);
    printf PROBLEM1 "%-13s %-30d %-11s\n",$key, $size, $WCSS;
}

print PROBLEM1 "\nDONE! Data points have convergenced in $iteration Iterations!\n";
print "Congratulations! You've completed Problem 1 of the Midterm!\n";

##################################################################################################

print "Now starting Problem 2 of the Midterm\n";

################## Problem 2 #####################################################################

my @k_clusters = (20, 40,60,80,100,150,200);
printf PROBLEM2 "Number of k Clusters\tWSS\tTSS\tBSS\tPercentage of variance explained\tNumber of Iterations\n";

foreach my $k (@k_clusters){
    # Step 1 of K-means
    open(F1, "gene_expression.txt") or die $!;
    $go = 1;
    %gene_centers = ();    
    while (<F1>){
	$_=~ s/\r//;
	chomp;
	my @cols = split(/\t/, $_);

	# HoA: assign the first k genes in the file as initial centers
	# must attack for genes with same gene profile
	if ($go <= $k) {
	    push(@{$gene_centers{$cols[0]}}, @cols[1..$#cols]);
	    $go++;
	}
	delete $gene_centers{'100008586'}; # removes gene 5

	if ($k == 200)  {
	    delete $gene_centers{'100129484'}; # removes gene 170
	}
    }
    close F1;

    %clusters = ();

    # Step 2 of K-means
    foreach my $key ( sort {$a <=> $b} (keys %all_genes) ) {
	push(@{$clusters{pearson_coefficient($all_genes{$key}, \%gene_centers)}}, $key);
    }

    %gene_centers = ();  # empty out this hash
    $go = 1; # reinitilize
    
    # Step 3: Compute the new centers of the data points
    foreach my $center ( keys %clusters ) {
	# returns a 114-dimensional vector of the gene expression values
	# for the new gene center
	push(@{$gene_centers{"C$go"}}, compute_centers($clusters{$center}, \%all_genes));
	$go++;  # iteratively generate new key for %gene_centers
    }
    
    $iteration = 1;
    print "Now on iteration: $iteration of k=$k means clustering\n";

    # Step 4
    $temp_ref = [];
    
    # determine the initial size of the each cluster
    foreach my $key ( keys %clusters ) {
	my $size = scalar( @{$clusters{$key}} );
	push(@{ $temp_ref }, $size );    
    }
    
    $iteration++;
    
    # check if the number of genes per cluster has changed
    while (1){
	unless (repeat($k, $temp_ref)){  
	    last;
	} else { # the genes have not converged
	    $iteration++;
	    $temp_ref = []; # reinitialize array
	    # determine the size of each of the new clusters
	    foreach my $key ( keys %clusters ) {
		my $size = scalar( @{$clusters{$key}} );
		push(@{ $temp_ref }, $size );    
	    }
	}
    }
    
    print "DONE! Data points have convergenced in $iteration Iterations!\n";
    
    my $TSS = total_sos(\%all_genes);
    my $WSS = 0;
    
    # loops through the centroid of each cluster
    foreach my $key ( sort { @{$clusters{$b}} <=> @{$clusters{$a}}} keys %clusters ) {
	my $size = scalar( @{$clusters{$key}} );
	# sum up the WSS for each centroid
	$WSS += within_sos(\@{$clusters{$key}}, $key);
    }
    
    # finally calculate the percentage of variarance explained
    my $BSS = ($TSS - $WSS);
    my $percent_explained = $BSS/$TSS;
    print PROBLEM2 "$k\t$WSS\t$TSS\t$BSS\t$percent_explained\t$iteration\n";
    my $new_k;
    if ($k < 100 ){  
	$new_k = $k + 20;
    } else {
	$new_k = $k + 50;
    }
    
    if ($k == 200){
	print "\nWoohoo!!! Completed clustering for k = $k means. Congratulations, you've completed the Midterm!\n";
	print "Go ahead, take a look at the output files. You've earned it!!!!\n";
    } else{
	print "\nCompleted clustering for k = $k means. Now performing the first TWO iterations for k = $new_k means.\n";
	print "This may take a while. Please wait.\n";
    }
}

############################## Subroutines #######################################

# compute the within cluster sum of squares for a single cluster
sub within_sos{
    # takes two params: (1) an array_ref of the genes in a cluster and
    # (2) the gene center which we use to immediately determine x-bar
    my $array_ref = $_[0];  # contains a list of all the genes in a cluster
    my @x_bar = @{$gene_centers{$_[1]}};
    my $sum;

    # for each gene in the cluster
    foreach my $gene (@{ $array_ref }){
	my @x_i = @{$all_genes{$gene}};  # gene profile for this gene
	my @out = map { $x_i[$_] - $x_bar[$_] } 0..$#x_i;
	$sum += dProduct(\@out, \@out);
    }
    return $sum; # return the sum of squares
}

# returns the total sum of squares for the entire data set
sub total_sos{
    # takes 1 param: a hash ref of all the gene profiles
    my $hash_ref = $_[0];  # %{$all_genes}
    my @arrays;
    
    # create an AoA of all the gene expression profiles in the file
    foreach (values %{ $hash_ref } ){
	push(@arrays, $_);
    }
    
    my $SIZE = scalar(@{ $arrays[0] });
     
    my @x_bar;
    foreach my $aref (@arrays) {
	for (my $i=0; $i<$SIZE; $i++) {
	    $x_bar[$i] += $aref->[$i];
	}
    }
    
    @x_bar = map { $_ / @arrays } @x_bar;
    
    my $running_sum;
    foreach my $aref (@arrays){
	my @a1 = @{ $aref };
	my @out = map { $a1[$_] - $x_bar[$_] } 0 .. $#a1;
	$running_sum += dProduct(\@out, \@out);
    }
    
    return $running_sum;
}


# return the dot product of two arrays
sub dProduct {
    my $x = $_[0];
    my $y = $_[1];
    my $sum;
    my $ct1 = $#{$x} + 1;   # items in $x

    for (my $i=0;$i<$ct1;$i++) {
	$sum += $$x[$i] * $$y[$i];
    }
    return $sum; # return the dot product
}

my $previous_cluster_ref = {};

# return true if the number of genes in each cluster has changed
sub repeat{
    # takes 2 params: k, the number of iterations and the array_ref
    $start = time();
    my ($k, $temp_ref) = @_;
    my $array_ref = [];

    print "Now on iteration: $iteration of k=$k means clustering\n";
    print "\tComputing... Please wait\n\n";

    # store the original set of clusters
    $previous_cluster_ref = \%clusters;

    %clusters = ();  #initialize empty out this hash

    # the new set of clusters
    foreach my $key ( sort {$a <=> $b} (keys %all_genes) ) {
	push(@{$clusters{pearson_coefficient($all_genes{$key}, \%gene_centers)}}, $key);
    }

    my $end = time();
    my $elapsed = $end - $start;

    print "\t\tCompleted iteration of algorithm in $elapsed seconds\n";

    # determine the new size of each cluster
    foreach my $key ( keys %clusters ) {
	push(@{ $array_ref }, scalar( @{$clusters{$key}}));
    }

    print "Size of each cluster is now: (@{ $array_ref })\n";
    print "Size of each cluster in previous iteration: (@{ $temp_ref })\n";

    # check if the size of each cluster is the same
    if (!([sort @{ $temp_ref }] ~~ [sort @{ $array_ref }])) {
	# if these two arrays are not identical, then there is fluctuation 
	# therefore repeat steps 2 and 3.

	# enter this loop if the cluster sizes have changed

	%gene_centers = ();  # empty out this hash
	$go = 1; # reinitilize

	# Step 3: Compute the new centers of the data points
	foreach my $center ( keys %clusters ) {
	    # returns a 114-dimensional vector of the gene expression values for the new gene center
	    push(@{$gene_centers{"C$go"}}, compute_centers($clusters{$center}, \%all_genes));
	    $go++;  # iteratively generate new key for %gene_centers
	}

	return 1;  # return true indicating that there was fluctuation(!)

    } else {
	# check that the same genes have in fact converged and have not switched clusters between clusters
	print "\nNow that the size each of each cluster has remained constant between iterations, determine that no genes are switching between clusters\n";
	
	# if the size of the cluster hasn't changed then
	# check if the cluster members are the same
	my ($s1, $s2);
	my %previous_cluster = %{ $previous_cluster_ref };
	
	foreach ( keys %gene_centers ){
	    # Compares two lines exactly
	    $s1 = join(',', @{$clusters{$_}});
	    $s2 = join(',', @{$previous_cluster{$_}});
	}

	if ($s1 eq $s2){
	    return 0; # false: the points have converged
	} else {
	    return 1; # true: repeat steps 2 and 3
	}
    }
}

# compute the new centers. returns an N-dimensional array of gene expression profile
sub compute_centers{
    # takes 1 param: an array_ref of all the genes in a cluster
    my $array_ref = $_[0];
    
    my $center = $array_ref->[0];
    my ($AoA_ref, $sums) = [];
    
    if (scalar@{$array_ref} == 1){
	return @{$all_genes{$center}};
    }
    
    # remove the cluster center from the cluster; usually it is
    # not part of the calculation for the new cluster center
    #@{ $array_ref } = grep { $_ != $center } @{$array_ref};
    
    # push the gene profile for every gene onto an array of arrays
    foreach (@{ $array_ref }){
	push ( @{ $AoA_ref }, \@{$all_genes{$_}} );
    }
    
    my $size_of_cluster = @{$AoA_ref};
    
    foreach my $column (0..$#{@{$AoA_ref}[0]}) {
	my $sum = 0;
	foreach my $aref (@{ $AoA_ref }){
	    $sum += $aref->[$column];
	}
	push ( @{ $sums } ,$sum/$size_of_cluster );
    }
    return @{$sums};
}

# calculates 10 Pearson correlation coefficients (PCC) between
# the gene of interest and each of the k = 10 clusters centers
# returns the gene center in %gene_center with the highest PCC
sub pearson_coefficient{
    # takes 2 params: array_ref of gene profile for single gene and an HoA_ref of centroids
    my ($gene_exp_ref, $centroids_ref) = @_;
    my ($numerator, $denominator) = 0;
    
    my ($prod_ref, $diff_x_ref, $diff_y_ref, $x_sq_ref, $y_sq_ref) = [];
    my %gene_center_pcc;   # where diff_x is gene of interest
    
    my $gene_exp_average = mean($gene_exp_ref);  # the gene profile average
    
    for my $gene_exp (@{$gene_exp_ref}) {  # loop through every gene expression value in profile
	push(@{ $diff_x_ref }, ($gene_exp - $gene_exp_average));
    }
    
    # possible bottleneck
    for my $centroid_gene_exp_ref (values %{$centroids_ref}){
	$diff_y_ref = [];  # initilize back to empty array
	my $mean = mean($centroid_gene_exp_ref);    # cache the average of the of centroid
	
	for my $index (@{$centroid_gene_exp_ref}) {
	    push(@{ $diff_y_ref }, ($index - $mean));
	}
	
	@{ $prod_ref } = map { @{ $diff_x_ref }[$_] * @{ $diff_y_ref }[$_] } 0..$#{ $diff_x_ref };
	
	$numerator = sum($prod_ref);
	
	@{ $x_sq_ref } = map {$_*$_}@$diff_x_ref;
	@{ $y_sq_ref } = map {$_*$_}@$diff_y_ref;
	
	$denominator = sqrt(sum($x_sq_ref)) * sqrt(sum($y_sq_ref));
	
	my $r = $numerator/$denominator;
	
	# matches a gene with it's gene expression profile
	my ($center) = grep { @{$gene_centers{$_}} ~~ @$centroid_gene_exp_ref } keys %gene_centers;
	$gene_center_pcc{$center} = $r; #store the PCC value for that gene
    }
    
    #return the centroid with the highest PCC value
    return (sort {$gene_center_pcc{$b} <=> $gene_center_pcc{$a}}
	    keys %gene_center_pcc)[0];
}

# sum returns the summmation of all the elements in the array.
# Dividing by @_ (which in scalar context gives the # of elements) returns the average
sub mean{
    # takes 1 param: an array_ref of the gene expression profile of a particular gene
    return sum($_[0])/scalar(@{$_[0]});
}

# returns the summation of all the elements in a numeric array
sub sum{
    # takes 1 param: an array_ref of the gene expression profile of a particular gene
    my $sum = 0;
    for ( @{$_[0]} ) {
	$sum += $_;
    }
    return $sum;
}
