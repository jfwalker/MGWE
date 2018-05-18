use Data::Dumper;

#This subroutine gets the different gene lengths
#This passes all tests
sub GetParts {
	
	@temp_array = ();
	@temp_array = @_;
	$PartFile = $temp_array[0];
	$IsItATest = $temp_array[1];
	$CountOfGenes = 0;
	open(Parts, "$PartFile")||die "Can't Find Parts file";
	while($line = <Parts>){

		chomp $line;
		$location = ($line =~ m/.*? = (.*)/)[0];
		if($location eq ""){
			print "Error program killed check outfile\n";
			print StatsOut "!!!!!!!!!!!!!!Error!!!!!!!!!!!!!!\n";
			print StatsOut "Need to fix parts file so it matches: ";
			print StatsOut "name_of_gene = start-stop\n";
			die;
		}
		#This can be used to specify only a few genes
		#CHANGE AFTER TESTING!!!!!!!!!!!!!!!!
		if($IsItATest eq "True"){
			if($CountOfGenes < 2){
				
				push @loc, $location;
		
			}else{}
		}else{
			push @loc, $location;
		}
		$CountOfGenes++;		
	}
	return @loc;
}

#This subroutine gets the bipartition information
#This passes all tests
sub GetBipartitions {
	
	$conflict = ""; local *conflict = $_[0];
	$clade = ""; @array = ();
	#This array will hold all relationships being examined
	push @array, $conflict;
	open(log_file, "bp.log") || die "No log File";
	while($line = <log_file>){
		
		chomp $line;
		if($line =~ /^CLADE/){
			
			$test = "false";
			$bipart = ($line =~ m/CLADE: (.*)/)[0];

			@array1 = split " ", $bipart;
			@array2 = split ",", $conflict;
			
			#print "ARRAY1\n";
			#print Dumper(\@array1);
			#print "ARRAY2\n";
			#print Dumper(\@array2);
			#this is really hacky and stupid but I can't find a better
			#way to match bipartitions and conflict
			if($#array2 == $#array1){
				#make an empty hash
				%comp = ();
				#use values in array1 as undefined
				@comp{@array1} = undef;
				#if its in array2 delete is
				delete @comp{@array2};
				#get whats left to see if it's all gone
				$comp_size = keys %comp;
				if($comp_size == 0){
					$test = "true";
				}
				
			}		
		}elsif($test eq "true"){
			if($line =~ /\tCOUNT/){
				if($line !~ /ICA/){
					#print "$line\n";
					$clade = ($line =~ m/ \t (.*?) \tCOUNT.*/)[0];
					#print "HERE: $clade\n";
					$clade =~ s/ /,/g;
					$clade =~ s/,\t,(.*?),$/$1/;
					push @array, $clade;
				}
			}elsif($line =~ /FREQ/){
				
				#$ICA = ($line =~ m/.*?ICA:\t(.*?)\t.*/)[0];
				#print StatsOut "The ICA for unique trees in this relationship is: $ICA\n";
			}
		}
	}
	if($verbose eq "True"){
		print StatsOut "##################Your Conflicts################\n";
		foreach $i (0..$#array){
		
			print StatsOut "Conflict $i: $array[$i]\n";
		}
	}
	return @array;		
}

#Get the bipartitions that each tree is associated with
#everything seems functional
sub TreeBiparts {
	
	@conflict = (); local *conflict = $_[0];
	$tree_count = 0; $clade = "";
	%TREE_HASH = (); $count = 0;
	@Conf = (); @Conf_And_Tree = ();
	if($verbose eq "True"){
		print StatsOut "###############Your Trees And their Conflicts##########\n";
	}
	#opens the tree file from the analysis
	open(trees, "Unique.tre")||die "No tree file\n";
	while($line = <trees>){
		
		chomp $line;
		open(tempout, ">temptesttre");
		print tempout "$line\n";
		system("$pxbp -t temptesttre > temp.log");
		open(temp, "temp.log");
		@Conf = ();
		while($temp = <temp>){
			chomp $temp;
			if($temp =~ /^CLADE/){
				$clade = ($temp =~ m/CLADE: (.*?)\t.*/)[0];
				#print "Here: $clade\n";
				@array1 = ();
				@array1 = split " ", $clade;
				foreach $i (0..$#conflict){
					
					@array2 = ();
					@array2 = split ",", $conflict[$i];
					if($#array1 == $#array2){
						%comp = ();
						@comp{@array1} = undef;
						delete @comp{@array2};
						$comp_size = keys %comp;
						if($comp_size == 0){
							if($verbose eq "True"){
								print StatsOut "Tree: $count Conflict: $i $conflict[$i]\n";
							}
							push @Conf, $i;
						}
					}
					
				}
			}
		}
		push @Conf_And_Tree, [@Conf];
		$count++;
	}
	return @Conf_And_Tree;
	system("rm temp.log temptesttre");
}

#Read in the Supermatrix
#This seems good
sub ReadSuperMatrix {
	
	
	$count = 0;
	$SuperName = ""; 
    local *SuperName = $_[0];
    %FastaHash = (); $name = "";
    $seq = "";
    #print "$SuperName\n";
    open(Supermatrix, "$SuperName");
	while($line = <Supermatrix>){
		
		chomp $line;
		if($line =~ /^>/){
			if($count != 0){
				
				$FastaHash{$name} = $seq;
			}
			$name = ($line =~ m/>(.*)/)[0];
			$seq = "";
		}else{
			$seq .= $line;
		}
		$count++;
	}
	$FastaHash{$name} = $seq;
	return %FastaHash;
	
	
}
sub top_o_da_conflict {
	
	$Likelihoods = ""; local *Likelihoods = $_[0];
	#Tree is the number and then theres the conflicts it matches
	@TreeAndCon = (); local *TreeAndCon = $_[1];
	@array = (); %HASH = ();
	@array = split " ", $Likelihoods;
	for $i (0..$#TreeAndCon){
		
		$ref = $TreeAndCon[$i];
		for $j (0..$#{$ref}){
			
			#print StatsOut "Tree $i = Conflict $TreeAndCon[$i][$j]\n";
			$HASH{$TreeAndCon[$i][$j]} .= "$array[$i]:$i,"
			
		}
		
	}
	$hash_size = keys %HASH;
	@array = (); @trees = (); @likes = ();
	$temp_like = 0; $temp_tree = "";
	foreach $i (0..($hash_size-1)){
		
		#print "$i: $HASH{$i}\n";
		@array = split ",", $HASH{$i};
		$temp_like = -9999999999999.999;
		foreach $j (0..$#array){
		
			($like, $tree) = split ":", $array[$j], 2;
			if($temp_like < $like){
				$temp_like = $like;
				$temp_tree = $tree;
			}
		}
		push @likes, $temp_like;
		push @trees, $temp_tree;
		#print "$i: $temp_tree\t$temp_like\n";
	}
	return (\@likes, \@trees);
}



open(Configure, "$ARGV[0]")||die "Please See Configure File\n";
while($line = <Configure>){
	
	if($line =~ /^pxrmt:/){
		$pxrmt = ($line =~ /.*?: (.*)/)[0];	
	}elsif($line =~ /^pxbp:/){
		$pxbp = ($line =~ /.*?: (.*)/)[0];	
	}elsif($line =~ /^raxml:/){
		$raxml = ($line =~ /.*?: (.*)/)[0];
	}elsif($line =~ /^Species:/){
		$conflicting_node = ($line =~ /.*?: (.*)/)[0];
	}elsif($line =~ /^outfile:/){
		$outfile = ($line =~ /.*?: (.*)/)[0];
	}elsif($line =~ /^Supermatrix:/){
		$SuperMatrix = ($line =~ /.*?: (.*)/)[0]; 
	}elsif($line =~ /^PartsFile:/){
		$PartFile = ($line =~ /.*?: (.*)/)[0]; 
	}elsif($line =~ /^Set:/){
		$TreeFile = ($line =~ /.*?: (.*)/)[0]; 
	}elsif($line =~ /^Threads:/){
		$threads = ($line =~ /.*?: (.*)/)[0]; 
	}elsif($line =~ /^Test:/){
		$IsItATest = ($line =~ /.*?: (.*)/)[0];
	}elsif($line =~ /^Verbose:/){
		$verbose  = ($line =~ /.*?: (.*)/)[0];	
	}elsif($line =~ /^Folder:/){
		$folder = ($line =~ /.*?: (.*)/)[0];	
	}elsif($line =~ /^secret:/){
		$secret = ($line =~ /.*?: (.*)/)[0];
	}elsif($line =~ /^Topologies:/){
		$Topos = ($line =~ /.*?: (.*)/)[0];
	}
}
open(StatsOut, ">$outfile")||die "In program give a name of the outfile\n";
print StatsOut "Information for the MGWE Analysis\n";
print StatsOut "################################################################\n";
print StatsOut "If any of this is wrong the analysis won't work, double check!!!!\n";
print StatsOut "You have set pxrmt to be in the path: $pxrmt\n";
print StatsOut "You have set pxbp to be in the path: $pxbp\n";
print StatsOut "You have set raxml to be in the path: $raxml\n";
print StatsOut "The species in the relationship you are looking at are: $conflicting_node\n";
print StatsOut "You are testing this many Topos: $Topos\n";
print StatsOut "Your supermatrix is file called: $SuperMatrix\n";
print StatsOut "Your partition file is called: $PartFile\n";
print StatsOut "Your tree set file is (newick hopefully?): $TreeFile\n";
print StatsOut "Running a test?: $IsItATest\n";
print StatsOut "You're output folder is: $folder\n";
print StatsOut "Running it in verbose?: $verbose\n";
print StatsOut "################################################################\n";
if($secret ne "True"){
	system("rm -Rf $folder && mkdir $folder");
}
print StatsOut "You're Edge Info is in bp.log\nYou're Unique Trees are in Unique.tre\n";

system("pxrr -u -t $TreeFile -o trees.unroot");
system("$pxbp -t trees.unroot -u | grep \"\(\" > Unique.tre");
system("$pxbp -t Unique.tre -v > bp.log");

#loc will store locations
@loc = ();
@loc = GetParts($PartFile, $IsItATest);

#CountOfGenes is total to be analyzed
$CountOfGenes = 0;
$CountOfGenes = ($#loc+1);

print StatsOut "You're total genes: $CountOfGenes\n";

#ConflictHash will have conflicting relationships
#element 0 will contain the bipart of interest the
#rest will be the ones that conflict with it
@Conflict = ();
@Conflict = GetBipartitions(\$conflicting_node);


#Array of Arrays where tree is the position and
#The array in that position contains the conflicts
#it matches
@Tree_Conflict = ();
@Tree_Conflict = TreeBiparts(\@Conflict);
print "Trees have been processed now Beginning Gene Likelihoods\n";
if($verbose eq "True"){
	print StatsOut "################################################################\n";
	print StatsOut "#HERE ARE THE TAXA AND MISSING ONES FROM EACH GENE:\n";
}
%FastaHash = (); $count = 0; %TotalTaxaHash = ();
%FastaHash = ReadSuperMatrix(\$SuperMatrix);
$all_seqs = 0;
$all_seqs = keys %FastaHash;


$sitecount = 0; $genecount = 0;
@like_sum = (); @tree_sum = ();
@parameter_sum = (); @Tree_total = ();
@gene_lengths = (); @gene_comp = ();
foreach $i (0..$#loc){
	#This block creates a temporary Fasta for the gene
	open(out, ">temp.fa");
	($start,$stop) = split "-", $loc[$i], 2;
	$dif = $stop - $start + 1;
	$to_remove = ""; $seq_count = 0;
	for $keys (sort keys %FastaHash){
		
		$seq = substr $FastaHash{$keys}, ($start-1), $dif;
		#check for sequences made of all missing data
		if($AMINO eq "TRUE"){
			$missing = $seq =~ tr/\-X/\-X/;
		}else{
			$missing = $seq =~ tr/N\-X/N\-X/;
		}
		if($missing != $dif){
			print out ">$keys\n$seq\n";
			#All the seqs going into the likelihood calc, helps with
			#Number of parameters
			$seq_count++;
		}else{
			#print out ">$keys\n$seq\n"; #Temp
			$to_remove .= "$keys,"
		}
	}
	#For parameters stores number of taxa
	if($verbose eq "True"){
		print StatsOut "Gene_$i: Has $seq_count in it\n";
	}
	push @parameter_sum, $seq_count;
	#Anything that was missing from the gene gets eliminated here
	if($verbose eq "True"){
		print StatsOut "These were removed from likelihood calc due to all missing data: $to_remove\n";
	}
	
	#Removes taxa from tree set to not interfere with likelihood calc
	if($to_remove ne ""){
		system("$pxrmt -t Unique.tre -n $to_remove > TempTree.tre");
	}else{
		system("cp Unique.tre TempTree.tre");
	}

	#Deletes RAxML Files
	system("rm RAxML_*");
	
	#calculate the SSLL for each tree
	print "(☞ﾟヮﾟ)☞ Processing Gene $i\n";
	system("$raxml -f g -T $threads -s temp.fa -m GTRGAMMA -z TempTree.tre -n EX_SSLL | grep \"Tree \" | grep \":\"");
	
	#read in SSLL and make a new file of them
	$t_count = 0;
	$gene1 = "";
	$ref_t = ""; $ref_l = "";
	@best_trees = (); @best_likes = ();
	open(RAXML, "RAxML_perSiteLLs.EX_SSLL");
	while($line = <RAXML>){
		
		if($t_count != 0){	
			chomp $line;
			($one,$two) = split " ", $line, 2;

			@temp = split " ", $two;
			$held = 0;
			foreach $i (0..$#temp){
				$held += $temp[$i];
			}

			$gene1 .= "$held ";
			if($t_count <= $Topos){
				@gene_comp[$t_count - 1] += $held;
			}
			$sitecount++;
		}
		$t_count++;
	}
	($ref_l, $ref_t) = top_o_da_conflict(\$gene1, \@Tree_Conflict);
	@best_trees = @$ref_t; @best_likes = @$ref_l;
	$very_best = -9999999999999999.9999999;
	$temp_tre = "";
	foreach $j (0..$#best_likes){
		
		$like_sum[$j] += $best_likes[$j];	
		$tree_sum[$j] .= "t$best_trees[$j],";
		if($very_best < $best_likes[$j]){
			
			$very_best = $best_likes[$j];
			$temp_tre = $j;
			
		}
		
	}
	

	push @gene_lengths, $dif;
	push @Tree_total, $temp_tre;
}


#This whole monstrosity of a section of code is about calculating
#The likelihood for the trees if you don't parameterize branch lengths
@MatrixLikes = (); $likelihood = "";
if($secret eq "True"){
	print "Secret Mode is chosen, not calculating SSLL's for matrix\n";
	open(File, "$folder/MatrixNoBrInfo.SSLL")||die "Did you delete the file MatrixNoBrInfo.SSLL?\nCannot run in secret mode without that file\n";
	while($line = <File>){
			
		chomp $line;
		if ($line =~ /Tree/){
			if($line =~ /Tree.*?:/){
				
				$likelihood = ($line =~ m/Tree.*?: (.*)/)[0];
				push @MatrixLikes, $likelihood;
			}
		}
			
	}
}else{
	if($Topos != 0){
		print "You have chosen to test $Topos topologies\n";
		print "Performing that analysis now, can be timely\n";
		print "( ﾟдﾟ)\n";
		$count = 0;
		open(trout, ">Questionable.tre");
		open(Trees, "Unique.tre")||die "No unique trees";
		while($line = <Trees>){
	
			chomp $line;
			if($count < $Topos){
			
				print trout "$line\n";
			}
			$count++;
		}
		system("$raxml -f g -T $threads -s $SuperMatrix -q $PartFile -m GTRGAMMA -z Questionable.tre -n Topologies_SSLL | grep \"Tree \" | grep \":\"");
		system("mv RAxML_perSiteLLs.Topologies_SSLL MatrixNoBrSSLLs.SSLL");
		system("mv RAxML_info.Topologies_SSLL MatrixNoBrInfo.SSLL");
		print "Done! ( ﾟヮﾟ)\n";
		open(File, "MatrixNoBrInfo.SSLL");
		while($line = <File>){
			
			chomp $line;
			if ($line =~ /Tree.*?:/){
				if($line =~ /:/){
				
					$likelihood = ($line =~ m/Tree.*?: (.*)/)[0];
					push @MatrixLikes, $likelihood;
				}
			}
			
		}
	}else{
		print "#################################################\n";
		print "Interesting you have Zero topologies?\n";
		print "This will give you some to look further into\n";
		print "But it will do it simply by finding common\n";
		print "Bipartitions so maybe worth re-running with that\n";
		print "In mind\n";
		print "#################################################\n";

	}
}
#End of the monstrosity of code about branch lengths


#This grabs the number of parameters, will be different between
#supermatrix and edges because supermatrix uses only one
#set of branch lengths
print "#################################################\n";
print "Went through data, now processing final bits\n";


#Parameters are 6 for GTR+G, 3 for empiracle sites and (2n-3) for the branch lengths
#Topology is not really a param because it creates the point in likelihood space
#this is what allows for Theobald's 2010 formal test of universal common ancestry to 
#make sense


$parameters_of_MGWE = 0; $parameters_of_matrix = 0;
foreach $i (0..$#parameter_sum){

	$parameters_of_MGWE += (9 + (2*$parameter_sum[$i] -3));
	
}
$parameters_of_matrix = ((2*$all_seqs)-3) + (9*($#parameter_sum+1));

print StatsOut "#######################Gene Counts###############################\n";
%HASH = (); %HASH_Sites = ();
foreach $i (0..$#Tree_total){

	#print "Best Edge is: $Tree_total[$i]\n";
	if(exists $HASH{$Tree_total[$i]}){
		
		$HASH{$Tree_total[$i]}++;
		$HASH_Sites{$Tree_total[$i]} += $gene_lengths[$i];
		
	}else{
		
		$HASH{$Tree_total[$i]} = 1;
		$HASH_Sites{$Tree_total[$i]} = $gene_lengths[$i];
	}	
}
#print Dumper(\%HASH);
foreach $i (0..$#Conflict){
	
	print StatsOut "Conflict $i $Conflict[$i]: $HASH{$i}\n";
	
}

if($verbose eq "True"){
	print StatsOut "#######################Site Counts###############################\n";
	
	foreach $i (0..$#Conflict){
	
		print StatsOut "Conflict $i $Conflict[$i]: $HASH_Sites{$i}\n";
	
	}
}


print StatsOut "#######################Paramater Info############################\n";
print StatsOut "Parameters of your regular old supermatrix: $parameters_of_matrix\n";
print StatsOut "Parameters used for the Edge analysis: $parameters_of_MGWE\n";

print StatsOut "#######################Your Likelihoods##########################\n";
@temp_sort = (); @sorted = ();
if($Topos != 0){
	
	foreach $i (0..$#MatrixLikes){
		
		print StatsOut "Tree $i: $MatrixLikes[$i]\n";
		push @temp_sort, $MatrixLikes[$i];
	}
	foreach $i (0..$#gene_comp){
		print StatsOut "Indv BrLen $i: $gene_comp[$i]\n";
		push @temp_sort, $gene_comp[$i];
	}
}
foreach $i (0..$#best_likes){

	print StatsOut "Edge $i $Conflict[$i]: $like_sum[$i]\n";
	push @temp_sort, $like_sum[$i];
	
}
@sorted = sort {$b <=> $a} @temp_sort;
print StatsOut "The Best Likelihood is: $sorted[0]\n";

print StatsOut "#####################Your AIC Scores##########################\n";
@temp_sort = (); $MatrixAIC = 0; @sorted = ();
if($Topos != 0){
	
	foreach $i (0..$#MatrixLikes){
		
		$MatrixAIC = (-2*$MatrixLikes[$i]) + (2*$parameters_of_matrix);
		print StatsOut "Tree $i: $MatrixAIC\n";
		push @temp_sort, $MatrixAIC;
	}
	foreach $i (0..$#gene_comp){
		$MatrixAIC = (-2*$gene_comp[$i]) + (2*$parameters_of_MGWE);
		print StatsOut "Indv BrLen $i: $MatrixAIC\n";
		push @temp_sort, $MatrixAIC;
	}

}
foreach $i (0..$#best_likes){

	$MatrixAIC = (-2*$like_sum[$i]) + (2*$parameters_of_MGWE);
	print StatsOut "Edge $i $Conflict[$i]: $MatrixAIC\n";
	push @temp_sort, $MatrixAIC;
	
}
@sorted = sort {$a <=> $b} @temp_sort;
print StatsOut "The Best AIC is: $sorted[0]\n";
$best_aic = 0;
$best_aic = $sorted[0];
#Correcting the AICc because why not
print StatsOut "#######################Your AICc Scores########################\n";
@temp_sort = (); $MatrixAICc = 0; @sorted = ();
($start,$length_of_supermatrix) = split "-", $loc[$#loc], 2;
$correction_for_AICc_Plain = (((2*($parameters_of_matrix**2)) + (2*$parameters_of_matrix)) / ($length_of_supermatrix - $parameters_of_matrix - 1));
$correction_for_AICc_cool = (((2*($parameters_of_MGWE**2)) + (2*$parameters_of_MGWE)) / ($length_of_supermatrix - $parameters_of_MGWE - 1));
if($Topos != 0){
	
	foreach $i (0..$#MatrixLikes){
		
		$MatrixAICc = (-2*$MatrixLikes[$i]) + (2*$parameters_of_matrix) + $correction_for_AICc_Plain;
		print StatsOut "Tree $i: $MatrixAICc\n";
		push @temp_sort, $MatrixAICc;
	}
	foreach $i (0..$#gene_comp){
		$MatrixAICc = (-2*$gene_comp[$i]) + (2*$parameters_of_MGWE) + $correction_for_AICc_cool;
		print StatsOut "Indv BrLen $i: $MatrixAICc\n";
		push @temp_sort, $MatrixAICc;
	}

}
foreach $i (0..$#best_likes){

	$MatrixAICc = (-2*$like_sum[$i]) + (2*$parameters_of_MGWE) + $correction_for_AICc_cool;
	print StatsOut "Edge $i $Conflict[$i]: $MatrixAICc\n";
	push @temp_sort, $MatrixAICc;
	
}
@sorted = sort {$a <=> $b} @temp_sort;
print StatsOut "The Best AICc is: $sorted[0]\n";
$best_AICc = 0;
$best_AICc = $sorted[0];

print StatsOut "###################Your Delta AIC Scores#######################\n";
$MatrixAIC = 0; $DeltaAIC = 0; $TotalDelta = 0;
if($Topos != 0){
	
	foreach $i (0..$#MatrixLikes){
		
		$MatrixAIC = (-2*$MatrixLikes[$i]) + (2*$parameters_of_matrix);
		$DeltaAIC = $MatrixAIC - $best_aic;
		$TotalDelta += exp(-0.5 * $DeltaAIC);
		print StatsOut "Tree $i: $DeltaAIC\n";
		
	}
	foreach $i (0..$#gene_comp){
		$MatrixAIC = (-2*$gene_comp[$i]) + (2*$parameters_of_MGWE);
		$DeltaAIC = $MatrixAIC - $best_aic;
		$TotalDelta += exp(-0.5 * $DeltaAIC);
		print StatsOut "Indv BrLen $i: $DeltaAIC\n";
	}

}
foreach $i (0..$#best_likes){

	$MatrixAIC = (-2*$like_sum[$i]) + (2*$parameters_of_MGWE);
	$DeltaAIC = $MatrixAIC - $best_aic;
	$TotalDelta += exp(-0.5 * $DeltaAIC);
	print StatsOut "Edge $i $Conflict[$i]: $DeltaAIC\n";
	
}
#Correcting the AICc because why not
print StatsOut "##################Your Delta AICc Scores########################\n";
@temp_sort = (); $MatrixAICc = 0; @sorted = (); $DeltaAICc = 0;
($start,$length_of_supermatrix) = split "-", $loc[$#loc], 2;
$correction_for_AICc_Plain = (((2*($parameters_of_matrix**2)) + (2*$parameters_of_matrix)) / ($length_of_supermatrix - $parameters_of_matrix - 1));
$correction_for_AICc_cool = (((2*($parameters_of_MGWE**2)) + (2*$parameters_of_MGWE)) / ($length_of_supermatrix - $parameters_of_MGWE - 1));
if($Topos != 0){
	
	foreach $i (0..$#MatrixLikes){
		
		$MatrixAICc = (-2*$MatrixLikes[$i]) + (2*$parameters_of_matrix) + $correction_for_AICc_Plain;
		$DeltaAICc = $MatrixAICc - $best_AICc;
		print StatsOut "Tree $i: $DeltaAICc\n";
	}
	foreach $i (0..$#gene_comp){
		$MatrixAICc = (-2*$gene_comp[$i]) + (2*$parameters_of_MGWE) + $correction_for_AICc_cool;
		$DeltaAICc = $MatrixAICc - $best_AICc;
		print StatsOut "Indv BrLen $i: $DeltaAICc\n";
	}

}
foreach $i (0..$#best_likes){

	$MatrixAICc = (-2*$like_sum[$i]) + (2*$parameters_of_MGWE) + $correction_for_AICc_cool;
	$DeltaAICc = $MatrixAICc - $best_AICc;
	print StatsOut "Edge $i $Conflict[$i]: $DeltaAICc\n";
	
}


#Jury seems to still be out about the utility of weights in this context
#print StatsOut "###################Your AIC Weigths#########################\n";
#$weight = 0; $MatrixAIC = 0; $DeltaAIC = 0;
#if($Topos != 0){
	
#	foreach $i (0..$#MatrixLikes){
		
#		$MatrixAIC = (-2*$MatrixLikes[$i]) + (2*$parameters_of_matrix);
#		$DeltaAIC = $MatrixAIC - $best_aic;
#		$weight = (exp(-0.5 * $DeltaAIC) / $TotalDelta);
#		print StatsOut "Tree $i: $weight\n";
		
#	}
#	foreach $i (0..$#gene_comp){
#		$MatrixAIC = (-2*$gene_comp[$i]) + (2*$parameters_of_MGWE);
#		$DeltaAIC = $MatrixAIC - $best_aic;
#		$weight = (exp(-0.5 * $DeltaAIC) / $TotalDelta);
#		print StatsOut "Indv BrLen $i: $weight\n";
#	}
#
#}
#foreach $i (0..$#best_likes){

#	$MatrixAIC = (-2*$like_sum[$i]) + (2*$parameters_of_MGWE);
#	$DeltaAIC = $MatrixAIC - $best_aic;
#	$weight = (exp(-0.5 * $DeltaAIC) / $TotalDelta);
#	print StatsOut "Edge $i $Conflict[$i]: $weight\n";
	
#}

if($secret ne "True"){
	system("mv Unique.tre bp.log trees.unroot MatrixNoBrInfo.SSLL MatrixNoBrSSLLs.SSLL phyx.logfile $folder");
}else{
	system("mv Unique.tre bp.log trees.unroot phyx.logfile $folder");
}
system("rm RAxML_info.EX_SSLL RAxML_perSiteLLs.EX_SSLL temp.log temp.fa temptesttre TempTree.tre");


print "################################################################\n";
print "(☞ﾟヮﾟ)☞ Program is done running please double check results ☜(ﾟヮﾟ☜)\n";
print "IN THE FOLDER $folder you should have the files:\n";
print "bp.log: This contains bipartition info for you tree set\n";
print "MatrixNoBrInfo: This is the info from the Supermatrix Analysis, without Brlengths to check parameters and blah\n";
print "phyx.logfile: Your phyx logfile to make sure everything ran good and dandy\n";
print "trees.unroot: Your trees unrooted\n";
print "Unique.tre: All the unique trees from your tree set\n";
print "The results of the analysis are in the file $outfile\n";
print "################################################################\n";
