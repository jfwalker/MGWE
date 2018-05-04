use Data::Dumper;

#This subroutine gets the different gene lengths
#This passes all tests
sub GetParts {
	
	@temp_array = ();
	@temp_array = @_;
	$PartFile = $temp_array[0];
	$IsItATest = "";
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

##############
#Things to mess with
$raxml = "raxmlHPC";
$threads = 3;
##############



$SuperMatrix = $ARGV[0] || die "perl TwoTop.pl Supermatrix Model Trees outfile\n";
$PartFile = $ARGV[1] || die "perl TwoTop.pl Supermatrix Model Trees outfile\n";
$Trees = $ARGV[2] || die "perl TwoTop.pl Supermatrix Model Trees outfile\n";
$Outfile = $ARGV[3] || die "perl TwoTop.pl Supermatrix Model Trees outfile\n";
open(out, ">$Outfile");
open(ForR, ">ForR");

@loc = ();
@loc = GetParts($PartFile);
system("$raxml -f g -T $threads -s $SuperMatrix -q $PartFile -m GTRGAMMA -z $Trees -n Topologies_SSLL");
$count = 0; @array = (); @array2 = (); @sum_array = ();
@array1 = (); $site1sum = 0; $site2sum = 0;
@site1_array = (); @site2_array = ();
open(file, "RAxML_perSiteLLs.Topologies_SSLL");
while($line = <file>){

	chomp $line;
	if($count != 0){
		
		($one,$two) = split " ", $line, 2;
		@array = split " ", $two;
		if($count == 1){
			
			@site1_array = split " ", $two;
		}else{
			
			@site2_array = split " ", $two;
		}
		$supersum = 0;
		foreach $i (0..$#loc){
			
			($one, $two) = split "-", $loc[$i];
			$one -= 1;
			$two -= 1;
			$sum = 0;
			
			foreach $j ($one..$two){
				
				$sum += $array[$j];
				$supersum += $array[$j];
			}
			if($count == 1){
				push @array1, $sum;
			}else{
				push @array2, $sum;
			}
		}
		push @sum_array, $supersum;
	}
	$count++;
}
@part_array = ();
open(Parts, "$PartFile")||die "Can't Find Parts file";
while($line = <Parts>){

	$line =~ m/(.*?) = .*/;
	push @part_array, $1;


}
$diff = abs($sum_array[0] - $sum_array[1]);
print out "######################################################\n";
print out "Difference in supermatrices: $diff\n";
print out "######################################################\n";
$SuperSum = 0; $CoalSum = 0;
foreach $i (0..$#array1){

	$diff = ($array1[$i] - $array2[$i]);
	print out "$part_array[$i]: $diff\n";
	print ForR "$diff,";
	if($array2[$i] < $array1[$i]){
		
		$SuperSum++;
	}else{
		
		$CoalSum++;	
	}
	
	
}
print ForR "\n";
print out "#######################Gene Counts#######################\n";
print out "Supermatrix: $SuperSum\nCoalescent: $CoalSum\n";
print out "#######################Site Counts#######################\n";
$super_site = 0; $coal_site = 0;
foreach $i (0..$#site1_array){
	
	if($site2_array[$i] < $site1_array[$i]){
		$super_site++;
	}else{
		$coal_site++;
	}
	
}
print out "Supermatrix: $super_site\nCoalescent: $coal_site\n";
