##!/usr/bin/perl


use Getopt::Long;
use Math::Complex;


my (%opts, %cov);
GetOptions(\%opts, 'h', "dir=s", "pat=s", "out=s");

my $help = "USAGE:
	perl Calculate_Coverage_Statistics.pl –dir LANDDIR –pat sufixLAND -out out.stat\n
PARAMETERS
	dir 		VL directory path
	pat 		VL files sufix
	out 		Output file\n
OUTPUT
	A file with coverage statistics:
	o	Total number of CS’s analyzed.
	o	The coverage mean
	o	The coverage standard deviation
	o	The maximum coverage
	o	The first quartile, the median, the third quartile and the IQR of the coverage\n";

die $help if $opts{h};

open(OUT, ">$opts{'out'}");

opendir(DIR, "$opts{'dir'}");
my @files = grep{/$opts{'pat'}$/} readdir(DIR);
closedir(DIR);

### Todos los datos de cobertura se guardan en el hash cov.
### $cov{30} = 1000, significa que hubo 1000 CS's que atrajeron 30 lecturas
$cont = 0;
$total = 0;
foreach $file(@files){
	$file_path = $opts{'dir'}.$file;
	#print "$file_path\n";
	open(IN, $file_path);
	while(<IN>){
		chomp($_);
		@data = split(/\t/, $_);
		if ($data[1] == 0){ next; }

		$cov{$data[1]}++;
		$cont++;
		$total = $total + $data[1];
	}
	close(IN);
	#last;
}

### int solo trunca la parte entera
print OUT "CONT - $cont\n";

my ($nmedian_all_ref, $class_all) = get_median_class_n($cont);
my @nmedian_all = @$nmedian_all_ref;

$quartile_cont = $nmedian_all[0]-1;
my ($quartile_median_ref, $class_quartile) = get_median_class_n($quartile_cont);
my @quartile_median = @$quartile_median_ref;

for ($i=0; $i<$class_quartile; $i++){
	$third_quartile_median[$i] = $quartile_median[$i] + $nmedian_all[$class_all-1];
}

print OUT "N-MEDIAN = @nmedian_all \n";
print OUT "N-QUARTILE = @quartile_median\n";
$mean = $total/$cont;
print OUT "MEAN = $mean\n";
#die;

### OBTENIENDO LA MEDIANA
### El dato de cobertura que se encuentre en la posicion media (relativa al numero de CS's)  es la mediana

$max = 0;
$sd = 0;
my (@median, @fst_quartile, @trd_quartile);

#### MEDIAN = @median
#### FIRST QUARTILE = @quartile_median
#### THIRD QUARTILE = @trd_quartile
####
#### 99%, 99.9%, 99.99%, 99.999% with coverage lower than ....
####

my ($ninety9_cont, $ninety99_cont, $ninety999_cont, $ninety9999_cont );
my ($ninety9, $ninety99, $ninety999, $ninety9999);

$ninety9_cont = int((0.99)*$cont);
$ninety99_cont = int((0.999)*$cont);
$ninety999_cont = int((0.9999)*$cont);
$ninety9999_cont = int((0.99999)*$cont);

foreach $key (sort {$a<=>$b} keys (%cov) ){  #### Ordeno los valores de cobertura
	
	$temp = $max + $cov{$key};
        $max_cov = $key;

	for ($i=0; $i<$class_all; $i++){
		if ($nmedian_all[$i] > $max && $nmedian_all[$i] <= $temp){
			push @median, $max_cov;
		}
	}
	for ($i=0; $i<$class_quartile; $i++){
		if ($quartile_median[$i] > $max && $quartile_median[$i] <= $temp){
			push @fst_quartile, $max_cov;
		}
		if ($third_quartile_median[$i] > $max && $third_quartile_median[$i] <= $temp){
			push @trd_quartile, $max_cov;
		}
		
	}
	
	if ($ninety9_cont > $max && $ninety9_cont <= $temp){ $ninety9 = $max_cov; }
        if ($ninety99_cont > $max && $ninety99_cont <= $temp){ $ninety99 = $max_cov; }
        if ($ninety999_cont > $max && $ninety999_cont <= $temp){ $ninety999 = $max_cov; }
        if ($ninety9999_cont > $max && $ninety9999_cont <= $temp){ $ninety9999 = $max_cov; }
 
	
	$max=$temp;
	$sd = (($key - $mean)**2)*$cov{$key};
	
}

$sd = $sd/($cont-1);
$sd = sqrt($sd);

if ($class_all == 1){
	$median = $median[0];
}else{
	$median = ($median[0] + $median[1])/2;
}

if ($class_quartile == 1){
	$fst_quartile = $fst_quartile[0];
	$trd_quartile = $trd_quartile[0];
}else{
	$fst_quartile = ($fst_quartile[0] + $fst_quartile[1])/2;
	$trd_quartile = ($trd_quartile[0] + $trd_quartile[1])/2;
}
print OUT "MAX_COV = $max_cov\n";
print OUT "SD = $sd\n";

$IQR  = $trd_quartile - $fst_quartile;

print OUT "FIRST QUARTILE - $fst_quartile\n";
print OUT "LA MEDIANA ES: $median\n";
print OUT "THIRD QUARTILE - $trd_quartile\n";
print OUT "EL IQR ES: $IQR\n\n";

my $frac_zero = ($cov{0}/$cont)*100;
print OUT "CS'S COV 0: $cov{0}\t$frac_zero%\n";
print OUT "99% READS COVERAGE LOWER THAN - $ninety9\n";
print OUT "99.9% READS COVERAGE LOWER THAN - $ninety99\n";
print OUT "99.99% READS COVERAGE LOWER THAN - $ninety999\n";
print OUT "99.999% READS COVERAGE LOWER THAN - $ninety9999\n";


sub get_median_class_n {
	$n = $_[0];
	if ($n%2 == 1){
		@nmedian = (int(($n/2)+0.5));
		$class = 1;
	}
	if($n%2 == 0){
		@nmedian = ($n/2, ($n/2)+1);
		$class = 2;
	}
	return(\@nmedian, $class);
}
