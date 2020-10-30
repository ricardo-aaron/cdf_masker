# cdf_masker.pl

cdf_masker.pl Perl script:

```perl
#!/usr/bin/perl

# ARGV:
# 0:cdf_file, 1:CELfile, 2:newchipname, 3:probeset list

use lib "/path/to/Bio-Affymetrix-0.5/lib/";
use Bio::Affymetrix::CDF;
use Getopt::Long qw/:config gnu_getopt/;

# use or modification of this software should reference: Hammond et al. Plant Methods 2005, 1:10. Using genomic DNA-based probe-selection to improve the sensitivity of high-density oligonucleotide arrays when applied to heterologous species.

# version 2.0 - Friday 13th April 2007


sub usage {
    my $fh=shift;

    print $fh "Usage: ".$0." {options} CDF_file_name CEL_file_name new_chip_type

options are:
-t, --threshold: CEL file intensity cut-off threshold (default:700)
-h, --thresholdH: CEL file HIGH intensity cut-off threshold (default:1000000)
-p, --percentage: FOR SIMULATIONS - random percent of probes to be removed (over-rides all other options)
-o, --output: File to output to (default standard out)
-v, --verbose: Be verbose.
--help, --usage: show this message
";

}

my $threshold=700;
my $thresholdH=1000000;
my $thresholdPer=0;
my $output;
my $usage;
my $verbose;

if (!GetOptions(
				"threshold|t=i"=>\$threshold,
				"output|o=s"=>\$output,
				"verbose|v"=>\$verbose,
				"help|usage"=>\$usage,
				"thresholdH|h=i"=>\$thresholdH,
				"thresholdPer|p=i"=>\$thresholdPer,)
	)
	{
	    usage(\*STDERR);
	    exit 0;
	};

if ($usage)
{
    usage(\*STDOUT);
    exit 1;
}

if (scalar @ARGV!=4)
{
    usage(\*STDERR);
    exit 0;
}

my $CDFfile = $ARGV[0];
my $CELfile = $ARGV[1];
my $newchipname = $ARGV[2];
my $file_affys = $ARGV[3];
# Build up coordinates of CEL file with intensity > threshold

if ($verbose)
{
    print STDERR "Reading CEL file and finding probes above threshold...";
}

open CEL, "<",$CELfile or die "Can't open CEL file $CELfile";
binmode CEL,":crlf";
# Get to the intensity section
{
    my $junk;
    while (defined($junk=<CEL>) && $junk !~/\[INTENSITY\]/) {;}
	# Two lines of trivia to skip
    $junk=<CEL>;
    $junk=<CEL>;
}

my $count=0;
my $throw=0;
my $keep=0;

# Make a grid of cells that are over the threshold
# The grid -to be used as a 2D array
my @grid;
while (my $line=<CEL>)
{
    chomp $line;
    # Stop at the [MASKS] section
    if ($line eq "[MASKS]")
    {
		last;
    }
    if ($line ne "")
    {
		$line=~/^\s*([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\.\d]+)\s*$/i or die "Can't understand CEL file format";
		if ($thresholdPer)
		{
			if (int(rand(100))>$thresholdPer)
			{
			    $grid[$1][$2]=1;
				$keep++;
		    }
		    else
		    {
				$throw++;
			}
		}
		else
		{
			if (($3 > $threshold)&&($3 < $thresholdH))
			{
			    $grid[$1][$2] = 1;
		    }
		}
		$count++;
		if ($verbose && ($count % 10000 == 0))
		{
		    print STDERR ".";
		}
    }
}
close CEL;
my $percent = int(100*($throw/$count));

if ($verbose)
{
	print STDERR "\n\nParsing CDF file...";
}

# So now we have @grid- a two dimensional array where $grid[$x][$y]
# returns 1 if the cel file has it over the threshold

#Now read the FILTER file, which contains all the affy ids to preserve in the CDF file

#### IMPORTANTE

my %keep_ids = ();
open(FIL, $file_affys) || die;   #Archivo con los affy_ids que deseamos analizar
while(<FIL>)
{
	next if ($_ =~ /^#/);
	chomp($_);
	$keep_ids{$_} = 1;
}
close(FIL);

# Parse existing CDF file
my $CDF=new Bio::Affymetrix::CDF({probemode=>1});
$CDF->parse_from_file($CDFfile);
if ($verbose)
{
	print STDERR "done\nRemoving probes...";
}

my $count2 = 0;

my $num_probesets_removed = 0;
my $num_probes_removed = 0;
my $num_probesets_kept = 0;
my $num_probes_kept = 0;

# Plan:
# Loop over every probeset on the array. For each probeset, loop over
# every probe_pair. For every probe, find the perfect match, and check
# if its X/Y coordinates is in the list above. If not, then delete
# it. Finally check that the probeset as a whole has some probes left-
# if not, delete it.

# Note- the Affymetrix software has the capability for probe_pairs
# with three mismatches and one perfect match. This software can
# handle this, although as far as we are aware, no chips of this type
# were ever produced.

my $probes=$CDF->probesets();

# Loop over every probeset in the CDF file
foreach my $probename (keys %$probes)
{
    # Need to filter out probes with unit numbers <1000- these are controls
    if ($probename >= 0)
    {
		my $nombre = $probes->{$probename}->name();
		if($keep_ids{$nombre})
		{
			# Loop over every probepair in the probeset
			my $probepairs = $probes->{$probename}->probe_pairs();
			# A new list of sucessful probe pairs
			my $newprobepairs = [];
			# Add every successful probeset to the new list
			foreach my $z (@$probepairs)
			{
				$num_probes_kept++;
				# Add to list
				push @$newprobepairs,$z;
			}
			# Replace old list with new list
			$probes->{$probename}->probe_pairs($newprobepairs);
			# If a probeset has no probes left, we have to remove the entire probeset
			$num_probesets_kept++;
        }
        else
        {
			$num_probesets_removed++;
			delete $probes->{$probename};
        }
	}
    $count2++;
    if ($verbose &&($count2 % 1000 ==0))
    {
		print STDERR ".";
	}
}

if ($verbose)
{
	print STDERR "done\nWriting CDF file...";
}

$CDF->name($newchipname);

if (!defined $output)
{
    print $CDF->write_to_filehandle(\*STDOUT,"MAS5");
}
else
{
    $CDF->write_to_file($output,"MAS5");
}

if ($verbose)
{
	my $metrics = "\nMetrics\n=======\n\n\n\ttail=$threshold\n\ttop=$thresholdH\n\trandom=$thresholdPer\%\n\n\tTotal original probes (PM + MM): $count\n\t*Random* Throw: $throw\n\t*Random* Keep $keep\n\tActual percentage of all *random* probes removed: $percent\n\n\tNumber of probe PAIRS removed: $num_probes_removed\n\tNumber of probe SETS removed: $num_probesets_removed.\n\n\tNumber of probe PAIRS kept: $num_probes_kept\n\tNumber of probe SETS kept: $num_probesets_kept.\n";
	print STDERR "done\n$metrics";

	my $metricname = $newchipname.".metrics.txt";
	open(MET, ">>$metricname") || die "Can't open output file";
	print MET $metrics;
	close(MET);
	print "\n..generating/appending-to metrics file ($metricname) done...\n";
}

1;
```
