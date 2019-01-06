#!/usr/bin/perl

#===============================================================================
#   Author: Alejandro Manzano-Marin (https://www.researchgate.net/profile/Alejandro_Manzano-Marin/).
#
#   File: primeScaff
#   Date: 01-09-2013
#   Version: 1.6.1
#
#   Usage:
#      perl primeScaff.pl -i inFile [-d outDir] [-p outFilePrefix] [options]
#
#      Check out 'perl primeScaff.pl -h' for short usage manual and info on the software.
#
#    Description: primeScaff is intended to automate the sometines tedious process of 
#                 manually designing specific primer pairs around gaps of genomic
#                 scaffolds and speedup the genome finishing stage of a genome sequencing
#                 project. It incorporates de-novo repeat finding using RECON to avoid
#                 as much as possible designing primers in repetitive regions and offers
#                 the possibility to easily fine-tune the primer design options using primer3.
#                 It outputs the repeat, gap and primer annotations in gff2 and gff3 to easily
#                 load and visualize the results using tools like UGENE (recommended)
#                 or artemis.
#
#    Contact: Contact the author at alejandro.manzano@uv.es using 'primeScaff: ' as
#             as begining for subject for bug reporting, feature request or whatever
#             (of course realted to the software).
#
#    COPYRIGHT: Copyright (C) 2012-2013  Alejandro Manzano-Marin (https://www.researchgate.net/profile/Alejandro_Manzano-Marin/).
#
#    LICENCE: This program is free software: you can redistribute it and/or modify it under the terms
#             of the GNU General Public License as published by the Free Software Foundation, either
#             version 3 of the License, or (at your option) any later version.
#             This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#             without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#             See the GNU General Public License for more details.
#             You should have received a copy of the GNU General Public License along with this program.
#             If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


# Load modules
use Getopt::Long;
use Pod::Usage;
use File::Copy;
use File::Basename;
use Cwd;


# Define subroutines
sub printVersion {
	# Capture variables
	my $software= $_[0];
	my $version= $_[1];
	my $fileHandle= $_[2];
	print $fileHandle "$software v$version\n";
	exit (0);
}

sub fasta2seqNameLst {
	########################################
	# fasta2seqNameLst v1.0
	#
	# This subroutine runs takes the path to a FASTA fromatted file and makes the seqNameLst file needed
	# to run RECON. On success, it returns 0, on failure to open the fasta file or creating the seqNameLst
	# file it dies.
	# Usage: &fasta2seqNameLst("pathTo/file.fasta", "prefixOutFile")
	#
	########################################
	
	# Define variables
	my $seqNum= 0;      # Keeps the number of FASTA records in the FASTa file
	my @tempArray= (); # Keeps the header of the FASTA records
	my $line= "";       # Used to store the current line being read from the file
	my $seqNum= 0;      # Used to store number of seqeunces in the FSTA file
	my @tempArray= ();  # Used to store the FASTA headers
	my $tempString= ""; # Used to store the splitted line
	local (*SEQNAMELST);
	local (*FASTA);
	
	
	# Capture arguments
	my $inFile= $_[0];    # Path to the summary/eles file to axtract reads from
	my $outFile= $_[1];
	
	# Subroutine
	
	## Read eles file and extract repeat families
	open (SEQNAMELST, ">$outFile") || die ("Unable to open file for writting: $outFile\n$!\n");
	open (FASTA, "<$inFile") || die ("Unable to open file for reading: $inFile\n$!\n");
	$seqNum=0;
	while ($line=<FASTA>){
		if ($line=~ m/^>(\S+)/){
			push (@tempArray, $1);
			$seqNum++;
		}
	}
	close (FASTA);
	print SEQNAMELST "$seqNum\n";
	@tempArray= sort @tempArray; # sort lexically as required by RECON
	foreach $tempString (@tempArray){
		print SEQNAMELST $tempString . "\n";
	}
	close (SEQNAMELST);

	return (0);
}

sub runRECON {
	########################################
	# subRECON v1.0
	#
	# This subroutine runs RECON. It returns 0 on success and dies when any given system call
	# returns a non-zero value.
	# Usage: &runRECON("seqFile.fasta", "file.seqNameLst", "outfilePrefix", int[0,inf), "pathTo/formatdb", "pathTo/megablast", "pathTo/RECON-x.xx.x", int[0,1])
	#
	########################################
	
	# Define variables
	my $cmd= "";           # Used to store a command to run system calls
	
	# Capture variables
	my $inFile= $_[0];	      # Captures fasta file to use for repeat detection
	my $seqNameLstFile= $_[1]; # Captures path to file.seqNameLst used by recon.pl
	my $prefix= $_[2];	      # Captures prefix for outfiles
	my $integer= $_[3];         # Captures integer to use with recon.pl
	my $formatdbBin= $_[4];    # Captures formatdb binary path
	my $megablastBin= $_[5];   # Capturesmegablast binary path
	my $reconPath= $_[6];      # Captures RECON-x.xx.x folder path
	my $verbose= $_[7];        # Capture the verbose level [0,1]
	
	# Subroutine
	
	## Run blast for RECON
	if ($verbose){
		print STDERR "Formatting database and running megablast for RECON.. ";
	}
	
	$cmd="$formatdbBin -i $inFile -p F";
	system ($cmd) == 0 || die "Unable to run formatdb: $?\n$!\n";

	$cmd="$megablastBin -i $inFile -d $inFile -o $prefix.blast.out -p 97 -W 16 -e 1e-05 -m 8";
	system ($cmd) == 0 || die "Unable to run megablast: $?\n$!\n";
	if ($verbose){
		print STDERR "Done\n";
	}
	
	## Run MSPCollect.pl and recon.pl
	if ($verbose){
		print STDERR "Running RECON.. ";
	}
	
	$cmd="perl $reconPath/scripts/MSPCollect.pl $prefix.blast.out > $prefix.blast.MSP";
	system ($cmd) == 0 || die "Unable to run MSPCollect.pl: $?\n$!\n";
	if (-s "$prefix.blast.MSP"){
		$cmd="perl $reconPath/scripts/recon.pl $seqNameLstFile $prefix.blast.MSP $integer \'$reconPath\'";
		system ($cmd) == 0 || die "Unable to run recon.pl: $?\n$!\n";
	}
	
	if ($verbose){
		print STDERR "Done\n";
	}
	
	return (0);
}

sub eles2GFFX {
	########################################
	# eles2GFFX v1.0
	#
	# This subroutine runs takes the summary/eles file from a RECON run and writes them to a given GFF
	# file. On success, it returns the number of repeats found and number of repeat families found in
	# a string separated by a "-" character. On failure it dies.
	# Usage: &eles2GFFX("pathTo/summary/eles", *GFFFILEHANDLE, int[2,3]|"art", int[0,inf))
	#
	########################################
	
	# Define variables
	my $foundRepeats= 0;  # Used to store the repeats found
	my $foundFamilies= 0; # Used to store the number of new families found in a RECON run
	my $line= "";          # Used to store the current line being read from the file
	my @tempArray= ();    # Used to store the splitted line
	local (*ELES);
	
	
	# Capture arguments
	my $elesPath= $_[0];    # Path to the summary/eles file to axtract reads from
	local (*GFF)= $_[1];     # Filehandle for GFF file
	my $gffVersion= $_[2];  # version of GFF to write (currently 2 or 3)
	my $famNumCorr= $_[3];  # Correction number for family number. (usefull if you run iterative searches)
	
	# Subroutine
	
	## Read eles file and extract repeat families
	open (ELES, "<$elesPath") || die ("Unable to open file for reading: $elesPath\n$?\n$!\n");
	while ($line=<ELES>){
		if ($line=~ m/^#/){
			next;
		}
		$foundRepeats++;
		$line=~ s/^\s+//;
		$line=~ s/\s+$//;
		@tempArray= split(/\s+/, $line);
		if ($tempArray[0] > $foundFamilies){
			$foundFamilies= $tempArray[0]; 
		}
		$tempArray[0]+= $famNumCorr;
		if ($tempArray[2] == -1){
			$tempArray[2]="-";
		}
		elsif ($tempArray[2] == 1){
			$tempArray[2]="+";
		}
		if ($gffVersion == 2){
			print GFF "$tempArray[3]\tRECON\trepeat_region\t$tempArray[4]\t$tempArray[5]\t.\t$tempArray[2]\t.\tlocus_tag=REP_F$tempArray[0]E$tempArray[1]; rpt_family=$tempArray[0]; note=\"element:$tempArray[1]\"\n";
		}
		elsif ($gffVersion == 3){
			print GFF "$tempArray[3]\tRECON\trepeat_region\t$tempArray[4]\t$tempArray[5]\t.\t$tempArray[2]\t.\tID=REP_F$tempArray[0]E$tempArray[1];Name=REP_F$tempArray[0]E$tempArray[1];rpt_family=$tempArray[0];element=$tempArray[1]\n";
		}
		elsif ($gffVersion eq "art"){
			print GFF "$tempArray[3]\tRECON\trepeat_region\t$tempArray[4]\t$tempArray[5]\t.\t$tempArray[2]\t.\tID=REP_F$tempArray[0]E$tempArray[1];Rpt_family=$tempArray[0];Note=\"element:$tempArray[1]\"\n";
		}
	}
	close (ELES);

	return ($foundFamilies . "-" . $foundRepeats);
}

sub eles2rptHash {
	########################################
	# eles2rptHash v1.0
	#
	# This subroutine runs takes the summary/eles file from a RECON run and stores it to a given Hash.
	# On success, it returns the number of repeats found and number of repeat families found in a string
	# separated by a "-" character. On failure it dies.
	# Usage: &eles2rptHash("pathTo/summary/eles", int[0,inf), \%hash)
	#
	########################################
	
	# Define variables
	my $foundRepeats= 0;  # Used to store the repeats found
	my $foundFamilies= 0; # Used to store the number of new families found in a RECON run
	my $line= "";          # Used to store the current line being read from the file
	my @tempArray= ();    # Used to store the splitted line
	local (*ELES);
	
	
	# Capture arguments
	my $elesPath= $_[0];    # Path to the summary/eles file to axtract reads from
	my $famNumCorr= $_[1];  # Correction number for family number. (usefull if you run iterative searches)
	my $refHash= $_[2];  # Hash reference
	
	# Subroutine
	
	## Read eles file and extract repeat families
	open (ELES, "<$elesPath") || die ("Unable to open file for reading: $elesPath\n$?\n$!\n");
	while ($line=<ELES>){
		if ($line=~ m/^#/){
			next;
		}
		$foundRepeats++;
		$line=~ s/^\s+//;
		$line=~ s/\s+$//;
		@tempArray= split(/\s+/, $line);
		if ($tempArray[0] > $foundFamilies){
			$foundFamilies= $tempArray[0]; 
		}
		$tempArray[0]+= $famNumCorr;
		$refHash->{$tempArray[3]}{"REP_F$tempArray[0]E$tempArray[1]"}[0]= $tempArray[4]-1; # -1 to save in 0-based index for future use in the script
		$refHash->{$tempArray[3]}{"REP_F$tempArray[0]E$tempArray[1]"}[1]= $tempArray[5]-1; # -1 to save in 0-based index for future use in the script
	}
	close (ELES);

	return ($foundFamilies . "-" . $foundRepeats);
}

sub maskFileFASTA{
	########################################
	# maskFileFASTA v1.0
	#
	# This subroutine runs takes a hash containg repeat elements positions in the format
	# $hash{"ctgID"}{""}[0]= startIndex(zero-based);
	# $hash{"ctgID"}{""}[1]= endIndex(zero-based);
	# On success, it returns 0.
	# Usage: &maskFileFASTA("pathTo/file.fasta", "pathTo/outfile.masked.fasta", \%hash, "maskChar", int[0,1])
	#
	########################################
	
	# Define variables
	my $line= "";          # Used to store the current line being read from the file
	my $sequence= "";     # Used to store the current sequence
	my %temp= ();           # Temporary variables hash
	$temp{'string'}= "";
	$temp{'seq'}= "";
	my $rptElement= "";   # Store repeat element ID
	my $i= 0;              # Counter
	local (*NEWFASTA);
	local (*FASTA);
	
	
	# Capture arguments
	my $inFile= $_[0];    # Path to infile
	my $outFile= $_[1];   # Path to outfile
	my $refHash= $_[2];   # Hash reference
	my $maskChar= $_[3];  # Character used for masking
	my $verbose= $_[4];   # Verbose [0,1]
	
	# Subroutine
	
	## Read eles file and extract repeat families
	open (NEWFASTA, ">$outFile") || die ("Unable to open file for writing: $outFile\n$!\n");
	open (FASTA, "<$inFile") || die ("Unable to open file for reading: $inFile\n$!\n");
	while ($line=<FASTA>){
		chomp $line;
		if ($line=~ /^>(\S+)/ || eof){
			if (eof){
				$sequence.=uc($line);
			}
			if (length($sequence)){
				if ($verbose){
					print STDERR "Masking repeats for seqeunce: $seqID (length=" . length($sequence) . ")\n";
				}
				foreach $rptElement (keys %{$refHash->{$seqID}}){
					$temp{'string'}= substr($sequence, $refHash->{$seqID}{$rptElement}[0], (($refHash->{$seqID}{$rptElement}[1])-($refHash->{$seqID}{$rptElement}[0])+1));
					if ($maskChar=~ m/^lower$/i){
						$temp{'string'}= lc($temp{'string'});
					}
					else {
						$temp{'string'}=~ s/\S{1}/$maskChar/g;
					}
					if ($verbose){
						print STDERR "repeat " . $rptElement . ": " . $refHash->{$seqID}{$rptElement}[0] . ", " . $refHash->{$seqID}{$rptElement}[1] . "\n";
					}
					substr($sequence, $refHash->{$seqID}{$rptElement}[0], (($refHash->{$seqID}{$rptElement}[1])-($refHash->{$seqID}{$rptElement}[0])+1), $temp{'string'})
				}
				print NEWFASTA ">$seqID\n";
				for ($i=0; $i<length($sequence); $i+=70){
					if (length(substr($sequence, $i, length($sequence)-$i))<70){
						print NEWFASTA substr($sequence, $i, length($sequence)-$i)."\n";
					}
					else {
						print NEWFASTA substr($sequence, $i, 70)."\n";
					}
				}
			}
			$seqID=$1;
			$sequence='';
			next;
		}
		$sequence.=uc($line);
	}
	close (FASTA);
	close (NEWFASTA);
	
	return (0);
}

sub revComp {
	my $sequence= $_[0];
	$sequence= scalar(reverse($sequence));
	$sequence=~ tr/atgcATGC/tacgTACG/;
	return ($sequence);
}

# Variable definition

## Define other variables
my $line='';
my %temp= ('num' => 0);
my $newRepeats= 1;
my %repeatSeq= ();
my @tempArray= ();
my $tempElem= 0;
my $tempString= '';
my $scfId= '';
my $i= 0;
my $j= 0;
my $ctgNum= -1;
my $sequence= '';
my $gapPrefix= '';
my $gapSuffix= '';
my $maxSize= 0;
my $tempSeq= '';
my %scf= ();
my $tempNum= 0;
my $cmd= '';
my $ioVersion= '4';
my $p3File= '';
my $prevSize= 0;
my $priNumLeft= 0;
my $priNumRight= 0;
my @scfEndSeq= ();
my $fixedGap= '';
my %scfBridgeEndSeq= ();

## General variables
my $PROGRAMNAME= 'primeScaff';
my $VERSION= '1.6.1';

## Prerequisite software
my $scriptPath= dirname(__FILE__);
my $primer3Bin= $scriptPath . '/bin/primer3/primer3_core';
my $formatdbBin= $scriptPath . '/bin/blast/formatdb';
my $megablastBin= $scriptPath . '/bin/blast/megablast';
my $reconPath= $scriptPath . '/bin/recon';
my $thParamPath= $scriptPath . '/bin/primer3/primer3_config/';


## Define options default values
my $opt_inFile= '';

my $opt_outDir= 'primeScaff_out';
my $opt_prefix= 'primeScaff_out';
my $opt_best= 0;
my $opt_cleanUp= 0;

my $opt_gapChar= 'N';
my $opt_minGap= 5;
my $opt_integer= 1;

my $opt_ampSize= '100-5000';
my $opt_priSize= '18-21-27';
my $opt_priGC= '40-50-60';
my $opt_priTm= '58-60-62';
my $opt_tmFormula= 1; # Default primer3 v2.3.4: SantaLucia JR (1998) "A unified view of polymer, dumbbell and oligonucleotide DNA nearest-neighbor thermodynamics", Proc Natl Acad Sci 95:1460-65 http://dx.doi.org/10.1073/pnas.95.4.1460
my $opt_gcClamp= 0;
my $opt_polyX= 4;
my $opt_maxN= 0;
my $opt_selfAny= 47; # Primer self complementarity	
my $opt_selfEnd= 47;
my $opt_dntpConc= '0.8'; # Default primer3
my $opt_saltCorr= 1; # Default by primer3: SantaLucia JR (1998) "A unified view of polymer, dumbbell and oligonucleotide DNA nearest-neighbor thermodynamics", Proc Natl Acad Sci 95:1460-65 http://dx.doi.org/10.1073/pnas.95.4.1460].
my $opt_saltMono= '50.0';
my $opt_saltDiva= '15.0';

my $opt_verbose= 0;
my $opt_man= 0;
my $opt_help= 0;
my $opt_printVersion= 0;

## Define options hash
GetOptions(\%opts, 
	'infile|i=s' => \$opt_inFile, 
	'dir|d:s' => \$opt_outDir, 
	'prefix|p:s' => \$opt_prefix, 
	'best!' => \$opt_best, 
	'clean!' => \$opt_cleanUp, 
	'scfEnds!' => \$opt_scfEnds, 
	'gapChar:s' => \$opt_gapChar, 
	'minGap:i' => \$opt_minGap, 
	'int:i' => \$opt_integer, 
	'ampSize:s' => \$opt_ampSize, 
	'priSize:s' => \$opt_priSize, 
	'priGC:s' => \$opt_priGC, 
	'priTm:s' => \$opt_priTm, 
	'tmFormula:i' => \$opt_tmFormula, 
	'gcClamp:i' => \$opt_gcClamp, 
	'polyX:i' => \$opt_polyX, 
	'maxN:i' => \$opt_maxN, 
	'selfAny:f' => \$opt_selfAny, 
	'selfEnd:f' => \$opt_selfEnd, 
	'dntpConc:f' => \$opt_dntpConc, 
	'saltCorr:i' => \$opt_saltCorr, 
	'saltMono:f' => \$opt_saltMono, 
	'saltDiva:f' => \$opt_saltDiva, 
	'verbose|v!' => \$opt_verbose, 
	'help|h!' => \$opt_help, 
	'man!' => \$opt_man, 
	'version!' => \$opt_printVersion) || pod2usage(-exitval => 1,  -verbose => 2);

if ($opt_help){
	pod2usage(-exitval => 1,  -verbose => 1);
}
if ($opt_man){
	pod2usage(-exitval => 0, -verbose => 2);
}
if ($opt_printVersion){
	&printVersion($PROGRAMNAME, $VERSION, \*STDERR);
}


## Script documetation
=pod

=head1 NAME

primeScaff - Platform-independent automated primer design around gaps in genomic scaffolds

=head1 VERSION

primeScaff v1.6.1

=head1 SYNOPSIS

perl primeScaff.pl -i|--infile inFile.fasta [-d|--dir outDir] [-p|--prefix outFilePrefix] [--best]
[--clean] [--scfEnds] [--gapChar 'N'] [--minGap integer] [--int integer] [--ampSize 'min-max']
[--priSize 'min-opt-max'] [--priGC 'min-opt-max'] [--priTm 'min-opt-max'] [--tmFormula integer]
[--gcClamp integer] [--polyX integer] [--maxN integer] [--selfAny float] [--selfEnd float] [--dntpConc float]
[--saltCorr integer] [--saltMono float] [--saltDiva float] [-h|--help] [--man] [--version]

=head1 DESCRIPTION

primeScaff is intended to automate the sometimes tedious process of manually designing specific primer pairs
around gaps of genomic scaffolds and speedup the genome closure stage of a genome sequencing project. It
incorporates de-novo repeat finding using RECON to avoid as much as possible designing primers in repetitive
regioins and offers the possibility to easily fine-tune the primer design options using primer3. It outputs
the repeat, gap and primer annotations in gff2 and gff3 to easily load and visualize the results using tools
like UGENE (recommended) or artemis.

=head1 OPTIONS

=head2 INPUT

=over 8

=item B<-i> | B<--infile> <string> (mandatory)

File containing scaffolded contigs in FASTA format.

=back

=head2 OUTPUT

=over 8

=item B<-d> | B<--dir> <string> (default: "primeScaff_out")

Name of directory to write outfiles to.

=item B<-p> | B<--prefix> <string> (default: "primeScaff_out")

Prefix for outfiles names.

=item B<--best> <boolean> (default: 0)

Only report the best primer pair.

=item B<--clean> <boolean> (default: 0)

Clean all intermediate files.

=item B<--scfEnds> <boolean> (default: 0)

If used, also predict primers for joining all combinations of contigs.
WARNING: If you have too many scaffold (for example >=50) it could take almost forever to
finish predicting the primers, It is recommended only t use it when one has few scaffolds.

=back

=head2 TUNNING PARAMETERS

=over 8

=item B<--gapChar> <string> (default: 'N')

Character considered to be a gap.

=item B<--minGap> <integer> (default: 5)

Minimum gap Size (used to avoid considering an undefined nucleotide N as a gap).

=back

=head3 RECON

=over 8

=item B<--int> <integer> (default: 1)

Integer to use in RECON (see RECON 00README file for details).

=back

=head3 PRIMER3

=over 8

=item B<--ampSize> <string> (default: '100-5000')

Range of accepted amplicon size(s) in the format 'min-max' (see primer3 manual:
PRIMER_PRODUCT_SIZE_RANGE).

=item B<--priSize> <string> (default: '18-21-27')

Minimum, optimal and maximum primer sizes in the format "min-opt-max" (see primer3 manual:
PRIMER_MIN_SIZE, PRIMER_OPT_SIZE, PRIMER_MAX_SIZE).

=item B<--priGC> <string> (default: '40-50-60')

Minimum, optimal and maximum primer GC content in the format "min-opt-max" (see primer3 manual:
PRIMER_MIN_GC, PRIMER_OPT_GC_PERCENT, PRIMER_MAX_GC).

=item B<--priTm> <string> (default: '58-60-62')

Minimum, optimal and maximum primer Tm in the format "min-opt-max" (see primer3 manual:
PRIMER_MIN_TM, PRIMER_OPT_TM, PRIMER_MAX_TM).

=item B<--tmFormula> <integer> (default: 1; [0,1])

Specifies details of melting temperature calculation as in primer3 (see primer3 manual:
PRIMER_TM_FORMULA).

=item B<--gcClamp> <integer> (default: 0; >=0)

Specify length of primers GC clamp (see primer3 manual: PRIMER_GC_CLAMP).

=item B<--polyX> <integer> (default: 4; >=0)

Maximum number of homopolymer bases in primer (see primer3 manual: PRIMER_MAX_POLY_X).

=item B<--maxN> <integer> (default: 0; >=0)

Maximum number of "N" characters (see primer3 manual: PRIMER_MAX_NS_ACCEPTED).

=item B<--selfAny> <float> (default: 47; >=0)

Tendency of a primer to bind to itself (see primer3 manual: PRIMER_MAX_SELF_ANY_TH).

=item B<--selfEnd> <float> (default: 47; >=0)

Score of the best binding it can find of a 3'-END to an identical primer (see primer3 manual:
PRIMER_MAX_SELF_END_TH).

=item B<--dntpConc> <float> (default: 0.8; >0)

Millimolar concentration of the sum of all deoxyribonucleotide triphosphates(see primer3 manual:
PRIMER_DNTP_CONC).

=item B<--saltCorr> <integer> (default: 1; [0,2])

Specifies the salt correction formula for the melting temperature calculation (see primer3 manual:
PRIMER_SALT_CORRECTIONS).

=item B<--saltMono> <float> (default: 50; >0)

Millimolar (mM) concentration of monovalent salt cations (usually KCl) in the PCR(see primer3 manual:
PRIMER_SALT_MONOVALENT).

=item B<--saltDiva> <float> (default: 15; >=0)

Millimolar (mM) concentration of divalent salt cations (usually MgCl^(2+)) in the PCR(see primer3 manual:
PRIMER_SALT_DIVALENT).

=back

=head2 INFO AND HELP

=over 8

=item B<-v> | B<--verbose> <boolean> (default: 0)

Prints status and info messages during processing.

=item B<-h> | B<--help> <boolean>

Print useful help on using this script.

=item B<--man> <boolean>

Print the full documentation.

=item B<--version> <boolean>

Print program version.

=back

=head1 AUTHOR

Alejandro Manzano-Marin, C<< <alejandro_dot_manzano_at_uv_dot_es> >>

=head1 BUGS

If you find a bug please email me at C<< <alejandro_dot_manzano_at_uv_dot_es> >> so that I can keep making primeScaff better.

=head1 COPYRIGHT

Copyright (C) 2012-2013  Alejandro Manzano-Marin (https://www.researchgate.net/profile/Alejandro_Manzano-Marin/).

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

# Assign values to variables dependant on options

## Input

## PRIMER 3
(my $minPriSize, my $optPriSize, my $maxPriSize)= split(/-/, $opt_priSize);
(my $minPriGC, my $optPriGC, my $maxPriGC)= split(/-/, $opt_priGC);
(my $minPriTm, my $optPriTm, my $maxPriTm)= split(/-/, $opt_priTm);

# Check necessary options and paths
if ($opt_help || !$opt_inFile){
	if (!$opt_help){
		print STDERR "ERROR:\n";
		if (!$opt_inFile){
			print STDERR "infile missing\n";
		}
		print STDERR "Please check usage manual\n";
	}
	print STDERR "\n";
	pod2usage(-exitval => 1,  -verbose => 1);
}

### Check binary and library dependencies
if (!(-e $primer3Bin) || !(-e $formatdbBin) || !(-e $megablastBin) || !(-d $thParamPath) || !(-d $reconPath)){
	print STDERR "ERROR: One or more dependencies missing or not found, please check that the paths of the following are real\n";
	print STDERR "primer3         : $primer3Bin\n";
	print STDERR "primer3_config  : $thParamPath\n";
	print STDERR "recon           : $reconPath\n";
	print STDERR "formatdb        : $formatdbBin\n";
	print STDERR "megablast       : $megablastBin\n";
	pod2usage(-exitval => 1,  -verbose => 1);
}


# Main program

## Print timestamp
if ($opt_verbose){
	print STDERR "primeScaff v" . $VERSION . " start run: " . scalar(localtime) . "\n";
}

## Check if outFile directory exists, if it does delete it and re-make it
if (-d $opt_outDir) {
	unlink glob ("$opt_outDir/* $opt_outDir/*.*");
	rmdir ("$opt_outDir");
}
if (!(-d $opt_outDir)){
	mkdir ("$opt_outDir")
}

## Identify repeat sequences

### Get sequences names and sort them to create seqNameLst file for use with RECON
if ($opt_verbose){
	print STDERR "Making $opt_prefix.seqNameLst file for RECON.. ";
}

&fasta2seqNameLst($opt_inFile, "$opt_outDir/$opt_prefix.seqNameLst");

if ($opt_verbose){
	print STDERR "DONE\n";
}

### Start gff3 File
open (GFFTHR, ">$opt_outDir/$opt_prefix.seq.gff3") || die ("Unable to open file for writing: $opt_outDir/$opt_prefix.seq.gff3\n$!\n");
print GFFTHR "##gff-version 3\n";


### Run RECON iteratively until no more repeats are found
copy ("$opt_inFile", "$opt_outDir/$opt_prefix.maskedN.temp");
copy ("$opt_inFile", "$opt_outDir/$opt_prefix.masked.fasta");
$temp{'num'}= 0;

if ($opt_verbose){
	print STDERR "Running iterative repeat serch using RECON.. ";
}

while ($newRepeats>0) {
	$newRepeats= 0;
	#### Run RECON 
	&runRECON("$opt_outDir/$opt_prefix.maskedN.temp", "$opt_outDir/$opt_prefix.seqNameLst", "$opt_outDir/$opt_prefix", $opt_integer, $formatdbBin, $megablastBin, $reconPath, $opt_verbose);
	
	#### Take results from recon from file summary/eles and write them to gff files
	if (-s "$opt_outDir/$opt_prefix.blast.MSP"){
		$tempString= &eles2GFFX ("summary/eles", *GFFTHR, 3, $temp{'num'});
		&eles2rptHash ("summary/eles", $temp{'num'}, \%repeatSeq);
		$tempString=~ m/(\d+)\-(\d+)/;
		$temp{'num'}+= $1;
		$newRepeats= $2;
	
		#### Create N masked files for RECON re-run
		&maskFileFASTA("$opt_outDir/$opt_prefix.maskedN.temp", "$opt_outDir/$opt_prefix.maskedN.new", \%repeatSeq, "N", $opt_verbose);
		unlink("$opt_outDir/$opt_prefix.maskedN.temp");
		copy ("$opt_outDir/$opt_prefix.maskedN.new", "$opt_outDir/$opt_prefix.maskedN.temp");
		unlink ("$opt_outDir/$opt_prefix.maskedN.new");
	
		#### Create lower case masked files for use with primer3
		&maskFileFASTA("$opt_outDir/$opt_prefix.masked.fasta", "$opt_outDir/$opt_prefix.masked.fasta.new", \%repeatSeq, "lower", $opt_verbose);
		unlink("$opt_outDir/$opt_prefix.masked.fasta");
		copy ("$opt_outDir/$opt_prefix.masked.fasta.new", "$opt_outDir/$opt_prefix.masked.fasta");
		unlink ("$opt_outDir/$opt_prefix.masked.fasta.new");
	}
	#### Cleanup RECON files
	if ($opt_cleanUp){
		unlink glob ("edge_redef_res/* edge_redef_res/*.* ele_def_res/* ele_def_res/*.* ele_redef_res/* ele_redef_res/*.* images/* images/*.* summary/* summary/*.*");
		rmdir ("edge_redef_res");
		rmdir ("ele_def_res");
		rmdir ("ele_redef_res");
		rmdir ("images");
		rmdir ("summary");
		unlink ("formatdb.log", "$opt_outDir/$opt_prefix.maskedN.temp.nhr", "$opt_outDir/$opt_prefix.maskedN.temp.nin", "$opt_outDir/$opt_prefix.maskedN.temp.nsq", "$opt_outDir/$opt_prefix.blast.MSP", "$opt_outDir/$opt_prefix.blast.out");
	}
}
unlink("$opt_outDir/$opt_prefix.maskedN.temp");
if ($opt_cleanUp){
	unlink ("$opt_outDir/$opt_prefix.seqNameLst");
}

if ($opt_verbose){
	print STDERR "Done\n";
}


## Begin primer checking for scaffolds
if ($opt_verbose){
	print STDERR "Starting primer desing\n";
}

open (PRIFASTA, ">$opt_outDir/$opt_prefix.primer.fasta") || die ("Unable to open file for writing: $opt_outDir/$opt_prefix.primer.fasta\n$!\n");
open (FASTA, "<$opt_outDir/$opt_prefix.masked.fasta") || die ("Unable to open file for reading: $opt_outDir/$opt_prefix.masked.fasta\n$!\n");
while ($line= <FASTA>){
	chomp $line;
	if ($line=~ /^>(\S+)/ || eof){
		if (eof){
			$sequence.= $line;
		}
		if (length($sequence)){
			print PRIFASTA "# Scaffold: $scfId\n";
			$prevSize= 0;
			$opt_ampSize=~ m/\d+\-(\d+)/;
			$maxSize= $1;
			if ($opt_scfEnds){
				if (length($sequence)-$maxSize >= 0){
					push (@scfEndSeq, {'ID' => $scfId, 'fivePrime' => substr($sequence, 0, $maxSize), 'threePrime' => substr($sequence, (length($sequence)-$maxSize), $maxSize)});
				}
				else {
					push (@scfEndSeq, {'ID' => $scfId, 'fivePrime' => substr($sequence, 0, length($sequence)), 'threePrime' => substr($sequence, 0, length($sequence))});
				}
			}

			### Detect and save gaps
			$ctgNum= 0;
			while ($sequence=~ m/(N{$opt_minGap,})/gi){
				$scf{$scfId}{'gap'}{$ctgNum}[0]= $1;          # Store gap string 
				$scf{$scfId}{'gap'}{$ctgNum}[1]= length ($`); # Store gap start in zero-based index 
#				$opt_ampSize=~ m/\d+\-(\d+)/;
#				$maxSize= $1;
				if ($scf{$scfId}{'gap'}{$ctgNum}[1]-$maxSize >= 0){
					$gapPrefix= substr($sequence, $scf{$scfId}{'gap'}{$ctgNum}[1]-$maxSize, $maxSize);
				}
				else {
					$gapPrefix= substr($sequence, 0, $scf{$scfId}{'gap'}{$ctgNum}[1]);
				}
				if (length($scf{$scfId}{'gap'}{$ctgNum}[0])+$scf{$scfId}{'gap'}{$ctgNum}[1]+$maxSize+1 <= length($sequence)){
					$gapSuffix= substr($sequence, (length($scf{$scfId}{'gap'}{$ctgNum}[0])+$scf{$scfId}{'gap'}{$ctgNum}[1]), $maxSize);
				}
				else {
					$gapSuffix= substr($sequence, (length($scf{$scfId}{'gap'}{$ctgNum}[0])+$scf{$scfId}{'gap'}{$ctgNum}[1]), (length($sequence)-(length($scf{$scfId}{'gap'}{$ctgNum}[0])+$scf{$scfId}{'gap'}{$ctgNum}[1])+1));
				}
				print GFFTHR "$scfId\tprimeScaff\tgap\t" . ($scf{$scfId}{'gap'}{$ctgNum}[1]+1) . "\t" . (length($scf{$scfId}{'gap'}{$ctgNum}[0])+$scf{$scfId}{'gap'}{$ctgNum}[1]) . "\t.\t+\t.\tID=" . $scfId . "_" . $ctgNum . "-" . ($ctgNum+1) . ";Name=" . $scfId . "_" . $ctgNum . "-" . ($ctgNum+1) . ";estimated_length=" . length($scf{$scfId}{'gap'}{$ctgNum}[0]) . "\n";
				if ($opt_verbose){
					print STDERR "Predicting primers in scaffold " . $scfId . " for gap between contig " . $ctgNum . " and contig " . ($ctgNum+1) . ".. ";
				}
				$tempSeq= $gapPrefix . $scf{$scfId}{'gap'}{$ctgNum}[0] . $gapSuffix;
				$p3File= $opt_prefix . "_" . $scfId . "_" . $ctgNum . "-" . ($ctgNum+1);
				open (SET, ">$opt_outDir/$p3File.setting") || die ("Unable to open file for writing: $opt_outDir/$p3File.setting\n$!\n");
				print SET "PRIMER_FIRST_BASE_INDEX=1\n";
				print SET "SEQUENCE_ID=" . $scfId . "_" . $ctgNum . "-" . ($ctgNum+1) . "\n";
				print SET "SEQUENCE_TEMPLATE=" . $tempSeq . "\n";
				print SET "SEQUENCE_TARGET=" . (length($gapPrefix)+1) . "," . length($scf{$scfId}{'gap'}{$ctgNum}[0]) ."\n";
				print SET "PRIMER_LOWERCASE_MASKING=1\n";
				print SET "PRIMER_TASK=pick_pcr_primers\n";
				print SET "PRIMER_NUM_RETURN=5\n";
				print SET "PRIMER_PRODUCT_SIZE_RANGE=" . $opt_ampSize . "\n"; # 100-5000
				print SET "PRIMER_MIN_SIZE=" . $minPriSize . "\n"; # 18
				print SET "PRIMER_OPT_SIZE=" . $optPriSize . "\n"; # 21
				print SET "PRIMER_MAX_SIZE=" . $maxPriSize . "\n"; # 27
				print SET "PRIMER_MIN_GC=" . $minPriGC . "\n"; # 40
				print SET "PRIMER_OPT_GC_PERCENT=" . $optPriGC . "\n"; # 50
				print SET "PRIMER_MAX_GC=" . $maxPriGC . "\n"; # 60
				print SET "PRIMER_MIN_TM=" . $minPriTm . "\n"; # 58
				print SET "PRIMER_OPT_TM=" . $optPriTm . "\n"; # 60
				print SET "PRIMER_MAX_TM=" . $maxPriTm . "\n"; # 62
				print SET "PRIMER_TM_FORMULA=" . $opt_tmFormula . "\n";
				print SET "PRIMER_GC_CLAMP=" . $opt_gcClamp ."\n"; # usually 1 or 2
				print SET "PRIMER_MAX_NS_ACCEPTED=0\n"; # Make sure you do not allow N characters in primers
				print SET "PRIMER_MAX_POLY_X=" . $opt_polyX . "\n"; # 5
				print SET "PRIMER_MAX_NS_ACCEPTED=" . $opt_maxN . "\n"; # 0
				print SET "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1\n";
				print SET "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=1\n";
				print SET "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" . $thParamPath . "\n";
				print SET "PRIMER_MAX_SELF_ANY_TH=" . $opt_selfAny . "\n";# default 47
				print SET "PRIMER_MAX_SELF_END_TH=" . $opt_selfEnd . "\n"; # default 47
				print SET "PRIMER_DNTP_CONC=" . $opt_dntpConc . "\n"; # 0.8
				print SET "PRIMER_SALT_CORRECTIONS=" . $opt_saltCorr . "\n";
				print SET "PRIMER_SALT_MONOVALENT=" . $opt_saltMono . "\n"; #KCl or NaCl 50
				print SET "PRIMER_SALT_DIVALENT=" . $opt_saltDiva . "\n"; #MgCl 15
				print SET "P3_FILE_FLAG=1\n";	
				print SET "=";
				close (SET);
				$cmd= "$primer3Bin -format_output -io_version=". $ioVersion . " -echo_settings_file -strict_tags -output=" . $opt_outDir . "/" . $p3File . ".out < " . $opt_outDir . "/" . $p3File . ".setting";
				system ($cmd) == 0 || die "Unable to run priner3_core: error code $?\n$!\n";
				open (PRIMER, "<$opt_outDir/$p3File.out") || die ("Unable to open file for writing: $opt_outDir/$p3File.out\n$!\n");
				$priNumLeft=0; # Start left primer count for gap
				$priNumRight=0; # Start right primer count for gap
				while ($line=<PRIMER>){
					if ($opt_best && $priNumLeft==1 && $priNumRight==1){
						last;
					}
					# Fields correspond to: primer_side[0] start[1] length[2] tm[3] gc%[4] any_th[5] 3'_th[6] hairpin[7] primer_seq[8]
					if ($line=~ m/(LEFT|RIGHT) PRIMER\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+([ATGCatgc]+)$/){
						@tempArray=($1, $2, $3, $4, $5, $6, $7, $8, $9);
						print PRIFASTA ">". $scfId . "_" . $ctgNum . "-" . ($ctgNum+1) . "_";
						if ($tempArray[0] eq "LEFT"){
							print PRIFASTA $priNumLeft . "_F\n$tempArray[8]\n";
						}
						else {
							print PRIFASTA $priNumRight . "_R\n$tempArray[8]\n\n";
						}
						if (!$tempHash{$tempArray[0]}{$tempArray[1] . "-" . $tempArray[2]}){ # Check if primer is not defined yet since many primer pairs can share one primer
							$tempHash{$tempArray[0]}{$tempArray[1] . "-" . $tempArray[2]}=1;
							$i++;
							print GFFTHR "$scfId\tprimeScaff\tprimer\t";
							if ($tempArray[0] eq "LEFT"){
								print GFFTHR ($scf{$scfId}{'gap'}{$ctgNum}[1]-length($gapPrefix)+$tempArray[1]) . "\t" . ($scf{$scfId}{'gap'}{$ctgNum}[1]-length($gapPrefix)+$tempArray[1]+$tempArray[2]-1) . "\t.\t+\t.\tID=" . $scfId . "_" . $ctgNum . "-" . ($ctgNum+1) . "_" . $priNumLeft . "_F" . ";Name=" . $scfId . "_" . $ctgNum . "-" . ($ctgNum+1) . "_" . $priNumLeft . "_F"; # If left primer, no need to calculate start position
								$priNumLeft++;
							}
							else {
								print GFFTHR ($scf{$scfId}{'gap'}{$ctgNum}[1]-length($gapPrefix)+$tempArray[1]-$tempArray[2]+1) . "\t". ($scf{$scfId}{'gap'}{$ctgNum}[1]-length($gapPrefix)+$tempArray[1]) . "\t.\t-\t.\tID=" . $scfId . "_" . $ctgNum . "-" . ($ctgNum+1) . "_" . $priNumRight . "_R" . ";Name=" . $scfId . "_" . $ctgNum . "-" . ($ctgNum+1) . "_" . $priNumRight . "_R";
								$priNumRight++;
							}
							print GFFTHR ";Tm=$tempArray[3];GC=$tempArray[4];any_th=$tempArray[5];3_th=$tempArray[6];hairpin=$tempArray[7];sequence=$tempArray[8]\n";
						}
						@tempArray=();
					}
				}
				%tempHash=();
				close (PRIMER);
				if ($opt_verbose){
					print STDERR "Done\n";
				}
				$ctgNum++;
			}
		} # Close if (length($sequence)) Line 306
		$scfId=$1;
		$ctgNum=0;
		$sequence='';
		### Cleanup primer3 files
		if ($opt_cleanUp){
			unlink glob ("*.for *.rev $opt_outDir/*.out $opt_outDir/*.setting");
		}
		next;
	} # Close if ($line=~ /^>(\S+)/ || eof) Line 302
	$sequence.= $line;
} # Close while ($line=<FASTA>) Line 300
close (FASTA);
close (PRIFASTA);

if ($opt_scfEnds && scalar(@scfEndSeq)>=2){
	for ($i=0; $i<1000; $i++){
		$fixedGap.= 'N';
	}
	open (SCFPRIFASTA, ">$opt_outDir/$opt_prefix.scfPrimer.fasta") || die ("Unable to open file for writing: $opt_outDir/$opt_prefix.scfPrimer.fasta\n$!\n");
	open (SCFPRITXT, ">$opt_outDir/$opt_prefix.scfPrimer.txt") || die ("Unable to open file for writing: $opt_outDir/$opt_prefix.scfPrimer.txt\n$!\n");
	for ($i=0; $i< scalar(@scfEndSeq); $i++){
		for ($j=$i+1; $j< scalar(@scfEndSeq); $j++){
			if ($scfEndSeq[$i]{'ID'} eq $scfEndSeq[$j]{'ID'}){
				next;
			}

			$scfBridgeEndSeq{($scfEndSeq[$i]{'ID'} . '-five_' . $scfEndSeq[$j]{'ID'} . '-five')}= $scfEndSeq[$i]{'fivePrime'} . '[' .  $fixedGap . ']' . $scfEndSeq[$j]{'fivePrime'};
			$scfBridgeEndSeq{($scfEndSeq[$i]{'ID'} . '-five_' . $scfEndSeq[$j]{'ID'} . '-three')}= $scfEndSeq[$i]{'fivePrime'} . '[' .  $fixedGap . ']' . &revComp($scfEndSeq[$j]{'threePrime'});
		}
	}
	foreach $scfId (sort {lc($a) cmp lc($b)} keys %scfBridgeEndSeq){
		
		if ($opt_verbose){
			print STDERR "Predicting primers between scaffolds: " . $scfId . ".. ";
		}
		
		$tempSeq= $scfBridgeEndSeq{$scfId};
		$tempSeq=~ m/^(\w+)\[N{1000}\]/;
		$tempSeq=~ s/\[//;
		$tempSeq=~ s/\]//;
		$gapPrefix=length($1);
		$p3File= $scfId;
		open (SET, ">$opt_outDir/$p3File.setting") || die ("Unable to open file for writing: $opt_outDir/$p3File.setting\n$!\n");
		print SET "PRIMER_FIRST_BASE_INDEX=1\n";
		print SET "SEQUENCE_ID=" . $scfId . "\n";
		print SET "SEQUENCE_TEMPLATE=" . $tempSeq . "\n";
		print SET "SEQUENCE_TARGET=" . ($gapPrefix+1) . ",1000\n";
		print SET "PRIMER_LOWERCASE_MASKING=1\n";
		print SET "PRIMER_TASK=pick_pcr_primers\n";
		print SET "PRIMER_NUM_RETURN=5\n";
		print SET "PRIMER_PRODUCT_SIZE_RANGE=" . $opt_ampSize . "\n"; # 100-5000
		print SET "PRIMER_MIN_SIZE=" . $minPriSize . "\n"; # 18
		print SET "PRIMER_OPT_SIZE=" . $optPriSize . "\n"; # 21
		print SET "PRIMER_MAX_SIZE=" . $maxPriSize . "\n"; # 27
		print SET "PRIMER_MIN_GC=" . $minPriGC . "\n"; # 40
		print SET "PRIMER_OPT_GC_PERCENT=" . $optPriGC . "\n"; # 50
		print SET "PRIMER_MAX_GC=" . $maxPriGC . "\n"; # 60
		print SET "PRIMER_MIN_TM=" . $minPriTm . "\n"; # 58
		print SET "PRIMER_OPT_TM=" . $optPriTm . "\n"; # 60
		print SET "PRIMER_MAX_TM=" . $maxPriTm . "\n"; # 62
		print SET "PRIMER_TM_FORMULA=" . $opt_tmFormula . "\n";
		print SET "PRIMER_GC_CLAMP=" . $opt_gcClamp ."\n"; # usually 1 or 2
		print SET "PRIMER_MAX_NS_ACCEPTED=0\n"; # Make sure you do not allow N characters in primers
		print SET "PRIMER_MAX_POLY_X=" . $opt_polyX . "\n"; # 5
		print SET "PRIMER_MAX_NS_ACCEPTED=" . $opt_maxN . "\n"; # 0
		print SET "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1\n";
		print SET "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=1\n";
		print SET "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" . $thParamPath . "\n";
		print SET "PRIMER_MAX_SELF_ANY_TH=" . $opt_selfAny . "\n";# default 47
		print SET "PRIMER_MAX_SELF_END_TH=" . $opt_selfEnd . "\n"; # default 47
		print SET "PRIMER_DNTP_CONC=" . $opt_dntpConc . "\n"; # 0.8
		print SET "PRIMER_SALT_CORRECTIONS=" . $opt_saltCorr . "\n";
		print SET "PRIMER_SALT_MONOVALENT=" . $opt_saltMono . "\n"; #KCl or NaCl 50
		print SET "PRIMER_SALT_DIVALENT=" . $opt_saltDiva . "\n"; #MgCl 15
		print SET "P3_FILE_FLAG=1\n";	
		print SET "=";
		close (SET);
		$cmd= "$primer3Bin -format_output -io_version=". $ioVersion . " -echo_settings_file -strict_tags -output=" . $opt_outDir . "/" . $p3File . ".out < " . $opt_outDir . "/" . $p3File . ".setting";
		system ($cmd) == 0 || die "Unable to run priner3_core: error code $?\n$!\n";
		open (PRIMER, "<$opt_outDir/$p3File.out") || die ("Unable to open file for writing: $opt_outDir/$p3File.out\n$!\n");
		$priNumLeft=0; # Start left primer count for gap
		$priNumRight=0; # Start right primer count for gap
		while ($line=<PRIMER>){
#			if ($opt_best && $priNumLeft==1 && $priNumRight==1){
			if ($priNumLeft==1 && $priNumRight==1){
				last;
			}
			# Fields correspond to: primer_side[0] start[1] length[2] tm[3] gc%[4] any_th[5] 3'_th[6] hairpin[7] primer_seq[8]
			if ($line=~ m/^(LEFT|RIGHT) PRIMER\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+([ATGCatgc]+)$/){
				#[1]LEFT/RIGHT [2]start [3]len [4]tm [5]gc% [6]any_th [7]3'_th [8]hairpin [9]seq
				@tempArray=($1, $2, $3, $4, $5, $6, $7, $8, $9);
				print SCFPRIFASTA ">". $scfId;
				if ($tempArray[0] eq "LEFT"){
					print SCFPRIFASTA $priNumLeft . "_F Tm=$tempArray[3] GC=$tempArray[4] any_th=$tempArray[5] 3_th=$tempArray[6] hairpin=$tempArray[7]\n$tempArray[8]\n";
				}
				else {
					print SCFPRIFASTA $priNumRight . "_R Tm=$tempArray[3] GC=$tempArray[4] any_th=$tempArray[5] 3_th=$tempArray[6] hairpin=$tempArray[7]\n$tempArray[8]\n\n";
				}
				if (!$tempHash{$tempArray[0]}{$tempArray[1] . "-" . $tempArray[2]}){ # Check if primer is not defined yet since many primer pairs can share one primer
					$tempHash{$tempArray[8]}{$tempArray[1] . "-" . $tempArray[2]}=1;
				}
#				if (!$tempHash{$tempArray[0]}{$tempArray[1] . "-" . $tempArray[2]}){ # Check if primer is not defined yet since many primer pairs can share one primer
#					$tempHash{$tempArray[0]}{$tempArray[1] . "-" . $tempArray[2]}=1;
#					$i++;
#					print GFFTHR "$scfId\tprimeScaff\tprimer\t";
#					if ($tempArray[0] eq "LEFT"){
#						print GFFTHR ($scf{$scfId}{'gap'}{$ctgNum}[1]-length($gapPrefix)+$tempArray[1]) . "\t" . ($scf{$scfId}{'gap'}{$ctgNum}[1]-length($gapPrefix)+$tempArray[1]+$tempArray[2]-1) . "\t.\t+\t.\tID=" . $scfId . "_" . $ctgNum . "-" . ($ctgNum+1) . "_" . $priNumLeft . "_F" . ";Name=" . $scfId . "_" . $ctgNum . "-" . ($ctgNum+1) . "_" . $priNumLeft . "_F"; # If left primer, no need to calculate start position
#						$priNumLeft++;
#					}
#					else {
#						print GFFTHR ($scf{$scfId}{'gap'}{$ctgNum}[1]-length($gapPrefix)+$tempArray[1]-$tempArray[2]+1) . "\t". ($scf{$scfId}{'gap'}{$ctgNum}[1]-length($gapPrefix)+$tempArray[1]) . "\t.\t-\t.\tID=" . $scfId . "_" . $ctgNum . "-" . ($ctgNum+1) . "_" . $priNumRight . "_R" . ";Name=" . $scfId . "_" . $ctgNum . "-" . ($ctgNum+1) . "_" . $priNumRight . "_R";
#						$priNumRight++;
#					}
#					print GFFTHR ";Tm=$tempArray[3];GC=$tempArray[4];any_th=$tempArray[5];3_th=$tempArray[6];hairpin=$tempArray[7];sequence=$tempArray[8]\n";
#				}
				@tempArray=();
			}
		}
		%tempHash=();
		close (PRIMER);
		if ($opt_verbose){
			print STDERR "Done\n";
		}
	}
	close (SCFPRIFASTA);
	close (SCFPRITXT);
}


## Create GFF2 and GFF3 artemis compatible and add sequences to GFF3 file
if ($opt_verbose){
	print STDERR "Making gff files with embeded sequences.. ";
}

### Create GFF2 and add headers
open (GFFTWO, ">$opt_outDir/$opt_prefix.seq.gff2") || die ("Unable to open file for writing: $opt_outDir/$opt_prefix.seq.gff2\n$!\n");
print GFFTWO "##gff-version 2\n##source-version $PROGRAMNAME v$VERSION\n##source-version RECON v1.07.1\n";

### Add Add sequences to GFF2 and GFF3
print GFFTHR "##FASTA\n";
open (NEWFASTA, "<$opt_outDir/$opt_prefix.masked.fasta") || die ("Unable to open file for reading: $opt_outDir/$opt_prefix.masked.fasta\n$!\n");
$i=0;
while ($line=<NEWFASTA>){
	print GFFTHR $line;
	if ($line=~ m/^>/){
		$line=~ s/^>/DNA /;
		if ($i){
			print GFFTWO "##end-DNA\n";
		}
		$i++;
	}
	$line= "##" . $line;
	print GFFTWO $line;
}
close (NEWFASTA);
#print GFFTHR "\n###\n";
print GFFTWO "##end-DNA\n";
close (GFFTHR);

### Add features to GFF2 file and write GFF3 artemis compatible file
open (GFFART, ">$opt_outDir/$opt_prefix.art.seq.gff3") || die ("Unable to open file for writing: $opt_outDir/$opt_prefix.art.seq.gff3\n$!\n");
open (GFFTHR, "<$opt_outDir/$opt_prefix.seq.gff3") || die ("Unable to open file for writing: $opt_outDir/$opt_prefix.seq.gff3\n$!\n");
while ($line=<GFFTHR>){
	chomp $line;
	if ($line=~ m/^[^\t]+\tRECON\trepeat_region/){
		$line=~ s/Name=REP_F\d+E\d+;//;
		$line=~ s/element=/Note=\"element:/;
		$line.= "\"";
		print GFFART $line . "\n";
		$line=~ s/ID=REP/locus_tag=REP/;
		$line=~ s/Note=/note=/;
		$line=~ s/;/; /g;
		print GFFTWO $line . "\n";
	}
	elsif ($line=~ m/^[^\t]+\tprimeScaff\tgap/){
		$line=~ s/Name=[^;]+;//;
		print GFFART $line . "\n";
		$line=~ s/ID=[^;]+;//;
		print GFFTWO $line . "\n";
	}
	elsif ($line=~ m/^[^\t]+\tprimeScaff\tprimer/){
		$line=~ s/\tprimer\t/\tprimer_bind\t/;
		$line=~ s/Name=[^;]+;/Note=\"/;
		$line=~ s/Tm=/Tm:/;
		$line=~ s/GC=/GC:/;
		$line=~ s/any_th=/any_th:/;
		$line=~ s/3_th=/3_th:/;
		$line=~ s/hairpin=/hairpin:/;
		$line=~ s/sequence=/sequence:/;
		$line.= "\"";
		print GFFART $line . "\n";
		$line=~ s/ID=/locus_tag=/;
		$line=~ s/Note=/note=/;
		$line=~ s/;/; /;
		print GFFTWO $line . "\n";
	}
	else{
		print GFFART $line . "\n";
	}
}
close (GFFTHR);
close (GFFART);
close (GFFTWO);					

if ($opt_verbose){
	print STDERR "Done\n";
	print STDERR "primeScaff v" . $VERSION . " end run: " . scalar(localtime) . "\n";
}


exit (0);
