#!/usr/bin/perl

#===============================================================================
#   Author: Alejandro Manzano-Marin (https://www.researchgate.net/profile/Alejandro_Manzano-Marin/).
#
#   File: primerCheck
#   Date: 08-04-2013
#   Version: 1.1.1
#
#   Usage:
#      perl primerCheck.pl -i inFile -s seqFile.fasta -p outFilePrefix [options]
#
#      Check out 'perl primerCheck.pl -h' for short usage manual and info on the software.
#
#    Description: This script extends primeScaff by using primersearch from the emboss package
#                 to check for the amplicons that can be produced using your primer pairs. It
#                 outputs a simple tabular file in the format:
#                 primerPairID\tprimerFseq\tprimerRseq\tsequenceHit\tstart\tsize
#
#    Contact: Contact the author at alejandro.manzano@uv.es using 'primeScaff: ' as
#             as begining for subject for bug reporting, feature request or whatever
#             (of course realted to the software).
#
#    COPYRIGHT: Copyright (C) 2012-2013 Alejandro Manzano-Marin (https://www.researchgate.net/profile/Alejandro_Manzano-Marin/).
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
use File::Basename;
use File::Copy;


# Define subroutines
sub printVersion {
	my $software= $_[0];
	my $version= $_[1];
	my $fileHandle= $_[2];
	print $fileHandle "$software v$version\n";
	exit (0);
}


# Variable definition

## Define other variables
my $line= '';
my $tempString= '';
my $tempString2= '';
my $cmd= '';
my %primers= ();
my %primerBind= ();
my $primerId= '';
my $ampNum= 0;

## General variables
my $PROGRAMNAME= 'primerCheck';
my $VERSION= '1.1.1';

## Prerequisite software (Defaults were done according to a local install)
my $scriptPath= dirname(__FILE__);
my $primersearchBin= $scriptPath . '/bin/EMBOSS/emboss/primersearch';

## Define options default values
my $opt_inFile= '';
my $opt_seqFile= '';
my $opt_mismatch= 5; #If primers 20 or over this is equivalent to around 1 bp

my $opt_prefix= 'primerCheck_out';

my $opt_verbose= 0;
my $opt_man= 0;
my $opt_help= 0;
my $opt_printVersion= 0;

## Define options hash
GetOptions(\%opts, 
	'infile|i=s', => \$opt_inFile, 
	'seq|s=s', => \$opt_seqFile, 
	'prefix|p:s', => \$opt_prefix, 
	'mismatch:i', => \$opt_mismatch, 
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

primerCheck - Checking all primer amplicons

=head1 VERSION

primerCheck v1.1.1

=head1 SYNOPSIS

primerCheck.pl -i|--infile inFile --seq seqFile.fasta --prefix outFilePrefix [--mismatch integer]
[-h|--help] [--man] [--version]

=head1 DESCRIPTION

primerCheck script extends primeScaff by using primersearch from the emboss package to check
for the amplicons that can be produced using your primer pairs. It outputs a simple tabular file

=head1 OPTIONS

=head2 INPUT

=over 8

=item B<-i> | B<--infile> <string> (mandatory)

File containing primer FASTA file (provided by primeScaff).

=item B<-s> | B<--seq> <string> (mandatory)

File containing sequences in FASTA format.

=back

=head2 OUTPUT

=over 8

=item B<-p> | B<--prefix> <string> (default: "priScaff_out")

Prefix for outfile name.

=back

=head2 TUNNING PARAMETERS

=over 8

=item B<--mismatch> <integer> (default: 5)

Mismatch percent allowed while running primersearch.

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


# Check necessary options and paths
if ($opt_help || !$opt_inFile || !$opt_seqFile){
	if (!$opt_help){
		print STDERR "ERROR:\n";
		if (!$opt_inFile){
			print STDERR "infile missing\n";
		}
		if (!$opt_seqFile){
			print STDERR "Sequence FASTA file path missing\n";
		}
		print STDERR "Please check usage manual\n";
	}
	print STDERR "\n";
	pod2usage(-exitval => 1,  -verbose => 1);
}


### Check binary and library dependencies
if (!(-e $primersearchBin)){
	print STDERR "ERROR: One or more dependencies not found, please check that the paths of the following are real\n";
	print STDERR "primersearch       : $primersearchBin\n";
	pod2usage(-exitval => 1,  -verbose => 1);;
}


# Main program

### Get sequences names and sort them to create seqNameLst file for use with RECON
if ($opt_verbose){
	print STDERR "Making $opt_prefix.primerSearch file for primersearch.. ";
}

open (PRIMERSEARCHIN, ">$opt_prefix.primerSearchIn") || die ("Unable to open file for writing: $opt_prefix.primerSearchIn\n$!\n");
open (PRIMERFASTA, "<$opt_inFile") || die ("Unable to open file for reading: $opt_inFile\n$!\n");
while ($line= <PRIMERFASTA>){
	chomp $line;
	if ($line=~ m/# Scaffold:/){
		print PRIMERSEARCHIN "\n$line\n";
		next;
	}
	if ($line=~ m/^>(\S+_\d+\-\d+)_(\d+)_([FR])/){
		if ($3 eq 'F'){
			print PRIMERSEARCHIN "$1_$2_$3_";
			$tempString = "$1_$2_$3_";
			$tempString2= <PRIMERFASTA>;
			chomp $tempString2;
			$tempString2.= " ";
		}
		if ($3 eq 'R'){
			print PRIMERSEARCHIN "$2_$3 ";
			$tempString.= "$2_$3 ";
			$tempString2.= <PRIMERFASTA>;
			chomp $tempString2;
			print PRIMERSEARCHIN "$tempString2\n";
			$tempString.= $tempString2;
			$tempString=~ m/^(\S+) (\S+) (\S+)/;
			$primers{$1}= "$2\t$3";
		}
		next;
	}
}
close (PRIMERFASTA);
close (PRIMERSEARCHIN);

if ($opt_verbose){
	print STDERR "Done\n";
}

if ($opt_verbose){
	print STDERR "Running primersearch.. ";
}

$cmd= "$primersearchBin -seqall $opt_seqFile -infile $opt_prefix.primerSearchIn -mismatchpercent $opt_mismatch -outfile $opt_prefix.primerSearchOut";
system ($cmd) == 0 || die "Unable to run primersearch: $?\n$!\n";
if ($opt_verbose){
	print STDERR "Done\n";
}

if ($opt_verbose){
	print STDERR "Making outFile $prefix.primerCheck.. ";
}
open (PRIMERCHECKOUT, ">$opt_prefix.primerCheck") || die ("Unable to open file for writing: $opt_prefix.primerCheck\n$!\n");
open (PRIMERSEARCHOUT, "<$opt_prefix.primerSearchOut") || die ("Unable to open file for reading: $opt_prefix.primerSearchOut\n$!\n");
print PRIMERCHECKOUT "#primerPair_ID\tprimerF_sequence\tprimerR_sequence\trefSequenceHit_ID\trefSequenceHit_start\tamplicon_size\n";
while ($line= <PRIMERSEARCHOUT>){
	chomp $line;
	if ($line=~ m/^Primer name (\S+)/){
		$primerId= $1;
		$line= <PRIMERSEARCHOUT>;
		while ($line=~ m/^Amplimer (\d+)/){
			$ampNum= $1;
			$line= <PRIMERSEARCHOUT>;
			$line=~ m/Sequence: (\S+)/;
			print PRIMERCHECKOUT $primerId . "\t" . $primers{$primerId} . "\t" . $1 . "\t";
			$line= <PRIMERSEARCHOUT>;
			$line= <PRIMERSEARCHOUT>;
			$line=~ m/hits forward strand at (\d+) with/;
			print PRIMERCHECKOUT $1 . "\t";
			$line= <PRIMERSEARCHOUT>;
			$line= <PRIMERSEARCHOUT>;
			$line=~ m/Amplimer length: (\d+) bp/;
			print PRIMERCHECKOUT $1 . "\n";
			$line= <PRIMERSEARCHOUT>;
		}
	}
}
close (PRIMERSEARCHOUT);
close (PRIMERCHECKOUT);

if ($opt_verbose){
	print STDERR "Done\n";
}
exit (0);
