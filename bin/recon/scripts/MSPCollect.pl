#!/usr/bin/perl

if (!@ARGV) {
    die "usage: MSPCollect BLAST_output_file\n";
}

$score_cutoff = 0;
$iden_cutoff = 0;

$score = -1;

open (BLAST, "<$ARGV[0]") || die "usage: MSPCollect BLAST_output_file\nCan not open the BLAST_output_file $ARGV[0]\n";
while ($line=<BLAST>) {
	# Megablast line def: QUERY	SBJCT	%Id	AlnLength	Mismatch	GAP	QRYstart	QRYend	SBJCTstart	SBJCTend	E-val	BITSCORE
	if ($line=~ m/^(\S+)\t(\S+)\t(\S+)\t\S+\t\S+\t\S+\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+\t(\S+)$/){
		$query=$1;
		$sbjct=$2;
		$iden=$3;
		$head_query=$4;
		$tail_query=$5;
		$head_sbjct=$6;
		$tail_sbjct=$7;
		$score=$8;
		if ($score =~ /(\S+)[eE]\+(\d+)$/) {
			$score = $1 * exp($2*log(10));
		}
		if ($score > 0 && $score >= $score_cutoff && $iden >= $iden_cutoff && !($query==$sbjct && $head_query==$head_sbjct && $tail_query==$tail_sbjct)){
			printf("%06d %03d %05d %05d %s %05d %05d %s \n", $score, $iden, $head_query, $tail_query, $query, $head_sbjct, $tail_sbjct, $sbjct);
		}
	}
}
close (BLAST);
