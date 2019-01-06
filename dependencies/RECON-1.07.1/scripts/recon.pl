#!/usr/bin/perl

$path = $ARGV[3] . "/bin";

if (!$path) {
    die "Please modifiy this script to specify the path to the binaries.\n";
}

if (@ARGV < 2) {
    die "usage: recon seq_name_list_file MSP_file integer\nSee 00README for details.\n";
}

open(SEQ, "$ARGV[0]") || die "usage: recon seq_name_list_file MSP_file integer\nCan not open the seq_name_list_file $ARGV[0].\nSee 00README for details.\n";
close (SEQ);

open(MSP, "$ARGV[1]") || die "usage: recon seq_name_list_file MSP_file integer\nCan not open the MSP_file $ARGV[1].\nSee 00README for details.\n";
close (MSP);


unlink ("core");

unlink glob ("summary/* summary/*.*");
rmdir ("summary");
mkdir ("summary");



# prepare images

unlink glob ("images/* images/*.*");
rmdir ("images");
mkdir ("images");

`$path/imagespread $ARGV[0] $ARGV[1] $ARGV[2]`;

if ($?) {die "imagespread failed. Exit code $?\n";}

if (@ARGV < 3) {
    $sect = 1;
} else {
    $sect = $ARGV[2];
}

for ($i=1; $i<=$sect; $i++) {
   $spread = "images/spread" . $i;   
   `sort -k 3,3 -k 4n,4n -k 5nr,5nr $spread >> images/images_sorted`;
   if ($?) {die "sort failed for $spread.\n";}
}

unlink glob ("images/spread*");


# initial definition of elements

unlink glob ("ele_def_res/* ele_def_res/*.*");
rmdir ("ele_def_res");
mkdir("ele_def_res");

`$path/eledef $ARGV[0] $ARGV[1] single`;
if ($?) {die "eledef failed. Exit code $?\n";}



# re-defining elements

unlink glob ("ele_redef_res/* ele_redef_res/*.*");
rmdir ("ele_redef_res");
mkdir ("ele_redef_res");

unlink ("tmp", "tmp2");
symlink ("ele_def_res", "tmp");
symlink ("ele_redef_res", "tmp2");

`$path/eleredef $ARGV[0]`;
if ($?) {die "eleredef failed. Exit code $?\n";}

unlink ("tmp", "tmp2");


# re-defining edges

unlink glob ("edge_redef_res/* edge_redef_res/*.*");
rmdir ("edge_redef_res");
mkdir ("edge_redef_res");

unlink ("tmp", "tmp2");
symlink ("ele_redef_res", "tmp");
symlink ("edge_redef_res", "tmp2");

`$path/edgeredef $ARGV[0]`;
if ($?) {die "edgeredef failed. Exit code $?\n";}



# famdef

unlink ("tmp", "tmp2");
symlink ("edge_redef_res", "tmp");

`$path/famdef $ARGV[0]`;
if ($?) {die "famdef failed. Exit code $?\n";}

unlink ("tmp");

