#!/usr/bin/perl
#-------------------------------------------------------
#
#      Dessins des HR de mes grilles de modèles 
#
#--------------------------------------------------
# D.C., janvier 2005.
#
# Initialisation :
$file_in = $ARGV[0];
if ( $file_in eq "" ){
    print "Usage: $0 file_in\n";
    exit(0);
}
$title = "'LMC models, overshooting 0.4 Hp'";
$label = "LMC_0.4_HR";
# Ouverture des fichiers :
open(FILE_IN, "< $file_in") 
    or die " *** le fichier : $file_in n'est pas lisible : $!";

while ( defined($ligne=<FILE_IN>) ) {
    $ligne=~ s/\n//;
    $filenames= $filenames." "."'$ligne' using 4:3 with lines,";
}
close(FILE_IN);
$long=length($filenames);
$filenames= substr($filenames,0,$long-1);
print "$long $filenames\n";

# Dessin avec GNUPLOT
open(PLOT, "| gnuplot") or goto PBGNUPLOT;
    print PLOT qq{set term png\n};
    print PLOT qq{set title $title\n};
    print PLOT qq{set xrange [] reverse\n};
    print PLOT qq{set output '$label.png'\n};
    print PLOT qq{set xlabel "Log Teff"\n};
    print PLOT qq{set ylabel "Log L/Lsun"\n};
    print PLOT qq{set nokey\n};
    print PLOT qq{plot $filenames\n};
    print PLOT qq{set size 0.5,0.5\n};
$label=$label."_small";
    print PLOT qq{set output '$label.png'\n};
    print PLOT qq{plot $filenames\n};
close(PLOT);

exit(0);

PBGNUPLOT:

    print " ==>> GNUplot n'est pas installé sur cette machine !\n";
exit(0);
