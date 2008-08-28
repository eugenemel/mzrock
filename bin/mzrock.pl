#This is a wrapper script normally excuted by monitorfolder program
#the monitor folder runs this script with --projectfolder option
#where projetfolder contains Thermo raw files

#This script finds raw files, converts them to mzCSV file format
#and executes R script to process files

#EXECUTION on COMMAND LINE
#   perl mzrock.pl --projectfolder YourFolderName

use strict;
use Getopt::Long;

#location of folders and executables,  NOTE USE of "\\" in pathnames under windows  OSes
my $INSTALLDIR="C:\\mzRock\\bin";
my $REXECUTABLE="C:\\Program Files\\R\\R-2.7.0\\bin\\R.exe";
my $RAW2MZCSV="$INSTALLDIR\\raw2mzCSV.exe";
my $RSCRIPT="$INSTALLDIR\\mzRock.auto.R";

#name of the folder with RAW files
my %OPT; &GetOptions(\%OPT, "projectfolder:s");     #--projectfolder option must be specified 
my $PROJECTFOLDER = $OPT{'projectfolder'};

#check if we have everything to start processing
if ( not -d $PROJECTFOLDER ) { die "Could not find project folder $PROJECTFOLDER\n"; }
if ( not -d $INSTALLDIR)  { die "Could not find install folder $INSTALLDIR\n"; }
if ( not -e $REXECUTABLE)  { die "Could not find R  in $REXECUTABLE\n"; }
if ( not -e $RAW2MZCSV )  { die "Could not find RAW2mzCSV exectuble $RAW2MZCSV\n"; }

#start log
open(LOG, ">$PROJECTFOLDER/log.txt") or warn "Can't write to $PROJECTFOLDER\log.txt file";
chdir($PROJECTFOLDER) or die "Can't change to folder $PROJECTFOLDER\n";
opendir(DIR, ".") or die "Can't read files from folder $PROJECTFOLDER\n";

#find Thermo .RAW files
foreach my $filename (readdir(DIR)) { 
	next unless ( $filename  =~ /raw$/i ); 
	
	my $id = $filename; 
	$id =~ s/\.raw$//i;
	
    #generate mzCSV file
	my $mzCSVFile = "$id.mzCSV";    
	print LOG "Writing $mzCSVFile\n";
	system "$RAW2MZCSV \"$filename\" >  $mzCSVFile";	#Xcalibar must be installed for this to work
}

#execute R
system("\"$REXECUTABLE\" BATCH --no-save < $RSCRIPT");
close(LOG);
