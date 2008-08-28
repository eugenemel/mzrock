#monitor Data Directory for new raw files

use strict;
use File::stat;
use Data::Dumper;

my $CHECK_ATTEMPTS=3;

my $dirToMonitor = '..\DataSets';
my $mzRockPL	=  'perl.exe mzrock.pl';

my %fileList1;
&getFileList($dirToMonitor,\%fileList1);
my %alteredFiles;

while(1) {
    print STDERR "Checking $dirToMonitor for new files\n";
    my %fileList2;
    &getFileList($dirToMonitor,\%fileList2);    
    &compareLists(\%fileList1, \%fileList2);
    &processAlteredFiles();
    sleep(10);
}

sub processAlteredFiles() {
	my %alteredDirectories;
    foreach my $file (  keys %alteredFiles ) {
        print "Processing: $file\n";
		if ( &monitorSingleFile($file) ) {
			my $sb = stat("$file");
			$fileList1{$file} = $sb->size . $sb->mtime;

			#windows looks for longest string ending with backslash \
			my $dir = undef;
			if ( $file =~ /(.*\\)/) { $dir = $1; }	
			$alteredDirectories{$dir}++;

		}
        delete $alteredFiles{$file};
    }

	foreach my $dir ( keys %alteredDirectories ) {
		if ( &monitorDirectory($dir) ) { 
			$dir =~ s/\\$//;
			print STDERR "Processing folder: $mzRockPL --projectfolder \"$dir\"\n";
			system("$mzRockPL --projectfolder \"$dir\"");
			&getFileList("$dir",\%fileList1);	#syncrozize filelists
		}
	}
}

sub compareLists() {
    my $list1 = shift;
    my $list2 = shift;

    foreach my $file ( keys %$list2 ) {
        if ( not exists $list1->{$file}) {
            print "New file $file\n";
            $alteredFiles{$file}= $list2->{$file};
        } elsif ( $list1->{$file} ne $list2->{$file} ) {
            print "Changed file $file $list1->{$file} vs $list2->{$file}\n";
            $alteredFiles{$file}= $list2->{$file};
        }
    }
}


sub getFileList() {
    my $dir = shift;
    my $fileList = shift;

    opendir(DIR, $dir) or warn "getFileList: Can't open directory $dir\n";
    my @contents = readdir(DIR);

#	print "getFileList: $dir\n";

    foreach my $f (@contents) {
        next if $f eq '.' or $f eq '..';
        if (-d "$dir\\$f" ) {
            &getFileList( "$dir\\$f", $fileList );
        } elsif ( -f "$dir\\$f" and $f =~ /raw$/i ) {
            my $sb = stat("$dir\\$f");
            $fileList->{"$dir\\$f"} = $sb->size . $sb->mtime;
			#print "\t\t $f\n";

        }
    }
	closedir(DIR);
}

sub monitorDirectory() {
	my $dir = shift;
	my $noChangesCount;
	my $lastDirSize;

	#monitor directory for fize size changes
	print STDERR "Monitoring $dir\n";
	while($noChangesCount < $CHECK_ATTEMPTS) { 
		return 0 if not -e $dir;
		opendir(DIR, $dir) or warn "monitorDirectory: Can't open directory $dir\n";

		my $dirSize;
		#calculate total size of the directory
		foreach my $file ( readdir(DIR) ) {
			next unless -e "$dir\\$file";
		   	my $sb = stat("$dir\\$file");
			$dirSize += $sb->size;
			sleep(1);
		}
		closedir(DIR);
	
		print STDERR "\t total dir size $dirSize\n";
		$noChangesCount ++ if ($lastDirSize == $dirSize );
		$lastDirSize = $dirSize;	#update to new directory size
	}
	return 1;
}

sub monitorSingleFile() {
	my $file = shift;
	my $lastsb = stat("$file"); sleep(1);
	my $noChangesCount;
	while($noChangesCount < 1) {
		return 0 if not -e $file;
        my $sb = stat("$file");
		if ($sb->size == $lastsb->size and $sb->mtime == $lastsb->mtime ) {
			$noChangesCount++;
			print "\t\t no changes count=$noChangesCount\n";
		} else {
			$noChangesCount = 0;
		}

		$lastsb = $sb;
		sleep(1);
	}
	return 1;
}
