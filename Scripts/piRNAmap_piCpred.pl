#!/usr/bin/perl

use strict;
use warnings;

# Global constants
# Global variables
# Options
my $no_rn = 0;
$|=1; #Autoflush

# Program names/paths
my $unitas   = 'unitas';
my $fq_dump  = 'fastq-dump';
my $gffread  = 'gffread';
my $protrac  = 'proTRAC';
my $fastqc   = 'fastqc';
my $bowtie   = 'bowtie';
my $samtools = 'samtools';
my $bedtools = 'bedtools';
my $pTscript = "proTRAC_2.4.4.pl";
my $p_val = 0.01;
my $p_opt = "-clsize 5000 -pimin 23 -pimax 31 -1Tor10A 0.5 -1Tand10A 0.33 -clstrand 0.0 -clsplit 1.0 -distr 5.0-90.0 -spike 50-1000 -nomotif -pdens $p_val";
# Test dependencies
my $depends = '';
foreach my $prog ($unitas,$fq_dump,$gffread,$protrac) {
	system("$prog -h 2> progtest_err > progtest_log");
	open(my $fh, '<', "progtest_err");
	my $frst_ln = <$fh>;
	if ($frst_ln && $frst_ln =~ /not found/) { $depends .= "$prog\n" }
	unlink('progtest_err','progtest_log');
}
if ($depends) { die("\nERROR: Cannot start program.\nFollowing dependencies need to be installed:\n$depends\n"); }

# Program name
print("\n--- $0 ---\n");
# Usage messaging
my $USAGE = "perl $0 genomes_table [species (optional)]\n";
unless ($ARGV[0]) {
   die("\nUsage: $USAGE\n");
}
# Collect command line arguments
my $genomes_table = $ARGV[0];

# Get genome list
my $gnm_list = get_tab_fields($genomes_table);

# Get species folder(s)
my @spec_folders = $ARGV[1] ? ($ARGV[1]) : `ls -1d *_*/`;

# Clean up species folder name(s)
foreach my $spec_folder (@spec_folders) {
	$spec_folder =~ s/\s+$//; #better chomp
	$spec_folder =~ s/\///;
}

# Go through each given species folder
foreach my $spec_folder (@spec_folders) {
	# Get species name
	my $spec = $spec_folder;
	$spec =~ s/_/ /;

	# Get genome link
	my $genome_link_gb = $gnm_list->{$spec}->[14] ? $gnm_list->{$spec}->[14] : '';
	my $genome_link_rs = $gnm_list->{$spec}->[15] ? $gnm_list->{$spec}->[15] : '';
	my $genome_link = $genome_link_rs ? $genome_link_rs : $genome_link_gb;

	# Get paths to fasta files
	my @noann_files = `find $spec_folder -type f | grep no-annotation.fas`;

	## Download and modify genome reference data
	# Get genome file
	my($gnm_file_gz) = ($genome_link =~ /\/([^\/]+)$/);
	$gnm_file_gz .= '_genomic.fna.gz';
	my($gnm_file) = ($gnm_file_gz =~ /(.*)\.gz$/);
	# Check if genome file exists, otherwise download and/or unzip
	unless (-e "Genomes/$gnm_file") {
	    unless (-e "Genomes/$gnm_file_gz") {
	        system("wget $genome_link/$gnm_file_gz -P Genomes >>ftp_log 2>&1");
	        unless (-e "Genomes/$gnm_file_gz") { next; }
	    }
	    system("gunzip -k Genomes/$gnm_file_gz >>sap_log 2>&1");
	}
	# Clean genome file (shorten fasta headers)
	clean_fasta_header("Genomes/$gnm_file") unless $no_rn;
	# Get repeatmasker file
	my($rep_file_gz) = ($genome_link =~ /\/([^\/]+)$/);
	$rep_file_gz .= '_rm.out.gz';
	my($rep_file) = ($rep_file_gz =~ /(.*)\.gz$/);
	# Check if repeatmasker file exists, otherwise download and/or unzip
	unless (-e "Genomes/$rep_file") {
	    unless (-e "Genomes/$rep_file_gz") {
	        system("wget $genome_link/$rep_file_gz -P Genomes >>ftp_log 2>&1");
	        #unless (-e "Genomes/$rep_file_gz") { next; }
	    }
		if (-e "Genomes/$rep_file_gz") {
			system("gunzip -k Genomes/$rep_file_gz >>sap_log 2>&1");
		}
	}
	# Get gff file
	my($gff_file_gz) = ($genome_link =~ /\/([^\/]+)$/);
	$gff_file_gz .= '_genomic.gff.gz';
	my($gff_file) = ($gff_file_gz =~ /(.*)\.gz$/);
	# Check if gff file exists, otherwise download and/or unzip
	unless (-e "Genomes/$gff_file") {
	    unless (-e "Genomes/$gff_file_gz") {
	        system("wget $genome_link/$gff_file_gz -P Genomes >>ftp_log 2>&1");
	    }
	    if (-e "Genomes/$gff_file_gz") {
	        system("gunzip -k Genomes/$gff_file_gz >>sap_log 2>&1");
	    } elsif (-e "Genomes/$gff_file") {
	        open(my $in, '>', "Genomes/$gff_file");
	        close($in);
	    }
	}
	# Transform gff to gtf
	my $gtf_file = $gff_file;
	$gtf_file =~ s/\.gff/\.gtf/;
	system("$gffread -T Genomes/$gff_file -o Genomes/$gtf_file >>sap_log 2>&1");
	# Get assembly report file
	my($inf_file) = ($genome_link =~ /\/([^\/]+)$/);
	$inf_file .= '_assembly_report.txt';
	# Check if assembly report file exists, otherwise download
	unless (-e "Genomes/$inf_file") {
	    system("wget $genome_link/$inf_file -P Genomes >>ftp_log 2>&1");
	}
	# Check presence of chromosome names
	my $chr_names_pres = check_chr_names("Genomes/$inf_file");
	# If chromosome names are present, change ids to names
	if ($chr_names_pres) {
	    # Extract chromosome ids & names
	    my $chr_names_gb = get_assembly_ids("Genomes/$inf_file",'gb');
	    my $chr_names_rs = get_assembly_ids("Genomes/$inf_file",'rs');
	    #my $chr_names = $genome_link_rs ? $chr_names_rs : $chr_names_gb;
	    my $chr_names = $chr_names_rs ? $chr_names_rs : $chr_names_gb;
	    # Change chromosome names
	    change_chr_gnm("Genomes/$gnm_file",$chr_names) unless $no_rn;
	    change_chr_rep("Genomes/$rep_file",$chr_names) unless $no_rn;
	    change_chr_gff("Genomes/$gff_file",$chr_names) unless $no_rn;
	    change_chr_gff("Genomes/$gtf_file",$chr_names) unless $no_rn;
	}
	# Create bowtie genome index if it does not yet exist
	my($gnm_idx) = ($gnm_file =~ /(.*)\.fna$/);
	unless (-e "Genomes/$gnm_idx.rev.2.ebwt" || -e "Genomes/$gnm_idx.rev.2.ebwtl") {
	    system("bowtie-build Genomes/$gnm_file Genomes/$gnm_idx");
	}

	# Go through each identified srna file
	foreach my $srna_fas (@noann_files) {
	    $srna_fas =~ s/\s+$//; #better chomp
	    print("$srna_fas\n");
		# All mappers
	    # Name bam file
	    my $bam_file_a = $srna_fas;
	    $bam_file_a =~ s/\.fas/\.all.bam/;
	    # Map filtered reads to genome
	    unless (-e "$bam_file_a") {
	        system("bowtie -k 1 -v 0 Genomes/$gnm_idx -f $srna_fas -S | samtools view -S -b -F 4 | samtools sort > $bam_file_a");
	    }
	    # Name sam file
	    my $sam_file_a = $bam_file_a;
	    $sam_file_a =~ s/\.bam/\.sam/;
	    # Convert bam to sam
	    unless (-e "$sam_file_a") {
	        system("$samtools view $bam_file_a > $sam_file_a");
	    }
	    # Unique mappers
	    # Name bam file
	    my $bam_file_u = $srna_fas;
	    $bam_file_u =~ s/\.fas/\.uni.bam/;
	    # Map filtered reads to genome
	    unless (-e "$bam_file_u") {
	        system("bowtie -m 1 -v 0 Genomes/$gnm_idx -f $srna_fas -S | samtools view -S -b -F 4 | samtools sort > $bam_file_u");
	    }
	    # Name sam file
	    my $sam_file_u = $bam_file_u;
	    $sam_file_u =~ s/\.bam/\.sam/;
	    # Convert bam to sam
	    unless (-e "$sam_file_u") {
	        system("$samtools view $bam_file_u > $sam_file_u");
	    }
	    # Predict pirna clusters
	    my $tisdir = $srna_fas;
	    $tisdir =~ s/\/[^\/]+.no-annotation.fas//;
	    my($sra_id) = ($srna_fas =~ /\/([^\/]+).no-annotation.fas/);
	    my $map_d = "$sam_file_u";
	    my $gnm_d = "Genomes/$gnm_file";
	    my $rep_d = "Genomes/$rep_file";
	    my $gen_d = -s "Genomes/$gtf_file" ? "Genomes/$gtf_file" : "Genomes/$gff_file";
	    # Check if pic prediction already exists
	    #print("$sra_id\n");
	    opendir(my $pro_dh, $tisdir);
	    unless (grep { "$tisdir/$_" =~ /proTRAC_$sra_id.no*/ } readdir $pro_dh) {
	        close($pro_dh);
	        system("perl $pTscript -m $map_d -g $gnm_d -repeatmasker $rep_d -geneset $gen_d $p_opt >>pT_log 2>&1");
	        # Move output folder to sra file directory
	        if (my @pt_dir = grep(-d, <*proTRAC_$sra_id.no*>)) {
	            system("mv $pt_dir[0] $tisdir >>sap_log 2>&1");
	        }
	    }
	}
	# Copy genome/annotation data
	unless (-e "$spec_folder/$gnm_file") {
		system("cp Genomes/$gnm_file $spec_folder/$gnm_file");
	}
	unless (-e "$spec_folder/$rep_file") {
		system("cp Genomes/$rep_file $spec_folder/$rep_file");
	}
	unless (-e "$spec_folder/$gff_file") {
		system("cp Genomes/$gff_file $spec_folder/$gff_file");
	}
}


exit;

################################# subroutines #################################

sub clean_fasta_header {
	# Take fasta file name
	my($file) = @_;
	# Open in and out files
	my $gnm = open_infile($file);
	my $new = open_outfile("$file.clean");
	# Clean fasta headers
	while (my $line = <$gnm>) {
		$line =~ s/\s+$//; #better chomp
		if ($line =~ /^>/) { $line =~ s/\s.*//; }
		print($new "$line\n");
	}
	# Close files
	close($gnm);
	close($new);
	# Delete old file and replace with cleaned file
	unlink($file);
	system("mv $file.clean $file >>sap_log 2>&1");
}

sub check_chr_names {
	# Take name of report file and source (gb|rs)
	my($infile) = @_;
	# Get file data
	my @in_data = get_file_data_array($infile);
	# Global tree hashes for each species
	my $chr_names_present = 1;
	# Go through file data
	foreach my $line (@in_data) {
		# Skip comment lines
		if ($line =~ /#/) { next; }
		# Get line data
		my @d = split(/\t/,$line);
		my $name = $d[0];
		if ($name eq 'na') {
			$chr_names_present = 0;
			last;
		}
	}
	# Return data fields
	return $chr_names_present;
}

sub get_assembly_ids {
	# Take name of report file and source (gb|rs)
	my($infile,$src) = @_;
	# Get file data
	my @in_data = get_file_data_array($infile);
	# Global tree hashes for each species
	my %seq_names = ();
	# Go through file data
	foreach my $line (@in_data) {
		# Skip comment lines
		if ($line =~ /#/) { next; }
		# Get line data
		my @d = split(/\t/,$line);
		my $name = $d[0];
		my $ncbi = '';
		if ($src eq 'gb') {
			$ncbi = $d[4];
		} elsif ($src eq 'rs') {
			$ncbi = $d[6];
		} else {
			$ncbi = $d[4];
		}
		# Save data fields
		$seq_names{$ncbi} = $name;
	}
	# Return data fields
	return \%seq_names;
}

sub change_chr_gff {
	# Take name of genome file
	my($infile,$chr_names) = @_;
	# Return 0 if input file does not exist
	unless (-e $infile) { return; }
	# Open input file
	my $in = open_infile($infile);
	# Open output file
	my $outfile = $infile.'.rnm';
	my $out = open_outfile($outfile);
	# Go through each line of file
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Entry line
		if ($line !~ /^#/) {
			# Get refseq id
			my @d = split(/\t/,$line);
			my $chr_id = $d[0];
			# Change refseq id to chromosome name
			if ($chr_names->{$chr_id}) {
				$line =~ s/$chr_id/$chr_names->{$chr_id}/;
			}
		}
		print($out "$line\n");
	}
	# Close files
	close($in);
	close($out);
	# Delete old file and replace with cleaned file
	unlink($infile);
	system("mv $outfile $infile >>sap_log 2>&1");
}

sub change_chr_rep {
	# Take name of genome file
	my($infile,$chr_names) = @_;
	# Return 0 if input file does not exist
	unless (-e $infile) { return; }
	# Open input file
	my $in = open_infile($infile);
	# Open output file
	my $outfile = $infile.'.rnm';
	my $out = open_outfile($outfile);
	# Go through each line of fasta file
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Line starts with sw score
		if ($line =~ /^\s*\d+/) {
			# Remove leading whitespace
			my $line_mod = $line;
			$line_mod =~ s/^\s*//;
			# Get line data
			my @d = split(/\s+/,$line_mod);
			# Get chromosome id
			my $chr_id = $d[4];
			# Substitute chromosome id with name
			if ($chr_names->{$chr_id}) {
				$line =~ s/$chr_id/$chr_names->{$chr_id}/;
			}
		}
		# Print to file
		print($out "$line\n");
	}
	# Close files
	close($in);
	close($out);
	# Delete old file and replace with cleaned file
	unlink($infile);
	system("mv $outfile $infile >>sap_log 2>&1");
}

sub change_chr_gnm {
	# Take name of genome file
	my($infile,$chr_names) = @_;
	# Open input file
	my $in = open_infile($infile);
	# Open output file
	my $outfile = $infile.'.rnm';
	my $out = open_outfile($outfile);
	# Go through each line of fasta file
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Fasta header line
		if ($line =~ /^>/) {
			# Substitute chromosome id with name
			my($chr_id) = ($line =~ />(.*)/);
			my $chr_nm = $chr_names->{$chr_id} ? $chr_names->{$chr_id} : $chr_id;
			# Print to file
			print($out ">$chr_nm\n");
		} else {
			# Print to file
			print($out "$line\n");
		}
	}
	# Close files
	close($in);
	close($out);
	# Delete old file and replace with cleaned file
	unlink($infile);
	system("mv $outfile $infile >>sap_log 2>&1");
}

sub get_tab_fields {
	# Take name of pic tree file
	my($infile, $id_i) = @_;
	# Set 0 as default index
	unless ($id_i) { $id_i = 0 }
	# Get file data
	my @in_data = get_file_data_array($infile);
	# Global tree hashes for each species
	my %data_fields = ();
	# Go through file data
	foreach my $line (@in_data) {
        # Skip title line
        if ($line eq $in_data[0]) { next; }
        # Get line data
        my @d = split(/\t/,$line);
		s/"//g for @d;
        # Get id
        my $id = $d[$id_i];
        # Save data fields
        @{$data_fields{$id}} = @d;
    }
	# Return data fields
	return \%data_fields;
}

# Open input file
# Usage: my $in = open_infile($infile);
sub open_infile {
	# Take input file name
    my($file) = @_;
    # Open input file
    my $fh;
    if ($file =~ /.gz$/) {
		open($fh, "gunzip -c $file |") or die("Cannot open file '$file': $!\n");
	} else {
    	open($fh, '<', $file) or die("Cannot open file '$file': $!\n");
    }
    # Return filehandle
    return $fh;
}

# Open output file
# Usage: my $out = open_outfile($outfile);
sub open_outfile {
	# Take output file name
    my($file) = @_;
    # Open output file
    open(my $fh, '>', $file) or die("Cannot open file '$file': $!\n");
    # Return filehandle
    return $fh;
}

# Extract file data and save in array
# Usage: my @filedata = get_file_data_array($file);
sub get_file_data_array {
	# Take input file name
    my($file,$ref_opt) = @_;
    my @filedata = ();
    $ref_opt = 0 unless $ref_opt;
	# Open input file
    my $fh = open_infile($file);
	# Extract lines and save in array
    while (my $line = <$fh>) {
    	$line =~ s/\s+$//; #better chomp
    	push(@filedata, $line);
    }
	# Close file
    close($fh) or die("Unable to close: $!\n");
	# Return array containing file data
    if ($ref_opt) {
    	return \@filedata;
    } else {
    	return @filedata;
    }
}
