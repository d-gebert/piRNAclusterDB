#!/usr/bin/perl

use strict;
use warnings;

# Global constants
my $merge_dist = 10000;
# Global variables
# Options
$|=1; #Autoflush

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
    my($genus) = ($spec =~ /^(\w{1})/);
	my($suffx) = ($spec =~ /\s(\w{3})/);
    unless ($suffx) {
        ($suffx) = ($spec =~ /\s(\w{2})/);
    }
	my $sid = "$genus$suffx";
    # Check species
	unless ($gnm_list->{$spec}->[0]) { next; }
    print("$spec_folder\n");
    # Storage variable for all pic gtf data for current species
    my @pic_gtf_data = ();
    my $idx1 = 0;
    # Get paths to gtf files
	my @pic_gtf_files = `find $spec_folder -type f | grep /clusters.gtf`;
    # Go through each identified srna file
	foreach my $pic_gtf_file (@pic_gtf_files) {
	    $pic_gtf_file =~ s/\s+$//; #better chomp
	    #print("$pic_gtf_file\n");
        # Get file data
        $pic_gtf_data[$idx1] = get_gtf_fields($pic_gtf_file);
        $idx1++;
    }
    # Combine all pic loci in one list
    my @all_pic_locs = ();
    for (my $i = 0; $i < $idx1; $i++) {
        foreach my $pic_id (sort {$a <=> $b} keys %{$pic_gtf_data[$i]}) {
            my $chr = $pic_gtf_data[$i]->{$pic_id}->[0];
            my $src = $pic_gtf_data[$i]->{$pic_id}->[1];
            my $typ = $pic_gtf_data[$i]->{$pic_id}->[2];
            my $beg = $pic_gtf_data[$i]->{$pic_id}->[3];
            my $end = $pic_gtf_data[$i]->{$pic_id}->[4];
            #my($rds) = ($pic_gtf_data[$i]->{$pic_id}->[8] =~ /mapped sequence reads: (\d+);/);
            push(@all_pic_locs,[$chr,$src,$typ,$beg,$end,'.','.','.','.']);
        }
    }
    # Sort list of pic locs
    @all_pic_locs = sort { $a->[0] cmp $b->[0] || $a->[3] <=> $b->[3] } @all_pic_locs;
    # Output list of pic locs
    my $outfile_all = "$spec_folder/$spec_folder.piRNAclusters.all.gtf";
    my $out_all = open_outfile($outfile_all);
    foreach my $p (@all_pic_locs) {
        print($out_all "$p->[0]\t$p->[1]\t$p->[2]\t$p->[3]\t$p->[4]\t$p->[5]\t$p->[6]\t$p->[7]\t$p->[8]\n");
    }
    close($out_all);
    # Merge all pic loci to get consensus loci
    system("bedtools merge -d $merge_dist -i $outfile_all > $outfile_all.merged.bed");
    my $cons_pic_data = get_gtf_fields("$outfile_all.merged.bed");
    my @cons_pic_locs = ();
    foreach my $loc_id (sort {$a <=> $b} keys %{$cons_pic_data}) {
        my $chr = $cons_pic_data->{$loc_id}->[0];
        my $beg = $cons_pic_data->{$loc_id}->[1];
        my $end = $cons_pic_data->{$loc_id}->[2];
        my $len = $end-$beg+1;
        push(@cons_pic_locs,[$chr,'proTRAC','piRNA_cluster',$beg,$end,'.','.','.','.',$len]);
    }
    unlink("$outfile_all");
    unlink("$outfile_all.merged.bed");

    # Output merged consensus pic loc gtf
    my $outfile_gtf = "$spec_folder/$spec_folder.piRNAclusters.gtf";
    my $out_gtf = open_outfile($outfile_gtf);
    # Sort by locus length
    foreach my $p (sort {$b->[9] <=> $a->[9]} @cons_pic_locs) {
        print($out_gtf "$p->[0]\t$p->[1]\t$p->[2]\t$p->[3]\t$p->[4]\t$p->[5]\t$p->[6]\t$p->[7]\t$p->[8]\n");
    }
    close($out_gtf);

    # Get genome link
	my $genome_link_gb = $gnm_list->{$spec}->[14] ? $gnm_list->{$spec}->[14] : '';
	my $genome_link_rs = $gnm_list->{$spec}->[15] ? $gnm_list->{$spec}->[15] : '';
	my $genome_link = $genome_link_rs ? $genome_link_rs : $genome_link_gb;
    # Get genome file
	my($gnm_file_gz) = ($genome_link =~ /\/([^\/]+)$/);
	$gnm_file_gz .= '_genomic.fna.gz';
	my($gnm_file) = ($gnm_file_gz =~ /(.*)\.gz$/);
    # Extract genomic sequences
    my $chr_seqs = get_fasta_seqs("Genomes/$gnm_file",1);

    # Ouput merged consensus pic loc sequences
    my $outfile_fas = "$spec_folder/$spec_folder.piRNAclusters.fasta";
    #unless (-e $outfile_fas) {
        my $out_fas = open_outfile($outfile_fas);
        # Initialize index
        my $idx2 = 0;
        foreach my $p (sort {$b->[9] <=> $a->[9]} @cons_pic_locs) {
            my $chr = $p->[0];
            my $beg = $p->[3];
            my $end = $p->[4];
            my $len = $end-$beg+1;
            # Get pic sequence
            my $pic_seq = substr($chr_seqs->{$chr},($beg-1),$len);
            # Make a pic id
            $idx2++;
            my $pre = '0' x (4-length($idx2));
            my $pid = "${sid}_${pre}${idx2}";
            # Output pic sequence
            print($out_fas ">$pid $chr:$beg..$end\n");
    		my $fmt_seq = format_sequence($pic_seq,60);
    		print($out_fas "${$fmt_seq}");
        }
        close($out_fas);
    #}

    # Ouput merged consensus pic loc pirnas (unique)
    my $outfile_fas2 = "$spec_folder/$spec_folder.piRNAs.uniq.fasta";
    my $outfile_tbl = "$spec_folder/$spec_folder.piRNAclusters.rpm.tbl";
    #unless (-e $outfile_fas2) {
        # Get paths to sam files
    	my @pir_sam_files = `find $spec_folder -type f | grep ".uni.sam\$"`;
        # Get clustered piRNAs
        my @clustered_seqs = ();
        # Get pic expression per dataset
        my %rds_per_pic_per_sra = ();
        my @sra_ids = ();
        my %total_reads = ();
        # Go through each identified sam file
    	foreach my $pir_sam_file (@pir_sam_files) {
    	    $pir_sam_file =~ s/\s+$//; #better chomp
            # Get sra id
            my $sra_id = $pir_sam_file;
			$sra_id =~ s/\..+$//;
            $sra_id =~ s/^.*\///;
            push(@sra_ids,$sra_id);
    	    # Get piRNA read hit positions
            my($pir_hit_poss,$pir_reads) = get_sam_seq_hitposs($pir_sam_file);
            # Get total reads for current sra dataset
            foreach my $seq (keys %{$pir_reads}) {
                $total_reads{$sra_id} += $pir_reads->{$seq};
            }
            # Initialize index
            my $idx3 = 0;
            # Get clustered piRNAs for current sam file
            foreach my $p (sort {$b->[9] <=> $a->[9]} @cons_pic_locs) {
                # Get cluster coordinates
                my $chr = $p->[0];
                my $beg = $p->[3];
                my $end = $p->[4];
                # Make a pic id
                $idx3++;
                my $pre = '0' x (4-length($idx3));
                my $pid = "${sid}_${pre}${idx3}";
                # Go through each cluster position
                foreach my $pos ($beg..$end) {
                    # Check if cluster position overlaps with a piRNA read pos
                    if ($pir_hit_poss->{$chr}->{$pos}) {
                        # Save clustered piRNA sequences
                        my @seqs = @{$pir_hit_poss->{$chr}->{$pos}};
                        push(@clustered_seqs,@seqs);
                        # Get number of reads per cluster per sra dataset
                        foreach my $seq (@seqs) {
                            $rds_per_pic_per_sra{$pid}{$sra_id} += $pir_reads->{$seq};
                        }
                    }
                }
            }
        }
        # Ouput merged consensus pic loc piRNA reads (unique)
        my $out_fas2 = open_outfile($outfile_fas2);
        # Initialize index
        my $idx4 = 0;
        # Output unique clustered piRNAs
        my %clustered_seqs = map {$_ => 1} @clustered_seqs;
        foreach my $seq (sort keys %clustered_seqs) {
            $idx4++;
            print($out_fas2 ">${sid}_$idx4\n$seq\n");
        }
        # Ouput merged consensus pic loc expression per dataset
        my $out_tbl = open_outfile($outfile_tbl);
        # Title line
        print($out_tbl "Cluster_ID chr:beg..end");
        foreach my $sra_id (@sra_ids) {
            print($out_tbl "\t$sra_id");
        }
        print($out_tbl "\n");
        # Initialize index
        my $idx5 = 0;
        foreach my $p (sort {$b->[9] <=> $a->[9]} @cons_pic_locs) {
            my $chr = $p->[0];
            my $beg = $p->[3];
            my $end = $p->[4];
            # Make a pic id
            $idx5++;
            my $pre = '0' x (4-length($idx5));
            my $pid = "${sid}_${pre}${idx5}";
            print($out_tbl "$pid $chr:$beg..$end");
            # Go through each sra dataset
            foreach my $sra_id (@sra_ids) {
                my $rpm = int((($rds_per_pic_per_sra{$pid}{$sra_id}/$total_reads{$sra_id})*1000000)+0.5);
                print($out_tbl "\t$rpm");
            }
            print($out_tbl "\n");
        }
    #}

    # Ouput merged consensus pic loc pirnas (all)
    my $outfile_fas2_a = "$spec_folder/$spec_folder.piRNAs.fasta";
    #unless (-e $outfile_fas2_a) {
        # Get paths to sam files
    	my @pir_sam_files_a = `find $spec_folder -type f | grep ".all.sam\$"`;
        # Get clustered piRNAs
        my @clustered_seqs_a = ();
        # Go through each identified sam file
    	foreach my $pir_sam_file_a (@pir_sam_files_a) {
    	    $pir_sam_file_a =~ s/\s+$//; #better chomp
    	    # Get piRNA read hit positions
            my($pir_hit_poss_a,$pir_reads_a) = get_sam_seq_hitposs($pir_sam_file_a);
            # Get clustered piRNAs for current sam file
            foreach my $p (@cons_pic_locs) {
                my $chr = $p->[0];
                my $beg = $p->[3];
                my $end = $p->[4];
                # Go through each cluster position
                foreach my $pos ($beg..$end) {
                    if ($pir_hit_poss_a->{$chr}->{$pos}) {
                        my @seqs = @{$pir_hit_poss_a->{$chr}->{$pos}};
                        push(@clustered_seqs_a,@seqs);

                    }
                }
            }
        }
        # Ouput merged consensus pic loc piRNA reads (all)
        my $out_fas2_a = open_outfile($outfile_fas2_a);
        # Initialize index
        my $idx6 = 0;
        # Output all clustered piRNAs
        my %clustered_seqs_a = map {$_ => 1} @clustered_seqs_a;
        foreach my $seq (sort keys %clustered_seqs_a) {
            $idx6++;
            print($out_fas2_a ">${sid}_$idx6\n$seq\n");
        }
    #}

}

exit;

################################# subroutines #################################

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

sub get_gtf_fields {
	# Take name of pic tree file
	my($infile) = @_;
	# Set 0 as initial id
	my $id = 0;
	# Get file data
	my @in_data = get_file_data_array($infile);
	# Global tree hashes for each species
	my %data_fields = ();
	# Go through file data
	foreach my $line (@in_data) {
        # Skip title line
        if ($line =~ /#/) { next; }
        # Get line data
        my @d = split(/\t/,$line);
        # Get id
        $id++;
        # Save data fields
        @{$data_fields{$id}} = @d;
    }
	# Return data fields
	return \%data_fields;
}

sub get_sam_seq_hitposs {
	# Take name of sam file
	my($samfile) = @_;
	# Initialize variables
	my %sam_seqs = ();
    my %sam_seq_rds = ();
	# Open all mappers sam file
	my $sam = open_infile($samfile);
	# Go through each line
	while (my $line = <$sam>) {
		$line =~ s/\s+$//; #better chomp
		# Get map data
		my @d = split(/\t/,$line);
		my $rds = $d[0];
		my $flg = $d[1];
		my $str = $flg == 16 ? '-' : '+';
		my $seq = $d[9];
        my $chr = $d[2];
        my $pos = $d[3];
		# Reverse sequence if on minus strand
		if ($str eq '-') {
			$seq = reverse($seq);
			$seq =~ tr/ATGC/TACG/;
		}
		# Save sequence and number of reads
		push(@{$sam_seqs{$chr}{$pos}},$seq);
        $sam_seq_rds{$seq} = $rds;
	}
	return \(%sam_seqs,%sam_seq_rds);
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

# Save fasta data as hash
# Usage: my $sequences = get_fasta_seqs($infile);
sub get_fasta_seqs {
	# Take fasta file name
	my($file,$short) = @_;
	# Variables
	my $name = '';
	my %sequences = ();
	# Open file
	my $in = open_infile($file);
	# Extract sequence and save with name
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		if ($line =~ /^>/) {
			($name) = ($line =~ />(.*)$/);
			if ($short) { $name =~ s/\s.*// }
		} else {
			$sequences{$name} .= $line;
		}
	}
	return \%sequences;
}

# Format sequence with given line length
# Usage: format_sequence($sequence,$length);
sub format_sequence {
	# Take sequence and line length
    my($sequence, $length) = @_;
    my $fmtd_seq = '';
    # Print sequence in lines of $length
    for (my $pos = 0 ; $pos < length($sequence) ; $pos += $length) {
        $fmtd_seq .= substr($sequence, $pos, $length)."\n";
    }
    # Return formatted sequence reference
    return \$fmtd_seq;
}
