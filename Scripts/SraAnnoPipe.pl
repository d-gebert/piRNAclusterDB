#!/usr/bin/perl
use strict;
use warnings;

# Program names/paths
my $unitas  = 'unitas';
my $fq_dump = 'fastq-dump';
my $gffread = 'gffread';
my $protrac = 'proTRAC';
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

# Global constants
my @tissues_kws = (
	"testis","sperm","ovary","egg","oocyte","embryo","gonad","germline","reproductive",
	"brain","muscle","thorax","head","gill","soma","hemolymph","blastula","trochophore",
	"mantle","labial palps","d-shape","umbo","pediveliger","pereon","pleon"
);
my $tis_runs_max = 30;
my $tis_pi_runs_max = 10;
my $data_dir = "SRA_data";
my $anno_dir = "UNI_data";
my $refseq_file = "ref_ncRNAs.fas";
my $pwd = `pwd`;
my $n_test_seqs = 100000;
my $min_len = 16;
my $max_len = 36;
my $no_pt = 0;
my $no_dl = 0;
my $no_rn = 0;
# Global variables
# Options
$|=1; #Autoflush

# Program name
print("\n--- $0 ---\n");
# Usage messaging
my $USAGE = "perl $0 <sraruntable_clean.txt> <info_clean.txt> <AnimalGenomes.txt> <genus_species>\n";
unless ($ARGV[0]&&$ARGV[1]&&$ARGV[2]) {
	die("\nUsage: $USAGE\n");
}
# Collect command line arguments
my $sra_run_table = $ARGV[0];
my $sra_inf_table = $ARGV[1];
my $genomes_table = $ARGV[2];
my $species = $ARGV[3] ? $ARGV[3] : '';

# Get index lists
my $idxs = get_field_indexes($sra_run_table);
# Get sra run list
my $sra_list = get_tab_fields($sra_run_table,$idxs->{'Run'});
# Get info list
my $inf_list = get_tab_fields($sra_inf_table);
# Get genome list
my $gnm_list = get_tab_fields($genomes_table);

##### Select sra datasets
# Save selected run ids
my %selected_rids = ();
# Go through each species
foreach my $spec (sort {$inf_list->{$a}->[1] cmp $inf_list->{$b}->[1]} keys %{$inf_list}) {
	# Go through each run id
	foreach my $rid (sort keys %{$sra_list}) {
		# Find species
		if ($sra_list->{$rid}->[$idxs->{'Organism'}] eq $spec) {
			# Check if tissue is wanted
			$sra_list->{$rid}->[$idxs->{'tissue'}] = '' unless $sra_list->{$rid}->[$idxs->{'tissue'}];
			my @tiss = ();
			if (@tiss = grep {$sra_list->{$rid}->[$idxs->{'tissue'}] =~ /$_/i} @tissues_kws) {
				# Save run id
				push(@{$selected_rids{$spec}{$tiss[0]}},$rid);
			}
		}
	}
}

##### Find sra datasets containing pirna sequences
# Get unitas species list
system("$unitas -supp_spec | grep _ | sed 's/ - //g' >unitas_specs");
my @unitas_specs = get_file_data_array('unitas_specs');
# Create directory for data storage
mkdir($data_dir);
my $i_spec = 0;
# Go through each species
foreach my $spec (sort keys %selected_rids) {
	# Edit species name
	my $spec_ = $spec;
	$spec_ =~ s/ /_/g;
	# Skip species other than determined by user
	if ($species) {
		unless ($spec_ eq $species) { next; }
	}
	# Count species
	$i_spec++;
	# Initialize number of pi runs for species
	my $spec_pi_runs = 0;
	# Get genus name
	my($genus) = ($spec =~ /^(\w+)\s/);
	# Check if species is supported by unitas
	my $unitas_spec = 'x';
	my @sp_match = ();
	if (@sp_match = grep {$_ =~ /^$spec_$/i} @unitas_specs) {
		$unitas_spec = $sp_match[0];
	} elsif (@sp_match = grep {$_ =~ /^$genus\_/i} @unitas_specs) {
		$unitas_spec = $sp_match[0];
	}
	# Add refseq option if no unitas species found
	my $refseq_opt = $unitas_spec eq 'x' ? '-refseq '.$refseq_file : '';
	# Get genome ftp links
	my $genome_link_gb = $gnm_list->{$spec}->[14] ? $gnm_list->{$spec}->[14] : '';
	my $genome_link_rs = $gnm_list->{$spec}->[15] ? $gnm_list->{$spec}->[15] : '';
	my $genome_link = $genome_link_rs ? $genome_link_rs : $genome_link_gb;
	# Skip species if no genome link is available
	unless ($genome_link) { next; }
	# Get cDNA seqs from ncbi if no unitas species found
	if ($unitas_spec eq 'x') {
		# Check if cds file is available in gb or rs
		my $cds_pres = 1;
		my $cds_link = '';
		my $cds_url_gb = '';
		my $cds_url_rs = '';
		if ($genome_link_gb) {
			my($cds_file_gb) = ($genome_link_gb =~ /\/([^\/]+)$/);
			$cds_file_gb .= '_cds_from_genomic.fna.gz';
			$cds_url_gb = $genome_link_gb.'/'.$cds_file_gb;
		}
		if ($genome_link_rs) {
			my($cds_file_rs) = ($genome_link_rs =~ /\/([^\/]+)$/);
			$cds_file_rs .= '_cds_from_genomic.fna.gz';
			$cds_url_rs = $genome_link_rs.'/'.$cds_file_rs;
		}
		# Get assembly name
		my($assembly_name) = ($genome_link =~ /\/\w+_[^\_]+_([^\/]+)$/);
		my $cds_file_end = "${assembly_name}_cds_from_genomic.fna.gz";
		my @cds_files = glob("UNITAS_refdump_x/*$cds_file_end*");
		unless (@cds_files && $cds_files[0]) {
			if (head($cds_url_rs)) {
				$cds_link = $genome_link_rs;
			} elsif (head($cds_url_gb)) {
				$cds_link = $genome_link_gb;
			} else {
				$cds_pres = 0;
			}
		}
		# Download cds file if it is available
		if ($cds_pres) {
			# Get name of cds file
			my($cds_file_gz) = ($cds_link =~ /\/([^\/]+)$/);
			$cds_file_gz .= '_cds_from_genomic.fna.gz';
			my($cds_file) = ($cds_file_gz =~ /(.*)\.gz$/);
			if (@cds_files && $cds_files[0]) {
				$cds_files[0] =~ s/UNITAS_refdump_x\///;
				$cds_file = $cds_files[0];
				$cds_file =~ s/\.gz//;
			}
			# Check if cds file exists, otherwise download and/or unzip
			unless (-e "UNITAS_refdump_x/$cds_file") {
				unless (-e "UNITAS_refdump_x/$cds_file_gz") {
					system("wget $cds_link/$cds_file_gz -P UNITAS_refdump_x >>ftp_log 2>&1");
				}
				system("gunzip -k UNITAS_refdump_x/$cds_file_gz >>sap_log 2>&1");
			}
			# If cds file downloaded successfully
			if (-e "UNITAS_refdump_x/$cds_file") {
				# Transform cds file to unitas file format
				cds_to_unitas_format("UNITAS_refdump_x/$cds_file") unless $no_rn;
				# Attach cds file to refseq option
				$refseq_opt .= ' -refseq UNITAS_refdump_x/'.$cds_file;
			}
		}
	}
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
			unless (-e "Genomes/$rep_file_gz") { next; }
		}
		system("gunzip -k Genomes/$rep_file_gz >>sap_log 2>&1");
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
		} else {
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
		my $chr_names = $genome_link_rs ? $chr_names_rs : $chr_names_gb;
		# Change chromosome names
		change_chr_gnm("Genomes/$gnm_file",$chr_names) unless $no_rn;
		change_chr_rep("Genomes/$rep_file",$chr_names) unless $no_rn;
		change_chr_gff("Genomes/$gff_file",$chr_names) unless $no_rn;
		change_chr_gff("Genomes/$gtf_file",$chr_names) unless $no_rn;
	}
	# Create bowtie genome index if it does not yet exist
	my($gnm_idx) = ($gnm_file =~ /(.*)\.fna$/);
	unless (-d "Genomes/idx/$gnm_idx") {
		system("bowtie-build Genomes/$gnm_file Genomes/idx/$gnm_idx/$gnm_idx");
	}
	# Create directory for species
	mkdir($data_dir.'/'.$spec_);
	# Initialize hit counts per mapped sequence for whole species
	my %spec_seq_hits = ();
	my $spec_seq_hits = \%spec_seq_hits;
	# Go through each tissue
	foreach my $tiss (sort keys %{$selected_rids{$spec}}) {
		# Create directory for tissue
		my $tisdir = $data_dir.'/'.$spec_.'/'.$tiss;
		mkdir($tisdir);
		# Do not include more than x runs for given tissue
		my $tis_runs = 0;
		my $tis_pi_runs = 0;
		# Save rids with pic prediction
		my @pt_rids = ();
		# Go through each run id
		foreach my $rid (@{$selected_rids{$spec}{$tiss}}) {
			# Count run ids for this tissue
			if ($tis_runs == $tis_runs_max) { last; }
			if ($tis_pi_runs == $tis_pi_runs_max) { last; }
			$tis_runs++;
			# Messaging
			my $message = '';
			# Name sra test file
			my $test_file = $rid.'_test.fasta';
			# Print information on runs
			my $n_test_seqs_mb = $n_test_seqs/1000000;
			$message .= "$spec\t$tiss\t$rid\t$n_test_seqs_mb/$sra_list->{$rid}->[$idxs->{'MBases'}]Mb";
			#### Download sra run test file unless it already exists
			unless (-e "$tisdir/$test_file") {
				system("$fq_dump $rid -X $n_test_seqs --fasta -O $tisdir >>sratool_log 2>&1");
				if (-e "$tisdir/$rid.fasta") {
					system("mv $tisdir/$rid.fasta $tisdir/$test_file >>sap_log 2>&1");
				} else {
					next;
				}
			}
			# Check if download was successful
			if (-e "$tisdir/$test_file") {
				$message .= " ok";
			} else {
				$message .= " fail";
			}
			#### Annotate sra test dataset
			# Adapter sequence
			my $adapter = '?';
			my $trimopt = '-trim '.$adapter.' -trim_minlength '.$min_len;
			# Check for possible adapter presence
			my $adapter_pres = check_adapter_presence("$tisdir/$test_file");
			# If no adapter present, no trimming needed
			if (!$adapter_pres) { $trimopt = '' }
			# Filter sequences by length
			if (!$adapter_pres) { filter_seq_length("$tisdir/$test_file",$min_len,$max_len) }
			# Check if unitas annotation already exists
			opendir(my $sra_dh1, $tisdir);
			unless (grep { "$tisdir/$_" =~ /_${test_file}_#/ } readdir $sra_dh1) {
				# Start unitas annotation
				my $infile = $tisdir.'/'.$test_file;
				system("$unitas -i $infile -s $unitas_spec -skip_piR $trimopt $refseq_opt >>unitas_log 2>&1");
				# Move output folder to sra file directory
				if (my @unidir = grep(-d, <*${test_file}_#*>)) {
					system("mv $unidir[0] $tisdir >>sap_log 2>&1");
				}
			}
			closedir($sra_dh1);
			#### Check whether no-annotation fraction represents pirnas
			my $pi_like = 0;
			# Initialize random adapter variables
			my($tag_bi_len,$tag_3p_len) = (0,0);
			# Get unitas output directory
			opendir(my $sra_dh2, $tisdir);
			my @unidir = grep { "$tisdir/$_" =~ /_${test_file}_#/ } readdir $sra_dh2;
			closedir($sra_dh2);
			# Check if unitas annotation exists	at tissue directory
			if (@unidir) {
				# Get relative unitas directory
				my $unidir = $tisdir.'/'.$unidir[0];
				# Check if no-annotation file exists
				my $noan_file = $unidir.'/fasta/unitas.no-annotation.fas';
				if (-s $noan_file) {
					# Get no-annotation sequences, being putative pirna sequences
					my $noan_seqs = get_fasta_seqs_collapsed($noan_file);
					### Get total counts of no-annotation reads/sequences
					my($noan_nsqs,$noan_nrds) = fasta_read_counts($noan_seqs);
					### Search for presence of random tags that need trimming
					if ($adapter_pres) {
						($tag_bi_len,$tag_3p_len) = check_random_tags($noan_seqs,$unidir,$test_file);
					}
					# Trim random tags if present
					if ($tag_bi_len == 4) {
						$message .= "\trm 4N-tags";
						$noan_seqs = trim_random_bi_tags($noan_seqs,$tag_bi_len);
					} elsif ($tag_3p_len == 3) {
						$message .= "\trm 3N-3tag";
						$noan_seqs = trim_random_3p_tags($noan_seqs,$tag_3p_len);
					} elsif ($adapter_pres) {
						$message .= "\ttrimmed ok";
					} else {
						$message .= "\tno adapter";
					}
					### Look for pirna trait 1: size fraction of 24-32 nt + sequence diversity
					my($pi_min,$pi_max) = (24,32);
					# Get length distribution
					my($nred_lens,$read_lens) = length_distribution($noan_seqs,16,36);
					# Get maximum length
					# Search for highest mir peak
					my @max_rlen = sort {$b <=> $a} keys %{$read_lens};
					# Get share of pirna size fraction
					my($pir_share,$pir_sdivt) = pi_size_fraction_share($nred_lens,$read_lens,$pi_min,$pi_max);
					### Look for pirna trait 2: 1T and 10A rates
					# Get positional nucleotide counts
					my $pi_pos_nucs = pos_nuc_counts($noan_seqs,$pi_min,$pi_max);
					# Get mean nucleotide shares
					my $mean_nucs = mean_pos_nuc_freqs($pi_pos_nucs,1,20);
					# Get 1T and 10A shares
					my($t1_rate,$a10_rate) = get_1T_10A_rates($pi_pos_nucs);
					### Decide whether no-annotation fraction contains pi-like rnas
					# Minimum pirna size share: 20%
					if ($pir_share >= 20) {
						# Minimum share of 1T or 10A: 35% (if mean 25%)
						if ($t1_rate >= ($mean_nucs->{'T'}*1.4) || $a10_rate >= ($mean_nucs->{'A'}*1.4)) {
							$pi_like = 1;
							$pi_like = 0 if $t1_rate < 10 || $a10_rate < 10;
						}
						# Alternative -> Minimum share of 1T and 10A: 30% (if mean 25%)
						if ($t1_rate >= ($mean_nucs->{'T'}*1.2) && $a10_rate >= ($mean_nucs->{'A'}*1.2)) {
							$pi_like = 1;
						}
						# Minimum total sequence diversity of pirna size fraction: 10%
						if ($pir_sdivt < 10) {
							$pi_like = 0;
						}
						# Cutoff for length maximum: 26nt
						if ($max_rlen[0] < 26) {
							$pi_like = 0;
						}
					}
					# Minimum no-annotation fraction share of all input reads: 0.1%
					if ($noan_nrds < $n_test_seqs*0.001) {
						$pi_like = 0;
					}
				} else {
					$message .= "\tbad data";
				}
			} else {
				$message .= "\tbad data";
			}
			#### Get complete sra dataset annotation and map to genome
			# Download whole sra dataset if pirnas present
			if ($pi_like) {
				$spec_pi_runs++;
				$tis_pi_runs++;
				# Info message
				$message .= "\tpiRs:Y";
				$message .= "\tpiCs:";
				# Download sra dataset
				unless (-e "$tisdir/$rid.fasta") {
					my $nseq_opt = $no_dl ? "-X $n_test_seqs" : "";
					system("$fq_dump $rid --fasta $nseq_opt -O $tisdir >>sratool_log 2>&1");
				}
				unless ($no_dl) {
					my $size = -s "$tisdir/$rid.fasta";
					if ($size/1000000 < 30) {
						system("$fq_dump $rid --fasta -O $tisdir >>sratool_log 2>&1");
					}
				}
				# Skip trimming and piC prediction if no no_pt option is true
				if ($no_pt) {
					$message .= "Skip pred.\n";
					print("$message");
				} else {
					# Name sra file
					my $sra_file = $rid.'.fasta';
					if ($adapter_pres) {
						# Trimmed sra data file name
						$sra_file = $sra_file.'.trimmed';
						# Start unitas trimming first only
						my $trim_file = "$tisdir/$sra_file";
						my @trimdir = ();
						unless (-e "$tisdir/$sra_file") {
							my $infile = $tisdir.'/'.$rid.'.fasta';
							system("$unitas -i $infile -s $unitas_spec $refseq_opt $trimopt -skip_mapping -skip_dust >>unitas_log 2>&1");
							# Get trimmed output file
							if (@trimdir = grep(-d, <*${rid}.fasta_#*>)) {
								$trim_file = $trimdir[0].'/'.$rid.'.fasta.trimmed.collapsed';
							}
						}
						# Extract sequences
						my $trim_seqs = get_fasta_seqs_collapsed($trim_file);
						# Remove 3' adapters and random tags if random tags are present
						if ($tag_bi_len == 4 || $tag_3p_len == 3) {
							# Remove random tags
							if ($tag_bi_len == 4) {
								$trim_seqs = trim_random_bi_tags($trim_seqs,$tag_bi_len);
							} elsif ($tag_3p_len == 3) {
								$trim_seqs = trim_random_3p_tags($trim_seqs,$tag_3p_len);
							}
						}
						# Print de-tagged sequences
						unless (-e "$tisdir/$sra_file") {
							my $out = open_outfile("$tisdir/$sra_file");
							foreach my $seq (sort {$trim_seqs->{$b} <=> $trim_seqs->{$a}} keys %{$trim_seqs}) {
								print($out ">$trim_seqs->{$seq}\n$seq\n");
							}
						}
						# Delete unitas output
						system("rm -rf $trimdir[0] >>sap_log 2>&1") if @trimdir;
					} else {
						# Collapse sra reads
						my $trim_seqs = collapse_fasta_seqs("$tisdir/$sra_file");
						# Trimmed sra data file name
						$sra_file = $sra_file.'.trimmed';
						# Print de-tagged sequences
						unless (-e "$tisdir/$sra_file") {
							my $out = open_outfile("$tisdir/$sra_file");
							foreach my $seq (sort {$trim_seqs->{$b} <=> $trim_seqs->{$a}} keys %{$trim_seqs}) {
								print($out ">$trim_seqs->{$seq}\n$seq\n");
							}
						}
					}
					# Filter sequences by length
					filter_seq_length("$tisdir/$sra_file",$min_len,$max_len);
					## Map trimmed sra file to genome
					my $map_file = $sra_file.'.sam';
					unless (-e "$tisdir/$map_file") {
						unless (-e "$tisdir/$map_file.gz") {
							my $sqm_opts = '/output_all_matches';
							system("bowtie Genomes/idx/$gnm_idx/$gnm_idx $tisdir/$sra_file -f -S $tisdir/$map_file >>bwt_log 2>&1");
						}
					}
					# Get unitas output directory
					opendir(my $sra_dh3, $tisdir);
					# Check if unitas annotation already exists
					unless (grep { "$tisdir/$_" =~ /_${map_file}_#/ } readdir $sra_dh3) {
						# Set up files for unitas annotation
						unless (-e "$tisdir/$map_file") {
							system("gunzip $tisdir/$map_file.gz -k");
						}
						my $infile = $tisdir.'/'.$map_file;
						# Start unitas annotation
						system("$unitas -i $infile -s $unitas_spec -skip_piR -skip_dust -pp $refseq_opt >>unitas_log 2>&1");
						# Move output folder to sra file directory
						if (my @unidir = grep(-d, <*${map_file}_#*>)) {
							system("mv $unidir[0] $tisdir >>sap_log 2>&1");
						}
					}
					closedir($sra_dh3);
					# Get unitas output directory
					opendir(my $sra_dh4, $tisdir);
					my @sra_unidir = grep { "$tisdir/$_" =~ /_${map_file}_#/ } readdir $sra_dh4;
					closedir($sra_dh4);
					# Initialize number of predicted pics
					my $n_pics = 0;
					# No-annotation file
					my $noan_file = $tisdir.'/'.$sra_unidir[0].'/fasta/unitas.non-miR.fas';
					### If no-annotation file exists, clean map file and run proTRAC
					if (-s $noan_file) {
						# Clean up map file
						my $clean_map_file = $map_file;
						#$clean_map_file =~ s/\.sam/\.noa\.sam/;
						unless (-e "$tisdir/$clean_map_file") {
							unless (-e "$tisdir/$map_file") {
								system("gunzip $tisdir/$map_file.gz -k");
							}
							#filter_mapped_seqs("$tisdir/$map_file",$noan_file,"$tisdir/$clean_map_file");
							system("gzip $tisdir/$map_file") unless "-e $tisdir/$map_file.gz";
						}
						# Predict pirna clusters
						my $map_d = "$tisdir/$clean_map_file";
						my $gnm_d = "Genomes/$gnm_file";
						my $rep_d = "Genomes/$rep_file";
						my $gen_d = -s "Genomes/$gtf_file" ? "Genomes/$gtf_file" : "Genomes/$gff_file";
						my $p_opt = "-1Tor10A 0.5 -1Tand10A 0.3 -clsize 5000";
						# Check if pic prediction already exists
						opendir(my $pro_dh, $tisdir);
						unless (grep { "$tisdir/$_" =~ /proTRAC_$clean_map_file/ } readdir $pro_dh) {
							close($pro_dh);
							system("$protrac -m $map_d -g $gnm_d -repeatmasker $rep_d -geneset $gen_d $p_opt >>pT_log 2>&1");
							# Move output folder to sra file directory
							if (my @pt_dir = grep(-d, <*proTRAC_$clean_map_file*>)) {
								system("mv $pt_dir[0] $tisdir >>sap_log 2>&1");
							}
						}
						opendir($pro_dh, $tisdir);
						# Get number of piCs
						if (my @new_pt_dir = grep { "$tisdir/$_" =~ /proTRAC_$clean_map_file/ } readdir $pro_dh) {
							close($pro_dh);
							opendir(my $ndh, "$tisdir/$new_pt_dir[0]");
							my @pic_files = grep { "$tisdir/$new_pt_dir[0]/$_" =~ /\.html/ } readdir $ndh;
							$n_pics = @pic_files;
							if ($n_pics > 0) { push(@pt_rids,$rid)}
						}
						# Get hit counts for each mapped sequence
						my $seq_hits = get_seq_hits($map_d);
						$spec_seq_hits = merge_seq_hits($spec_seq_hits,$seq_hits);
					}
					# Info message
					$message .= "$n_pics\n";
					print("$message");
				}
			} else {
				# Info message
				$message .= "\tpiRs:N";
				$message .= "\tNo piCs\n";
				print("$message");
			}
		}
		## Find robust/consens pics for each tissue if pirnas found
		if ($tis_pi_runs) {
			# Get robust pic locs for data set of tissue
			my $i = 0;
			my @pic_locs = ();
			# Go through each sra id with proTRAC annotation
			foreach my $rid (@pt_rids) {
				opendir(my $pro_dh, $tisdir);
				if (my @pt_dir = grep { "$tisdir/$_" =~ /proTRAC_$rid\.fasta/ } readdir $pro_dh) {
					close($pro_dh);
					# Get clusters gtf file
					my $pics_file = $tisdir.'/'.$pt_dir[0].'/clusters.gtf';
					# Extract gtf data
					$pic_locs[$i] = get_gtf_fields($pics_file);
					$i++;
				}
			}
			# Get all pic locs of all sra runs and count identicals
			my ($all_locs,$all_loc_counts) = get_all_locs(\@pic_locs);
			# Output robust clusters
			my $out = open_outfile($tisdir.'/'.$tiss.'.robust_clusters.gtf');
			foreach my $id (sort {$a <=> $b} keys %{$all_locs}) {
				my $chr = $all_locs->{$id}->[0];
				my $beg = $all_locs->{$id}->[3];
				my $end = $all_locs->{$id}->[4];
				if ($all_loc_counts->{$id}/(scalar @pic_locs) >= 0.5) {
					foreach my $field (@{$all_locs->{$id}}) {
						print($out "$field\t");
					}
					print($out "\n");
				}
			}
		}
	}
	## Find robust/consens pics for each species if pirnas found
	if ($spec_pi_runs) {
		# Get robust pic locs for each tissue
		my $i = 0;
		my @pic_locs = ();
		# Go through each tissue
		foreach my $tiss (sort keys %{$selected_rids{$spec}}) {
			# Get directory for tissue and tissue gtf name
			my $tisdir = $data_dir.'/'.$spec_.'/'.$tiss;
			my $tisgtf = $tisdir.'/'.$tiss.'.robust_clusters.gtf';
			# Check if gtf exists
			if (-s $tisgtf) {
				# Extract gtf data
				$pic_locs[$i] = get_gtf_fields($tisgtf);
				$i++;
			}
		}
		# Get all pic locs of all sra runs and count identicals
		my ($all_locs,$all_loc_counts) = get_all_locs(\@pic_locs);
		# Output robust clusters
		my $out = open_outfile($data_dir.'/'.$spec_.'/'.$spec_.'.robust_clusters.gtf');
		foreach my $id (sort {$a <=> $b} keys %{$all_locs}) {
			my $chr = $all_locs->{$id}->[0];
			my $beg = $all_locs->{$id}->[3];
			my $end = $all_locs->{$id}->[4];
			if ($all_loc_counts->{$id}/(scalar @pic_locs) >= 0.5) {
				foreach my $field (@{$all_locs->{$id}}) {
					print($out "$field\t");
				}
				print($out "\n");
			}
		}
	}
	# Output genomic hits for all mapped sequences
	my $hits_out = open_outfile($data_dir.'/'.$spec_.'/'.$spec_.'.hit_counts.txt');
	foreach my $seq (sort {$spec_seq_hits->{$b} <=> $spec_seq_hits->{$a}} keys %{$spec_seq_hits}) {
		print($hits_out "$seq\t$spec_seq_hits->{$seq}\n");
	}
}

exit;

################################# subroutines #################################

sub merge_seq_hits {
	# Take hashes with all and current hit counts
	my($all_hits,$seq_hits) = @_;
	# Go through each sequence of current hash
	foreach my $seq (keys %{$seq_hits}) {
		# If not present in all hits hash
		unless ($all_hits->{$seq}) {
			# Add to all hits hash
			$all_hits->{$seq} = $seq_hits->{$seq};
		}
	}
	return $all_hits;
}

sub get_seq_hits {
	# Take map file name
	my($file) = @_;
	# Variables
	my %seq_hits = ();
	# Open file
	my $in = open_infile($file);
	# Extract sequence and save with hit count
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		if ($line !~ /^@/) {
			my @d = split(/\t/,$line);
			my $seq = $d[9];
			$seq_hits{$seq}++;
		}
	}
	return \%seq_hits;
}

sub get_all_locs {
	# Take pic locs of each data set
	my($pic_locs) = @_;
	# Storage variables
	my %all_locs = ();
	my %all_loc_counts = ();
	# Go through each dataset
	foreach my $x (0..((scalar @{$pic_locs})-1)) {
		# Go through each pic loc
		foreach my $pic_x (sort {$a <=> $b} keys %{$pic_locs->[$x]}) {
			# Get loc coordinates
			my $chr_x = $pic_locs->[$x]->{$pic_x}->[0];
			my $beg_x = $pic_locs->[$x]->{$pic_x}->[3];
			my $end_x = $pic_locs->[$x]->{$pic_x}->[4];
			# Check if this pic is already in all pics list
			my $pic_in_list = 0;
			my $id_y = 0;
			# Go through each pic loc
			foreach my $pic_y (sort {$a <=> $b} keys %all_locs) {
				# Get loc coordinates
				my $chr_y = $all_locs{$pic_y}->[0];
				my $beg_y = $all_locs{$pic_y}->[3];
				my $end_y = $all_locs{$pic_y}->[4];
				# Check if locs overlap
				if ($chr_x eq $chr_y && $beg_y <= $end_x && $end_y >= $beg_x) {
					$pic_in_list = 1;
					$id_y = $pic_y;
				}
			}
			if ($pic_in_list) {
				# Get loc coordinates
				my $chr_y = $all_locs{$id_y}->[0];
				my $beg_y = $all_locs{$id_y}->[3];
				my $end_y = $all_locs{$id_y}->[4];
				# Get limit coordinates
				my $beg_n = $beg_x < $beg_y ? $beg_x : $beg_y;
				my $end_n = $end_x > $end_y ? $end_x : $end_y;
				# Update coordinates
				$all_locs{$id_y}->[3] = $beg_n;
				$all_locs{$id_y}->[4] = $end_n;
				# Increment count for loc
				$all_loc_counts{$id_y}++;
			} else {
				# Get new id
				my $id_n = 1;
				foreach my $id (sort {$a <=> $b} keys %all_locs) { $id_n = $id+1 }
				# Save new loc
				@{$all_locs{$id_n}} = @{$pic_locs->[$x]->{$pic_x}};
				# Initialize count for loc
				$all_loc_counts{$id_n} = 1;
			}
		}
	}
	return \(%all_locs,%all_loc_counts);
}

sub filter_mapped_seqs {
	# Take fasta file name
	my($map_file,$fas_file,$new_file) = @_;
	# Open in and out files
	my $map = open_infile($map_file);
	my $out = open_outfile($new_file);
	# Get fasta sequences
	my $fas_seqs = get_fasta_seqs_collapsed($fas_file);
	# Go through each line
	while (my $line = <$map>) {
		$line =~ s/\s+$//; #better chomp
		if ($line !~ /@/) {
			# Get line data
	        my @d = split(/\t/,$line);
			# Get probe sequence
			my $seq = $d[9];
			# Print to file if sequence is part of fas seqs
			if ($fas_seqs->{$seq}) {
				print($out "$line\n");
			}
		} else {
			print($out "$line\n");
		}
	}
}

sub filter_seq_length {
	# Take fasta file name
	my($infile,$min,$max) = @_;
	# Storage variable
	my $header = '';
	my $seq = '';
	# Open in and out files
	my $in = open_infile($infile);
	my $outfile = $infile.'.filter';
	my $out = open_outfile($outfile);
	# Go through each line
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Header line
		if ($line =~ /^>/) {
			if ($seq) {
				my $len = length($seq);
				if ($len >= $min && $len <= $max) {
					print($out "$header\n$seq\n");
				}
				$seq = '';
			}
			$header = $line;
		} else {
			$seq .= $line;
		}
	}
	if ($seq) {
		my $len = length($seq);
		if ($len >= $min && $len <= $max) {
			print($out "$header\n$seq\n");
			$seq = '';
		}
	}
	# Close files
	close($in);
	close($out);
	# Delete old file and replace with cleaned file
	unlink($infile);
	system("mv $outfile $infile >>sap_log 2>&1");
}

sub check_adapter_presence {
	# Take fasta file name
	my($file) = @_;
	# Storage variable
	my %lens = ();
	my $n_seqs = 0;
	my $seq = '';
	# Open in and out files
	my $fas = open_infile($file);
	# Clean fasta headers
	while (my $line = <$fas>) {
		$line =~ s/\s+$//; #better chomp
		# Seq line
		if ($line =~ /^>/) {
			if ($seq) {
				$lens{length($seq)}++;
				$n_seqs++;
				$seq = '';
			}
		} else {
			$seq .= $line;
		}
	}
	# Check if all seq lengths are equal
	my $equal_lens = 0;
	foreach my $len (sort keys %lens) {
		if ($lens{$len} eq $n_seqs) {
			$equal_lens = 1;
		}
	}
	my $adapter_presence = $equal_lens ? 1 : 0;
	# Check if unequal lengths are in srna size range
	unless ($equal_lens) {
		my $srna_len_seqs = 0;
		foreach my $len (16..36) {
			$srna_len_seqs += $lens{$len} if $lens{$len};
		}
		if ($srna_len_seqs/$n_seqs < 0.25) {
			$adapter_presence = 1;
		}
	}
	return $adapter_presence;
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

sub change_chr_gff {
	# Take name of genome file
	my($infile,$chr_names) = @_;
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

sub cds_to_unitas_format {
	# Take fasta file name
	my($file) = @_;
	my $newf = $file.'.new';
	# Open in and out files
	my $fas = open_infile($file);
	my $new = open_outfile($newf);
	# mRNA bool
	my $gid = '';
	# Clean fasta headers
	while (my $line = <$fas>) {
		$line =~ s/\s+$//; #better chomp
		if ($line =~ /^>/) {
			($gid) = ($line =~ /^lcl\|(\S+)/);
			unless ($gid) { $gid = 'unknown' }
			print($new ">protein_coding|$gid\n");
		} else {
			print($new "$line\n");
		}
	}
	# Close files
	close($fas);
	close($new);
	# Delete old file and replace with cleaned file
	unlink($file);
	system("mv $newf $file >>sap_log 2>&1");
}

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

sub check_random_tags {
	# Take no-annotation seqs, unitas directory and test file name
	my($noan_seqs,$unidir,$test_file) = @_;
	# Random tag lengths
	my $tag_bi_len = 0;
	my $tag_3p_len = 0;
	# Get pre-annotation sequences
	my $totl_file = $unidir.'/'.$test_file.'.trimmed.collapsed';
	my $totl_seqs = get_fasta_seqs_collapsed($totl_file);
	# Get no-annotation share
	my($noan_nsqs,$noan_nrds) = fasta_read_counts($noan_seqs);
	my($totl_nsqs,$totl_nrds) = fasta_read_counts($totl_seqs);
	my $noan_psqs = $totl_nsqs ? int(($noan_nsqs/$totl_nsqs*100)+0.5) : 0;
	# Share of no-annotation sequences > 80% may indicate presence of random tags
	if ($noan_psqs > 80) {
		# Get positional nucleotide counts
		my $all_pos_nucs = pos_nuc_counts($noan_seqs);
		# Get putative shifted t rate peak
		my %pos_t_rates = ();
		foreach my $pos (1,5) {
			# Get t rate
			my $pos_count = 0;
			for (values %{$all_pos_nucs->{$pos}}) { $pos_count += $_ }
			my $t_count = $all_pos_nucs->{$pos}->{'T'} ? $all_pos_nucs->{$pos}->{'T'} : 0;
			my $t_share  = $pos_count ? int(($t_count/$pos_count*100)+0.5) : 0;
			# Save t rate
			$pos_t_rates{$pos} = $t_share;
		}
		# Search for highest t rate
		my @max_t_pos = sort {$pos_t_rates{$b} <=> $pos_t_rates{$a}} keys %pos_t_rates;
		my $max_t_pos = $max_t_pos[0];
		# Check if highest t rate is t peak
		if ($pos_t_rates{$max_t_pos} >= 30) {
			# Get presumed tag length
			$tag_bi_len = $max_t_pos-1;
		}
		# Get mirna peak if maybe xN 3' tag
		if ($max_t_pos == 1) {
			# Extract mir fasta files
			my $fasdir = $unidir.'/fasta';
			opendir(my $fas_dh, $fasdir);
			my @mirfiles = grep { "$fasdir/$_" =~ /\.miR\./ } readdir $fas_dh;
			closedir($fas_dh);
			unless (@mirfiles) { @mirfiles = ('unitas.miRNA.fas') } # -s x
			if (@mirfiles) {
				# Get mir sequences
				my $mirfile = '';
				foreach my $file (@mirfiles) {
					if ($file =~ /\.non-/) { next }
					$mirfile = $fasdir.'/'.$file;
				}
				if (-e $mirfile) {
					my $mir_seqs = get_fasta_seqs_collapsed($mirfile);
					my($nred_lens,$read_lens) = length_distribution($mir_seqs,22,28);
					# Search for highest mir peak
					my @max_pos = sort {$read_lens->{$b} <=> $read_lens->{$a}} keys %{$read_lens};
					my $max_pos = $max_pos[0];
					$tag_3p_len = $max_pos-22 if $max_pos;
				}
			}
		} else {
			$tag_3p_len = 0;
		}
	}
	return ($tag_bi_len,$tag_3p_len);
}

sub fasta_read_counts {
	# Take fasta seqs
	my($fasta_seqs) = @_;
	# Storage variables
	my $n_seqs = 0;
	my $n_reds = 0;
	# Go through each sequence
	foreach my $seq (keys %{$fasta_seqs}) {
		# Count sequences
		$n_seqs++;
		$n_reds += $fasta_seqs->{$seq};
	}
	return ($n_seqs,$n_reds);
}

sub pi_size_fraction_share {
	# Take length distribution and pirna size limits
	my($nred_lens,$read_lens,$pi_min,$pi_max) = @_;
	# Set pirna size range
	my @maxlen = sort {$b <=> $a} keys %{$nred_lens};
	$pi_max = $maxlen[0] < $pi_max ? $maxlen[0] : $pi_max if $maxlen[0];
	# Initialize variables
	my $pir_nreds = 0;
	my $all_nreds = 0;
	my $pir_reads = 0;
	my $all_reads = 0;
	my $pir_share = 0;
	my $pir_sdivt = 0;
	foreach my $len (sort keys %{$nred_lens}) {
		$nred_lens->{$len} = 0 unless $nred_lens->{$len};
		$read_lens->{$len} = 0 unless $read_lens->{$len};
		$all_nreds += $nred_lens->{$len};
		$pir_nreds += $nred_lens->{$len} if grep {$len == $_} ($pi_min..$pi_max);
		$all_reads += $read_lens->{$len};
		$pir_reads += $read_lens->{$len} if grep {$len == $_} ($pi_min..$pi_max);
	}
	$pir_share = int(($pir_nreds/$all_nreds*100)+0.5) if $all_nreds;
	$pir_sdivt = int(($pir_nreds/$pir_reads*100)+0.5) if $pir_reads;
	return ($pir_share,$pir_sdivt);
}

sub get_1T_10A_rates {
	# Take positional nucleotide counts
	my($pos_nucs) = @_;
	# Initialize variables
	my $pos1_count = 0;
	my $pos10_count = 0;
	for (values %{$pos_nucs->{1}}) { $pos1_count += $_ }
	for (values %{$pos_nucs->{10}}) { $pos10_count += $_ }
	my $t1_count = $pos_nucs->{1}->{'T'} ? $pos_nucs->{1}->{'T'} : 0;
	my $a10_count = $pos_nucs->{10}->{'A'} ? $pos_nucs->{10}->{'A'} : 0;
	my $t1_rate = $pos1_count ? int($t1_count/$pos1_count*100) : 0;
	my $a10_rate = $pos10_count ? int($a10_count/$pos10_count*100) : 0;
	# Return 1T and 10A rates
	return ($t1_rate,$a10_rate);
}

sub mean_pos_nuc_freqs {
	# Take positional nucleotide counts
	my($pos_nucs,$min,$max) = @_;
	# Get positional shares for each nucleotide
	my %pos_nuc_shares = ();
	foreach my $pos ($min..$max) {
		# Get total positions count
		my $pos_count = 0;
		for (values %{$pos_nucs->{$pos}}) { $pos_count += $_ }
		# Get counts for each nucleotide
		foreach my $nuc (sort keys %{$pos_nucs->{$pos}}) {
			my $nuc_count = $pos_nucs->{$pos}->{$nuc} ? $pos_nucs->{$pos}->{$nuc} : 0;
			my $nuc_share = $pos_count ? int($nuc_count/$pos_count*100) : 0;
			# Save nucleotide share
			push(@{$pos_nuc_shares{$nuc}},$nuc_share);
		}
	}
	# Calculate mean shares for each nucleotide
	my %mean_nuc_shares = ();
	foreach my $nuc (sort keys %pos_nuc_shares) {
		$mean_nuc_shares{$nuc} = get_mean($pos_nuc_shares{$nuc});
	}
	return \%mean_nuc_shares;
}

sub trim_random_3p_tags {
	# Take fasta seqs
	my($fasta_seqs,$tag_3p_len) = @_;
	# Storage variables
	my %trimmed_seqs = ();
	# Go through each sequence
	foreach my $seq (keys %{$fasta_seqs}) {
		# Get read length
		my $seq_len = length($seq);
		# Trim sequence
		my $trim_seq = substr($seq,0,$seq_len-$tag_3p_len);
		# Save trimmed sequence
		$trimmed_seqs{$trim_seq} += $fasta_seqs->{$seq};
	}
	return \%trimmed_seqs;
}

sub trim_random_bi_tags {
	# Take fasta seqs
	my($fasta_seqs,$tag_bi_len) = @_;
	# Storage variables
	my %trimmed_seqs = ();
	# Go through each sequence
	foreach my $seq (keys %{$fasta_seqs}) {
		# Get read length
		my $seq_len = length($seq);
		# Trim sequence
		my $n_tags = 2;
		my $trim_seq = substr($seq,$tag_bi_len,$seq_len-($tag_bi_len*$n_tags));
		# Save trimmed sequence
		$trimmed_seqs{$trim_seq} += $fasta_seqs->{$seq};
	}
	return \%trimmed_seqs;
}

sub unitas_noannotation_share {
	# Take summary file name
	my($file) = @_;
	# Get annotation summary counts
	my $anno_counts = unitas_summary_counts($file);
	my $noann_count = $anno_counts->{'no annotation'} ? $anno_counts->{'no annotation'} : 0;
	my $total_count = 0;
	for (values %{$anno_counts}) { $total_count += $_ }
	my $noann_share = $noann_count/$total_count*100;
}

sub unitas_noannotation_count {
	# Take summary file name
	my($file) = @_;
	# Get annotation summary counts
	my $anno_counts = unitas_summary_counts($file);
	my $noann_count = $anno_counts->{'no annotation'} ? $anno_counts->{'no annotation'} : 0;
}

sub unitas_summary_counts {
	# Take summary file name
	my($file) = @_;
	# Storage variable
	my %anno_counts = ();
	# Get file data
	my @sumr_data = get_file_data_array($file);
	# Parse file and extract read counts
	foreach my $line (@sumr_data) {
		if ($line =~ /^\w+/) {
			my @d = split(/\t/,$line);
			$anno_counts{$d[0]} = $d[1];
		}
	}
	return \%anno_counts;
}

sub pos_nuc_counts {
	# Take fasta seqs
	my($fasta_seqs,$min,$max) = @_;
	# Storage variables
	my %pos_nucs = ();
	# Go through each sequence
	foreach my $seq (keys %{$fasta_seqs}) {
		# Get read length
		my $len = length($seq);
		# Apply cutoff if provided
		if ($min && $max) {
			if ($len < $min || $len > $max) {
				next;
			}
		}
		# Slide through sequence
		for (my $pos = 0; $pos < $len; $pos++) {
			# Get nucleotide
			my $nuc = substr($seq,$pos,1);
			$pos_nucs{$pos+1}{$nuc} += 1;
		}
	}
	return \%pos_nucs;
}

sub length_distribution {
	# Take fasta seqs
	my($fasta_seqs,$min,$max) = @_;
	# Storage variables
	my %read_len_dist = ();
	my %nred_len_dist = ();
	my %div_rate = ();
	# Go through each sequence
	foreach my $seq (keys %{$fasta_seqs}) {
		# Get read length
		my $len = length($seq);
		# Apply cutoff if provided
		if ($min && $max) {
			if ($len < $min || $len > $max) {
				next;
			}
		}
		# Count reads
		$read_len_dist{$len} += $fasta_seqs->{$seq};
		$nred_len_dist{$len} += 1;
	}
	# Get diversity rate
	my $all_nred = 0;
	my $all_read = 0;
	foreach my $len (keys %read_len_dist) {
		$div_rate{$len} = $nred_len_dist{$len}/$read_len_dist{$len};
		$all_nred += $nred_len_dist{$len};
		$all_read += $read_len_dist{$len};
	}
	my $div_totl = $all_read ? $all_nred/$all_read*100 : 0;
	return \(%nred_len_dist,%read_len_dist,%div_rate,$div_totl);
}

sub get_field_indexes {
	# Take name of pic tree file
	my($sra_run_table) = @_;
	# Initialize variable
	my %names_indexes = ();
	# Get file data
	my @srarun_data = get_file_data_array($sra_run_table);
	my @field_names = split(/\t/,$srarun_data[0]);
	# Loop through field names
	foreach my $i (0..$#field_names) {
	    $names_indexes{$field_names[$i]} = $i;
	}
	return \%names_indexes;
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

# Calculate mean of array values
sub get_mean {
	# Take array list
	my($array) = @_;
	# Number of values
	my $N = scalar(@{$array});
	if ($N == 0) { return 0 }
	# Sum values
	my $sum = get_sum($array);
	# Calculate mean value
	my $mean = $sum/$N;
	# Return mean
	return $mean;
}

sub get_fasta_seqs_collapsed {
	# Take fasta file name
	my($file,$short) = @_;
	# Variables
	my $reads = '';
	my %sequences = ();
	# Open file
	my $in = open_infile($file);
	# Extract sequence and save with name
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		if ($line =~ /^>/) {
			($reads) = ($line =~ />(.*)$/);
			if ($short) { $reads =~ s/\s.*// }
		} else {
			$sequences{$line} += $reads;
		}
	}
	return \%sequences;
}

sub collapse_fasta_seqs {
	# Take fasta file name
	my($file) = @_;
	# Variables
	my %sequences = ();
	# Open file
	my $in = open_infile($file);
	# Extract sequence and save with name
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		if ($line !~ /^>/) {
			$sequences{$line}++;
		}
	}
	return \%sequences;
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
