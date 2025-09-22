#! /usr/bin/env perl
#结果是未比对到rDNA的read
# By Nathan Sheffield, University of Virginia, 2018

# This is an incredibly fast Perl utility that re-pairs fastq files that have been 
# de-paired by running a single-end alignment on paired-end data. 
# It assumes that the the first fastq file (say, read1) contains a subset of reads
# found in the second fastq file (say, read2). 
# It also assumes that the files were in proper order to begin with.
# It will return a filtered version of the second file, which only has reads
# present in the first tile.

# Perl is the right language for this utility; because this is a high
# IO task (spitting out hundreds of millions of lines), and Perl is highly
# optimized for this type of IO process. A corresponding python
# program will take 10-100 fold longer to do the same thing.

# Setup
$file_filter = shift;
#outfolder+"prealignments"+${id}+'_'+genome+'_unmap.fq'(没有比对到rDNA上的read)
$file_fq1 = shift;
#${id}_raw_1_processed_noUMI.fastq(*_R1_processed_noUMI.fastq)
$file_fq2 = shift;
#${id}_raw_2_trimmed_noUMI.fastq(*_R2_trimmed_noUMI.fastq)
$file_fq1_filtered = shift;
#outfolder+"prealignments"+${id}+'_'+genome+'_unmap_R1.fq'
$file_fq2_filtered = shift;
#outfolder+"prealignments"+${id}+'_'+genome+'_unmap_R2.fq'
open(my $fh_filter, "<", $file_filter);
#outfolder+"prealignments"+${id}+'_'+genome+'_unmap.fq'
# We can't read from compressed input because it messes with the way the files
# are read in, but this is how you would do it.
# open(my $fh_fq1, "gunzip -c $file_fq1 |");
# open(my $fh_fq2, "gunzip -c $file_fq2 |");
open(my $fh_fq1, "<", $file_fq1);
#${id}_raw_1_processed_noUMI.fastq(*_R1_processed_noUMI.fastq)
open(my $fh_fq2, "<", $file_fq2);
#${id}_raw_2_trimmed_noUMI.fastq(*_R2_trimmed_noUMI.fastq)


# write output files here
if ($file_fq1_filtered =~ /\.gz$/i) {
	print STDERR "gzipping output\n";
	open(FH_FQ1_FILT, "| /bin/gzip -c > $file_fq1_filtered");
} else {
	print STDERR "not gzipping output\n";
	open(FH_FQ1_FILT, ">", $file_fq1_filtered);
}

if ($file_fq2_filtered =~ /\.gz$/i) {
	open(FH_FQ2_FILT, "| /bin/gzip -c > $file_fq2_filtered");
} else {
	open(FH_FQ2_FILT, ">", $file_fq2_filtered);
}



# Loop through reads
my $skipped = 0;
my %bhash; 

# load some read names into buffer
for ($r = 1; $r < 1000000; $r++) {
	$readnamef = <$fh_filter> or last;
	# 从文件outfolder+"prealignments"+${id}+'_'+genome+'_unmap.fq'读取一行,若读取失败(如文件结束)则退出循环.
	$readnamef =~ s/[\s\/].*$//;
	#提取outfolder+"prealignments"+${id}+'_'+genome+'_unmap.fq'的每行的首个字段,即序列名称
	<$fh_filter>;<$fh_filter>;<$fh_filter>;
	# 跳过后续3行
	chomp($readnamef);
	$bhash{$readnamef} = 1;
	# print $readname;
}

# print "$_\n" for keys %bhash;
#打印每个键‌
# print "BLAH\n\n";
while($readname2 = <$fh_fq2>) {
	#${id}_raw_2_trimmed_noUMI.fastq(*_R2_trimmed_noUMI.fastq)
	$readname1 = <$fh_fq1>;
	#${id}_raw_1_processed_noUMI.fastq(*_R1_processed_noUMI.fastq)
	#print STDERR "readname1:\t$readname1";
	$readname2_copy = $readname2;
	$readname2 =~ s/[\s\/].*$//;
	#print STDERR "readname2:\t$readname2";
	# if ($skipped < 50) { print STDERR ($readname2)};
	chomp($readname2);
	# if ("$readname" eq "$readname2") {
	if (exists $bhash{$readname2}) {
 		print FH_FQ2_FILT $readname2_copy;
		#向outfolder+"prealignments"+${id}+'_'+genome+'_unmap_R2.fq'输出
 		print FH_FQ2_FILT $line = <$fh_fq2>;
 		print FH_FQ2_FILT $line = <$fh_fq2>;
 		print FH_FQ2_FILT $line = <$fh_fq2>;
 		print FH_FQ1_FILT $readname1;
		#向outfolder+"prealignments"+${id}+'_'+genome+'_unmap_R1.fq'输出
 		print FH_FQ1_FILT $line = <$fh_fq1>;
 		print FH_FQ1_FILT $line = <$fh_fq1>;
 		print FH_FQ1_FILT $line = <$fh_fq1>;
 		
 		delete $bhash{$readname2};
 		# Parse in a new read from the filter:
 		if ($readnamef = <$fh_filter>) {
			#outfolder+"prealignments"+${id}+'_'+genome+'_unmap.fq'
			$readnamef =~ s/[\s\/].*$//;
			<$fh_filter>;<$fh_filter>;<$fh_filter>;
			chomp($readnamef);
			$bhash{$readnamef} = 1;
		}
	} else {
		# advance to next r2 read
 		$skipped++;
		<$fh_fq2>.<$fh_fq2>.<$fh_fq2>;
		<$fh_fq1>.<$fh_fq1>.<$fh_fq1>;
	}
}

close(FH_FQ1_FILT);
#关闭outfolder+"prealignments"+${id}+'_'+genome+'_unmap_R1.fq'
close(FH_FQ2_FILT);
#关闭outfolder+"prealignments"+${id}+'_'+genome+'_unmap_R2.fq'

$lost_reads_count = keys %bhash;
print STDERR $skipped." reads skipped\n";
print STDERR $lost_reads_count." reads lost\n";
#即没有比对上的reads的数目
if ($lost_reads_count < 200) {
	print STDERR "$_\n" for keys %bhash;
}
