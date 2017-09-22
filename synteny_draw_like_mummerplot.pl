#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::FastaKit qw/IndexFasta/;
use SVG;
use Data::Dumper qw/Dumper/;
use constant USAGE =><<EOH;

usage: $0 last.output fasta[.fai] output.svg

    Plot synteny to SVG in MUMmerplot
    
Requirements:
	perl Modules: SVG, Data::Dumper, FuhaoPerl5Lib::FastaKit
	samtools

v20170110

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



### Input and Output
die "Error: invalid input synteny file\n" unless (defined $ARGV[0] and -s $ARGV[0]);
my $inputconfig=$ARGV[0];
die "Error: invalid input fasta.fai file\n" unless (defined $ARGV[1] and -s $ARGV[1]);
my $inputfasta=$ARGV[1];
die "Error: invalid output file\n" unless (defined $ARGV[2]);
die "Error: output should have a extension name: svg/SVG\n" unless ($ARGV[2]=~/\.svg$/i);
unlink $ARGV[2] if (-e $ARGV[2]);
my $outputsvg=$ARGV[2];



### Default
my %confighash=();
my %seqlength=();
my $linenum=0;


### Configure
%confighash=(
	'plot_width' => 1000, #Default: 200
	'plot_height' => 1000, #Default: 200
	'plot_left_margin' => 10, #Default: 10
	'plot_right_margin' => 10, #Default: 10
	'plot_top_margin' => 10, #Default: 10
	'plot_bottom_margin' => 10, #Default: 10
	
	'x_line_height' => 2, #Default: 2
	'x_line_color' => 'black', #Default: 'black'
	'x_tick_height' => 5, #Default: 5
	'x_mark_font' => 'Arial', #Default: 'Arial'
	'x_mark_size' => 10, ### 0 = no marks #Default: 5
	'x_legend_font' => 'Arial', #Default: 'Arial'
	'x_legend_size' => 10, ### 0 = no marks #Default: 5
	'x_legend_color' => 'black', #Default: 'black'
	
	'y_line_width' => 2, #Default: 2
	'y_line_color' => 'black', #Default: 'black'
	'y_tick_height' => 5, #Default: 5
	'y_mark_font' => 'Arial', #Default: 'Arial'
	'y_mark_size' => 10, ### 0 = no marks #Default: 
	'y_legend_font' => 'Arial', #Default: 'Arial'
	'y_legend_size' => 10, ### 0 = no marks #Default: 5
	'y_legend_color' => 'black', #Default: 'black'
	
	'grid_line_width' => 1, #Default: 1
	'grid_line_color' => 'grey', #Default: 'grey'
	
	'synteny_line_width' => 5, #Default: 2
	'synteny_forward_line_color' => 'red', #Default: 'red'
	'synteny_reverse_line_color' => 'blue', # Default: 'blue'
);
### X axis order
my @xorder=();
### Y axis order
my @yorder=();



### Get seq length
if ($inputfasta=~/(\.fa$)|(\.fas$)|\.fasta$/i) {
	unless (-s "$inputfasta.fai") {
		unless (IndexFasta($inputfasta)) {
			die "Error: creating fasta index failed\n";
		}
	}
	$inputfasta.='.fai';
}
elsif ($inputfasta=~/\.fai$/i) {
	###do something
}
else {
	die "Error: unknown fasta format, should be fasta or fasta.fai\n";
}
open (FASTAINDEX, "< $inputfasta") || die "Error: invalid fasta or fasta index $inputfasta\n";
while (my $line=<FASTAINDEX>) {
	chomp $line;
	$linenum++;
	my @arr=split(/\t/, $line);
	unless (scalar(@arr)>=2) {
		die "Error: invalid fasta.fai line($linenum): $line\n";
	}
	$seqlength{$arr[0]}=$arr[1];
}
close FASTAINDEX;
print "\n\nLENGTH SUMMARY\n\tTotal sequences: $linenum\n\tTotal hash: ", scalar(keys %seqlength), "\n";
#print "Test: \%seqlength\n"; print Dumper \%seqlength;print "\n";### For test ###


open (SYNTENY, "< $inputconfig") || die "Error: can not open synteny config file\n";
$linenum=0;
my %hash_a=();
my %hash_b=();
while (my $line=<SYNTENY>) {
	chomp $line;
	next if ($line=~/^#/);
	$linenum++;
	my @arr=split(/\t/, $line);
	unless (exists $hash_a{$arr[0]}) {
		unless (exists $seqlength{$arr[0]}) {
			die "Error: sequence in X no length: $arr[0]\n";
		}
		push (@yorder, $arr[0]);
	}
	$hash_a{$arr[0]}++;
	unless (exists $hash_b{$arr[3]}) {
		unless (exists $seqlength{$arr[3]}) {
			die "Error: sequence in Y no length: $arr[3]\n";
		}
		push (@xorder, $arr[3]);
	}
	$hash_b{$arr[3]}++;
}
%hash_a=();
%hash_b=();
close SYNTENY;
if ($linenum==0) {
	print STDERR "Warnings: empty synteny file\n";
	exit 0;
}
die "Error: empty X axis sequences\n" unless (scalar(@xorder)>0);
die "Error: empty Y axis sequences\n" unless (scalar(@yorder)>0);




my $total_x_length=0;
my %starting_x=();
foreach my $xseq (@xorder) {
	$starting_x{$xseq}=$total_x_length+1;
	if (exists $seqlength{$xseq} and $seqlength{$xseq}=~/^\d+$/) {
		$total_x_length+=$seqlength{$xseq};
	}
	else {
		die "Error: X sequence $xseq NO length\n";
	}
}
my $x_factor=($confighash{'plot_width'}-$confighash{'plot_left_margin'}-$confighash{'plot_right_margin'})/$total_x_length;
print "\n\nX SUMMARY\n\tTotal X length: $total_x_length\n\tX factor: $x_factor\n";


my $total_y_length=0;
my %starting_y=();
foreach my $yseq (@yorder) {
	$starting_y{$yseq}=$total_y_length+1;
	if (exists $seqlength{$yseq} and $seqlength{$yseq}=~/^\d+$/) {
		$total_y_length+=$seqlength{$yseq};
	}
	else {
		die "Error: Y sequence $yseq NO length\n";
	}
}
my $y_factor=($confighash{'plot_height'}-$confighash{'plot_top_margin'}-$confighash{'plot_bottom_margin'})/$total_y_length;
print "\n\nY SUMMARY\n\tTotal Y length: $total_y_length\n\tY factor: $y_factor\n";



my $vectorout=SVG->new(width=>$confighash{'plot_width'}, height=>$confighash{'plot_height'});

###Draw axis
$vectorout->line (   id=> "xaxis",
					x1 => $confighash{'plot_left_margin'},
					y1 => $confighash{'plot_height'}-$confighash{'plot_bottom_margin'}, 
					x2 => $confighash{'plot_width'}-$confighash{'plot_right_margin'},
					y2 => $confighash{'plot_height'}-$confighash{'plot_bottom_margin'}, 
					stroke => $confighash{'x_line_color'}, 
					"stroke-width" => $confighash{'x_line_height'}
);
$vectorout->line (   id=> "yaxis",
							x1 => $confighash{'plot_left_margin'},
							y1 => $confighash{'plot_top_margin'}, 
							x2 => $confighash{'plot_left_margin'}, 
							y2 => $confighash{'plot_height'}-$confighash{'plot_bottom_margin'}, 
							stroke => $confighash{'y_line_color'}, 
							"stroke-width" => $confighash{'y_line_width'}
						
);

#print "Test: \%starting_x\n"; print Dumper \%starting_x;print "\n";### For test ###
#print "Test: \%starting_y\n"; print Dumper \%starting_y;print "\n";### For test ###

### Draw ticks



### Draw grid
print "Info: Drawing grid lines\n";
foreach (keys %starting_y) {
	$vectorout->text(x => $confighash{'plot_left_margin'}-2,
					y => $confighash{'plot_height'}-$confighash{'plot_bottom_margin'}-($starting_y{$_}+$seqlength{$_}-1)*$y_factor, 
					"width" => $confighash{'y_legend_size'},
					"height" => $confighash{'y_legend_size'},
					"font-family"=>$confighash{'y_legend_font'},
					"text-anchor"=>"start",
					"text-color"=>$confighash{'y_legend_color'},
					"font-size"=>$confighash{'y_legend_size'},
					"style"=> "writing-mode: tb; glyph-orientation-vertical: 1;",
					"-cdata" => "$_");
	next unless (exists $starting_y{$_} and $starting_y{$_}>1);
	$vectorout->line (   id=> "y_grid_$starting_y{$_}",
					x1 => $confighash{'plot_left_margin'},
					y1 => $confighash{'plot_height'}-$confighash{'plot_bottom_margin'}-($starting_y{$_}-1)*$y_factor, 
					x2 => $confighash{'plot_width'}-$confighash{'plot_right_margin'},
					y2 => $confighash{'plot_height'}-$confighash{'plot_bottom_margin'}-($starting_y{$_}-1)*$y_factor,, 
					stroke => $confighash{'grid_line_color'},
					"stroke-width" => $confighash{'grid_line_width'}
	);
}
foreach (keys %starting_x) {
	$vectorout->text(x => $starting_x{$_}*$x_factor+$confighash{'plot_left_margin'},
					y => $confighash{'plot_height'}-$confighash{'plot_bottom_margin'}+$confighash{'x_legend_size'}, 
					"width" => $confighash{'x_legend_size'},
					"height" => $confighash{'x_legend_size'},
					"font-family"=>$confighash{'x_legend_font'},
					"text-anchor"=>"start",
					"text-color"=>$confighash{'x_legend_color'},
					"font-size"=>$confighash{'x_legend_size'},
					"-cdata" => "$_");
	next unless (exists $starting_x{$_} and $starting_x{$_}>1);
	$vectorout->line (   id=> "x_grid_$starting_x{$_}",
						x1 => $confighash{'plot_left_margin'}+($starting_x{$_}-1)*$x_factor,
						y1 => $confighash{'plot_top_margin'}, 
						x2 => $confighash{'plot_left_margin'}+($starting_x{$_}-1)*$x_factor, 
						y2 => $confighash{'plot_height'}-$confighash{'plot_bottom_margin'}, 
						stroke => $confighash{'grid_line_color'}, 
						"stroke-width" => $confighash{'y_line_width'}
	);
}


### Draw synteny lines
print "Info: Draw synteny lines\n";
open (SYNTENY, "< $inputconfig") || die "Error: can not open synteny config file\n";
$linenum=0;
while (my $line=<SYNTENY>) {
	chomp $line;
	$linenum++;
	next if ($line=~/^#/);
	my @arr=split(/\t/, $line);
	if ($arr[7] eq '+') {
		$vectorout->line (   id => "$linenum-$arr[0]-$arr[3]",
						x1 => $confighash{'plot_left_margin'}+($starting_x{$arr[3]}+$arr[4])*$x_factor,
						y1 => $confighash{'plot_height'}-($starting_y{$arr[0]}+$arr[1])*$y_factor-$confighash{'plot_bottom_margin'},
						x2 => $confighash{'plot_left_margin'}+($starting_x{$arr[3]}+$arr[5])*$x_factor,
						y2 => $confighash{'plot_height'}-($starting_y{$arr[0]}+$arr[2])*$y_factor-$confighash{'plot_bottom_margin'},
						stroke => $confighash{'synteny_forward_line_color'}, 
						"stroke-width" => $confighash{'synteny_line_width'}
		);
	}
	elsif ($arr[7] eq '-') {
		$vectorout->line (   id => "$linenum-$arr[0]-$arr[3]",
						x1 => $confighash{'plot_left_margin'}+($starting_x{$arr[3]}+($seqlength{$arr[3]}-$arr[5]+1))*$x_factor,
						y1 => $confighash{'plot_height'}-($starting_y{$arr[0]}+$arr[2])*$y_factor-$confighash{'plot_bottom_margin'},
						x2 => $confighash{'plot_left_margin'}+($starting_x{$arr[3]}+($seqlength{$arr[3]}-$arr[4]+1))*$x_factor,
						y2 => $confighash{'plot_height'}-($starting_y{$arr[0]}+$arr[1])*$y_factor-$confighash{'plot_bottom_margin'},
						stroke => $confighash{'synteny_reverse_line_color'}, 
						"stroke-width" => $confighash{'synteny_line_width'}
		);
	}
	else{
		die "Error: invalid senteny line: $line\n";
	}
}
close SYNTENY;
my $finalout = $vectorout->xmlify;
open SVGFILE, "> $outputsvg";
print SVGFILE $finalout;
close SVGFILE;
