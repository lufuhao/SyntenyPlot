#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper qw/Dumper/;
use SVG;
use constant USAGE =><<EOH;

usage: $0 synteny.config output.svg

Descriptions:
	Produce publishable vector-format syntonic gene expression views

Requirements: 
	Perl Modules: SVG, Data::Dumper

synteny.config [2 column data set]
#Line start with # would be ignored
#At the end mark one more line with chromosome length (eg: 1001) and value 0
#    to indicate the precise ploting on chromsome
#gene positions are not necessarily unique
#gene_position	expression_value
100	1.25
200	8.00
254	-4.00
1001	0

Version:
	v20180419

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


#####################################################
#############  Configure start ######################
#####################################################
my $debug=0;
my %confighash=(
##[plot]
	'plot_height' => 200,
	'plot_width' => 500,
	'plot_margin_left' => 20,
	'plot_margin_right' => 10,
	'plot_margin_top' => 10,
	'plot_margin_bottom' => 10,
	'plot_line_width' => 1,
	'plot_background_color' => 'white',
	'plot_font_family' => 'Arial',
	'plot_track_space' => 10,
##[synteny]
	'synteny_chrom_height' => 10,
	'synteny_chrom_fillin_color' => 'grey',
	'synteny_chrom_line_color' => 'black',
	'synteny_chrom_line_width' => 1,
	'synteny_chrom_length' => 371759578, #chromsome length or will use largest position in config if not specified
##[ruler]
	'ruler_mark_interval' => 2,
	'ruler_tick_color' => 'black',
	'ruler_tick_length' => 5,
	'ruler_tick_width' => 1,
	'ruler_font_size' => 12,
	'ruler_font_distance' => 2,
	'ruler_text_color' => 'black',
	'ruler_x_line_color' => 'black',
	'ruler_x_line_width' => 2,
	'ruler_y_line_color' => 'black',
	'ruler_y_line_width' => 2,
	'ruler_y_top_scalar' => 10,  ### will auto detect when set to '#N/A'
	'ruler_y_down_scalar' => -10, ### will auto detect when set to '#N/A'
###[expression]
#	[start(>=), end(<=), line_color, marker_color]
#	will use synteny_chrom_line_color if not set
	'exprs_color_array' => [[-12, 12, 'grey', 'grey']],
	'exprs_up_line_size' => 0.3,
	'exprs_up_marker_size' => 0.6,
	'exprs_down_line_size' => 0.3,
	'exprs_down_marker_size' => 0.6,
);
#$confighash{'ruler_y_top_scalar'}
#####################################################
#############  Configure end ########################
#####################################################
# unless you understand the code, leave below alone #
#####################################################



### Input and output
my ($configfile, $outputsvg)=@ARGV;
unless (defined $configfile and -s $configfile) {
	die "Error: invalid config file\n";
}
unless (defined $outputsvg and $outputsvg=~/^\S+$/) {
	die "Error: invalid output.svg file\n";
}
unlink "$outputsvg" if (-e "$outputsvg");



### read configure file
my %pos2exp=();
my $linenum=0;
my $validnum=0;
my $tick_up_scale=0;
my $tick_down_scale=0;
my $tick_up_autodetect=0;
my $tick_down_autodetect=0;
my $scale_min=0;
my $scale_max=0;
if (exists $confighash{'ruler_y_top_scalar'} and $confighash{'ruler_y_top_scalar'}=~/^[-+]{0,1}\d+\.*\d*$/) {
	print "Info: Y-axis top scale was set to ", $confighash{'ruler_y_top_scalar'}, "\n";
	$tick_up_scale=$confighash{'ruler_y_top_scalar'};
}
else {
	$tick_up_autodetect=1;
}
if (exists $confighash{'ruler_y_down_scalar'} and $confighash{'ruler_y_down_scalar'}=~/^[-+]{0,1}\d+\.*\d*$/) {
	print "Info: Y-axis top scale was set to ", $confighash{'ruler_y_down_scalar'}, "\n";
	$tick_down_scale=$confighash{'ruler_y_down_scalar'};
}
else {
	$tick_down_autodetect=1;
}

close CONFIGUREINPUT if (defined fileno(CONFIGUREINPUT));
unless (open CONFIGUREINPUT, "<", $configfile) {
	die "Error: can not open configure file\n";
}
while (my $line=<CONFIGUREINPUT>) {
	chomp $line;
	$linenum++;
	next if ($line=~/^#/);
	my @arr=split(/\t/, $line);
	unless (scalar(@arr)==2) {
		die "Error: NumCol!=2 at line ($linenum): $line\n";
	}
	unless (defined $arr[0] and $arr[0]=~/^\d+$/ and $arr[0]>0) {
		die "Error: Col1 not INT at line ($linenum): $line\n";
	}
	unless (defined $arr[1] and $arr[1]=~/^[\-+]{0,1}\d+\.*\d*$/) {
		print STDERR "Warnings: ignored line due to Col2 not number at line ($linenum): $line\n";
		next;
	}
#	if (exists $pos2exp{$arr[0]}) {
#		die "Error: duplicated position at line ($linenum): $line\n"
#	}
	$pos2exp{$arr[0]}{$arr[1]}++;
	if ($tick_up_autodetect==1) {
		if ($validnum==0) {
			$scale_max=$arr[1];
		}
		elsif ($arr[1]>$scale_max) {
			$scale_max=$arr[1];
		}
	}
	if ($tick_down_autodetect==1) {
		if ($validnum==0) {
			$scale_min=$arr[1];
		}
		elsif ($arr[1]<$scale_min) {
			$scale_min=$arr[1];
		}
	}
	$validnum++;
}
close CONFIGUREINPUT;
print "Info: ##########SUMMARY#############\n";
print "Info: Total lines: $linenum\n";
print "Info: Valid lines: $validnum\n";



###################### Plot #############################
my $vectorout=SVG->new(width=>$confighash{'plot_width'}, height=>$confighash{'plot_height'});

### check space
my $chrom_x_length=$confighash{'plot_width'}-$confighash{'plot_margin_left'}-$confighash{'plot_margin_right'}-
2*$confighash{'ruler_tick_length'};
unless ($chrom_x_length>0) {
	die "Error: not enough width to plot, try to increase plot_width or decrease left or right margin\n";
}
my $scale_y_height=$confighash{'plot_height'}-$confighash{'plot_margin_top'}-$confighash{'plot_margin_bottom'}-2*$confighash{'synteny_chrom_height'}-2*$confighash{'plot_track_space'};
unless ($scale_y_height>0) {
	die "Error: not enough height to plot, try to increase plot_height or decrease top or bottom margin\n";
}

### background
$vectorout-> rectangle (	x => 0, y=> 0, 
							width  	=> $confighash{'plot_width'}, 
							height => $confighash{'plot_height'}, 
							id=> "background",
							style => {'fill' => $confighash{'plot_background_color'},
									'stroke'         => 'rgb(225, 225, 225)',
									'stroke-width'   =>  $confighash{'plot_line_width'},
									'stroke-opacity' =>  0,
									'fill-opacity'   =>  1
							},
						);

### chrome
my $chrom_x_start=$confighash{'plot_margin_left'}+2*$confighash{'ruler_tick_length'};
my $chrom_x_end=$confighash{'plot_width'}-$confighash{'plot_margin_right'};
my $chrom_y_top=$confighash{'plot_margin_top'};
my $chrom_y_bottom=$confighash{'plot_height'}-$confighash{'plot_margin_bottom'}-$confighash{'synteny_chrom_height'};

$vectorout-> rectangle (	x => $chrom_x_start, y=> $chrom_y_top, 
							width  	=> $chrom_x_length, 
							height => $confighash{'synteny_chrom_height'}, 
							id=> "chrom_top",
							style => {'fill' => $confighash{'synteny_chrom_fillin_color'},
									'stroke'         => $confighash{'synteny_chrom_line_color'},
									'stroke-width'   =>  $confighash{'synteny_chrom_line_width'},
									'stroke-opacity' =>  0,
									'fill-opacity'   =>  1
							},
						);

$vectorout-> rectangle (	x => $chrom_x_start, y=> $chrom_y_bottom, 
							width  	=> $chrom_x_length, 
							height => $confighash{'synteny_chrom_height'}, 
							id=> "chrom_bottom",
							style => {'fill' => $confighash{'synteny_chrom_fillin_color'},
									'stroke'         => $confighash{'synteny_chrom_line_color'},
									'stroke-width'   =>  $confighash{'synteny_chrom_line_width'},
									'stroke-opacity' =>  0,
									'fill-opacity'   =>  1
								},
					);


### y - axis
my $scale_y_top=$confighash{'plot_margin_top'}+$confighash{'synteny_chrom_height'}+$confighash{'plot_track_space'};
my $scale_y_bottom=$confighash{'plot_height'}-$confighash{'plot_margin_bottom'}-$confighash{'synteny_chrom_height'}-$confighash{'plot_track_space'};

my $tick_x_end=$chrom_x_start-$confighash{'ruler_y_line_width'}/2;
my $tick_x_start=$tick_x_end-$confighash{'ruler_tick_length'};

$vectorout->line (id=> "y-axis",
			x1 => $tick_x_end,
			y1 => $scale_y_top, 
			x2 => $tick_x_end, 
			y2 => $scale_y_bottom, 
			stroke => $confighash{'ruler_y_line_color'}, 
			"stroke-width" => $confighash{'ruler_y_line_width'}
			);

### $scale_min $scale_max
if ($tick_up_autodetect==1) {
	$tick_up_scale=(int($scale_max/$confighash{'ruler_mark_interval'})+1)*$confighash{'ruler_mark_interval'};
}
if ($tick_down_autodetect==1) {
	$tick_down_scale=(int($scale_min/$confighash{'ruler_mark_interval'})-1)*$confighash{'ruler_mark_interval'};
}

my $total_scales=($tick_up_scale-$tick_down_scale)/$confighash{'ruler_mark_interval'}+1;
if ($debug) {
	print "Top scale: $tick_up_scale\nBottom scale: $tick_down_scale\nTotal ticks: $total_scales\n";
}
my $scale_tick_interval_real_height=$scale_y_height/($total_scales-1);
my $xaxis_scale=0;
if ($scale_min<=0 and $scale_max>=0) {
	$xaxis_scale=0;
}
elsif ($scale_min>=0) {
	$xaxis_scale=$tick_down_scale;
}
elsif ($scale_max<=0) {
	$xaxis_scale=$tick_up_scale;
}
else {
	die "Error: do NOT know where to draw x-axis\n";
}

my $tick_ind_y=$scale_y_bottom;
my $tick_ind_cdata=$tick_down_scale;
my $exprs_y_start=0;
for (my $tick_ind_num=1; $tick_ind_num<=$total_scales; $tick_ind_num++) {
	if ($tick_ind_cdata==$xaxis_scale) {
		$vectorout->line (id=> "x-axis",
			x1 => $chrom_x_start,
			y1 => $tick_ind_y, 
			x2 => $chrom_x_end, 
			y2 => $tick_ind_y, 
			stroke => $confighash{'ruler_x_line_color'}, 
			"stroke-width" => $confighash{'ruler_x_line_width'}
		);
		$exprs_y_start=$tick_ind_y;
	}

	if ($tick_ind_num==1 or $tick_ind_num==$total_scales) {### ignore end scale 
#		next;
	}
	else {
#		draw ticks
		$vectorout->line (id=> "y-ticks-$tick_ind_num",
			x1 => $tick_x_start,
			y1 => $tick_ind_y, 
			x2 => $tick_x_end, 
			y2 => $tick_ind_y, 
			stroke => $confighash{'ruler_x_line_color'}, 
			"stroke-width" => $confighash{'ruler_x_line_width'}
			);
#		draw tick marks
		$vectorout->text(	x => $tick_x_start-$confighash{'ruler_font_distance'}-$confighash{'ruler_tick_length'}, 
				y => $tick_ind_y+$confighash{'ruler_font_size'}/4, 
				width => $confighash{'ruler_font_size'}, 
				height => $confighash{'ruler_font_size'}, 
				"font-family"=>$confighash{'plot_font_family'}, 
				"text-anchor"=>"end",
				"font-size"=>$confighash{'ruler_font_size'}, 
				"-cdata" => "$tick_ind_cdata");
	}
	print "Scale: $tick_ind_cdata\n" if ($debug);
	$tick_ind_y-=$scale_tick_interval_real_height;
	$tick_ind_cdata+=$confighash{'ruler_mark_interval'};
}

### draw expression
###[expression]
#	'exprs_color_array' => [[-1000, -4, 'blue', 'blue'], [-4, 4, 'grey', 'grey'], [4, 1000, 'red', 'red']];
#	'exprs_up_line_size' => 0.1,
#	'exprs_up_marker_size' => 0.2,
#	'exprs_down_line_size' => 0.1,
#	'exprs_down_marker_size' => 0.2,
#	'synteny_chrom_length' => 371759578, 
#);
#$confighash{'exprs_color_array'}
#$chrom_x_start
#$chrom_x_end
my @posarr=sort {$a<=>$b} keys %pos2exp;
my $chrom_length=$posarr[-1];
if (exists $confighash{'synteny_chrom_length'} and $confighash{'synteny_chrom_length'}=~/^\d+$/ and $confighash{'synteny_chrom_length'}>=$chrom_length) {
	$chrom_length=$confighash{'synteny_chrom_length'};
}


foreach my $indpos (sort {$a<=>$b} @posarr) {
	my $ind_exprs_x=$chrom_x_start+($chrom_x_end-$chrom_x_start)*$indpos/$chrom_length;
	my $ind_exprs_line_color='';
	my $ins_exprs_mark_color='';
	foreach my $ind_exp_value (sort {$a<=>$b} keys %{$pos2exp{$indpos}}) {
		if (exists $confighash{'exprs_color_array'}) {
			my $test_color_set=0;
			foreach my $color_array (@{$confighash{'exprs_color_array'}}) {
				if ($ind_exp_value >= $color_array->[0] and $ind_exp_value <= $color_array->[1]) {
					$ind_exprs_line_color=$color_array->[2];
					$ins_exprs_mark_color=$color_array->[3];
					$test_color_set=1;
	#				print "POS $indpos EXP ", $ind_exp_value, " Line $ind_exprs_line_color Marker $ins_exprs_mark_color\n" if ($debug);
					last;
				}
			}
			if ($test_color_set==0) {
				print STDERR "Warnings: color array not set: POSITION $indpos VALUE $ind_exp_value\n";
				print STDERR "          use default synteny_chrom_line_color : ", $confighash{'synteny_chrom_line_color'}, "\n";
				$ind_exprs_line_color=$confighash{'synteny_chrom_line_color'};
				$ins_exprs_mark_color=$confighash{'synteny_chrom_line_color'};
			}
		}
		else {
			$ind_exprs_line_color=$confighash{'synteny_chrom_line_color'};
			$ins_exprs_mark_color=$confighash{'synteny_chrom_line_color'};
		}
		my $ins_exprs_y_end=$exprs_y_start-($ind_exp_value-$xaxis_scale)*$scale_y_height/($tick_up_scale-$tick_down_scale);
		if ($ind_exp_value>$tick_up_scale) {### control extra high
			$ins_exprs_y_end=$exprs_y_start-($tick_up_scale-$xaxis_scale)*$scale_y_height/($tick_up_scale-$tick_down_scale)
		}
		if ($ind_exp_value<$tick_down_scale) {### control extra low
			$ins_exprs_y_end=$exprs_y_start-($tick_down_scale-$xaxis_scale)*$scale_y_height/($tick_up_scale-$tick_down_scale)
		}

		if ($ind_exp_value>$xaxis_scale) {
			$vectorout->line (id=> "exprs-line-pos$indpos-expr$ind_exp_value",
					x1 => $ind_exprs_x,
					y1 => $exprs_y_start,
					x2 => $ind_exprs_x, 
					y2 => $ins_exprs_y_end,
					stroke => $ind_exprs_line_color, 
					"stroke-width" => $confighash{'exprs_up_line_size'}
					);
			$vectorout->circle (id=> "exprs-marker-pos$indpos-expr$ind_exp_value",
					cx => $ind_exprs_x,
					cy => $ins_exprs_y_end,
					r => $confighash{'exprs_up_marker_size'},
					style => {
						'stroke' => $ins_exprs_mark_color,
						"stroke-width" => $confighash{'exprs_up_line_size'}/8,
						'fill' =>$ins_exprs_mark_color,
						'stroke-opacity' => 1,
						'fill-opacity' => 1
						}
					);
		}
		elsif ($ind_exp_value<$xaxis_scale) {
			$vectorout->line (id=> "exprs-line-pos$indpos-expr$ind_exp_value",
					x1 => $ind_exprs_x,
					y1 => $exprs_y_start, 
					x2 => $ind_exprs_x, 
					y2 => $ins_exprs_y_end, 
					stroke => $ind_exprs_line_color, 
					"stroke-width" => $confighash{'exprs_down_line_size'}
					);
			$vectorout->circle (id=> "exprs-marker-pos$indpos-expr$ind_exp_value",
					cx => $ind_exprs_x,
					cy => $ins_exprs_y_end,
					r => $confighash{'exprs_down_marker_size'},
					style => {
						'stroke' => $ins_exprs_mark_color,
						"stroke-width" => $confighash{'exprs_down_line_size'}/8,
						'fill' => $ins_exprs_mark_color,
						'stroke-opacity' => 1,
						'fill-opacity' => 1
						}
					);
		}
	}
}



my $finalout = $vectorout->xmlify;
open SVGFILE, "> $outputsvg";
print SVGFILE $finalout;
close SVGFILE;
