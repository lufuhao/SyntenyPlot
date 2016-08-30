#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper qw/Dumper/;
use SVG;
use constant USAGE =><<EOH;

usage: $0 global.config synteny.config annotation.config output.svg

Descriptions:
	Produce publishable vector-format synteny views using user-defined parameter and datasets

Requirements: 
	Perl Modules: SVG, Data::Dumper

Version:
	v20160830

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk

EOH
die USAGE if (scalar(@ARGV) !=4 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



my ($userconfig, $syntenyconfig, $annotationconfig, $outputsvg);
unless (defined $ARGV[0] and -s $ARGV[0]) {
	die "Error: invalid global.config file\n";
}
$userconfig=$ARGV[0];
unless (defined $ARGV[1] and -s $ARGV[1]) {
	die "Error: invalid synteny.config file\n";
}
$syntenyconfig=$ARGV[1];
unless (defined $ARGV[2] and -s $ARGV[2]) {
	die "Error: invalid annotation.config file\n";
}
$annotationconfig=$ARGV[2];
unless (defined $ARGV[3] and $ARGV[3]=~/^\S+$/) {
	die "Error: invalid output.svg file\n";
}
unlink "$ARGV[3]" if (-e $ARGV[3]);
$outputsvg=$ARGV[3];



### Defaults
my %userdefault=();
my %confighash=();
my %syntenyhash=();
my %annotationhash=();
my @blockorder=();
my @tissueorder=();
my @trackoders=();
my %seqlength=();
my %alltrackheight=();
my %express_location=();
my $primaryloci=[];
my %loci=();
my %tissue_lengend=();
my $use_defined_order=0;
my $use_defined_tissueorder=0;
my $use_defined_trackorder=0;
my $auto_calc_expression_height=0;
my $track_remaining_space_height=0;
my $expression_top_space=3;


### Default
#%confighash= {'plot_height' => 500, 'plot_width' => 2000};
%userdefault=(
##[plot]
	'plot_height' => 500,
	'plot_width' => 2000,
	'plot_margin_left' => 10,
	'plot_margin_right' => 10,
	'plot_margin_top' => 10,
	'plot_margin_bottom' => 10,
	'plot_line_size' => 1,
	'plot_background_color' => 'white',
#	'plot_font_family' => 'Arial',
##[synteny]
	'synteny_height' => 50,
	'synteny_fillin_color' => 'rgb(204,233,223)',
	'synteny_line_color' => 'black',
	'synteny_line_width' => 1,
	'synteny_order' => 'ToBeDetected',
##[ruler]
	'ruler_mark_interval' => 5000,
	'ruler_tick_color' => 'black',
	'ruler_tick_length' => 14,
	'ruler_font_size' => 14,
	'ruler_tick_width' => 1,
	'ruler_font_distance' => 2,
	'ruler_text_color' => 'black',
##[track]
	'track_height' => 5,
	'track_space' => 3,
	'track_color' => 'black',
	'track_arrow_add' => 2,
	'track_arrow_length' => 15,
	'track_order' => 'ToBeDetected',
	'track_subtrack_space_height' => 3,
##[expression]
### [expression][gridline]
	'gridline_width' => 1,
	'gridline_color' => 'grey',
	'gridline_interval' => 2,
	'gridline_min' => 0,
	'gridline_max' => 6,
	'gridline_left_margin' => 15,
	'gridline_right_margin' => 5,
	'gridtext_color' => 'black',
	'gridtext_size' => 12,
### [expression][barplot]
	'barplot_order' => 'ToBeDetected',
	'bar_line_width' => 1,
	'bar_alignment' => 'MIDDLE',
	'bar_fill_width' => 15,
### [expression][stdev]
	'error_line_width' => 1,
	'error_line_length' => 5,
	'error_line_color' => 'black',
### [expression][LEGEND]
	'legend_plot_position' => 'TOP',
	'legend_text_color' => 'black',
	'legend_text_size' => 14,
	'legend_plot_aplignment' => 'RIGHT',
	'legend_between_interval' => 9,
);



### Reading configs
### read user config
unless (&readUsrCfg($userconfig)) {
	die "(MAIN)Error: unable to read user.config: $userconfig\n";
}
if ($confighash{'synteny_order'} eq 'ToBeDetected') {
	$use_defined_order=0;
}
else {
	$use_defined_order=1;
	@blockorder=split(/,/, $confighash{'synteny_order'});
}
if ($confighash{'track_order'} eq 'ToBeDetected') {
	print "Info: TRACK order not user-defined\n";
	$use_defined_trackorder=0;
}
else {
	print "Info: TRACK order user-defined\n";
	$use_defined_trackorder=1;
	@trackoders=split(/,/, $confighash{'track_order'});
	foreach (my $i=0; $i<scalar(@trackoders); $i++) {$trackoders[$i]=uc($trackoders[$i]);}
}
if ($confighash{'barplot_order'} eq 'ToBeDetected') {
	$use_defined_tissueorder=0;
}
else {
	$use_defined_tissueorder=1;
	@tissueorder=split(/,/, $confighash{'barplot_order'});
}

### read synteny config
unless (&readSyntenyCfg($syntenyconfig)) {
	die "(MAIN)Error: unable to read synteny.config: $syntenyconfig\n";
}
### read annotations config
unless (&readAnnotationCfg($annotationconfig)) {
	die "(MAIN)Error: unable to read annotation.config: $annotationconfig\n";
}



### Print Input SUMMARY ###
print "\n\n\n###SUMMARY: ###\n";
print "BLOCK order (predefined: $use_defined_order TotalElement: ", scalar(@blockorder), "): ", join(',', @blockorder), "\n";
print "TRACK order (predefined: $use_defined_trackorder TotalElement: ", scalar(@trackoders), "): ", join(',', @trackoders), "\n";
print "TRACK order (predefined: $use_defined_tissueorder TotalElement: ", scalar(@tissueorder), "): ", join(',', @tissueorder), "\n";
print "\n\n\n";



### Estimate track height
foreach my $thisblock (keys %seqlength) {
	if (exists $confighash{"ruler_${thisblock}_length"}) {
		unless ($confighash{"ruler_${thisblock}_length"}=~/^\d+$/) {
			print STDERR "Warnings: invalid ruler_${thisblock}_length parameter in user.config file, use autodetected in synteny and annotation files: $seqlength{$thisblock}\n";
			$confighash{"ruler_${thisblock}_length"}=$seqlength{$thisblock};
		}
		else {
			if ($confighash{"ruler_${thisblock}_length"}< $seqlength{$thisblock}) {
				print STDERR "Warnings: There is BLOCK $thisblock POSITION $seqlength{$thisblock} larger than the length defined in user.config ruler_${thisblock}_length=", $confighash{"ruler_${thisblock}_length"}, "\n";
			}
		}
	}
	else {
		$confighash{"ruler_${thisblock}_length"}=$seqlength{$thisblock};
	}
}
%seqlength=();

foreach my $idvblock (@blockorder) {
	foreach my $idvtrack (@trackoders) {
		my $subtracknum=0;
		if ($idvtrack eq 'EXPRESSION') {
			if (exists $annotationhash{$idvblock}{'EXPRESSION'}) {
				$alltrackheight{$idvblock}{$idvtrack}{'subtracks'}=1;
			}
			else {
				$alltrackheight{$idvblock}{$idvtrack}{'subtracks'}=0;
				$alltrackheight{$idvblock}{$idvtrack}{'height'}=0;
				next;
			}
			if (exists $confighash{"track_${idvblock}_EXPRESSION_height"}) {
				$alltrackheight{$idvblock}{$idvtrack}{'height'}=$confighash{"track_${idvblock}_EXPRESSION_height"};
			}
			elsif (exists $confighash{'track_EXPRESSION_height'}) {
				$alltrackheight{$idvblock}{$idvtrack}{'height'}=$confighash{'track_EXPRESSION_height'};
				$confighash{"track_${idvblock}_EXPRESSION_height"}=$confighash{'track_EXPRESSION_height'};
			}
			else {
				#$auto_calc_expression_height++;
			}
			next;
		}
		next unless (exists $annotationhash{$idvblock} and exists $annotationhash{$idvblock}{$idvtrack});
		my $layouttype='SEPARATE';
		if (exists $confighash{"track_${idvblock}_${idvtrack}_layout"}) {
			if ($confighash{"track_${idvblock}_${idvtrack}_layout"} =~/(^overlap$)|(^separate$)/i) {
				$layouttype=uc($confighash{"track_${idvblock}_${idvtrack}_layout"});
			}
			else {
				die "Error: invalid track_${idvblock}_${idvtrack}_layout value: ", $confighash{"track_${idvblock}_${idvtrack}_layout"}, " for BLCOK $idvblock TRACK $idvtrack\n";
			}
		}
		elsif (exists $confighash{"track_${idvtrack}_layout"}) {
			if ($confighash{"track_${idvtrack}_layout"} =~/(^overlap$)|(^separate$)/i) {
				$layouttype=uc($confighash{"track_${idvtrack}_layout"});
				
			}
			else {
				die "Error: invalid track_${idvtrack}_layout value: ", $confighash{"track_${idvtrack}_layout"}, " for BLCOK $idvblock TRACK $idvtrack\n";
			}
		}
#		print "Info: BLOCK ${idvblock} TRACK ${idvtrack} layout: $layouttype\n";### For test ###
		if ($layouttype =~/^overlap$/i) {
			$subtracknum=1;
		}
		elsif ($layouttype =~/^separate$/i) {
#			print "Info: BLOCK ${idvblock} TRACK ${idvtrack} layout: assign subtrack\n";### For test ###
			(my $test, $subtracknum)=&AssignSubtrack($idvblock, $idvtrack);
			unless ($test) {
				die "Error: unable to assign subtracks for BLCOK $idvblock TRACK $idvtrack\n";
			}
		}
		else {
			die "Error: invalid layouttype: $layouttype for BLCOK $idvblock TRACK $idvtrack";
		}
		$confighash{"track_${idvblock}_${idvtrack}_layout"}=$layouttype;
		print "Info: BLOCK ${idvblock} TRACK ${idvtrack} layout: $layouttype Subtracks: $subtracknum\n";
		if (exists $confighash{"track_${idvblock}_${idvtrack}_height"}) {
			if ($confighash{"track_${idvblock}_${idvtrack}_height"} =~/^\d+$/) {
				$alltrackheight{$idvblock}{$idvtrack}{'height'}=$confighash{"track_${idvblock}_${idvtrack}_height"};
			}
			else {
				die "Error: invalid track_${idvblock}_${idvtrack}_height parameter in global.config: $userconfig\n";
			}
		}
		elsif (exists $confighash{"track_${idvtrack}_height"}) {
			if ($confighash{"track_${idvtrack}_height"} =~/^\d+$/) {
				$alltrackheight{$idvblock}{$idvtrack}{'height'}=$confighash{"track_${idvtrack}_height"};
				$confighash{"track_${idvblock}_${idvtrack}_height"}=$confighash{"track_${idvtrack}_height"};
			}
			else {
				die "Error: invalid track_${idvtrack}_height parameter in global.config: $userconfig\n";
			}
		}
		else {
			$alltrackheight{$idvblock}{$idvtrack}{'height'}=$confighash{'track_height'};
			$confighash{"track_${idvblock}_${idvtrack}_height"}=$confighash{'track_height'};
		}
		$alltrackheight{$idvblock}{$idvtrack}{'subtracks'}=$subtracknum;
	}
}
if (0) {
	print "\n";
	print Dumper \%annotationhash;
	print "\n\n\n";
}



### Total Y Drawing block height
my $trackdrawspace=$confighash{'plot_height'}-$confighash{'synteny_height'}*(scalar(@blockorder)-1);
$trackdrawspace=$trackdrawspace-$confighash{'plot_margin_top'};
$trackdrawspace=$trackdrawspace-$confighash{'plot_margin_bottom'};
$trackdrawspace=$trackdrawspace-$confighash{'ruler_tick_length'}*(scalar(@blockorder))-$confighash{'ruler_tick_length'}*(scalar(@blockorder))/2;
$trackdrawspace=$trackdrawspace-$confighash{'track_space'}*(scalar(@blockorder))*2;
$trackdrawspace=$trackdrawspace-$confighash{'plot_line_size'}*(scalar(@blockorder));
#print "Test: Total Y writable space: $confighash{'plot_height'}-$confighash{'plot_line_size'}*(scalar(@blockorder))-$confighash{'synteny_height'}*(scalar(@blockorder)-1)-$confighash{'ruler_tick_length'}*(scalar(@blockorder))-$confighash{'plot_margin_top'}-$confighash{'plot_margin_bottom'}-$confighash{'track_space'}*(scalar(@blockorder))*2 = $trackdrawspace\n"; ### for test ###
print "Test: Total Y writable space: $trackdrawspace\n"; ### for test ###


###Double Check if every track gets its height
my %undeterminedtrackheight=();
my $total_trackheight=0;
foreach my $idvblock (@blockorder) {
	foreach my $idvtrack (@trackoders) {
		unless (exists $alltrackheight{$idvblock} and exists $alltrackheight{$idvblock}{$idvtrack}) {
			die "Error: non-existing alltrackheight for BLOCK $idvblock TRACK $idvtrack\n";
		}
		unless(exists $alltrackheight{$idvblock}{$idvtrack}{'subtracks'} and $alltrackheight{$idvblock}{$idvtrack}{'subtracks'}=~/^\d+$/) {
			die "Error: invalid track number for BLOCK $idvblock TRACK $idvtrack\n";
		}
		unless (exists $alltrackheight{$idvblock}{$idvtrack}{'height'} and $alltrackheight{$idvblock}{$idvtrack}{'height'}=~/^\d+$/) {
			if ($idvtrack eq 'EXPRESSION') {
				$undeterminedtrackheight{$idvblock}{$idvtrack}++;
				$auto_calc_expression_height++;
			}
			else {
				die "Error: invalid track height for BLOCK $idvblock TRACK $idvtrack\n";
			}
		}
		else {
			$total_trackheight += $alltrackheight{$idvblock}{$idvtrack}{'subtracks'} * $alltrackheight{$idvblock}{$idvtrack}{'height'} + $confighash{'track_subtrack_space_height'} * ($alltrackheight{$idvblock}{$idvtrack}{'subtracks'}-1);
		}
	}
	$total_trackheight+=$expression_top_space+$confighash{'track_space'}*(scalar(@trackoders)-1);
}

if ($total_trackheight > $trackdrawspace) {
	die "Error: not enough space to draw tracks. \nConsider to reduce synteny_height, or track height\n";
}
else {print "Test: total configure height: $total_trackheight\n";} ### For test ###



### Auto detect expression height
if ($auto_calc_expression_height) {
	print "\n\n\nInfo: automatically calculate undefined TRACK height\n";
	print "Info: total remaining height = ", ($trackdrawspace-$total_trackheight), "\n";
	print "Info: total remaining TRACK = $auto_calc_expression_height\n\n\n";
	my $autodefined_expression_height=($trackdrawspace-$total_trackheight)/$auto_calc_expression_height;
	foreach my $idvblock (keys %undeterminedtrackheight) {
		$alltrackheight{$idvblock}{'EXPRESSION'}{'height'}=$autodefined_expression_height;
	}
}
else {
	$track_remaining_space_height=($trackdrawspace-$total_trackheight)/(scalar(@blockorder)-1);
	print "Test: remaining height: $track_remaining_space_height\n";
}
%undeterminedtrackheight=();



###Final check \%alltrackheight
if (0) {
	print "\n\%alltrackheight\n";
	print Dumper \%alltrackheight;
	print "\n\n\n";
}
foreach my $idvblock (@blockorder) {
	foreach my $idvtrack (@trackoders) {
		unless (exists $alltrackheight{$idvblock}{'EXPRESSION'}{'subtracks'} and $alltrackheight{$idvblock}{'EXPRESSION'}{'subtracks'}=~/^\d+$/) {
			die "Error: no subtracks number BLOCK $idvblock TRACK $idvtrack\n";
		}
		unless (exists $alltrackheight{$idvblock}{'EXPRESSION'}{'height'} and $alltrackheight{$idvblock}{'EXPRESSION'}{'height'}=~/^\d+\.?\d*$/) {
			print $alltrackheight{$idvblock}{'EXPRESSION'}{'height'}, "\n";
			die "Error: no subtracks height BLOCK $idvblock TRACK $idvtrack\n";
		}
	}
}



### Drawable Y space
my $horizontal_drawing_space=$confighash{'plot_width'}- 2 * $confighash{'plot_line_size'}-$confighash{'plot_margin_left'}-$confighash{'plot_margin_right'};



###starting drawing
my $vectorout=SVG->new(width=>$confighash{'plot_width'}, height=>$confighash{'plot_height'});
### draw plot
$vectorout->rectangle(x => 0, y => 0, 
				   width  	=> $confighash{'plot_width'}, 
				   height => $confighash{'plot_height'}, 
				   id=> "background",
				   style => {'fill' => $confighash{'plot_background_color'},
                             'stroke'         => 'rgb(225, 225, 225)',
                             'stroke-width'   =>  0,
                             'stroke-opacity' =>  0,
                             'fill-opacity'   =>  1,
                    		},
                );
### draw blocks
my $starttingX=$confighash{'plot_margin_left'};
my $starttingY=$confighash{'plot_margin_top'};
my %block_coordinates=();
for (my $i=0; $i<scalar(@blockorder); $i++) {
	my $idvblock=$blockorder[$i];
	my $block_total_height=$confighash{'plot_line_size'}+$track_remaining_space_height+$confighash{'track_space'}+ $confighash{'ruler_tick_length'}+$confighash{'ruler_tick_length'}/2;
	for (my $j=(scalar(@trackoders)-1); $j>=0; $j--) {
		my $idvtrack=$trackoders[$j];
		unless (exists $alltrackheight{$idvblock}{$idvtrack}{'subtracks'}) {
			die "Error: no subtracks for BLOCK $idvblock TRACK $idvtrack\n";
		}
		next unless ($alltrackheight{$idvblock}{$idvtrack}{'subtracks'}>0);
		$block_total_height += $alltrackheight{$idvblock}{$idvtrack}{'height'}*$alltrackheight{$idvblock}{$idvtrack}{'subtracks'} + ($alltrackheight{$idvblock}{$idvtrack}{'subtracks'}-1)*$confighash{'track_subtrack_space_height'};
		$block_total_height +=$confighash{'track_space'};
	}
###Drawing block frame
	print "BLOCK $idvblock X $starttingX Y $starttingY Total Height $block_total_height\n";
	$block_coordinates{$idvblock}{'top'}=$starttingY;#For Synteny
	$block_coordinates{$idvblock}{'bottom'}=$starttingY+$block_total_height;### For synteny
	my $xv=[]; my $yv=[];
	my $points={};
#	$xv=[$confighash{'plot_margin_left'}, $confighash{'plot_width'}-$confighash{'plot_margin_right'}, $confighash{'plot_width'}-$confighash{'plot_margin_right'}, $confighash{'plot_margin_left'}];
#	$yv=[$starttingY, $starttingY, $starttingY+$block_total_height, $starttingY+$block_total_height];
#	$points=$vectorout->get_path( x => $xv, y => $yv, -type =>'polygon');
#	$vectorout->polygon (%$points,
#           			id=>"$idvblock-frame",
#						style=> { 'fill'=> 'rgb(240, 248, 255)',
#						'fill-opacity'=>1,
#						'stroke'=>'black',
#						'stroke-width'   =>  $confighash{'plot_line_size'},#$confighash{'plot_line_size'},
#						'stroke-opacity' => 1,
#						});
	$vectorout->rectangle(x => $starttingX, y => $starttingY, 
						width  	=> $confighash{'plot_width'}-$confighash{'plot_margin_right'}-$confighash{'plot_margin_right'}, 
						height => $block_total_height,
						id=> "$idvblock-frame",
						style => {'fill' => 'rgb(240, 248, 255)',
								'fill-opacity'   =>  1,
								'stroke'         => 'black',
								'stroke-width'   =>  $confighash{'plot_line_size'},
								'stroke-opacity' =>  1,
							},
						);
	my $inner1topXleft=$starttingX+$confighash{'plot_line_size'}/2;
	my $inner1topXright=$confighash{'plot_width'}-$confighash{'plot_line_size'}/2-$confighash{'plot_margin_right'};
	my $inner1topY=$starttingY+$confighash{'plot_line_size'}/2+$confighash{'plot_inner1_height'}/2;
	unless ($confighash{'plot_inner1_height'}>0) {
		$vectorout->line (   id=> "$idvblock-top-inner1",
							x1 => $inner1topXleft,
							y1 => $inner1topY, 
							x2 => $inner1topXright, 
							y2 => $inner1topY, 
							stroke => $confighash{'plot_inner1_color'}, 
							"stroke-width" => $confighash{'plot_inner1_height'}
						
						);
	}
	my $inner2topXleft=$starttingX+$confighash{'plot_line_size'}/2;
	my $inner2topXright=$confighash{'plot_width'}-$confighash{'plot_line_size'}/2-$confighash{'plot_margin_right'};
	my $inner2topY=$inner1topY+$confighash{'plot_inner1_height'}/2+$confighash{'plot_inner2_height'}/2;
	unless ($confighash{'plot_inner2_height'}>0) {
	$vectorout->line (   id=> "$idvblock-top-inner2",
						x1 => $inner2topXleft,
						y1 => $inner2topY, 
						x2 => $inner2topXright, 
						y2 => $inner2topY, 
						stroke => $confighash{'plot_inner2_color'}, 
						"stroke-width" => $confighash{'plot_inner2_height'}
					);
	}
	my $inner1bottomXleft=$starttingX+$confighash{'plot_line_size'}/2;
	my $inner1bottomXright=$confighash{'plot_width'}-$confighash{'plot_line_size'}/2-$confighash{'plot_margin_right'};
	my $inner1bottomY=$starttingY+$block_total_height-$confighash{'plot_line_size'}/2-$confighash{'plot_inner1_height'}/2;
	$vectorout->line (   id=> "$idvblock-bottom-inner1",
						x1 => $inner1bottomXleft,
						y1 => $inner1bottomY, 
						x2 => $inner1bottomXright, 
						y2 => $inner1bottomY,
						stroke => $confighash{'plot_inner1_color'}, 
						"stroke-width" => $confighash{'plot_inner1_height'}
					);
	my $inner2bottomXleft=$starttingX+$confighash{'plot_line_size'}/2;
	my $inner2bottomXright=$confighash{'plot_width'}-$confighash{'plot_line_size'}/2-$confighash{'plot_margin_right'};
	my $inner2bottomY=$inner1bottomY-$confighash{'plot_inner1_height'}/2-$confighash{'plot_inner2_height'}/2;
	$vectorout->line (   id=> "$idvblock-bottom-inner2",
						x1 => $inner2bottomXleft,
						y1 => $inner2bottomY, 
						x2 => $inner2bottomXright, 
						y2 => $inner2bottomY,
						stroke => $confighash{'plot_inner2_color'}, 
						"stroke-width" => $confighash{'plot_inner2_height'}
					);
###Drawing tiks and scale
	foreach (my $x=0; $x<$confighash{"ruler_${idvblock}_length"}; $x+=$confighash{'ruler_mark_interval'}) {
		my $scale_textX=$x*$horizontal_drawing_space/$confighash{"ruler_${idvblock}_length"}+$confighash{'plot_line_size'}+$confighash{'plot_margin_left'};
		if ($x==0) {
			my $scalar_block_Y1=$starttingY+$block_total_height-$confighash{'plot_line_size'}/2;
			$vectorout->text(x => $scale_textX+$confighash{'ruler_font_distance'}, 
						y => $scalar_block_Y1-$confighash{'ruler_font_distance'}, 
						width => $confighash{'ruler_font_size'}, 
						height => $confighash{'ruler_font_size'}, 
						"font-family"=>$confighash{'plot_font_family'}, 
						"text-anchor"=>"start",
						"font-size"=>$confighash{'ruler_font_size'}, 
						"-cdata" => "${idvblock}");
			next;
		}
		
		my $scalar_top_Y1=$starttingY+$confighash{'plot_line_size'}/2;
		my $scalar_top_Y2=$scalar_top_Y1+$confighash{'ruler_tick_length'}/2;
		$vectorout->line (id=> "$idvblock-top-tick-$x",
						x1 => $scale_textX,
						y1 => $scalar_top_Y1, 
						x2 => $scale_textX, 
						y2 => $scalar_top_Y2,
						stroke => $confighash{'ruler_tick_color'}, 
						"stroke-width" => $confighash{'ruler_tick_width'}
					);
		my $scalar_bottom_Y1=$starttingY+$block_total_height-$confighash{'plot_line_size'}/2;
		my $scalar_bottom_Y2=$scalar_bottom_Y1-$confighash{'ruler_tick_length'};
		$vectorout->line (id=> "$idvblock-bottom-tick-$x",
						x1 => $scale_textX,
						y1 => $scalar_bottom_Y1, 
						x2 => $scale_textX, 
						y2 => $scalar_bottom_Y2,
						stroke => $confighash{'ruler_tick_color'}, 
						"stroke-width" => $confighash{'ruler_tick_width'}
					);
		my $scale_text=$x;
		$scale_text=~s/000000000$/G/;$scale_text=~s/000000$/M/;$scale_text=~s/000$/K/;
		$vectorout->text(	x => $scale_textX+$confighash{'ruler_font_distance'}, 
						y => $scalar_bottom_Y1-$confighash{'ruler_font_distance'}, 
						width => $confighash{'ruler_font_size'}, 
						height => $confighash{'ruler_font_size'}, 
						"font-family"=>$confighash{'plot_font_family'}, 
						"text-anchor"=>"start",
						"font-size"=>$confighash{'ruler_font_size'}, 
						"-cdata" => "$scale_text");
	}
###Drawing tracks
	my $gridline_drawed=0;
	my $non_expression_startY=$starttingY+$block_total_height-$confighash{'plot_line_size'}/2-$confighash{'ruler_tick_length'};
	for (my $j=(scalar(@trackoders)-1); $j>=0; $j--) {
		my $idvtrack=$trackoders[$j];
		next unless ($alltrackheight{$idvblock}{$idvtrack}{'subtracks'}>0);
		
		$non_expression_startY=$non_expression_startY-$confighash{'track_space'};
#		my $svgtag=$vectorout->tag("TAG$idvblock-$idvtrack", id=>"TAG$idvblock-$idvtrack");
#		my $svganchor=$svgtag->anchor (id=>"Anchor$idvblock-$idvtrack");
		my $gridlineX1=$confighash{'plot_margin_left'}+$confighash{'plot_line_size'}/2+$confighash{'gridline_left_margin'};
		my $gridlineX2=$confighash{'plot_width'}-$confighash{'plot_line_size'}/2-$confighash{'plot_margin_right'}-$confighash{'gridline_right_margin'};
		
		if ($idvtrack eq 'EXPRESSION') {### For Expression ###
			###Drawing gridline
			if (($confighash{'bar_fill_width'}*scalar(@tissueorder)*(scalar(keys %{$annotationhash{$idvblock}{$idvtrack}})))>($gridlineX2-$gridlineX1)) {
			die "Error: not enought space to draw barplot, consider to reduce bar_fill_width BLOCK $idvblock TRACK $idvtrack\n";
		}
			my $barplot_inter_feature_space=($gridlineX2-$gridlineX1-$confighash{'bar_fill_width'}*scalar(@tissueorder)*(scalar(keys %{$annotationhash{$idvblock}{$idvtrack}})))/(scalar(keys %{$annotationhash{$idvblock}{$idvtrack}})+1);
			unless ($gridline_drawed) {
				unless ($gridlineX1<$gridlineX2) {
					die "Error: no space to draw gridline BLOCK TRACK \n";
				}
				for (my $y=$confighash{'gridline_min'}; $y<=$confighash{'gridline_max'}; $y+=$confighash{'gridline_interval'}) {
					my $gridlineY1=$non_expression_startY-$y*$alltrackheight{$idvblock}{$idvtrack}{'height'}/($confighash{'gridline_max'}-$confighash{'gridline_min'});
					$vectorout->line (id=> "$idvblock-$idvtrack$y",
						x1 => $gridlineX1,
						y1 => $gridlineY1, 
						x2 => $gridlineX2, 
						y2 => $gridlineY1,
						stroke => $confighash{'gridline_color'}, 
						"stroke-width" => $confighash{'gridline_width'}
					);
					### drawing scale
					$vectorout->text(	x => $confighash{'plot_margin_left'}+$confighash{'plot_line_size'}/2+$confighash{'gridline_left_margin'}/4, 
										y => $gridlineY1+$confighash{'gridtext_size'}/4, 
										width => $confighash{'gridline_left_margin'}/2, 
										height => $confighash{'gridline_left_margin'}/2, 
										"font-family"=>$confighash{'plot_font_family'}, 
										"text-anchor"=>"start",
										"text-color"=>$confighash{'gridtext_color'}, 
										"font-size"=>$confighash{'gridtext_size'}, 
										"-cdata" => "$y");
				}
				$gridline_drawed++;
			}
			### Draw barplot
			my %thistrackfeatures=();
			foreach my $featureid (keys %{$annotationhash{$idvblock}{$idvtrack}}) {
				next unless (exists $express_location{$featureid});
				$thistrackfeatures{$featureid}=$express_location{$featureid};
			}
			if (0) {### for test ###
				print "Test: \%thistrackfeatures\n";
				print Dumper \%thistrackfeatures;
				print "\n";
			}
			my $thisfeaturearr=[];
			$thisfeaturearr=&OrderFeatures(\%thistrackfeatures);
			%thistrackfeatures=();
			if (0) {### for test ###
				print "Test: \$thisfeaturearr\n";
				print Dumper $thisfeaturearr;
				print "\n";
			}
			### transform gene location from nucleotide to image
			$primaryloci=[];
			for (my $e=0; $e<scalar(@{$thisfeaturearr}); $e++) {
				my $left=$confighash{'plot_margin_left'}+$confighash{'plot_line_size'}/2+${$thisfeaturearr}[$e][1]*$horizontal_drawing_space/$confighash{"ruler_${idvblock}_length"};
				my $right=$confighash{'plot_margin_left'}+$confighash{'plot_line_size'}/2+${$thisfeaturearr}[$e][2]*$horizontal_drawing_space/$confighash{"ruler_${idvblock}_length"};
				my $middle=$confighash{'plot_margin_left'}+$confighash{'plot_line_size'}/2+((${$thisfeaturearr}[$e][1]+${$thisfeaturearr}[$e][2])/2)*$horizontal_drawing_space/$confighash{"ruler_${idvblock}_length"};
#				print "Test: LEFT ", ${$thisfeaturearr}[$e][1], " RIGHT ", ${$thisfeaturearr}[$e][2], " MIDDLE ", ((${$thisfeaturearr}[$e][1]+${$thisfeaturearr}[$e][2])/2), " for $idvblock TRACK $idvtrack FEATURE ", ${$thisfeaturearr}[$e][0], "\n"; ### For test ###
#				print "Test: LEFT $left RIGHT $right MIDDLE $middle for $idvblock TRACK $idvtrack FEATURE ", ${$thisfeaturearr}[$e][0], "\n"; ### For test ###
				if ($confighash{'bar_alignment'} =~/^middle$/i) {
					${$primaryloci}[$e]=[${$thisfeaturearr}[$e][0], $middle-scalar(@tissueorder)*$confighash{'bar_fill_width'}/2,  $middle+scalar(@tissueorder)*$confighash{'bar_fill_width'}/2];
				}
				elsif ($confighash{'bar_alignment'} =~/^left$/i) {
					${$primaryloci}[$e]=[${$thisfeaturearr}[$e][0], $left,  $left+scalar(@tissueorder)*$confighash{'bar_fill_width'} ];
				}
				elsif ($confighash{'bar_alignment'} =~/^right$/i) {
					${$primaryloci}[$e]=[${$thisfeaturearr}[$e][0], $right-scalar(@tissueorder)*$confighash{'bar_fill_width'} ];
				}
				else {
					die "Error: unknown bar_alignment parameter in user.config\n";
				}
			}
			$thisfeaturearr=[];
			if (0) {### for test ###
				print "Test: before ReArrangeLoci \$primaryloci\n";
				print Dumper $primaryloci;
				print "\n";
			}
			my $improvedloci=&ReArrangeLoci($gridlineX1+$confighash{'bar_fill_width'},$confighash{'bar_fill_width'}, $gridlineX2-$confighash{'bar_fill_width'}, 10000);
			if (0) {### for test ###
				print "Test: After  ReArrangeLoci \$primaryloci\n";
				print Dumper $primaryloci;
				print "\n";
			}
			for (my $f=0; $f<scalar(@{$primaryloci}); $f++) {
				my $thatfeatureid=${$primaryloci}[$f][0];
#				die "Error: ReArrangeError\n" unless ((${$primaryloci}[$f][2]-${$primaryloci}[$f][1]) == scalar(@tissueorder)*$confighash{'bar_fill_width'});
				die "Error: ReArrangeError2\n" unless (exists $annotationhash{$idvblock}{$idvtrack}{$thatfeatureid});
				my $thatblock=$loci{$thatfeatureid}{'ref'};
				my $thattrack=$loci{$thatfeatureid}{'track'};
				my $barstartX=${$primaryloci}[$f][1];
				for my $idv_tissue (@tissueorder) {
					die "Error: Expression value error BLOCK $thatblock TRACK $thattrack FEATURE_ID $thatfeatureid TISSUE $idv_tissue\n" unless (exists $annotationhash{$idvblock} and exists $annotationhash{$idvblock}{$idvtrack} and exists $annotationhash{$idvblock}{$idvtrack}{$thatfeatureid} and exists $annotationhash{$idvblock}{$idvtrack}{$thatfeatureid}{$idv_tissue} and exists $annotationhash{$idvblock}{$idvtrack}{$thatfeatureid}{$idv_tissue}{'value'} and $annotationhash{$idvblock}{$idvtrack}{$thatfeatureid}{$idv_tissue}{'value'}=~/^\d+\.*\d*$/);
					my $thisbarheight=$annotationhash{$idvblock}{$idvtrack}{$thatfeatureid}{$idv_tissue}{'value'}*$alltrackheight{$idvblock}{$idvtrack}{'height'}/($confighash{'gridline_max'}-$confighash{'gridline_min'});
					### Draw Bar
					unless ($annotationhash{$idvblock}{$idvtrack}{$thatfeatureid}{$idv_tissue}{'value'}==0) {
						$vectorout->rectangle(x => $barstartX, y => $non_expression_startY-$thisbarheight, 
							width  	=> $confighash{'bar_fill_width'}, 
							height => $thisbarheight,
							id=> "$thatblock-$thattrack-$thatfeatureid-$idv_tissue",
							style => {'fill' => $annotationhash{$idvblock}{$idvtrack}{$thatfeatureid}{$idv_tissue}{'color'},
									'fill-opacity'   =>  1,
									'stroke'         => 'black',
									'stroke-width'   =>  $confighash{'bar_line_width'},
									'stroke-opacity' =>  1,
								},
							);
					}
					### Draw Error
					unless (exists $annotationhash{$idvblock}{$idvtrack}{$thatfeatureid}{$idv_tissue}{'stdev'} and $annotationhash{$idvblock}{$idvtrack}{$thatfeatureid}{$idv_tissue}{'stdev'}==0) {
						my $thiserrheight=$annotationhash{$idvblock}{$idvtrack}{$thatfeatureid}{$idv_tissue}{'stdev'}*$alltrackheight{$idvblock}{$idvtrack}{'height'}/($confighash{'gridline_max'}-$confighash{'gridline_min'});
						my $errorX1=$barstartX+$confighash{'bar_fill_width'}/2+$confighash{'error_line_length'}/2;
						my $errorX2=$barstartX+$confighash{'bar_fill_width'}/2-$confighash{'error_line_length'}/2;
						my $errorY1=$non_expression_startY-$thisbarheight-$thiserrheight;
						my $errorY2=$non_expression_startY-$thisbarheight+$thiserrheight;
						$vectorout->line (id=> "$thatblock-$thattrack-$thatfeatureid-$idv_tissue-toperr",
										x1 => $errorX1,
										y1 => $errorY1, 
										x2 => $errorX2, 
										y2 => $errorY1,
										stroke => $confighash{'error_line_color'}, 
										"stroke-width" => $confighash{'ruler_tick_width'}
									);
						$vectorout->line (id=> "$thatblock-$thattrack-$thatfeatureid-$idv_tissue-miderr",
										x1 => ($errorX1+$errorX2)/2,
										y1 => $errorY1, 
										x2 => ($errorX1+$errorX2)/2, 
										y2 => $errorY2,
										stroke => $confighash{'error_line_color'}, 
										"stroke-width" => $confighash{'ruler_tick_width'}
									);
						$vectorout->line (id=> "$thatblock-$thattrack-$thatfeatureid-$idv_tissue-boterr",
										x1 => $errorX1,
										y1 => $errorY2, 
										x2 => $errorX2, 
										y2 => $errorY2,
										stroke => $confighash{'error_line_color'}, 
										"stroke-width" => $confighash{'ruler_tick_width'}
									);
					}
					$barstartX+=$confighash{'bar_fill_width'};
				}
			}
			$primaryloci=[];

			$non_expression_startY=$non_expression_startY-$confighash{'track_subtrack_space_height'}*($alltrackheight{$idvblock}{$idvtrack}{'subtracks'}-1)-$alltrackheight{$idvblock}{$idvtrack}{'subtracks'}*$alltrackheight{$idvblock}{$idvtrack}{'height'};
		}
### Draw non EXPRESSION tracks
		else {
			$non_expression_startY=$non_expression_startY-$confighash{'track_subtrack_space_height'}*($alltrackheight{$idvblock}{$idvtrack}{'subtracks'}-1)-$alltrackheight{$idvblock}{$idvtrack}{'subtracks'}*$alltrackheight{$idvblock}{$idvtrack}{'height'};
			foreach my $idvid (sort keys %{$annotationhash{$idvblock}{$idvtrack}}) {
				die "Error: invalid positions $annotationhash{$idvblock}{$idvtrack}{$idvid}{'position'} for BLOCK $idvblock TRACK $idvtrack ID $idvid\n" unless ($annotationhash{$idvblock}{$idvtrack}{$idvid}{'position'} =~ /^\d+-\d+$/);
				my ($this_trackX1, $this_trackX2)=split(/-/, $annotationhash{$idvblock}{$idvtrack}{$idvid}{'position'});
				my $thissubtrack=1;
				if (exists $annotationhash{$idvblock}{$idvtrack}{$idvid}{'subtrack'} and $annotationhash{$idvblock}{$idvtrack}{$idvid}{'subtrack'}=~/^\d+$/) {
					$thissubtrack=$annotationhash{$idvblock}{$idvtrack}{$idvid}{'subtrack'};
				}
				$this_trackX1=$this_trackX1*$horizontal_drawing_space/$confighash{"ruler_${idvblock}_length"}+$confighash{'plot_line_size'}+$confighash{'plot_margin_left'};
				$this_trackX2=$this_trackX2*$horizontal_drawing_space/$confighash{"ruler_${idvblock}_length"}+$confighash{'plot_line_size'}+$confighash{'plot_margin_left'};
				my $this_feature_width=$this_trackX2-$this_trackX1;
				my $this_trackY1=$non_expression_startY+($confighash{'track_subtrack_space_height'}+$alltrackheight{$idvblock}{$idvtrack}{'height'})*($thissubtrack-1);
				my $this_trackY2=$this_trackY1+$alltrackheight{$idvblock}{$idvtrack}{'height'};
				unless (exists $annotationhash{$idvblock}{$idvtrack}{$idvid}{'plot'}) {
					die "Error: unknown plot type for BLOCK $idvblock TRACK $idvtrack ID $idvid\n";
				}
				unless (exists $annotationhash{$idvblock}{$idvtrack}{$idvid}{'color'}) {
					die "Error: unknown plot color for BLOCK $idvblock TRACK $idvtrack ID $idvid\n";
				}
				
				if ($annotationhash{$idvblock}{$idvtrack}{$idvid}{'plot'} =~/^box$/) {
					$vectorout->rectangle(x => $this_trackX1, 
						y => $this_trackY1, 
						width  	=> $this_feature_width, 
						height => $alltrackheight{$idvblock}{$idvtrack}{'height'}, 
						id=> "$idvblock-$idvtrack-$idvid",
						style => {'fill' => $annotationhash{$idvblock}{$idvtrack}{$idvid}{'color'},
								'fill-opacity'   =>  1,
								'stroke'         => 'black',
								'stroke-width'   =>  $confighash{'plot_line_size'},
								'stroke-opacity' =>  1,
							},
					);
				}
				elsif ($annotationhash{$idvblock}{$idvtrack}{$idvid}{'plot'} =~/^arrow$/) {
					unless (exists $annotationhash{$idvblock}{$idvtrack}{$idvid}{'strand'} and $annotationhash{$idvblock}{$idvtrack}{$idvid}{'strand'}=~/^[+-]{1}$/) {
						die "Error: unknown strand of plot type 'ARROW' for BLOCK $idvblock TRACK $idvtrack ID $idvid\n";
					}
					$xv=[]; $yv=[];
#					print "Test: BLOCK $idvblock TRACK $idvtrack ID $idvid STRAND $annotationhash{$idvblock}{$idvtrack}{$idvid}{'strand'}\n"; ### for Test ###
#					print "Test: FEATURE width $this_feature_width ARROW $confighash{'track_arrow_length'} HEIGHT $alltrackheight{$idvblock}{$idvtrack}{'height'}\n"; ### For test ###
					if ($annotationhash{$idvblock}{$idvtrack}{$idvid}{'strand'} eq '+') {
						if ($this_feature_width<=$confighash{'track_arrow_length'} or $alltrackheight{$idvblock}{$idvtrack}{'height'}<=$confighash{'track_arrow_add'}*2) {
							$xv=[$this_trackX1, $this_trackX1, $this_trackX2];
							$yv=[$this_trackY1, $this_trackY2, ($this_trackY1+$this_trackY2)/2];
						}
						else {
							$xv=[$this_trackX1, $this_trackX2-$confighash{'track_arrow_length'}, $this_trackX2-$confighash{'track_arrow_length'}, $this_trackX2, $this_trackX2-$confighash{'track_arrow_length'}, $this_trackX2-$confighash{'track_arrow_length'}, $this_trackX1];
							$yv=[$this_trackY1+$confighash{'track_arrow_add'}, $this_trackY1+$confighash{'track_arrow_add'}, $this_trackY1, ($this_trackY1+$this_trackY2)/2, $this_trackY2, $this_trackY2-$confighash{'track_arrow_add'}, $this_trackY2-$confighash{'track_arrow_add'}];
						}
					}
					elsif ($annotationhash{$idvblock}{$idvtrack}{$idvid}{'strand'} eq '-') {
						if ($this_feature_width<=$confighash{'track_arrow_length'} or $alltrackheight{$idvblock}{$idvtrack}{'height'}<=$confighash{'track_arrow_add'}*2) {
							$xv=[$this_trackX1, $this_trackX2, $this_trackX2];
							$yv=[($this_trackY1+$this_trackY2)/2, $this_trackY1, $this_trackY2];
						}
						else {
							$xv=[$this_trackX1, $this_trackX1+$confighash{'track_arrow_length'}, $this_trackX1+$confighash{'track_arrow_length'}, $this_trackX2, $this_trackX2, $this_trackX1+$confighash{'track_arrow_length'}, $this_trackX1+$confighash{'track_arrow_length'}];
							$yv=[($this_trackY1+$this_trackY2)/2, $this_trackY1, $this_trackY1+$confighash{'track_arrow_add'}, $this_trackY1+$confighash{'track_arrow_add'}, $this_trackY2-$confighash{'track_arrow_add'}, $this_trackY2-$confighash{'track_arrow_add'}, $this_trackY2];
						}
					}
					else {
						die "Error: unknown strand for BLOCK $idvblock TRACK $idvtrack ID $idvid\n";
					}
					if (0) {
						print "Test: BLOCK $idvblock TRACK $idvtrack ID $idvid: \n";
						print Dumper $xv;
						print Dumper $yv;
						print "\n\n\n";
					}
					$points={};
					$points=$vectorout->get_path( x => $xv, y => $yv, -type =>'polygon');
					$vectorout->polygon (%$points,
#								id=>"$idvblock-$idvtrack-$idvid",
								style=> { 'fill'=> $annotationhash{$idvblock}{$idvtrack}{$idvid}{'color'},
								'fill-opacity'=>1,
								'stroke'=>'black',
								'stroke-width'   =>  1,#$confighash{'plot_line_size'},
								'stroke-opacity' => 0,
										});
				}
			}
		}
	}
	$starttingY+=$block_total_height+$confighash{'synteny_height'}+$confighash{'plot_line_size'};
}



### Draw synteny
#print "Test: @blockorder\n"; ### For test ###
#print "Test: \%syntenyhash\n"; print Dumper \%syntenyhash; print "\n"; ### For test ###
#print "Test: \%block_coordinates\n"; print Dumper \%block_coordinates; print "\n"; ### For test ###
foreach (my $g=0; $g<scalar(@blockorder); $g++) {
	for (my $h=$g; $h<scalar(@blockorder); $h++) {
		my $block1=$blockorder[$g];
		my $block2=$blockorder[$h];
		unless (exists $syntenyhash{"$block1-$block2"}) {
			next;
		}
		
#		print "Test: Drawing Synteny $block1=>$block2\n";### For test ###
		my $j=0;
		foreach my $i (@{$syntenyhash{"$block1-$block2"}}) {
#			print "Test: \$i\n"; print Dumper $i; print "\n"; ### For test ###
			my $block1_X1=$confighash{'plot_margin_left'}+$confighash{'plot_line_size'}+$i->[0] * $horizontal_drawing_space/$confighash{"ruler_${block1}_length"};
			my $block1_X2=$confighash{'plot_margin_left'}+$confighash{'plot_line_size'}+$i->[1] * $horizontal_drawing_space/$confighash{"ruler_${block1}_length"};
			my $block2_X1=$confighash{'plot_margin_left'}+$confighash{'plot_line_size'}+$i->[2] * $horizontal_drawing_space/$confighash{"ruler_${block2}_length"};
			my $block2_X2=$confighash{'plot_margin_left'}+$confighash{'plot_line_size'}+$i->[3] * $horizontal_drawing_space/$confighash{"ruler_${block2}_length"};
			my $block1_Y=$block_coordinates{$block1}{'bottom'}+$confighash{'plot_line_size'}/2;
			my $block2_Y=$block_coordinates{$block2}{'top'}-$confighash{'plot_line_size'}/2;
			my $xv=[]; my $yv=[];
			my $points={};
			$xv=[ $block1_X1, $block1_X2, $block2_X2, $block2_X1 ];
			$yv=[ $block1_Y, $block1_Y, $block2_Y, $block2_Y ];
			$points=$vectorout->get_path( x => $xv, y => $yv, -type =>'polygon');
			if (0) { ### For test ###
				print "Test: \$points\n";
				print Dumper $points;
				print "\n";
			}
			my $this_fillincolor;
			if (exists $confighash{'synteny_fillin_color'}) {
				$this_fillincolor=$confighash{'synteny_fillin_color'};
			}
			else {
				$this_fillincolor='rgb('.int(rand(256)).', '.int(rand(256)).', '.int(rand(256)).')';
			}
			$vectorout->polygon (%$points,
						id=>"$block1-$block2-synteny-$j",
						style=> { 'fill'=> $this_fillincolor,
						'fill-opacity'=>1,
						'stroke'=>'black',
						'stroke-width'   =>  0,
						'stroke-opacity' => 0,
						});
		$j++;
		}
		$j=0;
		foreach my $i (@{$syntenyhash{"$block1-$block2"}}) {
#			print "Test: \$i\n"; print Dumper $i; print "\n"; ### For test ###
			my $block1_X1=$confighash{'plot_margin_left'}+$confighash{'plot_line_size'}+$i->[0] * $horizontal_drawing_space/$confighash{"ruler_${block1}_length"};
			my $block1_X2=$confighash{'plot_margin_left'}+$confighash{'plot_line_size'}+$i->[1] * $horizontal_drawing_space/$confighash{"ruler_${block1}_length"};
			my $block2_X1=$confighash{'plot_margin_left'}+$confighash{'plot_line_size'}+$i->[2] * $horizontal_drawing_space/$confighash{"ruler_${block2}_length"};
			my $block2_X2=$confighash{'plot_margin_left'}+$confighash{'plot_line_size'}+$i->[3] * $horizontal_drawing_space/$confighash{"ruler_${block2}_length"};
			my $block1_Y=$block_coordinates{$block1}{'bottom'}+$confighash{'plot_line_size'}/2;
			my $block2_Y=$block_coordinates{$block2}{'top'}-$confighash{'plot_line_size'}/2;
			$vectorout->line (id=> "$block1-$block2-synteny-$j-left",
								x1 => $block1_X1,
								y1 => $block1_Y, 
								x2 => $block2_X1, 
								y2 => $block2_Y,
								stroke => $confighash{'synteny_line_color'}, 
								"stroke-width" => $confighash{'synteny_line_width'}
							);
			$vectorout->line (id=> "$block1-$block2-synteny-$j-right",
								x1 => $block1_X2,
								y1 => $block1_Y, 
								x2 => $block2_X2, 
								y2 => $block2_Y,
								stroke => $confighash{'synteny_line_color'}, 
								"stroke-width" => $confighash{'synteny_line_width'}
							);
		$j++;
		}
	}
}



### Drawing plot legend
my ($plotYmin, $plotYmax)=('', '');
foreach (keys %block_coordinates) {
	if ($plotYmin eq '' or $plotYmin>$block_coordinates{$_}{'top'}) {$plotYmin=$block_coordinates{$_}{'top'};next;}
	if ($plotYmax eq '' or $plotYmax<$block_coordinates{$_}{'bottom'}) {$plotYmax=$block_coordinates{$_}{'bottom'};next;}
}
print "Test: Plot MIN $plotYmin MAX $plotYmin\n";
my $legend_textY=0;
if ($confighash{'legend_plot_position'} =~/^TOP$/i) {
	$legend_textY=$plotYmin-$confighash{'plot_line_size'};
	if ($legend_textY<$confighash{'legend_text_size'}) {
		die "Error: no Y space to draw legend\n\tConsider to reduce legend_text_size or increase plot_margin_top in user.config\n";
	}
}
elsif ($confighash{'legend_plot_position'} =~/^BOTTOM$/i) {
	$legend_textY=$plotYmax+$confighash{'plot_line_size'}+$confighash{'legend_text_size'};
	if ($legend_textY> $confighash{'plot_height'}) {
		die "Error: no Y space to draw legend\n\tConsider to reduce legend_text_size or increase plot_margin_bottom in user.config\n";
	}
}
else {
	die "Error: invalid legend_plot_position parameter in user.config\n";
}
my $tissueXinterval=$confighash{'legend_between_interval'}*$confighash{'legend_text_size'};
my $total_lengthX=scalar(@tissueorder)*$tissueXinterval;
if ($total_lengthX>($confighash{'plot_width'}-$confighash{'plot_margin_left'}-$confighash{'plot_margin_right'})) {
	die "Error: no X space to draw legend\n\tConsider to reduce plot_margin_left/plot_margin_right/legend_text_size or increase plot_width in user.config\n";
}
my $legend_startX=0;
if ($confighash{'legend_plot_aplignment'}=~/^left$/i) {
	$legend_startX=$confighash{'plot_margin_left'};
}
elsif ($confighash{'legend_plot_aplignment'}=~/^middle$/i) {
	$legend_startX=$confighash{'plot_width'}/2-$total_lengthX/2;
}
elsif ($confighash{'legend_plot_aplignment'}=~/^right$/i) {
	$legend_startX=$confighash{'plot_width'}-$confighash{'plot_margin_right'}-$total_lengthX;
}
else {
	die "Error: invalid legend_plot_aplignment parameter in user.config\n";
}

foreach my $tissue (@tissueorder) {
	die "Error: $tissue no color\n" unless (exists $tissue_lengend{$tissue});
	my @tissuecolors=keys %{$tissue_lengend{$tissue}};
	if (scalar(@tissuecolors)<1) {
		die "Error: tissue $tissue have no color\n";
	}
	elsif (scalar(@tissuecolors)>1) {
		die "Error: tissue $tissue have multiple color: @tissuecolors\n" 
	}
	my $thistissuecolor=$tissuecolors[0];
	
	
	$vectorout->rectangle(x => $legend_startX, y => $legend_textY-$confighash{'legend_text_size'}, 
						width  	=> $confighash{'legend_text_size'}, 
						height => $confighash{'legend_text_size'},
						id=> "LEGEND-$tissue-rectangular",
						style => {'fill' => $thistissuecolor,
								'fill-opacity'   =>  1,
								'stroke'         => 'black',
								'stroke-width'   =>  1,
								'stroke-opacity' =>  1,
							},
						);
	$vectorout->text(	x => $legend_startX+$confighash{'plot_line_size'}+$confighash{'legend_text_size'},
						y => $legend_textY, 
						width => $confighash{'legend_text_size'}*$confighash{'legend_between_interval'}, 
						height => $confighash{'legend_text_size'}, 
						"font-family"=>$confighash{'plot_font_family'}, 
						"text-anchor"=>"start",
						"text-color"=>$confighash{'legend_text_color'}, 
						"font-size"=>$confighash{'legend_text_size'}, 
						"-cdata" => "$tissue");

	$legend_startX+=$tissueXinterval;
}


my $finalout = $vectorout->xmlify;
open SVGFILE, "> $outputsvg";
print SVGFILE $finalout;
close SVGFILE;



#####################################################################
######################### Sub Functions #############################
#####################################################################



### Load user.config parameters
### readUsrCfg ($user_configure_file)
### Global: %confighash
### Dependencies:
### Note:
sub readUsrCfg {
	my $RUCfile=shift;

	my $RUCsubinfo="(readUsrCfg)";
	local *RUCCONFIG;
	my $RUCnumline=0;
	
	unless (defined $RUCfile and -s $RUCfile) {
		print STDERR $RUCsubinfo, "Error: invalid user configure file\n";
		return 0;
	}
	close RUCCONFIG if (defined fileno(RUCCONFIG));
	unless (open (RUCCONFIG, "< $RUCfile")) {
		print STDERR $RUCsubinfo, "Error: can not open user.config file: $RUCfile\n";
		return 0;
	}
	
	while (my $RUCline=<RUCCONFIG>) {
		chomp $RUCline;
		$RUCnumline++;
		next if ($RUCline=~/(^#)|(^\[)|(^\s*$)/);
		my @RUCarr=split(/=/,$RUCline);
		unless (scalar(@RUCarr)==2) {### only one = allowed perl line
			print STDERR $RUCsubinfo, "Error: invalid user.config line($RUCnumline): more than two '=': $RUCline\n";
			return 0;
		}
		#remove spaces before and after strings
		$RUCarr[0]=~s/^\s+//;$RUCarr[0]=~s/\s+$//;
		$RUCarr[1]=~s/^\s+//;$RUCarr[1]=~s/\s+$//;
		$confighash{$RUCarr[0]}=$RUCarr[1];
	}
	close RUCCONFIG;


### check default [plot]
	unless (exists $confighash{'plot_height'} and $confighash{'plot_height'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default plot_height=", $userdefault{'plot_height'}, "\n";
		$confighash{'plot_height'}=$userdefault{"plot_height"};
	}
	unless (exists $confighash{'plot_width'} and $confighash{'plot_width'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default plot_width=", $userdefault{'plot_width'}, "\n";
		$confighash{'plot_width'}={'plot_width'};
	}
	unless (exists $confighash{'plot_margin_left'} and $confighash{'plot_margin_left'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default plot_margin_left=", $userdefault{'plot_margin_left'}, "\n";
		$confighash{"plot_margin_left"}=$userdefault{'plot_margin_left'};
	}
	unless (exists $confighash{'plot_margin_right'} and $confighash{'plot_margin_right'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default plot_margin_right=", $userdefault{'plot_margin_right'}, "\n";
		$confighash{'plot_margin_right'}=$userdefault{'plot_margin_right'};
	}
	unless (exists $confighash{'plot_margin_top'} and $confighash{'plot_margin_top'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default plot_margin_top=", $userdefault{'plot_margin_top'}, "\n";
		$confighash{'plot_margin_top'}=$userdefault{'plot_margin_top'};
	}
	unless (exists $confighash{'plot_margin_bottom'} and $confighash{'plot_margin_bottom'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default plot_margin_bottom=", $userdefault{'plot_margin_bottom'}, "\n";
		$confighash{'plot_margin_bottom'}=$userdefault{'plot_margin_bottom'};
	}
	unless (exists $confighash{'plot_line_size'} and $confighash{'plot_line_size'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default plot_line_size=", $userdefault{'plot_line_size'}, "\n";
		$confighash{'plot_line_size'}=$userdefault{'plot_line_size'};
	}
	unless (exists $confighash{'plot_background_color'} and $confighash{'plot_background_color'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default plot_background_color=", $userdefault{'plot_background_color'}, "\n";
		$confighash{'plot_background_color'}=$userdefault{'plot_background_color'};
	}


### check default [synteny]
	unless (exists $confighash{'synteny_height'} and $confighash{'synteny_height'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default synteny_height=", $userdefault{'synteny_height'}, "\n";
		$confighash{"synteny_height"}=$userdefault{'synteny_height'};
	}
	unless (exists $confighash{'synteny_fillin_color'} and $confighash{'synteny_fillin_color'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default synteny_fillin_color=", $userdefault{'synteny_fillin_color'}, "\n";
		$confighash{"synteny_fillin_color"}=$userdefault{'synteny_fillin_color'};
	}
	unless (exists $confighash{'synteny_line_color'} and $confighash{'synteny_line_color'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default synteny_line_color=", $userdefault{'synteny_line_color'}, "\n";
		$confighash{"synteny_line_color"}=$userdefault{'synteny_line_color'};
	}
	unless (exists $confighash{'synteny_line_width'} and $confighash{'synteny_line_width'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default synteny_line_width=", $userdefault{'synteny_line_width'}, "\n";
		$confighash{"synteny_line_width"}=$userdefault{'synteny_line_width'};
	}
	unless (exists $confighash{'synteny_order'} and $confighash{'synteny_order'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default synteny_order=", $userdefault{'synteny_order'}, "\n";
		$confighash{'synteny_order'}=$userdefault{'synteny_order'};
	}



### check default [ruler]
	unless (exists $confighash{'ruler_mark_interval'} and $confighash{'ruler_mark_interval'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default ruler_mark_interval=", $userdefault{"ruler_mark_interval"}, "\n";
		$confighash{'ruler_mark_interval'}=$userdefault{'ruler_mark_interval'};
	}
	unless (exists $confighash{'ruler_tick_color'} and $confighash{'ruler_tick_color'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default ruler_tick_color=", $userdefault{'ruler_tick_color'}, "\n";
		$confighash{'ruler_tick_color'}=$userdefault{'ruler_tick_color'};
	}
	unless (exists $confighash{'ruler_tick_length'} and $confighash{'ruler_tick_length'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default ruler_tick_length=", $userdefault{'ruler_tick_length'}, "\n";
		$confighash{'ruler_tick_length'}=$userdefault{'ruler_tick_length'};
	}
	unless (exists $confighash{'ruler_font_size'} and $confighash{'ruler_font_size'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default ruler_font_size=", $confighash{'ruler_font_size'}, "\n";
		$confighash{'ruler_font_size'}=$userdefault{'ruler_font_size'};
	}
	unless (exists $confighash{'ruler_tick_width'} and $confighash{'ruler_tick_width'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default ruler_tick_width=", $userdefault{'ruler_tick_width'}, "\n";
		$confighash{'ruler_tick_width'}=$userdefault{'ruler_tick_width'};
	}
	unless (exists $confighash{'ruler_font_distance'} and $confighash{'ruler_font_distance'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, 'Warnings: use default ruler_font_distance=', $userdefault{'ruler_font_distance'}, "\n";
		$confighash{'ruler_font_distance'}=$userdefault{'ruler_font_distance'};
	}
	unless (exists $confighash{'ruler_text_color'} and $confighash{'ruler_text_color'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default ruler_text_color=", $userdefault{'ruler_text_color'}, "\n";
		$confighash{'ruler_text_color'}=$userdefault{'ruler_text_color'};
	}



### check default [track]
	unless (exists $confighash{'track_height'} and $confighash{'track_height'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default track_height=", $userdefault{'track_height'}, "\n";
		$confighash{'track_height'}=$userdefault{'track_height'};
	}
	unless (exists $confighash{'track_space'} and $confighash{'track_space'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default track_space=", $userdefault{'track_space'}, "\n";
		$confighash{'track_space'}=$userdefault{'track_space'};
	}
	unless (exists $confighash{'track_color'} and $confighash{'track_color'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default track_color=", $userdefault{'track_color'}, "\n";
		$confighash{'track_color'}=$userdefault{'track_color'};
	}
	unless (exists $confighash{'track_arrow_add'} and $confighash{'track_arrow_add'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default track_arrow_add=", $userdefault{'track_arrow_add'}, "\n";
		$confighash{'track_arrow_add'}=$userdefault{'track_arrow_add'};
	}
	unless (exists $confighash{'track_arrow_length'} and $confighash{'track_arrow_length'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default track_arrow_length=", $userdefault{'track_arrow_length'}, "\n";
		$confighash{'track_arrow_length'}=$userdefault{'track_arrow_length'};
	}
	unless (exists $confighash{'track_order'} and $confighash{'track_order'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default track_order=", $userdefault{'track_order'}, "\n";
		$confighash{'track_order'}=$userdefault{'track_order'};
	}
	unless (exists $confighash{'track_subtrack_space_height'} and $confighash{'track_subtrack_space_height'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default track_subtrack_space_height=", $userdefault{'track_subtrack_space_height'}, "\n";
		$confighash{'track_subtrack_space_height'}=$userdefault{'track_subtrack_space_height'};
	}

### check default [expression]
### check default [expression][gridline]
	unless (exists $confighash{'gridline_width'} and $confighash{'gridline_width'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default gridline_width=", $userdefault{'gridline_width'}, "\n";
		$confighash{'gridline_width'}=$userdefault{'gridline_width'};
	}
	unless (exists $confighash{'gridline_color'} and $confighash{'gridline_color'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default gridline_color=", $userdefault{'gridline_color'}, "\n";
		$confighash{'gridline_color'}=$userdefault{'gridline_color'};
	}
	unless (exists $confighash{'gridline_interval'} and $confighash{'gridline_interval'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default gridline_interval=", $userdefault{'gridline_interval'}, "\n";
		$confighash{'gridline_interval'}=$userdefault{'gridline_interval'};
	}
	unless (exists $confighash{'gridline_min'} and $confighash{'gridline_min'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default gridline_min=", $userdefault{'gridline_min'}, "\n";
		$confighash{'gridline_min'}=$userdefault{'gridline_min'};
	}
	unless (exists $confighash{'gridline_max'} and $confighash{'gridline_max'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default gridline_max=", $userdefault{'gridline_max'}, "\n";
		$confighash{'gridline_max'}=$userdefault{'gridline_max'};
	}
	unless (exists $confighash{'gridline_left_margin'} and $confighash{'gridline_left_margin'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default gridline_left_margin=", $userdefault{'gridline_left_margin'}, "\n";
		$confighash{'gridline_left_margin'}=$userdefault{'gridline_left_margin'};
	}
	unless (exists $confighash{'gridline_right_margin'} and $confighash{'gridline_right_margin'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default gridline_right_margin=", $userdefault{'gridline_right_margin'}, "\n";
		$confighash{'gridline_right_margin'}=$userdefault{'gridline_right_margin'};
	}
	unless (exists $confighash{'gridtext_color'} and $confighash{'gridtext_color'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default gridtext_color=", $userdefault{'gridtext_color'}, "\n";
		$confighash{'gridtext_color'}=$userdefault{'gridtext_color'};
	}
	unless (exists $confighash{'gridtext_size'} and $confighash{'gridtext_size'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default gridtext_size=", $userdefault{'gridtext_size'}, "\n";
		$confighash{'gridtext_size'}=$userdefault{'gridtext_size'};
	}
### check default [expression][barplot]
	unless (exists $confighash{'barplot_order'} and $confighash{'barplot_order'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default barplot_order=", $userdefault{'barplot_order'}, "\n";
		$confighash{'barplot_order'}=$userdefault{'barplot_order'};
	}
	unless (exists $confighash{'bar_line_width'} and $confighash{'bar_line_width'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default bar_line_width=", $userdefault{'bar_line_width'}, "\n";
		$confighash{'bar_line_width'}=$userdefault{'bar_line_width'};
	}
	unless (exists $confighash{'bar_alignment'} and $confighash{'bar_alignment'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default bar_alignment=", $userdefault{'bar_alignment'}, "\n";
		$confighash{'bar_alignment'}=$userdefault{'bar_alignment'};
	}
	unless (exists $confighash{'bar_fill_width'} and $confighash{'bar_fill_width'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default bar_fill_width=", $userdefault{'bar_fill_width'}, "\n";
		$confighash{'bar_fill_width'}=$userdefault{'bar_fill_width'};
	}
### check default [expression][stdev]
	unless (exists $confighash{'error_line_width'} and $confighash{'error_line_width'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default error_line_width=", $userdefault{'error_line_width'}, "\n";
		$confighash{"error_line_width"}=$userdefault{'error_line_width'};
	}
	unless (exists $confighash{'error_line_length'} and $confighash{'error_line_length'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default error_line_length=", $userdefault{'error_line_length'}, "\n";
		$confighash{'error_line_length'}=$userdefault{'error_line_length'};
	}
	unless (exists $confighash{'bar_line_width'} and $confighash{'bar_line_width'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default bar_line_width=", $userdefault{'bar_line_width'}, "\n";
		$confighash{'bar_line_width'}=$userdefault{'bar_line_width'};
	}
### [expression][LEGEND]
	unless (exists $confighash{'legend_plot_position'} and $confighash{'legend_plot_position'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default legend_plot_position=", $userdefault{'legend_plot_position'}, "\n";
		$confighash{'legend_plot_position'}=$userdefault{'legend_plot_position'};
	}
	unless (exists $confighash{'legend_text_color'} and $confighash{'legend_text_color'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default legend_text_color=", $userdefault{'legend_text_color'}, "\n";
		$confighash{'legend_text_color'}=$userdefault{'legend_text_color'};
	}
	unless (exists $confighash{'legend_text_size'} and $confighash{'legend_text_size'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default legend_text_size=", $userdefault{'legend_text_size'}, "\n";
		$confighash{'legend_text_size'}=$userdefault{'legend_text_size'};
	}
	unless (exists $confighash{'legend_plot_aplignment'} and $confighash{'legend_plot_aplignment'}=~/^\S+$/) {
		print STDERR $RUCsubinfo, "Warnings: use default legend_plot_aplignment=", $userdefault{'legend_plot_aplignment'}, "\n";
		$confighash{'legend_plot_aplignment'}=$userdefault{'legend_plot_aplignment'};
	}
	unless (exists $confighash{'legend_between_interval'} and $confighash{'legend_between_interval'}=~/^\d+\.*\d*$/) {
		print STDERR $RUCsubinfo, "Warnings: use default legend_between_interval=", $userdefault{'legend_between_interval'}, "\n";
		$confighash{'legend_between_interval'}=$userdefault{'legend_between_interval'};
	}
	return 1;
}



### read synteny config
### &readSyntenyCfg($synteny_config)
### Global: %syntenyhash, @blockorder, %seqlength; $use_defined_order
### Dependencies: 
### Note:
sub readSyntenyCfg {
	my $RSCfile=shift;

	my $RSCsubinfo="(readSyntenyCfg)";
	local *RSCCONFIG;
	my $RSCnumline=0;
	my %RSCorder_hash=();

	unless (defined $RSCfile and -s $RSCfile) {
		print STDERR $RSCsubinfo, "Error: invalid synteny configure file\n";
		return 0;
	}
	close RSCCONFIG if (defined fileno(RSCCONFIG));
	unless (open (RSCCONFIG, "< $RSCfile")) {
		print STDERR $RSCsubinfo, "Error: can not open synteny.config file: $RSCfile\n";
		return 0;
	}
	while (my $RSCline=<RSCCONFIG>) {
		chomp $RSCline;
		$RSCnumline++;
		next if ($RSCline=~/(^#)|(^\s*$)/);
		my @RSCarr=split(/\t/, $RSCline);
		###check format
		unless (scalar(@RSCarr)>=6) {
			print STDERR $RSCsubinfo, "Error: invalid synteny line($RSCnumline): less than 6 column: $RSCline\n";
			return 0;
		}
		$RSCarr[1]=~s/^\s+//;$RSCarr[1]=~s/\s+$//;
		$RSCarr[2]=~s/^\s+//;$RSCarr[2]=~s/\s+$//;
		unless ($RSCarr[1]=~/^\d+$/ and $RSCarr[2]=~/^\d+$/ and $RSCarr[1] < $RSCarr[2]) {
			print STDERR $RSCsubinfo, "Error: invalid synteny line($RSCnumline): non-number: $RSCline\n";
			return 0;
		}
		$RSCarr[4]=~s/^\s+//;$RSCarr[4]=~s/\s+$//;
		$RSCarr[5]=~s/^\s+//;$RSCarr[5]=~s/\s+$//;
		unless ($RSCarr[4]=~/^\d+$/ and $RSCarr[5]=~/^\d+$/ and $RSCarr[4] < $RSCarr[5]) {
			print STDERR $RSCsubinfo, "Error: invalid synteny line($RSCnumline): non-number2: $RSCline\n";
			return 0;
		}
		### load seqlength
		if (exists $seqlength{$RSCarr[0]}) {
			$seqlength{$RSCarr[0]}=$RSCarr[2] if ($RSCarr[2]>$seqlength{$RSCarr[0]});
		}
		else {
			$seqlength{$RSCarr[0]}=$RSCarr[2];
		}
		if (exists $seqlength{$RSCarr[3]}) {
			$seqlength{$RSCarr[3]}=$RSCarr[5] if ($RSCarr[5]>$seqlength{$RSCarr[3]});
		}
		else {
			$seqlength{$RSCarr[3]}=$RSCarr[5];
		}
		###load order
		if ($use_defined_order==0) {
			push (@blockorder, $RSCarr[0]) unless (exists $RSCorder_hash{$RSCarr[0]});
			$RSCorder_hash{$RSCarr[0]}++;
			push (@blockorder, $RSCarr[3]) unless (exists $RSCorder_hash{$RSCarr[3]});
			$RSCorder_hash{$RSCarr[3]}++;
		}
		push (@{$syntenyhash{"$RSCarr[3]-$RSCarr[0]"}}, [$RSCarr[4], $RSCarr[5], $RSCarr[1], $RSCarr[2]]);
		push (@{$syntenyhash{"$RSCarr[0]-$RSCarr[3]"}}, [$RSCarr[1], $RSCarr[2], $RSCarr[4], $RSCarr[5]]);
	}
	close RSCCONFIG;
	
	return 1;
}


### read annotations config
### &readAnnotationCfg($annotation_config)
### Global: %annotationhash, @blockorder
### Dependencyies: 
### Note: 
sub readAnnotationCfg {
	my $RACfile=shift;
	
	my $RACsubinfo="(readAnnotationCfg)";
	local *RACCONFIG;
	my $RACnumline=0;
	my %RACtissueorders=();
	my %RACtracks=();
	unless (defined $RACfile and -s $RACfile) {
		print STDERR $RACsubinfo, "Error: invalid annotation configure file\n";
		return 0;
	}
	close RACCONFIG if (defined fileno(RACCONFIG));
	unless (open (RACCONFIG, "< $RACfile")) {
		print STDERR $RACsubinfo, "Error: can not open annotation.config file: $RACfile\n";
		return 0;
	}
	while (my $RACline=<RACCONFIG>) {
		chomp $RACline;
		$RACnumline++;
		next if ($RACline=~/(^#)|(^\s*$)/);
		my @RACarr=split(/\t/, $RACline);
		unless (scalar(@RACarr) ==9) {### check if none columns
			print STDERR $RACsubinfo, "Error: not 9 columns at line ($RACnumline): $RACline\n";
			return 0;
		}
		$RACarr[6]=uc($RACarr[6]);
		if ($use_defined_trackorder==0) {
			$RACtracks{$RACarr[6]}++;
		}
		if ($RACarr[6] =~/^expression$/i) {
			unless ($RACarr[7] eq 'xyplot') {
				print STDERR $RACsubinfo, "Error: invalid plot type at line ($RACnumline): expression currently only supports xyplot: $RACline\n";
				return 0;
			}
			$annotationhash{$RACarr[0]}{'EXPRESSION'}{$RACarr[1]}{$RACarr[2]}{'value'}=$RACarr[5];
			$annotationhash{$RACarr[0]}{'EXPRESSION'}{$RACarr[1]}{$RACarr[2]}{'stdev'}=$RACarr[4];
			$annotationhash{$RACarr[0]}{'EXPRESSION'}{$RACarr[1]}{$RACarr[2]}{'color'}=$RACarr[8];
			$tissue_lengend{$RACarr[2]}{$RACarr[8]}++;
			$loci{$RACarr[1]}{'ref'}=$RACarr[0];
			$loci{$RACarr[1]}{'track'}='EXPRESSION';
			if ($use_defined_tissueorder==0) {
				push (@tissueorder, $RACarr[2]) unless (exists $RACtissueorders{$RACarr[2]});
			}
			$RACtissueorders{$RACarr[2]}++;
			my %RACtemp1=();
			for (my $RACi=0; $RACi<scalar(@tissueorder); $RACi++) {
				$RACtemp1{$tissueorder[$RACi]}=$RACi;
			}
			unless (exists $RACtemp1{$RACarr[2]}) {
				print STDERR $RACsubinfo, "Error: invalid tissues at col3 at line ($RACnumline): expression currently only supports xyplot: $RACline\n";
				return 0;
			}
		}
		else {
#DD	23601	26290	+	AA0860670	.	gene	arrow	purple
#AA	11932	12009	.	.	.	repeat	box	black
			if (exists $annotationhash{$RACarr[0]} and exists $annotationhash{$RACarr[0]}{$RACarr[6]} and exists $annotationhash{$RACarr[0]}{$RACarr[6]}{$RACarr[4]}) {
				die $RACsubinfo, "Error: duplicated feature ID at line ($RACnumline): $RACarr[4]\n";
			}
			else {
				$annotationhash{$RACarr[0]}{$RACarr[6]}{$RACarr[4]}{'position'}="$RACarr[1]-$RACarr[2]";
				$annotationhash{$RACarr[0]}{$RACarr[6]}{$RACarr[4]}{'strand'}="$RACarr[3]";
				$annotationhash{$RACarr[0]}{$RACarr[6]}{$RACarr[4]}{'plot'}="$RACarr[7]";
				$annotationhash{$RACarr[0]}{$RACarr[6]}{$RACarr[4]}{'color'}="$RACarr[8]";
			}
			if (exists $express_location{$RACarr[4]}) {
				die $RACsubinfo, "Error: duplicated feature ID at line ($RACnumline): $RACarr[4]\nEach feature ID should be unique\n";
			}
			else {
				$express_location{$RACarr[4]}=[$RACarr[1], $RACarr[2]];
			}
		}
	}
	close RACCONFIG;
	@trackoders=sort (keys %RACtracks) if ($use_defined_trackorder==0);

	return 1;
}



### Assign subtrack if position overlaps
### &AssignSubtrack($blockname, $trackname)
### Global: %annotationhash %confighash
### Dependencies:
### Note:
sub AssignSubtrack {
	my ($ASblockname, $AStrackname)=@_;

	my $ASsubinfo='(AssignSubtrack)';
	my $AStracknum=0;
	my %ASpositions=();
	my %ASfinaltracks=();
	my %ASblockmax=();

	unless (exists $annotationhash{$ASblockname} and exists $annotationhash{$ASblockname}{$AStrackname}) {
		print STDERR $ASsubinfo, "Warnings: TRACK $AStrackname in BLOCK $ASblockname have no annotations\n";
		return (1, $AStracknum);
	}
	foreach my $idvid (sort keys %{$annotationhash{$ASblockname}{$AStrackname}}) {
#		print $ASsubinfo, "Test: BLOCK $ASblockname TRACK $AStrackname FeatureID $idvid\n"; ### For test ###
		if (exists $annotationhash{$ASblockname}{$AStrackname}{$idvid}{'position'}) {
			if ($annotationhash{$ASblockname}{$AStrackname}{$idvid}{'position'} =~/^(\d+)-(\d+)$/) {
				$ASpositions{$1}{$2}{$idvid}++;
			}
			else {
				print STDERR $ASsubinfo, "Warnings: invalid position ID: $idvid TRACK $AStrackname in BLOCK $ASblockname\n";
				return 0;
			}
		}
		else {
			print STDERR $ASsubinfo, "Test: no positions for BLOCK $ASblockname TRACK $AStrackname FeatureID $idvid\n"; ### For test ###
			return 0;
		}
	}
	if (0) {### For test ###
		print "\n";
		print $ASsubinfo, "Test: \%ASpositions\n";
		print Dumper \%ASpositions;
		print "\n\n\n";
	}
	foreach my $ASstart (sort {$a<=>$b} keys %ASpositions) {
		foreach my $ASend (sort {$a<=>$b} keys %{$ASpositions{$ASstart}}) {
			foreach my $ASthisid (keys %{$ASpositions{$ASstart}{$ASend}}) {
				foreach (sort {$a<=>$b} keys %ASblockmax) {
					if ($ASstart>$ASblockmax{$_}) {
						$ASfinaltracks{$ASthisid}=$_;
						$ASblockmax{$_}=$ASend;
						last;
					}
				}
				unless (exists $ASfinaltracks{$ASthisid}) {
					$AStracknum++;
					$ASblockmax{$AStracknum}=$ASend;
					$ASfinaltracks{$ASthisid}=$AStracknum;
				}
				print "Test: Feature ID $ASthisid TRACK $AStracknum\n"; ### For test ###
			}
		}
	}
	if (0) {### For test ###
		print "\n";
		print $ASsubinfo, "Test: \%ASfinaltracks\n";
		print Dumper \%ASfinaltracks;
		print "\n\n\n";
	}
	if (0) {### For test ###
		print "\n";
		print $ASsubinfo, "Test: \%ASblockmax\n";
		print Dumper \%ASblockmax;
		print "\n\n\n";
	}
	%ASpositions=();
	%ASblockmax=();
	foreach my $idvid (sort keys %{$annotationhash{$ASblockname}{$AStrackname}}) {
		if (exists $ASfinaltracks{$idvid} and $ASfinaltracks{$idvid}=~/^\d+$/) {
#			print $ASsubinfo, "Test: BLOCK $ASblockname TRACK $AStrackname FeatureID $idvid Subtrack $ASfinaltracks{$idvid}\n"; ### For test ###
			$annotationhash{$ASblockname}{$AStrackname}{$idvid}{'subtrack'}=$ASfinaltracks{$idvid};
		}
		else {
			print STDERR $ASsubinfo, "Warnings: No subtrack number for ID: $idvid TRACK $AStrackname in BLOCK $ASblockname\n";
			return 0;
		}
	}
	return (1, $AStracknum);
}



#%hash=( $feature_key1 => [$start1, $end1], $genename2 => [$start2, $end2], ...)
### Order feature by it's location, start > end> feature_key
### &OrderFeatures (\%hash)
### Global:
### Dependency:
### Note:
### Return: $arr=[[$feature_key1, $start1, $end1], ...]
sub OrderFeatures {
	my $OFhash=shift;
	
	my $OFsubinfo='(OrderFeatures)';
	my $OForderarr=[];
	my %OFtemphash=();
	my $OFd=0;
	
	foreach my $OFi (keys %{$OFhash}) {
		$OFtemphash{${$OFhash}{$OFi}[0]}{${$OFhash}{$OFi}[1]}{$OFi}++;
	}
	foreach my $OFa (sort {$a<=>$b} keys %OFtemphash) {
		foreach my $OFb (sort {$a<=>$b} keys %{$OFtemphash{$OFa}}) {
			foreach my $OFc (sort keys %{$OFtemphash{$OFa}{$OFb}}) {
				${$OForderarr}[$OFd++]=[$OFc, $OFa, $OFb];
			}
		}
	}
	
	if (0) {### For test ###
		print $OFsubinfo, "Test: \$OForderarr\n";
		print Dumper $OForderarr;
		print "\n";
	}
	
	return $OForderarr;
}



### rearrange loci if feature overlaps
### &ReArrangeLoci($primaryloci, $min, $max);
### Global: $primaryloci
### Dependency:
### Note:
sub ReArrangeLoci {
	my ($RAFmin, $RAFinterval, $RAFmax, $RAFcount)=@_;
	
	my $RAFsinfo='(ReArrangeLoci)';
	my $RAFoverlap=0;
	my $RAFtimesmove=0;
	#test overlap
	
	if (${$primaryloci}[0][1]<$RAFmin) {
		my $RAFtemp1=${$primaryloci}[0][1];
		${$primaryloci}[0][1]=$RAFmin;
		${$primaryloci}[0][2]+=$RAFmin-$RAFtemp1;
	}
	if (${$primaryloci}[scalar(@{$primaryloci})-1][2]>$RAFmax) {
		my $RAFtemp1=${$primaryloci}[scalar(@{$primaryloci})-1][2];
		${$primaryloci}[scalar(@{$primaryloci})-1][2]=$RAFmax;
		${$primaryloci}[scalar(@{$primaryloci})-1][1]=${$primaryloci}[scalar(@{$primaryloci})-1][1]-(${$primaryloci}[scalar(@{$primaryloci})-1][2]-$RAFmax);
	}
	for (my $RAFi=0; $RAFi<scalar(@{$primaryloci}); $RAFi++) {
		next if ($RAFi==0);

		if ((${$primaryloci}[$RAFi][1]-${$primaryloci}[$RAFi-1][2]) < $RAFinterval) {
			$RAFoverlap++;
			if ($RAFi==1 and (${$primaryloci}[$RAFi-1][1]-$RAFinterval/2)<$RAFmin) {
				my $RAFtemp1=${$primaryloci}[$RAFi-1][2]+$RAFinterval-${$primaryloci}[$RAFi][1];
				${$primaryloci}[$RAFi][1]+=$RAFtemp1;
				${$primaryloci}[$RAFi][2]+=$RAFtemp1;
			}
			elsif ($RAFi==(scalar(@{$primaryloci})-1) and (${$primaryloci}[$RAFi][2]+$RAFinterval/2)>$RAFmax) {
				my $RAFtemp1=$RAFinterval-(${$primaryloci}[$RAFi][1]-${$primaryloci}[$RAFi-1][2]);
				${$primaryloci}[$RAFi-1][1]-=$RAFtemp1;
				${$primaryloci}[$RAFi-1][2]-=$RAFtemp1;
			}
			else {
				my $RAFtemp1=(${$primaryloci}[$RAFi][1]+${$primaryloci}[$RAFi-1][2])/2;
				my $RAFtemp2=$RAFtemp1+$RAFinterval/2-${$primaryloci}[$RAFi][1];
				${$primaryloci}[$RAFi][1]+=$RAFtemp2;
				${$primaryloci}[$RAFi][2]+=$RAFtemp2;
				
				my $RAFtemp3=${$primaryloci}[$RAFi-1][2]-($RAFtemp1-$RAFinterval/2);
				${$primaryloci}[$RAFi-1][1]-=$RAFtemp3;
				${$primaryloci}[$RAFi-1][2]-=$RAFtemp3;
			}
		}
		else {
			next;
		}
	}
	if (--$RAFcount<1) {
		die "Error: can not ReArrange EXPRESSION locations because of overlap\n";
	}
	if ($RAFoverlap) {
		&ReArrangeLoci($RAFmin, $RAFinterval, $RAFmax, $RAFcount);
	}
	else {
		return 1;
	}
}
