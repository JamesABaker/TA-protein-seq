#!/usr/bin/perl -w

# KDcalc.pl	Simply KD calculator, window length specified at top, scale -4.5 to +4.5

# input		KDcalc_in.txt	list of fasta sequences

# output	KDcalc_out.txt	KDwin values for each sequence, plus extra o/p on minima and dist to Ct
#				keep the disorder prediction, and indicate whether a + or - charge at site
#		Also files specifically to visualise the CT

# As for the added CT files, we have further utility to focus on TA proteins.
# Limit the CT40 file to just proteins with one predicted TM.
# Algorithm for this is to predict first TM (KD >= 1.6), then block out an area of winlen2 around this.
# Then look for next KD >= 1.6 that is at least winlen2 (possibly plus a couple) away from blocked area.
# If found then increment TM counter and continue until done.  Put the TM count number in all o/p.

# With regard to the CT files, note that we have a dist_cut variable for specifying CTcut output - currently at 40,
# this matches (when add on winlen2) the value of 50 aas within CT for TA id from the 2007 Arabidopsis paper.

# Modified Oct 2014 to calculate KD average and std deviations at each location from the beginning.

if ((exists $ENV{TM_GOODBAD}) and ($ENV{TM_GOODBAD} eq "yes")) {	# allow for HCD or not in first word of data records
  $TM_goodbad		= "yes";
} else {								# default
  $TM_goodbad		= "no";
}

# ?? AS OF OCTOBER 2014, EXTRA analysis, relevant for picking out TAs I think is DEFAULT OFF - turn on as below ??

if ((exists $ENV{EXTRA_TA}) and ($ENV{EXTRA_TA} eq "yes")) {	# ?? see some errors on occasion when yes ?? BOLLOCKS
  $extra_TA		= "yes";
} else {								# default
  $extra_TA		= "no";
}

# initiliase various for mean and stddev calculation - include per site and win around site sets
$ntop		= 100000;
for ($i=1; $i<=$ntop; $i++) {
  $nx[$i]		= 0;
  $xsum[$i]		= 0;
  $xsqsum[$i]		= 0;
  $nx_win[$i]		= 0;
  $xsum_win[$i]		= 0;
  $xsqsum_win[$i]	= 0;
}

$winlen		= 5;		# window for calculation of KD nb 19 advised for KD TM (1.6 threshold) detection
$winlen2	=  2;		# = (winlen-1)/2
#$winlen	= 21;		# window for calculation of KD nb 19 advised for KD TM (1.6 threshold) detection
#$winlen2	= 10;		# = (winlen-1)/2
if ((2*$winlen2+1) != $winlen) {
	print "window or half window wrong\n";
	exit;
}

# List the amino acids, for seq entropy (will convert to uc here at least ??)

#@aa_single	= ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");

# SET UP THE AMINO ACID PROPENSITIES FOR THE FOLLOWING PROPERTIES:
#	dis hash is Rune and Linding disorder propensity : GlobPlot : NAR (2003) 31:3701
#	chr hash gives just -1 to (D,E) and +1 to (K,R)
#		use chr to get chr_mean_abs over window length winlen = ABS (sum(chr)/winlen)
#	KD hash is Kyte-Doolittle hydropathy index : JMB (1982) 157:105
#		scale/norm following Uversky : Proteins (2000) 41:415, (-4.5,+4.5) to (0,1)
#		use KDnorm to get KDnorm_mean = sum(KDnorm) / winlen

$dis {'A'}	= -0.26154;
$dis {'C'}	= -0.01515;
$dis {'D'}	=  0.22763;
$dis {'E'}	= -0.20469;
$dis {'F'}	= -0.22557;
$dis {'G'}	=  0.43323;
$dis {'H'}	= -0.00122;
$dis {'I'}	= -0.42224;
$dis {'K'}	= -0.10009;
$dis {'L'}	= -0.33793;
$dis {'M'}	= -0.22590;
$dis {'N'}	=  0.22989;
$dis {'P'}	=  0.55232;
$dis {'Q'}	= -0.18768;
$dis {'R'}	= -0.17659;
$dis {'S'}	=  0.14288;
$dis {'T'}	=  0.00888;
$dis {'V'}	= -0.38618;
$dis {'W'}	= -0.24338;
$dis {'Y'}	= -0.20751;
$dis {'X'}	=  0.0;

$chr {'A'}	=  0.0;
$chr {'C'}	=  0.0;
$chr {'D'}	= -1.0;
$chr {'E'}	= -1.0;
$chr {'F'}	=  0.0;
$chr {'G'}	=  0.0;
$chr {'H'}	=  0.0;
$chr {'I'}	=  0.0;
$chr {'K'}	=  1.0;
$chr {'L'}	=  0.0;
$chr {'M'}	=  0.0;
$chr {'N'}	=  0.0;
$chr {'P'}	=  0.0;
$chr {'Q'}	=  0.0;
$chr {'R'}	=  1.0;
$chr {'S'}	=  0.0;
$chr {'T'}	=  0.0;
$chr {'V'}	=  0.0;
$chr {'W'}	=  0.0;
$chr {'Y'}	=  0.0;
$chr {'X'}	=  0.0;

$dis {'a'}	= -0.26154;
$dis {'c'}	= -0.01515;
$dis {'d'}	=  0.22763;
$dis {'e'}	= -0.20469;
$dis {'f'}	= -0.22557;
$dis {'g'}	=  0.43323;
$dis {'h'}	= -0.00122;
$dis {'i'}	= -0.42224;
$dis {'k'}	= -0.10009;
$dis {'l'}	= -0.33793;
$dis {'m'}	= -0.22590;
$dis {'n'}	=  0.22989;
$dis {'p'}	=  0.55232;
$dis {'q'}	= -0.18768;
$dis {'r'}	= -0.17659;
$dis {'s'}	=  0.14288;
$dis {'t'}	=  0.00888;
$dis {'v'}	= -0.38618;
$dis {'w'}	= -0.24338;
$dis {'y'}	= -0.20751;
$dis {'x'}	=  0.0;

$chr {'a'}	=  0.0;
$chr {'c'}	=  0.0;
$chr {'d'}	= -1.0;
$chr {'e'}	= -1.0;
$chr {'f'}	=  0.0;
$chr {'g'}	=  0.0;
$chr {'h'}	=  0.0;
$chr {'i'}	=  0.0;
$chr {'k'}	=  1.0;
$chr {'l'}	=  0.0;
$chr {'m'}	=  0.0;
$chr {'n'}	=  0.0;
$chr {'p'}	=  0.0;
$chr {'q'}	=  0.0;
$chr {'r'}	=  1.0;
$chr {'s'}	=  0.0;
$chr {'t'}	=  0.0;
$chr {'v'}	=  0.0;
$chr {'w'}	=  0.0;
$chr {'y'}	=  0.0;
$chr {'x'}	=  0.0;

$KD {'A'}	=  0.11;
$KD {'C'}	=  -0.13;
$KD {'D'}	=  3.49;
$KD {'E'}	=  2.68;
$KD {'F'}	=  -0.32;
$KD {'G'}	=  0.74;
$KD {'H'}	= -2.06;
$KD {'I'}	= -0.60;
$KD {'K'}	= 2.71;
$KD {'L'}	= -0.55;
$KD {'M'}	= -0.10;
$KD {'N'}	= 2.05;
$KD {'P'}	= 2.23;
$KD {'Q'}	= 2.36;
$KD {'R'}	= 2.58;
$KD {'S'}	= 0.84;
$KD {'T'}	= 0.52;
$KD {'V'}	= -0.31;
$KD {'W'}	=  0.30;
$KD {'Y'}	=  0.68;
$KD {'X'}	=  0.0;

# Normalise the KD hydropathy values.

#@KD_keys	= keys %KD;
#foreach $keyhere (@KD_keys) {
#  $KDtmp		= $KD{$keyhere} + 4.5;
#  $KDnorm{$keyhere}	= $KDtmp / 9.0;
#  print "KDnorm for $keyhere is $KDnorm{$keyhere}\n";
#}

# get the sequences

open (SEQS_IN, "KD_calc_in.txt") or die "cannot open KD_calc_in.txt\n";
my $nseq		= 0;
while ($line 		= <SEQS_IN>) {
	chomp $line;
	my @words 	= split ("",$line);
	if (exists $words[0]) {
	    if ($words[0] eq ">") {			# new protein
		$nseq++;
		$seqid[$nseq]	= $line;
#		$seqtmp{$seqid[$nseq]}	= "";
		$seqtmpB[$nseq]		= "";		# swap to array, not hash, since multiple IDs could now be the same (e.g. TMs)
	    }
	    else {
		if ($nseq == 0) {
			print "did not see first seq id\n";
			exit;
		}
#		$seqtmp{$seqid[$nseq]}	.= $line;
		$seqtmpB[$nseq]		.= $line;
	    }
	}
}
$nseqs			= $nseq;
close (SEQS_IN);

for ($nseq=1; $nseq<=$nseqs; $nseq++) {		# remove * sequence end flags
#	$seqhere	= $seqtmp{$seqid[$nseq]};
	$seqhere	= $seqtmpB[$nseq];
	$lentmp		= length $seqhere;
	@words		= split ("",$seqhere);
	$naa		= 0;
	for ($n=1; $n<=$lentmp; $n++) {
	  if ($words[$n-1] ne " ") {
	    if ($words[$n-1] eq "*") {
	      if ($naa == 0) {
	        print "STOPPING - first seq char is a * seq end flag\n";
		exit;
	      }
	    } else {
	      $naa++;
	      $sequpd[$naa]	= $words[$n-1];
	    }
	  }
	}
	$seqstr		= "";
	for ($n=1; $n<=$naa; $n++) {
	  $seqstr	.= $sequpd[$n];
	}
#	$seq{$seqid[$nseq]}	= $seqstr;	# move to array to distinuisgh same ID (e.g. with TMs)
	$seqB[$nseq]		= $seqstr;
}

# calculate the disorder (by windowed propensity) values
# also get aa site charge and KD window values, output with overall seq info

if ($TM_goodbad eq "yes") { $nTMseqs = 0; }

open (SEQS_OUT, ">KDcalc_out.txt") or die "cannot open KDcalc_out.txt\n";
for ($nseq=1; $nseq<=$nseqs; $nseq++) {				# OPEN SEQUENCES LOOP
#	$seqhere	= $seq{$seqid[$nseq]};
	$seqhere	= $seqB[$nseq];
	$seqlen[$nseq]	= length $seqhere;
	@words		= split ("",$seqhere);
	$winmax		= -999;
	$nwinmax	=  999;
	$winmaxA[$nseq]		= $winmax;
	$nwinmaxA[$nseq]	= $nwinmax;

	if ($TM_goodbad eq "yes") {
	  $KD_TMsum	= 0;
	  $nKD_TMsum	= 0;
	}

	for ($naa=1; $naa<=$seqlen[$nseq]; $naa++) {		# OPEN AAs in SEQs LOOP
		$wlow			= $naa - $winlen2;
		$wtop			= $naa + $winlen2;
		if (exists $chr{$words[$naa-1]}) {
		  if    ($chr{$words[$naa-1]} > 0) {
		    $chrthis		= "+";
		  } elsif ($chr{$words[$naa-1]} < 0) {
		    $chrthis		= "-";
		  } else {
		    $chrthis		= " ";
		  }
		} else {
		  $chrthis		= " ";
		}
		$propsum		= 0;
		$nprop			= 0;
		$KDwin			= 0;
		$nKDwin			= 0;

		if (exists $KD{$words[$naa-1]}) {	# get the per site avg, std dev, as well as for windows
		  $nx[$naa]++;
		  $xsum[$naa]		+= $KD{$words[$naa-1]};
		  $xsqsum[$naa]		+= $KD{$words[$naa-1]} * $KD{$words[$naa-1]};
		}

		if ($TM_goodbad eq "yes") {		# set from TM_GOODBAD env variable, assume we have TMs 1-20 (say)
		  if ($naa <= 20) {			# and then rank (to file ??) biggest KD differences within a protein
		    if (exists $KD{$words[$naa-1]}) {
		      $KD_TMsum		+= $KD{$words[$naa-1]};
		      $nKD_TMsum++;
		    }
		    if ($naa == 20) {			# record vaerage KD for this starting 20 aas (assumed TM in this option)
		      $nTMseqs++;
		      $KD_TM[$nTMseqs]		= $KD_TMsum / $nKD_TMsum;
		      $KD_TM_ID[$nTMseqs]	= $seqid[$nseq];

print "KD TM value stored = $KD_TM[$nTMseqs] for seq ID = $KD_TM_ID[$nTMseqs]\n";	# bollocks

		    }
		  }
		}

		for ($nw=$wlow; $nw<=$wtop; $nw++){		# second the disorder
			if (($nw>=1) and ($nw<=$seqlen[$nseq])) {
			    if (exists $dis{$words[$nw-1]}) {
				$propsum += $dis{$words[$nw-1]};
				$nprop++;
			    }
	   	            if (exists $KD{$words[$nw-1]}) {
			      $nKDwin++;
		  	      $KDwin	+= $KD{$words[$nw-1]};
	        	    } else {
		  	      print " KD PROBLEM naa words = $naa $words[$nw-1]\n";
	        	    }
			}
		}
		if ($nprop != 0) {
			$windis[$naa]		= $propsum/$nprop;
			if ($windis[$naa] > 0) {
				$winprop[$naa]	= "*";
			} else {
				$winprop[$naa]	= " ";
			}
		} else {
			$windis[$naa]		= 0;
			$winprop[$naa]		= " ";
		}
		$windis[$nseq][$naa]	= $windis[$naa];
		$winpropA[$nseq][$naa]	= $winprop[$naa];
		if ($nKDwin != 0) {
			$winKD[$naa]		= $KDwin/$nKDwin;
			$winKDstore[$nseq][$naa] = $winKD[$naa];
			if ($nKDwin == $winlen) {
			  if ($winKD[$naa] > $winmax) {
			    $winmax		= $winKD[$naa];
			    $nwinmax		= $naa;
			    $winmaxA[$nseq]	= $winmax;
			    $nwinmaxA[$nseq]	= $nwinmax;
			  }
			}
		} else {
			$winKD[$naa]		= 0;
			$winKDstore[$nseq][$naa] = $winKD[$naa];
		}
		$winchr[$naa]			= $chrthis;
		$winchrA[$nseq][$naa]		= $winchr[$naa];
	}							# CLOSE AAs in SEQs LOOP

	printf SEQS_OUT "$seqid[$nseq]\n";
	$dist		= $seqlen[$nseq] - $nwinmax;
	printf SEQS_OUT "Highest complete window KD = ";
	printf SEQS_OUT "%5.2f", $winmax;
	printf SEQS_OUT " at aa number $nwinmax, at $dist from CT\n";

	printf SEQS_OUT "SEQ=";
	for ($naa=1; $naa<=$seqlen[$nseq]; $naa++) {
	  printf SEQS_OUT "  $words[$naa-1]  ";
	}

	printf SEQS_OUT "\n";
	printf SEQS_OUT "DIS=";
	for ($naa=1; $naa<=$seqlen[$nseq]; $naa++) {
	  for ($m=1; $m<=5; $m++) {
	    printf SEQS_OUT "$winprop[$naa]";
	  }
	}

	printf SEQS_OUT "\n";
	printf SEQS_OUT " KD=";
	for ($naa=1; $naa<=$seqlen[$nseq]; $naa++) {
	  printf SEQS_OUT "%5.1f", $winKD[$naa];
	  if ($naa < $ntop) {
	    $nx_win[$naa]++;				# increment stats for this location in this segment
	    $xsum_win[$naa]	+= $winKD[$naa];
	    $xsqsum_win[$naa]	+= $winKD[$naa]*$winKD[$naa];
	  }
	}

	printf SEQS_OUT "\n";
	printf SEQS_OUT "CHR=";
	for ($naa=1; $naa<=$seqlen[$nseq]; $naa++) {
	  for ($m=1; $m<=5; $m++) {
	    printf SEQS_OUT "$winchr[$naa]";
	  }
	}
	printf SEQS_OUT "\n\n";
}
# close (SEQS_OUT);

if ($extra_TA ne "yes") { goto CIRCUMVENT; }

for ($nseq=1; $nseq<=$nseqs; $nseq++) {				# OPEN SEQUENCES LOOP
  $nTM[$nseq]		= 0;
  if ($seqlen[$nseq] > $winlen) {
#    $seqhere	= $seq{$seqid[$nseq]};
    $seqhere	= $seqB[$nseq];
    @words		= split ("",$seqhere);
    $wlow		= 1 + $winlen2;
    $wtop		= $seqlen[$nseq] - $winlen2;
    $full		= 'no';
    for ($naa=$wlow; $naa<=$wtop; $naa++) {			# OPEN AAs in SEQs LOOP
	$covered[$naa]	= "no";
    }
    while ($full eq "no") {					# record whether we have anywhere else to fit a TM
	$winMAX		= -5;
	$nwinMAX	= -999;
	for ($naa=$wlow; $naa<=$wtop; $naa++) {			# OPEN AAs in SEQs LOOP
	    if ( ($winKDstore[$nseq][$naa] > $winMAX) and ($covered[$naa] eq "no") ) {
	        $winMAX		= $winKDstore[$nseq][$naa];
		$nwinMAX	= $naa;
	    }
	}
	if ($winMAX >= 1.6) {
	    if ($covered[$nwinMAX] ne 'no') {
		print "STOPPING - covered wrong for nwinMAX\n";
		exit;
	    } else {
		$nTM[$nseq]++;
		$TMcurr			= $nTM[$nseq];
		if ($nwinMAX == -999) {
		  print "STOPPING - nwinMAX still at -999\n";
		  exit;
		}
		$TMcen[$nseq][$TMcurr]	= $nwinMAX;

#??sat
#print "nseq nTM cenTM = $nseq $TMcurr $TMcen[$nseq][$TMcurr]\n";

		$mlow		= $nwinMAX - $winlen - 1;
		if ($mlow < 1) { $mlow = 1; }
		$mtop		= $nwinMAX + $winlen + 1;
		if ($mtop > $seqlen[$nseq]) { $mtop	= $seqlen[$nseq]; }
		for ($maa=$mlow; $maa<=$mtop; $maa++) {		# block out to where next TM centre would be allowed
		    $covered[$maa]	= "yes";		# i.e. winlen (+ 1 for good measure) - see how well works
		}
	    }
	} else {
	    $full		= "yes";
	}
    }
  }
}

open (CT_OUT, ">KDcalc_CT.txt") or die "cannot open KDcalc_CT.txt\n";
open (CT_OUTCUT, ">KDcalc_CT40_1TM.txt") or die "cannot open KDcalc_CT40_1TM.txt\n";
open (CT_FASTA, ">KDcalc_CT40_1TM.fasta") or die "cannot open KDcalc_CT40_1TM.fasta\n";
for ($nseq=1; $nseq<=$nseqs; $nseq++) {				# OPEN SEQUENCES LOOP
  $written[$nseq]		= "no";
}
$dist_cut		= 40;
for ($nseq=1; $nseq<=$nseqs; $nseq++) {				# OPEN SEQUENCES LOOP
  $currmax		= -999;
  $ncurrmax		=  999;
  for ($mseq=1; $mseq<=$nseqs; $mseq++) {			# OPEN SEQUENCES LOOP
    if ($written[$mseq] eq "no") {
      if ($winmaxA[$mseq] > $currmax) {
        $currmax	= $winmaxA[$mseq];
	$ncurrmax	= $mseq;
      }
    }
  }
  $nextseq		= $ncurrmax;
  $written[$nextseq]	= "yes";

#	$seqhere	= $seq{$seqid[$nextseq]};
	$seqhere	= $seqB[$nextseq];

# ?? we have some off error cropping up occasionally, have not sorted it yet ??
if (!defined $seqhere) { print "BOLLOCKS seqhere not defined for nseq = $nseq $seqid[$nextseq]\n"; }

	@words		= split ("",$seqhere);
	@names		= split (" ",$seqid[$nextseq]);
	$dist		= $seqlen[$nextseq] - $nwinmaxA[$nextseq];
        if ( ($dist <=$dist_cut) and ($nTM[$nextseq] == 1) ) {
	  $CT_oneTM	= "yes";
	  $tailanchor[$nextseq]	= "yes";
	} else {
	  $CT_oneTM	= "no";
	  $tailanchor[$nextseq]	= "no";
	}
	printf CT_OUT "$names[0] KDmax= ";
	printf CT_OUT "%5.2f", $winmaxA[$nextseq];
	printf CT_OUT " $dist from CT";
	printf CT_OUT "  and $nTM[$nextseq] predicted TMs\n";
	if ($CT_oneTM eq "yes") {
	  printf CT_OUTCUT "$names[0] KDmax= ";
	  printf CT_OUTCUT "%5.2f", $winmaxA[$nextseq];
	  printf CT_OUTCUT " $dist from CT";
	  printf CT_OUTCUT "  and $nTM[$nextseq] predicted TMs\n";
	  printf CT_FASTA "$names[0]\n";
	  for ($nnaa=1; $nnaa<=$seqlen[$nextseq]; $nnaa++) {
	    printf CT_FASTA "$words[$nnaa-1]";
	  }
	  printf CT_FASTA "\n";
	}

	$low		= $seqlen[$nextseq] - 59;
	printf CT_OUT "CTSEQ=";
	if ( $CT_oneTM eq "yes" ) { printf CT_OUTCUT "CTSEQ="; }
	for ($naa=$low; $naa<=$seqlen[$nextseq]; $naa++) {
	  if ($naa >= 1) {
	    printf CT_OUT "$words[$naa-1]";
	    if ($CT_oneTM eq "yes") { printf CT_OUTCUT "$words[$naa-1]" };
	  } else {
	    printf CT_OUT " ";
	    if ($CT_oneTM eq "yes") { printf CT_OUTCUT " "; }
	  }
	}
	printf CT_OUT "\n";
	printf CT_OUT "CTDIS=";
	if ($CT_oneTM eq "yes") {
	  printf CT_OUTCUT "\n";
	  printf CT_OUTCUT "CTDIS=";
	}
	for ($naa=$low; $naa<=$seqlen[$nextseq]; $naa++) {
	  if ($naa >= 1) {
	    printf CT_OUT "$winpropA[$nextseq][$naa]";
	    if ($CT_oneTM eq "yes") { printf CT_OUTCUT "$winpropA[$nextseq][$naa]"; }
	  } else {
	    printf CT_OUT " ";
	    if ($CT_oneTM eq "yes") { printf CT_OUTCUT " "; }
	  }
	}
	printf CT_OUT "\n";
	printf CT_OUT "MAXKD=";
	if ($CT_oneTM eq "yes") { printf CT_OUTCUT "\nMAXKD="; }
	for ($naa=$low; $naa<=$seqlen[$nextseq]; $naa++) {
	  if ($naa >= 1) {
	    if ($naa == $nwinmaxA[$nextseq]) {
	      printf CT_OUT "M";
	      if ($CT_oneTM eq "yes") { printf CT_OUTCUT "M"; }
	    } else {
	      printf CT_OUT " ";
	      if ($CT_oneTM eq "yes") { printf CT_OUTCUT " "; }
	    }
	  } else {
	    printf CT_OUT " ";
	    if ($CT_oneTM eq "yes") { printf CT_OUTCUT " "; }
	  }
	}
	printf CT_OUT "\n";
	printf CT_OUT "CTCHR=";
	if ($CT_oneTM eq "yes") { printf CT_OUTCUT "\nCTCHR="; }
	for ($naa=$low; $naa<=$seqlen[$nextseq]; $naa++) {
	  if ($naa >= 1) {
	    printf CT_OUT "$winchrA[$nextseq][$naa]";
	    if ($CT_oneTM eq "yes") { printf CT_OUTCUT "$winchrA[$nextseq][$naa]"; }
	  } else {
	    printf CT_OUT " ";
	    if ($CT_oneTM eq "yes") { printf CT_OUTCUT " "; }
	  }
	}
	printf CT_OUT "\n\n";
	if ($CT_oneTM eq "yes") { printf CT_OUTCUT "\n\n"; }
}
close (CT_OUT);
close (CT_OUTCUT);
close (CT_FASTA);

#??thurs below
# one more analaysis - get profiles of various props over a TM segment:
# aa KD [1] : win KD [2] : charge [3]
# averaged over:
# ALL predicted TMs : TA protein TMs : a specified subset (if this file exists)
# ??thurs - SHOULD ALSO PUT IN AA COMPOSITIONS

$nprof_all		= 0;
$nprof_ta		= 0;
$nprof_sel		= 0;
if (-e "./TMselect.txt") {
    $select		= "yes";
    open (SEL, "<TMselect.txt") or die "cannot open TMselect.txt\n";
    for ($nseq=1; $nseq<=$nseqs; $nseq++) {	# OPEN SEQUENCES LOOP
	$selected[$nseq]	= "no";
    }
    while ($line 	= <SEL>) {
	chomp $line;
	@words		= split (" ",$line);
	if (exists $words[0]) {				# check against the sequence list
	    for ($nseq=1; $nseq<=$nseqs; $nseq++) {	# OPEN SEQUENCES LOOP
	        @names	= split (" ",$seqid[$nseq]);
		if ($names[0]	eq $words[0]) {
		    $selected[$nseq]	= "yes";
		}
	    }
	}
    }
    close (SEL);
} else {
    $select		= "no";
}
for ($maa=1; $maa<=$winlen; $maa++) {
  for ($what=1; $what<=3; $what++) {
    $prof_all[$what][$maa]	= 0;
    $prof_ta[$what][$maa]	= 0;
    $prof_sel[$what][$maa]	= 0;
  }
}
for ($nf=1; $nf<=10; $nf++) {					# zero the window flank counters
  for ($what=1; $what<=3; $what++) {
   for ($end=1; $end<=2; $end++) {
    $flank_all[$end][$what][$nf]	= 0;
    $nflank_all[$end][$what][$nf]	= 0;			# have separate counters - beyond window - maybe null
    $flank_ta[$end][$what][$nf]		= 0;
    $nflank_ta[$end][$what][$nf]	= 0;
    $flank_sel[$end][$what][$nf]	= 0;
    $nflank_sel[$end][$what][$nf]	= 0;
   }
  }
}
for ($nseq=1; $nseq<=$nseqs; $nseq++) {				# OPEN SEQUENCES LOOP
  if ($nTM[$nseq] > 0) {
#    $seqhere		= $seq{$seqid[$nseq]};
    $seqhere		= $seqB[$nseq];
    @words		= split ("",$seqhere);
    for ($iTM=1; $iTM<=$nTM[$nseq]; $iTM++) {
        $cenhere	= $TMcen[$nseq][$iTM];
        $wlow		= $cenhere - $winlen2;
        $wtop		= $cenhere + $winlen2;
	$nprof_all++;
	if (exists $tailanchor[$nseq]) {
	    if ($tailanchor[$nseq] eq "yes") {$nprof_ta++; }
	} else {
	    print "STOPPING - tailanchor unset for nseq = $nseq\n";
	    exit;
	}
	if (($select eq "yes") and ($selected[$nseq] eq "yes")) {$nprof_sel++; }
        for ($maa=$wlow; $maa<=$wtop; $maa++) {			# OPEN AAs in TM LOOP
	    $paa		= $maa - $wlow + 1;
	    $prof_all[1][$paa]	+= $KD{$words[$maa-1]};
	    $prof_all[2][$paa]	+= $winKDstore[$nseq][$maa];
	    $prof_all[3][$paa]	+= $chr{$words[$maa-1]};
	    if ($tailanchor[$nseq] eq "yes") {
	      $prof_ta[1][$paa]	+= $KD{$words[$maa-1]};
	      $prof_ta[2][$paa]	+= $winKDstore[$nseq][$maa];
	      $prof_ta[3][$paa]	+= $chr{$words[$maa-1]};
	    }
	    if (($select eq "yes") and ($selected[$nseq] eq "yes")) {
	     $prof_sel[1][$paa]	+= $KD{$words[$maa-1]};
	     $prof_sel[2][$paa]	+= $winKDstore[$nseq][$maa];
	     $prof_sel[3][$paa]	+= $chr{$words[$maa-1]};
	    }
	}

	for ($end=1; $end<=2; $end++) {
	  if ($end == 1) {
	    $first		= $cenhere - $winlen2 - 10;
	  } else {
	    $first		= $cenhere + $winlen2 + 1;
	  }
	  $last			= $first + 9;
          for ($maa=$first; $maa<=$last; $maa++) {			# OPEN AAs in TM LOOP
	   $paa		= $maa - $first + 1;
	   if (($maa >= 1) and ($maa <= $seqlen[$nseq])) {
	    $flank_all[$end][1][$paa]	+= $KD{$words[$maa-1]};
	    $flank_all[$end][3][$paa]	+= $chr{$words[$maa-1]};
	    $nflank_all[$end][1][$paa]++;
	    $nflank_all[$end][3][$paa]++;
	    if ($tailanchor[$nseq] eq "yes") {
	      $flank_ta[$end][1][$paa]	+= $KD{$words[$maa-1]};
	      $flank_ta[$end][3][$paa]	+= $chr{$words[$maa-1]};
	      $nflank_ta[$end][1][$paa]++;
	      $nflank_ta[$end][3][$paa]++;
	    }
	    if (($select eq "yes") and ($selected[$nseq] eq "yes")) {
	     $flank_sel[$end][1][$paa]	+= $KD{$words[$maa-1]};
	     $flank_sel[$end][3][$paa]	+= $chr{$words[$maa-1]};
	     $nflank_sel[$end][1][$paa]++;
	     $nflank_sel[$end][3][$paa]++;
	    }
	  }
	}
      }
    }
  }
}

open (PROF, ">KDcalc_profiles.txt") or die "cannot open KDcalc_profiles.txt\n";
printf PROF "total number of TMs found for (all, TA, sel) = ( $nprof_all $nprof_ta $nprof_sel )\n";
printf PROF "window    all   all   all    TA    TA    TA   sel   sel   sel\n";
printf PROF "   loc     KD winKD     Q    KD winKD     Q    KD winKD     Q\n";

$end		= 1;
for ($maa=1; $maa<=10; $maa++) {
  printf PROF "    %2.0f ", $maa-10;
  for ($what=1; $what<=3; $what++) {
    if ($nflank_all[$end][$what][$maa] != 0) { $flank_all[$end][$what][$maa] /= $nflank_all[$end][$what][$maa]; }
    printf PROF "%6.2f", $flank_all[$end][$what][$maa];
  }
  for ($what=1; $what<=3; $what++) {
    if ($nflank_ta[$end][$what][$maa] != 0) { $flank_ta[$end][$what][$maa] /= $nflank_ta[$end][$what][$maa]; }
    printf PROF "%6.2f", $flank_ta[$end][$what][$maa];
  }
  for ($what=1; $what<=3; $what++) {
    if (($select eq "yes") and ($nflank_sel[$end][$what][$maa] != 0)) { $flank_sel[$end][$what][$maa] /= $nflank_sel[$end][$what][$maa]; }
    printf PROF "%6.2f", $flank_sel[$end][$what][$maa];
  }
  printf PROF "\n";
}

for ($maa=1; $maa<=$winlen; $maa++) {
  printf PROF "    %2.0f ", $maa;
  for ($what=1; $what<=3; $what++) {
    if ($nprof_all != 0) { $prof_all[$what][$maa]	= $prof_all[$what][$maa] / $nprof_all; }
    printf PROF "%6.2f", $prof_all[$what][$maa];
  }
  for ($what=1; $what<=3; $what++) {
    if ($nprof_ta != 0) { $prof_ta[$what][$maa]	= $prof_ta[$what][$maa] / $nprof_ta; }
    printf PROF "%6.2f", $prof_ta[$what][$maa];
  }
  for ($what=1; $what<=3; $what++) {
    if (($select eq "yes") and ($nprof_sel != 0)) { $prof_sel[$what][$maa]	= $prof_sel[$what][$maa] / $nprof_sel; }
    printf PROF "%6.2f", $prof_sel[$what][$maa];
  }
  printf PROF "\n";
}

$end		= 2;
for ($maa=1; $maa<=10; $maa++) {
  printf PROF "    %2.0f ", $maa+$winlen;
  for ($what=1; $what<=3; $what++) {
    if ($nflank_all[$end][$what][$maa] != 0) { $flank_all[$end][$what][$maa] /= $nflank_all[$end][$what][$maa]; }
    printf PROF "%6.2f", $flank_all[$end][$what][$maa];
  }
  for ($what=1; $what<=3; $what++) {
    if ($nflank_ta[$end][$what][$maa] != 0) { $flank_ta[$end][$what][$maa] /= $nflank_ta[$end][$what][$maa]; }
    printf PROF "%6.2f", $flank_ta[$end][$what][$maa];
  }
  for ($what=1; $what<=3; $what++) {
    if (($select eq "yes") and ($nflank_sel[$end][$what][$maa] != 0)) { $flank_sel[$end][$what][$maa] /= $nflank_sel[$end][$what][$maa]; }
    printf PROF "%6.2f", $flank_sel[$end][$what][$maa];
  }
  printf PROF "\n";
}

close (PROF);
#??thurs above - put in aa compositions

CIRCUMVENT:

# output the stats for KD variation at each location

printf SEQS_OUT "PER WINDOW section of KD stats over these WINDOWS, aligned by starting aa\n\n";
printf SEQS_OUT "WINavg=";
for ($i=1; $i<=$ntop; $i++) {
  if ($nx_win[$i] != 0) {
    $x_avg_win[$i]		= $xsum_win[$i] / $nx_win[$i];
    printf SEQS_OUT "%6.2f", $x_avg_win[$i];
  }
}
printf SEQS_OUT "\n\nWINdev=";
for ($i=1; $i<=$ntop; $i++) {
  if ($nx_win[$i] > 1) {
    $tmp_win		= $xsqsum_win[$i] + $nx_win[$i]*$x_avg_win[$i]*$x_avg_win[$i] - 2*$x_avg_win[$i]*$xsum_win[$i];
    if ($tmp_win > 0) {
      $std_dev_win[$i]	= sqrt (($tmp_win / $nx_win[$i]));
      printf SEQS_OUT "%6.2f", $std_dev_win[$i];
    }
  }
}
printf SEQS_OUT "\n\n";

printf SEQS_OUT "PER AASITE section of KD stats over these AASITES, aligned by starting aa\n\n";
printf SEQS_OUT "AA_avg=";
for ($i=1; $i<=$ntop; $i++) {
  if ($nx[$i] != 0) {
    $x_avg[$i]		= $xsum[$i] / $nx[$i];
    printf SEQS_OUT "%6.2f", $x_avg[$i];
  }
}
printf SEQS_OUT "\n\nAA_dev=";
for ($i=1; $i<=$ntop; $i++) {
  if ($nx[$i] > 1) {
    $tmp		= $xsqsum[$i] + $nx[$i]*$x_avg[$i]*$x_avg[$i] - 2*$x_avg[$i]*$xsum[$i];
    if ($tmp > 0) {
      $std_dev[$i]	= sqrt (($tmp / $nx[$i]));
      printf SEQS_OUT "%6.2f", $std_dev[$i];
    }
  }
}
printf SEQS_OUT "\n";

close (SEQS_OUT);

# AT PRESENT, just take contiguous TM pairings for comparison
# ?? COULD CONSIDER ADDING taking all TM combinations within a single protein ??

if ($TM_goodbad eq "yes") {		# now for ranking (and output) by biggest within seq TM KD difference - functional ??
  open (RANKING, ">TM_goodbad_ranking.txt") or die "cannot open TM ranking file\n";
  $nIDs			= 0;
  $ID_last		= $KD_TM_ID[1];
  $KD_last		= $KD_TM[1];
  $KD_diff_max		= 0;
  for ($iTM=2; $iTM<=$nTMseqs; $iTM++) {	# assume that TMs from same ID/protein are contiguous in our arrays
    $ID_here		= $KD_TM_ID[$iTM];
    $KD_here		= $KD_TM[$iTM];
    if ($ID_here eq $ID_last) {
      $KD_diff		= abs ($KD_here - $KD_last);	# and in fact we are ONLY TESTING consecutive TMs - ?? COULD BE UPDATED
      if ($KD_diff > $KD_diff_max) { $KD_diff_max = $KD_diff; }
      $KD_last		= $KD_here;
    } else {					# process the last protein
      $nIDs++;
      $ID_store[$nIDs]		= $ID_last;
      $diff_store[$nIDs]	= $KD_diff_max;
      printf RANKING "$ID_last ";
      printf RANKING "%6.2f\n", $KD_diff_max;		# ?? FOR NOW WRITE - SHOULD SORT AND THEN WRITE ??
      $ID_last			= $ID_here;
      $KD_last			= $KD_here;
      $KD_diff_max		= 0;
    }
  }
  $nIDs++;			# that last one to do
  $ID_store[$nIDs]	= $ID_here;
  $diff_store[$nIDs]	= $KD_diff_max;
  printf RANKING "$ID_here ";
  printf RANKING "%6.2f\n", $KD_diff_max;		# ?? FOR NOW WRITE - SHOULD SORT AND THEN WRITE ??
}

# ?? SORT - and consider beyond differencing just contiguous TMs ??

close (RANKING);

exit;
