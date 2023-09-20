#!/usr/bin/perl -w
use warnings;
use strict;
use FindBin;                # where was script installed?
use lib $FindBin::Bin;      # use that dir for libs, too
use Getopt::Long;
use Data::Dumper;

################
# TENTAKKL
# Taxonomic Extraction Normalization Transformation and Accounting for Kraken/bracKen fiLes
# S. Conlan (2023)

my $version="0.3";

#options
my $opt = {config => 'minitax.config', comment=>'#', sep=>"\t", verbose=>1, outfmt=>"list", normalize=>"n", precision=>"4" };
if (-t STDIN && ! @ARGV){usage()};
GetOptions ($opt, 'config=s', 'kraken=s@{1,}', 'bracken=s@{1,}', 'verbose=i', 'outfmt=s', 'normalize=s', 'add_lost', 'add_unclassified', 'out=s',
            'precision=i','logfile=s','help');
if ($opt->{help}){usage()};
my $verbose=$opt->{'verbose'};
my $ofh;
if ($opt->{out}){open($ofh,">".$opt->{out}) or die "can not open ".$opt->{out}."\n"}
else {$ofh=\*STDOUT};
my $lvlauto="NA"; #automatically detect bracken level

#print STDERR Dumper($opt);

#########
# Main

# Load taxonomic targets from file
my ($targ,$drop)=load_targets($opt->{config});

my %counts;      #counts table           - target counts and root, unclassified if requested
my %fcount;      #total counts from file - functionally, this is the deminator for bracken noramlization
my %dropped;     #reads in dropped taxa  - needed to adjust kraken counts later
my %root_totals; #bracken counts at root - calculating lost reads

foreach my $f (@{$opt->{'bracken'}})
  {
    print STDERR "INFO: Working on $f...\n" if ($verbose>0);
    my @lvl_stack;   #keep track of level
    my @name_stack;  #keep track of target name
    my $sample=sample_from_file($f);

    my %atlvlsums;   #store sums of reads at a given level (G, S, etc...) for error checking
    open(INF,"<$f") or die "can't open bracken file $f\n";
    while (my $newlin=<INF>)
     {
       chomp($newlin);
       my ($pct,$atbelow,$atlvl,$taxlvl,$taxid,$sciname)=split(/\t/,$newlin);
       $sciname=~s/^\s+//g; #strip off spaces, not reliable for computing on

       #tracking what level taxonomic counts are at, should be just 1 for bracken
       if ($atlvl>0){$lvlauto=$taxlvl;$atlvlsums{$taxlvl}+=$atlvl};

       if ($sciname eq "root")
         {
           #initialize target stack
           push(@lvl_stack,$taxlvl); # R
           push(@name_stack,$sciname); # root
           if ($atlvl>0){die "FATAL: root has atlvl counts\n"};
           $root_totals{$sample}=$atbelow;
           next;
         }

       while(is_less_or_equal($taxlvl,$lvl_stack[$#lvl_stack])==1)
        {
          #current taxonomic level is above (less) than current level
          #pop targets until we get to a parent of the current level
          my $q=pop(@lvl_stack);
          $q=pop(@name_stack);
        }

       #drop taxa user requested dropped, keep track of those reads
       if (exists($drop->{$sciname}))
         {
           #taxon to drop
           $root_totals{$sample}-=$atlvl; #modify root count
           $dropped{$sample}+=$atlvl;     #save for later
           print STDERR join("\t",$lvl_stack[$#lvl_stack],$name_stack[$#name_stack],"<-X-",$sciname,tax_to_num($taxlvl),$pct,$atbelow,$atlvl,$taxlvl,$taxid)."\n" if ($verbose>1);
           next;
         }

       if (exists($targ->{$sciname}))
         {
           #if this is a user-requested target, append it to the stack
           push(@lvl_stack,$taxlvl);
           push(@name_stack,$sciname);
         }
       if ($atlvl>0){print STDERR join("\t",$lvl_stack[$#lvl_stack],$name_stack[$#name_stack],"<---",$sciname,tax_to_num($taxlvl),$pct,$atbelow,$atlvl,$taxlvl,$taxid)."\n" if ($verbose>1)};

       #now that the current target level is the last item in @name_stack, count reads
       $counts{$sample}{$name_stack[$#name_stack]}+=$atlvl;
       $fcount{$sample}+=$atlvl
     }
    close(INF);

    #Error check, bracken should only have atlvl reads at 1 taxonomic level
    if (scalar(keys(%atlvlsums)) !=1)
      {
        #error checking, right now we expect bracken files that only have counts at one taxonomic level
        print STDERR Dumper(\%atlvlsums);
        die "FATAL: more than one level is populated\n";
      }

    # If requested, this adds the difference between what bracken reports at root
    # and the sum of the reads at the bracken assignment level, fcount (e.g., S)
    # add to root and increase the fcount sum (since we are now counting them)
    my $rootdiff=0; #difference between root and assigned
    if ($opt->{'add_lost'})
      {
        $targ->{'root'}=0;
        $rootdiff=$root_totals{$sample}-$fcount{$sample};
        $counts{$sample}{'root'}+=$rootdiff;
        $fcount{$sample}+=$rootdiff;
      }
    print STDERR "INFO: $f ".$atlvlsums{$lvlauto}." reads assigned to $lvlauto of ".$root_totals{$sample}." under root; ".$fcount{$sample}." in table ($rootdiff assigned to root)\n" if ($verbose>0);
  }

##########
# kraken2
# if raw kraken files are supplied, grab unclassified reads

my $kdat;
if (exists($opt->{'kraken'})){$kdat=get_kraken($opt->{'kraken'})}
elsif (exists($opt->{'add_unclassified'})){die "FATAL: user requested unclassified reads but didn't supply kraken files\n"};

# validate kraken data and add reads if requested. kraken reports can supply
# 1. reads bracken dropped (e.g., Reads not distributed)
# 2. unclassified reads
foreach my $k (sort keys %root_totals)
  {
    if (!exists($kdat->{$k})){die "FATAL: $k is in bracken files but not paired kraken files, check names\n";next}
    if (exists($dropped{$k})){$kdat->{$k}->{'root'}-=$dropped{$k}}; #if we dropped reads, adjust those here

    if ( (abs($root_totals{$k} - $kdat->{$k}->{'root'})/$root_totals{$k}) > 0.005 )
      {
        # because of bracken's Reads not distributed, the kraken classified won't exactly match bracken root
        # but if they differ by a lot that's an issue
        print STDERR "WARN: $k, bracken root and kraken classified differ by >0.5% : ".$root_totals{$k}.",".$kdat->{$k}->{'root'}."\n";
      }
    #report those reads and recover if requested...
    my $rootdiff=0;
    if ($opt->{'add_lost'})
      {
        $rootdiff=$kdat->{$k}->{'root'} - $root_totals{$k};
        $counts{$k}{'root'}+=$rootdiff;
        $fcount{$k}+=$rootdiff;
      }
    print STDERR $k." lost some reads during redistribution: ".($kdat->{$k}->{'root'} - $root_totals{$k})." ($rootdiff assigned to root)\n" if ($verbose>0);
    if ($opt->{'add_unclassified'})
      {
        $targ->{'unclassified'}="-1";
        $counts{$k}{'unclassified'}=$kdat->{$k}->{'unclassified'};
      }
  }

##########
# output

if ($opt->{'outfmt'} eq "list")
  {
    #tidy list format
    print $ofh join($opt->{'sep'},"sample","order","taxon","count","b_fraction","k_fraction")."\n";
    foreach my $s (@{$opt->{'bracken'}})
      {
        my $sample=sample_from_file($s);
        foreach my $t (sort {$targ->{$a}<=>$targ->{$b}} keys %{$targ})
          {
            print $ofh join($opt->{'sep'},$sample,$targ->{$t},$t).$opt->{'sep'};
            if ( exists($counts{$sample}{$t}) )
              {
                #counts and bracken normalized
                print $ofh join($opt->{'sep'},$counts{$sample}{$t},
                                              sprintf("%.".$opt->{'precision'}."f",($counts{$sample}{$t}/$fcount{$sample}) ) );
                #if kraken data supplied, report kraken normalized
                if (exists($kdat->{$sample}))
                  {print $ofh $opt->{'sep'}.sprintf("%.".$opt->{'precision'}."f",($counts{$sample}{$t}/($fcount{$sample} + $kdat->{$sample}->{'unclassified'})) )}
                else {print $ofh $opt->{'sep'}."NA"};
                print $ofh "\n";
              }
            else {print $ofh join($opt->{'sep'},"0","0","0")."\n"};
          }
      }
  }
elsif ($opt->{'outfmt'} eq "wide")
  {
    #spreadsheet format
    foreach my $s (@{$opt->{'bracken'}}){print $ofh $opt->{'sep'}.sample_from_file($s)};
    print "\n";

    foreach my $t (sort {$targ->{$a}<=>$targ->{$b}} keys %{$targ})
      {
        print $ofh $t;
        foreach my $s (@{$opt->{'bracken'}})
         {
            my $sample=sample_from_file($s);
            if ( exists($counts{$sample}{$t}) )
              {
                if ($opt->{'normalize'} eq 'n'){print $ofh $opt->{'sep'}.$counts{$sample}{$t}};
                if ($opt->{'normalize'} eq 'b'){print $ofh $opt->{'sep'}.sprintf("%.".$opt->{'precision'}."f",($counts{$sample}{$t}/$fcount{$sample}))};
                if ($opt->{'normalize'} eq 'k'){print $ofh $opt->{'sep'}.sprintf("%.".$opt->{'precision'}."f",($counts{$sample}{$t}/($fcount{$sample} + $kdat->{$sample}->{'unclassified'})) )};
              }
            else {print $ofh $opt->{'sep'}."0"};
          }
        print "\n";
      }
  }
else {die "FATAL: unrecognized outfmt\n"};

#########
# log - TODO, finish this

if ($opt->{'logfile'})
  {
    open(LOG,">$opt->{'logfile'}") or die "can't open log\n";
    print LOG join($opt->{'sep'},"sample","k_root","k_unclass","b_root","b_assigned","b_dropped","b_lost")."\n";
    foreach my $s (@{$opt->{'bracken'}})
      {
        my $sample=sample_from_file($s);
        print LOG $sample;
        if (exists($kdat->{$sample})){print LOG $opt->{'sep'}.$kdat->{$sample}->{'root'}} else {print LOG $opt->{'sep'}."NA"};
        if (exists($kdat->{$sample})){print LOG $opt->{'sep'}.$kdat->{$sample}->{'unclassified'}} else {print LOG $opt->{'sep'}."NA"};
        print LOG $opt->{'sep'}.$root_totals{$sample};
        print LOG $opt->{'sep'}.$fcount{$sample};
        if (exists($dropped{$sample})){print LOG $opt->{'sep'}.$dropped{$sample}} else {print LOG $opt->{'sep'}."NA"};
        if ($opt->{'add_lost'}){print LOG $opt->{'sep'}.($kdat->{$sample}->{'root'} - $root_totals{$sample})} else {print LOG $opt->{'sep'}."NA"};
        print LOG "\n";        
      }
    close(LOG);
  }

#########
# SUB

sub usage
  {
    print STDERR "TENTAKKL v$version\nTaxonomic Extraction Normalization Transformation and Accounting for Kraken/bracKen fiLes\n\n";
    print STDERR "usage perl -w tentakkl.pl [options] -config minitax.cfg -bracken kraken_output/*_bracken.kreport\n\n";
    print STDERR "-config file       list of target taxa, text file, one taxon per line\n";
    print STDERR "                   preface taxon with '-' to have reads at that level removed. Child taxa unaffected\n";
    print STDERR "-bracken files     bracken-corrected kraken taxonomy reports\n";
    print STDERR "-kraken files      [optional] kraken taxonomy reports, used to extract unclassified reads\n";
    print STDERR "-out file          output filename (STDOUT by default)\n";
    print STDERR "-outfmt list|wide  list=tidy or wide=matrix; list has counts and normalized; wide uses normalize option\n";
    print STDERR "-normalize n|b|k   (outfmt=wide) n=counts\n";
    print STDERR "                                 b=normalize to 1 using bracken total\n";
    print STDERR "                                 k=normalize to 1 using kraken total (incl. unclassified)\n";
    print STDERR "-add_lost          add reads lost during bracken read reassignment to a root taxon\n";
    print STDERR "-add_unclassified  add unclassified reads category\n";
    print STDERR "-precision int     number of decimal places when normalizing\n";
    print STDERR "-logfile file      write a logfile (still under development)\n";
    print STDERR "-verbose int       0=no info, 1=file-level info, 2=taxa-level info\n";
    exit(0);
  }

#sub is_less
#  {
#    #return 1 if first value is less than second
#    #R<D<K<P<C<O<F<G<S
#    #also S<S1<S2
#    my %order=qw( R 1 D 2 K 3 P 4 C 5 O 6 F 7 G 8 S 9 );
#    my ($v1,$v2)=@_;
#    if (tax_to_num($v1)<tax_to_num($v2)){return(1)} else {return(0)};
#  }

sub is_less_or_equal
  {
    #return 1 if first value is less than second
    #R<D<K<P<C<O<F<G<S
    #also S<S1<S2
    my %order=qw( R 1 D 2 K 3 P 4 C 5 O 6 F 7 G 8 S 9 );
    my ($v1,$v2)=@_;
    if (tax_to_num($v1)<=tax_to_num($v2)){return(1)} else {return(0)};
  }

sub tax_to_num
  {
    my %order=qw( R 1 D 2 K 3 P 4 C 5 O 6 F 7 G 8 S 9 );
    my $x=$_[0];
    if (exists($order{$x})){return($order{$x})};
    if ($x=~/^([RDKPCOFGS])(\d)$/){return($order{$1}).".".$2};
    die "can't parse $x to a numeric value\n";    
  }

sub sample_from_file
  {
    my $file=$_[0];
    $file=~s/^.+\///g; #remove path
    $file=~s/\..+$//g; #remove extensions
    if ($file=~/^([A-Z0-9]+)/gi){return($1)}; #special case of grabbing letters/numbers 
    return($file);
  }

sub load_targets
  {
    my %dat;
    my %del;
    my $cfg=$_[0];
    my $n=0;
    open(INF,"<$cfg") or die "can't open target taxa\n";
    while (my $newlin=<INF>)
     {
       chomp($newlin);
       if ($newlin=~/^\s*\#/ || $newlin=~/^\s*$/){next};
       $newlin=~s/^\s+//g;
       $newlin=~s/\s+$//g;
       if ($newlin=~/^\-/)
         {
           #taxa to drop
           $newlin=~s/^\-//;
           $del{$newlin}=1;
         }
       else
         {
           $n++;
           $dat{$newlin}=$n;
         }
     }
    close(INF);
    return(\%dat,\%del);
  }

sub get_kraken
  {
    my %dat;
    my $kf=$_[0];
    foreach my $file (@{$kf})
      {
         open(KRA,"<$file") or die "can't open $file\n";
         while (my $newlin=<KRA>)
           {
             # 18.11  848465  848465  U       0       unclassified
             # 81.89  3835311 1636    R       1       root

             chomp($newlin);
             $newlin=~s/^\s+//g;
             my @f=split(/\s+/,$newlin);
             if ($f[3] eq 'U' && $f[5] eq 'unclassified'){$dat{sample_from_file($file)}{'unclassified'}=$f[1]};
             if ($f[3] eq 'R' && $f[5] eq 'root'){$dat{sample_from_file($file)}{'root'}=$f[1]};
           }
         close(KRA); 
      }
    return(\%dat);
  }
