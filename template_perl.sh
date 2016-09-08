#!/usr/bin/env perl

use English;
$WARNING = 1;
use strict;
use Getopt::Long;
use Data::Dumper;
use YAML qw(LoadFile);
use Cwd qw(abs_path);
use List::Util qw( min max );

# # custom modules in local lib/
# use File::Basename;
# use lib dirname (__FILE__);
# use lib::Levenshtein qw(distance);


my ($fasta,$yaml) = @{read_options()};
my $SETTINGS = LoadFile($yaml);
    print Dumper $SETTINGS; exit;




sub read_options{

    # defaults
    my $fasta = '';
    abs_path($0) =~ /(.*)\//;
    my $yaml = $1 . '/match_fastaIDs_by_keywords.yml';


    GetOptions('fasta=s'  => \$fasta,
               'yaml=s' => \$yaml
               );

    return [$fasta,$yaml];
}
