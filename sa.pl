#!/usr/bin/env perl
###############################################################################
#
#    sa.pl
#
#    Reads in scaffold data and summarizes into table
#
#    Copyright (C) Inka VanWonterghem
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;
use Carp;
#use Math;

#CPAN modules
use Bio::SeqIO;
use Data::Dumper;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# edit here to log all external commands
my $global_log_commands = 0;

# ext command failure levels
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

# get input params and print copyright
printAtStart();
my $global_options = checkParams();

######################################################################
# CODE HERE
######################################################################

#print(Dumper($global_options));
# default options
createDefault('read_cutoff',1000000);

#print($global_options->{'read_cutoff'}."\n");

# read file
my $ss_fileHandle = openRead($global_options->{'sspace_summary'});
foreach my $line (<$ss_fileHandle>) {
    #print($line);
}
close($ss_fileHandle);

# write to a file
my $out_fh = openWrite("bob");
print $out_fh "freddy\n", plus2(8), "\n";
close($out_fh);

# bioperl
my $seqio = Bio::SeqIO->new(-file => $global_options->{'sspace_summary'}, '-format' => 'Fasta');
while(my $seq = $seqio->next_seq)
{
    my $seq_string = $seq->seq;
    my $seq_id = $seq->id;
    my $seq_len = length($seq_string);
    #print $seq_id, $seq_len, "\n";
}

# hash of hashes
my %hoh = ();
my @data = ([5,3],[3,3],[3,2],[3,3],[5,6],[5,3],[3,3]);

print "\n\n\n\n\n\n\nBEFORE:\n";
print Dumper @data;
print Dumper %hoh;

my %read_count_store = ();
my $dumb_read_amount = 8;

foreach my $dat (@data) {
  my $from_bin = @{$dat}[0];
  my $to_bin = @{$dat}[1];

  if($from_bin eq $to_bin) {
    # we should store the read count
    if(not exists $read_count_store{$to_bin}) {
      my @tmp_array = ($dumb_read_amount);
       $read_count_store{$to_bin} = \@tmp_array;
    } else {
      push(@{$read_count_store{$to_bin}}, $dumb_read_amount);
    }
  }

  foreach my $nuffin ((1..2)) {
    if(not exists $hoh{$from_bin}) {
      # never seen before
      my %tmp_hash = ();
        $tmp_hash{$to_bin} = 1;
        $hoh{$from_bin} = \%tmp_hash;
    } else {
      # we have an existing hash
      if (exists $hoh{$from_bin}{$to_bin}){
        $hoh{$from_bin}{$to_bin} += 1;
      } else {
        $hoh{$from_bin}{$to_bin} = 1;
      }
    }
    $from_bin = @{$dat}[1];
    $to_bin = @{$dat}[0];
  }
}

print "AFTER:\n";

#print mean($read_count_store[3]);
print Dumper %hoh;
print Dumper %read_count_store;
######################################################################
# CUSTOM SUBS
######################################################################

sub plus2 {
  my ($in) = @_;
  return $in + 2;

}

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    # "long|single<type>"
    #
    # :i number (integer)
    # :f number (decimal)
    # :s string
    # +  flag
    #
    my @standard_options = ( "help|h+", "sspace_summary|s:s", "read_cutoff|r:i");
    my %options;

    # Add any other command line options, and the code to handle them
    #
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    #if(!exists $options{''} ) { printParamError (""); }
    if(!exists $options{'sspace_summary'} ) { printParamError ("You need to supply the sspace summary CSV file"); }

    return \%options;
}

sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #
    my ($error) = @_;
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}

sub createDefault
{
    #-----
    # Set default values for parameters
    #
    my ($option_name, $default_value) = @_;
    if(not exists $global_options->{$option_name})
    {
      $global_options->{$option_name} = $default_value;
    }
}

######################################################################
# FILE IO

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

######################################################################
# EXTERNAL COMMANDS
#
# checkAndRunCommand("ls", {
#                          -a => ""
#                          },
#                          WARN_ON_FAILURE);

sub checkFileExists {
    #-----
    # Does a file exists?
    #
    my ($file) = @_;
    unless(-e $file) {
        croak "**ERROR: $0 : Cannot find:\n$file\n";
    }
}

sub logExternalCommand
{
    #-----
    # Log a command line command to the command line!
    #
    if(1 == $global_log_commands) {
        print $_[0], "\n";
    }
}

sub isCommandInPath
{
    #-----
    # Is this command in the path?
    #
    my ($cmd, $failure_type) = @_;
    if (system("which $cmd |> /dev/null")) {
        handleCommandFailure($cmd, $failure_type);
    }
}

sub runExternalCommand
{
    #-----
    # Run a command line command on the command line!
    #
    my ($cmd) = @_;
    logExternalCommand($cmd);
    system($cmd);
}

sub checkAndRunCommand
{
    #-----
    # Run external commands more sanelier
    #
    my ($cmd, $params, $failure_type) = @_;

    isCommandInPath($cmd, $failure_type);

    # join the parameters to the command
    my $param_str = join " ", map {formatParams($_)} @{$params};

    my $cmd_str = $cmd . " " . $param_str;

    logExternalCommand($cmd_str);

    # make sure that all went well
    if (system($cmd_str)) {
         handleCommandFailure($cmd_str, $failure_type)
    }
}

sub formatParams {

    #---------
    # Handles and formats the different ways of passing parameters to
    # checkAndRunCommand
    #
    my $ref = shift;

    if (ref($ref) eq "ARRAY") {
        return join(" ", @{$ref});
    } elsif (ref($ref) eq "HASH") {
        return join(" ", map { $_ . " " . $ref->{$_}} keys %{$ref});
    }
    croak 'The elements of the $params argument in checkAndRunCommand can ' .
        'only contain references to arrays or hashes\n';
}


sub handleCommandFailure {
    #-----
    # What to do when all goes bad!
    #
    my ($cmd, $failure_type) = @_;
    if (defined($failure_type)) {
        if ($failure_type == DIE_ON_FAILURE) {
            croak "**ERROR: $0 : " . $! . "\n";
        } elsif ($failure_type == WARN_ON_FAILURE) {
            carp "**WARNING: $0 : " . $! . "\n";
        }
    }
}


######################################################################
# MISC

sub printAtStart {
print<<"EOF";
----------------------------------------------------------------
 $0
 Copyright (C) Inka VanWonterghem

 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
----------------------------------------------------------------
EOF
}

__DATA__

=head1 NAME

    sa.pl

=head1 COPYRIGHT

   copyright (C) Inka VanWonterghem

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

   Insert detailed description here

=head1 SYNOPSIS

    sa.pl -s <SSPACE_CSV> [-r <INT>] [-help|h]

      -sspace_summary -s           sspace CSV from Jason's scripts
      [-read_cutoff -r]            reject all links with less than this many links
      [-help -h]                   Displays basic usage information

=cut

