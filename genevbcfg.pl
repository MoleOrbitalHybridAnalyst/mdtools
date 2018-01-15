#!/bin/perl

use File::Basename;
use warnings;

$debug=1;
sub say_something_in_comments
{
    my $something=shift;
    my $fh=shift;
    print $fh "\n";
    print $fh ":" x 27;
    print $fh "\n";
    print $fh "::$something\n";
    print $fh ":" x 27;
    print $fh "\n\n";
}

$where_is_this_script=dirname($0);
$input_count=0;
print "write system specific setups based on evb.cfg.head\n";
my @all_the_lines_of_evbcfg_head=qw{};
print "which evb.cfg.head?\n" if $debug;
while(<STDIN>)
{
    if($input_count==0) {#get to know the filename of evb.cfg.head
        chomp;
        if(-e $_)   {$fn_evbcfg_head=$_;}
        else        {$fn_evbcfg_head=$where_is_this_script."/".$_;}
        print "$fn_evbcfg_head\n" if $debug;
        open IN,"<",$fn_evbcfg_head;
        while(<IN>) {
            chomp;
            push @all_the_lines_of_evbcfg_head,$_;
        }
        close IN;
        print "what name of evb.cfg output?\n" if $debug;
    }
    elsif($input_count==1) {#get to know the file name of output
        print if $debug;
        chomp;
        open OUT,">",$_;
        for(@all_the_lines_of_evbcfg_head) {print OUT $_."\n";}
        print "what about model selection?\n" if $debug;
    }
    elsif($input_count==2) {#get to know what models 
        print if $debug;
        chomp;
        &say_something_in_comments("Model Selection",\*OUT);
        my @models=split ' ';
        for(@models) {
            print OUT "#define $_\n";
        }
        print "what about molecules to declare?\n" if $debug;
    }
    elsif($input_count==3 or $input_count==4) {#get to know the molecules
        print if $debug;
        chomp;
        &say_something_in_comments("Declare Molecules",\*OUT) if $input_count==3;
        print "HERE\n" if $debug==2;
        my $molecule_count=0;
        my @molecules=split ' ';
        for(@molecules) {
            print OUT sprintf "settype $_\t%d\n",++$molecule_count if $input_count==4;
        }
        if($input_count==3) {
            for(split ' ') { print OUT "#define $_\n";}
        }
        print "what about reactions?\n" if ($debug and $input_count==4);
    }
    elsif($input_count==5) {
        print if $debug;
        chomp;
        &say_something_in_comments("Reaction Types",\*OUT);
        my @reactions=split ' ';
        for(@reactions) {
            print OUT "#define $_\n";
        }
        print "what about LAMMPS types?\n" if $debug;
    }
    elsif($input_count==6) {
        print if $debug;
        chomp;
        &say_something_in_comments("Declare LAMMPS Types",\*OUT);
        print OUT "#include \"$_\"\n";
        print "what about evb parameters?\n";
    }
    elsif($input_count==7) {
        print if $debug;
        chomp;
        &say_something_in_comments("Declare EVB parameters",\*OUT);
        print OUT "#include \"$_\"\n";
    }
    $input_count++;
    last if $input_count > 7;
}
print OUT "\n";
print OUT ":" x 78;
print OUT "\n";
print OUT ":" x 78;
close OUT;
