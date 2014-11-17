#!/usr/bin/perl
#===============================================================================
#
#         FILE: bioDbnet.pl
#
#        USAGE: ./bioDbnet.pl inputtype inputtaxid inputfileaddress outputfile
#
#  DESCRIPTION: This script is used to query biodbnet to match genesymbols and ENST IDs to Refseq ID
#      OPTIONS: ./bioDbnet.pl genesymbolandsynonyms 9606 missingGenename.9606_uniq  9606_rnaGeneName2refseqs.hash
# REQUIREMENTS: use file name and input type (2rd parameter) as input
#         BUGS: ---
#        NOTES: the REST service can't accept more than 500 input at one time
#       AUTHOR: Jiali Wang
# ORGANIZATION: University luxembourg
#      VERSION: 1.0
#      CREATED: 11/10/2014 02:17:32 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Storable;
use List::MoreUtils 'uniq';
use LWP::Simple;
use local::lib;
use Data::Dumper qw(Dumper);
use JSON::Parse 'parse_json';

my $project_dir = "/work/projects/pdgfr_kit/miRNA/wj/";
my $data_dir = $project_dir . "data/";
my $inputtype=$ARGV[0];#inputtype could be genesymbolandsynonyms,ensembltranscriptid, ensemblgeneid
my $inputtaxid=$ARGV[1];#10090 or 9606
my $inputfile=$ARGV[2];#uniq files of uniq input ids
my $outputfile=$ARGV[3];#the file that to be read in and open as hash
my %target_hash=%{retrieve($data_dir.'annotations/oldhash/'.$outputfile)};
my $baseurl= 'http://biodbnet.abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&input=';
#my $url = 'http://biodbnet.abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&input=genesymbolandsynonyms&inputValues=KIT,PDGFRA&outputs=refseqmrnaaccession&taxonId=9606&format=row';
my @all=();
my $counter=1;
my $temp="";
print $inputfile;
if(open (INPUT,$data_dir.'missing/'.$inputfile)){
    @all=<INPUT>;
    close INPUT;
     $counter=1;
    
     $temp="";
}else{print "file can't be open!";}
    while(@all){
        my $line =shift(@all);
        #print scalar(@all)."\n";
        chomp ($line);
        if ($counter>=400){  


            $temp=$temp.','.$line;
            #call and excute one rest
            print "start RESTing \n";
            &getrest($temp);

            $counter=1;
            $temp="";
            sleep(20);
        }
        elsif(!@all){
            $temp=$temp.','.$line;
            #call and excute one rest
            print "start RESTing \n";
            &getrest($temp);
            
            } 
        else{
            $counter++;
            $temp=$temp.','.$line;
        }#else


    }#while(@all)


#print Dumper \%target_hash;
store \%target_hash,$data_dir.'annotations/'.$outputfile.'_new';     


sub getrest{

    my ($temp)= @_;
    my $url = "http://biodbnet.abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&input=$inputtype&inputValues=$temp&outputs=refseqmrnaaccession&taxonId=$inputtaxid&format=row";
    my $response = get $url;
    die 'Error getting $url' unless defined $response;
    print $response,"\n";
    my $obj= parse_json($response);
#print ref $obj,"\n";
#print ref @{$obj}[1],"\n";
    foreach my $current (@{$obj}){
        my %output=%{$current};
        if ($output{"RefSeq mRNA Accession"} eq "-"){
            next;        
        }#this is empty output
        else{
            my $currentInput=$output{"InputValue"};
            my @refseqs = split('//',$output{"RefSeq mRNA Accession"});
            print "@refseqs\n";
            $target_hash{$currentInput}=\@refseqs;            
        }
        
    
    }
}



