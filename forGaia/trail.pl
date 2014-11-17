
#===============================================================================
#
#         FILE: bioDbnet.pl
#
#        USAGE: ./bioDbnet.pl inputtype inputtaxid inputfileaddress  
#
#  DESCRIPTION: This script is used to query biodbnet to match genesymbols and ENST IDs to Refseq ID
#      OPTIONS: ---
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
my $inputtype=$ARGV[1];
my $inputtaxid=$ARGV[2];#10090 or 9606
my $inputfile=$ARGV[3];#uniq files of uniq input ids
my $outputfile=$ARGV[4];#the file that to be read in and open as hash
my %target_hash=%{retrieve($data_dir.'annotations/'.$inputtaxid.$outputfile)};
my $baseurl= 'http://biodbnet.abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&input=';
#my $url = 'http://biodbnet.abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&input=genesymbolandsynonyms&inputValues=KIT,PDGFRA&outputs=refseqmrnaaccession&taxonId=9606&format=row';
if(open (INPUT,$inputfile)){
    my @all=<INPUT>;
    close INPUT;
    my $counter=1;
    my $temp="";
    while(@all){
        $line =shift(@all);
        chomp ($line);
        if ($counter>=400){
            $temp=$temp.','.$line;
            #call and excute one rest
            &getrest($temp);

            $counter=1;
            $temp="";
        } else{
            $counter++;
            $temp=$temp.','.$line;
        }#else


    }#while(@all)

}
print Dumper \%target_hash;
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
        %output=%{$current};
        if ($output{"RefSeq mRNA Accession"} eq "-"){
            next;        
        }#this is empty output
        else{
            my $currentInput=$output{"InputValue"};
            my @refseqs=split("\\",$output{"RefSeq mRNA Accession"});
            print @refseqs;
            %target_hash{$currentInput}=@refseqs;            
        }
        
    
    }
}



