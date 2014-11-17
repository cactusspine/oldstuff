#!/usr/bin/perl -w

use File::Basename;
use warnings;
use utf8;
use strict;
use Storable;
use Data::Dumper qw(Dumper);
use List::MoreUtils 'uniq';
use Time::localtime;
###########################################################
my $project_dir = "/home/jiali/workspace/projects/vs/";
my $data_dir = $project_dir . "data/";
my $result_dir = $project_dir . "result/";
my %all_resulted_uni_ids_hash = ();
my $result_file = $result_dir .&timestamp()."mirna_targets.txt";
close(RESULT_FILE);
open(RESULT_FILE,">>$result_file");
#close(RESULT_FILE);
my $tax_ID="9606";
my $line = "";
my @all = ();
###################################################################
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
###############################################################
#### comparative analysis ............#############################################
my @hash_files= glob "$data_dir/*/*.hash";
my @hash_to_check;
my @hash_to_check_name;
my %weight_hash;
my $counter=0;
print RESULT_FILE "MIRBASE_ACCESSION\tREFSEQ_ID\t";
for(0..$#hash_files){
    if($hash_files[$_]!~/annotations|miRBase/){
        push @hash_to_check,retrieve($hash_files[$_]);#@hash_to_check contains the reference of each hash
        
        print "$hash_files[$_]\n";
        my $basename=basename($hash_files[$_]);
        (my $title=$basename)=~s/\.hash//;
        print RESULT_FILE "$title\t";
        print "$basename\t";
        push @hash_to_check_name,$basename;
        if($hash_files[$_]=~/tarbase/i){
            $weight_hash{$counter}=3;
            print "$basename  weight =".$weight_hash{$counter}."\n";
        }else{
            $weight_hash{$counter}=1;
            print "$basename weight = $weight_hash{$counter}\n";
        }
        $counter++;
    }
}

print RESULT_FILE "priority\n";
####################################################

=begin  BlockComment  # BlockCommentNo_1

my %mirnaAcc2mirAlias_hash=%{retrieve($data_dir.'miRBase/mir_alias.hash')};
my @all_result_mirna=keys %mirnaAcc2mirAlias_hash;
print "the number of mirnaID:".scalar(@all_result_mirna)."\n";
my %union_refseq;

my %refseq2ensg_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_refseq2ensg.hash')};
my %refseq2entrez_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_refseq2enterz.hash')};
my %refseq2ensembl_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_refseqID2ensembl_transcriptID.hash')};
my %refseq2geneName_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_refseq2rnaGeneNames.hash')};
foreach my $key (keys %refseq2entrez_hash){
    $union_refseq{$key}=1;
}

foreach my $key (keys %refseq2ensg_hash){
    $union_refseq{$key}=1;
}

foreach my $key (keys %refseq2ensg_hash){
    $union_refseq{$key}=1;
}

foreach my $key (keys %refseq2geneName_hash){
    $union_refseq{$key}=1;
}
my @all_refseq=keys %union_refseq;
print "size of refeq array".scalar(@all_refseq)."\n";
my %result_hash;

=end    BlockComment  # BlockCommentNo_1

=cut

#########################################################
my %all_result_mirna;
my %union_refseq;
#my $counter=0;
my %result_hash;
foreach my $ref_hash(@hash_to_check){
    
    #print "$ref_hash\n";
    my %current_hash=%{$ref_hash};
    # print Dumper \%current_hash;
    foreach my $mirna_acc(sort keys %current_hash){
        foreach my $refseq_acc(keys %{$current_hash{$mirna_acc}}){
            #print "$mirna_acc,$refseq_acc\n";
            if($mirna_acc=~'ARRAY' or $refseq_acc=~'ARRAY'){print "the $hash_to_check_name[$counter] is wrong\n";}

            if($refseq_acc=~'ARRAY' or $refseq_acc=~'ENST'){print "the $hash_to_check_name[$counter] refseq containing ENST\n";}
            
            if($mirna_acc=~'ARRAY' or $mirna_acc=~'ENST'){print "the $hash_to_check_name[$counter] mirna_acc is containing ENST\n";}
            $all_result_mirna{$mirna_acc}=1;
            $union_refseq{$refseq_acc}=1;
        }    
    }
#$counter++;
}
#print Dumper \%all_result_mirna;

#print Dumper \%union_refseq;
my @all_result_mirna=keys %all_result_mirna;
print "the number of mirnaID:".scalar(@all_result_mirna)."\n";
my @all_refseq=keys %union_refseq;
print "size of refeq array".scalar(@all_refseq)."\n";
sleep(30);
$counter=0;#reset counter
foreach my $mirbase_id(@all_result_mirna){

    foreach my $refseq_id(@all_refseq){
        my $currentline="$mirbase_id\t$refseq_id";
        my $at_least_one_target=0;
        my $priority=0;
        

        for my $ref_hash(@hash_to_check){
            #print "$ref_hash\n";
            my %current_hash=%{$ref_hash};
            if (exists $current_hash{$mirbase_id}{$refseq_id}){
                $currentline.="\tyes";
                $priority=$priority+$weight_hash{$counter};
                $at_least_one_target=1;
            }else{
                $currentline.="\tno";
            }
            $counter++;
        }
        print "$currentline\t$priority\n";
        $counter=0;#reset counter before checking next miRNA&refseq pair after gone through all the hashes
        #sleep(5);
        if($at_least_one_target==1){
            print RESULT_FILE "$currentline\t$priority\n";
            $result_hash{$mirbase_id}{$refseq_id}=$priority;
            
        }        
    
    }

}
close(RESULT_FILE);
my $result_hash_file=$result_dir.$tax_ID.'result.hash';
store\%result_hash,$result_hash_file;

