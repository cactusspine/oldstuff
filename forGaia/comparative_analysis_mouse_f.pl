#!/usr/bin/perl -w

use File::Basename;
use warnings;
use utf8;
use strict;
use Storable;
use Data::Dumper qw(Dumper);
use List::MoreUtils 'uniq';
use Time::localtime;
use v5.16;
use feature "switch";
###########################################################
my $project_dir = "/work/projects/pdgfr_kit/miRNA/wj/";
my $data_dir = $project_dir . "data/";
#my $data_dir=$project_dir."data/trail/";
my $result_dir = $project_dir . "result/";
my %all_resulted_uni_ids_hash = ();
my $result_file = $result_dir .&timestamp()."mirna_targets.txt";
my $result_file_d = $result_file."_detail";
close(RESULT_FILE);
open(RESULT_FILE,">>$result_file");
open (RESULT_FILE_D,">>$result_file_d");
#close(RESULT_FILE);
my $tax_ID="10090";
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
#my @hash_files= glob "$data_dir/*/*.hash";
my @hash_files= glob "$data_dir/*/*10090*.hash";
my @hash_to_check;
my @hash_to_check_name;
my %weight_hash;
my $counter=0;
print RESULT_FILE "MIRBASE_ACCESSION\tREFSEQ_ID\t";
for(0..$#hash_files){
    if($hash_files[$_]!~/annotations|miRBase/){
        push @hash_to_check,retrieve($hash_files[$_]);#@hash_to_check contains the reference of each hash
        
        #print "$hash_files[$_]\n";
        my $basename= basename($hash_files[$_]);
        (my $title=$basename)=~s/\.hash//;
        print RESULT_FILE "$title\t";
        print "$basename\t";
        print RESULT_FILE_D "$title\t";
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
my $counter=0;
my %result_hash;
foreach my $ref_hash(@hash_to_check){
    
    #print "$ref_hash\n";
    my %current_hash=%{$ref_hash};
    # print Dumper \%current_hash;
    print "current underchecking hash is $hash_to_check_name[$counter]\n";
    my $current_mirna_nr=scalar keys %current_hash;
    print "current hash contrains miRNA accession Nr is $current_mirna_nr\n";
    foreach my $mirna_acc(sort keys %current_hash){
        my $current_target_nr = scalar keys %{$current_hash{$mirna_acc}};
        print "current mirna $mirna_acc has $current_target_nr targets\n";
	
        if($mirna_acc=~'MIMAT'){$all_result_mirna{$mirna_acc}=1;}else{print "it is a wrong entry\n";next;}
        foreach my $refseq_acc(keys %{$current_hash{$mirna_acc}}){
            #print "$mirna_acc,$refseq_acc\n";
            if($mirna_acc=~'ARRAY' or $refseq_acc=~'ARRAY'){print "the $hash_to_check_name[$counter] is wrong\n";next;}

            if($refseq_acc=~'ARRAY' or $refseq_acc=~'ENST'){print "the $hash_to_check_name[$counter] refseq containing ENST\n";next;}
            
            if($mirna_acc=~'ARRAY' or $mirna_acc=~'ENST'){print "the $hash_to_check_name[$counter] mirna_acc is containing ENST\n";next;}
            $union_refseq{$refseq_acc}=1;
            $result_hash{$mirna_acc}{$refseq_acc}{'priority'}=0;


        }    
    }
$counter++;
}
#print Dumper \%all_result_mirna;

#print Dumper \%union_refseq;
my @all_result_mirna=keys %all_result_mirna;
print "the number of mirnaID:".scalar(@all_result_mirna)."\n";
my @all_refseq=keys %union_refseq;
print "size of refeq array".scalar(@all_refseq)."\n";
sleep(30);
##################the below part need improvement###############################################################################################
=begin
my %result_detail_hash;
foreach my $mirbase_id(sort keys %result_hash){
    my $counter_target_miR=0;#set counter of number of targets for this miR

    foreach my $refseq_id(keys %{$result_hash{$mirbase_id}}){
        my $currentline="$mirbase_id\t$refseq_id";
        my $at_least_one_target=0;
        my $priority=0;       
        $counter=0;#reset counter
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
        #print "$currentline\t$priority\n";
        #$counter=0;#reset counter before checking next miRNA&refseq pair after gone through all the hashes
        #sleep(5);
        if($at_least_one_target==1){
            print RESULT_FILE "$currentline\t$priority\n";
            $result_detail_hash{$mirbase_id}{$refseq_id}{"priority"}=$priority;
            $counter_target_miR=1;

            
        }
        
    
    }
if($counter_target_miR==0){
    print "there is something wrong\n";

}
print "$mirbase_id has already been checked\n";
}
close(RESULT_FILE);
my $result_hash_file=$result_dir.$tax_ID.'result.hash';
store\%result_detail_hash,$result_hash_file;
=end

=cut

#################################################################################################################
my %result_detail_hash;
$counter=0;#this is the counter to keep track of ref_hashes
for my $ref_hash(@hash_to_check){
    print "current underchecking hash is $hash_to_check_name[$counter]\n";
    my %current_hash=%{$ref_hash};
    foreach my $mirbase_id(sort keys %result_hash){
       
        foreach my $refseq_id(keys %{$result_hash{$mirbase_id}}){
            #my $currentline="$mirbase_id\t$refseq_id";
            if(exists $current_hash{$mirbase_id}{$refseq_id}){
                $result_hash{$mirbase_id}{$refseq_id}{'priority'}=$result_hash{$mirbase_id}{$refseq_id}{'priority'}+$weight_hash{$counter};#update the priority
                $result_hash{$mirbase_id}{$refseq_id}{"summary"}[$counter]="yes"; 
            ###########################################################################################################
            #this part was added to keep track of the score and evidences ,will be store in {'details'}
            given($hash_to_check_name[$counter]){
            when (/eimmo/) {$result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]=$current_hash{$mirbase_id}{$refseq_id} ;print "found $hash_to_check_name[$counter]\n";}#eimmo P value, bigger better
            when (/microcosm/) {$result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]=$current_hash{$mirbase_id}{$refseq_id}{"score"} ;print "found $hash_to_check_name[$counter]\n";}#microcosm score , bigger better
            when (/targetScan/) {$result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]=$current_hash{$mirbase_id}{$refseq_id}{"total_Context_Score"}.';'.$current_hash{$mirbase_id}{$refseq_id}{"aggregatePCT"} ;print "found $hash_to_check_name[$counter]\n";}#eimmo P value, bigger better
            when (/_pictar_/) {$result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]=$current_hash{$mirbase_id}{$refseq_id}{'pictar_score'} ;print "found $hash_to_check_name[$counter]\n";}#pictar score, higher better
            when (/_pita_/) {$result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]=$current_hash{$mirbase_id}{$refseq_id}{'pita_score_per_gene'} ;print "found $hash_to_check_name[$counter]\n";}#pita value, smaller better
            when (/microrna_org/) {$result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]=$current_hash{$mirbase_id}{$refseq_id}{"miRanda_align_score"}.';'.$current_hash{$mirbase_id}{$refseq_id}{"misvr_score"} ;print "found $hash_to_check_name[$counter]\n";}
            when (/mirdb/) {$result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]=$current_hash{$mirbase_id}{$refseq_id} ;print "found $hash_to_check_name[$counter]\n";}#mirdb score, bigger better
            when (/microT/) {$result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]=$current_hash{$mirbase_id}{$refseq_id} ;print "found $hash_to_check_name[$counter]\n";}#eimmo P value, bigger better
            when (/_mirtarbase/) {$result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]=$current_hash{$mirbase_id}{$refseq_id} ;print "found $hash_to_check_name[$counter]\n";}#targetgene:experiment:reference PMID
            when (/targetspy/) {$result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]=$current_hash{$mirbase_id}{$refseq_id}{"targetspyScore"} ;print "found $hash_to_check_name[$counter]\n";}#eimmo P value, bigger better
            when (/_paccmit\.hash/) {$result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]=$current_hash{$mirbase_id}{$refseq_id}{'logP'} ;print "found $hash_to_check_name[$counter]\n";}#

            when (/paccmit_cds/) {$result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]=$current_hash{$mirbase_id}{$refseq_id}{'P_SH'} ;print "found $hash_to_check_name[$counter]\n";}#
            when (/_tarbase/) {$result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]=$current_hash{$mirbase_id}{$refseq_id}{'evidence'}.':'.$current_hash{$mirbase_id}{$refseq_id}{'regulation'} ;print "found $hash_to_check_name[$counter]\n";}#eimmo P value, bigger better
            default{
            $result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]="yes";#eimmo P value, bigger better
        }

            }

            ##########################################################################################################
            }#if (exists $current_hash{$mirbase_id}{$refseq_id}){this line exist in the current under checking hash
            else{#this pair does not exist in the current under checking hash
                $result_hash{$mirbase_id}{$refseq_id}{"summary"}[$counter]="no"; 
                $result_hash{$mirbase_id}{$refseq_id}{'details'}[$counter]="no";
            }#else{#this pair does not exist in the current under checking hash
        }# foreach my $refseq_id(keys %{$result_hash{$mirbase_id}}){
    }#foreach my $mirbase_id(sort keys %result_hash)
    $counter++;#move on to next hash
}

my $result_hash_file = $result_dir.$tax_ID.'result.hash';
%result_detail_hash = %result_hash;
store\%result_detail_hash,$result_hash_file;

foreach my $mirbase_id(sort keys %result_hash){
    my $counter_target_miR=0;#set counter of number of targets for this miR
    foreach my $refseq_id(keys %{$result_hash{$mirbase_id}}){
        my $currentline="$mirbase_id\t$refseq_id";
        my $currentline_d="$mirbase_id\t$refseq_id";
        my $at_least_one_target=0;
        $counter_target_miR=1;
        foreach my $current_stat(@{$result_hash{$mirbase_id}{$refseq_id}{"summary"}}){
            $currentline=$currentline."\t".$current_stat;
        }#foreach my $current_stat(@{$result_hash{$mirbase_id}{$refseq_id}{"summary"}})
        
        foreach my $current_sta(@{$result_hash{$mirbase_id}{$refseq_id}{'details'}}){
            $currentline_d=$currentline_d."\t".$current_sta;
        }#foreach my $current_stat(@{$result_hash{$mirbase_id}{$refseq_id}{"summary"}})
        $currentline = $currentline."\t".$result_hash{$mirbase_id}{$refseq_id}{'priority'};
        print RESULT_FILE "$currentline\n";
        print RESULT_FILE_D "$currentline_d\n";
        print "$currentline\n";
        }#foreach my $refseq_id(keys %{$result_hash{$mirbase_id}}){
    
        if($counter_target_miR==0){
            print "there is something wrong\n";
        }
print "$mirbase_id has already been checked\n";
}

close(RESULT_FILE);
close(RESULT_FILE_D);

