#!/usr/bin/perl -w
#===============================================================================
#
#         FILE: mirna2target_wj_v1.pl
#
#        USAGE: ./mirna2target_wj_v1.pl  
#  DESCRIPTION: this the script to process all the different miRNA-target prediction databases in the /work/projects/pdgfr_kit/miRNA/
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS:
#        NOTES: ---
#       AUTHOR: Jiali Wang, 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 04/24/2014 01:27:32 PM
#     REVISION:
#===============================================================================

use feature "switch";
use warnings;
use utf8;
use strict;
use Storable;
use Data::Dumper qw(Dumper);
use List::MoreUtils 'uniq';
use Time::localtime;
###########################################################
my $project_dir = "/work/projects/pdgfr_kit/miRNA/wj/";
my $data_dir = $project_dir . "data/";
my $result_dir = $project_dir . "result/";
my %all_resulted_uni_ids_hash = ();
my $result_file = $result_dir .&timestamp()."mirna_targets.txt";
#close(RESULT_FILE);

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
###miRBASE alias matching related==miRBASE alias===#########################################################################################
sub mirBaseAlias(){
    my $miRBase_alias_file=$data_dir."miRBase/aliases.txt";
    close(MIRALIAS);
    my %mir_alias_hash=();#MIRBASE accession as key
    my %mirAlias2Accession_hash=();#all the documented mir alias as key to the accession
    if (open (MIRALIAS,$miRBase_alias_file)){
        @all=<MIRALIAS>;
        close(MIRALIAS);
        while(@all){
            $line = shift @all;
            chomp ($line);#every line looks like :MI0000001	cel-let-7L;cel-let-7;
            chop ($line);#remove the last ;
            my ($mirbase_accession,$mir_alias_list)=split("\t",$line);
            my @mir_aliases = split(";",$mir_alias_list);
            @{$mir_alias_hash{$mirbase_accession}}=@mir_aliases;
            foreach my $i(@mir_aliases){
                $mirAlias2Accession_hash{$i}=$mirbase_accession;
            }
        }
    }
    #print $mir_alias_hash{"MI0000001"};#test
    #print $mirAlias2Accession_hash{"cel-let-7"};#test
    store \%mir_alias_hash, $data_dir.'miRBase/mir_alias.hash';
    store \%mirAlias2Accession_hash, $data_dir.'miRBase/mirAlias2Accession.hash';
}#mirBaseAlias()
##END===SECTION_miRbase_Alias==################################################################################################################## 

##################################################################################################################################################
###################################################################################
#NM_001207023	mRNA linear	IQCF3
#NM_001164434	mRNA linear	KRTAP22-2 KRTAP22 2
#this file contains all the RefSeq(which are not DNAs by grep) and their corresponding HGNC symbols and RNA type
#attention! symbols which containing - and . need special treatment to eliminate their pieces
############################################################
sub refseq2rnaGeneName(){
    my ($species) = @_;#$species='human','mouse';
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my $refseq2rna_geneName_file = $data_dir.'annotations/'.$tax_ID.'_refseq2RNA_GeneNames.txt';
    print "try read in refseq2rnaGeneNames FILE:$refseq2rna_geneName_file\n";
    if(open(REFSEQ2RNAGENENAME,$refseq2rna_geneName_file)){
        @all = <REFSEQ2RNAGENENAME>;
        close(REFSEQ2RNAGENENAME);
    }else{print "Can't open refseq2rnaGeneName file : $refseq2rna_geneName_file\n";}
    my %rnaGeneName2refseqs_hash=();#the hash containing key:rnaGeneName=>@non redundent refseqs ids
    my %refseqs2rnaGeneNames_hash=();#the hash containing refseqs=>@ non redundent rnaGeneNames(still possible one transcript multiple miRNAs)
    while(@all){
        $line = shift(@all);
        chomp($line);
        if($line){
            my($refseq_ID,$molecule_type,$geneSymbols_string)= split("\t",$line);          
            my @geneSymbols_temp= uniq(split(" ",$geneSymbols_string));#it will print a array even without token
            my @geneSymbols_temp_detailed= split(/[.|\-| ]/,$geneSymbols_string);#split based on multiple delimer including _. and \s
            if (scalar(@geneSymbols_temp)<scalar(@geneSymbols_temp_detailed)){
                my %temp=();
                foreach (@geneSymbols_temp_detailed){                                               
                    $temp{$_}++;                    
                }#  foreach (@geneSymbols_temp_detailed){                
                foreach(@geneSymbols_temp){
                    if((exists $temp{$_}) and $temp{$_}>1){#this key only showed once in the detailed hash. eg. fami.ly fami ly=>fami.ly,fami,ly=detailed>fami,ly,fami,ly
                        next;#this is a entity need to be deleted
                    }else{
                        push @{$refseqs2rnaGeneNames_hash{$refseq_ID}},$_;
                        push @{$rnaGeneName2refseqs_hash{$_}},$refseq_ID;
                    }
                }#foreach(@geneSymbols_temp){
            }#if (scalar(@geneSymbols_temp)<scalar(@geneSymbols_temp_detailed)){#strategy:counting the duplicated string with ._delimer
            else{#there is no hyphen, dot or exception
                push @{$refseqs2rnaGeneNames_hash{$refseq_ID}},@geneSymbols_temp;
                foreach(@geneSymbols_temp){
                    push @{$rnaGeneName2refseqs_hash{$_}},$refseq_ID;
                }# foreach(@geneSymbols_temp){
            }#else#there is no hyphen, dot or exception
        }#if($line){
    }#while(@all){
    my $rnaGeneName2refseqs_hash_file=$data_dir.'annotations/'.$tax_ID.'_rnaGeneName2refseqs.hash';
    my $refseq2_rnaGeneNames_hash_file=$data_dir.'annotations/'.$tax_ID.'_refseq2rnaGeneNames.hash';
    print Dumper \%refseqs2rnaGeneNames_hash;
    print Dumper \%rnaGeneName2refseqs_hash;
    store \%rnaGeneName2refseqs_hash,$rnaGeneName2refseqs_hash_file;     
    store \%refseqs2rnaGeneNames_hash,$refseq2_rnaGeneNames_hash_file;
    return ($rnaGeneName2refseqs_hash_file,$refseq2_rnaGeneNames_hash_file); 
}#sub refseq2rnaGeneName   
##==END_ABOUT_REFSEQ2RNA_GENE_NAMES#######################################################################################

##==ABOUT REFSEQ2ENSMBL== ###########################################################################################
#NM_032333	ENST00000372187
#NM_015900	ENST00000273371
#this file have a one to one relationship
############################################################
sub refseq2ensembl(){
    my ($species) = @_;#$species='human','mouse';
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my $refseq2ensembl_tr_file = $data_dir.'annotations/'.$tax_ID.'_refseq2ensembl_transcript.txt';
    print "try read in refseq2ensembl FILE:$refseq2ensembl_tr_file\n";
    if(open(REFSEQ2ENSEMBL,$refseq2ensembl_tr_file)){
        @all = <REFSEQ2ENSEMBL>;
        close(REFSEQ2ENSEMBL);
    }else{print "Can't open refseq2ensembl file : $refseq2ensembl_tr_file\n";}
    my %ensembl_tr2refseq_hash=();#the hash containing key:rna ensembl transcript ID=>$refseqs id
    my %refseq2ensembl_tr_hash=();#
    while(@all){
        $line= shift(@all);
        chomp($line);
        if($line){
            my($refseq_ID,$ensembl_tr)= split("\t",$line);
            $ensembl_tr2refseq_hash{$ensembl_tr}=$refseq_ID;
            $refseq2ensembl_tr_hash{$refseq_ID}=$ensembl_tr;        
        }#if($line){
    }#while(@all){



    my $ensembl_tr2refseq_hash_file=$data_dir.'annotations/'.$tax_ID.'_ensembl_transcriptID2refseqID.hash';
    my $refseq2ensembl_tr_hash_file=$data_dir.'annotations/'.$tax_ID.'_refseqID2ensembl_transcriptID.hash';
    store \%ensembl_tr2refseq_hash,$ensembl_tr2refseq_hash_file;     
    store \%refseq2ensembl_tr_hash,$refseq2ensembl_tr_hash_file;
    return ($ensembl_tr2refseq_hash_file,$refseq2ensembl_tr_hash_file); 
}#sub refseq2ensembl
##==END_ABOUT_REFSEQ2ENSEMBL#######################################################################################

##==ABOUT biomart_ensemblENSG2refseq== ###########################################################################################
#this file was created from ensembbiomart servicehttp://www.ensembl.org/biomart/martview/
#Dataset Genome
#no filter
#attributes: ensembl gene ID;Ensembl transcript ID; Refseq miRNA predicted ;refseqmRNA NM_001195597
#attention the downloaded file need treatment some time the order of tab can change so first change the order of column by awk
#awk -F "\t" '{print $4,$3,$1,$2}' OFS="\t" mouse_ENSG2refseq.txt>mouse_ENSG2refseq.txt_2
#then delete those lines without matching by grep
#grep -v "^\t\t" mouse_ENSG2refseq.txt_2>mouse_ENSG2refseq.txt_final
#by the way check if there is any difference between refseq SRS result and biomart add not added relationship onto the file
############################################################
sub ensg2refseqs(){
    my ($species) = @_;#$species='human','mouse';
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my $ensg2refseq_file = $data_dir.'annotations/'.$species.'_ENSG2refseq.txt_final';
    print "try read in ensg2refseqs FILE:$ensg2refseq_file\n";
    if(open(ENSG2REFSEQ,$ensg2refseq_file)){
        @all = <ENSG2REFSEQ>;
        shift(@all);#this file have head line
        close(ENSG2REFSEQ);
    }else{print "Can't open ensembl geneID 2 refseq file : $ensg2refseq_file\n";}
    ###############this part is newly added ###################################
    my $ensembl_tr2refseq_hash_file=$data_dir.'annotations/'.$tax_ID.'_ensembl_transcriptID2refseqID.hash';
    my $refseq2ensembl_tr_hash_file=$data_dir.'annotations/'.$tax_ID.'_refseqID2ensembl_transcriptID.hash';
    my %refseq2ensembl_tr_hash;
    my %ensembl_tr2refseq_hash;
    while(not -e $ensembl_tr2refseq_hash_file){
        &refseq2ensembl($species);
    }
    while(not -e $refseq2ensembl_tr_hash_file){
        &refseq2ensembl($species);
    }    
    %ensembl_tr2refseq_hash=%{retrieve($ensembl_tr2refseq_hash_file)};
    %refseq2ensembl_tr_hash=%{retrieve($refseq2ensembl_tr_hash_file)};  
    #################this part is newly added################################
    my %ensg2refseqs_hash=();#the hash containing key:ensembl gene IDs=>@non redundent refseqs ids
    my %refseq2ensg_hash=();#the hash containing refseqs=> ensembl gene IDs
    while(@all){
        $line = shift(@all);
        chomp($line);
        if($line){
            #$RefSeq mRNA [e.g. NM_001195597]	$RefSeq mRNA predicted [e.g. XM_001125684]	$Ensembl Gene ID	$Ensembl Transcript ID
            my($refseq_mRNA,$refseq_predicted,$ensg,$enst)= split("\t",$line);            
            if (!$refseq_mRNA eq ''){
                push @{$ensg2refseqs_hash{$ensg}},$refseq_mRNA;
                $refseq2ensg_hash{$refseq_mRNA}=$ensg;
            }
            if (!$refseq_predicted eq ''){
                push @{$ensg2refseqs_hash{$ensg}},$refseq_predicted;
                $refseq2ensg_hash{$refseq_predicted}=$ensg;
            }
            #######newly added #########################################
            if (not exists $ensembl_tr2refseq_hash{$enst}){
                print "$enst has not been mapped ! :$refseq_mRNA: $refseq_predicted\n";
                #rule: if there are two match, then chose the$refseq_mRNA,ditch the predicted 
                if(!$refseq_mRNA eq ''){
                    $ensembl_tr2refseq_hash{$enst}=$refseq_mRNA;
                    $refseq2ensembl_tr_hash{$refseq_mRNA}=$enst;
                }else{
                    $ensembl_tr2refseq_hash{$enst}=$refseq_predicted;
                    $refseq2ensembl_tr_hash{$refseq_predicted}=$enst;
                }
            }
            ################newly added_END###################################
            @{$ensg2refseqs_hash{$ensg}} = uniq(@{$ensg2refseqs_hash{$ensg}});  
            #######################this part is newly added 
        }#if($line){
    }#while(@all){
    my $ensg2refseqs_hash_file=$data_dir.'annotations/'.$tax_ID.'_ensg2refseqs.hash';
    my $refseq2ensg_hash_file=$data_dir.'annotations/'.$tax_ID.'_refseq2ensg.hash';
    print Dumper \%ensg2refseqs_hash;
    store \%ensg2refseqs_hash,$ensg2refseqs_hash_file;     
    store \%refseq2ensg_hash,$refseq2ensg_hash_file;
########newly added##################################################
    store \%ensembl_tr2refseq_hash,$ensembl_tr2refseq_hash_file;     
    store \%refseq2ensembl_tr_hash,$refseq2ensembl_tr_hash_file;
###########newly added################################################
    return ($ensg2refseqs_hash_file,$refseq2ensg_hash_file); 
}#sub ensg2refseqs   
##==END_ABOUT_REFSEQ2ENSEMBLGENEs#######################################################################################

#########refseq2ENTREZ############################################################
sub refseq2entrez(){
    my ($species) = @_;#$species='human','mouse';
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my $refseq2enterzgene_file = $data_dir.'annotations/'.$tax_ID.'.Refseq2enterzgene_id.txt';
    print "try read in refseq2entrez FILE:$refseq2enterzgene_file\n";
    if(open(REFSEQ2ENTERZ,$refseq2enterzgene_file)){
        @all = <REFSEQ2ENTERZ>;
        close(REFSEQ2ENTERZ);
    }else{print "Can't open refseq2entrez file : $refseq2enterzgene_file\n";}
    my %enterzgene2refseqs_hash=();#the hash containing key:enterz=>@non redundent refseqs ids
    my %refseq2enterzgene_hash=();#the hash containing refseq=>enterz
    while(@all){
        $line = shift(@all);
        chomp($line);
        if($line){     
            my($refseq_ID,$enterz)= split("\t",$line);
            push @{$enterzgene2refseqs_hash{$enterz}},$refseq_ID;
            $refseq2enterzgene_hash{$refseq_ID}=$enterz;
        }#if($line){
    }#while(@all){
    my $enterz2refseqs_hash_file = $data_dir.'annotations/'.$tax_ID.'_enterz2refseqs.hash';
    my $refseq2enterz_hash_file = $data_dir.'annotations/'.$tax_ID.'_refseq2enterz.hash';
    print Dumper \%refseq2enterzgene_hash;
    print Dumper \%enterzgene2refseqs_hash;
    store \%enterzgene2refseqs_hash,$enterz2refseqs_hash_file;     
    store \%refseq2enterzgene_hash,$refseq2enterz_hash_file;
    return ($enterz2refseqs_hash_file,$refseq2enterz_hash_file); 
}#sub refseq2enterz
##==END_ABOUT_REFSEQENTERZ#######################################################################################

### ==BEGIN_MICROCOSM==#########################
#&microcosm('human');
#&microcosm('mouse');
#Microcosm_V5 no updating is using miRanda algorithem to identify potential BS find highly complementary at 5 end(0-100)->Vienna routine for thermodynamic stability->check ortholog conservation.
#About MicroCosm Targetshttp://www.ebi.ac.uk/enright-srv/microcosm/htdocs/targets/v5/info.html
#this downloading version was created on 2007-11-1 and there is no update ever since
##this scoring seems don't take into consider the multiple BS effect it simplely chose the highest scored BS
sub microcosm(){ #input the name of file to be parsed
    #my ($microcosm_mirna2target_file) =shift(@_);   #
    my $species = shift(@_);
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my $miRNA_Family_file = $data_dir."targetscan/".$tax_ID."_miR_Family_Info.txt";
    my $microcosm_mirna2target_file=$data_dir."microcosm/".$species.'_microcosm_mirna2target.txt';#    
    my %ens_transcript2uni_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_ensembl_transcriptID2refseqID.hash')};
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    my $missing_GeneName_file = $data_dir.'microcosm/microcosm_'.$species.'_missingENSEMBL.txt';
    close(MICROCOSM_FILE);
    my %microcosm_hash = ();
    my %microcosm_detail_hash = ();
    if (open(MICROCOSM_FILE, $microcosm_mirna2target_file)){
        @all = <MICROCOSM_FILE>;
        close(MICROCOSM_FILE);
        $line = shift @all; # header info deleted
    }# if (open(MICROCOSM_FILE, $mouse_microcosm_mirna2target_file)){
    else{
    print"Can't open file $microcosm_mirna2target_file\n";
    }
    while(@all){
        $line = shift @all;
        chomp($line);
        if ($line){            
            ##GROUP SEQ METHOD  FEATURE CHR START END STRAND  PHASE SCORE PVALUE_OG TRANSCRIPT_ID\ENSEMBLTRANSCRIPT EXTERNAL_NAME\GENESYMBOL	
            my ($group, $mirna_id, $method,  $feature, $chr, $start, $end, $strand,  $phase, $score, $pvalue_og, $ensembl_transcript_id, $external_name) = split("\t", $line);
            #$mirna_id=~s/\*|star|_star//;$keep mirna_id as it is ...it is the old mir base ID like hsa-miR-497*, it could be matched by the miRbase dictionary
            my @uni_id = (); #use uni_id as uni id
            if ($ens_transcript2uni_hash{$ensembl_transcript_id}) {#use ensembl_transcript_id for mapping
                my $tmp_uni_id = $ens_transcript2uni_hash{$ensembl_transcript_id};
                push @uni_id, $tmp_uni_id;
                }
                else{
                    print"can't find corresponding UNI_ID for $ensembl_transcript_id\n";
                    open (MISSING_GENE,'>>'.$missing_GeneName_file)||die "Could not open file $!";
                    print MISSING_GENE "$ensembl_transcript_id\n";
                    close(MISSING_GENE);
                    }#else{
            #if ($genename2uni_hash{$external_name}) {#====to revise===
            #my $tmp_uni_id = $genename2uni_hash{$external_name};
            #push @uni_id, $tmp_uni_id;
            # }#this "if" might be omitted since every ensembl_transcript could be better matched
            if (@uni_id) {#if there is id which could be matched at this line
                @uni_id = &nonRedundantList(\@uni_id);	#tide up in case of multiple match one ensemblTranscript->multiple uni_ID
            foreach my $uni_id (@uni_id) {
                $all_resulted_uni_ids_hash{$uni_id}=1;# this uniID gene have match
                if ($mirAlias2Accession_hash{$mirna_id}){
                    my $mirna_accession = $mirAlias2Accession_hash{$mirna_id};
                    if ($microcosm_hash{$mirna_id}{$uni_id}) {#if there is already a entry with same mirna_accession and uni_id update it,else do nothing # the higher the better?
                        my $old_score = $microcosm_detail_hash{$mirna_accession}{$uni_id}{"score"};#the microcosm score seem will keep the highest for multiple BS within the same pair
                        if ($score > $old_score) {
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"chr"} = $chr;
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"start"} = $start;
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"end"} = $end;
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"strand"} = $strand;
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"score"} = $score;
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"pvalue_og"} = $pvalue_og;
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"gene_symbol"} = $external_name;
                            # $microcosm_hash{$mirna_accession}{$uni_id} = 1;
                            $microcosm_hash{$mirna_accession}{$uni_id} = 1;
                        }# if ($score > $old_score) {#if there is already a entry with same mirna_accession and uni_id update it,else do nothing # the higher the better?
                    } # if ($microcosm_hash{$mirna_accession}{$uni_id}) {
                    else {
                        $microcosm_detail_hash{$mirna_accession}{$uni_id}{"chr"} = $chr;
                        $microcosm_detail_hash{$mirna_accession}{$uni_id}{"start"} = $start;
                        $microcosm_detail_hash{$mirna_accession}{$uni_id}{"end"} = $end;
                        $microcosm_detail_hash{$mirna_accession}{$uni_id}{"strand"} = $strand;
                        $microcosm_detail_hash{$mirna_accession}{$uni_id}{"score"} = $score;
                        $microcosm_detail_hash{$mirna_accession}{$uni_id}{"pvalue_og"} = $pvalue_og;
                        $microcosm_detail_hash{$mirna_accession}{$uni_id}{"external_name"} = $external_name;
                        $microcosm_hash{$mirna_accession}{$uni_id} = 1;
                        #$microcosm_hash{$uni_id}{$mirna_accession} = 1;
                    } # else {# if there is no entry of this miRNA and uni_id pair yet
                }#if this mir_id find match in mir_dictionary if ($mirAlias2Accession_has{$mirna_id}){                   
            } # foreach my $uni_id (@uni_id) { 
         } # if (@uni_id) {
        else { print "problem - microcosm - no uni_id - line : $line\n";} # else {#no uni_id match for this gene
    } # if ($line) {
}# while(@all){
#my $microcosm_detail_hash_file=$data_dir."microcosm/".$tax_ID.'_microcosm_mirna2target.detail_hash';    
#store \%microcosm_detail_hash,$microcosm_detail_hash_file;
my $microcosm_hash_file=$data_dir."microcosm/".$tax_ID.'_microcosm_mirna2target.hash';
store \%microcosm_detail_hash,$microcosm_hash_file;
return ($microcosm_hash_file);
}
#genes could be matched to ensembl_transcript id

###END==microcosm==##########################################
##############################################################

###==BEGIN_TARGETSCAN==############################################
#it might need a subroutine to link the stem loop pre miRNA to the mature type to extend the matching of miRNA(don't need to ,it's website doesn't do this either)
##need to a dictionary between the miRfamily&MIRbaseID(Attention, this file contains the relationship between miRfamily vs all the mature miR)#ATTENTION:MATURE_MIR ONLY
#better to choose thresholds for targets (total context scores or total Pct) rather than thresholds for sites, because multiple weak sites to the same miRNA can add up to more repression than a single strong site.
##use summary file################################################################
##push all the miRNAs in the miR family into the hash
#target_Scan_hash{RefseqTranscriptID}{miRNAAccession}{"Context_Score"}
#---------------------------------------------{"Aggregate_PCT"}
#---------------------------------------------{"miRNAfamily"} @ push in the name of miRNA names
#standard: Species id =10090(Mouse) or 9606(Human)
#Total num conserved sites>1 need to be recoded
#Total context score <-0.3
#Agggegate PCT>0.8
#NM_001001130	Zfp85-rs1	AAAUUCG	10116	0	0	0	0	1	0	0	1	rno-miR-10a-3p	-0.011	NULL
####################################################################################
#&targetscan('human','broadly');
#&targetscan('human','conserved');
#&targetscan('human','nonconserved');
#&targetscan('mouse','broadly');
#&targetscan('mouse','conserved');
#&targetscan('mouse','nonconserved');
sub targetscan(){

    my ($species,$family_conservation) = @_;
     my $targetscan_hash_file=$data_dir.'targetscan/'.$species.'_'.$family_conservation.'_targetScan.hash';
    if ($family_conservation eq "broadly"){$family_conservation=2;}elsif($family_conservation eq "conserved"){$family_conservation=1;}else{$family_conservation=0;}
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my $miRNA_Family_file = $data_dir."targetscan/".$tax_ID."_miR_Family_Info.txt";
    #targetScan/miR_Family_Info.txt
    my %mirFamily2mirAccession=();# hash mirFamilySeed_sequence->@miRAccessions in target Taxonomy
    my %mirFamily2mirFamilyString=();
    if (open(MIRFAMILY_FILE, $miRNA_Family_file)){
        @all = <MIRFAMILY_FILE>;
        close(MIRFAMILY_FILE);
        $line = shift @all; # header info deleted
    }# if (open(MIRFAMILY_FILE, $miRNA_Family_file)){
    else{
    print"Can't open file $miRNA_Family_file\n";
    }#else{
    while(@all){
        $line = shift @all;
        chomp($line);#remove new line
        if ($line){
            #miR family	Seed+m8	Species ID	MiRBase ID	Mature sequence	Family Conservation?	MiRBase Accession
            #let-7/98/4458/4500	GAGGUAG	9606	hsa-let-7a	UGAGGUAGUAGGUUGUAUAGUU 2	MIMAT0000062
            my ($miR_family_string, $seed_Seq, $species_ID, $mirBaseID, $matureSeq, $conservation_Nr,$mirBaseAccession) = split("\t", $line);
            if($species_ID eq $tax_ID and $conservation_Nr>=$family_conservation){
                push @{$mirFamily2mirAccession{$seed_Seq}},$mirBaseAccession;#array containing all the mirFamily members MirBase accession(species limited)
                if(not exists $mirFamily2mirFamilyString{$seed_Seq}){
                    $mirFamily2mirFamilyString{$seed_Seq}=$miR_family_string;          
                    }#if(not exists $mirFamily2mirFamilyString{$seed_Seq}){#
                }#if($tax_ID eq $species_ID){#if it is a target species miRNA
            }#if($line){
        }#while(@all){

    my $targetScan_summary_file = $data_dir."targetscan/mouse_targetscan_mirnaFamily2RefSeq.txt";
    if($tax_ID eq "9606"){$targetScan_summary_file =$data_dir."targetscan/human_targetscan_mirnaFamily2RefSeq.txt";}#if($tax_ID eq "9606"){#
    my %targetscan_hash=();
    if(open(TARGETSCAN_SUMMARY,$targetScan_summary_file)){
        @all = <TARGETSCAN_SUMMARY>;
        close(TARGETSCAN_SUMMARY);
        $line =shift @all;#delete head line
    }else{print "Can't open file $targetScan_summary_file\n";}
    while(@all){
        $line = shift(@all);
        chomp($line);
        if($line){
            #$Transcript_ID	GeneSymbol	$miRNAfamilySeed $Species_ID	$Total_num_conserved_sites	Number_of_conserved_8mer_sites	Number_of_conserved_7mer-m8_sites	Number_of_conserved_7mer-1a_sites	Total_num_nonconserved_sites	Number_of_nonconserved_8mer_sites	Number_of_nonconserved_7mer-m8_sites	Number_of_nonconserved 7mer-1a sites	Representative_miRNA	$Total_context_score	$Aggregate_PCT
            #NM_001001130	Zfp85-rs1	AAAGUGC	10090	0	0	0	0	1	0	0	1	mmu-miR-106a	-0.041	0.000
            my($refseq_Transcript_ID,$gene_symbol,$miR_family,$species_ID,$tncs,$nc8s,$nc78s,$nc7s,$tnncs,$nnc8s,$nnc78s,$nnc7s,$representative_miR,$total_context_score,$aggregatePCT)= split("\t",$line);
            #next if ($total_context_score eq "NULL" or $aggregatePCT eq "NULL" );
            #if($tncs>=1 and ($total_context_score +0.3)<=0 and $aggregatePCT>=0.75){#at least one consevative site match, tcns<-0.3,aPCT>=0.75
            if($tncs>=1 ){#at least one consevative site match,else delete
                if(exists $mirFamily2mirAccession{$miR_family}){
                    foreach my $currentmiRNA (@{$mirFamily2mirAccession{$miR_family}}){
                        $targetscan_hash{$currentmiRNA}{$refseq_Transcript_ID}{"total_Context_Score"}=$total_context_score;
                        $targetscan_hash{$currentmiRNA}{$refseq_Transcript_ID}{"aggregatePCT"}=$aggregatePCT;
                        #==TO__ADD a tracker or lable that there is result in targetScan for this miRNA and refSeq_Transcript_ID & change cutoff to variablesa
                    }#foreach my $currentmiRNA (@{$mirFamily2mirAccession{$miR_family}}){#take in every miRNA accession which have high match and exist in target organism
                }#if(exists $mirFamily2mirAccession{$miR_family}){#
            }#if($tncs>=1 and ($total_context_score +0.3)<=0 and $aggregatePCT>=0.75){#
            #else{print "$refseq_Transcript_ID\t$representative_miR\t$total_context_score\t$aggregatePCT\n";}#else{# TO__DELETE       
        }#if($line){#    
    }#while(@all){#
    ###FOR_TESTING###
    #print Dumper \%targetscan_hash;
    #print "------------------------------------------\n";
    #foreach my $mirna (sort keys %targetscan_hash){
    # foreach my $refSeq_T_ID (keys %{$targetscan_hash{$mirna}}){
    #        print "$mirna,$refSeq_T_ID: $targetscan_hash{$mirna}{$refSeq_T_ID}{'total_Context_Score'}\n";
    #    }#foreach my $refSeq_T_ID (keys %{$targetscan_hash{$mirna}}){
    #}#foreach my $mirna (sort keys %targetscan_hash){
    ###END_FOR_TESTING###
   
    store \%targetscan_hash,$targetscan_hash_file;
    return ($targetscan_hash_file);
}#sub targetScan
#END_SECTION_TARGETSCAN############################################################################

#####==BEGIN_PICTAR##################################################################################
#PICTAR offer the all the target site between each pair of miRNA(miRNA_name) and target(GeneSymbol),and the total score was calculated by each pair
#&pictar('mouse','most_conserved');
#&pictar('mouse','most_conserved');
#&pictar('human','mammal_conserved');
#&pictar('human','mammal_conserved');

sub pictar(){
    my ($species,$conservation) = @_;#$species='human','mouse';$conservation='most_conserved','mammal_conserved'
    #$species ||= 'mouse';
    #$conservation ||= 'most_conserved';$defined the inputfile
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my $pictar_target_file = $data_dir.'pictar/'.$tax_ID."_pictar_".$conservation.'.csv';
    my $missing_GeneName_file=  $data_dir.'pictar/'.$tax_ID."_pictar_".$conservation.'_missingGeneNames.txt';
    my %gene_symbol2uni_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_rnaGeneName2refseqs.hash')};
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')}; 
    print "try read in PICTAR file : $pictar_target_file\n";
    if(open(PICTAR_TARGET,$pictar_target_file)){
        @all = <PICTAR_TARGET>;
        close(PICTAR_TARGET);
        shift @all for 1 .. 4;#delete headers 4 lines in all
    }else{print "Can't open PICTAR file $pictar_target_file\n";}
    my %pictar_hash=();#$pictar_hash{mirna_accession}{target_uni_id}[{'strand+target_sites'}=@array{'pictar_score'}]
    my %pictar_tracking_hash=();# this is the hash to keep track of the mirna-genesymbol pair, to check if it already exists $pictar_tracking_hash{miRNA_name}{target_gene_symbol}
    while(@all){
        $line = shift(@all);
        chomp($line);
        if($line){
            #track_name,gene_symbol/NCBI RefSeq,gene_symbol/NCBI_RefSeq_location,$data_source,$score,target_site,target_site_location,genomic_strand,top-percent_value/////csv
            #hsa-miR-9-5p,ONECUT2,chr18:55102916-55158530,PICTAR,35,hsa-miR-9-5p,chr18:55153447-55153453,+,N/A
            #multiple lines to one mir-target pair, with different bidning site
            my($miRNA_name,$gene_symbol,$gene_location,$data_source,$pictar_score,$target_site,$target_site_location,$genomic_strand,$top_percent_value)= split(",",$line);
            my $mirna_accession;
            if(exists $mirAlias2Accession_hash{$miRNA_name}){
                $mirna_accession = $mirAlias2Accession_hash{$miRNA_name};
            }# if(exists $mirAlias2Accession_hash{$miRNA_name}){
            else{
                print "WRONG miRNA name: $miRNA_name\n";
            }#else           
            my @uni_ids=();            
            if(exists $gene_symbol2uni_hash{$gene_symbol}){
                @uni_ids=@{$gene_symbol2uni_hash{$gene_symbol}};
            }else{
                print "Can't find corresponding uni_ids for gene symbol $gene_symbol\n";
                open (MISSING_GENE,'>>'.$missing_GeneName_file)||die "Could not open file $!";
                print MISSING_GENE "$gene_symbol\n";
                close(MISSING_GENE);
            }#check if there is coresponding uni_id to mentioned gene symbol
            if(exists $pictar_tracking_hash{$miRNA_name}{$gene_symbol}){
                foreach my $uni_id(@uni_ids){
                    push @{$pictar_hash{$mirna_accession}{$uni_id}{'target_sites'}},$genomic_strand.$target_site_location;
                }#foreach my $uni_id(@uni_ids){
            }#if(exists $pictar_tracking_hash{$miRNA_name}{$gene_symbol})
            else{
                $pictar_tracking_hash{$miRNA_name}{$gene_symbol}=$pictar_score;
                foreach my $uni_id(@uni_ids){
                    push @{$pictar_hash{$mirna_accession}{$uni_id}{'target_sites'}},$genomic_strand.$target_site_location;
                    $pictar_hash{$mirna_accession}{$uni_id}{'pictar_score'}=$pictar_score;
                }#foreach my $uni_id(@uni_ids){                
            }#else this is a new pair
        }#if($line){
    }#while(@all){
    my $pictar_hash_file=$data_dir.'pictar/'.$tax_ID.'_pictar_'.$conservation.'.hash';
    store \%pictar_hash,$pictar_hash_file;
    return ($pictar_hash_file);
}#sub pictar
#END_SECTION_PICTAR##################################################################################

#==BEGIN_PITA=####################################################################################
#&pita('human','no_flank','TOP');
#&pita('human','flank','TOP');
#&pita('human','no_flank','ALL');
#&pita('human','flank','ALL');
#&pita('mouse','no_flank','TOP');
#&pita('mouse','flank','TOP');
#&pita('mouse','no_flank','ALL');
#&pita('mouse','flank','ALL');
sub pita(){
    my ($species,$flank,$stringency) = @_;
    #$species='human','mouse';$flank='no_flank','flank';$stringency='TOP','ALL'
    #flank will allow certain shift of matching: 3upstream, 15downstream
    #already replaced the soft links
    ##stringency:TOP=>the top predictions (having a full match 7- or 8-mer seed and a conservation score of 0.9 or higher)All=>the complete list of predictions 
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my $pita_target_file = $data_dir.'pita/'.$tax_ID."_pita_".$flank.'_'.$stringency.'.tab';
    my $missing_GeneName_file=  $data_dir.'pita/'.$tax_ID."_pita_".$flank.'_'.$stringency.'.missingGeneNames.txt';
    my %gene_symbol2uni_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_rnaGeneName2refseqs.hash')};
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')}; 
    print "try read in PITA target file : $pita_target_file\n";
    if(open(PITA_TARGET,$pita_target_file)){
        @all = <PITA_TARGET>;
        close(PITA_TARGET);
        shift @all;#delete header 1 line
    }else{print "Can't open PITA_TARGET file $pita_target_file\n";}
    my %pita_hash=();#$pita_hash{mirna_accession}{target_uni_id}[{'strand+target_sites'}=@,${'pictar_score'}]
    my %pita_tracking_hash=();# this is the hash to keep track of the mirna-genesymbol pair, to check if it already exists $pictar_tracking_hash{miRNA_name}{target_gene_symbol}
    while(@all){
        $line = shift(@all);
        chomp($line);
        if($line){
            #RefSeq	Names	microRNA	Nr_Sites	Score
            #NM_006599;NM_138713;NM_138714;NM_173214	NFAT5	hsa-miR-1207-5p	6	-37.69
            #NM_001012426;NM_001012427;NM_138457	FOXP4	hsa-miR-939	19	-37.46
            my($refseqs_string,$gene_symbol,$miRNA_name,$nr_sites,$pita_score)= split("\t",$line);
            my @pita_refseqs=split(";",$refseqs_string);
            my $mirna_accession;
            if(exists $mirAlias2Accession_hash{$miRNA_name}){
                $mirna_accession = $mirAlias2Accession_hash{$miRNA_name};
            }# if(exists $mirAlias2Accession_hash{$miRNA_name}){
            else{
                print "WRONG miRNA name: $miRNA_name\n";
            }#else           
            my @uni_ids=();            
            if(exists $gene_symbol2uni_hash{$gene_symbol}){#check if the genesymbol correspondent NMs 
                @uni_ids=@{$gene_symbol2uni_hash{$gene_symbol}};
                #push @pita_refseqs,@uni_ids;
                @pita_refseqs=uniq(@pita_refseqs);                
            }else{
                print "Can't find corresponding PITA uni_ids for gene symbol $gene_symbol\n";
                open (MISSING_GENE,'>>'.$missing_GeneName_file)||die "Could not open file $!";
                print MISSING_GENE "$gene_symbol\n";
                close(MISSING_GENE);
            }#else#check if there is coresponding uni_id to mentioned gene symbol
            foreach(@pita_refseqs){
                $pita_hash{$mirna_accession}{$_}{"nr_pita_binding_sites_per_gene"}=$nr_sites;
                $pita_hash{$mirna_accession}{$_}{"pita_score_per_gene"}=$pita_score;
            }# foreach(@pita_refseqs){
        }#if($line){
    }#while(@all){
    my $pita_hash_file= $data_dir.'pita/'.$tax_ID."_pita_".$flank.'_'.$stringency.'.hash';
    store \%pita_hash,$pita_hash_file;
    return ($pita_hash_file);
}#sub pictar
##END_SECTION_PITA###################################################################################

sub nonRedundantList{
    my $val = shift;
    my @val=@$val;                                                
    my %cpt;
    foreach (@val){                                               
        $cpt{$_}++;                                                 
    }                                                             
    return keys %cpt;
} #sub redundantList{

####==BEGIN_MICRORNA.ORG===#####################################
    #microRNA.org : http://www.microrna.org/microrna/getDownloads.do
    #it used a modified miRanda algorithm used an alignment score threshold of 120 instead of the previous threshold of 140
    #mirSVR is a regression model that computes a weighted sum of a number of sequence and context features f the predicted miRNA::mRNA duplex.mirSVR is a regression model that computes a  weighted sum of a number of sequence and context features of the predicted miRNA::mRNA duplex.
#microrna_org('human','conserved');
#microrna_org('human','non_conserved');

#microrna_org('mouse','conserved');
#microrna_org('mouse','non_conserved');
sub microrna_org(){
    my ($species,$conservation) = @_;
    #$species='human','mouse',$conservation='conserved','non_conserved'
    ##conservation indicated if the miRNA under intestigation is conserved across species
    #names of file need redefined    
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my $microrna_org_file = $data_dir.'microrna_org/'.$tax_ID."_microrna_org_".$conservation.'.txt';
    print "try read in  mocrorna_org file : $microrna_org_file\n";
    if(open(MICRO_ORG,$microrna_org_file)){
        @all = <MICRO_ORG>;
        close(MICRO_ORG);
        shift @all;#delete header 1 line
    }else{print "Can't open micro_org file $microrna_org_file\n";}
    my %microrna_org_hash = ();#$microrna_hash{mirna_accession}{target_uni_id}   
    while(@all){
        $line = shift(@all);
        chomp($line);
        if($line){
            #mirbase_acc	mirna_name	gene_id	gene_symbol	transcript_id	ext_transcript_id	mirna_alignment	alignment	gene_alignment	mirna_start	mirna_end	gene_start	gene_end	genome_coordinates	conservation	align_score	seed_cat	energy	mirsvr_score
            #MIMAT0000062	hsa-let-7a	5270	SERPINE2	uc002vnu.2	NM_006216	uuGAUAUGUUGGAUGAU-GGAGu	  | :|: ||:|| ||| |||| 	aaCGGUGAAAUCU-CUAGCCUCu	2	21	495	516	[hg19:2:224840068-224840089:-]	0.5684	122	0	-14.73	-0.7269
            my($mirbase_acc,$mirna_name,$gene_id,$gene_symbol,$transcript_id,$ext_transcript_id,$mirna_alignment,$alignment,$gene_alignment,$mirna_start,$mirna_end,$gene_start,$gene_end,$genome_coordinates,$conservation,$align_score,$seed_cat,$energy,$mirsvr_score)= split("\t",$line);
            push @{$microrna_org_hash{$mirbase_acc}{$ext_transcript_id}{"miRanda_align_score"}},$align_score;
            push @{$microrna_org_hash{$mirbase_acc}{$ext_transcript_id}{"misvr_score"}},$mirsvr_score;
        }#if($line){
    }#while(@all){
    my $microrna_org_hash_file= $data_dir.'microrna_org/'.$tax_ID."_microrna_org_".$conservation.'.hash';
    store \%microrna_org_hash,$microrna_org_hash_file;
    return ($microrna_org_hash_file);
}#sub microrna_org

###==END_SECTION_MICRORNA.ORG==###############################

###BEGIN_SECTION_MIRDB=====########################################
#&mirdb('human');
#&mirdb(mouse);
sub mirdb(){
    my ($species) = @_;
    #$species='human','mouse'
    my $mirdb_file = $data_dir.'mirdb/'.$species.'_mirdb.txt';
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')}; 
    print "try read in  mirdb file : $mirdb_file\n";
    if(open(MIRDB,$mirdb_file)){
        @all = <MIRDB>;
        close(MIRDB);       
    }else{print "Can't open miRDB file $mirdb_file\n";}
    my %miRDB_hash = ();#$miRDB_hash{mirna_accession}{target_uni_id}   
    while(@all){
        $line = shift(@all);
        chomp($line);
        if($line){
            #hsa-let-7a-2-3p	XM_003403774	70.5358226315502
            #hsa-let-7a-2-3p	NM_001099733	75.4047713654
            ##this file have a one(miR-refseq pair)-one score relationship
            my($miRNA_name,$refseq_id,$miRDB_score)= split("\t",$line);
            my $mirna_accession;
            if(exists $mirAlias2Accession_hash{$miRNA_name}){
                $mirna_accession = $mirAlias2Accession_hash{$miRNA_name};
                $miRDB_hash{$mirna_accession}{$refseq_id}=$miRDB_score;
            }# if(exists $mirAlias2Accession_hash{$miRNA_name}){
            else{
                print "WRONG miRNA name: $miRNA_name\n";
            }#else           
           
        }#if($line){
    }#while(@all){
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my $miRDB_hash_file= $data_dir.'mirdb/'.$tax_ID.'_mirdb.hash';
    store \%miRDB_hash,$miRDB_hash_file;
    return ($miRDB_hash_file);
}#sub mirdb
###END_SECTION_MIRDB===############################################
###MIRGEN#############################################################
#MIRGEN is now swiched to webserver..it is a database of microRNA genomic information and regulation.
###this database has been updated and the previous prediction was done with microT.
###==BEGIN_MIRNAMAP==########################################################
#miRNAMap 2.0 http://mirnamap.mbc.nctu.edu.tw/html/about.html
#miRNAMap collected the know varified and miRNAtargets(Tarbase and Surveying literature) in human,mouse rat, In addition targetScan,miranda rnahybird mirTar were used in prediction. the creteria: 1)predicted by two or moretools;2)target gene containd multiple sites, 3) target sites is accessible by Sfold
#mirnamap('human',1,1,1);
#mirnamap('human',1,2,1);
#mirnamap('human',2,1,1);
#mirnamap('human',2,2,1);
#mirnamap('mouse',1,1,1);
#mirnamap('mouse',1,2,1);
#mirnamap('mouse',2,1,1);
#mirnamap('mouse',2,2,1);
sub mirnamap(){
    my ($species,$creterion1,$creterion2,$creterion3) = @_;
    #$species='human','mouse'$creterion1:Nr, the minimum of reported tools;$creterion2_minimum sites nr between target and miR;$accessible or not (1 0)is this site
    $species//='human';
    $creterion1//=2;
    $creterion2//=1;
    $creterion3//=1;    
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my $mirnamap_target_file = $data_dir.'mirnamap/'.$tax_ID.'_target_mirnamap.txt';
    my $missing_ensembl_file=  $data_dir.'mirnamap/'.$tax_ID.'_mirnamap_'.$creterion1.$creterion2.$creterion3.'missingEnsemblID.txt';
    my %ensembl_tr2refseq_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_ensembl_transcriptID2refseqID.hash')};
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};     
    print "try read in miRNAmap target file : $mirnamap_target_file\n";
    if(open(MIRNAMAP_TARGET,$mirnamap_target_file)){
        @all = <MIRNAMAP_TARGET>;
        close(MIRNAMAP_TARGET);
        shift @all;#delete header 1 line
    }else{print "Can't open miRNAmap TARGET file $mirnamap_target_file\n";}
    #only the matches passed the three creterion will be added to the corresponding hash
    my %mirnamap_hash=();#$pita_hash{mirna_accession}{target_uni_id}{location}/{reported tools}{}
    while(@all){
        $line = shift(@all);
        chomp($line);
        if($line){
            #mature_miRNA	Ensembl_transcript_ID	target_start	tsrget_end	miRNA_3-5	alignment	target_5-3	tool_name	criterion_1	criterion_2	criterion_3
            #hsa-let-7a	ENST00000373266	55	62	U----------------UGAUAUGUUGGAU-----GAUGGAGU	                 |||:  |::||       |||||||	AGGGCAGGGCCCAGCA ACUGCCCGGCCAGUCCUCCUACCUCC	TargetScan	2	1	1
            my($mature_miR_name,$ensembl_transcript_ID,$target_start,$target_end,$miRNA_3_5,$alignment,$target_5_3,$tool_name,$criter_1,$criter_2,$criter_3)= split("\t",$line);
            my $mirna_accession;
            my $refseq_ID;
            if(exists $mirAlias2Accession_hash{$mature_miR_name}){
                $mirna_accession = $mirAlias2Accession_hash{$mature_miR_name};
            }# if(exists $mirAlias2Accession_hash{$miRNA_name}){
            else{
                print "WRONG miRNA name: $mature_miR_name\n";
                last;#exit th cycle, this line is abnormaline like 
            }#else 
            if(exists $ensembl_tr2refseq_hash{$ensembl_transcript_ID}){
                #######################====BUG:there are many old ensembl transcript ID need to be transformed to check if they retired or have successor by ensembl API.##########################################
                $refseq_ID = $ensembl_tr2refseq_hash{$ensembl_transcript_ID};
                if($criter_1>=$creterion1 and $criter_2>=$creterion2 and $criter_3>=$creterion3){
                    push @{$mirnamap_hash{$mirna_accession}{$refseq_ID}{'site_location'}},$target_start.'_'.$target_end;
                    push @{$mirnamap_hash{$mirna_accession}{$refseq_ID}{'prediction_tool'}},$tool_name;
                    push @{$mirnamap_hash{$mirna_accession}{$refseq_ID}{'common_tool_nr'}},$criter_1;
                }# if($criter_1>=$creterion1 and $criter_2>){
            } else{
                    print "$line\n";
                    open (MISSING_GENE,'>>'.$missing_ensembl_file)||die "Could not open file $!";
                    print MISSING_GENE "$ensembl_transcript_ID\n";
                    close(MISSING_GENE);
            }#else#check if there is coresponding uni_id to mentioned gene ensembl transcript id
            
        }#if($line){
    }#while(@all){
    my $mirnamap_hash_file= $data_dir.'mirnamap/'.$tax_ID."_mirnamap_".$creterion1.'_'.$creterion2.'_'.$creterion3.'.hash';
    store \%mirnamap_hash,$mirnamap_hash_file;
    return ($mirnamap_hash_file);
}#sub mirnamap
#####==END_SECTION_MIRNAMAP==#####################################################

######==BEGIN_REPTAR##################################################################
#&reptar('human');
#&reptar('mouse');
sub reptar(){
    my ($species,$mfe_upper,$nmfeb_lower,$rf_lower) = @_;#define the species(the profile to readin), the minimal free energy of inding , #the normalized minimal free energy of binding, and the repeating motifs minimum
    $mfe_upper//= -(15);
    $nmfeb_lower//=0.1;
    $rf_lower//=2;#BUG:currently still don't know where to use this parameter for cut
    #$species='human','mouse',set default
    my $reptar_file = $data_dir.'reptar/'.$species.'_reptar.txt';
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my %refseq2_rnaGeneNames_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_refseq2rnaGeneNames.hash')};
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    print "try read in  reptar target file : $reptar_file\n";
    if(open(REPTAR,$reptar_file)){
        @all = <REPTAR>;#this file has no headline
        close(REPTAR);       
    }else{print "Can't open reptar file $reptar_file\n";}
    my %reptar_hash = ();#$miRDB_hash{mirna_accession}{target_uni_id} 
   
    while(@all){
        $line = shift(@all);
        chomp($line);
        if($line){
            #gene symbol #gnene accesssion #microRNA symbol # binding ste start position # binding site end position#minimal free energy<=-10, default to be less than -15 #normallized miminul free energydefault >=0.1 #g-U based pairs #binding sites pattern#3 ustr conservation score[0-1]the higher the better, default 0 #know repeting elements #repeated motifs default >2# algorithm used
            #ABCC2:::NM_000392	hsa-miR-193bstar	beg:171	end:187	mfe:-16.50nm:0.374	gu:2/12	prof:10,m8	pic:3' AGTAGAGCGGGAGTTTTGGGGC 5'&||.|     .|||||||    &5' ----CTTGA----GAAACCCC- 3'	b_cons:0.00	t_cons:0.02	rep:LINE/L1-L1ME3B,	Rep
            #ABCC2:::NM_000392	hsa-miR-609	beg:112	end:131	mfe:-16.90	nm:0.476 gu:2/14	prof:10,m8	pic:3' TCTCTACTCTCTTTGTGGGA 5'&     .||| ||  .|| ||||    &5' --GGATAAGT-GAACACCC- 3'	b_cons:0.02	t_cons:0.02	rep:LINE/L1-L1ME3B,	Rep
            my ($gene_name_refseq_ids, $mirna_id, $start, $end, $mfe_sc, $nm_sc, $gu, $prof, $pic, $b_cons, $t_cons, $rep, $alg) = split("\t", $line);
            #print "$mirna_id\n";
            $mirna_id=~s/star/\*/;
            #print "$mirna_id\n";
            my $mirna_accession;     
            my ($gene_symbol,$refseq_id)=split(":::",$gene_name_refseq_ids);
            my ($label,$mfe)=split(":",$mfe_sc);
            my ($label,$nm)=split(":",$nm_sc);
            if(exists $mirAlias2Accession_hash{$mirna_id} and exists  $refseq2_rnaGeneNames_hash{$refseq_id}){
                $mirna_accession = $mirAlias2Accession_hash{$mirna_id};
                if($mfe<$mfe_upper and $nm>=$nmfeb_lower ){#if this line passed the quality control
                    push @{$reptar_hash{$mirna_accession}{$refseq_id}},$start.$end.$prof;
                } #if($mfe<$mfe_upper and $nm>=$nmfeb_lower ){
            }# if(exists $mirAlias2Accession_hash{$miRNA_name}){
            else{
                print "WRONG miRNA name:$mirna_id\n" if(! exists $mirAlias2Accession_hash{$mirna_id});
                print "WRONG refseq ID:$refseq_id\n" if (! exists $refseq2_rnaGeneNames_hash{$refseq_id});
            }#else
        }#if($line){
    }#while(@all){
    my $reptar_hash_file= $data_dir.'reptar/'.$tax_ID.'_nfe_'.$mfe_upper.'_nm_'.$nmfeb_lower.'_reptar.hash';
    store \%reptar_hash,$reptar_hash_file;
    return ($reptar_hash_file);
}#sub reptar
######==END_SECTION_REPTAR#############################################################

######==BEGIN_EIMMO===##################################################################
#eimmo can set the ct off of P>0.5 more dedium confidence and P>0.8 for high confidence miRNA target sites besides, the cons can be sent if this site is perfectly conserved in the calde
#the output hash key will be mirna_accession, refseq_id, the value is a Pvalue
#&eimmo('human',undef,undef,'');# my ($species,$P_lower,$conserve_min,$cds)= @_;#define the species(the profile to readin), the minimal free energy of inding , #the normalized minimal free energy of binding, and the repeating motifs minimum#
#&eimmo('mouse',undef,undef,'');
#&eimmo('human',undef,undef,'CDS');
#&eimmo('mouse',undef,undef,'CDS');
sub eimmo(){
    my ($species,$P_lower,$conserve_min,$cds)= @_;#define the species(the profile to readin), the minimal free energy of inding , #the normalized minimal free energy of binding, and the repeating motifs minimum
    $species//= 'mouse';
    print "species: $species\n";
    $P_lower//=0.5;#by default midium confidence, can be increased to 0.8 for high confidence
    print "P_lower: $P_lower\n";
    $conserve_min//=0;#by default there is no conservation cut ,if there is request then increase conservation to 1
    print "conserve_min: $conserve_min\n";
    $cds//='';# if wanna see the matchs in CDS region, the input need to be 'CDS'
    print "cds: $cds\n";
    
    #$species='human','mouse',set default
    my $eimmo_file = $data_dir.'eimmo/'.$species.$cds.'_eimmo.txt';
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my %refseq2_rnaGeneNames_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_refseq2rnaGeneNames.hash')};
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    print "try read in  eimmo target file : $eimmo_file\n";
    if(open(EIMMO,$eimmo_file)){
        @all = <EIMMO>;#this file has no headline
        shift(@all);
        close(EIMMO);       
    }else{print "Can't open eimmo file $eimmo_file\n";}
    my %eimmo_hash = ();#$miRDB_hash{mirna_accession}{target_uni_id} 
   
    while(@all){
        $line = shift(@all);
        chomp($line);
        if($line){
            #chrom	strand	chrFrom	chrTo	3UtrFrom	3UtrTo	Annot	Cons	p
            #chr1	+	76228456	76228463	8	15	NM_000016:hsa-miR-139-5p	00000110	0.136
            #chr1	+	76228514	76228521	66	73	NM_000016:hsa-miR-129*	10011110	0.575  
            my ($chrom,$strand,$chrFrom,$chrTo,$Utr3From,$Utr3To,$refseq_id_mirna_id,$Cons,$p) = split("\t", $line);               
            my ($refseq_id,$mirna_id)=split(":",$refseq_id_mirna_id);
            if(exists $mirAlias2Accession_hash{$mirna_id}
                # and exists  $refseq2_rnaGeneNames_hash{$refseq_id}
            ){
                my $mirna_accession = $mirAlias2Accession_hash{$mirna_id};
                if($p>$P_lower){#if this line passed the quality control
                    
                    #my @current;
                    #  push @current,$Utr3From."\t".$Utr3To."\t".$p."\t".$Cons;
                    if (exists $eimmo_hash{$mirna_accession}{$refseq_id} and $eimmo_hash{$mirna_accession}{$refseq_id}>$p){
                        #do nothing 
                    }else{
                        #$eimmo_hash{$mirna_accession}{$refseq_id}=\@current;
                        $eimmo_hash{$mirna_accession}{$refseq_id}=$p;}
                    print "mirna:ACC:$mirna_accession\t";
                    print "mirna:ACC:$refseq_id\n";

                } #if($mfe<$mfe_upper and $nm>=$nmfeb_lower ){
            }# if(exists $mirAlias2Accession_hash{$miRNA_name}){
            else{
                print "WRONG miRNA name:$mirna_id\n" if(! exists $mirAlias2Accession_hash{$mirna_id});
                print "WRONG refseq ID:$refseq_id\n" if (! exists $refseq2_rnaGeneNames_hash{$refseq_id});
            }#else
        }#if($line){
    }#while(@all){
    my $eimmo_hash_file= $data_dir.'eimmo/'.$tax_ID.'_P_'.$P_lower.$cds.'_eimmo.hash';
    # print Dumper \%eimmo_hash;
    
    store \%eimmo_hash,$eimmo_hash_file;
    return ($eimmo_hash_file);
}#sub eimmo
######==END_SECTION_EIMMO#############################################################

######==BEGIN_microT4===##################################################################
#http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=site/help&topic=microtv4
#the output hash key will be mirna_accession, refseq_id, the value is a miTG score it is multiple line to one of the ensembl gene id and miRNA pair, documented each biding location.by default the cutoff filter for CDS is 0.7
#perl -pne 'if(not $.%2){s/\n/,/;}' microtv4_data_trail.csv -p can be canceled 
#the files need to be treated first since it splitted one lines into several we only need the ones that are headlines
#grep ",mmu-" microT_CDS_brief.csv >mouse_microT_CDS_brief.csv
# grep ",hsa-" microtv4_brief.csv  >human_microtv4_brief.csv
#&microT4();#mouse, microt4,0.7
#&microT4('human',);#human,micro4,0.7
#&microT4('mouse','cds');#mouse,cds microt, 0.7
#&microT4('human','cds');

sub microT4(){
    my ($species,$cds,$cutoff)= @_;#define the species
    $species//= 'mouse';
    $cds//='microtv4';#by default take the microv4 file/$cds=CDS then take the microCDS file
    $cutoff//=0.7;#by default the cutoff for miRTg score is 0.7,could be adjusted to higher
    my $microT_file = $data_dir."microt/".$species.'_'.$cds."_microT.tsv";
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}

    my $missing_GeneName_file=  $data_dir.'microt/microT_'.$tax_ID.$cds.$cutoff.'_missingENST.txt';
    my %gene_symbol2uni_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_rnaGeneName2refseqs.hash')};
    my %ensembl_tr2refseq_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_ensembl_transcriptID2refseqID.hash')};
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    print "try read in  microT target file : $microT_file\n";
    if(open(MICROT,$microT_file)){
        @all = <MICROT>;#this file has no headline
        #shift(@all);#there is no header
        close(MICROT);       
    }else{print "Can't open microT file $microT_file\n";}
    my %microT_hash = ();#$microT_hash{mirna_accession}{target_uni_id} 
    while(@all){
        $line = shift(@all);
        chomp($line);
        if($line){     
            #TranscriptId,GeneId(name),Mirna-Name(miRBase-version),miTG-score
            #F52F10.2,F52F10.2(F52F10.2),cel-miR-62(18),0.488            #
            my ($transcriptid,$geneid_name,$mirnaName_version,$miTGscore) = split(",", $line);               
            my ($mirna_id,$rest)=split('\\(', $mirnaName_version);
            my ($gene_ID,$gene_name)=split(/[()]+/,$geneid_name);#need checking  
            my $mirna_accession;
            my $refseq_id;
            if(exists $mirAlias2Accession_hash{$mirna_id}){
                $mirna_accession = $mirAlias2Accession_hash{$mirna_id};}else{
                print "WRONG miRNA name:$mirna_id\n";               
            }#else
            if(exists $ensembl_tr2refseq_hash{$transcriptid}){
                $refseq_id= $ensembl_tr2refseq_hash{$transcriptid};#this line still need review to ensure that every ensembl_tr can find a refseq id##==BUG
            }else{
                print "$transcriptid not exist!\n";
                open (MISSING_GENE,'>>'.$missing_GeneName_file)||die "Could not open file $!";
                print MISSING_GENE "$transcriptid:$gene_ID:$gene_name\n";
                close(MISSING_GENE);
            }
            if($miTGscore>=$cutoff and $refseq_id and $mirna_accession){#if this line passed the quality control
                $microT_hash{$mirna_accession}{$refseq_id}=$miTGscore;
            } #if($miTGscore>=$cutoff ){
           
            
        }#if($line){
    }#while(@all){
    my $microT_hash_file= $data_dir.'microt/'.$tax_ID.'_TG'.$cutoff.'_'.$cds.'_microT.hash';
    store \%microT_hash,$microT_hash_file;
    return ($microT_hash_file);
}#sub microt4
######==END_SECTION_microT4#############################################################




######==BEGIN_mirTARBASE===##################################################################
#http://mirtarbase.mbc.nctu.edu.tw/php/download.php
# the experimentally validated microRNA-target interactions database
#&mirtarbase();
#&mirtarbase('human');
#&mirtarbase(,'SE');
#&mirtarbase('human','SE');
sub mirtarbase(){
    my ($species,$evidence)= @_;#define the species(human, mouse) and the support by evidence:(SE or WK) SE strong evidence(Reporter assay and Western blot),WK including every type of evidence
    $species//= 'mouse';
    $evidence//='WK';# by default,take in any evidence     if only want higher level, then set to SE
    my $mirtarbase_file = $data_dir.'mirtarbase/'.$species.'_MTI.tsv';
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my %gene_symbol2uni_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_rnaGeneName2refseqs.hash')};
    my %entrez2uni_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_enterz2refseqs.hash')};
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    print "try read in  miRTarBase target file : $mirtarbase_file\n";
    if(open(MIRTARBASE,$mirtarbase_file)){
        @all = <MIRTARBASE>;#this file has no headline
        shift(@all);#delete head line
        close(MIRTARBASE);       
    }else{print "Can't open miRTarBase file $mirtarbase_file\n";}
    my %mirtarbase_hash = ();#$microT4_hash{mirna_accession}{target_uni_id} 
    while(@all){
        $line = shift(@all);
        chomp($line);        
        if($line){
            #miRTarBase ID	miRNA	Species (miRNA)	Target Gene	Target Gene (Entrez ID)Species (Target Gene)	Experiments	Support Type	References (PMID)
            #MIRT015952	mmu-miR-124-3p	Mus musculus	Aatk	11302	Mus musculus	Microarray	Functional MTI (Weak)	18619591
            #MIRT015602	mmu-miR-9-5p	Mus musculus	Aatk	11302	Mus musculus	Sequencing	Functional MTI 	19536157
            #MIRT015489	mmu-miR-9-5p	Mus musculus	Abca1	11303	Mus musculus	Sequencing	Functional MTI (Weak)	19536157           
            my ($miRTarBase_ID,$miRNA,$species_string,$target_Gene,$target_Entrez,$species,$experiment,$support_Type,$references)= split("\t", $line);            
            my $mirna_accession;
            $target_Gene = uc($target_Gene);
            my @refseq_ids;
            if(exists $mirAlias2Accession_hash{$miRNA}){
                $mirna_accession = $mirAlias2Accession_hash{$miRNA};
            }else{ print "WRONG miRNA name:$miRNA\n";next;}            
            if($support_Type =~ /Weak/ and $evidence eq 'SE'){#if this line dont fit the required evidence strengh
                 next;
                } # if($experiment_Support_Type=~/Weak/ and $evidence eq 'SE'){#if this line dont fit the required evidence strengh
            if((!exists $gene_symbol2uni_hash{$target_Gene}) and (! exists $entrez2uni_hash{$target_Entrez})){
                 print "can't find refseqs : $target_Gene, $target_Entrez\n";
                 next;
           }else{
               push @refseq_ids,@{$gene_symbol2uni_hash{$target_Gene}} if (exists $gene_symbol2uni_hash{$target_Gene});
               push @refseq_ids, @{$entrez2uni_hash{$target_Entrez}} if (exists $entrez2uni_hash{$target_Entrez});
               uniq(@refseq_ids);
               foreach my $refseq_id(@refseq_ids){
                   $mirtarbase_hash{$mirna_accession}{$refseq_id} = $target_Gene.':'.$experiment.':'.$references;
               }            
            }#else        
        }#if($line){
    }#while(@all){
    my $mirtarbase_hash_file= $data_dir.'mirtarbase/'.$tax_ID.'_Evidence_'.$evidence.'_mirtarbase.hash';
    store \%mirtarbase_hash,$mirtarbase_hash_file;
    return ($mirtarbase_hash_file);
}#sub mirtarbase

######==END_SECTION_mirTARBASE#############################################################

######==BEGIN_TARGETRANK===##################################################################
#this part still need improvement ,since the scoring system build in targetrank was not taken into consideration 
#&targetrank('human',1);
#&targetrank('human',0);
#&targetrank('mouse',1);
#&targetrank('mouse',0);
sub targetrank(){
    my ($species,$conservation)= @_;#define the species(human, mouse) and the level of conservation0(not)1 yes
    $species//= 'mouse';
    $conservation//=0;# by default,take in any evidence     
    my $targetrank_file = $data_dir.'targetrank/'.$species.'_mirbase_targetrank.txt';
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my %gene_symbol2uni_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_rnaGeneName2refseqs.hash')};    
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    print "try read in  TARGETRANK target file : $targetrank_file\n";
    if(open(TARGETRANK,$targetrank_file)){
        @all = <TARGETRANK>;#this file has no headline
        shift(@all);#delete head line
        close(TARGETRANK);       
    }else{print "Can't open targetrank file $targetrank_file\n";}
    my %targetrank_hash = ();#$microT4_hash{mirna_accession}{target_uni_id} 
    while(@all){
        $line = shift(@all);
        chomp($line);        
        if($line){
            #miRBase ID	Gene Symbol	Number of Refseq UTRs	UTR Number	Refseq ID	Seed Match Type	Seed Match Conservation#multiple line correspond to one NM and mirna pair, different match site
            #hsa-miR-323-3p	UBE2Q1	1	1	NM_017582	8mer	1
            #hsa-miR-323-3p	UBE2Q1	1	1	NM_017582	M8 7mer	0
            #hsa-miR-323-3p	UBE2Q1	1	1	NM_017582	A1 7mer	1  
            #For each match S(m) = SSeedMatchType(m)+R5'conservation(m)+R3'AU(m),current don't know how to calculate the score schema, but so temperally collect every match
            my ($miRBase_ID,$Gene_Symbol,$Number_of_Refseq_UTRs,$UTR_Number,$Refseq_ID,$Seed_Match_Type,$Seed_Match_Conservation)= split("\t+", $line);            
            my $mirna_accession;
            print "seed match conservation $Seed_Match_Conservation\n";
           
            my @refseq_ids;
            if(exists $mirAlias2Accession_hash{$miRBase_ID}){
                $mirna_accession = $mirAlias2Accession_hash{$miRBase_ID};
            }else{ print "WRONG miRNA name: $miRBase_ID\n";next;}            
            if($Seed_Match_Conservation < $conservation){#if this line dont fit the required evidence strengh
                 next;
                } # if($Seed_Match_Conservation < $conservation){#if this line dont fit the required conservation eg it is 0 while conservation cutoff is 
            else{
                push @{$targetrank_hash{$mirna_accession}{$Refseq_ID}}, $Seed_Match_Type.':'.$Seed_Match_Conservation.':'.$Number_of_Refseq_UTRs.':'.$UTR_Number;               
            }#else        
        }#if($line){
    }#while(@all){
    my $targetrank_hash_file= $data_dir.'targetrank/'.$tax_ID.'_Conservation_'.$conservation.'_targetrank.hash';
    store \%targetrank_hash,$targetrank_hash_file;
    return ($targetrank_hash_file);
}#sub targetrank

######==END_SECTION_TARGETRANK#############################################################

######==BEGIN_MIRTAR===##################################################################
#Putative interactions between miRNA and genes (default constraints MFE<-14, and score >140)
#&mirtar('3utr');
#&mirtar('5utr');
#&mirtar('all');
#&mirtar('cds');
sub mirtar(){
    my ($region)= @_;#currently there is only mirtar base for human!!!there is no data for other oganism
    my $species//= 'human';
    $region//='3utr';# by default,take the 3.utr, however it can also be (all/5utr/cds)     
    my $mirtar_file = $data_dir.'/mirtar/'.$species.'_'.$region.'_mirtar.txt';
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my %gene_symbol2uni_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_rnaGeneName2refseqs.hash')};    
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    my %entrez2uni_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_enterz2refseqs.hash')};
    print "try read in  mirtar target file : $mirtar_file\n";
    if(open(MIRTAR,$mirtar_file)){
        @all = <MIRTAR>;#this file has no headline
        close(MIRTAR);       
    }else{print "Can't open mirtar file $mirtar_file\n";}
    my %mirtar_hash = ();#$microT4_hash{mirna_accession}{target_uni_id} 
    while(@all){
        $line = shift(@all);
        chomp($line);        
        if($line){
            #mature mirnaname, entrezidlist, gene symbol list
            my ($miRBase_ID,$entrezid_list,$Gene_Symbol_list)= split("\t", $line);            
            my $mirna_accession;
            my @gene_symbols=split(';',$Gene_Symbol_list);
            my @entrezids=split(';',$entrezid_list);
            my @refseq_ids;
            if(exists $mirAlias2Accession_hash{$miRBase_ID}){
                $mirna_accession = $mirAlias2Accession_hash{$miRBase_ID};
            }else{
                print "WRONG miRNA name: $miRBase_ID\n";next;
            } 
            foreach my $gene_symbol (@gene_symbols){
                if (exists $gene_symbol2uni_hash{$gene_symbol}){push @refseq_ids, @{$gene_symbol2uni_hash{$gene_symbol}};}                          
            }
            foreach my $entrezid (@entrezids){
                if (exists $entrez2uni_hash{$entrezid}){push @refseq_ids, @{$entrez2uni_hash{$entrezid}};}                          
            }
            @refseq_ids=uniq(@refseq_ids);
            foreach my $refseq_ID(@refseq_ids){
            $mirtar_hash{$mirna_accession}{$refseq_ID}=1;
            }               
        }#if($line){
    }#while(@all){
    my $mirtar_hash_file= $data_dir.'mirtar/'.$tax_ID.'_'.$region.'_mirtar.hash';
    store \%mirtar_hash,$mirtar_hash_file;
    return ($mirtar_hash_file);
}#sub mirtar

######==END_SECTION_mirtar#############################################################

######==BEGIN_TARGETSPY===##################################################################
#http://www.targetspy.org/index.php?down=true
#score cut off-sensitive or more strict, which is specific
#have seed mode and non-seed mode, which require a perfect seed match
#&targetspy('human','sens','seed');
#&targetspy('human','sens','all');
#&targetspy('human','spec','seed');
#&targetspy('mouse','sens','seed');
#&targetspy('mouse','sens','all');
#&targetspy('mouse','spec','seed');

sub targetspy(){
    my ($species,$sensitivity,$seed)= @_;#define the species(human, mouse) 
    $species//= 'mouse';
    $sensitivity//='sens';# by default,senc,can also be altered to spec   
    $seed//='seed';#or all
    my $targetspy_file;
    if($seed eq 'seed'){
        $targetspy_file = $data_dir.'targetspy/'.$species.'_'.$sensitivity.'_'.$seed.'_targetspy.tsv';
    }else{
         $targetspy_file = $data_dir.'targetspy/'.$species.'_all_targetspy.tsv';
    }
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    print "try read in  targetspy target file : $targetspy_file\n";
    if(open(TARGETSPY,$targetspy_file)){
        @all = <TARGETSPY>;#this file has no headline
        close(TARGETSPY);       
    }else{print "Can't open targetspy file $targetspy_file\n";}
    my %targetspy_hash = ();#$microT4_hash{mirna_accession}{target_uni_id} 
    while(@all){
        $line = shift(@all);
        chomp($line);        
        if($line){
            #let-7a	NM_177398	897	921	0	22	ggctaatagtatgtacctacctca((((.(((..((.((.((((((((&)))))))))).))..)))))))	-15.5	0.9912094304815581
            #let-7a	NM_145257	160	177	1	19	TATTGGTCACCTACCTC	((((..((..(((((((&)))))))..))..).)))	-13.5	0.9980730598432394
            #attention, this file does not contain the complete mirna name, need add mmu, or hsa acording to species
            my ($miRNA,$refseq_id,$starting,$stop_pos,$strand,$chrome,$match,$detail,$energy,$targetspy_score)= split("\t", $line);            
            if ($species eq "human"){$miRNA='hsa-'.$miRNA;}else{$miRNA='mmu-'.$miRNA;}
            my $mirna_accession;
            if(exists $mirAlias2Accession_hash{$miRNA}){
                $mirna_accession = $mirAlias2Accession_hash{$miRNA};
            }else{ print "WRONG miRNA name:$miRNA\n";next;}
            $targetspy_hash{$mirna_accession}{$refseq_id}{'location'}=$stop_pos.':'.$strand.':'.$chrome.':'.$match.':'.$detail;
            $targetspy_hash{$mirna_accession}{$refseq_id}{'energy'}=$energy;
            $targetspy_hash{$mirna_accession}{$refseq_id}{'targetspyScore'}= $targetspy_score;
        }#if($line){
    }#while(@all){
    my $targetspy_hash_file= $data_dir.'targetspy/'.$tax_ID.'_'.$sensitivity.'_'.$seed.'_targetspy.hash';
    store \%targetspy_hash,$targetspy_hash_file;
    return ($targetspy_hash_file);
}#sub targetspy

######==END_SECTION_TARGETSPY#############################################################
######################HGNC_history###########################################################
    #this part is for the mapping of HGNC old symbol to the newest symbol because of TARBASE, since it contain many old Gene symbol
    #the mapping file was generated from the BIOMART
    #http://www.genenames.org/biomart/martview/dab0043cad05015e7371f7ee48182b39
    #stored as annotation/HGNC_history.txt
sub hgncHistory(){
    my $hgncHistory_file=$data_dir.'annotations/HGNC_history.txt';
    my %hgnc_history_hash=();#all the documented mir alias as key to the accession
    if (open (HGNCHISTORY,$hgncHistory_file)){
        @all=<HGNCHISTORY>;
        shift(@all);
        close(HGNCHISTORY);
    }else{
        print "Can't open file $hgncHistory_file\n";
    }
    while(@all){
        $line = shift @all;
        print ($line);
        chomp ($line);
        #EGFL2~withdrawn	
        #NBN	NBS, NBS1
        my ($current_HGNC,$previous_HGNCs)=split("\t",$line);
        next if ($current_HGNC=~/withdrawn/);#check if this id already retired
        my @previous_HGNC_ids = split(", ",$previous_HGNCs);
        push @previous_HGNC_ids,$current_HGNC;
        foreach my $i(@previous_HGNC_ids){
            $hgnc_history_hash{$i} = $current_HGNC;
            print "$i: $current_HGNC\n";
        }
    }
    store \%hgnc_history_hash, $data_dir.'annotations/hgnc_history.hash';
    return (\%hgnc_history_hash);
}#sub hgnc_history
##############################END_HGNC_history##################################################

=begin  BlockComment  # BlockCommentNo_1

###############BEGIN_TARBASE########################################################
#This part temperary skipped that there is no reply from diana yet
#However the code was generated based on older version
sub tarbase(){
    my ($species)= @_;#define the species(human, mouse)
    $species//= 'mouse';
    my $tarbase_file = $data_dir.'tarbase/TarBase_V5.0.txt';
    my $tax_ID = "10090";# by default genes are from mouse
    my $organism_pre="mmu-";
    if ($species eq "human"){$tax_ID = "9606";$organism_pre="hsa-";}
    my %hgnc_history_hash=%{retrieve($data_dir.'annotations/hgnc_history.hash')};
    my %genename2refseqs_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_rnaGeneName2refseqs.hash')};
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    my %ensg2refseqs_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_ensg2refseqs.hash')};
    my $debug_file =$data_dir.'tarbase/debug_'.$species.'.file';
    open (DEBUG,'>>'.$debug_file)||die "Could not open file $!";
    print "try read in  tarbase target file : $tarbase_file\n";
    if(open(TARBASE,$tarbase_file)){
        @all = <TARBASE>;
        shift(@all);
        close(TARBASE);       
    }else{print "Can't open tarbase file $tarbase_file\n";}
    my %tarbase_hash = ();#$tarbase_hash{mirna_accession}{target_uni_id} 
    while(@all){
        $line = shift(@all);
        chomp($line);        
        if($line){
            #Id	Id_V4	Data_Type	Support_Type	Organism	miRNA	HGNC_Symbol	Gene	Isoform	Ensembl	Chr_loc	MRE	S_S_S	I_S	D_S	Validation	Paper	Target_seq	miRNA_seq	Seq_location	PMID	KEGG	Protein_type	Differentially_expressed_in	Pathology_or_Event	Mis_Regulation	Gene_Expression	Tumour_Involvement	Bibliographic_Notes	Cell_Line_Used	HGNC_ID	SwissProt
            #106	2	Unknown	FALSE	Human	let-7b	n_a	RPIA	_	ENSG00000153574 	chr2	88830465-88889707	Given	No	n_a	in vitro reporter gene assay (Luciferase)	Class 1a	Kiriakidou et al, 2004	n_a	UGAGGUAGUAGGUUGUGUGGUU	n_a	15131085			n_a	n_a	n_a	n_a	n_a		n_a	n_a	n_a
            #107	4	mRNA repression	TRUE	Human	let-7b	LIN28	Lin28	_	ENSG00000131914	chr1	26421410-26440350	Given	Yes	n_a	in vitro reporter gene assay (Luciferase)	Class 1a	Kiriakidou et al, 2004	n_a	UGAGGUAGUAGGUUGUGUGGUU	n_a	15131085	hsa+79727		n_a	n_a	n_a	n_a	carcinoma		n_a	15986	Q9H9Z2        

            my ($id, $id_v4, $data_type, $support_type, $org, $mirna_family, $hgnc_symbol, $gene_symbol, $isoform, $ensembl_gene_id, $chr_loc, $mre, $s_s_s, $i_s, $d_s, $validation, $paper, $target_seq, $mirna_seq, $seq_location, $pmid, $kegg, $protein_type, $differentially_expressed_in, $pathology_or_event, $mis_regulation, $gene_expression, $tumor_involvement, $bibliographic_notes, $cell_line_used, $hgnc_id, $swissprot)= split("\t", $line);            
            if ($org=~/$species/i){
                #this line is evidence in target organism
            }else{
                print "this $org is not target organism\n";
                next;#don't further process this line
            }
            my @mirna_accessions;
            my @refseqs;
            my $mirnaAliasCheck;

            foreach my $mirAlias (sort keys %mirAlias2Accession_hash){
                $mirnaAliasCheck.="$mirAlias:$mirna_family\n";
                if ($mirAlias =~/$mirna_family/i and (not $mirAlias=~/$mirna_family\d/i) and $mirAlias=~/$organism_pre/){#make sure that it contains pattern
                    print "there is a match: $mirAlias:$mirna_family\n";#this mirAlias contain the family partbut not lie let-75# may need further confirmation
                    push @mirna_accessions,$mirAlias2Accession_hash{$mirAlias};
                }# if ($mirAlias =~/$mirna_family/i){#this mirAlias contain the family part# may need further confirmation
            }
            @mirna_accessions = uniq(@mirna_accessions);            
            print DEBUG "$line\n";
            print DEBUG (Dumper @mirna_accessions);
            #if (not @mirna_accessions){print $mirnaAliasCheck;}
            #need a match between the gene symbol and refseqs with ensembl gene ID
            if(exists $ensg2refseqs_hash{$ensembl_gene_id}){
                push @refseqs, @{$ensg2refseqs_hash{$ensembl_gene_id}};
            }
            if(exists $hgnc_history_hash{$hgnc_symbol}){
                if(exists $genename2refseqs_hash{$hgnc_history_hash{$hgnc_symbol}}){
                push @refseqs, @{$genename2refseqs_hash{$hgnc_history_hash{$hgnc_symbol}}};
            }
            }else{
                print "$hgnc_symbol is not in history\n";
            }
            @refseqs= uniq(@refseqs);
            next if(not @refseqs);
            print DEBUG Dumper @refseqs;
            foreach my $mirna_accession(@mirna_accessions){
                foreach my $refseq_id(@refseqs){
                    $tarbase_hash{$mirna_accession}{$refseq_id}{'Validation'}=$validation;
                    $tarbase_hash{$mirna_accession}{$refseq_id}{'PMID'}=$pmid;
                    $tarbase_hash{$mirna_accession}{$refseq_id}{'detail'}=$line;                    
                }            
            }
        }#if($line){
    }#while(@all){
    #print Dumper \%tarbase_hash;
    close (DEBUG);
    my $tarbase_hash_file= $data_dir.'tarbase/'.$tax_ID.'_tarbase.hash';
    store \%tarbase_hash,$tarbase_hash_file;
    return ($tarbase_hash_file);
}#sub tarbase
################END_TARBASE######################################################

=end    BlockComment  # BlockCommentNo_1

=cut

###############BEGIN_TARBASE########################################################
#This part is based on the new diana TARBASE 6.0!!!
sub tarbase(){
    my ($species)= @_;#define the species(human, mouse)
    $species//= 'mouse';
    my $tarbase_file = $data_dir.'tarbase/TarBase_V6.0.csv';
    my $tax_ID = "10090";# by default genes are from mouse
    my $organism_pre="mmu";
    if ($species eq "human"){$tax_ID = "9606";$organism_pre="hsa";}
    my %hgnc_history_hash=%{retrieve($data_dir.'annotations/hgnc_history.hash')};
    my %genename2refseqs_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_rnaGeneName2refseqs.hash')};
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    my %ensg2refseqs_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_ensg2refseqs.hash')};
    my $debug_file =$data_dir.'tarbase/debug_'.$species.'.file';
    open (DEBUG,'>'.$debug_file)||die "Could not open file $!";
    print "try read in  tarbase target file : $tarbase_file\n";
    if(open(TARBASE,$tarbase_file)){
        @all = <TARBASE>;
        shift(@all);
        close(TARBASE);       
    }else{print "Can't open tarbase file $tarbase_file\n";}
    my %tarbase_hash = ();#$tarbase_hash{mirna_accession}{target_uni_id} 
    while(@all){
        $line = shift(@all);
        chomp($line);        
        if($line){
#miRNA-mimat,ensgid,gene-name,reporter_gene,nothern_blot,western_blot,qPCR,proteomics,microarray,sequencing,degradome_seq,other
#MIMAT0000001,F11A1.3,daf-12,POSITIVE,\N,\N,\N,\N,\N,\N,\N,\N
#MIMAT0000001,F11A1.3,daf-12,\N,\N,\N,\N,\N,\N,\N,\N,UNKNOWN
#MIMAT0000001,NM_063774(cel),NM_063774(cel),\N,\N,\N,\N,\N,\N,\N,\N,UNKNOWN
#MIMAT0000001,ZK792.6,let-60,POSITIVE,\N,\N,\N,\N,\N,\N,\N,\N
            my ($mirna_accession,$ensgid,$genename,$reporter_gene,$nothern_blot,$western_blot,$qPCR,$proteomics,$microarray,$sequencing,$degradome_seq,$other)= split(",", $line);   
            # $miRNAmimat is fine, always mature miRNA accession, $ENSGID could be ENSG, or genename, or refseq id ,or id+(organism)
            my @refseqs;    
            my ($ensembl_gene_id, $org)= split(/\(/, $ensgid );
            #chop $org;
            print "$genename\n";
            my @genenames=split(/\(/,$genename);
            my $hgnc_symbol=shift(@genenames);
            print "$hgnc_symbol \n";
            if(exists $ensg2refseqs_hash{$ensembl_gene_id}){
                push @refseqs, @{$ensg2refseqs_hash{$ensembl_gene_id}};
            }
            elsif(exists $genename2refseqs_hash{$hgnc_symbol} and $org =~ m/$organism_pre/){
                push @refseqs, @{$genename2refseqs_hash{$hgnc_symbol}};
            }
            elsif(exists $hgnc_history_hash{$hgnc_symbol}){
                if(exists $genename2refseqs_hash{$hgnc_history_hash{$hgnc_symbol}}){
                push @refseqs, @{$genename2refseqs_hash{$hgnc_history_hash{$hgnc_symbol}}};
            }
            }elsif($org =~m/$organism_pre/ and $ensembl_gene_id=~ m/NM_|XM_/){

                push @refseqs,$ensembl_gene_id;#in such case the id is acturally a refseq id
                print "$ensembl_gene_id,this is a refseq id \n";
            }
            
            else{
                if($org =~ m/$organism_pre/){
                print "can't find match $ensgid,$genename\n";                
                print DEBUG "$line\n";}
                
            }
            @refseqs= uniq(@refseqs);
            next if(not @refseqs);
            foreach my $refseq_id(@refseqs){
                $tarbase_hash{$mirna_accession}{$refseq_id}{'detail'}=$ensgid."\t".$genename;
                my @list=($reporter_gene,$nothern_blot,$western_blot,$qPCR,$proteomics,$microarray,$sequencing,$degradome_seq,$other);
                if(not my @matched= grep{/NEGATIVE|POSITIVE|UNKNOWN/} @list){#there is no details about the regulation direction
                     $tarbase_hash{$mirna_accession}{$refseq_id}{'regulation'}='NA';
                     $tarbase_hash{$mirna_accession}{$refseq_id}{'evidence'}='NA';
                }else{
                    my @index=grep{$list[$_]=~/NEGATIVE|POSITIVE|UNKNOWN/} 0 ..$#list;
                    my $index=shift @index;
                    my $matched=shift @matched;
                    $tarbase_hash{$mirna_accession}{$refseq_id}{'regulation'}=$matched;
                    given($index) {
                       when (0) {$tarbase_hash{$mirna_accession}{$refseq_id}{'evidence'}='Reporter Assay';}
                       when (1) {$tarbase_hash{$mirna_accession}{$refseq_id}{'evidence'}='Nothern Blot';}
                       when (2) {$tarbase_hash{$mirna_accession}{$refseq_id}{'evidence'}='Western Blot';}
                       when (3) {$tarbase_hash{$mirna_accession}{$refseq_id}{'evidence'}='qPCR';}
                       when (4) {$tarbase_hash{$mirna_accession}{$refseq_id}{'evidence'}='proteomics';}
                       when (5) {$tarbase_hash{$mirna_accession}{$refseq_id}{'evidence'}='microarray';}
                       when (6) {$tarbase_hash{$mirna_accession}{$refseq_id}{'evidence'}='sequenceing';}
                       when (7) {$tarbase_hash{$mirna_accession}{$refseq_id}{'evidence'}='degradome_seq';}
                       when (8) {$tarbase_hash{$mirna_accession}{$refseq_id}{'evidence'}='otherMethods';}
                       default {print "You are just making it up\n";}                        
                    }
                }                   
            }            
            
        }#if($line){
    }#while(@all){
    print Dumper \%tarbase_hash;
    close (DEBUG);
    my $tarbase_hash_file= $data_dir.'tarbase/'.$tax_ID.'_tarbase.hash';
    store \%tarbase_hash,$tarbase_hash_file;
    return ($tarbase_hash_file);
}#sub tarbase
################END_TARBASE######################################################
######==BEGIN_PACCMIT===##################################################################
#paccmit('','_access','_2_5')#
#paccmit('','','')#no filter
#paccmit('','_access','')#access only filter
#paccmit('_cons','','')#conservation filter only
#paccmit('_cons','_access','')#conservation+acess
#paccmit()#most strict, with acessibility, restriction and conservation filter
#http://lcpt.epfl.ch/MicroRNA_target_predictions
#http://lcpt.epfl.ch/PACCMIT-CDS
#this programm didn't offer prediction result for mouse so need a otholog conversion or download the ensembl mouse to run the program,program was also downloaded
#matching between the ENSG and 
#this predictions has been obtained by restricting the location of the mucleation region also optimized Pcutoff of 0.2 was used 
sub paccmit(){#this program only have human data available , and different format for CDS an orgion
    my ($con,$accessibility,$restriction)= @_;#$con if it is conservatived, $considertiong accessiblity or not $location restriction 2-5
    my $species='human';
    $accessibility//='_access';
    $restriction//= '_2_5';#restrcted region or ''
    $con//='_cons';# by default,take conservative or ''   
    my $paccmit_file = $data_dir.'paccmit/paccmit'.$accessibility.$con.$restriction.'.txt';
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my %ensg2refseqs_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_ensg2refseqs.hash')};
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    print "try read in  paccmit target file : $paccmit_file\n";
    if(open(PACCMIT,$paccmit_file)){
        @all = <PACCMIT>;#this file has no headline
        shift(@all);#delete head line
        shift(@all);
        close(PACCMIT);       
    }else{print "Can't open PACCMIT file $paccmit_file\n";}
    my %paccmit_hash = ();#$microT4_hash{mirna_accession}{target_uni_id} 
    while(@all){
        $line = shift(@all);
        chomp($line);        
        if($line){
            #Gene            miRNA           c_filter log(P_SH)
            #--------------------------------------------------
            #ENSG00000169181	hsa-miR-3134	68	-81.754
            #ENSG00000148735	hsa-miR-4747-5p	41	-56.963
            #ENSG00000148735	hsa-miR-5196-5p	41	-56.963
                  
            my ($ensg_id,$miRNA,$c_filter,$logP)= split("\t", $line);            
            my $mirna_accession;
            my @refseq_ids;
            if(exists $mirAlias2Accession_hash{$miRNA}){
                $mirna_accession = $mirAlias2Accession_hash{$miRNA};
            }else{ print "WRONG miRNA name:$miRNA\n";next;}            
            if(!exists $ensg2refseqs_hash{$ensg_id}){
                 print "can't find refseqs : $ensg_id\n";
                 next;
           }else{
               push @refseq_ids,@{$ensg2refseqs_hash{$ensg_id}};
               uniq(@refseq_ids);
               foreach my $refseq_id(@refseq_ids){
                   $paccmit_hash{$mirna_accession}{$refseq_id}{'c_filter'}=$c_filter;
                   $paccmit_hash{$mirna_accession}{$refseq_id}{'logP'}=$logP;
               }            
            }#else        
        }#if($line){
    }#while(@all){
    my $paccmit_hash_file= $data_dir.'paccmit/'.$tax_ID.$accessibility.$con.$restriction.'_paccmit.hash';
    store \%paccmit_hash,$paccmit_hash_file;
    return ($paccmit_hash_file);
}#sub paccmit
######==END_SECTION_PACCMIT#############################################################
######==BEGIN_PACCMIT_CDS===##################################################################

#http://lcpt.epfl.ch/PACCMIT-CDS
#this programm didn't offer prediction result for mouse so need a otholog conversion or download the ensembl mouse to run the program,program was also downloaded
#matching between the ENSG and 
#this predictions has been obtained by restricting the location of the mucleation region also optimized Pcutoff of 0.2 was used 
#paccmit_cds();
#paccmit_cds('');
sub paccmit_cds(){#this program only have human data available , and different format for CDS an orgion
    my ($con)= @_;#$con if it is conservatived, $considertiong accessiblity or not $location restriction 2-5
    my $species='human';
    $con//='_cons';# by default,take conservative or ''   
    my $paccmit_cds_file = $data_dir.'paccmit/paccmit_cds_'.$species.$con.'.txt';
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}
    my %ensembl_tr2refseq_hash=%{retrieve($data_dir.'annotations/'.$tax_ID.'_ensembl_transcriptID2refseqID.hash')};
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    my $missing_GeneName_file=  $data_dir.'paccmit/paccmit_cds_'.$species.$con.'_missingGeneNames.txt';
    print "try read in  paccmit_cds target file : $paccmit_cds_file\n";
    if(open(PACCMITCDS,$paccmit_cds_file)){
        @all = <PACCMITCDS>;#this file has no headline
        shift(@all);#delete head line
        shift(@all);
        close(PACCMITCDS);       
    }else{print "Can't open paccmit_cds file $paccmit_cds_file\n";}
    my %paccmit_cds_hash = ();#$microT4_hash{mirna_accession}{target_uni_id} 
    while(@all){
        $line = shift(@all);
        chomp($line);        
        if($line){
            #ENST00000225964	17	0.0000000000	hsa-miR-3144-5p
            #ENST00000225964	18	0.0000000000	hsa-miR-634
            #ENST00000225964	22	0.0000000000	hsa-miR-29a-3p,hsa-miR-29b-3p,hsa-miR-29c-3p
            my ($enst_id,$c_cons,$P_SH,$miRNAs)= split("\t", $line); 
            my @mirnas =split(",",$miRNAs);
            my @mirna_accession;
            my $refseq_id;
            foreach my $miRNA(@mirnas){
                if(exists $mirAlias2Accession_hash{$miRNA}){
                    push @mirna_accession,$mirAlias2Accession_hash{$miRNA};
                }else{ print "WRONG miRNA name:$miRNA\n";next;} 
            }           
            if(!exists $ensembl_tr2refseq_hash{$enst_id}){
                 print "can't find ensembl transcript : $enst_id\n";#this part still have some problem! first some id retired, nolonger in DB#second some could not be mapped to refseq
                 open (MISSING_GENE,'>>'.$missing_GeneName_file)||die "Could not open file $!";
                 print MISSING_GENE "$enst_id\n";
                 close(MISSING_GENE);
                 next;
           }else{
               $refseq_id=$ensembl_tr2refseq_hash{$enst_id};          
            }#else
            foreach(@mirna_accession){
                $paccmit_cds_hash{$_}{$refseq_id}{'c_cons'}=$c_cons;
                $paccmit_cds_hash{$_}{$refseq_id}{'P_SH'}=$P_SH;
            }
        }#if($line){
    }#while(@all){
    my $paccmit_cds_hash_file= $data_dir.'paccmit/'.$tax_ID.$con.'_paccmit_cds.hash';
    store \%paccmit_cds_hash,$paccmit_cds_hash_file;
    return ($paccmit_cds_hash_file);
}#sub paccmit_cds
######==END_SECTION_PACCMIT_CDS#############################################################
######==BEGIN_TARGETMINER===##################################################################
#TargetMiner was used to search genome-wide potential conserved targets as in miRanda, TargetScan and PicTar for a phastcons cutoff value of 0.57
#observed that the average number of targets predicted per miRNA (TPM) for TargetMiner is moderate (1084.5) as compared with several other methods

sub targetminer(){#this program only have human data available , and different format for CDS an orgion
    my $species='human';
    my $targetminer_file = $data_dir.'targetminer/targetminer_'.$species.'.txt';
    my $tax_ID = "10090";# by default genes are from mouse
    if ($species eq "human"){$tax_ID = "9606";}   
    my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    print "try read in  targetminer target file : $targetminer_file\n";
    if(open(TARGETMINER,$targetminer_file)){
        @all = <TARGETMINER>;#this file has no headline
        shift(@all);#delete head line
        close(TARGETMINER);       
    }else{print "Can't open targetminer file $targetminer_file\n";}
    my %targetminer_hash = ();#$microT4_hash{mirna_accession}{target_uni_id} 
    while(@all){
        $line = shift(@all);
        chomp($line);        
        if($line){
            #miRNA_id	miRNA_sequence	mRNA	Chromosome	6mer (count)	6mer (position)	7mer-A1 (count)	7mer-A1 (position)	7mer-m8 (count)	7mer-m8 (position)	8mer (count)	8mer (position)
            #hsa-let-7c	 UGAGGUAGUAGGUUGUAUGGUU	NM_000025	chr8	1	456	-	-	1	1158	-	-            
            my ($miRNA_id,$miRNA_sequence,$mRNA,$Chromosome,$mer6_count,$mer6_pos,$mer7_A1_count,$mer7_A1_pos,$mer78_count,$mer78_pos,$mer8_count,$mer8_pos)= split("\t", $line); 
            $miRNA_id=~s/^\s+|\s+$//g;
            my $mirna_accession;
            if(exists $mirAlias2Accession_hash{$miRNA_id}){
                $mirna_accession = $mirAlias2Accession_hash{$miRNA_id};
            }else{ print "WRONG miRNA name:$miRNA_id\n";
                print "$line\n";
                next;}

            $targetminer_hash{$mirna_accession}{$mRNA}=1;
                       
        }#if($line){
    }#while(@all){
    my $targetminer_hash_file= $data_dir.'targetminer/targetminer_'.$tax_ID.'_targetminer.hash';
    store \%targetminer_hash,$targetminer_hash_file;
    return ($targetminer_hash_file);
}#sub targetminer 
############END_SECTION_TARGETMINER###################################################################
    #my(%human_microcosm_hash,%human_microcosm_detail_hash)=&microcosm('human_microcosm_mirna2target.txt');
#my(%mouse_microcosm_hash,%mouse_microcosm_detail_hash)=&microcosm('mouse_microcosm_mirna2target.txt');

    # my %microcosm_hash=%{retrieve($data_dir.'annotations/human_microcosm_mirna2target.txt.hash')};
#my %microcosm_detail_hash=%{retrieve($data_dir.'annotations/human_microcosm_mirna2target.txt.detail_hash')};
#my $pita_hash_file = &pita('human','no_flank','TOP');#$species='human','mouse';$flank='no_flank','flank';$stringency='TOP','ALL'
#my $human_microrna_hash_file= &microrna_org('human','conserved');#$species='human','mouse',$conservation='conserved','non_conserved'

#my %human_targetscan_hash = %{ retrieve( &targetScan('human') )} ;
#my $human_miRDB_hash_file= &mirdb('human');
#my $human_mirnamap_hash_file=&mirnamap('human',2,2,1);#my ($species,$creterion1,$creterion2,$creterion3) = @_;$species='human','mouse'$creterion1:Nr, the minimum of reported tools;$creterion2_minimum sites nr between target and miR;$accessible or not (1 0)is this site
#my $human_reptar_hash_file= &reptar('human');#($species,$mfe_upper,$nmfeb_lower,$rf_lower) #define the species(the profile to readin), the minimal free energy of inding , #the normalized minimal free energy of binding, and the repeating motifs minimum
#my $human_eimmo_hash_file=&eimmo('human');



#################testing zone##############################
#
# my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
# my %ens_transcript2uni_hash = %{retrieve($data_dir.'annotations/ens_transcript2uni.hash')};
#  my ($rnaGeneName2refseqs_hash_file,$refseq2_rnaGeneNames_hash_file)=&refseq2rnaGeneName('human');
#my $pictar_hash_file= &pictar('human','most_conserved'); 

    #my(%human_microcosm_hash,%human_microcosm_detail_hash)=&microcosm('human_microcosm_mirna2target.txt');
#my(%mouse_microcosm_hash,%mouse_microcosm_detail_hash)=&microcosm('mouse_microcosm_mirna2target.txt');

    # my %microcosm_hash=%{retrieve($data_dir.'annotations/human_microcosm_mirna2target.txt.hash')};
#my %microcosm_detail_hash=%{retrieve($data_dir.'annotations/human_microcosm_mirna2target.txt.detail_hash')};
#my $pita_hash_file = &pita('human','no_flank','TOP');#$species='human','mouse';$flank='no_flank','flank';$stringency='TOP','ALL'
#my $human_microrna_hash_file= &microrna_org('human','conserved');#$species='human','mouse',$conservation='conserved','non_conserved'
#my ($ensembl_tr2refseq_hash_file,$refseq2ensembl_tr_hash_file)=&refseq2ensembl('human'); 


#my %human_targetscan_hash = %{ retrieve( &targetScan('human') )} ;
#my $human_miRDB_hash_file= &mirdb('human');
#my $human_mirnamap_hash_file=&mirnamap('human',2,2,1);#my ($species,$creterion1,$creterion2,$creterion3) = @_;$species='human','mouse'$creterion1:Nr, the minimum of reported tools;$creterion2_minimum sites nr between target and miR;$accessible or not (1 0)is this site
#my $human_reptar_hash_file= &reptar('human');#($species,$mfe_upper,$nmfeb_lower,$rf_lower) #define the species(the profile to readin), the minimal free energy of inding , #the normalized minimal free energy of binding, and the repeating motifs minimum
#my $human_eimmo_hash_file=&eimmo('human');
#my ($entrez2refseqs_hash_file,$refseq2entrez_hash_file)=&refseq2entrez('human');
#my ($ensg2refseqs_hash_file,$refseq2ensg_hash_file)=&ensg2refseqs('human');
#my ($ensg2refseqs_hash_file,$refseq2ensg_hash_file)=&ensg2refseqs('mouse');

#&mirtar();
##############################################################
#&mirtarbase();
#&mirtarbase('human');
#&mirtarbase(,'SE');
#&mirtarbase('human','SE');
#########################################################
#&microT4();#mouse, microt4,0.7
#&microT4('human',);#human,micro4,0.7
#&microT4('mouse','CDS');#mouse,cds microt, 0.7
#&microT4('human','CDS');#human,cds microT,0.7
#####################################################################
#&hgncHistory();#must be excuted before tarbase
#my ($tarbase_hash_file)=&tarbase('human');
#my ($tarbase_hash_file)=&tarbase('mouse');
#################################################################
#&paccmit('','_access','_2_5');#access+restriction
#&paccmit('','','');#no filter
#&paccmit('','_access','');#access only filter
#&paccmit('_cons','','');#conservation filter only
#&paccmit('_cons','_access','');#conservation+acess
#&paccmit();#most strict, with acessibility, restriction and conservation filter
######targetminer###############################################
#&paccmit_cds();
#&paccmit_cds('');
########################
#&targetminer();
#########################################################################
   
####################END_of_testing zone#####################################################
#&mirBaseAlias();
#&refseq2rnaGeneName('human');
#&refseq2rnaGeneName('mouse');
#&refseq2ensembl('human');
#&refseq2ensembl('mouse');
#&ensg2refseqs('human');
#&ensg2refseqs('mouse'); 
#my ($ensembl_tr2refseq_hash_file,$refseq2ensembl_tr_hash_file)=&refseq2ensembl('human'); 
#my ($entrez2refseqs_hash_file,$refseq2entrez_hash_file)=&refseq2entrez('human');
#my ($ensg2refseqs_hash_file,$refseq2ensg_hash_file)=&ensg2refseqs('human');
#my ($entrez2refseqs_hash_file,$refseq2entrez_hash_file)=&refseq2entrez('mouse');
#&refseq2entrez('human');
#&refseq2entrez('mouse');
#&microcosm('human');
#&microcosm('mouse');
#&targetscan('human','broadly');
#&targetscan('human','conserved');
#&targetscan('human','nonconserved');
#&targetscan('mouse','broadly');
#&targetscan('mouse','conserved');
#&targetscan('mouse','nonconserved');
    #############################################
#&pictar('mouse','most_conserved');
#&pictar('mouse','mammal_conserved');
#&pictar('human','most_conserved');
#&pictar('human','mammal_conserved');
    ###################finished however there are many problem about the pimatching of gene symbols 
#&pita('human','no_flank','TOP');
#&pita('human','flank','TOP');
#&pita('human','no_flank','ALL');
#&pita('human','flank','ALL');
#&pita('mouse','no_flank','TOP');
#&pita('mouse','flank','TOP');
#&pita('mouse','no_flank','ALL');
#&pita('mouse','flank','ALL');
    ###################finished however there are many problem about hte matching of gene symbols################
#&microrna_org('human','conserved');
#&microrna_org('human','non_conserved');
#&microrna_org('mouse','conserved');
#&microrna_org('mouse','non_conserved');
    #########################finished
#&mirdb('human');
#&mirdb('mouse');
    ############finished###########
##&mirnamap('human',1,1,1);
##&mirnamap('human',1,2,1);
#&mirnamap('human',2,1,1);
##&mirnamap('human',2,2,1);
##&mirnamap('mouse',1,1,1);
##&mirnamap('mouse',1,2,1);
#&mirnamap('mouse',2,1,1);
##&mirnamap('mouse',2,2,1);
#&reptar('human');
#&reptar('mouse');
#&eimmo('human',undef,undef,'CDS');
&eimmo('mouse',undef,undef,'');
#&eimmo('human',undef,undef,'');
#&eimmo('human',0,undef,'');
&eimmo('mouse',0,undef,'');
######################
&microT4();#mouse, microt4,0.7
&microT4('human',undef,undef);#human,micro4,0.7
&microT4('mouse','cds',undef);#mouse,cds microt, 0.7
&microT4('human','cds',undef);
&mirtarbase();
&mirtarbase('human',undef);
#&mirtarbase(undef,'SE');
#&mirtarbase('human','SE');
&mirtar('3utr');
#&mirtar('5utr');
#&mirtar('all');
#&mirtar('cds');
##&targetrank('human',1);
##&targetrank('human',0);
##&targetrank('mouse',1);
##&targetrank('mouse',0);
##&targetspy('human','sens','seed');
&targetspy('human','sens','all');
&targetspy('human','spec','seed');
&targetspy('mouse','sens','seed');
&targetspy('mouse','sens','all');
&targetspy('mouse','spec','seed');
&tarbase('human');
&tarbase('mouse');
&paccmit('','_access','_2_5');#
&paccmit('','','');#no filter
&paccmit('','_access','');#access only filter
&paccmit('_cons','','');#conservation filter only
&paccmit('_cons','_access','');#conservation+acess
&paccmit();#most strict, with acessibility, restriction and conservation filter

&paccmit_cds();
&paccmit_cds('');

&targetminer();#this one is human only

