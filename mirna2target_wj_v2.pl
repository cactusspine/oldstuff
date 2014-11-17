#!/usr/bin/perl -w

use strict;
use Storable;
use Data::Dumper qw(Dumper);

my $project_dir = "/home/jiali/workspace/projects/vs/";
my $data_dir = $project_dir . "data/";
my $result_dir = $project_dir . "result/";
my %all_resulted_uni_ids_hash = ();
my $result_file = $result_dir . "mirna_targets.txt";
#close(RESULT_FILE);

my $line = "";
my @all = ();



###miRBASE alias matching related==miRBASE alias===
=begin
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

##END===miRbase_Alias== 

##################################################################################################################################################
## ens_transcript2uni_hash file creation##############
my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
my $ensemblTranscriptID = $data_dir . "annotations/EnsmblTranscript_id.txt";
close(TRAIL_ENSEMBLTID_FILE);
my %ens_transcript2uni_hash = ();
if (open(TRAIL_ENSEMBLTID_FILE, $ensemblTranscriptID)){
    @all = <TRAIL_ENSEMBLTID_FILE>;
    close(TRAIL_ENSEMBLTID_FILE);
}
while(@all){
    $line = shift @all;
    chomp($line);
    if ($line) {
        my $ensg=$line;
        #my ($mouse_ensg, $mouse_refseq) = split("\t", $line);
        push @{$ens_transcript2uni_hash{$ensg}}, $ensg;#last one is value, key here is uni_id=ensg

    } # if ($line){
}# while(@all){
foreach my $uni_id (sort keys %ens_transcript2uni_hash){
    foreach my $cell(@{$ens_transcript2uni_hash{$uni_id}}){
        print "$uni_id,$cell\n";    
    }
}#foreach %ens_transcript2uni_hash to print#testing
store\%ens_transcript2uni_hash,$data_dir.'annotations/ens_transcript2uni.hash';
###END_ens_ranscript2uni_hash file creation##### 


### Microcosm hash part ==microcosm==#########################
#Microcosm_V5 no updating is using miRanda algorithem to identify potential BS find highly complementary at 5 end(0-100)->Vienna routine for thermodynamic stability->check ortholog conservation.
##%this part need to be packed and reused for human, the place format are the same date.dir/human_microcosm_mirna2target.txt
sub microcosm(){ #input the name of file to be parsed
    my ($microcosm_mirna2target_file) =shift(@_);   
    #==need revise, add input for different organism
    $microcosm_mirna2target_file=$data_dir."microcosm/".$microcosm_mirna2target_file;#
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
                my @tmp_uni_id = @{$ens_transcript2uni_hash{$ensembl_transcript_id}};
                push @uni_id, @tmp_uni_id;}
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
                        my $old_score = $microcosm_detail_hash{$mirna_accession}{$uni_id}{"score"};
                        if ($score > $old_score) {
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"chr"} = $chr;
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"start"} = $start;
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"end"} = $end;
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"strand"} = $strand;
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"score"} = $score;
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"pvalue_og"} = $pvalue_og;
                            $microcosm_detail_hash{$mirna_accession}{$uni_id}{"gene_symbol"} = $external_name;
                            $microcosm_hash{$mirna_accession}{$uni_id} = 1;
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
                        $microcosm_hash{$uni_id}{$mirna_accession} = 1;
                    } # else {# if there is no entry of this miRNA and uni_id pair yet
                }#if this mir_id find match in mir_dictionary if ($mirAlias2Accession_has{$mirna_id}){                   
            } # foreach my $uni_id (@uni_id) { 
         } # if (@uni_id) {
        else { print "problem - microcosm - no uni_id - line : $line\n";} # else {#no uni_id match for this gene
    } # if ($line) {
}# while(@all){
store \%microcosm_detail_hash, $microcosm_mirna2target_file.'.detail_hash';
store \%microcosm_hash,$microcosm_mirna2target_file.'.hash';
return (\%microcosm_hash,\%microcosm_detail_hash)
}
#genes could be matched to ensembl_transcript id

###END==microcosm==
=end

=cut
sub nonRedundantList{
    my $val = shift;
    my @val=@$val;                                                
    my %cpt;
    foreach (@val){                                               
        $cpt{$_}++;                                                 
    }                                                             
    return keys %cpt;
} #sub redundantList{
##############################################################
### ==BEGIN_TARGETSCAN======
##need to a dictionary between the miRfamily&MIRbaseID
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
sub targetScan(){
    my $species = shift(@_);
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
            if($species_ID eq $tax_ID){
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
            next if ($total_context_score eq "NULL" or $aggregatePCT eq "NULL" );
            if($tncs>=1 and ($total_context_score +0.3)<=0 and $aggregatePCT>=0.75){#at least one consevative site match, tcns<-0.3,aPCT>=0.75
                if(exists $mirFamily2mirAccession{$miR_family}){
                    foreach my $currentmiRNA (@{$mirFamily2mirAccession{$miR_family}}){
                        $targetscan_hash{$currentmiRNA}{$refseq_Transcript_ID}{"total_Context_Score"}=$total_context_score;
                        $targetscan_hash{$currentmiRNA}{$refseq_Transcript_ID}{"aggregatePCT"}=$aggregatePCT;
                        #==TO__ADD a tracker or lable that there is result in targetScan for this miRNA and refSeq_Transcript_ID & change cutoff to variables
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
    my $targetscan_hash_file=$data_dir.'targetscan/'.$species.'_targetScan_detail_hash';
    store \%targetscan_hash,$targetscan_hash_file;
    return ($targetscan_hash_file);
}#sub targetScan
#END_SECTION_ENDSCTION############################################################################
=begin
#################testing zone##############################
  my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};
    my %ens_transcript2uni_hash = %{retrieve($data_dir.'annotations/ens_transcript2uni.hash')};

    my(%human_microcosm_hash,%human_microcosm_detail_hash)=&microcosm('human_microcosm_mirna2target.txt');
#my(%mouse_microcosm_hash,%mouse_microcosm_detail_hash)=&microcosm('mouse_microcosm_mirna2target.txt');

    my %microcosm_hash=%{retrieve($data_dir.'annotations/human_microcosm_mirna2target.txt.hash')};
#my %microcosm_detail_hash=%{retrieve($data_dir.'annotations/human_microcosm_mirna2target.txt.detail_hash')};
=end

=cut

my %human_targetscan_hash = %{ retrieve( &targetScan('human') )} ;
####################END_of_testing zone#####################################################
