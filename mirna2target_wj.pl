#!/usr/bin/perl -w

use strict;
use Storable;

my $project_dir = "/home/jiali/workspace/projects/vs/";
my $data_dir = $project_dir . "data/";
my $result_dir = $project_dir . "result/";
my %all_resulted_uni_ids_hash = ();
my $result_file = $result_dir . "mirna_targets.txt";
#close(RESULT_FILE);

my $line = "";
my @all = ();

=begin  BlockComment  # BlockCommentNo_2

open(RESULT_FILE, ">$result_file");#open Result matrix to write .. all the scoring from different methods will be stored here


### Kiran pain related ###

### Pain_related_genes
my $pain_related_genes_file = $data_dir . "kiran_exp/Pain_related_genes.txt";#need to be changed at every differerent input 
close(PAIN_GENES_FILE);
my %pain_genes_hash = ();
if (open(PAIN_GENES_FILE, $pain_related_genes_file)){
    @all = <PAIN_GENES_FILE>;
    close(PAIN_GENES_FILE);
    while(@all){
        $line = shift @all;
        chomp($line);
        if ($line) {
            my ($category, $gene_name, $mgi_id) = split("\t", $line);
#print "mgi_id : $mgi_id\n";				
            my $gene_name__category = $gene_name . "____" . $category;
            push @{$pain_genes_hash{$mgi_id}}, $gene_name__category#autoviv key,val pair
        } # if ($line) {
    }# while(@all){
} # if (open(PAIN_GENES_FILE, $pain_related_genes_file)){
##pain_gene_hash {mgi_id,gene_name_category}

### Pain_related mRNAs after profiling ....#%this part is not to tag genes which are involved in pain(in GIST)
my $pain_related_mrna_file = $data_dir . "kiran_exp/mrna_fc1_entrezgene2mgi.txt";
close(PAIN_MRNA_FILE);
my %pain_mrna_hash = ();
if (open(PAIN_MRNA_FILE, $pain_related_mrna_file)){
    @all = <PAIN_MRNA_FILE>;
    close(PAIN_MRNA_FILE);
    while(@all){
        $line = shift @all;
        chomp($line);
        if ($line) {
            my ($mouse_ensembl_gene, $mgi_id, $mgi_symbol, $entrezgene_id) = split("\t", $line);
            $pain_mrna_hash{$mgi_id} = $mgi_symbol;
#print "pain_mrna_hash mgi_id : $mgi_id, mgi_symbol : $mgi_symbol\n";				
        } # if ($line) {
    }# while(@all){
} # if (open(PAIN_MRNA_FILE, $pain_related_mrna_file)){
=end    BlockComment  # BlockCommentNo_2

=cut

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
=end

=cut
##END===miRbase_Alias== 
=begin
### MGI related ....
my $mgi_coordinate_file = $data_dir . "mgi/MGI_Coordinate.rpt"; 
close(MGI);
my %mgi_coordinate_hash = ();
my %mgi2entresgene_from_coordinate = ();
my %genename2uni_hash = ();

if (open(MGI, $mgi_coordinate_file)){
    @all = <MGI>;
    close(MGI);
    $line = shift @all; # get rid off header info

    while(@all){
        $line = shift @all;
        chomp($line);
        if ($line) {
            my ($mgi_id, $marker_type, $gene_symbol, $gene_name, $representative_genome_id, $representative_genome_chromosome, $representative_genome_start, $representative_genome_end, $representative_genome_strand, $representative_genome_build, $entrezgene_id, $NCBI_gene_chromosome, $NCBI_gene_start, $NCBI_gene_end, $NCBI_gene_strand, $ensembl_gene_id, $Ensembl_gene_chromosome, $Ensembl_gene_start, $Ensembl_gene_end, $Ensembl_gene_strand, $vega_gene_id, $VEGA_gene_chromosome, $VEGA_gene_start, $VEGA_gene_end, $VEGA_gene_strand, $UniSTS_gene_chromosome, $UniSTS_gene_start, $UniSTS_gene_end, $MGI_QTL_gene_chromosome, $MGI_QTL_gene_start, $MGI_QTL_gene_end, $miRBase_gene_id, $miRBase_gene_chromosome, $miRBase_gene_start, $miRBase_gene_end, $miRBase_gene_strand, $Roopenian_STS_gene_start, $Roopenian_STS_gene_end) = split("\t", $line);
            $mgi_coordinate_hash{$mgi_id}{'gene_symbol'} = $gene_symbol;
            $mgi_coordinate_hash{$mgi_id}{'gene_name'} = $gene_name;
            $mgi_coordinate_hash{$mgi_id}{'entrezgene_id'} = $entrezgene_id;
            $mgi_coordinate_hash{$mgi_id}{'ensembl_gene_id'} = $ensembl_gene_id;
            $mgi_coordinate_hash{$mgi_id}{'vega_id'} = $vega_gene_id;
            $mgi2entresgene_from_coordinate{$mgi_id}{$entrezgene_id} = 1;
            $mgi2entresgene_from_coordinate{$entrezgene_id}{$mgi_id} = 1;
            $genename2uni_hash{$gene_symbol} = $mgi_id;
        }
    }
}


my $mrk_sequence_file = $data_dir . "mgi/MRK_Sequence.rpt";
close(MGI);
my %mrk_sequence_hash = ();
my %mgi2refseq_mrk_sequence = ();
my %refseq2mgi_mrk_sequence = ();

if (open(MGI, $mrk_sequence_file)){
    @all = <MGI>;
    close(MGI);
    $line = shift @all; # get rid off header info

    while(@all){
        $line = shift @all;
        chomp($line);
        if ($line) {
            my ($mgi_id, $gene_symbol, $status, $type, $gene_name, $cm_position, $chr, $genbank_acc_ids, $unigene_id, $refseq_ids, $vega_transcript_id, $ensemble_transcript_id) = split("\t", $line);
            $mrk_sequence_hash{$mgi_id}{'gene_symbol'} = $gene_symbol;
            $mrk_sequence_hash{$mgi_id}{'gene_name'} = $gene_name;
            my @refseq_ids = split(" ", $refseq_ids);
            foreach my $refseq_id (@refseq_ids) {
                push @{$mgi2refseq_mrk_sequence{$mgi_id}}, $refseq_id;
                push @{$refseq2mgi_mrk_sequence{$refseq_id}}, $mgi_id;
                $genename2uni_hash{$gene_symbol} = $mgi_id;
            }	
        }
    }
}


### mouse ensembl 2 MGI mappings .....
my $ensembl2mgi_file = $data_dir . "ensembl/mouse_ensembl2mgi.txt";
close(ENS2MGI_FILE);
my %ensembl2mgi_hash = ();
my %ens_transcript2mgi_hash = ();
my %ens_gene2mgi_hash = ();
my %mgi2ens_transcript_hash = ();
my %mgi_info_hash = ();
if (open(ENS2MGI_FILE, $ensembl2mgi_file)){
    @all = <ENS2MGI_FILE>;
    close(ENS2MGI_FILE);
    $line = shift @all; # get rid off header info

    while(@all){
        $line = shift @all;
        chomp($line);
        if ($line) {
            my ($ens_gene, $ens_transcript, $ens_protein, $mgi_id, $mgi_symbol, $mgi_des) = split("\t", $line);
            if ($ens_transcript && $mgi_id) {
                push @{$ens_transcript2mgi_hash{$ens_transcript}}, $mgi_id;
                push @{$ens_gene2mgi_hash{$ens_gene}}, $mgi_id;
                push @{$mgi2ens_transcript_hash{$mgi_id}}, $ens_transcript;
            }		

            $ensembl2mgi_hash{$ens_transcript}{"ens_gene"} = $ens_gene;
            push @{$ensembl2mgi_hash{$ens_transcript}{"mgi_id"}}, $mgi_id;
            $ensembl2mgi_hash{$ens_transcript}{"ens_protein"} = $ens_protein;
            $ensembl2mgi_hash{$ens_transcript}{$mgi_id}{"mgi_symbol"} = $mgi_symbol;
            $ensembl2mgi_hash{$ens_transcript}{$mgi_id}{"mgi_des"} = $mgi_des;

            $mgi_info_hash{$mgi_id}{"mgi_symbol"} = $mgi_symbol;
            $mgi_info_hash{$mgi_id}{"mgi_des"} = $mgi_des;
        } # if ($line) {
    }# while(@all){
} # if (open(ENS2MGI_FILE, $ensembl2mgi_file)){


### human entrezgene2ensembl ....
my $human_entrezgene2ensembl_file =  $data_dir . "ensembl/human_ensembl2entrez.txt";#ENTREZGENEID => ENSEMBL ID
close(HUMAN_ENTZ2ENSG);
my %human_entrezgene2ensg_hash = ();
if (open(HUMAN_ENTZ2ENSG, $human_entrezgene2ensembl_file)){
    @all = <HUMAN_ENTZ2ENSG>;
    close(HUMAN_ENTZ2ENSG);
    while(@all){
        $line = shift @all;
        chomp($line);
        if ($line) {
            my ($human_ensg, $human_entrezgene) = split("\t", $line);
            push @{$human_entrezgene2ensg_hash{$human_entrezgene}}, $human_ensg;
        } # if ($line){
    }# while(@all){
} # if (open(HUMAN_ENTZ2ENSG, $human_entrezgene2ensembl_file)){


### human2mouse orthologs  ....
my $human2mouse_orthologs_file = $data_dir . "ensembl/human2mouse_orthologs";
close(HUMAN2MOUSE_ORTHO);
my %human2mouse_ortholog_hash = ();
if (open(HUMAN2MOUSE_ORTHO, $human2mouse_orthologs_file)){
    @all = <HUMAN2MOUSE_ORTHO>;
    close(HUMAN2MOUSE_ORTHO);
    while(@all){
        $line = shift @all;
        chomp($line);
        if ($line) {
            my ($human_ensg, $mouse_ensg) = split("\t", $line);
            push @{$human2mouse_ortholog_hash{$human_ensg}}, $mouse_ensg;
        } # if ($line){
    }# while(@all){
} # if (open(HUMAN_ENTZ2ENSG, $human2mouse_orthologs_file)){


### mouse refseq_dna to ensembl_gene  ....
my $mouse_refseq_dna_file = $data_dir . "ensembl/mouse_ensg2refseq_combined.txt";
close(MOUSE_REFSEQ_FILE);
my %mouse_ensg2refseq_hash = ();sdrr
if (open(MOUSE_REFSEQ_FILE, $mouse_refseq_dna_file)){
    @all = <MOUSE_REFSEQ_FILE>;
    close(MOUSE_REFSEQ_FILE);
    while(@all){
        $line = shift @all;
        chomp($line);
        if ($line) {
            my ($mouse_ensg, $mouse_refseq) = split("\t", $line);
            push @{$mouse_ensg2refseq_hash{$mouse_refseq}}, $mouse_ensg;
        } # if ($line){
    }# while(@all){
} # if (open(MOUSE_REFSEQ_FILE, $mouse_refseq_dna_file)){
=end

=cut

##################################################################################################################################################
my %mirAlias2Accession_hash=%{retrieve($data_dir.'miRBase/mirAlias2Accession.hash')};

my $ensemblTranscriptID = $data_dir . "annotations/EnsmblTranscript_id.txt";
close(TRAIL_ENSEMBLTID_FILE);
my %ens_transcript2uni_hash = ();
if (open(TRAIL_ENSEMBLTID_FILE, $ensemblTranscriptID)){
    @all = <TRAIL_ENSEMBLTID_FILE>;
    close(TRAIL_ENSEMBLTID_FILE);
}
    while(@all){
        #$line = shift @all;
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
=begin
##############################################################
my(%human_microcosm_hash,%human_microcosm_detail_hash)=&microcosm('human_microcosm_mirna2target.txt');
my(%mouse_microcosm_hash,%mouse_microcosm_detail_hash)=&microcosm('mouse_microcosm_mirna2target.txt');
### mus_musculus Microcosm hash part ==microcosm==#########################
#Microcosm_V5 no updating is using miRanda algorithem to identify potential BS find highly complementary at 5 end(0-100)->Vienna routine for thermodynamic stability->check ortholog conservation.
##%this part need to be packed and reused for human, the place format are the same date.dir/human_microcosm_mirna2target.txt
sub microcosm(){ #input the name of file to be parsed
    my ($microcosm_mirna2target_file) =shift(@_);
    $microcosm_mirna2target_file=$data_dir."microcosm/".$microcosm_mirna2target_file;#
    close(MICROCOSM_FILE);
    my %microcosm_hash = ();
    my %microcosm_detail_hash = ();
    if (open(MICROCOSM_FILE, $microcosm_mirna2target_file)){
        @all = <MICROCOSM_FILE>;
        close(MICROCOSM_FILE);
        $line = shift @all; # header info deleted
    }# if (open(MICROCOSM_FILE, $mouse_microcosm_mirna2target_file)){
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
=end
=cut
###END==microcosm==

=begin
### mus_musculus mouse_targetscan_mirna2human_entrezgene hash part
#==targetscan==
##%it is currently using the Conserved_Site_Context_Scores.txt possible usage summary which taking into consider all 7mer and 8mer
#
my $mouse_targetscan_mirna2human_file = $data_dir . "targetscan/mouse_targetscan_mirna2human_entrezgene.txt";
close(TARGETSCAN_FILE);
my %targetscan_hash = ();
if (open(TARGETSCAN_FILE, $mouse_targetscan_mirna2human_file)){
    @all = <TARGETSCAN_FILE>;
    close(TARGETSCAN_FILE);
    $line = shift @all; # header info

    while(@all){
        $line = shift @all;
        chomp($line);
        if ($line) {
# Gene_ID Gene_Symbol Transcript_ID Gene_Tax_ID miRNA Site_Type UTR_start UTR_end 3prime_pairing  local_AU  position  TA  SPS context_score context_score_percentile
# Gene_ID  Gene_Symbol	Transcript-ID	Gene_Tax-ID	miRNA	Site_Type	UTR_start	UTR_end	3prime_pairing	local_AU	position	TA	SPS	context+ score	context+_score_percentile
#10090,10116,13616,8364,9031,9258,9544,9598,9606,9615,9796,9913->possible Gene Tax ID
#71667	0610007L01Rik	NM_001081394	8364	xtr-miR-9a	3	117	124	0.045	-0.006	-0.091	0.019 0.032	-0.248	88
#71667	0610007L01Rik	NM_001081394	8364	xtr-miR-9	3	117	124	0.045	-0.006	-0.091	0.019 0.032	-0.248	88
            my ($entrezgene_id, $gene_symbol, $refseq_id, $species_id, $mirna_id, $site_type, $utr_start, $utr_end, $three_prime_pairing, $local_au, $position, $TA, $SPS, $context_score, $context_percentile) = split("\t", $line);
            $mirna_id=~s/\*|star|_star//;
            if ($species_id == 10090) {
                my @temp_mgi_id = keys %{$mgi2entresgene_from_coordinate{$entrezgene_id}};
                my @mgi_id = ();
                if ($refseq2mgi_mrk_sequence{$refseq_id}) {
                    @mgi_id = @{$refseq2mgi_mrk_sequence{$refseq_id}};
                }
                push @mgi_id, @temp_mgi_id;

                if (my $tmp_mgi_id = $genename2mgiid_hash{$gene_symbol}) {
                    push @mgi_id, $tmp_mgi_id;
                }

                if (@mgi_id) {
                    @mgi_id = &nonRedundantList(\@mgi_id);	
                foreach my $mgi_id (@mgi_id) {
                    $all_resulted_mgi_ids_hash{$mgi_id}=1;
                    $targetscan_hash{$mirna_id}{$mgi_id} = 1;
                    $targetscan_hash{$mirna_id}{$mgi_id} = 1;
                } # foreach my $mgi_id (@mgi_id) {
            } # if (@mgi_id) {
            else {
                print "problem - targetscan - no mgi - line : $line\n";
            }
        } # if ($species_id == 10090) {
    } # if ($line) {
}# while(@all){
} # if (open(TARGETSCAN_FILE, $mouse_targetscan_mirna2human_file)){

### pictar part
# pictar results with 2 stringencies ... one with conservation between 7 species and other one with 13 species ...

## 7 species conservation ...==pictar==
my $pictar_7spec_file = $data_dir . "pictar/mouse_pictar_7spec_conservation.txt";

close(PT7);
my %pictar7_hash = ();
if (open(PT7, $pictar_7spec_file)){
    @all = <PT7>;
    close(PT7);

    while(@all){
        $line = shift @all;
        chomp($line);
        if ($line) {
            my ($refseq_id, $mirna_id) = split("_", $line);
            $mirna_id=~s/\*|star|_star//;
            my @mgi_id = ();
            if ($refseq2mgi_mrk_sequence{$refseq_id}) {
                @mgi_id = @{$refseq2mgi_mrk_sequence{$refseq_id}};
            }
            if (@mgi_id) {
                @mgi_id = &nonRedundantList(\@mgi_id);	
            foreach my $mgi_id (@mgi_id) {
                $all_resulted_mgi_ids_hash{$mgi_id}=1;
                $pictar7_hash{$mirna_id}{$mgi_id} = 1;
                $pictar7_hash{$mgi_id}{$mirna_id} = 1;
            }
        }	
        else {
            print "problem - pictar7 - no mgi - line : $line\n";
        }
    }
}
    }	

## 13 species conservation ...==pictar==
    #this part source file will be in format as below
#track_name,gene symbol/NCBI_RefSeq,gene symbol/NCBI RefSeq_location,data_source,PICTAR_score,target_site,target_site_location,genomic_strand,top-percent_value, csv human_hg19
    #add a seperationfor the Refseq ID
#hsa-miR-9-5p,ONECUT2,chr18:55102916-55158530,PICTAR,35,hsa-miR-9-5p,chr18:55153447-55153453,+,N/A
#hsa-miR-9-5p,ONECUT2,chr18:55102916-55158530,PICTAR,35,hsa-miR-9-5p,chr18:55146406-55146412,+,N/A
    #target files:/work/projects/pdgfr_kit/miRNA/mirna/picTar/human_pictar_prediction_hg19.csv
    #/work/projects/pdgfr_kit/miRNA/mirna/picTar/mouse_pictar_prediction_mm9.csv    #
    my $pictar_13spec_file = $data_dir . "pictar/mouse_pictar_13spec_conservation.txt";

    close(PT13);
    my %pictar13_hash = ();
    if (open(PT13, $pictar_13spec_file)){
        @all = <PT13>;
        close(PT13);

        while(@all){
            $line = shift @all;
            chomp($line);
            if ($line) {
                my ($refseq_id, $mirna_id) = split("_", $line);
                $mirna_id=~s/\*|star|_star//;
                my @mgi_id = ();
                if($refseq2mgi_mrk_sequence{$refseq_id}) {
                    @mgi_id = @{$refseq2mgi_mrk_sequence{$refseq_id}};
                }
                if (@mgi_id) {
                    @mgi_id = &nonRedundantList(\@mgi_id);
                foreach my $mgi_id (@mgi_id) {
                    $all_resulted_mgi_ids_hash{$mgi_id}=1;
                    $pictar13_hash{$mirna_id}{$mgi_id} = 1;
                    $pictar13_hash{$mgi_id}{$mirna_id} = 1;
                }
            }
            else {
                print "problem - pictar13 - no mgi - line : $line\n";
            }
        }
    }
}
#==pictar==
### ==microRNA.org== there is no updte, however there is the human prediction version (miranda/svr)
my $microrna_org_file = $data_dir . "microrna_org/mouse_predictions_S_C_aug2010.txt";

close(FH);
my %microrna_org_hash = ();
if (open(FH, $microrna_org_file)){
    @all = <FH>;
    close(FH);

    while(@all){
        $line = shift @all;
        chomp($line);
        if ($line) {
#mirbase_acc  mirna_name  gene_id gene_symbol transcript_id ext_transcript_id mirna_alignment alignment gene_alignment  mirna_start mirna_end gene_start  gene_end  genome_coordinates  conservation  align_score seed_cat  energy  mirsvr_score
            my ($mirbase_acc, $mirna_id,  $entrezgene_id, $gene_symbol, $ucsc_transcript_id, $ext_transcript_id, $mirna_alignment, $alignment, $gene_alignment, $mirna_start, $mirna_end, $gene_start, $gene_end, $genome_coordinates, $conservation, $align_score, $seed_cat, $energy, $mirsvr_score) = split("\t", $line);
            $mirna_id=~s/\*|star|_star//;
            my @mgi_id = ();
            @mgi_id = keys %{$mgi2entresgene_from_coordinate{$entrezgene_id}};
            if (my $tmp_mgi_id = $genename2mgiid_hash{$gene_symbol}) {
                push @mgi_id, $tmp_mgi_id;
            }
            if (@mgi_id) {
                @mgi_id = &nonRedundantList(\@mgi_id);	
            foreach my $mgi_id (@mgi_id) {
                $all_resulted_mgi_ids_hash{$mgi_id}=1;
                $microrna_org_hash{$mirna_id}{$mgi_id} = 1;
                $microrna_org_hash{$mgi_id}{$mirna_id} = 1;
            }
        }
        else {
            print "problem - microrna_org_hash - no mgi - line : $line\n";
        }
    }
}
  }

#### mouse miRDB 
  my $mouse_mirdb_file = $data_dir . "mirdb/mouse_mirdb.txt";
  close(MIRDB);
  my %mirdb_hash = ();
  if (open(MIRDB, $mouse_mirdb_file)){
      @all = <MIRDB>;
      close(MIRDB);
#    $line = shift @all; # header info
#name geneName  refGene readNum targetScanSites picTarSites RNA22Sites  PITASites miRandaSites
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
              my ($mirna_id, $refseq_id, $score) = split("\t", $line);
              $mirna_id=~s/\*|star|_star//;
              my @mgi_id = ();
              if ($refseq2mgi_mrk_sequence{$refseq_id}) {
                  @mgi_id = @{$refseq2mgi_mrk_sequence{$refseq_id}};
              }
              if (@mgi_id) {
                  @mgi_id = &nonRedundantList(\@mgi_id);	
              foreach my $mgi_id (@mgi_id) {
                  $all_resulted_mgi_ids_hash{$mgi_id}=1;
                  $mirdb_hash{$mirna_id}{$mgi_id} = 1;
                  $mirdb_hash{$mgi_id}{$mirna_id} = 1;
              }
          }
          else {
              print "problem - mirdb - no mgi - line : $line\n";
          }
      } # if ($line) {
  }# while(@all){
  } # if (open(MIRDB, $mouse_mirdb_file)){

#### mouse miRGen 
### it has 3 files ....
  my %mirgen_hash = ();

  my $mouse_mirgen_refseq_file = $data_dir . "mirgen/Refseq_ID.table";
  close(MIRGEN);
  if (open(MIRGEN, $mouse_mirgen_refseq_file)){
      @all = <MIRGEN>;
      close(MIRGEN);
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
              my ($refseq_id, $mirna_ids_length) = split("\t", $line);
              my @mirna_ids_length = split(", ", $mirna_ids_length);
              my @mirnas = ();
              foreach my $mirna_id_length (@mirna_ids_length) {
                  my ($mirna_id, $length) = split(":", $mirna_id_length);
                  push @mirnas, $mirna_id;
              }

              my @mgi_id = ();
              if ($refseq2mgi_mrk_sequence{$refseq_id}) {
                  @mgi_id = @{$refseq2mgi_mrk_sequence{$refseq_id}};
              }
              if (@mgi_id) {
                  @mgi_id = &nonRedundantList(\@mgi_id);
              foreach my $mgi_id (@mgi_id) {
                  $all_resulted_mgi_ids_hash{$mgi_id}=1;
                  foreach my $mirna_id (@mirnas) {
                      $mirna_id=~s/\*|star|_star//;
                      $mirgen_hash{$mirna_id}{$mgi_id} = 1;
                      $mirgen_hash{$mgi_id}{$mirna_id} = 1;
                  }
              }
          }
          else {
              print "problem - mirgen refseq - no mgi - line : $line\n";
          }
      } # if ($line) {
  }# while(@all){
  } # if (open(MIRGEN, $mouse_mirgen_refseq_file)){

  my $mouse_mirgen_entrez_file = $data_dir . "mirgen/ENTREZ_ID.table";

  close(MIRGEN);
  if (open(MIRGEN, $mouse_mirgen_entrez_file)){
      @all = <MIRGEN>;
      close(MIRGEN);
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
              my ($entrezgene_id, $mirna_ids_length) = split("\t", $line);
              my @mirna_ids_length = split(", ", $mirna_ids_length);
              my @mirnas = ();
              foreach my $mirna_id_length (@mirna_ids_length) {
                  my ($mirna_id, $length) = split(":", $mirna_id_length);
                  push @mirnas, $mirna_id;
              }

              my @mgi_id = ();
              @mgi_id = keys %{$mgi2entresgene_from_coordinate{$entrezgene_id}};
              if (@mgi_id) {
                  foreach my $mgi_id (@mgi_id) {
                      $all_resulted_mgi_ids_hash{$mgi_id}=1;
                      foreach my $mirna_id (@mirnas) {
                          $mirna_id=~s/\*|star|_star//;
                          $mirgen_hash{$mirna_id}{$mgi_id} = 1;
                          $mirgen_hash{$mgi_id}{$mirna_id} = 1;
                      }
                  }
              }
              else {
                  print "problem - mirgen entrezgene - no mgi - line : $line\n";
              }
          } # if ($line) {
      }# while(@all){
  } # if (open(MIRGEN, $mouse_mirgen_entrez_file)){

  my $mouse_mirgen_ensg_file   = $data_dir . "mirgen/ENSG_ID.table";

  close(MIRGEN);
  if (open(MIRGEN, $mouse_mirgen_ensg_file)){
      @all = <MIRGEN>;
      close(MIRGEN);
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
              my ($ensembl_gene_id, $mirna_ids_length) = split("\t", $line);
              my @mirna_ids_length = split(", ", $mirna_ids_length);
              my @mirnas = ();
              foreach my $mirna_id_length (@mirna_ids_length) {
                  my ($mirna_id, $length) = split(":", $mirna_id_length);
                  push @mirnas, $mirna_id;
              }

              my @mgi_id = ();
              if ($ens_gene2mgi_hash{$ensembl_gene_id}) {
                  @mgi_id = @{$ens_gene2mgi_hash{$ensembl_gene_id}};
              }
              if (@mgi_id) {
                  foreach my $mgi_id (@mgi_id) {
                      $all_resulted_mgi_ids_hash{$mgi_id}=1;
                      foreach my $mirna_id (@mirnas) {
                          $mirna_id=~s/\*|star|_star//;
                          $mirgen_hash{$mirna_id}{$mgi_id} = 1;
                          $mirgen_hash{$mgi_id}{$mirna_id} = 1;
                      }
                  }
              }
              else {
                  print "problem - mirgen ensembl - no mgi - line : $line\n";
              }
          } # if ($line) {
      }# while(@all){
  } # if (open(MIRGEN, $mouse_mirgen_ensg_file)){


### miRNAMap ... part ....
  my $mirnamap_file   = $data_dir . "mirnamap/miRNA_targets_mmu.txt";
  my %mirnamap_hash = ();
  close(MIRNAMAP);
  if (open(MIRNAMAP, $mirnamap_file)){
      @all = <MIRNAMAP>;
      $line = shift @all; # header info
      close(MIRNAMAP);
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
#mature miRNA    Ensembl transcript ID   target start    tsrget end      miRNA 3-5       alinment      target 5-3       tool name       criterion 1     criterion 2     criterion 3
              my ($mirna_id, $ensembl_transcript_id, $target_start, $target_end, $mirna_3_5, $alinment, $target_5_3, $tool_name, $criterion_1, $criterion_2, $criterion_3) = split("\t", $line);

              $mirna_id=~s/\*|star|_star//;
              my @mgi_id = ();
              if ($ens_transcript2mgi_hash{$ensembl_transcript_id}) {
                  @mgi_id = @{$ens_transcript2mgi_hash{$ensembl_transcript_id}};
              }
              if (@mgi_id) {
                  foreach my $mgi_id (@mgi_id) {
                      $all_resulted_mgi_ids_hash{$mgi_id}=1;
                      $mirnamap_hash{$mirna_id}{$mgi_id} = 1;
                      $mirnamap_hash{$mgi_id}{$mirna_id} = 1;
                  } 
              }
              else {
                  print "problem - mirnamap - no mgi - line : $line\n";
              } 
          } # if ($line) { 
      }# while(@all){
  } # if (open(MIRNAMAP, $mirnamap_file)){


### miRTarBase ... part ....
  my $mirtarbase_file   = $data_dir . "mirtarbase/mmu_MTI.txt";
  my %mirtarbase_hash = ();
  close(MIRTARBASE);
  if (open(MIRTARBASE, $mirtarbase_file)){
      @all = <MIRTARBASE>;
      $line = shift @all; # header info
      close(MIRTARBASE);
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
#miRTarBase ID   miRNA   Species (miRNA) Target Gene     Target Gene (Entrez ID) Species (Target Gene) Experiments      Support Type    References (PMID)
              my ($mirtarbase_id, $mirna_id, $species_mirna, $gene_symbol, $entrezgene_id, $species_target, $experiments, $support_type, $references) = split("\t", $line);
              $mirna_id=~s/\*|star|_star//;
              if (($species_mirna eq "Mus musculus") && ($species_target eq "Mus musculus")) {
                  my @mgi_id = ();
                  @mgi_id = keys %{$mgi2entresgene_from_coordinate{$entrezgene_id}};
                  if (my $tmp_mgi_id = $genename2mgiid_hash{$gene_symbol}) {
                      push @mgi_id, $tmp_mgi_id;
                  }
                  if (@mgi_id) {
                      @mgi_id = &nonRedundantList(\@mgi_id);	
                  foreach my $mgi_id (@mgi_id) {
                      $all_resulted_mgi_ids_hash{$mgi_id}=1;
                      $mirtarbase_hash{$mirna_id}{$mgi_id} = 1;
                      $mirtarbase_hash{$mgi_id}{$mirna_id} = 1;
                  }
              }
              else {
                  print "problem - mirtarbase mmu_mti - no mgi - line : $line\n";
              }
          } # if (($species_mirna eq "Mus musculus") && ($species_target eq "Mus musculus")) {
      } # if ($line) { 
  }# while(@all){
  } # if (open(MIRTARBASE, $mirtarbase_file)){

### miRTarBase ... part ...Supported by strong experimental evidences (Reporter assay or Western blot)
  my $mirtarbase_wr_file   = $data_dir . "mirtarbase/miRTarBase_SE_WR.txt";
  my %mirtarbase_wr_hash = ();
  close(MIRTARBASE);
  if (open(MIRTARBASE, $mirtarbase_wr_file)){
      @all = <MIRTARBASE>;
      $line = shift @all; # header info
      close(MIRTARBASE);
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
#miRTarBase ID   miRNA   Species (miRNA) Target Gene     Target Gene (Entrez ID) Species (Target Gene) Experiments      Support Type    References (PMID)
              my ($mirtarbase_id, $mirna_id, $species_mirna, $gene_symbol, $entrezgene_id, $species_target, $experiments, $support_type, $references) = split("\t", $line);
              $mirna_id=~s/\*|star|_star//;
              if (($species_mirna eq "Mus musculus") && ($species_target eq "Mus musculus")) {
                  my @mgi_id = ();
                  @mgi_id = keys %{$mgi2entresgene_from_coordinate{$entrezgene_id}};
                  if (my $tmp_mgi_id = $genename2mgiid_hash{$gene_symbol}) {
                      push @mgi_id, $tmp_mgi_id;
                  }
                  if (@mgi_id) {
                      @mgi_id = &nonRedundantList(\@mgi_id);	
                  foreach my $mgi_id (@mgi_id) {
                      $all_resulted_mgi_ids_hash{$mgi_id}=1;
                      $mirtarbase_wr_hash{$mirna_id}{$mgi_id} = 1;
                      $mirtarbase_wr_hash{$mgi_id}{$mirna_id} = 1;
                  }
              }
              else {
                  print "problem - mirtarbase wr - no mgi - line : $line\n";
              }
          } # if (($species_mirna eq "Mus musculus") && ($species_target eq "Mus musculus")) {
      } # if ($line) {
  }# while(@all){
  } # if (open(MIRTARBASE, $mirtarbase_wr_file)){

### miRTarBase ... part ...Supported by strong experimental evidences (Western blot)
my $mirtarbase_w_file   = $data_dir . "mirtarbase/miRTarBase_SE_W.txt";
  my %mirtarbase_w_hash = ();
  close(MIRTARBASE);
  if (open(MIRTARBASE, $mirtarbase_w_file)){
      @all = <MIRTARBASE>;
      $line = shift @all; # header info
      close(MIRTARBASE);
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
#miRTarBase ID   miRNA   Species (miRNA) Target Gene     Target Gene (Entrez ID) Species (Target Gene) Experiments      Support Type    References (PMID)
              my ($mirtarbase_id, $mirna_id, $species_mirna, $gene_symbol, $entrezgene_id, $species_target, $experiments, $support_type, $references) = split("\t", $line);
              $mirna_id=~s/\*|star|_star//;
              if (($species_mirna eq "Mus musculus") && ($species_target eq "Mus musculus")) {
                  my @mgi_id = ();
                  @mgi_id = keys %{$mgi2entresgene_from_coordinate{$entrezgene_id}};
                  if (my $tmp_mgi_id = $genename2mgiid_hash{$gene_symbol}) {
                      push @mgi_id, $tmp_mgi_id;
                  }
                  if (@mgi_id) {
                      @mgi_id = &nonRedundantList(\@mgi_id);	
                  foreach my $mgi_id (@mgi_id) {
                      $all_resulted_mgi_ids_hash{$mgi_id}=1;
                      $mirtarbase_w_hash{$mirna_id}{$mgi_id} = 1;
                      $mirtarbase_w_hash{$mgi_id}{$mirna_id} = 1;
                  }
              }
              else {
                  print "problem - mirtarbase w - no mgi - line : $line\n";
              }
          } # if (($species_mirna eq "Mus musculus") && ($species_target eq "Mus musculus")) {
      } # if ($line) {
  }# while(@all){
  } # if (open(MIRTARBASE, $mirtarbase_w_file)){

### miRTarBase ... part ...Supported by strong experimental evidences (Reporter assay)
  my $mirtarbase_r_file   = $data_dir . "mirtarbase/miRTarBase_SE_R.txt";
  my %mirtarbase_r_hash = ();
  close(MIRTARBASE);
  if (open(MIRTARBASE, $mirtarbase_r_file)){
      @all = <MIRTARBASE>;
      $line = shift @all; # header info
      close(MIRTARBASE);
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
#miRTarBase ID   miRNA   Species (miRNA) Target Gene     Target Gene (Entrez ID) Species (Target Gene) Experiments      Support Type    References (PMID)
              my ($mirtarbase_id, $mirna_id, $species_mirna, $gene_symbol, $entrezgene_id, $species_target, $experiments, $support_type, $references) = split("\t", $line);
              $mirna_id=~s/\*|star|_star//;
              if (($species_mirna eq "Mus musculus") && ($species_target eq "Mus musculus")) {
                  my @mgi_id = ();
                  @mgi_id = keys %{$mgi2entresgene_from_coordinate{$entrezgene_id}};
                  if (my $tmp_mgi_id = $genename2mgiid_hash{$gene_symbol}) {
                      push @mgi_id, $tmp_mgi_id;
                  }
                  if (@mgi_id) {
                      @mgi_id = &nonRedundantList(\@mgi_id);	
                  foreach my $mgi_id (@mgi_id) {
                      $all_resulted_mgi_ids_hash{$mgi_id}=1;
                      $mirtarbase_r_hash{$mirna_id}{$mgi_id} = 1;
                      $mirtarbase_r_hash{$mgi_id}{$mirna_id} = 1;
                  }
              }
              else {
                  print "problem - mirtarbase r - no mgi - line : $line\n";
              }
          } # if (($species_mirna eq "Mus musculus") && ($species_target eq "Mus musculus")) {
      } # if ($line) {
  }# while(@all){
  } # if (open(MIRTARBASE, $mirtarbase_r_file)){

### miRTarBase ... part ...Supported by weak experimental evidences (Microarray or pSILAC)
  my $mirtarbase_mp_file   = $data_dir . "mirtarbase/miRTarBase_WE_MP.txt";
  my %mirtarbase_mp_hash = ();
  close(MIRTARBASE);
  if (open(MIRTARBASE, $mirtarbase_mp_file)){
      @all = <MIRTARBASE>;
      $line = shift @all; # header info
      close(MIRTARBASE);
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
#miRTarBase ID   miRNA   Species (miRNA) Target Gene     Target Gene (Entrez ID) Species (Target Gene) Experiments      Support Type    References (PMID)
              my ($mirtarbase_id, $mirna_id, $species_mirna, $gene_symbol, $entrezgene_id, $species_target, $experiments, $support_type, $references) = split("\t", $line);
              $mirna_id=~s/\*|star|_star//;
              if (($species_mirna eq "Mus musculus") && ($species_target eq "Mus musculus")) {
                  my @mgi_id = ();
                  @mgi_id = keys %{$mgi2entresgene_from_coordinate{$entrezgene_id}};
                  if (my $tmp_mgi_id = $genename2mgiid_hash{$gene_symbol}) {
                      push @mgi_id, $tmp_mgi_id;
                  }
                  if (@mgi_id) {
                      @mgi_id = &nonRedundantList(\@mgi_id);	
                  foreach my $mgi_id (@mgi_id) {
                      $all_resulted_mgi_ids_hash{$mgi_id}=1;
                      $mirtarbase_mp_hash{$mirna_id}{$mgi_id} = 1;
                      $mirtarbase_mp_hash{$mgi_id}{$mirna_id} = 1;
                  }
              }
              else {
                  print "problem - mirtarbase mp - no mgi - line : $line\n";
              }
          } # if (($species_mirna eq "Mus musculus") && ($species_target eq "Mus musculus")) {
      } # if ($line) {
  }# while(@all){
  } # if (open(MIRTARBASE, $mirtarbase_mp_file)){

### PITA ... part ...
  my $pita_top_file   = $data_dir . "pita/PITA_targets_mm9_TOP.tab";
  my %pita_top_hash = ();
  close(PITA);
  if (open(PITA, $pita_top_file)){
      @all = <PITA>;
      $line = shift @all; # header info
      close(PITA);
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
#RefSeq	Name	microRNA
              my ($refseq_ids, $gene_symbol, $mirna_id) = split("\t", $line);
              #$mirna_id=~s/\*|star|_star//;
              my @refseq_ids = split(";", $refseq_ids);
              my @mgi_id = ();
              foreach my $refseq_id (@refseq_ids) {
                  my @temp_mgi_id = ();
                  if ($refseq2mgi_mrk_sequence{$refseq_id}) {
                      @temp_mgi_id = @{$refseq2mgi_mrk_sequence{$refseq_id}};
                      push @mgi_id, @temp_mgi_id;
                  }
              }
              if (my $tmp_mgi_id = $genename2mgiid_hash{$gene_symbol}) {
                  push @mgi_id, $tmp_mgi_id;
              }

              if (@mgi_id) {
                  @mgi_id = &nonRedundantList(\@mgi_id);
              foreach my $mgi_id (@mgi_id) {
                  $all_resulted_mgi_ids_hash{$mgi_id}=1;
                  $pita_top_hash{$mirna_id}{$mgi_id} = 1;
                  $pita_top_hash{$mgi_id}{$mirna_id} = 1;
              }
          }
          else {
              print "problem - pita top - no mgi - line : $line\n";
          }
      } # if ($line) {
  }# while(@all){
  } # if (open(PITA, $pita_top_file)){

### PITA ... part ...ALL
  my $pita_all_file   = $data_dir . "pita/PITA_targets_mm9_ALL.tab";
  my %pita_all_hash = ();
  close(PITA);
  if (open(PITA, $pita_all_file)){
      @all = <PITA>;
      $line = shift @all; # header info
      close(PITA);
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
#RefSeq	Name	microRNA
              my ($refseq_ids, $gene_symbol, $mirna_id) = split("\t", $line);
              $mirna_id=~s/\*|star|_star//;
              my @refseq_ids = split(";", $refseq_ids);
              my @mgi_id = ();
              foreach my $refseq_id (@refseq_ids) {
                  my @temp_mgi_id = ();
                  if ($refseq2mgi_mrk_sequence{$refseq_id}) {
                      @temp_mgi_id = @{$refseq2mgi_mrk_sequence{$refseq_id}};
                      push @mgi_id, @temp_mgi_id;
                  }
              }
              if (my $tmp_mgi_id = $genename2mgiid_hash{$gene_symbol}) {
                  push @mgi_id, $tmp_mgi_id;
              }

              if (@mgi_id) {
                  @mgi_id = &nonRedundantList(\@mgi_id);
              foreach my $mgi_id (@mgi_id) {
                  $all_resulted_mgi_ids_hash{$mgi_id}=1;
                  $pita_all_hash{$mirna_id}{$mgi_id} = 1;
                  $pita_all_hash{$mgi_id}{$mirna_id} = 1;
              }
          }
          else {
              print "problem - pita all - no mgi - line : $line\n";
          }
      } # if ($line) {
  }# while(@all){
  } # if (open(PITA, $pita_all_file)){


### RepTar ... part ...
  my $reptar_file   = $data_dir . "reptar/mouse-predictions.txt";
  my %reptar_hash = ();
  close(REPTAR);
  if (open(REPTAR, $reptar_file)){
      @all = <REPTAR>;
      close(REPTAR);
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
              my ($gene_name_refseq_ids, $mirna_id, $start, $end, $mfe, $nm, $gu, $prof, $pic, $b_cons, $t_cons, $rep, $crep) = split("\t", $line);
              $mirna_id=~s/\*|star|_star//;
              my ($gene_symbol, $refseq_id) = split(":::", $gene_name_refseq_ids);
              my @mgi_id = ();
              if ($refseq2mgi_mrk_sequence{$refseq_id}) {
                  @mgi_id = @{$refseq2mgi_mrk_sequence{$refseq_id}};
              }
              if (my $tmp_mgi_id = $genename2mgiid_hash{$gene_symbol}) {
                  push @mgi_id, $tmp_mgi_id;
              }
              if (@mgi_id) {
                  @mgi_id = &nonRedundantList(\@mgi_id);
              foreach my $mgi_id (@mgi_id) {
                  $all_resulted_mgi_ids_hash{$mgi_id}=1;
                  $reptar_hash{$mirna_id}{$mgi_id} = 1;
                  $reptar_hash{$mgi_id}{$mirna_id} = 1;
              }
          }
          else {
              print "problem - reptar - no mgi - line : $line\n";
          }
      } # if ($line) {
  }# while(@all){
  } # if (open(REPTAR, $reptar_file)){

#### mouse starbase
  my $mouse_starbase_file = $data_dir . "starbase/starBase_Mouse_miRNA-target_interactions2012-05-07_22-40.csv";
  close(STARBASE);
  my %starbase_hash = ();
  my %starbase_detail_hash = ();
  if (open(STARBASE, $mouse_starbase_file)){
      @all = <STARBASE>;
      close(STARBASE);
      $line = shift @all; # header info
#name geneName  refGene readNum targetScanSites picTarSites RNA22Sites  PITASites miRandaSites
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
              my ($mirna_id, $gene_symbol, $refseq_id, $read_num, $targetScanSites, $picTarSites, $RNA22Sites, $PITASites, $miRandaSites) = split(",", $line);
              $mirna_id=~s/\*|star|_star//;

              my @mgi_id = ();
              if ($refseq2mgi_mrk_sequence{$refseq_id}) {
                  @mgi_id = @{$refseq2mgi_mrk_sequence{$refseq_id}};
              }
              if (my $tmp_mgi_id = $genename2mgiid_hash{$gene_symbol}) {
                  push @mgi_id, $tmp_mgi_id;
              }
              if (@mgi_id) {
                  @mgi_id = &nonRedundantList(\@mgi_id);
              foreach my $mgi_id (@mgi_id) {
                  $all_resulted_mgi_ids_hash{$mgi_id}=1;
                  $starbase_hash{$mirna_id}{$mgi_id} = 1;
                  $starbase_hash{$mgi_id}{$mirna_id} = 1;

                  $starbase_detail_hash{$mirna_id}{$mgi_id}{"gene_name"} = $gene_symbol;
                  $starbase_detail_hash{$mirna_id}{$mgi_id}{"targetScanSites"} = $targetScanSites;
                  $starbase_detail_hash{$mirna_id}{$mgi_id}{"picTarSites"} = $picTarSites;
                  $starbase_detail_hash{$mirna_id}{$mgi_id}{"RNA22Sites"} = $RNA22Sites;
                  $starbase_detail_hash{$mirna_id}{$mgi_id}{"PITASites"} = $PITASites;
                  $starbase_detail_hash{$mirna_id}{$mgi_id}{"miRandaSites"} = $miRandaSites;
              }
          }
          else {
              print "problem - starbase - no mgi - line : $line\n";
          }
      } # if ($line) {
  }# while(@all){
  } # if (open(STARBASE, $mouse_starbase_file)){

#### mouse tarbase
  my $mouse_tarbase_file = $data_dir . "tarbase/TarBase_V5.0.txt";
  close(TARBASE);
  my %tarbase_hash = ();
  my %tarbase_detail_hash = ();
  if (open(TARBASE, $mouse_tarbase_file)){
      @all = <TARBASE>;
      close(TARBASE);
      $line = shift @all; # header info
#Id      Id_V4   Data_Type       Support_Type    Organism        miRNA   HGNC_Symbol     Gene    Isoform Ensembl Chr_loc     MRE     S_S_S   I_S     D_S     Validation      Paper   Target_seq      miRNA_seq       Seq_location    PMID	KEGG    Protein_type    Differentially_expressed_in     Pathology_or_Event      Mis_Regulation  Gene_Expression Tumour_Involvement  Bibliographic_Notes     Cell_Line_Used  HGNC_ID SwissProt
      while(@all){
          $line = shift @all;
          chomp($line);
          if ($line) {
              my ($id, $id_v4, $data_type, $support_type, $org, $mirna, $hgnc_symbol, $gene_symbol, $isoform, $ensembl_gene_id, $chr_loc, $mre, $s_s_s, $i_s, $d_s, $validation, $paper, $target_seq, $mirna_seq, $seq_location, $pmid, $kegg, $protein_type, $differentially_expressed_in, $pathology_or_event, $mis_regulation, $gene_expression, $tumor_involvement, $bibliographic_notes, $cell_line_used, $hgnc_id, $swissprot) = split(",", $line);

              if ($org eq "Mouse") {
                  my $mirna_id = "mmu-" . $mirna;
                  $mirna_id=~s/\*|star|_star//;
                  my @mgi_id = ();
                  if ($ens_gene2mgi_hash{$ensembl_gene_id}) {
                      @mgi_id = @{$ens_gene2mgi_hash{$ensembl_gene_id}};
                  }
                  if (my $tmp_mgi_id = $genename2mgiid_hash{$gene_symbol}) {
                      push @mgi_id, $tmp_mgi_id;
                  }
                  if (@mgi_id) {
                      @mgi_id = &nonRedundantList(\@mgi_id);
                  foreach my $mgi_id (@mgi_id) {
                      $all_resulted_mgi_ids_hash{$mgi_id}=1;
                      $tarbase_hash{$mirna_id}{$mgi_id} = 1;
                      $tarbase_hash{$mgi_id}{$mirna_id} = 1;

                      $tarbase_detail_hash{$mirna_id}{$mgi_id}{"data_type"} = $data_type;
                      $tarbase_detail_hash{$mirna_id}{$mgi_id}{"gene_name"} = $gene_symbol;
                      $tarbase_detail_hash{$mirna_id}{$mgi_id}{"paper"} = $paper;
                      $tarbase_detail_hash{$mirna_id}{$mgi_id}{"pmid"} = $pmid;
                      $tarbase_detail_hash{$mirna_id}{$mgi_id}{"pathology_or_event"} = $pathology_or_event;
                      $tarbase_detail_hash{$mirna_id}{$mgi_id}{"mis_regulation"} = $mis_regulation;
                      $tarbase_detail_hash{$mirna_id}{$mgi_id}{"bibliographic_notes"} = $bibliographic_notes;
                      $tarbase_detail_hash{$mirna_id}{$mgi_id}{"cell_line_used"} = $cell_line_used;
                  }
              }
              else {
                  print "problem - tarbase - no mgi - line : $line\n";
              }
          } # if ($org eq "Mouse") {
      } # if ($line) {
  }# while(@all){
  } # if (open(TARBASE, $mouse_tarbase_file)){


### process master file

#	print MICROCOSM_RESULT_FILE "pain_related_gene\tpain_related_mRNA\tmirbase_id\tscore\tpvalue_og\tmgi_id\tmgi_symbol\tmgi_des\n";				
#	print TARGETSCAN_RESULT_FILE "pain_related_gene\tpain_related_mRNA\tmirbase_id\tcontext_score\tcontext_percentile\tmgi_id\tmgi_symbol\tmgi_des\n";				
#	print PICTAR_RESULT_FILE "pain_related_gene\tpain_related_mRNA\tmirbase_id\tscore\trank\tmgi_id\tmgi_symbol\tmgi_des\n";				
  my %microcosm_result_hash = ();
  my %targetscan_result_hash = ();
  my %pictar_result_hash = ();
  my @mirbase_ids = ();
  my %mirna_info_hash = ();
  my @all_resulted_mgi_ids = ();
  my $master_file = $data_dir . "kiran_exp/mirbase_ids__fold_change.txt";
#print "master_file : $master_file\n";	
  close(MASTER_FILE);
  if (open(MASTER_FILE, $master_file)){
      @all = <MASTER_FILE>;
      close(MASTER_FILE);
      my $i = 0;
      while(@all){
          $i++;
          $line = shift @all;
          chomp($line);
          my ($mirbase_id, $mirna_fold_change) = split("\t", $line);
          $mirbase_id=~s/\*//;
#print "----------------i:$i   mirbase_id : $mirbase_id\n";			
          push @mirbase_ids, $mirbase_id;
          my $mirna_regulation = "miRNA up";
          if ($mirna_fold_change =~ '-') {
              $mirna_regulation = "miRNA down";
          }
          $mirna_info_hash{$mirbase_id}{'mirna_fold_change'} = $mirna_fold_change;
          $mirna_info_hash{$mirbase_id}{'mirna_regulation'} = $mirna_regulation;
      }# while(@all){

#### comparative analysis ............
      print RESULT_FILE "priority\tpain_related\tmicrocosm\ttargetscan\tpictar7\tpictar13\tmicrorna_org\tmirdb\tmirgen\tmirnamap\tmirtarbase\tmirtarbase_wr(V)\tpita_top\tpita_all\treptar\tstarbase\ttarbase(V)\tmirbase_id\tmgi_id\tmgi_symbol\tmgi_des\n";				
      @all_resulted_mgi_ids = keys %all_resulted_mgi_ids_hash;
      @mirbase_ids = &nonRedundantList(\@mirbase_ids);

  foreach my $mirbase_id (@mirbase_ids) {
      chomp($mirbase_id);
      my $mirna_fold_change = $mirna_info_hash{$mirbase_id}{'mirna_fold_change'};
      my $mirna_regulation = $mirna_info_hash{$mirbase_id}{'mirna_regulation'};
      foreach my $mgi_id (@all_resulted_mgi_ids) {
#print "mirbase_id : $mirbase_id, mgi_id : $mgi_id\n";			
          my $pain_related = "";
          my $pain_related_mrna = "";
          my $pain_related_mrna_des = "";
          my $pain_related_mrna_symbol = "";
          my $pain_related_mrna_fold_change = "";
          my $pain_related_mrna_p_value = "";
          my $mgi_symbol   = "";
          my $mgi_des      = "";
          my $priority = 0;
          my $at_least_one_target = 0;
          if ($pain_genes_hash{$mgi_id}) {
              $pain_related = "pain_yes";
              $priority++;
          }

          if ($pain_mrna_hash{$mgi_id}) {
              $pain_related_mrna = "yes";
          }

          $mgi_symbol = $mgi_coordinate_hash{$mgi_id}{"gene_symbol"};
          $mgi_des = $mgi_coordinate_hash{$mgi_id}{"gene_name"};

### microcosm (prediction)
          my $microcosm;
          if ($microcosm_hash{$mirbase_id}{$mgi_id}) {
              $microcosm = "yes";
              $at_least_one_target = 1;
              $priority++;
          }

### targetscan (prediction)
          my $targetscan;
          if ($targetscan_hash{$mirbase_id}{$mgi_id}) {
              $targetscan = "yes";
              $at_least_one_target = 1;
              $priority++;
          }

### pictar 7 species part (prediction)
          my $pictar7;
          if ($pictar7_hash{$mirbase_id}{$mgi_id}) {
              $pictar7 = "yes";
              $at_least_one_target = 1;
              $priority++;
          }

### pictar 13 species part (prediction)
          my $pictar13;
          if ($pictar13_hash{$mirbase_id}{$mgi_id}) {
              $pictar13 = "yes";
              $at_least_one_target = 1;
              $priority++;
          }

###  microrna_org (prediction)
          my $microrna_org;
          if ($microrna_org_hash{$mirbase_id}{$mgi_id}) {
              $microrna_org = "yes";
              $at_least_one_target = 1;
              $priority++;
          }

###  mirdb (prediction)
          my $mirdb;
          if ($mirdb_hash{$mirbase_id}{$mgi_id}) {
              $mirdb = "yes";
              $at_least_one_target = 1;
              $priority++;
          }

###  mirgen (prediction)
          my $mirgen;
          if ($mirgen_hash{$mirbase_id}{$mgi_id}) {
              $mirgen = "yes";
              $at_least_one_target = 1;
              $priority++;
          }

###  mirnamap (prediction)
          my $mirnamap;
          if ($mirnamap_hash{$mirbase_id}{$mgi_id}) {
              $mirnamap = "yes";
              $at_least_one_target = 1;
              $priority++;
          }

###  mirtarbase (prediction)
          my $mirtarbase;
          if ($mirtarbase_hash{$mirbase_id}{$mgi_id}) {
              $mirtarbase = "yes";
              $at_least_one_target = 1;
              $priority++;
          }

###  mirtarbase (Supported by strong experimental evidences (Reporter assay or Western blot))
          my $mirtarbase_wr;
          if ($mirtarbase_wr_hash{$mirbase_id}{$mgi_id}) {
              $mirtarbase_wr = "yes";
              $at_least_one_target = 1;
              $priority = $priority+3;
          }
          elsif ($mirtarbase_r_hash{$mirbase_id}{$mgi_id}) {
              $mirtarbase_wr = "yes";
              $at_least_one_target = 1;
              $priority = $priority+3;
          }
          elsif ($mirtarbase_w_hash{$mirbase_id}{$mgi_id}) {
              $mirtarbase_wr = "yes";
              $at_least_one_target = 1;
              $priority = $priority+3;
          }
          elsif ($mirtarbase_mp_hash{$mirbase_id}{$mgi_id}) {
              $mirtarbase_wr = "yes";
              $at_least_one_target = 1;
              $priority = $priority+2;
          }

###  pita top(prediction)
          my $pita_top;
          if ($pita_top_hash{$mirbase_id}{$mgi_id}) {
              $pita_top = "yes";
              $at_least_one_target = 1;
              $priority++;
          }

###  pita all(prediction)
          my $pita_all;
          if ($pita_all_hash{$mirbase_id}{$mgi_id}) {
              $pita_all = "yes";
              $at_least_one_target = 1;
              $priority++;
          }

###  reptar(prediction)
          my $reptar;
          if ($reptar_hash{$mirbase_id}{$mgi_id}) {
              $at_least_one_target = 1;
              $reptar = "yes";
              $priority++;
          }

###  starbase(prediction)
          my $starbase;
          if ($starbase_hash{$mirbase_id}{$mgi_id}) {
              $starbase = "yes";
              $at_least_one_target = 1;
              $priority++;
          }

###  tarbase(obtained from literature with experimental evidance)
          my $tarbase;
          if ($tarbase_hash{$mirbase_id}{$mgi_id}) {
              $tarbase = "yes";
              $at_least_one_target = 1;
              $priority = $priority+3;
          }



          # print only pain related mRNAs
#				if ($pain_related_mrna eq "yes") {
#print "srk srk srk microcosm : $microcosm, targetscan : $targetscan, pictar : $pictar\n";				
          if ($at_least_one_target) {
#print "--Radhe-Shyam Radhe-Shyam Radhe-Shyam\n";					
              print RESULT_FILE "$priority\t$pain_related\t$microcosm\t$targetscan\t$pictar7\t$pictar13\t$microrna_org\t$mirdb\t$mirgen\t$mirnamap\t$mirtarbase\t$mirtarbase_wr\t$pita_top\t$pita_all\t$reptar\t$starbase\t$tarbase\t$mirbase_id\t$mgi_id\t$mgi_symbol\t$mgi_des\n";				
          } # if ($at_least_one_target) {
#				} # if ($pain_related_mrna eq "yes") {

      } # foreach my $mgi_id (@all_resulted_mgi_ids) {
  }	# foreach my $mirbase_id (@mirbase_ids) {

    } # if (open(MASTER_FILE, $master_file)){
    close(RESULT_FILE);
#	close(MICROCOSM_RESULT_FILE);
#	close(TARGETSCAN_RESULT_FILE);
#	close(PICTAR_RESULT_FILE);
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
