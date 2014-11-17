#!/usr/bin/perl -w
use strict;
use warnings;

use Data::Dumper qw(Dumper);
use Bio::EnsEMBL::Registry;

my $regsitry='Bio::EnsEMBL::Registry';
$regsitry->load_registry_from_db(-HOST => 'ensembldb.ensembl.org',-PORT => 5306,-USER => 'anonymous',-VERBOSE => 0);
my $archiveStableIdAdaptor = $regsitry->get_adaptor( 'Human', 'Core', 'ArchiveStableId' );

my $oldensemblfile='/home/jiali/workspace/projects/vs/9606_mirnamap_221missingEnsemblID.txt.uniq';
if (open (MIRALIAS,$oldensemblfile)){
    my @all=<MIRALIAS>;
    close(MIRALIAS);
    while(@all){
        my $original_id = shift @all;
        chomp ($original_id);#every line contains a id
        my $archive_id =$archiveStableIdAdaptor->fetch_by_stable_id($original_id);       
        
        #print Dumper $archive_id;
        my $successor = $archive_id->get_latest_incarnation();
        if($original_id eq $successor->stable_id){print "Original ID:".$original_id.", Latest version of ".$archive_id->stable_id." is ".$successor->stable_id."\n";}


        # Then get the associated archive info which tells us what Transcripts it had at retirement 
#        my $archived_info = $archive_id->get_all_associated_archived();
#        foreach my $associated_ref(@{$archived_info}){
#           print Dumper $associated_ref;
#           print "--------------------------------------------------\n";    #    }
        # print "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    }
}



=begin  BlockComment  # BlockCommentNo_1

foreach my $associated (@{$archived_info}) {
  my ($arch_gene, $arch_tr, $arch_tl, $pep_seq) = @{$associated};
  my $successor = $arch_tr->get_latest_incarnation();
  if($successor->release() == $release) {
    my $live_transcript = $ta->fetch_by_stable_id($successor->stable_id());
    my $gene = $live_transcript->get_Gene();
    $possible_new_genes{$gene->stable_id()} = $gene;
  }
}

foreach my $id (keys %possible_new_genes) {
  printf("%d | %s -> %s\n", $release, $original_id, $id);
}

=end    BlockComment  # BlockCommentNo_1

=cut

