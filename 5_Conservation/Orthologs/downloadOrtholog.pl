#!/usr/bin/perl

# An example script demonstrating the use of BioMart API.
# This perl API representation is only available for configuration versions >=  0.5 
use strict; use warnings;
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;
use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input transcript ID bedfiles (e.g. hg19_gencode.bed)>\n" unless @ARGV == 1;
my $org1 = "hsapiens_gene_ensembl";
my $org2 = "mmusculus";
my ($folder, $name) = mitochy::getFilename($input, "folder");
my $confFile = "/usr/local/bin/Perl/biomart/registry.xml";
#
# NB: change action to 'clean' if you wish to start a fresh configuration  
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
#

my $action='cached';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
my $registry = $initializer->getRegistry;

my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');

my @gene = @{process_input($input)};

$query->setDataset("$org1");
#$query->addFilter("ensembl_transcript_id", \@gene);
$query->addFilter("ensembl_transcript_id", ["ENST00000335137"]);
$query->setDataset("hsapiens_gene_ensembl");
$query->addAttribute("ensembl_transcript_id");
$query->addAttribute("mmusculus_homolog_canonical_transcript_protein");
$query->addAttribute("mmusculus_homolog_perc_id");
$query->addAttribute("mmusculus_homolog_perc_id_r1");
$query->addAttribute("mmusculus_homolog_ensembl_gene");
$query->addAttribute("mmusculus_homolog_orthology_confidence");

$query->formatter("TSV");

my $query_runner = BioMart::QueryRunner->new();
############################## GET COUNT ############################
# $query->count(1);
# $query_runner->execute($query);
# print $query_runner->getCount();
#####################################################################


############################## GET RESULTS ##########################
# to obtain unique rows only
# $query_runner->uniqueRowsOnly(1);

$query_runner->execute($query);
open (OUT, ">", "$org2\_$name.out") or die "Cannot write to $org2\_$name.out: $!\n";
$query_runner->printResults(\*OUT);
#####################################################################


sub process_input {
	my ($input) = @_;
	my @gene;
	open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
	while (my $line = <$in>) {
		chomp($line);
		if ($input =~ /.bed/) {
			my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
			push(@gene, $name);
		}
	}
	return(\@gene);
}
