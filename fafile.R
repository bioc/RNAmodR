fafile <- "sacCer3_R64-2-1.fasta"
gfffile <- "sacCer3_R64-2-1.gff3"

# download sacCer3 information
utils::download.file(paste0("http://downloads.yeastgenome.org/sequence/",
                            "S288C_reference/genome_releases/",
                            "S288C_reference_genome_R64-2-1_20150113.tgz"),
                     paste0("S288C_reference_genome_R64-2-1_20150113.tgz"))
utils::untar(paste0("S288C_reference_genome_R64-2-1_20150113.tgz"),
             exdir = ".")

# Remove fasta sequences attached to gff file
gff <- rtracklayer::import.gff3(paste0("S288C_reference_genome_R64-2-1_20150113/",
                                       "saccharomyces_cerevisiae_R64-2-1_20150113.gff"))
rtracklayer::export.gff3(gff, gfffile)

# move fasta file
file.rename(paste0("S288C_reference_genome_R64-2-1_20150113/",
                   "S288C_reference_sequence_R64-2-1_20150113.fsa"),
            fafile)

# fix fasta names
gff <- rtracklayer::import.gff3(gfffile)
fsa <- Biostrings::readDNAStringSet(fafile)
names(fsa) <- unique(as.character(rtracklayer::chrom(gff)))
Biostrings::writeXStringSet(fsa, fafile)