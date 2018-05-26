#!/usr/bin/env perl
####!/usr/local/perls/perl-5.16.3/bin/perl
use Modern::Perl '2013';
use autodie;

# parse SNPs from multiple alignment fasta file.

die "Usage: $0 alignment_file.fasta\n" unless $ARGV[0];

my( $ref_name, $ref_seq ); # reference (query) sequence
my $ref_seq2; # ungapped reference sequence
my $flag = 0;

local $/ = "\n>";
my $sbjct_count = 0; # strain count
my $snps = {}; # sparse array ( HoH ) to hold SNPs
my $snps2 ={}; # tally position in original ungapped reference
my @strain_names;
my %snp_cnt; # hold strain index and count SNPs in each strain
my %snps3; # key-value = position in alignment-position in reference gene.

while( <> ) {
    chomp; # it chomps "\n>", not just \n because of $/
    next if ! length($_);
    
    if( $flag == 1 ) {
        # get name and sequence for each subject strain
        my( $sbjct_name, $sbjct_seq ) = get_data( $_ );
        push @strain_names, $sbjct_name;
        
        # populate sparse array with SNPs.
        # rows are SNP positions, cols are different strains
        for( my $i = 0; $i < length( $ref_seq ); ++$i ) {
            my $ref_base = substr( $ref_seq, $i, 1 );
			my $aln_pos = $i+1; # base position in alignment (including gaps)

			if( $ref_base ne substr( $sbjct_seq, $i, 1 ) ) {
				$snps->{ $aln_pos }{ $sbjct_count } = substr( $sbjct_seq, $i, 1 ); # populate sparse array with SNPs
				$snp_cnt{ $sbjct_count }++; # strain number => total SNPs in strain
				$snps2->{ $aln_pos }{ substr( $sbjct_seq, $i, 1 ) }++; # count SNPs at each position
			}
		}
        ++$sbjct_count; # increment strain number
    } 
    else {
        # get name and sequence for query strain
        ( $ref_name, $ref_seq ) = get_data( $_ );
        $flag = 1;
        push @strain_names, ( 'reference_' . $ref_name );
        $snp_cnt{ $sbjct_count } = ''; 
        ++$sbjct_count; # increment strain count
        
		# count position in original ungapped reference
		my $gap = 0; # substract this from align position to get position in original ungapped reference
		my $ref_pos; # base position in original reference (excluding gaps)
		for( my $i = 0; $i < length( $ref_seq ); ++$i ) {
            my $ref_base = substr( $ref_seq, $i, 1 );
			my $aln_pos = $i+1; # base position in alignment (including gaps)
			if( $ref_base ne '-' ) {
				$ref_pos = $aln_pos - $gap;
				$snps3{ $aln_pos } = $ref_pos;
			} 
			else {
				$snps3{ $aln_pos } = '-';
				++$gap;
			}
		}
		# remove gaps to get original reference sequence.
		( $ref_seq2 = $ref_seq ) =~ s/\-//g;

    }
}

# insert query sequence nt at position that has SNP
for my $pos ( sort{ $a <=> $b } keys %$snps ) {
    if( exists( $snps->{ $pos } ) ) {
        $snps->{ $pos }{0} = substr( $ref_seq, $pos-1, 1 );
    }
}

# evaluate synonymous/non-synonymous substitutions.
for my $i ( sort{ $a<=>$b } keys %$snps2 ) {
	my $ref_pos = $snps3{ $i };
	if( $ref_pos ne '-' )  { # skip insertions in reference
		for my $base( keys %{ $snps2->{ $i } } ) {
	        if( $base ne '-' ) {    # skip deletions in subject
            	my( $ref_cod, $sbjct_cod );
				if( ! (($ref_pos) % 3) ) { # 3rd position in codon
					$ref_cod = substr( $ref_seq2, ($ref_pos-3), 3 );
					$sbjct_cod = $ref_cod;
					substr( $sbjct_cod, 2, 1, $base );
				} 
				elsif( (($ref_pos) % 3 ) == 2 ) { # 2nd position in codon
					$ref_cod = substr( $ref_seq2, ($ref_pos-2), 3 );
					$sbjct_cod = $ref_cod;
					substr( $sbjct_cod, 1, 1, $base );
				} 
				elsif( (($ref_pos) % 3 ) == 1 ) { # 1st position in codon
					$ref_cod = substr( $ref_seq2, ($ref_pos-1), 3 );
					$sbjct_cod = $ref_cod;
					substr( $sbjct_cod, 0, 1, $base );
				}
				
				#if( $ref_cod !~ /\-/ ) {   # skip indels
					my $substitution = compare_codons( $ref_cod, $sbjct_cod );
					$snps2->{ $i }{ $base . "sub" } = $substitution;
				#} 
			} 
			else {
				next;
			}
		}
	} 
	else {
		next;
	}
}


# print reference name and names of strains that have snps.
say "nt_pos_in_ref\t",join( "\t", 
    ( ( map $strain_names[$_], 
    ( sort{ $a<=>$b } keys %snp_cnt) ), 
    "A", "C", "G", "T", "-", "A_sub", "C_sub", "G_sub", "T_sub" )
    );

# print SNP results
for my $pos( sort{ $a<=>$b } keys %$snps ) {
    if( exists( $snps->{ $pos } ) ) {
        print $snps3{ $pos },"\t";
        for my $strain( sort{ $a <=> $b } keys %snp_cnt ) {
            if( exists( $snps->{ $pos }{ $strain } ) ) {
                print $snps->{ $pos }{ $strain }, "\t";
            } 
            else {
                print "\t";
            }
        }
        print $snps2->{ $pos }{A} // '', "\t";
        print $snps2->{ $pos }{C} // '', "\t";
        print $snps2->{ $pos }{G} // '', "\t";
        print $snps2->{ $pos }{T} // '', "\t";
        print $snps2->{ $pos }{"-"} // '', "\t";
        print $snps2->{ $pos }{Asub} // '', "\t";
        print $snps2->{ $pos }{Csub} // '', "\t";
        print $snps2->{ $pos }{Gsub} // '', "\t";
        say $snps2->{ $pos }{Tsub} // '';
    }
}

# print SNPs per strain
say "\t", map{ "$snp_cnt{$_}\t" } sort{ $a<=>$b } keys %snp_cnt;
    


###########################
sub get_data {
# parse header and DNA 
# for seqeunce entries.
    
    my( $header, $seq );
    
    # use '?' because first fasta record has '>' in front
    # while it's chomped in all remaining.
    ( $header ) = /^>?(.*)\n/;
    $header =~ s/\s*//g;
    s/^>?.*\n//;
    s/\n//g;
    s/\s*//g;
    $seq = uc($_);
    return( $header, $seq );
}

sub compare_codons {
# compare codons

    my( $cod1, $cod2 ) = @_;

    # persistent reference to hash using state, 
    # cannot initialize state variable in list context.
    state $transl_table_11 = { 
                        TTT => 'F', # Phe      
                        TCT => 'S', # Ser      
                        TAT => 'Y', # Tyr      
                        TGT => 'C', # Cys  
                        TTC => 'F', # Phe      
                        TCC => 'S', # Ser      
                        TAC => 'Y', # Tyr      
                        TGC => 'C', # Cys  
                        TTA => 'L', # Leu      
                        TCA => 'S', # Ser      
                        TAA => '*', # STOP     
                        TGA => '*', # STOP 
                        TTG => 'L', # Leu START    
                        TCG => 'S', # Ser      
                        TAG => '*', # STOP      
                        TGG => 'W', # Trp  
                        CTT => 'L', # Leu      
                        CCT => 'P', # Pro      
                        CAT => 'H', # His      
                        CGT => 'R', # Arg  
                        CTC => 'L', # Leu      
                        CCC => 'P', # Pro      
                        CAC => 'H', # His      
                        CGC => 'R', # Arg  
                        CTA => 'L', # Leu      
                        CCA => 'P', # Pro      
                        CAA => 'Q', # Gln     
                        CGA => 'R', # Arg  
                        CTG => 'L', # Leu START    
                        CCG => 'P', # Pro      
                        CAG => 'Q', # Gln      
                        CGG => 'R', # Arg  
                        ATT => 'I', # Ile START  
                        ACT => 'T', # Thr      
                        AAT => 'N', # Asn      
                        AGT => 'S', # Ser  
                        ATC => 'I', # Ile START    
                        ACC => 'T', # Thr      
                        AAC => 'N', # Asn      
                        AGC => 'S', # Ser  
                        ATA => 'I', # Ile START    
                        ACA => 'T', # Thr      
                        AAA => 'K', # Lys      
                        AGA => 'R', # Arg  
                        ATG => 'M', # Met START    
                        ACG => 'T', # Thr      
                        AAG => 'K', # Lys      
                        AGG => 'R', # Arg  
                        GTT => 'V', # Val      
                        GCT => 'A', # Ala      
                        GAT => 'D', # Asp      
                        GGT => 'G', # Gly  
                        GTC => 'V', # Val      
                        GCC => 'A', # Ala      
                        GAC => 'D', # Asp      
                        GGC => 'G', # Gly  
                        GTA => 'V', # Val      
                        GCA => 'A', # Ala      
                        GAA => 'E', # Glu      
                        GGA => 'G', # Gly  
                        GTG => 'V', # Val START    
                        GCG => 'A', # Ala      
                        GAG => 'E', # Glu      
                        GGG => 'G', # Gly
                    };
    if( $transl_table_11->{ $cod1 } eq $transl_table_11->{ $cod2 } ) {
        my $sub = "synonym";
    } 
    else {
        my $sub = $$transl_table_11{ $cod1 } . "->" . $$transl_table_11{ $cod2 };
    }
}


