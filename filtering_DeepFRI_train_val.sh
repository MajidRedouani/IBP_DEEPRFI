#first step: make dataset that contains only training and validation sequences used in DeepFRI model learing process


#we take the entry names of the file that contains training validation and test data
grep "^>" nrPDB-GO_2019.06.18_sequences.fasta | sed 's/>//' > main_headers.txt

#take entry names of file that contains only test data and remove the suffix "_test"
grep "^>" nrPDB-GO_2019.06.18_test_sequences.fasta | sed 's/>//; s/_test//' > formatted_test_headers.txt

# create new file that has only training and validation data:  entry name with corresponding sequences
comm -23 <(sort main_headers.txt) <(sort formatted_test_headers.txt) > headers_to_keep.txt
 awk '
    BEGIN {
        # Load headers to keep from headers_to_keep.txt into an array
        while ((getline line < "headers_to_keep.txt") > 0) {
            keep[line] = 1
        }
        RS=">"       # Set the record separator to ">"
        ORS=""       # Output separator is empty to concatenate records properly
    }
    NR > 1 {         # Skip the first empty record
        split($0, lines, "\n")  # Split each record by newline into an array
        header = lines[1]       # The first line is the header

        # Remove the header from the remaining lines to get only the sequence
        if (header in keep) {   # If header is in the keep list
            print ">" header "\n"  # Print header
            for (i = 2; i <= length(lines); i++) {
                if (lines[i] ~ /^>/) break  # Stop if the line starts with ">"
                print lines[i] "\n"         # Print each line of the sequence
            }
        }
    }
' nrPDB-GO_2019.06.18_sequences.fasta > train_val_sequences.fasta

#delete whitelines between sequences
sed '/^$/d' train_val_sequences.fasta > clean_train_val_sequences.fasta

#second step: blast sequences from our dataset against the sequences of the training and validation data used in DeepFRI


#extract sequence from original dataset as fasta
awk 'NR > 1 {print ">" $1 "\n" $3}' "remaining_sequences (1).tsv" > remaining_sequences.fasta

#prepare dataset with training and validation data to be used as a blast database
makeblastdb -in clean_train_val_sequences.fasta -dbtype prot -out train_val_db

#blast our sequences against training and validation data
blastp -query remaining_sequences.fasta -db train_val_db -outfmt "6 qseqid sseqid pident length evalue bitscore" -evalue 1e-2 > blast_results.tsv

#select id's that had a similarity score higher than 70% 
awk '$3 >= 70 && $4 >= (0.5 * length($3)) {print $1}' blast_results.tsv | sort | uniq > filtered_ids.txt

#remove these entries from our original dataset
awk 'NR==FNR {ids[$1]; next} NR==1 || !($1 in ids)' filtered_ids.txt remaining_sequences.tsv > dataset_filtered.tsv


 
