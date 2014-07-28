CC=gcc

all:
	$(CC) src/kmer_refset_process.c -o bin/kmer_refset_process -lm
	$(CC) src/kmer_jaccard_index.c -o bin/kmer_jaccard_index -lm
	$(CC) src/kmer_reads_process_stdin.c -o bin/kmer_reads_process_stdin -lm
	$(CC) src/intersect_kmer_lists_filelist.c -o bin/intersect_kmer_lists_filelist -lm
clean:
	rm bin/*

