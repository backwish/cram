CC = gcc
CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O9 -fPIC
cram_test: cram.c cram.h memory.c memory.h
	$(CC)  $(CFLAGS) -o cram_test cram_test.c cram.c memory.c huffman.c buddy.c mmap.c heap.c -lm -pthread
huffman_test: 
	$(CC)  $(CFLAGS) -o huffman_test huffman_test.c huffman.c memory.c mmap.c heap.c -lm -pthread
huffman_test2: 
	$(CC)  $(CFLAGS) -o huffman_test2 huffman_test2.c huffman.c memory.c mmap.c heap.c -lm -pthread
buddy_test:
	$(CC)  $(CFLAGS) -o buddy_test buddy_test.c buddy.c -lm -pthread
clean:
	rm cram_test huffman_test
