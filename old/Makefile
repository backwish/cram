all: mkcram cram_test crambg_test

#CC = gcc43
CC = gcc

#CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -fPIC -g
CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O9 -fPIC

#CSALIB_DIR = ../csalib

mkcram: cram.c cram.h memory.c memory.h
	$(CC)  $(CFLAGS) -o mkcram cram.c memory.c huffman.c mmap.c heap.c -lm -pthread -DMKIDX_MAIN 
cram_test: cram.c cram.h memory.c memory.h
	$(CC)  $(CFLAGS) -o cram_test cram.c memory.c huffman.c mmap.c heap.c -lm -pthread -DTEST_MAIN
crambg_test: cram_bg.c cram_bg.h memory.c memory.h
	$(CC)  $(CFLAGS) -o crambg_test cram_bg.c memory.c huffman.c mmap.c heap.c -lz

clean:
	rm -f mkcram cram_test crambg_test
