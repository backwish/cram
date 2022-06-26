#include "memory.h"

typedef struct {
  i64 n; // number of bytes
  int bs; // block size
  darray *da;

} CRAMBG;

CRAMBG *crambg_initialize(char *filename, int bs);
i64 crambg_usedmemory(CRAMBG *cram);
void crambg_read(CRAMBG *cram, i64 s, i64 t, uchar *buf);
void crambg_write(CRAMBG *cram, i64 s, i64 t, uchar *buf);
void crambg_save(CRAMBG *cram, uchar *filename, uchar *memfilename);
CRAMBG *crambg_load(uchar *filename, uchar *memfilename);
void crambg_free(CRAMBG *cram);
