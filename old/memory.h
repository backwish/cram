#include "typedef.h"
#include "mman.h"

typedef byte mtype;
//typedef i64 ptype;
typedef u16 ofs_t;

typedef struct {
  int a, b; // min/max size of blocks (in bytes)

  i64 *List; // heads of lists of segments
  i64 Free; // head of list of free segments

  size_t memsize; // current size of memory (#segments)
  mtype *mem;
//  i64 freesize;
  int k;
  MMAP *map;
} dblock;

typedef struct {
  mtype *addr[2];
  ofs_t len[2];
} dblock_entry;

dblock *dblock_initialize(int a, int b, i64 m, int k);
dblock_entry dblock_allocate_block(dblock *db, int s);
mtype *dblock_addr(dblock_entry *blk, int i);
i64 dblock_free_block(dblock *db, dblock_entry *blk);
i64 dblock_usedmemory(dblock *db);
i64 block_seg(dblock *db, dblock_entry *blk);
i64 block_pos(dblock *db, dblock_entry *blk);
i64 block_len(dblock_entry *blk);

void dblock_write(dblock *db, FILE *out, FILE *mout);
dblock *dblock_read(uchar **map, uchar *memfilename);
void dblock_free(dblock *db);


typedef struct {
  i64 n; // number of blocks
  int k; // number of bytes to represent n

//  i64 *Seg; // segment storing i-th block
  uchar *Seg; // segment storing i-th block
  ofs_t *Pos; // relative location of i-th block in Seg[i]
//  ofs_t *Len; // length (in bytes) of i-th block

  dblock *db;
} darray;

darray *darray_initialize(i64 n, int a, int b, i64 m, int k);
//int darray_blocksize(darray *da, i64 i);
dblock_entry darray_address(darray *da, i64 i);
void darray_change(darray *da, i64 i, int new_size);
i64 darray_usedmemory(darray *da);
void darray_check(darray *da);

void darray_write(darray *da, FILE *out, FILE *mout);
darray *darray_read(uchar **map, uchar *memfilename);
void darray_free(darray *da);


/**********************************************************
A segment consists of
 size: size (bytes) of blocks in the segment
 offset: position of the first block in the segment
   Note.  [0..offset-1] may be used for another block, which begins in the prev segment,
          if the segment is not the first one in the list
 prev: pointer to its previous segment
 next: pointer to its next segment
 block: block_data

All live segments are stored in [1..last]


A block_data consists of
  id: index of the block
  data: block itself

 **********************************************************/
