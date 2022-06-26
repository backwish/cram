#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "memory.h"

#define ID_DBLOCK     0x08
#define ID_DARRAY     0x09

#define USE_MMAP 1


#define mymalloc(p,n) {p = malloc((n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory in line %d\n",__LINE__); exit(1);};}

#ifndef min
 #define min(x,y) (((x)<(y))?(x):(y))
#endif



static u64 getuint(uchar *s, i64 w)
{
  u64 x;
  i64 j;
  x = 0;
  for (j=0; j<w; j++) {
    x += ((u64)(*s++)) << (j*8);
  }
  return x;
}

static void putuint(uchar *s, u64 x, i64 w)
{
  i64 j;
  for (j=0; j<w; j++) {
    *s++ = x & 0xff;
    x >>= 8;
  }
}

static u64 getuint_block(dblock_entry *s, int i, i64 w)
{
  u64 x;
  i64 j;
  x = 0;
  for (j=0; j<w; j++) {
    x += ((u64)(*dblock_addr(s, i+j))) << (j*8);
  }
  return x;
}

static void putuint_block(dblock_entry *s, int i, u64 x, i64 w)
{
  i64 j;
  for (j=0; j<w; j++) {
    *dblock_addr(s, i+j) = x & 0xff;
    x >>= 8;
  }
}

static u64 readuint(uchar *s, i64 w)
{
  u64 x;
  i64 j;
  x = 0;
  for (j=0; j<w; j++) {
    x += ((u64)(*s++)) << (j*8); // little endian
  }
  return x;
}

static void writeuint(int k, u64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc((int)(x & 0xff),f); // little endian
    x >>= 8;
  }
}


#define GETSEG(da, i) getuint((da)->Seg + (i)*da->k, da->k)
#define SETSEG(da, i, x) putuint((da)->Seg + (i)*da->k, (x), da->k)
#define GETPOS(da, i) (da->Pos[i])
#define SETPOS(da, i, x) (da->Pos[i] = x)

#define SSIZ 0
#define SOFS (SSIZ+sizeof(ofs_t))
#define SPREV (SOFS+sizeof(ofs_t))
#define SNEXT (SPREV+sizeof(i64))
#define SBLK (SNEXT+sizeof(i64))

#define GETSIZE(mp) getuint(&(mp)[SSIZ], sizeof(ofs_t))
#define GETOFFSET(mp) getuint(&(mp)[SOFS], sizeof(ofs_t))
#define GETPREV(mp) getuint(&(mp)[SPREV], sizeof(i64))
#define GETNEXT(mp) getuint(&(mp)[SNEXT], sizeof(i64))
#define SETSIZE(mp,x) putuint(&(mp)[SSIZ], (x), sizeof(ofs_t))
#define SETOFFSET(mp,x) putuint(&(mp)[SOFS], (x), sizeof(ofs_t))
#define SETPREV(mp,x) putuint(&(mp)[SPREV], (x), sizeof(i64))
#define SETNEXT(mp,x) putuint(&(mp)[SNEXT], (x), sizeof(i64))

//#define BLOCK_SIZE(db) ((db)->b*3)
#define BLOCK_SIZE(db) ((db)->b)
#define SEGMENT_SIZE(db) (SBLK+BLOCK_SIZE(db))
//#define SEGMENT_SIZE(db) (1<<(8*sizeof(ofs_t)))
//#define BLOCK_SIZE(db) (SEGMENT_SIZE(db)-SBLK)
#define SEGMENT_ADDR(db,i) (&(db)->mem[i*SEGMENT_SIZE(db)])

#define NIL 0

#define GETIDX(b,k) getuint_block((b), 0, (k))
#define SETIDX(b,x,k) putuint_block((b), 0, (x), (k))

static int blog(i64 x) // blog(n)+1 bits are necessary and sufficient for storing a number in [0,n]
{
int l;
  l = -1;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
}

dblock *dblock_initialize(int a, int b, i64 initial_size, int k)
{
  dblock *db;
  i64 i, ss, ms;
  mtype *mp;

  mymalloc(db, 1);
  db->a = a;  db->b = b;

  mymalloc(db->List, b-a+1);
  for (i=a; i<=b; i++) db->List[i-a] = NIL;

  ss = SEGMENT_SIZE(db);
  db->memsize = ms = initial_size / ss;
  db->k = k/ss + 1;
#if USE_MMAP
  db->map = mymmap_w("mmaptmp", ms*ss);
//  db->map = mymmap_w("/dev/zero", ms*ss);
//  db->map = mymmap_anom(ms*ss);
  db->mem = db->map->addr;
#else
  mymalloc(db->mem, ms * ss);
  db->map = NULL;
#endif
#if 0
  for (i=0; i<ms*ss; i++) db->mem[i] = 0;
#else
//  memset(db->mem, 0, ms*ss*sizeof(db->mem[0]));
  for (i=0; i<ss; i++) {
    printf("write %ld \r",i);  fflush(stdout);
    memset(&db->mem[ms*i], 0, ms*sizeof(db->mem[0]));
  }
#endif
  db->Free = 1; // the first free segment

  return db;
}

void dblock_free(dblock *db)
{
  free(db->List);
  if (db->map != NULL) {
    mymunmap(db->map);
  } else {
    free(db->mem);
  }
  free(db);
}

void dblock_write(dblock *db, FILE *out, FILE *mout)
{
  i64 i, s, m;

  writeuint(1, ID_DBLOCK, out);
  writeuint(sizeof(db->a), db->a, out);
  writeuint(sizeof(db->b), db->b, out);
  writeuint(sizeof(db->k), db->k, out);
////  writeuint(sizeof(db->memsize), db->memsize, out);
  writeuint(sizeof(db->memsize), db->Free, out);
//  writeuint(sizeof(db->Free), db->Free, out);
  for (i=db->a; i<=db->b; i++) writeuint(sizeof(db->List[0]), db->List[i-db->a], out);
#if 0
  fwrite(SEGMENT_ADDR(db,0), SEGMENT_SIZE(db), db->Free, mout);
#else
  s = 1 << 16;
  for (i=0; i<db->Free; i+=s) {
    fprintf(stderr, "write %ld/%ld \r", i/s, db->Free/s);  fflush(stderr);
    m = s;  if (i+m > db->Free) m = db->Free - i;
    fwrite(SEGMENT_ADDR(db,i), SEGMENT_SIZE(db), m, mout);
  }
#endif
}

dblock *dblock_read(uchar **map, uchar *memfilename)
{
  dblock *db;
  uchar *p;
  int id;
  i64 i;
  
  mymalloc(db, 1);
  p = *map;
  
  if ((id = readuint(p,1)) != ID_DBLOCK) {
    printf("dblock_read: id = %d\n",id);
    exit(1);
  }
  p += 1;
  db->a = readuint(p,sizeof(db->a));  p += sizeof(db->a);
  db->b = readuint(p,sizeof(db->b));  p += sizeof(db->b);
  db->k = readuint(p,sizeof(db->k));  p += sizeof(db->k);
  db->memsize = readuint(p,sizeof(db->memsize));  p += sizeof(db->memsize);
  db->Free = db->memsize;

  mymalloc(db->List, db->b - db->a + 1);
  for (i=db->a; i<=db->b; i++) {
    db->List[i-db->a] = readuint(p, sizeof(db->List[0]));  p += sizeof(db->List[0]);
  }

  db->map = mymmap_w(memfilename, db->memsize * SEGMENT_SIZE(db));
  db->mem = db->map->addr;

  *map = p;
  return db;
}

static i64 allocate_new_segment(dblock *db)
{
  i64 free_seg;
  size_t new_memsize;
  mtype *mp;
  i64 i, s;
  
  free_seg = db->Free;
  if (free_seg == db->memsize) {
    new_memsize = db->memsize + db->k;
//    fprintf(stderr,"allocate_new_segment: memsize %ld -> %ld (%ld bytes)\n", db->memsize, new_memsize, new_memsize*SEGMENT_SIZE(db));
    if (db->map != NULL) {
      mp = mymremap(db->map, new_memsize * SEGMENT_SIZE(db));
    } else {
      mp = realloc(db->mem, new_memsize * SEGMENT_SIZE(db));
    }
    if (mp == NULL) {
      printf("allocate_new_segment: cannot realloc\n");
      exit(1);
    }
    if (mp != db->mem) {
//      fprintf(stderr,"allocate_new_segment: memory has moved.\n"); // not fatal
      db->mem = mp;
    }
    mp = SEGMENT_ADDR(db, db->memsize);
    for (i=0; i<SEGMENT_SIZE(db)*db->k; i++) *mp++ = 0;
    db->Free = db->memsize;
    db->memsize = new_memsize;
  }
  db->Free++;
  return free_seg;
}

static i64 free_segment(dblock *db, i64 i)
{
  i64 j, last, k;
  i64 p, n;
  mtype *mi, *mj;
  i64 freesize;

  j = db->Free - 1;
  db->Free--;
  if (i != j) { // move the last live segment j to i
    mj = SEGMENT_ADDR(db, j);
    p = GETPREV(mj);
    n = GETNEXT(mj);
    if (p != NIL) {
      SETNEXT(SEGMENT_ADDR(db, p), i);
      SETPREV(SEGMENT_ADDR(db, i), p); // !
    } else { // j was the head of List[sj];
      int sj;
      sj = GETSIZE(mj);
      db->List[sj - db->a] = i;
    }
    if (n != NIL) {
      SETPREV(SEGMENT_ADDR(db, n), i);
      SETNEXT(SEGMENT_ADDR(db, i), n); // !
    }

    mi = SEGMENT_ADDR(db, i);
    for (k = 0; k < SEGMENT_SIZE(db); k++) {
      *mi++ = *mj++;
    }

#if 0 // for debug
    mj = SEGMENT_ADDR(db, j) + SBLK;
    for (k = 0; k < db->b; k++) {
      *mj++ = 0xee;
    }
#endif
  } else {
    i = NIL; // !
  }

  // shrink the memory
  freesize = db->memsize - db->Free;
  if (freesize >= 2*db->k) {
    mtype *mp;
    size_t new_memsize;
    new_memsize = db->memsize - db->k;
//    fprintf(stderr,"free_segment: memsize %ld -> %ld\n", db->memsize, new_memsize);
    if (db->map != NULL) {
      mp = mymremap(db->map, new_memsize * SEGMENT_SIZE(db));
    } else {
      mp = realloc(db->mem, new_memsize * SEGMENT_SIZE(db));
    }
    if (mp == NULL) {
      printf("free_segment: cannot realloc\n");
      exit(1);
    }
    if (mp != db->mem) {
//      fprintf(stderr,"free_segment: memory has moved.\n"); // not fatal
      db->mem = mp;
    }
    db->memsize = new_memsize;
  }
  return i;
}

dblock_entry dblock_allocate_block(dblock *db, int s)
{
  mtype *mp;
  i64 i, j, ss;
  int ofs;
  dblock_entry new_block;

  if (s < db->a || s > db->b) {
    printf("dblock_allocate_block: size %d is out of [%d, %d]\n", s, db->a, db->b);
    exit(1);
  }

  i = db->List[s - db->a]; // the first segment for size s
  ofs = GETOFFSET(SEGMENT_ADDR(db,i)); // local address of the first block in i-th segment
    // if the list is NULL, it points to segment 0, whose offset is 0.  So a new segment will be allocated.

  if (s <= ofs) { // new block can be stored in the segment
    new_block.addr[0] = SEGMENT_ADDR(db, i) + SBLK + ofs-s;
    new_block.len[0] = s;
    new_block.addr[1] = NULL;
    new_block.len[1] = 0;
    SETOFFSET(SEGMENT_ADDR(db, i), ofs-s);
  } else { // new block cannot be stored
    j = allocate_new_segment(db);
    SETSIZE(SEGMENT_ADDR(db, j), s);
    db->List[s - db->a] = j; // new segment is the head of List[s]
    SETPREV(SEGMENT_ADDR(db, j), NIL);
    SETNEXT(SEGMENT_ADDR(db, j), i);
    SETPREV(SEGMENT_ADDR(db, i), j);

    new_block.addr[1] = (ofs > 0) ? SEGMENT_ADDR(db, i) + SBLK : NULL;
    new_block.len[1] = ofs;
//    SETOFFSET(SEGMENT_ADDR(db, i), 0); // !
//    new_block.addr[0] = SEGMENT_ADDR(db, j) + SBLK + (db->b - (s-ofs));
    new_block.addr[0] = SEGMENT_ADDR(db, j) + SEGMENT_SIZE(db) - (s-ofs);
    new_block.len[0] = s - ofs;
//    SETOFFSET(SEGMENT_ADDR(db, j), db->b - (s-ofs));
    SETOFFSET(SEGMENT_ADDR(db, j), BLOCK_SIZE(db) - (s-ofs));
  }
  return new_block;
}

i64 dblock_usedmemory(dblock *db)
{
  i64 size;
  int i, s;
  size = sizeof(dblock);
  size += (db->b - db->a + 1) * sizeof(*db->List);
#if 0
  printf("dblock: List %ld bytes\n", (db->b - db->a + 1) * sizeof(*db->List));
  s = 0;
  for (i=db->a; i<=db->b; i++) {
    if (db->List[i-db->a] != NIL) s++;
  }
  printf("number of lists %d redundancy %d bytes\n", s, s*SEGMENT_SIZE(db));
#endif
  size += SEGMENT_SIZE(db) * db->Free;
  return size;
}

i64 block_seg(dblock *db, dblock_entry *blk)
{
  return (blk->addr[0] - db->mem) / SEGMENT_SIZE(db);
}

i64 block_pos(dblock *db, dblock_entry *blk)
{
  return (blk->addr[0] - db->mem) % SEGMENT_SIZE(db);
}

i64 block_len(dblock_entry *blk)
{
  return blk->len[0] + blk->len[1];
}

mtype *dblock_addr(dblock_entry *blk, int i)
{
  if (i < 0) {
    printf("dblock_addr: i = %d\n", i);
    exit(1);
  }
  if (i < blk->len[0]) return blk->addr[0]+i;
  if (i >= blk->len[0] + blk->len[1]) {
    printf("dblock_addr: i = %d size = %d\n", i, blk->len[0] + blk->len[1]);
    exit(1);
  }
  return blk->addr[1] + i - blk->len[0];
}

dblock_entry block_addr_ofs(dblock_entry *blk, int ofs)
{
  dblock_entry blk2;

  if (ofs < blk->len[0]) {
    blk2.addr[0] = blk->addr[0] + ofs;
    blk2.len[0] = blk->len[0] - ofs;
    blk2.addr[1] = blk->addr[1];
    blk2.len[1] = blk->len[1];
  } else {
    blk2.addr[0] = blk->addr[1] + (ofs - blk->len[0]);
    blk2.len[0] = blk->len[1] - (ofs - blk->len[0]);
    blk2.addr[1] = NULL;
    blk2.len[1] = 0;
  }
  return blk2;
}


static dblock_entry first_block(dblock *db, i64 i) // first block in segment i
{
  i64 j;
  int ofs;
  dblock_entry blk;
  mtype *mp;
  int s;

  s = GETSIZE(SEGMENT_ADDR(db,i));
  ofs = GETOFFSET(SEGMENT_ADDR(db,i));
  if (s <= BLOCK_SIZE(db) - ofs) {
    blk.addr[0] = SEGMENT_ADDR(db, i) + SBLK + ofs;
    blk.len[0] = s;
    blk.addr[1] = NULL;
    blk.len[1] = 0;
  } else {
    blk.addr[0] = SEGMENT_ADDR(db, i) + SBLK + ofs;
    blk.len[0] = BLOCK_SIZE(db) - ofs;
    j = GETNEXT(SEGMENT_ADDR(db,i));
    mp = SEGMENT_ADDR(db,j);
    blk.addr[1] = mp + SBLK;
    blk.len[1] = s - blk.len[0];
  }
  return blk;
}

static dblock_entry next_block(dblock *db, dblock_entry *blk0)
{
  i64 i, j;
  int ofs;
  dblock_entry blk;
  mtype *mp;
  int s;

  i = block_seg(db, blk0);
  ofs = block_pos(db, blk0);
  s = block_len(blk0);

  if (ofs + s < SEGMENT_SIZE(db)) { // next block starts in the same segment
    blk.addr[0] = blk0->addr[0] + s;
    if (ofs + s + s <= SEGMENT_SIZE(db)) { // next block fits in the same segment
      blk.len[0] = s;
      blk.addr[1] = NULL;
      blk.len[1] = 0;
    } else {
      blk.len[0] = SEGMENT_SIZE(db) - (ofs + s);
      j = GETNEXT(SEGMENT_ADDR(db,i));
      mp = SEGMENT_ADDR(db,j);
      blk.addr[1] = mp + SBLK;
      blk.len[1] = s - blk.len[0];
    }
  } else {
    j = GETNEXT(SEGMENT_ADDR(db,i));
    mp = SEGMENT_ADDR(db,j);
    blk.addr[0] = mp + SBLK + blk0->len[1];
    blk.len[0] = s;
    blk.addr[1] = NULL;
    blk.len[1] = 0;
  }
  return blk;
}



i64 dblock_free_block(dblock *db, dblock_entry *blk) // blk may change
{
  int s, ss;
  i64 i,j, f;
  dblock_entry first_blk;
  mtype *mp, *q, *memtmp;
  int ofs;

  s = blk->len[0] + blk->len[1];
  i = db->List[s - db->a];
  first_blk = first_block(db, i);

  // move the first block
  if (first_blk.addr[0] != blk->addr[0]) {
    for (j=0; j<s; j++) {
      *dblock_addr(blk, j) = *dblock_addr(&first_blk, j);
//      *dblock_addr(&first_blk, j) = 0xff; // for debug
    }
  } else {
    blk->addr[0] = blk->addr[1] = NULL;
    blk->len[0] = blk->len[1] = 0;
  }

  memtmp = db->mem;
  
  // free the first block
  f = NIL;
  if (first_blk.len[1] == 0) { // the first block fits in the first segment
    mp = SEGMENT_ADDR(db,i);
    ofs = GETOFFSET(mp);
    if (ofs + s < BLOCK_SIZE(db)) { // the first segment does not become empty
      SETOFFSET(mp, ofs+s);
    } else {
      j = GETNEXT(SEGMENT_ADDR(db,i));
      db->List[s - db->a] = j;
      SETPREV(SEGMENT_ADDR(db,j), NIL);
      f = free_segment(db, i);
    }
  } else { // first block lies on two segments
    j = GETNEXT(SEGMENT_ADDR(db,i));
    mp = SEGMENT_ADDR(db,j);
    SETOFFSET(mp, first_blk.len[1]); // !
    db->List[s - db->a] = j;
    SETPREV(SEGMENT_ADDR(db,j), NIL);
    f = free_segment(db, i);
  }

  if (db->mem != memtmp) {
    if (blk->addr[0] != NULL) blk->addr[0] += db->mem - memtmp;
    if (blk->addr[1] != NULL) blk->addr[1] += db->mem - memtmp;
  }

  return f;
}

static void dump_memory(dblock *db)
{
  i64 i,j,ss;
  mtype *mp;

  ss = SEGMENT_SIZE(db);
  for (i=0; i<db->memsize; i++) {
    mp = SEGMENT_ADDR(db,i);
    printf("segment %ld  addr = %p\n", i, mp);
    printf("  size %d offset %d prev %ld next %ld\n", GETSIZE(mp), GETOFFSET(mp), GETPREV(mp), GETNEXT(mp));
    printf("  data:");
    for (j=SBLK; j<ss; j++) {
      printf(" %02x(%c)", mp[j], mp[j]);
    }
    printf("\n\n");
  }
  for (i=db->a; i<=db->b; i++) {
    printf("List %d:", i);
    j = db->List[i - db->a];
    while (j != NIL) {
      printf(" %d", j);
      j = GETNEXT(SEGMENT_ADDR(db,j));
    }
    printf("\n");
  }
  printf("Free %ld  memsize %ld  \n", db->Free, db->memsize);
}





darray *darray_initialize(i64 n, int a, int b, i64 m, int k)
{
  darray *da;
  dblock *db;
  i64 i;
  mtype *mp;

  mymalloc(da, 1);
  da->n = n;
//  da->k = sizeof(da->n);
  da->k = (blog(n+1)+1+8-1)/8;

  da->db = dblock_initialize(a+da->k, b+da->k, m, k);

//  mymalloc(da->Seg, n);
  mymalloc(da->Seg, n * da->k);
  mymalloc(da->Pos, n);

  for (i=0; i<n; i++) {
//    if (i % (n/100) == 0) {
//      fprintf(stderr, "%ld \r", i/(n/100));  fflush(stderr);
//    }
//    da->Seg[i] = NIL;
//    da->Pos[i] = 0;
    SETSEG(da, i, NIL);
    SETPOS(da, i, 0);
  }
  return da;
}

void darray_free(darray *da)
{
  dblock_free(da->db);
  free(da->Seg);
  free(da->Pos);
  free(da);
}

void darray_write(darray *da, FILE *out, FILE *mout)
{
  i64 i;
  int k;

  writeuint(1, ID_DARRAY, out);
  k = da->k;
  writeuint(1, k, out);
  writeuint(k, da->n, out);
  for (i=0; i<da->n; i++) {
    writeuint(k, GETSEG(da, i), out);
  }
  for (i=0; i<da->n; i++) {
    writeuint(sizeof(da->Pos[i]), GETPOS(da, i), out);
  }
  dblock_write(da->db, out, mout);
}

darray *darray_read(uchar **map, uchar *memfilename)
{
  darray *da;
  int k, id;
  i64 i, n;
  uchar *p;
  
  mymalloc(da, 1);
  p = *map;
  
  if ((id = readuint(p,1)) != ID_DARRAY) {
    printf("dblock_read: id = %d\n",id);
    exit(1);
  }
  p += 1;
  da->k = k = readuint(p,1);  p += 1;
  da->n = n = readuint(p,k);  p += k;

  mymalloc(da->Seg, n * da->k);
  mymalloc(da->Pos, n);

  for (i=0; i<da->n; i++) {
    SETSEG(da, i, readuint(p, k));  p += k;
  }
  for (i=0; i<da->n; i++) {
    SETPOS(da, i, readuint(p, sizeof(da->Pos[0])));  p += sizeof(da->Pos[0]);
  }

  da->db = dblock_read(&p, memfilename);
  *map = p;
  return da;
}

#if 0
int darray_blocksize(darray *da, i64 i)
{
  i64 seg;

//  seg = da->Seg[i];
  seg = GETSEG(da, i);
  if (seg == NIL) return 0;
  return GETSIZE(SEGMENT_ADDR(da->db, seg)) - da->k;
}
#endif

static dblock_entry darray_address_internal(darray *da, i64 i)
{
  dblock_entry blk;
  i64 seg, ofs, j;
  int s;
  mtype *mp;

  if (i < 0 || i >= da->n) {
    printf("darray_address: i = %ld n = %ld\n", i, da->n);
    exit(1);
  }

//  seg = da->Seg[i];
  seg = GETSEG(da, i);

  if (seg == NIL) {
    printf("darray_address: block %ld is not allocated.\n", i);
    exit(1);
  }

  ofs = GETPOS(da, i);
  mp = SEGMENT_ADDR(da->db, seg);
  s = GETSIZE(mp);

  if (ofs + s <= SEGMENT_SIZE(da->db)) {
    blk.addr[0] = mp + ofs;
    blk.len[0] = s;
    blk.addr[1] = NULL;
    blk.len[1] = 0;
  } else {
    blk.addr[0] = mp + ofs;
    blk.len[0] = SEGMENT_SIZE(da->db) - ofs;
    j = GETNEXT(SEGMENT_ADDR(da->db,seg));
    mp = SEGMENT_ADDR(da->db, j);
    blk.addr[1] = mp + SBLK;
    blk.len[1] = s - blk.len[0];
  }

  return blk;
}

dblock_entry darray_address(darray *da, i64 i)
{
  dblock_entry blk;

  blk = darray_address_internal(da, i);
  return block_addr_ofs(&blk, da->k);
}

#if 1
static void darray_check_i(darray *da, i64 i)
{
  int j,s;
  dblock_entry blk;

//  if (da->Seg[i] != NIL) {
  if (GETSEG(da, i) != NIL) {
    blk = darray_address_internal(da, i);
    j = GETIDX(&blk, da->k);
    if (i != j) {
      printf("darray_check_i: i = %ld inverse = %ld\n", i, j);
    }
  }
}

void darray_check(darray *da)
{
  i64 i;
  int j,s;
  dblock_entry blk;
  mtype tmp[sizeof(mtype *)];

  for (i=0; i<da->n; i++) {
//    if (da->Seg[i] != NIL) {
    if (GETSEG(da, i) != NIL) {
      blk = darray_address_internal(da, i);
      j = GETIDX(&blk, da->k);
      if (i != j) {
        printf("darray_check: i = %ld inverse = %ld\n", i, j);
      }
    }
  }
}
#endif

void darray_change(darray *da, i64 i, int new_size)
{
  dblock_entry blk1, blk2, blk3;
  i64 moved_segment;
  i64 j;
  int s, ofs;
  mtype tmp[sizeof(mtype *)];
  mtype *mp;

  if (GETSEG(da, i) == NIL) { // first allocation
    blk2 = dblock_allocate_block(da->db, new_size+da->k);
    putuint(tmp, i, da->k);
    for (j=0; j<da->k; j++) {
      *dblock_addr(&blk2, j) = tmp[j];  // write the index (i)
    }
    SETSEG(da, i, block_seg(da->db, &blk2));
    SETPOS(da, i, block_pos(da->db, &blk2));
  } else { // change the size
    blk2 = dblock_allocate_block(da->db, new_size+da->k); // new address
    blk1 = darray_address_internal(da, i); // current address

#if 0
    s = min(da->k + new_size, block_len(&blk1));
    for (j=0; j<s; j++) {
      *dblock_addr(&blk2, j) = *dblock_addr(&blk1, j);  // copy the index and content
    }
    for (j=s; j<da->k + new_size; j++) {
      *dblock_addr(&blk2, j) = 0; // for debug
    }
#else
//    for (j=0; j<da->k; j++) {
//      *dblock_addr(&blk2, j) = *dblock_addr(&blk1, j);  // copy the index
//    }
      SETIDX(&blk2, GETIDX(&blk1, da->k), da->k);  // copy the index
#endif
    SETSEG(da, i, block_seg(da->db, &blk2));
    SETPOS(da, i, block_pos(da->db, &blk2));

    moved_segment = dblock_free_block(da->db, &blk1); // free the old block.  blk1 may be broken

    // fix the reference to the block moved from the head of the list
    if (blk1.addr[0] != NULL) {
      j = GETIDX(&blk1, da->k);
      SETSEG(da, j, block_seg(da->db, &blk1));
      SETPOS(da, j, block_pos(da->db, &blk1));
    }

    if (moved_segment != NIL) { // fix the references to moved segment
      blk3 = first_block(da->db, moved_segment);
      while (1) {
        j = GETIDX(&blk3, da->k);
        SETSEG(da, j, block_seg(da->db, &blk3));
        SETPOS(da, j, block_pos(da->db, &blk3));
        blk3 = next_block(da->db, &blk3);
        if (block_seg(da->db, &blk3) != moved_segment) break;
      }
    }
  }

//  darray_check(da);
}

i64 darray_usedmemory(darray *da)
{
  i64 size;
  size = sizeof(darray);
  size += da->n * da->k * sizeof(*da->Seg);
//  printf("darray: Seg %ld bytes\n", da->n * da->k * sizeof(*da->Seg));
  size += da->n * sizeof(*da->Pos);
//  printf("darray: Pos %ld bytes\n", da->n * sizeof(*da->Pos));
  size +=dblock_usedmemory(da->db);
  return size;
}

#if 1
static void dump_darray(darray *da)
{
  i64 i;
  int j,s;
  dblock_entry blk;
  mtype tmp[sizeof(mtype *)];

  printf("darray n = %ld\n", da->n);
  for (i=0; i<da->n; i++) {
//    if (da->Seg[i] != NIL) {
    if (GETSEG(da, i) != NIL) {
      blk = darray_address_internal(da, i);

      for (j=0; j<da->k; j++) {
        tmp[j] = *dblock_addr(&blk, j);
      }
      j = getuint(tmp, da->k);

      s = block_len(&blk) - da->k;
      printf("block %ld[%ld] (len %d, seg %ld pos %d):", i, j, s, GETSEG(da, i), GETPOS(da, i));
      for (j=0; j<s; j++) {
        printf(" %02x", *dblock_addr(&blk, da->k + j));
      }
      printf("\n");
    }
  }
}
#endif

#if 0
int main(int argc, char *argv[])
{
  dblock *db;
  dblock_entry blk, blk2, blk3;
  i64 i,j;

  db = dblock_initialize(1, 10, 100, 1);
  dump_memory(db);

  blk = dblock_allocate_block(db, 1);
  *dblock_addr(&blk,0) = '1';
  blk = dblock_allocate_block(db, 1);
  *dblock_addr(&blk,0) = '2';
  blk = dblock_allocate_block(db, 1);
  *dblock_addr(&blk,0) = '3';

  blk = dblock_allocate_block(db, 2);
  *dblock_addr(&blk,0) = '4';
  *dblock_addr(&blk,1) = '5';
  blk2 = blk;

  for (i=0; i<4; i++) {
    blk = dblock_allocate_block(db, 3);
    for (j=0; j<3; j++) {
      *dblock_addr(&blk,j) = 'a'+i*3+j;
    }
    if (i == 0) blk3 = blk;
  }

  dump_memory(db);

  dblock_free_block(db, blk2);
  dump_memory(db);

  dblock_free_block(db, blk3);
  dump_memory(db);

  blk = dblock_allocate_block(db, 4);
  *dblock_addr(&blk,0) = '0';
  *dblock_addr(&blk,1) = '1';
  *dblock_addr(&blk,2) = '2';
  *dblock_addr(&blk,3) = '3';
  dump_memory(db);

  blk = dblock_allocate_block(db, 3);
  *dblock_addr(&blk,0) = '3';
  *dblock_addr(&blk,1) = '4';
  *dblock_addr(&blk,2) = '5';
  dump_memory(db);

}
#endif

#if 0
int main(int argc, char *argv[])
{
  darray *da;
  dblock_entry blk, blk2, blk3;
  i64 i,j;

  da = darray_initialize(10, 1, 10, 100, 1);
  dump_darray(da);
  dump_memory(da->db);

  darray_change(da, 0, 1);
  blk = darray_address(da, 0);
  *dblock_addr(&blk, 0) = '0';

  darray_change(da, 1, 1);
  blk = darray_address(da, 1);
  *dblock_addr(&blk, 0) = '1';

  darray_change(da, 2, 1);
  blk = darray_address(da, 2);
  *dblock_addr(&blk, 0) = '2';

  darray_change(da, 3, 2);
  blk = darray_address(da, 3);
  *dblock_addr(&blk, 0) = '4';
  *dblock_addr(&blk, 1) = '5';

  for (i=4; i<8; i++) {
    darray_change(da, i, 3);
    blk = darray_address(da, i);
    for (j=0; j<3; j++) {
      *dblock_addr(&blk, j) = 'a'+i*3+j;
    }
  }

  dump_darray(da);
  dump_memory(da->db);

  darray_change(da, 0, 5);

  dump_darray(da);
  dump_memory(da->db);

}
#endif
