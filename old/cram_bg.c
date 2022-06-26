#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include "cram_bg.h"
#include "mman.h"

#define ID_CRAMBG     0x0b


#define mymalloc(p,n) {p = malloc((n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory in line %d\n",__LINE__); exit(1);};}

#ifndef min
 #define min(x,y) (((x)<(y))?(x):(y))
#endif
#ifndef max
 #define max(x,y) (((x)>(y))?(x):(y))
#endif

#define D (8*sizeof(uchar))

static u64 getuint(uchar *s, i64 i, i64 n, i64 w)
{
  u64 x,c;
  i64 j;

  x = 0;
  for (j=0; j<w; j++) {
    if (i+j < n) c = s[i+j]; else c = 0;
    //x += c << (j*8); // little endian
    x += c << ((w-1-j)*8); // big endian
  }
  return x;
}

static int setbit(dblock_entry *blk, i64 i,int x)
{
  i64 j,l;

  j = i / D;
  l = i % D;
  switch (x) {
  case 0:
    *dblock_addr(blk, j) &= (~(1L<<(D-1-l)));
    break;
  case 1:
    *dblock_addr(blk, j) |= (1L<<(D-1-l));
    break;
  default:
    printf("error setbit x=%d\n",x);
    exit(1);
  }
  return x;
}

static u64 setbits0(dblock_entry *blk, i64 i, int d, u64 x)
{
  int j;

  for (j=0; j<d; j++) {
    setbit(blk,i+j,(x>>(d-j-1))&1);
  }
  return x;
}

static u64 setbits(dblock_entry *blk, i64 i, int d, u64 x)
{
  u64 y,m;
  int d2;
  i64 iq, ir;
  mtype *mp;

  iq = i / D;
  ir = i % D;

  while (ir+d > D) {
    d2 = D-ir;
    y = x >> (d-d2);
    m = (1<<d2)-1;
//    *B = (*B & (~m)) | y;
    mp = dblock_addr(blk, iq);
    *mp = (*mp & (~m)) | y;
    iq++;  ir=0;
    d -= d2;
    x &= (1<<d)-1;
  }
  m = (1<<d)-1;
  y = x << (D-ir-d);
  m <<= (D-ir-d);
//  *B = (*B & (~m)) | y;
  mp = dblock_addr(blk, iq);
  *mp = (*mp & (~m)) | y;

  return x;
}


static u64 getbit32(dblock_entry *blk, i64 i)
{
  i64 j,l,x;
  int ii, bs;

  bs = block_len(blk);
  j = i / D;
  l = i % D;
  x = 0;
  for (ii=0; ii<4+1; ii++) {
    x <<= 8;
    if (j+ii < bs) {
      x += *dblock_addr(blk, j+ii);
    }
  }
  return (x >> (D-l)) & 0xffffffff;
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


CRAMBG *crambg_initialize(char *filename, int bs)
{
  CRAMBG *cram;
  MMAP *map;
  i64 i,j,n,nb, s, t;
  int k, sigma, c, c2;
  uchar *text;
  double *f2;
  i64 csize, csize2;
  darray *da;
  int len, r;
  u64 x;
  uchar *compbuf;
  dblock_entry blk;
  unsigned long compsize, origsize;
  
  mymalloc(cram,1);

  cram->bs = bs;

  map = mymmap(filename);
  cram->n = n = map->len;
  text = map->addr;

  nb = (n+bs-1)/bs; // number of blocks

  cram->da = da = darray_initialize(nb, 1, bs*2, n/4, n/100);

  mymalloc(compbuf, bs*10);

  csize = 0;
  for (i=0; i<nb; i++) { // compress i-th block
    if (i % 1000 == 0) {
      fprintf(stderr, "%ld/%ld \r", i, nb);  fflush(stderr);
    }
    s = i*bs;  t = min((i+1)*bs-1, n-1);
    origsize = t - s + 1;
    compsize = bs*10;
    if (compress2(compbuf, &compsize, &text[s], origsize, 1) != Z_OK) {
      error("compress()");
    }
    csize += compsize;
    darray_change(da, i, compsize); // allocate a block
    blk = darray_address(da, i);

    for (j=0; j<compsize; j++) {
      *dblock_addr(&blk, j) = compbuf[j];
    }
  }
  free(compbuf);

  fprintf(stderr,"n = %ld csize = %ld (%1.2f bpc)\n", n, csize, (double)csize*8/n);
  fprintf(stderr,"darray used memory %ld bytes (%1.2f bpc)\n", darray_usedmemory(cram->da), (double)darray_usedmemory(cram->da)*8/n);

  return cram;
}

void crambg_save(CRAMBG *cram, uchar *filename, uchar *memfilename)
{
  FILE *out, *mout;

  out = fopen(filename, "w");
  mout = fopen(memfilename, "w");
  if (out == NULL || mout == NULL) {
    perror("cram_save: fopen\n");
    exit(1);
  }

  writeuint(1, ID_CRAMBG, out);
  writeuint(sizeof(cram->n), cram->n, out);
  writeuint(sizeof(cram->bs), cram->bs, out);

  darray_write(cram->da, out, mout);

  fclose(out);
  fclose(mout);
}

CRAMBG *crambg_load(uchar *filename, uchar *memfilename)
{
  CRAMBG *cram;
  int k,i,j,c, id;
  MMAP *map;
  uchar *p;
  
  mymalloc(cram, 1);

  map = mymmap(filename);
  p = map->addr;
  
  if ((id = readuint(p,1)) != ID_CRAMBG) {
    printf("crambg_load: id = %d\n",id);
    exit(1);
  }
  p += 1;
  cram->n = readuint(p,sizeof(cram->n));  p += sizeof(cram->n);
  cram->bs = readuint(p,sizeof(cram->bs));  p += sizeof(cram->bs);

  cram->da = darray_read(&p, memfilename);

  return cram;
}

void crambg_free(CRAMBG *cram)
{
  darray_free(cram->da);
  free(cram);
}

i64 crambg_usedmemory(CRAMBG *cram)
{
  i64 size;
  size = sizeof(CRAMBG);
  size += darray_usedmemory(cram->da);
  return size;
}

static void crambg_read_block(CRAMBG *cram, i64 b, uchar *buf) // decode b-th block into buf
{
  i64 i,j,b1,b2,n;
  i64 i1, i2;
  i64 bs;
  int k, r, c;
  dblock_entry blk;
  u64 x;
  uchar *compbuf;
  unsigned long compsize, origsize;

  n = cram->n;
  bs = cram->bs;
  if (b < 0 || b >= (n+bs-1)/bs) {
    printf("cram_read_block: b = %ld #blocks = %ld\n", b, (n+bs-1)/bs);
    exit(1);
  }

  blk = darray_address(cram->da, b);
  compsize = block_len(&blk);
  mymalloc(compbuf, compsize);

  for (i=0; i<compsize; i++) compbuf[i] = *dblock_addr(&blk, i);

  origsize = bs;
  if (uncompress(buf, &origsize, compbuf, compsize) != Z_OK) {
    error("uncompress()");
  }

  free(compbuf);
}

static void crambg_write_block(CRAMBG *cram, i64 b, uchar *buf) // write buf to b-th block
{
  i64 i,j,b1,b2,n;
  i64 i1, i2, s, t;
  i64 bs;
  int k, r, c, len;
  dblock_entry blk;
  u64 x;
  i64 csize2;
  uchar *compbuf;
  unsigned long compsize, origsize;

  n = cram->n;
  bs = cram->bs;
  if (b < 0 || b >= (n+bs-1)/bs) {
    printf("cram_read_block: b = %ld #blocks = %ld\n", b, (n+bs-1)/bs);
    exit(1);
  }

  mymalloc(compbuf, bs*10);

  s = b*bs;  t = min((b+1)*bs-1, n-1);
  origsize = t - s + 1;
  compsize = bs*10;
  if (compress2(compbuf, &compsize, buf, origsize, 1) != Z_OK) {
    error("compress()");
  }
  darray_change(cram->da, b, compsize); // allocate a block
  blk = darray_address(cram->da, b);

  for (j=0; j<compsize; j++) {
    *dblock_addr(&blk, j) = compbuf[j];
  }

  free(compbuf);
}



void crambg_read(CRAMBG *cram, i64 s, i64 t, uchar *buf) // decode T[s..t] into buf
{
  i64 i,j,b1,b2,b;
  i64 i1, i2;
  i64 bs;
  int k, r, c;
  dblock_entry blk;
  uchar *tmpbuf;

  if (s < 0 || t >= cram->n) {
    printf("cram_read: s = %ld t = %ld n = %ld\n", s, t, cram->n);
    exit(1);
  }
  bs = cram->bs;

  mymalloc(tmpbuf, bs);

  b1 = s / bs;  b2 = t / bs;
  for (b = b1;  b <= b2; b++) { // decode each block
    crambg_read_block(cram, b, tmpbuf);
    i1 = 0;     if (b == b1) i1 = s % bs;
    i2 = bs-1;  if (b == b2) i2 = t % bs;
    for (i=i1; i<=i2; i++) *buf++ = tmpbuf[i];
  }
  free(tmpbuf);
}

void crambg_write(CRAMBG *cram, i64 s, i64 t, uchar *buf) // write buf into T[s..t]
{
  i64 i,j,b1,b2,b,n;
  i64 i1, i2;
  i64 bs;
  int k, r, c;
  dblock_entry blk;
  uchar *tmpbuf;
  i64 nb;

  if (s < 0 || t >= cram->n) {
    printf("cram_write: s = %ld t = %ld n = %ld\n", s, t, cram->n);
    exit(1);
  }
  bs = cram->bs;
  n = cram->n;

  mymalloc(tmpbuf, bs);

  b1 = s / bs;  b2 = t / bs;
  for (b = b1;  b <= b2; b++) { // change each block
    i1 = 0;     if (b == b1) i1 = s % bs;
    i2 = bs-1;  if (b == b2) i2 = t % bs;

    // read current block, which will be overwritten
    crambg_read_block(cram, b, tmpbuf);

    // rewrite
    for (i=i1; i<=i2; i++) tmpbuf[i] = *buf++;
    crambg_write_block(cram, b, tmpbuf);
  }

  free(tmpbuf);
}

#if 1

#include <sys/timeb.h>
#ifndef _MYTIMESTRUCT_
#define _MYTIMESTRUCT_
typedef struct timeb mytimestruct;
void mygettime(mytimestruct *t)
{
  ftime(t);
}
double mylaptime(mytimestruct *before,mytimestruct *after)
{
  double t;
  t = after->time - before->time;
  t += (double)(after->millitm - before->millitm)/1000;
  return t;
}
#endif

#define RANDOM 0


int main(int argc, char *argv[])
{
  CRAMBG *cram, *cram2;
  i64 s,t,i,n, n2, ss;
  uchar *buf;
  MMAP *map;
  uchar *text2;
  int b, bu;
  mytimestruct before,after;
  
  b = 1024;
//  if (argc>=4) b = atoi(argv[3]);
  fprintf(stderr, "block size = %d\n", b);

  mygettime(&before);
  cram = crambg_initialize(argv[1], b);
  mygettime(&after);
  fprintf(stderr, "initialize %f sec\n", mylaptime(&before,&after));

  n = cram->n;
  fprintf(stderr,"used memory %ld bytes (%1.2f bpc)\n", crambg_usedmemory(cram), (double)crambg_usedmemory(cram)*8/cram->n);

  buf = malloc(n);
#if 0
  // 圧縮されたメモリを読む時間を計る
  
  for (b = 16; b <= 1024; b <<= 1) {
    srand(0);
    mygettime(&before);
    for (s=0; s<n; s+=b) {
      if (s % ((n/b/100)*b) == 0) {
        fprintf(stderr, "%d %lf\n", s / ((n/b/100)*b), (double)crambg_usedmemory(cram)*8/n);
      }
#if RANDOM
      ss = rand() % (n-b);
#else
      ss = s;
#endif
      t = ss+b-1;  if (t >= n) t = n-1;
      crambg_read(cram, ss, t, buf);
    }
    mygettime(&after);
    fprintf(stdout, "read (%d-byte blocks) %f sec\n", b, mylaptime(&before,&after));
  }
#endif

#if 0
  // 圧縮されたメモリに同じものを上書きする時間を計る

  map = mymmap(argv[1]);
  text2 = map->addr;
  n2 = map->len;
  if (n2 > n) n2 = n;

  for (b = 16; b <= 1024; b <<= 1) {
    mygettime(&before);
    for (s=0; s<n; s+=b) {
      if (s % ((n2/b/100)*b) == 0) {
        fprintf(stderr, "%d %lf\n", s / ((n2/b/100)*b), (double)crambg_usedmemory(cram)*8/n);
      }
#if RANDOM
      ss = rand() % (n-b);
#else
      ss = s;
#endif
      t = ss+b-1;  if (t >= n) t = n-1;
      crambg_write(cram, ss, t, &text2[ss]);
#if 1
      crambg_read(cram, ss, t, &buf[ss]);
      for (i=ss; i<=t; i++) {
        if (buf[i] != text2[i] || 0) {
          printf("??? s = %ld i = %ld  %c  %c\n", ss, i, buf[i], text2[i]);
        }
      }
#endif
    }
    mygettime(&after);
    fprintf(stdout, "write (%d-byte blocks) %f sec\n", b, mylaptime(&before,&after));
  }
#endif

//  fprintf(stderr,"used memory %ld bytes (%1.2f bpc)\n", cram_usedmemory(cram), (double)cram_usedmemory(cram)*8/n);

#if 1
  // 圧縮されたメモリに違うものを上書きする時間を計る
  map = mymmap(argv[2]);
  text2 = map->addr;
  n2 = map->len;
  if (n2 > n) n2 = n;


  b = 1024;
  if (argc>=4) b = atoi(argv[3]);
  mygettime(&before);
  for (s=0; s<n2; s+=b) {
    if (s % ((n2/b/100)*b) == 0) {
      fprintf(stderr, "%ld \r", s / ((n2/b/100)*b));  fflush(stderr);
      printf("%d %lf\n", s / ((n2/b/100)*b), (double)crambg_usedmemory(cram)*8/n);
    }
    t = s+b-1;  if (t >= n) t = n-1;
    crambg_write(cram, s, t, &text2[s]);
#if 0
    crambg_read(cram, s, t, &buf[s]);
    for (i=s; i<=t; i++) {
      if (buf[i] != text2[i]) {
        printf("??? s = %ld i = %ld %c %c\n", s, i, buf[i], text2[i]);
      }
    }
#endif
  }
  mygettime(&after);
  fprintf(stderr, "rewrite %f sec b = %d\n", mylaptime(&before,&after), b);
  fprintf(stderr, "used memory %ld bytes (%1.2f bpc)\n", crambg_usedmemory(cram), (double)crambg_usedmemory(cram)*8/cram->n);
  printf("%d %lf\n", 100, (double)crambg_usedmemory(cram)*8/n);
#endif

#if 0
  // ディスクに正しく保存されているか調べる
  crambg_save(cram, "crambg1.dat", "crambg2.dat");
  crambg_free(cram);
  cram2 = crambg_load("crambg1.dat", "crambg2.dat");

  mygettime(&before);
  for (s=0; s<n2; s+=b) {
    if (s % (b*100) == 0) {
      fprintf(stderr, "%ld/%ld \r", s, n2);  fflush(stderr);
    }
    t = s+b-1;  if (t >= n2) t = n2-1;
    crambg_read(cram2, s, t, buf);
    //printf("block %ld: ", s/128);
    for (i=s; i<=t; i++) {
      if (buf[i-s] != text2[i]) {
        printf("??? s = %ld i = %ld\n", s, i);
      }
    }
//    for (i=0; i<t-s+1; i++) printf("%c", buf[i]);
    //printf("\n");
  }
  mygettime(&after);
  fprintf(stderr, "read %f sec\n", mylaptime(&before,&after));
#endif
  return 0;
}
#endif
