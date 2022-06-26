#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cram.h"
#include "mman.h"

#define ID_CRAM     0x0a


#define mymalloc(p,n) {p = malloc((n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory in line %d\n",__LINE__); exit(1);};}

#ifndef min
 #define min(x,y) (((x)<(y))?(x):(y))
#endif
#ifndef max
 #define max(x,y) (((x)>(y))?(x):(y))
#endif

#define D (8*sizeof(uchar))

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

#if 0
static u64 setbits0(dblock_entry *blk, i64 i, int d, u64 x)
{
  int j;

  for (j=0; j<d; j++) {
    setbit(blk,i+j,(x>>(d-j-1))&1);
  }
  return x;
}
#endif

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

static u64 getbits(dblock_entry *blk, i64 i, int d)
{
  u64 x;
  if (d > 32) {
    printf("getbits: d = %d\n", d);
    exit(1);
  }
  x = getbit32(blk, i);
  x >>= (32 - d);
  return x;
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


static void compute_rank(int sigma, rank_code *code)
{
  int i;
  for (i = 0; i < sigma; i++) code[i].rank = 0;
}

CRAM *cram_initialize(char *filename, int bs, int sbs, int bu, u16 opt)
{
  CRAM *cram;
  MMAP *map;
  i64 i,j,n,nb;
  int k, sigma, c, c2, nc;
  i64 *freq;
  uchar *text;
  Huffman2 *huf;
//  double *f2;
  i64 csize, csize2;
  darray *da;
  dblock_entry blk;
  int len, r, r2, lsb;
  u64 x;
  
  uchar *tmpbuf;
  tmpbuf = malloc(bs);


  mymalloc(cram,1);

  cram->opt = opt;

  cram->k = k = 2;

  cram->bs = bs;
  cram->sbs = sbs/k; // number of super-characters in a small block

  map = mymmap(filename);
  cram->n = n = map->len;
  text = map->addr;

  cram->logn = (blog(n)+1+7)/8;

  nb = (n+bs-1)/bs; // number of blocks

  cram->phase = 1;
//  cram->step = 0;
  cram->step = nb*(bu-1)/(bu);
  cram->bu = bu;

  cram->counter = 0;
  cram->nbytes = 0;

//  cram->divergence = 0.0;

  cram->da = da = darray_initialize(nb, 1, bs*2, n/4, n/100);


  cram->sigma = sigma = 1 << (k*8);

  nc = (cram->opt & OPT_DEAMORTIZE) ? 3 : 2;

#if USE_HUFFMAN

  for (i=0; i<nc; i++) {
    mymalloc(cram->freq[i], sigma);
  }
  freq = cram->freq[0];
  for (c=0; c<sigma; c++) freq[c] = 1;

  for (i=0; i<(n+k-1)/k; i++) {
    if (i % (n/k/100) == 0) {
      fprintf(stderr, "counting freq %ld \r", i / (n/k/100));  fflush(stderr);
    }
    c = getuint(text, i*k, n, k);
    freq[c]++;
  }

  for (i=1; i<nc; i++) {
    for (c=0; c<sigma; c++) {
      cram->freq[i][c] = freq[c];
    }
  }
  for (i=0; i<nc; i++) {
    cram->huf[i] = huf = MakeHuffman2Tree2(sigma, freq, 8*k);
  }
#else

  for (i=0; i<nc; i++) {
    mymalloc(cram->code[i], sigma);
  }

  for (c=0; c<sigma; c++) cram->code[0][c].freq = 1;

  for (i=0; i<(n+k-1)/k; i++) {
    c = getuint(text, i*k, n, k);
    cram->code[0][c].freq++;
  }

  for (i=1; i<nc; i++) {
    for (c=0; c<sigma; c++) {
      cram->code[i][c].freq = cram->code[0][c].freq;
    }
  }

  for (i=0; i<nc; i++) {
    compute_rank(sigma, cram->code[i]);
  }
#endif

  if (bs % k != 0) {
    printf("cram_initialize: bs = %d k = %d\n", bs, k);
    exit(1);
  }

  csize = 0;
  for (i=0; i<nb; i++) { // compress i-th block
    if (i % (nb/100) == 0) {
      fprintf(stderr, "encoding %ld \r", i / (nb/100));  fflush(stderr);
    }
//    if (i % 1000 == 0) {
//      fprintf(stderr, "%ld/%ld \r", i, nb);  fflush(stderr);
//    }
    csize2 = 0;
    for (j=0; j<bs; j+=k) {
      c = getuint(text, i*bs+j, n, k);
      csize += huf->clen[c];
      csize2 += huf->clen[c];
    }
    if (cram->sbs > 0) {
      lsb = blog(csize2-1)+1;
      while (1) {
        int b;
        b = (csize2+ lsb * (bs/k/cram->sbs) + 7)/8;
        if (blog(b*8-1)+1 == lsb) break;
        lsb++;
        if (lsb > 20) {
          printf("??? c = %ld\n", lsb);
        }
      }
      csize2 += lsb * (bs/k/cram->sbs);
      csize  += lsb * (bs/k/cram->sbs);
    }
    darray_change(da, i, (csize2+7)/8); // allocate a block
    blk = darray_address(da, i);

    r = 0; // bit address
    if (cram->sbs > 0) {
      r2 = lsb * (bs/k/cram->sbs);
    } else {
      r2 = 0;
    }
    for (j=0; j<bs; j+=k) {
      c = getuint(text, i*bs+j, n, k);
      len = huf->clen[c];
      if (len == 0) {
        printf("clen[%d] = %d\n", c, len);
        exit(1);
      }
      if (cram->sbs > 0) {
        if ((j/k) % cram->sbs == 0) {
          setbits(&blk, (j/k/cram->sbs)*lsb, lsb, r);
        }
      }
      setbits(&blk, r2 + r, len, huf->code[c] >> (64-len));
      r += len;
    }
#if 0
{
  i64 i2;
  i2 = i + bs - 1;
  if (i2 > n-1) i2 = n-1;
  cram_read(cram, i, i2, tmpbuf);
  for (j=i; j<=i2; j++) {
    if (tmpbuf[j-i] != text[j]) {
      printf("error i = %ld j = %ld tmp = [%c] text = [%c]\n", i, j, tmpbuf[j-i], text[j]);
    }
  }
}
#endif
  }

  fprintf(stderr,"n = %ld csize = %ld (%1.2f bpc)\n", n, csize, (double)csize/n);
  fprintf(stderr,"darray used memory %ld bytes (%1.2f bpc)\n", darray_usedmemory(cram->da), (double)darray_usedmemory(cram->da)*8/n);

//  darray_check(da);

  return cram;
}

void cram_free(CRAM *cram)
{
  int i, nc;

  nc = (cram->opt & OPT_DEAMORTIZE) ? 3 : 2;
  for (i=0; i<nc; i++) {
#if USE_HUFFMAN
    free(cram->freq[i]);
    freeHuffman2(cram->huf[i]);
#else
#endif
  }
  darray_free(cram->da);
  free(cram);
}

i64 cram_usedmemory(CRAM *cram)
{
  i64 size;
  size = sizeof(CRAM);
  size += darray_usedmemory(cram->da);
  size += cram->sigma * sizeof(*cram->freq[0]) * 3;
  if (cram->opt & OPT_DEAMORTIZE) {
    size += Huffman2_usedmemory(cram->huf[0]) * 3;
  } else {
    size += Huffman2_usedmemory(cram->huf[0]) * 2;
  }
  return size;
}

static Huffman2 *current_huf(CRAM *cram, i64 b)
{
  Huffman2 *huf;
  if (cram->opt & OPT_DEAMORTIZE) {
    if (b < cram->step) { // already rewritten with new code
      huf = cram->huf[(cram->phase-1) % 3];
    } else { // encoded with old code
      huf = cram->huf[(cram->phase+1) % 3];
    }
  } else {
    if (b < cram->step) { // already rewritten with new code
      huf = cram->huf[(cram->phase-1) % 2];
    } else { // encoded with old code
      huf = cram->huf[cram->phase % 2];
    }
  }
  return huf;
}

static void cram_read_block(CRAM *cram, i64 b, uchar *buf) // decode b-th block into buf
{
  i64 i,j,b1,b2,n;
  i64 i1, i2;
  i64 bs;
  int k, r, c;
  dblock_entry blk;
  Huffman2 *huf;
  u64 x;
  int lsb;

  n = cram->n;
  bs = cram->bs;
  if (b < 0 || b >= (n+bs-1)/bs) {
    printf("cram_read_block: b = %ld #blocks = %ld\n", b, (n+bs-1)/bs);
    exit(1);
  }
  k = cram->k;

  huf = current_huf(cram, b);

  blk = darray_address(cram->da, b);
  i2 = bs;  if (b*i2 >= n) i2 = n % bs;
  if (cram->sbs > 0) {
    lsb = blog(block_len(&blk)*8-1)+1;
    r = lsb * (bs/k/cram->sbs);
  } else {
    r = 0;
  }
  for (i=0; i<i2; i+=k) {
    x = getbit32(&blk, r) << 32;
//    c = DecodeHuffman2(huf, x);
    c = DecodeHuffman2_tbl(huf, x);
    for (j=0; j<k; j++) {
      *buf++ = c >> (8*(k-1-j));
    }
    r += huf->clen[c];
  }
}

static void cram_read_sb(CRAM *cram, i64 b, i64 s, i64 t, uchar *buf) // decode a part [s,t] of b-th block into buf
{
  i64 i,j,b1,b2,n;
  i64 i1, i2;
  i64 bs;
  int k, r, c;
  dblock_entry blk;
  Huffman2 *huf;
  u64 x;
  int lsb;
  int sb, tb;

  n = cram->n;
  bs = cram->bs;
  if (b < 0 || b >= (n+bs-1)/bs) {
    printf("cram_read_block: b = %ld #blocks = %ld\n", b, (n+bs-1)/bs);
    exit(1);
  }
  k = cram->k;

  huf = current_huf(cram, b);

  blk = darray_address(cram->da, b);

  sb = s / k / cram->sbs;
//  tb = ((t+k-1) / k  + cram->sbs - 1)/ cram->sbs;

  lsb = blog(block_len(&blk)*8-1)+1;
  r = getbit32(&blk, sb * lsb) >> (32-lsb);  // offset
  r += lsb * (bs/k/cram->sbs);

  i1 = sb * cram->sbs * k;
//  i2 = tb * cram->sbs * k;
  i2 = t;
  for (i=i1; i<=i2; i+=k) {
    x = getbit32(&blk, r) << 32;
//    c = DecodeHuffman2(huf, x);
    c = DecodeHuffman2_tbl(huf, x);
    for (j=0; j<k; j++) {
      if (s <= i+j && i+j <= t) {
        *buf++ = c >> (8*(k-1-j));
      }
    }
    r += huf->clen[c];
  }
}

static void cram_write_block(CRAM *cram, i64 b, uchar *buf) // write buf to b-th block
{
  i64 i,j,b1,b2,n;
  i64 i1, i2;
  i64 bs;
  int k, r, c, len, r2;
  dblock_entry blk;
  Huffman2 *huf;
  u64 x;
  i64 csize2;
  int lsb;

  n = cram->n;
  bs = cram->bs;
  if (b < 0 || b >= (n+bs-1)/bs) {
    printf("cram_write_block: b = %ld #blocks = %ld\n", b, (n+bs-1)/bs);
    exit(1);
  }
  k = cram->k;

  huf = current_huf(cram, b);

  i2 = bs;  if (b*i2 >= n) i2 = n % bs;
  csize2 = 0;
  for (i=0; i<i2; i+=k) {
    c = getuint(buf, i, i2, k);
    csize2 += huf->clen[c];
  }

  if (cram->sbs > 0) {
    lsb = blog(csize2-1)+1;
    while (1) {
      int b;
      b = (csize2+ lsb * (bs/k/cram->sbs) + 7)/8;
      if (blog(b*8-1)+1 == lsb) break;
      lsb++;
      if (lsb > 20) {
        printf("??? c = %ld\n", lsb);
      }
    }
    csize2 += lsb * (bs/k/cram->sbs);
  }

  darray_change(cram->da, b, (csize2+7)/8); // change the block size
  blk = darray_address(cram->da, b);

  if (cram->sbs > 0) {
    r2 = lsb * (bs/k/cram->sbs);
  } else {
    r2 = 0;
  }
  r = 0; // bit address
  for (i=0; i<i2; i+=k) {
    c = getuint(buf, i, i2, k);
    len = huf->clen[c];
    if (len == 0) {
      printf("clen[%d] = %d\n", c, len);
      exit(1);
    }
    if (cram->sbs > 0) {
      if ((i/k) % cram->sbs == 0) {
        setbits(&blk, (i/k/cram->sbs)*lsb, lsb, r);
      }
    }
    setbits(&blk, r2 + r, len, huf->code[c] >> (64-len));
    r += len;
  }
}

static void cram_write_sb(CRAM *cram, i64 b, i64 s, i64 t, uchar *buf) // write buf to a part [s,t] of b-th block
{
  i64 i,j,b1,b2,sb1, sb2, n;
  i64 i1, i2;
  i64 bs;
  int k, r, c, len, r2, r0;
  dblock_entry blk, tmpblk;
  Huffman2 *huf;
  u64 x;
  i64 csize2;
  int lsb, lsb0;
  int sbs;
  uchar *tmpbuf;
  int p1, p2, p3, p, ofs;

  n = cram->n;
  bs = cram->bs;
  if (b < 0 || b >= (n+bs-1)/bs) {
    printf("cram_write_sb: b = %ld #blocks = %ld\n", b, (n+bs-1)/bs);
    exit(1);
  }
  k = cram->k;
  sbs = cram->sbs;

  huf = current_huf(cram, b);

  b1 = s / k;  b2 = t / k; // first and last super-characters
  sb1 = b1 / sbs;  sb2 = b2 / sbs; // first and last small blocks

  blk = darray_address(cram->da, b);
  r = block_len(&blk);
  mymalloc(tmpbuf, r);
  for (i=0; i<r; i++) tmpbuf[i] = *dblock_addr(&blk, i);
  tmpblk.addr[0] = tmpbuf;  tmpblk.addr[1] = NULL;
  tmpblk.len[0] = r;  tmpblk.len[1] = 0;

  lsb0 = blog(r*8-1)+1;
  p1 = getbit32(&tmpblk, sb1 * lsb0) >> (32-lsb0); // first bit of the small block to be changed
  if (sb2+1 < bs/k/sbs) {
    p2 = getbit32(&tmpblk, (sb2+1) * lsb0) >> (32-lsb0); // first bit of the next of last small block to be changed
  } else {
    p2 = r*8 - lsb0*(bs/k/sbs);
  }
  p3 = r*8 - lsb0*(bs/k/sbs);

//  for (i=0; i<bs/k/sbs; i++) {
//    printf("sb[%d] %d\n", i, getbit32(&tmpblk, i * lsb0) >> (32-lsb0));
//  }



  csize2 = 0;
  for (i=0; i<=(sb2-sb1+1)*sbs-1; i++) {
    c = getuint(buf, i*k, t+1, k);
    csize2 += huf->clen[c];
  }
  csize2 += p1;
  if (sb2+1 < bs/k/sbs) {
    ofs = csize2 - p2;
  } else {
    ofs = 0; //?
  }
  csize2 += p3 - p2;
//  printf("ofs = %d\n", ofs);

  lsb = blog(csize2-1)+1;
  while (1) {
    int b;
    b = (csize2+ lsb * (bs/k/sbs) + 7)/8;
    if (blog(b*8-1)+1 == lsb) break;
    lsb++;
    if (lsb > 20) {
      printf("??? c = %ld\n", lsb);
    }
  }
  csize2 += lsb * (bs/k/sbs);

  darray_change(cram->da, b, (csize2+7)/8); // change the block size
  blk = darray_address(cram->da, b);

  r2 = lsb * (bs/k/sbs);
  // copy data
  r0 = lsb0 * (bs/k/sbs);
  for (p=0; p+32<=p1; p+=32) {
    x = getbit32(&tmpblk, r0 + p);
    setbits(&blk, r2 + p, 32, x);
  }
  if (p < p1) {
    x = getbit32(&tmpblk, r0 + p) >> (32-(p1-p));
    setbits(&blk, r2 + p, p1-p, x);
  }
  // copy pointers
  for (i=0; i<sb1; i++) {
    x = getbit32(&tmpblk, i*lsb0) >> (32-lsb0);
    setbits(&blk, i*lsb, lsb, x);
  }

  r = p1; // bit address
  for (i=0; i<=(sb2-sb1+1)*sbs-1; i++) {
    c = getuint(buf, i*k, bs, k);
    len = huf->clen[c];
    if (len == 0) {
      printf("clen[%d] = %d\n", c, len);
      exit(1);
    }
    if (i % sbs == 0) {
      setbits(&blk, (i/sbs+sb1)*lsb, lsb, r);
    }
    setbits(&blk, r2 + r, len, huf->code[c] >> (64-len));
    r += len;
  }

  // copy data
  r0 = lsb0 * (bs/k/sbs);
  for (p=p2; p+32<=p3; p+=32) {
    x = getbit32(&tmpblk, r0 + p);
    setbits(&blk, r2 + r, 32, x);
    r += 32;
  }
  if (p < p3) {
    x = getbit32(&tmpblk, r0 + p) >> (32-(p3-p));
    setbits(&blk, r2 + r, p3-p, x);
  }
  // copy pointers
  for (i=sb2+1; i<bs/k/sbs; i++) {
    x = getbit32(&tmpblk, i*lsb0) >> (32-lsb0);
    setbits(&blk, i*lsb, lsb, x + ofs);
  }

//  for (i=0; i<bs/k/sbs; i++) {
//    printf("sb[%d] %d\n", i, getbit32(&blk, i * lsb0) >> (32-lsb0));
//  }


  free(tmpbuf);
}



void cram_read(CRAM *cram, i64 s, i64 t, uchar *buf) // decode T[s..t] into buf
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
  k = cram->k;

  mymalloc(tmpbuf, bs);

  b1 = s / bs;  b2 = t / bs;
  for (b = b1;  b <= b2; b++) { // decode each block
    i1 = 0;     if (b == b1) i1 = s % bs;
    i2 = bs-1;  if (b == b2) i2 = t % bs;
    if (i1 == 0 && i2 == bs-1) {
      cram_read_block(cram, b, tmpbuf);
      for (i=i1; i<=i2; i++) *buf++ = tmpbuf[i];
    } else {
      if (cram->sbs > 0) {
        cram_read_sb(cram, b, i1, i2, tmpbuf);
        for (i=i1; i<=i2; i++) *buf++ = tmpbuf[i-i1];
      } else {
        cram_read_block(cram, b, tmpbuf);
        for (i=i1; i<=i2; i++) *buf++ = tmpbuf[i];
      }
    }
  }
  free(tmpbuf);
}

void cram_write(CRAM *cram, i64 s, i64 t, uchar *buf) // write buf into T[s..t]
{
  i64 i,j,b1,b2,b,n, sb1, sb2;
  i64 i1, i2;
  i64 bs;
  int k, r, c;
  dblock_entry blk;
  uchar *tmpbuf;
  i64 nb;
  double *f2;
  double f;
  int sbs;

  if (s < 0 || t >= cram->n) {
    printf("cram_write: s = %ld t = %ld n = %ld\n", s, t, cram->n);
    exit(1);
  }
  bs = cram->bs;
  k = cram->k;
  n = cram->n/k;
  sbs = cram->sbs;

  mymalloc(tmpbuf, bs);

  b1 = s / bs;  b2 = t / bs;
  for (b = b1;  b <= b2; b++) { // change each block
    i1 = 0;     if (b == b1) i1 = s % bs;
    i2 = bs-1;  if (b == b2) i2 = t % bs;
    if ((i1 == 0 && i2 == bs-1) || (sbs == 0) || 0) {
      // read current block, which will be overwritten
      cram_read_block(cram, b, tmpbuf);
      // update frequency table
      for (i=0; i<(bs+k-1)/k; i++) {
        c = getuint(tmpbuf, i*k, bs, k);
        if (cram->opt & OPT_DEAMORTIZE) {
          cram->freq[cram->phase % 3][c]--;
        } else {
          cram->freq[cram->phase % 2][c]--;
        }
      }

      // rewrite
      for (i=i1; i<=i2; i++) tmpbuf[i] = *buf++;
      // update frequency table
      for (i=0; i<(bs+k-1)/k; i++) {
        c = getuint(tmpbuf, i*k, bs, k);
        if (cram->opt & OPT_DEAMORTIZE) {
          cram->freq[cram->phase % 3][c]++;
        } else {
          cram->freq[cram->phase % 2][c]++;
        }
      }
      cram_write_block(cram, b, tmpbuf);
    } else { // change a part of a block
      i64 b1, b2;
      b1 = i1 / k;  b2 = i2 / k; // first and last super-characters to change
      sb1 = b1 / sbs;  sb2 = b2 / sbs;
      cram_read_sb(cram, b, sb1*sbs*k, (sb2+1)*sbs*k-1, &tmpbuf[sb1*sbs*k]);
      // update frequency table
      for (i=b1; i<=b2; i++) {
        c = getuint(tmpbuf, i*k, bs, k);
        if (cram->opt & OPT_DEAMORTIZE) {
          cram->freq[cram->phase % 3][c]--;
        } else {
          cram->freq[cram->phase % 2][c]--;
        }
      }

      // rewrite
      for (i=i1; i<=i2; i++) tmpbuf[i] = *buf++;
      // update frequency table
      for (i=b1; i<=b2; i++) {
        c = getuint(tmpbuf, i*k, bs, k);
        if (cram->opt & OPT_DEAMORTIZE) {
          cram->freq[cram->phase % 3][c]++;
        } else {
          cram->freq[cram->phase % 2][c]++;
        }
      }
      cram_write_sb(cram, b, sb1*sbs*k, (sb2+1)*sbs*k-1, &tmpbuf[sb1*sbs*k]);
    }
  }

  // re-encode
  nb = (cram->n+cram->bs-1)/cram->bs; // number of blocks

  cram->nbytes += (t-s+1) * cram->bu;
  while (cram->nbytes >= cram->bs && cram->step < nb) {
    b = cram->step;
    cram_read_block(cram, b, tmpbuf);
    cram->step += 1; // to cheat (to use new code) the write_block routine ???
    cram_write_block(cram, b, tmpbuf);
    cram->step -= 1;
    if (cram->opt & OPT_DEAMORTIZE) {
      if (cram->counter < cram->sigma) {
        cram->freq[(cram->phase-1) % 3][cram->counter] = cram->freq[(cram->phase+1) % 3][cram->counter];
        cram->counter++;
      }
    }
    cram->step += 1;
    cram->nbytes -= cram->bs;
  }

  if (cram->step >= nb) { // end of phase
//    printf("end of phase %ld              \n", cram->phase);
    if (cram->opt & OPT_DEAMORTIZE) {
      if (cram->counter < cram->sigma) {
//      fprintf(stderr, "counter = %ld sigma = %ld\n", cram->counter, cram->sigma);
        for (c=cram->counter; c<cram->sigma; c++) {
          cram->freq[(cram->phase-1) % 3][c] = cram->freq[(cram->phase+1) % 3][c];
        }
      }
      cram->counter = 0;

      freeHuffman2(cram->huf[cram->phase % 3]);
      cram->huf[cram->phase % 3] = MakeHuffman2Tree2(cram->sigma, cram->freq[(cram->phase-1) % 3], 8*k);
    } else {
      for (c=0; c<cram->sigma; c++) {
        cram->freq[(cram->phase+1) % 2][c] = cram->freq[cram->phase % 2][c];
      }
      freeHuffman2(cram->huf[cram->phase % 2]);
      cram->huf[cram->phase % 2] = MakeHuffman2Tree2(cram->sigma, cram->freq[(cram->phase) % 2], 8*k);
    }
    cram->step = 0;
    cram->phase++;
  }

  free(tmpbuf);
}

void cram_save(CRAM *cram, uchar *filename, uchar *memfilename)
{
  int i, j, c, k;
  FILE *out, *mout;

  out = fopen(filename, "w");
  mout = fopen(memfilename, "w");
  if (out == NULL || mout == NULL) {
    perror("cram_save: fopen\n");
    exit(1);
  }

  writeuint(1, ID_CRAM, out);
  k = cram->logn;
  writeuint(1, k, out);
  writeuint(k, cram->n, out);
  writeuint(k, cram->bs, out);
  writeuint(k, cram->sbs, out);
  writeuint(k, cram->bu, out);

  writeuint(k, cram->phase, out);
  writeuint(k, cram->step, out);
  writeuint(k, cram->counter, out);
  writeuint(k, cram->nbytes, out);

  writeuint(sizeof(cram->opt), cram->opt, out);

  writeuint(k, cram->k, out);
  writeuint(k, cram->sigma, out);

  c = (cram->opt & OPT_DEAMORTIZE) ? 3 : 2;
  for (i=0; i<c; i++) {
    for (j=0; j<cram->sigma; j++) {
      writeuint(k, cram->freq[i][j], out);
    }
  }
  for (i=0; i<c; i++) {
    Huffman2_write(cram->huf[i], out);
  }

  darray_write(cram->da, out, mout);

  fclose(out);
  fclose(mout);
}

CRAM *cram_load(uchar *filename, uchar *memfilename)
{
  CRAM *cram;
  int k,i,j,c, id;
  MMAP *map;
  uchar *p;
  
  mymalloc(cram, 1);

  map = mymmap(filename);
  p = map->addr;
  
  if ((id = readuint(p,1)) != ID_CRAM) {
    printf("cram_load: id = %d\n",id);
    exit(1);
  }
  p += 1;
  cram->logn = k = readuint(p,1);  p += 1;
  cram->n = readuint(p,k);  p += k;
  cram->bs = readuint(p,k);  p += k;
  cram->sbs = readuint(p,k);  p += k;
  cram->bu = readuint(p,k);  p += k;

  cram->phase = readuint(p,k);  p += k;
  cram->step = readuint(p,k);  p += k;
  cram->counter = readuint(p,k);  p += k;
  cram->nbytes = readuint(p,k);  p += k;
  

  cram->opt = readuint(p,sizeof(cram->opt));  p += sizeof(cram->opt);

  cram->k = readuint(p,k);  p += k;
  cram->sigma = readuint(p,k);  p += k;

  c = (cram->opt & OPT_DEAMORTIZE) ? 3 : 2;
  for (i=0; i<c; i++) {
    mymalloc(cram->freq[i], cram->sigma);
    for (j=0; j<cram->sigma; j++) {
      cram->freq[i][j] = readuint(p,k);  p += k;
    }
  }
  for (i=0; i<c; i++) {
    cram->huf[i] = Huffman2_read(&p);
    cram->huf[i]->tbl_width = 8*cram->k;
    mkhufdecodetable2(cram->huf[i]);
  }

  cram->da = darray_read(&p, memfilename);

  return cram;
}

double cram_entropy(CRAM *cram)
{
  double e, p;
  i64 *freq;
  i64 c, n;

  freq = cram->freq[cram->phase % 3];
  n = cram->n / cram->k;
  e = 0.0;
  for (c=0; c<cram->sigma; c++) {
    p = (double)freq[c] / n;
    e += p * log(1/p);
  }
  e /= cram->k;
  e /= log(2.0);
  return e;
}

double cram_divergence(CRAM *cram)
{
  double e, p, q;
  i64 *freq, *freq2;
  i64 c, n;
#if 1
  freq = cram->freq[cram->phase % 3];
  freq2 = cram->freq[(cram->phase+1) % 3];
  n = cram->n / cram->k;
  e = 0.0;
  for (c=0; c<cram->sigma; c++) {
    p = (double)freq[c] / n;
    q = (double)freq2[c] / n;
    e += p * log(p/q);
  }
  e /= cram->k;
  e /= log(2.0);
  return e;
#else
  return cram->divergence / cram->k / log(2.0);
#endif
}


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



#ifdef MKIDX_MAIN

int main(int argc, char *argv[])
{
  CRAM *cram, *cram2;
  i64 s,t,i,n, n2, ss;
  uchar *buf;
  MMAP *map;
  uchar *text2;
  int b, bu, sbs;
  mytimestruct before,after;
  
  b = 1024;
  sbs = 64;
  bu = 4;
  if (argc>=4) b = atoi(argv[2]);
  if (argc>=5) sbs = atoi(argv[3]);
  if (argc>=6) bu = atoi(argv[4]);
  fprintf(stderr, "block size = %d small-block size = %d update = %d\n", b, sbs, bu);

  mygettime(&before);
  cram = cram_initialize(argv[1], b, sbs, bu, 0);
  mygettime(&after);
  fprintf(stderr, "initialize %f sec\n", mylaptime(&before,&after));

  cram_save(cram, "cram1.dat", "cram2.dat");

  n = cram->n;
  fprintf(stderr,"used memory %ld bytes (%1.2f bpc)\n", cram_usedmemory(cram), (double)cram_usedmemory(cram)*8/cram->n);
}
#endif

#ifdef TEST_MAIN

#include <string.h>


#define RANDOM 0

#if RANDOM
unsigned long long genrand64_int64(void);
void init_by_array64(unsigned long long init_key[],
		     unsigned long long key_length);
#endif

int main(int argc, char *argv[])
{
  CRAM *cram, *cram2;
  i64 s,t,i,n, n2, ss;
  uchar *buf;
  MMAP *map;
  uchar *text2;
  int b, bu, sbs;
  mytimestruct before,after;
  unsigned long long init[4]={0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL};

  
#if 1
  // 圧縮されたデータを b バイトずつ読み込む時間を測定
  // (正しく復元されたかどうか text2 と比較)
  cram = cram_load(argv[2], argv[3]);
  n = cram->n;

  map = mymmap(argv[1]);
  text2 = map->addr;
  mymalloc(buf, 1024);

//  srand(0);
#if RANDOM
  init_by_array64(init, 4);
#endif
  for (b = 1024; b >= 16; b >>= 1) {
//    srand(0);
    mygettime(&before);
    for (s=0; s<n; s+=b) {
      if (s % ((n/b/100)*b) == 0) {
        fprintf(stderr, "%d\r", s / ((n/b/100)*b));
      }
#if RANDOM
//      ss = rand() % (n-b);
      ss = genrand64_int64() % (n-b);
#else
      ss = s;
#endif
      t = ss+b-1;  if (t >= n) t = n-1;
      cram_read(cram, ss, t, buf);
#if 0
      for (i=ss; i<=t; i++) {
        if (buf[i-ss] != text2[i] || 0) {
          printf("??? s = %ld i = %ld  %c  %c\n", ss, i, buf[i-ss], text2[i]);
        }
      }
#endif
    }
    mygettime(&after);
    fprintf(stdout, "read (%d-byte blocks) %f sec\n", b, mylaptime(&before,&after));
  }
  cram_free(cram);
#endif

#if 0
  // 圧縮されていないデータを読む時間を測定
  
  map = mymmap(argv[1]);
  text2 = map->addr;
//  n2 = map->len;
//  if (n2 > n) n2 = n;
  n = map->len;
  mymalloc(buf, 1024);

//  srand(0);
#if RANDOM
  init_by_array64(init, 4);
#endif
  for (b = 1024; b >= 16; b >>= 1) {
//    srand(0);
    mygettime(&before);
    for (s=0; s<n; s+=b) {
      if (s % ((n/b/100)*b) == 0) {
        fprintf(stderr, "%d\r", s / ((n/b/100)*b));
      }
#if RANDOM
//      ss = rand() % (n-b);
      ss = genrand64_int64() % (n-b);

#else
      ss = s;
#endif
      t = ss+b-1;  if (t >= n) t = n-1;
      memcpy(buf, &text2[ss], t-ss+1);
    }
    mygettime(&after);
    fprintf(stdout, "read memory (%d-byte blocks) %f sec\n", b,
                     mylaptime(&before,&after));
  }
  mymunmap(map);
#endif



#if 0
  // text2 を圧縮する時間を測定
  // (正しく圧縮できたか確認)
  cram = cram_load(argv[2], argv[3]);
  n = cram->n;

  map = mymmap(argv[1]);
  text2 = map->addr;
  mymalloc(buf, 1024);

  n2 = map->len;
  if (n2 > n) n2 = n;

  for (b = 16; b <= 1024; b <<= 1) {
    mygettime(&before);
    for (s=0; s<n; s+=b) {
      if (s % ((n2/b/100)*b) == 0) {
        fprintf(stderr, "%d %lf\n", s / ((n2/b/100)*b), (double)cram_usedmemory(cram)*8/n);
      }
#if RANDOM
      ss = rand() % (n-b);
#else
      ss = s;
#endif
      t = ss+b-1;  if (t >= n) t = n-1;
//      printf("write %ld\n", s);
      cram_write(cram, ss, t, &text2[ss]);
#if 1
      cram_read(cram, s*1, t, &buf[s*1]);
      for (i=s*1; i<=t; i++) {
        if (buf[i] != text2[i] || 0) {
          printf("??? s = %ld i = %ld  %c  %c\n", s, i, buf[i], text2[i]);
        }
      }
#endif
    }
    mygettime(&after);
    fprintf(stderr, "write (%d-byte blocks) %f sec\n", b, mylaptime(&before,&after));
  }
#endif

#if 0
  // データを圧縮し，別のデータを上書きしていく時間を測定
  // (正しく圧縮・復元できるか確認)
  b = 1024;
  sbs = 64;
  bu = 4;

  cram = cram_initialize(argv[1], b, sbs, bu, 0);
  n = cram->n;
  fprintf(stderr,"used memory %ld bytes (%1.2f bpc)\n", cram_usedmemory(cram), (double)cram_usedmemory(cram)*8/cram->n);

  map = mymmap(argv[2]);
  text2 = map->addr;
  n2 = map->len;
  if (n2 > n) n2 = n;
  mymalloc(buf, 1024);

  b = 1024;
  if (argc>=4) b = atoi(argv[3]);
  mygettime(&before);
  for (s=0; s<n2; s+=b) {
    if (s % ((n2/b/100)*b) == 0) {
      fprintf(stderr, "%ld \r", s / ((n2/b/100)*b));  fflush(stderr);
      printf("%d %lf\n", s / ((n2/b/100)*b), (double)cram_usedmemory(cram)*8/n);
//      printf("%d %lf\n", s / ((n2/b/100)*b), cram_entropy(cram));
//      printf("%d %lf\n", s / ((n2/b/100)*b), cram_divergence(cram));
    }
    t = s+b-1;  if (t >= n) t = n-1;
    cram_write(cram, s, t, &text2[s]);
#if 0
    for (i=s; i<=t; i+=b) {
      i64 j;
      cram_read(cram, i, i+b-1, &buf[i-s]);
//      for (j=i; j<i+b; j++) printf("%c", buf[j]);
//      printf("\n");
    }
//    cram_read(cram, s*1, t, buf);
//    for (i=0; i<t-s+1; i++) printf("%c", buf[i]);
    for (i=s; i<=t; i++) {
      if (buf[i-s] != text2[i]) {
        printf("??? s = %ld i = %ld [%c] [%c]\n", s, i, buf[i], text2[i]);
      }
    }
#endif
  }
  mygettime(&after);
  fprintf(stderr, "rewrite %f sec b = %d\n", mylaptime(&before,&after), b);
  fprintf(stderr, "used memory %ld bytes (%1.2f bpc)\n", cram_usedmemory(cram), (double)cram_usedmemory(cram)*8/cram->n);
  printf("%d %lf\n", 100, (double)cram_usedmemory(cram)*8/n);
//  printf("%d %lf\n", 100, cram_entropy(cram));
//  printf("%d %lf\n", 100, cram_divergence(cram));
#endif

//  cram_save(cram, "cram1.dat", "cram2.dat");
//  cram_free(cram);

#if 0
  // ディスクに正しく保存されているか確認
  {
    cram = cram_load("cram1.dat", "cram2.dat");

//    map = mymmap(argv[1]);
//    text2 = map->addr;

  }
#endif

#if 0
  map = mymmap(argv[1]);
  text2 = map->addr;
  n2 = cram->n;
  b = 1024;
  mymalloc(buf, b);
  mygettime(&before);
  for (s=0; s<n2; s+=b) {
    if (s % (b*100) == 0) {
      fprintf(stderr, "%ld/%ld \r", s, n2);  fflush(stderr);
    }
    t = s+b-1;  if (t >= n2) t = n2-1;
    cram_read(cram, s, t, buf);
    //printf("block %ld: ", s/128);
    for (i=s; i<=t; i++) {
      if (buf[i-s] != text2[i]) {
        printf("??? s = %ld i = %ld [%c] [%c]\n", s, i, buf[i-s], text2[i]);
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
