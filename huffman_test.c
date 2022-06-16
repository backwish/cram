#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "huffman.h"
#include "heap.h"
#include "mman.h"
#include "memory.h"

#define mymalloc(p,n) {p = malloc((n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory in line %d\n",__LINE__); exit(1);};}
#define SIGMA 65536
#define D (8*sizeof(uchar))
static u64 getuint(uchar *s, i64 i, i64 n, i64 w)
{
  u64 x,c;
  i64 j;

  x = 0;
  for (j=0; j<w; j++) {
    if (i+j < n) c = s[i+j]; else c = 0;
    x += c << ((w-1-j)*8); // big endian
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
    mp = dblock_addr(blk, iq);
    *mp = (*mp & (~m)) | y;
    iq++;  ir=0;
    d -= d2;
    x &= (1<<d)-1;
  }
  m = (1<<d)-1;
  y = x << (D-ir-d);
  m <<= (D-ir-d);
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
i64 *buildFreq(uchar *text,i64 n,int k){
    i64 *freq;
    int c,i;
    mymalloc(freq,SIGMA);
    for(c=0;c<SIGMA;c++) freq[c] = 1;
    for(i=0;i<(n+k-1)/k;i++){
        c = getuint(text,i*k,n,k);
        freq[c]++;
    }
    return freq;
}
void compressBlock(Huffman2 *huf,darray *da,uchar *buf,i64 b,i64 bs){
    int i,c,k,csize,r,len;
    dblock_entry blk;
    k=2;
    for(i=0;i<bs;i+=k){
        c = getuint(buf,i,bs,k);
        csize+=huf->clen[c];
    }
    darray_change(da,b,(csize+7)/8);
    blk = darray_address(da,b);
    r = 0;
    for(i=0;i<bs;i+=k){
        c = getuint(buf,i,bs,k);
        len = huf->clen[c];
        setbits(&blk,r,len,huf->code[c] >> (64-len));
        r+=len;
    }
}
void compress(Huffman2 *huf,darray *da,uchar *text,i64 nb,i64 bs){    
    int i,j;
    uchar *tmpbuf;
    tmpbuf = malloc(bs);
    for(i=0;i<nb;i++){
        for(j=0;j<bs;j++) tmpbuf[j] = text[i*bs+j];
        compressBlock(huf,da,tmpbuf,i,bs);
    }    
    free(tmpbuf);
}
int decompressBlock(Huffman2 *huf,darray *da,uchar *text,i64 b,i64 bs){
    int i,j,r,k,c;
    u64 x;
    dblock_entry blk;
    uchar *buf;
    buf = malloc(bs);
    blk = darray_address(da,b);
    r = 0;k = 2;
    for(i=0;i<bs;i+=k){
        x = getbit32(&blk,r)<<32;
        c = DecodeHuffman2_tbl(huf,x);
        r+=huf->clen[c];        
        for(j=0;j<k;j++){
            *buf++ = c >> (8*(k-1-j));
        }
    }
    for(i=0;i<bs;i++){
        assert(text[b*bs+i]==buf[i]);
        if(text[b*bs+i]!=buf[i]){
            printf("FAIL DECOMPRESS AT %d %d\n",b,i);
            return 0;
        } 
    }
    return 1;
}
int decompress(Huffman2 *huf,darray *da,uchar *text,i64 nb,i64 bs){
    int i,j;
    for(i=0;i<nb;i++){
        if(decompressBlock(huf,da,text,i,bs)==0) return 0;
    }
    return 1;
}


int testHuffman(uchar* src,int mode){
    MMAP *map;
    i64 n,bs,nb;
    uchar *text;
    Huffman2 *huf;
    i64 *freq;
    int k,c,i,j,csize,ret;
    darray *da;
    dblock_entry blk;
    
    bs = 1024;
    k = 2;
    assert(n%bs==0);
    nb = (n+bs-1)/bs;    
    map = mymmap(src);
    text = map->addr;
    n = map->len;
    freq = buildFreq(text,n,k);
    huf = MakeHuffman2Tree2(SIGMA,freq,8*k,mode);    
    if(mode==1){
        int r = huf->right[huf->rootIdx];
        assert(huf->right[r]==huf->right[0] && huf->left[r]==huf->left[0]);
    }else if(mode==2){        
        int r = huf->right[huf->right[huf->rootIdx]];
        printf("%d %d %d\n",huf->rootIdx,huf->right[huf->rootIdx],r);
        assert(huf->right[r]==huf->right[0] && huf->left[r]==huf->left[0]);
    }
    free(freq);   

    printf("HUFFMAN TREE BUILD FINISHED MODE : %d\n",mode);

    da = darray_initialize(nb, 1, bs*2, n/4, n/100); 
    printf("COMPRESS START\n");   
    compress(huf,da,text,nb,bs);
    printf("DECOMPRESS START\n");
    ret = decompress(huf,da,text,nb,bs);
    printf("DECOMPRESS END %d\n",ret);
    freeHuffman2(huf);
    free(text);
    return ret;
}


int main(int argc, char *argv[]){
    Huffman2 *huf;
    int i,ret;
    for(i=2;i<=2;i++){
        ret = testHuffman(argv[1],i);
        if(ret==1){
            printf("HUFFMAN TEST SUCCESS MODE %d\n",i);
        }else{
            printf("HUFFMAN TEST FAIL MODE %d\n",i);
        }
    }
    return 0;
}