#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "huffman.h"
#include "heap.h"

#define ID_HUFFMAN     0x05
#define mymalloc(p,n) {p = malloc((n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory in line %d\n",__LINE__); exit(1);};}


struct freqs2 {
  i64 freq;
  int idx;
};

static int tmpf_cmp2(int *s, int *t, i64 w)
{
  if (*s < *t) return -1;
  if (*s > *t) return  1;
  return 0;
}


static u64 getuint(uchar *s, i64 i, i64 w)
{
  u64 x;
  i64 j;
  s += i*w;
  x = 0;
  for (j=0; j<w; j++) {
    x += ((u64)(*s++)) << (j*8);
  }
  return x;
}
static int NODENUM(int n){
  return 2*n+3;
}
static int NILNODE(int n){
  return 2*n;
}



void freeHuffman(Huffman *p)
{
  free(p->left);
  free(p->right);
  free(p->clen);
  free(p->code);
  if (p->tbl_width > 0) free(p->tbl);
  free(p);
}
static int decodeNaive(Huffman *h, u64 x)
{
  unsigned int i = sizeof(x)*8 - 1;
  int p = h->rootIdx;
  while (1) {
    if ((x >> i) & 1) p = h->right[p]; else p = h->left[p];
    i--;    
    if (p < h->n) return p;
  }
}
void makeDecodeTable(Huffman *h)
{
  u64 i;
  u64 v,w;
  int width;

  width = h->tbl_width;
  mymalloc(h->tbl, 1<<width);

  for (i = 0; i < (1<<width); i++) {
    w = i << (sizeof(u64)*8-width); // bit stream
    v = decodeNaive(h, w);
    if (h->clen[v] <= width) {
      h->tbl[i] = v;
    } else {
      h->tbl[i] = -1;
    }
  }
}
int decodeFromHuffman(Huffman *h, u64 x)
{
  u64 w;
  
  if (h->tbl_width > 0) {
    w = x >> (sizeof(u64)*8-h->tbl_width);
    if (h->tbl[w] >= 0) return h->tbl[w];
  }
  return decodeNaive(h,x);
}



static Huffman *initVariables(int n,int tbl_width)
{
  int i;
  Huffman *p;
  mymalloc(p,1);
  p->n = n;    
  p->tbl_width = tbl_width;
  mymalloc(p->left,NODENUM(n));
  mymalloc(p->right,NODENUM(n));
  mymalloc(p->clen,NODENUM(n));  
  mymalloc(p->code,NODENUM(n));
  for (i=0; i<NODENUM(n); i++) p->clen[i] = 0;
  for (i=0; i<NODENUM(n); i++) p->left[i] = p->right[i] = NILNODE(n);
  return p;
}
static i64 getFreqSum(i64 *freq,int n){
  i64 freqSum = 0;
  for(int i=0;i<n;i++){
    freqSum+=freq[i];
  }
  return freqSum;
}
static int getLeftTreeLeafNum(i64 *freq,int n,double ratio){    
  i64 freqSum = getFreqSum(freq,n),pSum = 0;
  int ret = 0;
  for(int i=0;i<n;i++){
    pSum+=freq[i];
    if(pSum-freq[i] < freqSum*ratio && pSum>freqSum*ratio) return i;
  }
  return -1;
}
void setChildren(Huffman *h,int p,int l,int r){
  h->left[p] = l;
  h->right[p] = r;
}
int buildEdgePart(int leafStart,int leafEnd,int nonLeafStart,Huffman *h,i64 *freq){
  HEAP H;
  int i,n,l,r;
  n = leafEnd - leafStart;  
  struct freqs2 *tmpf, tmin1, tmin2, tnew;
  mymalloc(tmpf,NODENUM(n));
  for(i=0;i<n;i++){
    tmpf[i+1].idx = i+leafStart;
    tmpf[i+1].freq = (freq[i+leafStart]+3)/4;
  }
  heap_build(&H, n, (uchar *)tmpf, n*2, sizeof(tmpf[0]), sizeof(tmpf[0].freq),
            (void *)&(tmpf[0].freq) - (void *)&tmpf[0], (int (*)(uchar *, uchar *, i64))tmpf_cmp2);
  r = nonLeafStart-1;
  for(i=0;i<n-1;i++){
    heap_extract(&H, (uchar *)&tmin1);
    heap_extract(&H, (uchar *)&tmin2);
    r++;
    setChildren(h,r,tmin1.idx,tmin2.idx);    
    tnew.idx = r;
    tnew.freq = tmin1.freq + tmin2.freq;
    heap_insert(&H, (uchar *)&tnew);     
  }
  free(tmpf);
  return r;
}

int buildEdge(Huffman *h,i64 *freq,int paddingLen,int leftTreeLeafNum){
  int root,l,r,rl,rr;
  if(paddingLen==0){
    root = buildEdgePart(0,h->n,h->n,h,freq);
  }else if(paddingLen==1){
    l = buildEdgePart(0,h->n,h->n,h,freq);
    r = l+1;
    root = l+2;
    setChildren(h,root,l,r);    
  }else if(paddingLen==2){
    l = buildEdgePart(0,leftTreeLeafNum,h->n,h,freq);
    rl = buildEdgePart(leftTreeLeafNum,h->n,l+1,h,freq);
    rr = rl+1;
    r = rl+2;
    root = rl+3;
    setChildren(h,root,l,r);
    setChildren(h,r,rl,rr);    
  }
  return root;
}
static void makeCode(int n, int r, int d, Huffman *h, u64 c)
{
  h->clen[r] = d;
  if (h->left[r] < NILNODE(n)) {
    makeCode(n,h->left[r],d+1,h, c);
  } else {    
    h->code[r] = c;
  }
  if (h->right[r] < NILNODE(n)) {
    makeCode(n,h->right[r],d+1,h, c + (1L<<(sizeof(c)*8-1-d)));
  } else {
    h->code[r] = c;
  }
}

Huffman *makeHuffmanTree(int n, i64 *freq, int tbl_width,int paddingLen)
{  
  int i,j;
  Huffman *h;    
  int l,r,r2;
  int m1,m2,leftTreeLeafNum,freqSum,pSum;    
  double ratio = (double)2/3;
  h = initVariables(n,tbl_width);  
  leftTreeLeafNum = getLeftTreeLeafNum(freq,n,ratio);  
  assert(leftTreeLeafNum>=0); 
  
  h->rootIdx = buildEdge(h,freq,paddingLen,leftTreeLeafNum);
  makeCode(n,h->rootIdx,0,h,0);  
  makeDecodeTable(h);
  return h;
}

i64 Huffman_usedmemory(Huffman *h)
{
  i64 size;
  size = sizeof(*h);
  size += (NODENUM(h->n)) * sizeof(*h->left);
  size += (NODENUM(h->n)) * sizeof(*h->right);
  size += (NODENUM(h->n)) * sizeof(*h->clen);
  size += (NODENUM(h->n)) * sizeof(*h->code);
  if (h->tbl_width > 0) size += (1 << h->tbl_width) * sizeof(*h->tbl);
  return size;
}
