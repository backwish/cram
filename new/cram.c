#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include <assert.h>
#include "cram.h"
#include "mman.h"

#define mymalloc(p,n) {p = malloc((n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory in line %d\n",__LINE__); exit(1);};}
#define D (8*sizeof(uchar))
#define TBLWIDTH 16
#define min(a, b) (((a)<(b))?(a):(b))
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
static u64 setCode(dblock_entry *blk, i64 i, int d, u64 x)
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
static u64 getCode(dblock_entry *blk, i64 i)
{
  u64 code;
  i64 j,l;
  int ii, bs;

  bs = block_len(blk);
  j = i / D;
  l = i % D;
  code = 0;
  for (ii=0; ii<4+1; ii++) {
    code <<= 8;
    if (j+ii < bs) {
      code += *dblock_addr(blk, j+ii);
    }
  }
  code = ((code >> (D-l)) & 0xffffffff);  
  code<<=32;
  return code;
}
static u64 PADDING(CRAM *cram){
  if(cram->paddingLen==0) return 0;
  return (3L<<(64-cram->paddingLen));
}
static int getEncodedTreeIndex(CRAM *cram, i64 b){
  return cram->encodedTree[b];
}
static void setEncodedTreeIndex(CRAM *cram,i64 b){
  cram->encodedTree[b] = cram->currentTreeIndex;
}
static void getChangeMinFreq(CRAM *cram,int c){
  int targetLen = TBLWIDTH+cram->paddingLen;
  int minFreq = 2147483647;
  int i;
  for(i=0;i<cram->sigma;i++){
    if(cram->huf[c]->clen[i] == targetLen && cram->freq[c][i] < minFreq){
      minFreq = cram->freq[c][i];
    }
  }
  assert(minFreq < 2147483647);
  cram->changeMinFreq[c] = minFreq;
}

static void catchReadWriteException(CRAM *cram,i64 b){
  if (b < 0 || b >= cram->nb) {
    printf("cram_write_block: b = %ld #blocks = %ld\n", b, cram->nb);
    exit(1);
  }
}
static int getCodelenFromHuffman(CRAM *cram,int ch){
  return cram->writeBlockUsingSaved ? cram->currentLen[cram->currentTreeIndex][ch] : cram->huf[cram->currentTreeIndex]->clen[ch];
}
static u64 getRightAlignedCodeFromHuffman(CRAM *cram,int ch,int len){
  return (cram->writeBlockUsingSaved ? cram->currentCode[cram->currentTreeIndex][ch] : cram->huf[cram->currentTreeIndex]->code[ch]) >> (64-len);
}
static int getCodeSize(CRAM *cram,int codeNum,uchar *buf){
  int ch,codeSize;  
  codeSize = 0;
  for (int i=0; i<codeNum; i+=cram->k) {
    ch = getuint(buf, i, codeNum, cram->k);
    codeSize += getCodelenFromHuffman(cram,ch);
  }
  return codeSize;
}

static void cram_write_block(CRAM *cram,int hufIndex,i64 b, uchar *buf)
{  
  i64 i,j,codeNum,codeSize;
  int codePos, ch, len;
  dblock_entry blk;      

  catchReadWriteException(cram,b);
  setEncodedTreeIndex(cram,b);  
  codeNum = cram->bs;  if (b*codeNum >= cram->n) codeNum = cram->n % cram->bs;
  codeSize = getCodeSize(cram,codeNum,buf);
  
  darray_change(cram->da, b, (codeSize+7)/8); // change the block size
  blk = darray_address(cram->da, b);  
  codePos = 0;
  for (i=0; i<codeNum; i+=cram->k) {
    ch = getuint(buf, i, codeNum, cram->k);    
    len = getCodelenFromHuffman(cram,ch);
    assert(len>=0);
    setCode(&blk, codePos, len, getRightAlignedCodeFromHuffman(cram,ch,len));      
    codePos += len;
  }
}
static i64* freqInit(CRAM *cram,uchar *text,int nc){
  int c,i;
  i64 *freq;  
  for (i=0; i<nc; i++) {
    mymalloc(cram->freq[i], cram->sigma);
  }      
  freq = cram->freq[0];
  for (c=0; c<cram->sigma; c++) freq[c] = 1;  
  for (i=0; i<(cram->n+cram->k-1)/cram->k; i++) {    
    c = getuint(text, i*cram->k, cram->n, cram->k);
    freq[c]++;
  }    
  for (i=1; i<nc; i++) {
    for (c=0; c<cram->sigma; c++) {
      cram->freq[i][c] = freq[c];
    }
  }
  return freq;
}
static void initBlocks(CRAM *cram,uchar* text){
  uchar *tmpbuf;
  int i,j;
  tmpbuf = malloc(cram->bs);
  for (i=0; i<cram->nb; i++) { // compress i-th block
    for(j=0;j<cram->bs;j++) tmpbuf[j] = text[i*cram->bs+j];
    cram_write_block(cram,0,i,tmpbuf);    
  }
  free(tmpbuf);
}
static void initHuffmanTree(CRAM *cram,int nc){
  int i;
  for (i=0; i<nc; i++) {
    cram->huf[i]  = makeHuffmanTree(cram->sigma, cram->freq[0], 8*cram->k,cram->paddingLen);
  }
}
void catchInitException(CRAM *cram){
  if (cram->bs % cram->k != 0) {
    printf("cram_initialize: bs = %d k = %d\n", cram->bs, cram->k);
    exit(1);
  }      
}

CRAM *cram_init(char *filename, int bs, int bu, int paddingLen)
{
  CRAM *cram;
  MMAP *map;  
  int nc,i;  
  uchar *text;
  Huffman *huf;        
  
  mymalloc(cram,1);
  cram->paddingLen = paddingLen;
  cram->k = 2;  
  cram->bs = bs;  

  map = mymmap(filename);
  cram->n = map->len;
  text = map->addr;  
  cram->nb = (cram->n+bs-1)/bs; // number of blocks
  cram->currentTreeIndex = 0;
  cram->bu = bu;
  cram->step = bu==0 ? 0 : cram->nb*(cram->bu-1)/(cram->bu); 
  cram->da = darray_initialize(cram->nb, 1, cram->bs*2, cram->n/4, cram->n/100);
  cram->writeBlockUsingSaved = 0;
  cram->sigma = 1 << (cram->k*8);
  nc = 2;  
  catchInitException(cram);  

  //mymalloc(cram->encodedTree,(nb+63)/64);
  //for(int i=0;i<(nb+63)/64;i++) cram->encodedTree[i] = 0;
  mymalloc(cram->encodedTree,cram->nb);
  freqInit(cram,text,nc);  
  initHuffmanTree(cram,nc);  
  initBlocks(cram,text);  
  
  return cram;
}
void buildGlobalAdditional(CRAM *cram){
  cram->codeAllocator = makeCodeAllocator(TBLWIDTH);
  mymalloc(cram->pivotFreq,cram->sigma);
}
void freeGlobalAdditional(CRAM *cram){
  freeCodeAllocator(cram->codeAllocator);
  free(cram->pivotFreq);
}
void buildLocalAdditional(CRAM *cram,int hufIndex){
  int c;
  cram->additionalCodeTable[hufIndex] = makeAdditionalCodeTable(TBLWIDTH);
  getChangeMinFreq(cram,hufIndex);  
  mymalloc(cram->currentCode[hufIndex],cram->sigma);
  mymalloc(cram->currentLen[hufIndex],cram->sigma);      
  for(c=0;c<cram->sigma;c++){
    cram->pivotFreq[c] = cram->freq[hufIndex][c];
    cram->currentLen[hufIndex][c] = cram->huf[hufIndex]->clen[c];
    cram->currentCode[hufIndex][c] = cram->huf[hufIndex]->code[c];
  }
}
void freeLocalAdditional(CRAM *cram,int hufIndex){
  freeAdditionalCodeTable(cram->additionalCodeTable[hufIndex]);
  free(cram->currentCode[hufIndex]);
  free(cram->currentLen[hufIndex]);  
}
void initAdditional(CRAM *cram,int changeMultiple){
  int c;
  cram->canAddCode = 1;
  cram->newCodeIdx = 0;
  cram->changeMultiple = changeMultiple;    
  cram->writeBlockUsingSaved = 1;
  buildGlobalAdditional(cram);
  for(c=0;c<2;c++) buildLocalAdditional(cram,c);
}
void freeAdditional(CRAM *cram){
  int c;
  freeGlobalAdditional(cram);
  for(c=0;c<2;c++) freeLocalAdditional(cram,c);
}
void replaceAdditional(CRAM *cram,int hufIndex){
  freeGlobalAdditional(cram);
  freeLocalAdditional(cram,hufIndex);
  buildGlobalAdditional(cram);
  buildLocalAdditional(cram,hufIndex);
  cram->canAddCode = 1;
}
CRAM *cram_initialize(char *filename, int bs, int bu, int paddingLen,int changeMultiple){
  CRAM *cram;  
  assert(paddingLen>=0);
  cram = cram_init(filename,bs,bu,paddingLen);  
  if(paddingLen>0){
    initAdditional(cram,changeMultiple);
  }
  return cram;
}

void cram_free(CRAM *cram)
{
  int i, nc,c;
  nc = 2;
  for (i=0; i<nc; i++) {
    free(cram->freq[i]);
    freeHuffman(cram->huf[i]);
  }
  free(cram->encodedTree);
  if(cram->paddingLen>0){
    freeAdditional(cram);    
  }
  darray_free(cram->da);
  free(cram);
}

i64 cram_usedmemory(CRAM *cram)
{
  i64 size;
  size = sizeof(CRAM);
  fprintf(stderr, "cram size : %ld\n", size);
  size += darray_usedmemory(cram->da);
  //fprintf(stderr,"darray used memory %ld bytes (%1.2f bpc)\n", darray_usedmemory(cram->da), (double)darray_usedmemory(cram->da)*8/cram->n);
  size += cram->sigma * sizeof(*cram->freq[0]) * 3;
  size += Huffman_usedmemory(cram->huf[0]) * 2;  
  //return size;
  return darray_usedmemory(cram->da);
}

static int checkUsingTable(CRAM *cram,u64 code){
  if(cram->paddingLen > 0 && (code&PADDING(cram))==PADDING(cram)) return 1;
  return 0;
}
static void writeToBuffer(CRAM *cram,i64 i,int ch,uchar *buf){  
  for (int j=0; j<cram->k; j++) {
    buf[i+j] = ch >> (8*(cram->k-1-j));
  }   
}

int cram_read_block(CRAM *cram, i64 b, uchar *buf) // decode b-th block into buf
{
  i64 i,j,codeNum;  
  int codePos, ch,encodedTreeIndex;
  dblock_entry blk;
  Huffman *huf;
  u64 code;  
  TableElement tableElement;
  catchReadWriteException(cram,b);

  encodedTreeIndex = getEncodedTreeIndex(cram, b);
  huf = cram->huf[encodedTreeIndex];
  blk = darray_address(cram->da, b);
  codeNum = cram->bs;  
  if(b*codeNum >= cram->n) codeNum = cram->n % cram->bs;
  codePos = 0;  
    
  for (i=0; i<codeNum; i+=cram->k) {
    code = getCode(&blk, codePos);
    if(checkUsingTable(cram,code)) {      
      tableElement = decodeFromAdditionalCodeTable(cram->additionalCodeTable[encodedTreeIndex],code<<cram->paddingLen);      
      ch = tableElement.ch;
      codePos += tableElement.len+cram->paddingLen;
    }else{      
      ch = decodeFromHuffman(huf, code);
      codePos += huf->clen[ch];
    }
    writeToBuffer(cram,i,ch,buf);
  }
  return codePos;
}


static u64 getNewCodeWithPadding(CRAM *cram,int totalLen){
  u64 newCode;  
  newCode = codeAllocate(cram->codeAllocator,totalLen - cram->paddingLen);  
  if(newCode==1) return newCode;
  newCode = (newCode>>(cram->paddingLen));
  newCode |= PADDING(cram);  
  return newCode;
}
static int checkNotInsertCode(CRAM *cram){
  if(cram->paddingLen>0 && cram->canAddCode) return 0;
  return 1;
}
static int checkNotMakeNewCode(CRAM *cram,int currFreq,int c){
  if(cram->paddingLen==0 || !cram->canAddCode) return 1;
  if(currFreq > cram->changeMinFreq[cram->currentTreeIndex] && 
  currFreq > cram->pivotFreq[c]*cram->changeMultiple) return 0;
  return 1;
}
static void makeNewCode(CRAM *cram,int c){
  u64 newCode;
  int currFreq = cram->freq[cram->currentTreeIndex][c];    
  if(checkNotMakeNewCode(cram,currFreq,c)) return;
  cram->pivotFreq[c] = currFreq;
  cram->newCodeLen[cram->newCodeIdx] = min(cram->currentLen[cram->currentTreeIndex][c]-1,TBLWIDTH+(cram->paddingLen));
  newCode = getNewCodeWithPadding(cram,cram->newCodeLen[cram->newCodeIdx]);     
  if(newCode==1){      
    cram->canAddCode = 0;
    return;
  }            
  assert((newCode & ((1L<<(64-(cram->newCodeLen[cram->newCodeIdx])))-1)) == 0 );
  cram->newCodeCh[cram->newCodeIdx] = c;
  cram->newCode[cram->newCodeIdx] = newCode;
  cram->newCodeIdx++;  
}
static void insertCode(CRAM *cram){
  int i,len,ch;
  u64 code;    
  if(checkNotInsertCode(cram)) return;    
  for(i=0;i<cram->newCodeIdx;i++){
    code = cram->newCode[i];
    len = cram->newCodeLen[i];
    ch = cram->newCodeCh[i];
    cram->currentCode[cram->currentTreeIndex][ch] = code;
    cram->currentLen[cram->currentTreeIndex][ch] = len;
    code<<=(cram->paddingLen);
    assert((code & ( (1L<< (64-(len-cram->paddingLen)) ) -1)) == 0);      
    insertToAdditionalCodeTable(cram->additionalCodeTable[cram->currentTreeIndex],code,ch,len-cram->paddingLen);
  }
  cram->newCodeIdx = 0;
}
static int changeFreqUnit(CRAM *cram,int i,uchar *tmpbuf,int delta){
  int c = getuint(tmpbuf, i*cram->k, cram->bs, cram->k);
  cram->freq[cram->currentTreeIndex][c]+=delta;
  return c;
}
static void changeFreq(CRAM *cram,uchar *tmpbuf,uchar *buf){
  int i,c;
  for (i=0; i<(cram->bs+cram->k-1)/cram->k; i++){
    changeFreqUnit(cram,i,tmpbuf,-1);
  }
  for (i=0; i<cram->bs; i++) tmpbuf[i] = *buf++;
  for (i=0; i<(cram->bs+cram->k-1)/cram->k; i++){
    c = changeFreqUnit(cram,i,tmpbuf,1);
    makeNewCode(cram,c);
  }
}
static void changeBlock(CRAM *cram,i64 b,uchar *buf){  
  int i,c;
  uchar *tmpbuf;  
  mymalloc(tmpbuf, cram->bs);
  cram_read_block(cram, b, tmpbuf);      
  changeFreq(cram,tmpbuf,buf);  
  cram_write_block(cram, cram->currentTreeIndex, b, tmpbuf);      
  free(tmpbuf);
}
static void reencodeOneBlock(CRAM *cram,int b){
  uchar *tmpbuf;
  mymalloc(tmpbuf, cram->bs);
  cram_read_block(cram, b, tmpbuf);
  cram_write_block(cram,cram->currentTreeIndex,b, tmpbuf);
  free(tmpbuf);
}
static void reencodeBlocks(CRAM *cram){
  int i;  
  for(i=0;i<cram->bu;i++){
    if(cram->step >= cram->nb) break;    
    reencodeOneBlock(cram,cram->step);
    cram->step++;
  }  
}
static void reencodeAll(CRAM *cram){
  int i;
  for(i=0;i<cram->nb;i++){
    reencodeOneBlock(cram,i);
  }
}
static void rebuildHuffmanTree(CRAM *cram){
  int c;    
  for (c=0; c<cram->sigma; c++) {
    cram->freq[(cram->currentTreeIndex+1) % 2][c] = cram->freq[cram->currentTreeIndex][c];
  }
  cram->currentTreeIndex++;
  cram->currentTreeIndex%=2;
  freeHuffman(cram->huf[cram->currentTreeIndex]);
  cram->huf[cram->currentTreeIndex] = makeHuffmanTree(cram->sigma, cram->freq[(cram->currentTreeIndex+1) % 2], 8*cram->k,cram->paddingLen);    

  if(cram->paddingLen>0){
    replaceAdditional(cram,cram->currentTreeIndex);
  }
  cram->step = 0;  
}
int checkRebuildHuffmanTree(CRAM *cram){
  if(cram->bu > 0){
    return cram->step>=cram->nb ? 1 : 0;
  }else{
    return !cram->canAddCode;
  }
}
void cram_write_replacePart(CRAM *cram,i64 b,uchar *buf){  
  changeBlock(cram,b,buf);
  reencodeBlocks(cram);
  insertCode(cram);
  if(checkRebuildHuffmanTree(cram)){
    rebuildHuffmanTree(cram);
  }
}
void cram_write_replaceAll(CRAM *cram,i64 b,uchar *buf){
  changeBlock(cram,b,buf);
  insertCode(cram);
  if(checkRebuildHuffmanTree(cram)){    
    rebuildHuffmanTree(cram);
    reencodeAll(cram);
  }
}
void cram_write(CRAM *cram,i64 b,uchar *buf){
  if(cram->bu>0){
    cram_write_replacePart(cram,b,buf);
  }else if(cram->bu==0){
    cram_write_replaceAll(cram,b,buf);
  }
}