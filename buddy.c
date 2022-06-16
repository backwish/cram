#include "typedef.h"
#include "buddy.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define mymalloc(p,n) {p = malloc((n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory in line %d\n",__LINE__); exit(1);};}
#define INF 100
#define max(a, b) (((a)>(b))?(a):(b))
#define min(a, b) (((a)<(b))?(a):(b))
static int left_child(int idx){
    return idx*2;
}
static int right_child(int idx){
    return idx*2+1;
}
static int parent(int idx){
    return idx/2;
}
static int checkEndOfCodeLen(int i){
    return (i&(i+1))==0 ? 1 : 0;
}
static void initCodeLen(CodeAllocator *self,int nodeNum){
    int len = 1;
    mymalloc(self->requiredMinCodeLen,nodeNum);
    self->requiredMinCodeLen[0] = self->requiredMinCodeLen[1] = 1;
    for(int i=2;i<nodeNum;i++){
        self->requiredMinCodeLen[i] = len;        
        if(checkEndOfCodeLen(i)) len++;
    }
}
CodeAllocator *makeCodeAllocator(int maxDepth){
    CodeAllocator *p;
    int nodeNum = (1<<(maxDepth+1));        
    mymalloc(p,1);
    initCodeLen(p,nodeNum);    
    return p;
}
void freeCodeAllocator(CodeAllocator *self){
    free(self->requiredMinCodeLen);
    free(self);
}
static int checkCannotAllocate(CodeAllocator *self,int size){
    return size < self->requiredMinCodeLen[1] ? 1 : 0;
}
static void recalcRequirement(CodeAllocator *self,int pos){
    int l,r;
    self->requiredMinCodeLen[pos] = INF;
    pos = parent(pos);
    while(pos){
        l = left_child(pos);
        r = right_child(pos);
        self->requiredMinCodeLen[pos] = min(self->requiredMinCodeLen[l],self->requiredMinCodeLen[r]);
        pos = parent(pos);
    }
}
static int findCodePos(CodeAllocator *self,u64 *code,int size){
    int pos=1,l,r;
    for(int i=1;i<=size;i++){
        l = left_child(pos);
        r = right_child(pos);
        if(self->requiredMinCodeLen[l] > size){
            pos = r;
            (*code)|=(1L<<(sizeof(u64)*8-i));
        }else if(self->requiredMinCodeLen[r] > size){
            pos = l;
        }else{
            pos = (self->requiredMinCodeLen[l] >= self->requiredMinCodeLen[r]) ? l : r;
        }
    }
    return pos;
}
u64 codeAllocate(CodeAllocator *self,int size){
    int pos,i,l,r,go;
    u64 code = 0;    
    if(checkCannotAllocate(self,size)) return 1;    
    pos = findCodePos(self,&code,size);    
    recalcRequirement(self,pos);    
    return code;
}


AdditionalCodeTable *makeAdditionalCodeTable(int maxDepth){
    AdditionalCodeTable *p;
    int nodeNum = (1<<maxDepth),i;
    mymalloc(p,1);
    mymalloc(p->tbl,nodeNum);
    p->maxDepth = maxDepth;
    for(i=0;i<nodeNum;i++){
        (p->tbl[i]).ch = -1;
        (p->tbl[i]).len = -1;
    }
    return p;
}
void freeAdditionalCodeTable(AdditionalCodeTable *self){
    free(self->tbl);
    free(self);
}
void insertToAdditionalCodeTable(AdditionalCodeTable *self,u64 code,int ch,int len){    
    int tableIndex = (code>>64-self->maxDepth);
    assert(tableIndex < (1<<(self->maxDepth)));
    assert(self->tbl[tableIndex].ch==-1);//prefix free
    self->tbl[tableIndex].ch = ch;
    self->tbl[tableIndex].len = len;
}
TableElement decodeFromAdditionalCodeTable(AdditionalCodeTable *self,u64 code){    
    int tableIndex = (code>>64-self->maxDepth);
    assert(tableIndex < (1<<(self->maxDepth)));

    for(int i=0;i<self->maxDepth;i++){
        if(self->tbl[tableIndex].ch>=0) return self->tbl[tableIndex];
        tableIndex&=~(1<<i);        
    }
    assert(self->tbl[tableIndex].ch>=0);//assure exists
    return self->tbl[tableIndex];
}