#include <stdio.h>
#include <stdlib.h>
#include "buddy.h"
#include "typedef.h"
#define TABLECASES 6
#define CASES 8
#define CASESIZE 4
#define CASEMAXDEPTH 16
int input[CASES][CASESIZE] = 
{{1,2,3,3},{1,3,2,3},{1,3,3,2},{2,1,3,3},
{2,3,1,3},{2,3,3,1},{1,1,2,2},{1,2,1,2}};
int output[CASES][CASESIZE] = 
{{0,2,6,7},{0,4,3,5},{0,4,5,3},{0,1,2,3},
{0,2,1,3},{0,2,3,1},{0,1,-1,-1},{0,2,-1,3}};
int decoded[TABLECASES][CASESIZE] = 
{{5,10,15,20},{6,11,16,21},{55,44,33,22},
{777,666,555,444},{0,0,0,0},{65535,65535,65535,65535}};

int shiftCode(u64 code,int len){
    if(code==1) return -1;
    return (code>>(sizeof(u64)*8-len));
}
void testAdditionalTable(){
    TableElement ele;
    int i,j,len,ch;
    AdditionalCodeTable *additionalCodeTable;
    u64 code;
    for(i=0;i<TABLECASES;i++){
        additionalCodeTable = makeAdditionalCodeTable(CASEMAXDEPTH);
        for(int j=0;j<CASESIZE;j++){            
            len = input[i][j];
            code = (u64)output[i][j];
            code<<=(sizeof(u64)*8-len);
            ch = decoded[i][j];
            insertToAdditionalCodeTable(additionalCodeTable,code,ch,len);            
        }
        for(int j=0;j<CASESIZE;j++){
            len = input[i][j];
            code = (u64)output[i][j];
            code<<=(sizeof(u64)*8-len);            
            ele = decodeFromAdditionalCodeTable(additionalCodeTable,code);
            if(ele.ch!=decoded[i][j] || ele.len!=len){
                printf("case : %d, idx : %d, output : %d answer : %d\n",i,j,ch,decoded[i][j]);
                printf("ADDITIONAL TABLE TEST FAILED\n");
                return;
            }
        }
        freeAdditionalCodeTable(additionalCodeTable);
    }
    printf("ADDITIONAL TABLE TEST SUCCESS\n");;
}
void testCodeAllocate(){
    int i,j,shiftedCode;
    CodeAllocator *codeAllocator;
    u64 code;
    for(i=0;i<CASES;i++){
        codeAllocator = makeCodeAllocator(CASEMAXDEPTH);
        for(j=0;j<CASESIZE;j++){
            code = codeAllocate(codeAllocator,input[i][j]);
            shiftedCode = shiftCode(code,input[i][j]);
            if(shiftedCode!=output[i][j]){
                printf("case : %d, idx : %d, output : %d answer : %d\n",i,j,shiftedCode,output[i][j]);
                printf("CODE ALLOCATE TEST FAILED\n");
                return;
            }            
        }
        freeCodeAllocator(codeAllocator);
    }
    printf("CODE ALLOCATE TEST SUCCESS\n");;
}

int main(){
    int maxDepth = 4;
    testCodeAllocate();
    testAdditionalTable();
}