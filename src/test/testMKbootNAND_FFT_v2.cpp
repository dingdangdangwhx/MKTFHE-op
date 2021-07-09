#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <fstream>


#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"


#include "mkTFHEparams.h"
#include "mkTFHEkeys.h"
#include "mkTFHEkeygen.h"
#include "mkTFHEsamples.h"
#include "mkTFHEfunctions.h"


using namespace std;


// **********************************************************************************
// ********************************* MAIN *******************************************
// **********************************************************************************


void dieDramatically(string message) {
    cerr << message << endl;
    abort();
} 

// Guaker 20210521 实验 - 拼接密文与门
void MKAND(MKLweSample *result,                         // 用于存放结果
            const MKLweSample *ca,                      // 输入密文 ca
            const MKLweSample *cb,                      // 输入密文 cb
            const MKLweBootstrappingKeyFFT_v2 *bkFFT,   // 输入 MK 自举密钥
            const LweParams* LWEparams,                 // 输入 LWE 参数
            const LweParams *extractedLWEparams,        // 输入 LWE 拓展参数
            const TLweParams* RLWEparams,               // 输入 RLWE 参数
            const MKTFHEParams *MKparams,               // 输入 MK 参数
            const MKRLweKey *MKrlwekey)                 // 输入 MKRLWE 密钥（仅使用公钥部分）
{
    clock_t begin_MKAND = clock();
    MKbootsNAND_FFT_v2m2(result, ca, cb, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsNAND_FFT_v2m2(result, result, result, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    clock_t end_MKAND = clock();
    double time_MKAND = ((double) end_MKAND - begin_MKAND)/CLOCKS_PER_SEC;
    cout << "拼接 MKAND 完成，耗时： " << time_MKAND << " 秒" << endl;
}

// Guaker 20210521 实验 - 拼接密文或门
void MKOR(MKLweSample *result,                         // 用于存放结果
            const MKLweSample *ca,                      // 输入密文 ca
            const MKLweSample *cb,                      // 输入密文 cb
            const MKLweBootstrappingKeyFFT_v2 *bkFFT,   // 输入 MK 自举密钥
            const LweParams* LWEparams,                 // 输入 LWE 参数
            const LweParams *extractedLWEparams,        // 输入 LWE 拓展参数
            const TLweParams* RLWEparams,               // 输入 RLWE 参数
            const MKTFHEParams *MKparams,               // 输入 MK 参数
            const MKRLweKey *MKrlwekey)                 // 输入 MKRLWE 密钥（仅使用公钥部分）
{
    clock_t begin_MKOR = clock();
    MKLweSample *temp1 = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp2= new_MKLweSample(LWEparams, MKparams);
    MKbootsNAND_FFT_v2m2(temp1, ca, ca, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsNAND_FFT_v2m2(temp2, cb, cb, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsNAND_FFT_v2m2(result, temp1, temp2, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    clock_t end_MKOR = clock();
    double time_MKOR = ((double) end_MKOR - begin_MKOR)/CLOCKS_PER_SEC;
    cout << "拼接 MKOR 完成，耗时： " << time_MKOR << " 秒" << endl;
}

//密文非门（取反）
void MKNOT(MKLweSample *result,                         // 用于存放结果
            const MKLweSample *ca,                      // 输入密文 ca
            const MKLweBootstrappingKeyFFT_v2 *bkFFT,   // 输入 MK 自举密钥
            const LweParams* LWEparams,                 // 输入 LWE 参数
            const LweParams *extractedLWEparams,        // 输入 LWE 拓展参数
            const TLweParams* RLWEparams,               // 输入 RLWE 参数
            const MKTFHEParams *MKparams,               // 输入 MK 参数
            const MKRLweKey *MKrlwekey)                 // 输入 MKRLWE 密钥（仅使用公钥部分）
{
    clock_t begin_MKNOT = clock();
    MKbootsNAND_FFT_v2m2(result, ca, ca, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    clock_t end_MKNOT = clock();
    double time_MKNOT = ((double) end_MKNOT - begin_MKNOT)/CLOCKS_PER_SEC;
    cout << "拼接 MKNOT 完成，耗时： " << time_MKNOT << " 秒" << endl;
}

//密文异或门（加法）
void MKXOR(MKLweSample *result,                         // 用于存放结果
            const MKLweSample *ca,                      // 输入密文 ca
            const MKLweSample *cb,                      // 输入密文 cb
            const MKLweBootstrappingKeyFFT_v2 *bkFFT,   // 输入 MK 自举密钥
            const LweParams* LWEparams,                 // 输入 LWE 参数
            const LweParams *extractedLWEparams,        // 输入 LWE 拓展参数
            const TLweParams* RLWEparams,               // 输入 RLWE 参数
            const MKTFHEParams *MKparams,               // 输入 MK 参数
            const MKRLweKey *MKrlwekey)                 // 输入 MKRLWE 密钥（仅使用公钥部分）
{
    clock_t begin_MKXOR = clock();
    MKLweSample *temp1 = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp2= new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp3 = new_MKLweSample(LWEparams, MKparams);

    MKbootsNAND_FFT_v2m2(temp1, ca, cb, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsNAND_FFT_v2m2(temp2, temp1, ca, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsNAND_FFT_v2m2(temp3, temp1, cb, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsNAND_FFT_v2m2(result, temp2, temp3, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    clock_t end_MKXOR = clock();
    double time_MKXOR = ((double) end_MKXOR - begin_MKXOR)/CLOCKS_PER_SEC;
    cout << "拼接 MKXOR 完成，耗时： " << time_MKXOR << " 秒" << endl;
}


// // 函数功能：实现自举的与门
// EXPORT void MKbootsAND_FFT_v2m2(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb,
//                                  const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams *LWEparams, const LweParams *extractedLWEparams,
//                                  const TLweParams *RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey);

// // 函数功能：实现自举的或门
// EXPORT void MKbootsOR_FFT_v2m2(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb,
//                                  const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams *LWEparams, const LweParams *extractedLWEparams,
//                                  const TLweParams *RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey);

// // 函数功能：实现自举的非门
// EXPORT void MKlweNegate(MKLweSample *result, const MKLweSample *sample, const MKTFHEParams *MKparams);

// // 函数功能：实现自举的异或门
// EXPORT void MKbootsXOR_FFT_v2m2(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb,
//                                  const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams *LWEparams, const LweParams *extractedLWEparams,
//                                  const TLweParams *RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey);


//密文全加器 
vector<MKLweSample *> MKFullAdder(const MKLweSample *a, const MKLweSample *b, const MKLweSample *ci, 
        const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) 
{
    MKLweSample *temp1 = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp2 = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp3 = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp4 = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp5 = new_MKLweSample(LWEparams, MKparams);

    MKbootsXOR_FFT_v2m2(temp1, a, b, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsAND_FFT_v2m2(temp2, temp1, ci, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsAND_FFT_v2m2(temp3, a, b, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);


    MKbootsOR_FFT_v2m2(temp4, temp2, temp3, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsXOR_FFT_v2m2(temp5, temp1, ci, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);

    vector<MKLweSample *> result;
    result.push_back(temp4);//进位
    result.push_back(temp5);//结果

    delete_MKLweSample(temp1);
    delete_MKLweSample(temp2);
    delete_MKLweSample(temp3);

    return result;//返回一个数组;
}


//密文八位全加器 
vector<MKLweSample *> MKEightAdder(MKLweSample *a[],MKLweSample *b[],const int32_t bitnum,
        const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey){

    MKLweSample *res[bitnum+1]={NULL};

    MKLweSample *ci = new_MKLweSample(LWEparams, MKparams);  //初始为0的密文
    Torus32 Add8Const = modSwitchToTorus32(-1, 8);
    MKlweNoiselessTrivial(ci, Add8Const, MKparams);

    for(int32_t i=bitnum-1;i>=0;i--){
        vector<MKLweSample *> temp=MKFullAdder(a[i],b[i],ci,bkFFT,LWEparams,extractedLWEparams,RLWEparams,MKparams,MKrlwekey);
        ci=temp[0];
        res[i+1]=temp[1];
    }
    res[0]=ci;

    vector<MKLweSample *> result;
    for(int32_t i=0;i<=bitnum;i++){
        result.push_back(res[i]);
    }
    return result;//返回9个密文组成的数组
}





//一位大于器
MKLweSample* gT(const MKLweSample *a, const MKLweSample *b, const MKLweSample *c,  
    const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
    const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) 
{
    MKLweSample *result= new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp1 = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp2= new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp3 = new_MKLweSample(LWEparams, MKparams);

    MKbootsXOR_FFT_v2m2(temp1, a, c, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsXOR_FFT_v2m2(temp2, b, c, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsAND_FFT_v2m2(temp3, temp1, temp2, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsXOR_FFT_v2m2(result, a, temp3, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);

    delete_MKLweSample(temp1);
    delete_MKLweSample(temp2);
    delete_MKLweSample(temp3);
    return result;
}

//八位大于器
MKLweSample* eightGT(MKLweSample *a[], MKLweSample *b[], const int32_t bitnum,
    const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
    const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) 
{
    MKLweSample *result=new_MKLweSample(LWEparams, MKparams);;

    MKLweSample *ci = new_MKLweSample(LWEparams, MKparams);  //初始为0的密文
    Torus32 Add8Const = modSwitchToTorus32(-1, 8);
    MKlweNoiselessTrivial(ci, Add8Const, MKparams);

    for(int32_t i=bitnum-1;i>=0;i--){
        result=gT(a[i],b[i],ci,bkFFT,LWEparams,extractedLWEparams,RLWEparams,MKparams,MKrlwekey);
        ci=result;
    }
    return result;
}


// 实现一个半加器
void simiAdder(MKLweSample *result,                         // 用于存放结果
                const MKLweSample *ca,                      // 输入加数 a
                const MKLweSample *cb,                      // 输入加数 b
                const MKLweSample *cc,                      // 输入进位 c
                const MKLweBootstrappingKeyFFT_v2 *bkFFT,   // 输入 MK 自举密钥
                const LweParams *LWEparams,                 // 输入 LWE 参数
                const LweParams *extractedLWEparams,        // 输入 LWE 拓展参数
                const TLweParams *RLWEparams,               // 输入 RLWE 参数
                const MKTFHEParams *MKparams,               // 输入 MK 参数
                const MKRLweKey *MKrlwekey)                 // 输入 MKRLWE 密钥（仅适用公钥部分）
{
    MKLweSample *temp_result = new_MKLweSample(LWEparams, MKparams);
    MKbootsXOR_FFT_v2m2(temp_result, ca, cb, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsXOR_FFT_v2m2(result, temp_result, cc, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    delete_MKLweSample(temp_result);
    return;
}


// 实现一个进位器
void simiCounter(MKLweSample *result,                       // 用于存放结果
                const MKLweSample *ca,                      // 输入加数 a
                const MKLweSample *cb,                      // 输入加数 b
                const MKLweSample *cc,                      // 输入进位 c
                const MKLweBootstrappingKeyFFT_v2 *bkFFT,   // 输入 MK 自举密钥
                const LweParams *LWEparams,                 // 输入 LWE 参数
                const LweParams *extractedLWEparams,        // 输入 LWE 拓展参数
                const TLweParams *RLWEparams,               // 输入 RLWE 参数
                const MKTFHEParams *MKparams,               // 输入 MK 参数
                const MKRLweKey *MKrlwekey)                 // 输入 MKRLWE 密钥（仅适用公钥部分）
{
    MKLweSample *temp_result2 = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp_result3 = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp_result4 = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp_result5 = new_MKLweSample(LWEparams, MKparams);
    MKbootsAND_FFT_v2m2(temp_result2, ca, cb, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsAND_FFT_v2m2(temp_result3, ca, cc, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsAND_FFT_v2m2(temp_result4, cb, cc, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsOR_FFT_v2m2(temp_result5, temp_result2, temp_result3, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsOR_FFT_v2m2(result, temp_result5, temp_result4, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    delete_MKLweSample(temp_result2);
    delete_MKLweSample(temp_result3);
    delete_MKLweSample(temp_result4);
    delete_MKLweSample(temp_result5);
    return;
}


// 实现一位加法器（全加器）
void oneBitAdder(MKLweSample *result1,                      // 用于存放结果
                MKLweSample *result2,                       // 用于存放进位
                const MKLweSample *ca,                      // 输入加数 a
                const MKLweSample *cb,                      // 输入加数 b
                const MKLweSample *cc,                      // 输入进位 c
                const MKLweBootstrappingKeyFFT_v2 *bkFFT,   // 输入 MK 自举密钥
                const LweParams *LWEparams,                 // 输入 LWE 参数
                const LweParams *extractedLWEparams,        // 输入 LWE 拓展参数
                const TLweParams *RLWEparams,               // 输入 RLWE 参数
                const MKTFHEParams *MKparams,               // 输入 MK 参数
                const MKRLweKey *MKrlwekey)                 // 输入 MKRLWE 密钥（仅适用公钥部分）
{
    simiAdder(result1, ca, cb, cc, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    simiCounter(result2, ca, cb, cc, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    return;
}


// 实现原码八位加法器
void eightBitAdder(MKLweSample *result1[],                      // 用于存放结果
                MKLweSample *result2,                       // 用于存放进位
                MKLweSample *ca[],                      // 输入加数 a
                MKLweSample *cb[],                      // 输入加数 b
                const MKLweBootstrappingKeyFFT_v2 *bkFFT,   // 输入 MK 自举密钥
                const LweParams *LWEparams,                 // 输入 LWE 参数
                const LweParams *extractedLWEparams,        // 输入 LWE 拓展参数
                const TLweParams *RLWEparams,               // 输入 RLWE 参数
                const MKTFHEParams *MKparams,               // 输入 MK 参数
                const MKRLweKey *MKrlwekey)                 // 输入 MKRLWE 密钥（仅适用公钥部分）
{
    // 构造 c0
    MKLweSample *temp_adder = new_MKLweSample(LWEparams, MKparams);
    static const Torus32 Add8Const = modSwitchToTorus32(-1, 8);
    MKlweNoiselessTrivial(temp_adder, Add8Const, MKparams);

    // 8位加法器
    for (int32_t i = 0; i < 8; i++){
        oneBitAdder(result1[i] ,temp_adder, ca[i], cb[i], temp_adder, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    }
    result2 = temp_adder;
    return;
}

// 实现补码八位加法器
void eightBitAdder_com(MKLweSample *result1[],                      // 用于存放结果
                MKLweSample *ca[],                      // 输入加数 a
                MKLweSample *cb[],                      // 输入加数 b
                const MKLweBootstrappingKeyFFT_v2 *bkFFT,   // 输入 MK 自举密钥
                const LweParams *LWEparams,                 // 输入 LWE 参数
                const LweParams *extractedLWEparams,        // 输入 LWE 拓展参数
                const TLweParams *RLWEparams,               // 输入 RLWE 参数
                const MKTFHEParams *MKparams,               // 输入 MK 参数
                const MKRLweKey *MKrlwekey)                 // 输入 MKRLWE 密钥（仅适用公钥部分）
{
    // 构造 c0
    MKLweSample *temp_adder = new_MKLweSample(LWEparams, MKparams);
    static const Torus32 Add8Const = modSwitchToTorus32(-1, 8);
    MKlweNoiselessTrivial(temp_adder, Add8Const, MKparams);

    // 8位加法器
    for (int32_t i = 7; i >= 0; i--){
        oneBitAdder(result1[i] ,temp_adder, ca[i], cb[i], temp_adder, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    }
    return;
}

// GUaker 20210526 实现原码4位乘法器
void fourBitMulter_com(MKLweSample *result1[],                      // 用于存放结果
                MKLweSample *ca[],                      // 输入乘数 a
                MKLweSample *cb[],                      // 输入乘数 b
                const MKLweBootstrappingKeyFFT_v2 *bkFFT,   // 输入 MK 自举密钥
                const LweParams *LWEparams,                 // 输入 LWE 参数
                const LweParams *extractedLWEparams,        // 输入 LWE 拓展参数
                const TLweParams *RLWEparams,               // 输入 RLWE 参数
                const MKTFHEParams *MKparams,               // 输入 MK 参数
                const MKRLweKey *MKrlwekey)                 // 输入 MKRLWE 密钥（仅使用公钥部分）
{
    // 构造 c0
    static const Torus32 Mult8Const = modSwitchToTorus32(-1, 8);
    const int32_t bitnum = 8;

    MKLweSample *temp_adder = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp_and = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp_mul[bitnum];
    for (int k = 0; k < bitnum; ++k){
        temp_mul[k] = new_MKLweSample(LWEparams, MKparams);
        MKlweNoiselessTrivial(temp_mul[k], Mult8Const, MKparams);
    }
    for (int j = bitnum - 1; j >= 0; --j){
        MKlweNoiselessTrivial(temp_adder, Mult8Const, MKparams);
        MKbootsAND_FFT_v2m2(temp_and, ca[bitnum - 1], cb[j], bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        oneBitAdder(result1[j + bitnum] ,temp_adder, temp_mul[bitnum - 1], temp_and, temp_adder, bkFFT, LWEparams, 
                        extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        for(int i = bitnum - 2; i >= 0; --i){
            MKbootsAND_FFT_v2m2(temp_and, ca[i], cb[j], bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
            oneBitAdder(temp_mul[i + 1] ,temp_adder, temp_mul[i], temp_and, temp_adder, bkFFT, LWEparams, 
                        extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        }
        temp_mul[0] = temp_adder;
    }
    for (int j = 0; j < bitnum; ++j){
        result1[j] = temp_mul[j];
    }
    return;
}

// 测试补码 8 位加法器
void test_eightBitAdder_com(){
    // 测试轮数
    const int32_t nb_trials = 10;

    // 生成参数 
    static const int32_t k = 1;
    static const double ks_stdev = 2.44e-5;// 2.44e-5; //standard deviation
    static const double bk_stdev = 3.29e-10; // 3.29e-10; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space
    static const int32_t n = 560; //500;            // LWE modulus
    static const int32_t n_extract = 1024;    // LWE extract modulus (used in bootstrapping)
    static const int32_t hLWE = 0;         // HW secret key LWE --> not used
    static const double stdevLWE = 0.012467;      // LWE ciphertexts standard deviation
    static const int32_t Bksbit = 2;       // Base bit key switching
    static const int32_t dks = 8;          // dimension key switching
    static const double stdevKS = ks_stdev; // 2.44e-5;       // KS key standard deviation
    static const int32_t N = 1024;            // RLWE,RGSW modulus
    static const int32_t hRLWE = 0;        // HW secret key RLWE,RGSW --> not used
    static const double stdevRLWEkey = bk_stdev; // 3.29e-10; // 0; // 0.012467;  // RLWE key standard deviation
    static const double stdevRLWE = bk_stdev; // 3.29e-10; // 0; // 0.012467;     // RLWE ciphertexts standard deviation
    static const double stdevRGSW = bk_stdev; // 3.29e-10;     // RGSW ciphertexts standard deviation 
    static const int32_t Bgbit = 9;        // Base bit gadget
    static const int32_t dg = 3;           // dimension gadget
    static const double stdevBK = bk_stdev; // 3.29e-10;       // BK standard deviation
    static const int32_t parties = 2;      // 参与方数量


    // new parameters 
    // 2 parties, B=2^9, d=3 -> works
    // 4 parties, B=2^8, d=4 -> works
    // 8 parties, B=2^6, d=5 -> works 

    // params
    LweParams *extractedLWEparams = new_LweParams(n_extract, ks_stdev, max_stdev);
    LweParams *LWEparams = new_LweParams(n, ks_stdev, max_stdev);
    TLweParams *RLWEparams = new_TLweParams(N, k, bk_stdev, max_stdev);
    MKTFHEParams *MKparams = new_MKTFHEParams(n, n_extract, hLWE, stdevLWE, Bksbit, dks, stdevKS, N, 
                            hRLWE, stdevRLWEkey, stdevRLWE, stdevRGSW, Bgbit, dg, stdevBK, parties);

    cout << "参数生成完成。" << endl;

    cout << "密钥生成阶段开始。" << endl;
    ofstream OutFile("Test.txt"); //利用构造函数创建txt文本，并且打开该文本
    // Key generation 
    clock_t begin_KG = clock();

    // LWE 密钥       
    MKLweKey* MKlwekey = new_MKLweKey(LWEparams, MKparams);
    MKLweKeyGen(MKlwekey);
    OutFile << "MKLWE 密钥生成完成：MKlwekey" << endl;
    for (int i = 0;i < 2; ++i){
        for (int j = 0;j < n; ++j){
            OutFile << MKlwekey->key[i].key[j];
        }
        OutFile << endl;
    }
    OutFile << endl;

    // RLWE 密钥
    MKRLweKey* MKrlwekey = new_MKRLweKey(RLWEparams, MKparams);
    MKRLweKeyGen(MKrlwekey);
    OutFile << "MKRLWE 密钥生成完成：MKrlwekey" << endl;
    for (int i = 0;i < 2; ++i){
        for (int j = 0;j < N; ++j){
            OutFile << MKrlwekey->key->key->coefs[i*N + j];
        }
        OutFile << endl;
    }
    // OutFile << "下面输出公钥：" << endl;
    // for (int i = 0;i < 2; ++i){
    //     for (int j = 0;j < N; ++j){
    //         OutFile << MKrlwekey->Pkey->coefsT[i*N + j] << endl;
    //     }
    //     OutFile << endl;
    // }
    // OutFile << endl;

    // LWE key extracted 
    MKLweKey* MKextractedlwekey = new_MKLweKey(extractedLWEparams, MKparams);
    MKtLweExtractKey(MKextractedlwekey, MKrlwekey);
    OutFile << "提取 MKLWE 密钥生成完成：MKextractedlwekey" << endl;
    for (int i = 0;i < 2; ++i){
        for (int j = 0;j < n; ++j){
            OutFile << MKextractedlwekey->key[i].key[j];
        }
        OutFile << endl;
    }
    OutFile << endl;
    OutFile.close();            //关闭Test.txt文件

    // bootstrapping + key switching keys
    MKLweBootstrappingKey_v2* MKlweBK = new_MKLweBootstrappingKey_v2(LWEparams, RLWEparams, MKparams);
    MKlweCreateBootstrappingKey_v2(MKlweBK, MKlwekey, MKrlwekey, MKextractedlwekey, 
                                extractedLWEparams, LWEparams, RLWEparams, MKparams);
    cout << "MKLWE 自举密钥生成完成：MKlweBK" << endl;

    // bootstrapping FFT + key switching keys
    MKLweBootstrappingKeyFFT_v2* MKlweBK_FFT = new_MKLweBootstrappingKeyFFT_v2(MKlweBK, LWEparams, RLWEparams, MKparams);
    cout << "MKLWE FFT 自举密钥生成完成：MKlweBK_FFT" << endl;   

    clock_t end_KG = clock();
    double time_KG = ((double) end_KG - begin_KG)/CLOCKS_PER_SEC;
    cout << "密钥生成阶段完成，耗时： " << time_KG << " 秒" << endl;


    int32_t error_count_v2m2 = 0;	// 记录错误案例个数
    double argv_time_Add8 = 0.0;	// 记录平均时间

    for (int trial = 0; trial < nb_trials; ++trial)
    {
        cout << "****************" << endl;
        cout << "第 " << trial + 1 << " 次试验：" << endl;
        cout << "****************" << endl;
        srand(time(0));

        // // 生成两个输入（原码）
        // const int32_t bitnum = 8;     // 尝试 8 位 加法器
        // int32_t mod = pow(2, bitnum);
        // int32_t addnum1 = rand() % (mod - 1);
        // int32_t addnum2 = rand() % (mod - 1);
        // int32_t answer_add8 = addnum1 + addnum2;
        // cout << "原码加数 1 = " << addnum1 << endl;
        // cout << "原码加数 2 = " << addnum2 << endl;
        // cout << "原码正确答案 = " << answer_add8 << endl;

        // // 进制转换：存二进制输入加数
        // int32_t addnumarr1[bitnum];
        // int32_t addnumarr2[bitnum];

        // // 进制转换：存二进制输出结果
        // clock_t begin_conv = clock();
        // for (int32_t i = 0; i < bitnum; i++){
        //     addnumarr1[i] = addnum1 % 2;
        //     addnumarr2[i] = addnum2 % 2;
        //     addnum1 = addnum1>>1;
        //     addnum2 = addnum2>>1;
        // }
        // // 输出二进制转换结果
        // for (int32_t i = 0; i < bitnum; i++){
        //     cout << addnumarr1[bitnum - i - 1];
        // }
        // cout << endl;
        // for (int32_t i = 0; i < bitnum; i++){
        //     cout << addnumarr2[bitnum - i - 1];
        // }
        // cout << endl;
        // clock_t end_conv = clock();
        // double time_conv = ((double) end_conv - begin_conv)/CLOCKS_PER_SEC;
        // cout << "进制转换时间： " << time_conv << "秒" << endl;

    
        // Guaker 20210519 补码加法 加数生成
        const int32_t bitnum = 8;     // 尝试 8 位 加法器
        int32_t mod = pow(2, bitnum);
        int32_t addnum1_com = rand() % (mod / 2) - mod / 4;   // 补码实验所使用的的加数1
        int32_t addnum2_com = rand() % (mod / 2) - mod / 4;   // 补码实验所使用的的加数2
        int32_t answer_add8_com = addnum1_com + addnum2_com;
        cout << "补码加数 1 = " << addnum1_com << endl;
        cout << "补码加数 2 = " << addnum2_com << endl;
        cout << "补码正确答案 = " << answer_add8_com << endl;

        // Guaker 20210519 补码加法 进制转换
        int32_t addnumarr1_com[bitnum];
        int32_t addnumarr2_com[bitnum];

        if (addnum1_com >= 0) {addnumarr1_com[0] = 0;}
        else {
            addnumarr1_com[0] = 1;
            addnum1_com += mod;
        }
        for (int32_t i = bitnum - 1; i > 0; i--){
            addnumarr1_com[i] = addnum1_com % 2;
            addnum1_com /= 2;
        }
        if (addnum2_com >= 0) {addnumarr2_com[0] = 0;}
        else {
            addnumarr2_com[0] = 1;
            addnum2_com += mod;
        }
        for (int32_t i = bitnum - 1; i > 0; i--){
            addnumarr2_com[i] = addnum2_com % 2;
            addnum2_com /= 2;
        }
        cout << "补码加数 1 二进制形式 : ";
        for (int32_t i = 0; i < bitnum; i++){
            cout <<addnumarr1_com[i];
        }
        cout << endl;
        cout << "补码加数 2 二进制形式 : ";
        for (int32_t i = 0; i < bitnum; i++){
            cout <<addnumarr2_com[i];
        }
        cout << endl;


        // 输入加数
        MKLweSample *Add_in1[bitnum];
        MKLweSample *Add_in2[bitnum];
        // 存放结果
        MKLweSample *Add_out[bitnum + 1];
        // 加密
        clock_t begin_enc = clock();
        for (int32_t i = 0; i < bitnum; i++){
            Add_in1[i] = new_MKLweSample(LWEparams, MKparams);
            Add_in2[i] = new_MKLweSample(LWEparams, MKparams);
            Add_out[i] = new_MKLweSample(LWEparams, MKparams);
            MKbootsSymEncrypt(Add_in1[i], addnumarr1_com[i], MKlwekey);
            MKbootsSymEncrypt(Add_in2[i], addnumarr2_com[i], MKlwekey);
        }
        Add_out[bitnum] = new_MKLweSample(LWEparams, MKparams);

        // 构造 c0
        // MKLweSample *temp_adder = new_MKLweSample(LWEparams, MKparams);
        // static const Torus32 Add8Const = modSwitchToTorus32(-1, 8);
        // MKlweNoiselessTrivial(temp_adder, Add8Const, MKparams);
        clock_t end_enc = clock();
        double time_enc = ((double) end_enc - begin_enc)/CLOCKS_PER_SEC;
        cout << "加密时间： " << time_enc << "秒" << endl;

        // // 8位加法器
        // clock_t begin_Add8 = clock();
        // for (int32_t i = 0; i < bitnum; i++){
        //     oneBitAdder(Add_out[i] ,temp_adder, Add_in1[i], Add_in2[i], temp_adder, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        // }
        // Add_out[bitnum] = temp_adder;

        clock_t begin_Add8 = clock();
        // eightBitAdder(Add_out, temp_adder, Add_in1, Add_in2, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        eightBitAdder_com(Add_out, Add_in1, Add_in2, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        
        clock_t end_Add8 = clock();
        double time_Add8 = ((double) end_Add8 - begin_Add8)/CLOCKS_PER_SEC;
        cout << "一个 8 位加法器运行时间： " << time_Add8 << "秒" << endl;
        argv_time_Add8 += time_Add8;

        // // 原码解密
        // clock_t begin_dec = clock();
        // int32_t addoutnum = 0;
        // for (int32_t i = 0; i < bitnum; i++){
        //     cout << MKbootsSymDecrypt(Add_out[bitnum - i], MKlwekey);
        //     addoutnum += MKbootsSymDecrypt(Add_out[bitnum - i], MKlwekey);
        //     addoutnum = addoutnum<<1;
        // }
        // cout << MKbootsSymDecrypt(Add_out[0], MKlwekey) << endl;
        // addoutnum += MKbootsSymDecrypt(Add_out[0], MKlwekey);
        // cout << "结果 = " << addoutnum << endl;

        // // 验证原码 8 位加法器
        // if (addoutnum != answer_add8) {
        //     // 计算 计算错误
        //     error_count_v2m2 +=1;
        //     cout << "ERROR!!! " << endl;
        // }
        // clock_t end_dec = clock();
        // double time_dec = ((double) end_dec - begin_dec)/CLOCKS_PER_SEC;
        // cout << "一个 8 位加法器运行时间： " << time_dec << "秒" << endl;


        // 补码解密
        clock_t begin_dec = clock();
        int32_t addoutnum = 0;
        for (int32_t i = 0; i < bitnum - 1; i++){
            cout << MKbootsSymDecrypt(Add_out[i], MKlwekey);
            addoutnum += MKbootsSymDecrypt(Add_out[i], MKlwekey);
            addoutnum = addoutnum<<1;
        }
        cout << MKbootsSymDecrypt(Add_out[bitnum - 1], MKlwekey) << endl;
        addoutnum += MKbootsSymDecrypt(Add_out[bitnum - 1], MKlwekey);
        cout << "结果 = " << addoutnum << endl;

        if (addoutnum > (mod / 2) - 1) addoutnum -= mod;

        // 验证补码 8 位加法器
        if (addoutnum != answer_add8_com) {
            // 计算 计算错误
            error_count_v2m2 +=1;
            cout << "ERROR!!! " << endl;
        }
        clock_t end_dec = clock();
        double time_dec = ((double) end_dec - begin_dec)/CLOCKS_PER_SEC;
        cout << "一个 8 位加法器运行时间： " << time_dec << "秒" << endl;

        // // TODO : 释放空间
        // for (int32_t i = 0; i < bitnum; i++){
        //     delete_MKLweSample(Add_in1[i]);
        //     delete_MKLweSample(Add_in2[i]);
        //     delete_MKLweSample(Add_out[i]);
        // }
        // delete_MKLweSample(Add_out[bitnum]);
        // delete_MKLweSample(temp_adder);
    }

    cout << endl;
    cout << "计算结果错误数量: " << error_count_v2m2 << " 错误率：" << error_count_v2m2/nb_trials << endl;
    cout << "8 位加法器 的平均时间: " << argv_time_Add8/nb_trials << " 秒" << endl;

    // 释放密钥空间
    delete_MKLweBootstrappingKeyFFT_v2(MKlweBK_FFT);
    delete_MKLweBootstrappingKey_v2(MKlweBK);
    delete_MKLweKey(MKextractedlwekey);
    delete_MKRLweKey(MKrlwekey);
    delete_MKLweKey(MKlwekey);
    // 释放参数空间
    delete_MKTFHEParams(MKparams);
    delete_TLweParams(RLWEparams);
    delete_LweParams(LWEparams);
    delete_LweParams(extractedLWEparams);

    return;
}

// 测试门电路两种实现区别
void test_boolgate(){
    // 测试轮数
    const int32_t nb_trials = 10;

    // 生成参数 
    static const int32_t k = 1;
    static const double ks_stdev = 2.44e-5;// 2.44e-5; //standard deviation
    static const double bk_stdev = 3.29e-10; // 3.29e-10; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space
    static const int32_t n = 560; //500;            // LWE modulus
    static const int32_t n_extract = 1024;    // LWE extract modulus (used in bootstrapping)
    static const int32_t hLWE = 0;         // HW secret key LWE --> not used
    static const double stdevLWE = 0.012467;      // LWE ciphertexts standard deviation
    static const int32_t Bksbit = 2;       // Base bit key switching
    static const int32_t dks = 8;          // dimension key switching
    static const double stdevKS = ks_stdev; // 2.44e-5;       // KS key standard deviation
    static const int32_t N = 1024;            // RLWE,RGSW modulus
    static const int32_t hRLWE = 0;        // HW secret key RLWE,RGSW --> not used
    static const double stdevRLWEkey = bk_stdev; // 3.29e-10; // 0; // 0.012467;  // RLWE key standard deviation
    static const double stdevRLWE = bk_stdev; // 3.29e-10; // 0; // 0.012467;     // RLWE ciphertexts standard deviation
    static const double stdevRGSW = bk_stdev; // 3.29e-10;     // RGSW ciphertexts standard deviation 
    static const int32_t Bgbit = 9;        // Base bit gadget
    static const int32_t dg = 3;           // dimension gadget
    static const double stdevBK = bk_stdev; // 3.29e-10;       // BK standard deviation
    static const int32_t parties = 2;      // 参与方数量

    // params
    LweParams *extractedLWEparams = new_LweParams(n_extract, ks_stdev, max_stdev);
    LweParams *LWEparams = new_LweParams(n, ks_stdev, max_stdev);
    TLweParams *RLWEparams = new_TLweParams(N, k, bk_stdev, max_stdev);
    MKTFHEParams *MKparams = new_MKTFHEParams(n, n_extract, hLWE, stdevLWE, Bksbit, dks, stdevKS, N, 
                            hRLWE, stdevRLWEkey, stdevRLWE, stdevRGSW, Bgbit, dg, stdevBK, parties);

    cout << "参数生成完成。" << endl;

    cout << "密钥生成阶段开始。" << endl;
    ofstream OutFile("Test.txt"); //利用构造函数创建txt文本，并且打开该文本
    // Key generation 
    clock_t begin_KG = clock();

    // LWE 密钥       
    MKLweKey* MKlwekey = new_MKLweKey(LWEparams, MKparams);
    MKLweKeyGen(MKlwekey);
    OutFile << "MKLWE 密钥生成完成：MKlwekey" << endl;
    for (int i = 0;i < parties; ++i){
        for (int j = 0;j < n; ++j){
            OutFile << MKlwekey->key[i].key[j];
        }
        OutFile << endl;
    }
    OutFile << endl;

    // RLWE 密钥
    MKRLweKey* MKrlwekey = new_MKRLweKey(RLWEparams, MKparams);
    MKRLweKeyGen(MKrlwekey);
    OutFile << "MKRLWE 密钥生成完成：MKrlwekey" << endl;
    for (int i = 0;i < parties; ++i){
        for (int j = 0;j < N; ++j){
            OutFile << MKrlwekey->key->key->coefs[i*N + j];
        }
        OutFile << endl;
    }

    // LWE key extracted 
    MKLweKey* MKextractedlwekey = new_MKLweKey(extractedLWEparams, MKparams);
    MKtLweExtractKey(MKextractedlwekey, MKrlwekey);
    OutFile << "提取 MKLWE 密钥生成完成：MKextractedlwekey" << endl;
    for (int i = 0;i < parties; ++i){
        for (int j = 0;j < n; ++j){
            OutFile << MKextractedlwekey->key[i].key[j];
        }
        OutFile << endl;
    }
    OutFile << endl;
    OutFile.close();            //关闭Test.txt文件

    // bootstrapping + key switching keys
    MKLweBootstrappingKey_v2* MKlweBK = new_MKLweBootstrappingKey_v2(LWEparams, RLWEparams, MKparams);
    MKlweCreateBootstrappingKey_v2(MKlweBK, MKlwekey, MKrlwekey, MKextractedlwekey, 
                                extractedLWEparams, LWEparams, RLWEparams, MKparams);
    cout << "MKLWE 自举密钥生成完成：MKlweBK" << endl;

    // bootstrapping FFT + key switching keys
    MKLweBootstrappingKeyFFT_v2* MKlweBK_FFT = new_MKLweBootstrappingKeyFFT_v2(MKlweBK, LWEparams, RLWEparams, MKparams);
    cout << "MKLWE FFT 自举密钥生成完成：MKlweBK_FFT" << endl;   

    clock_t end_KG = clock();
    double time_KG = ((double) end_KG - begin_KG)/CLOCKS_PER_SEC;
    cout << "密钥生成阶段完成，耗时： " << time_KG << " 秒" << endl;


    for (int trial = 0; trial < nb_trials; ++trial)
    {
        cout << "****************" << endl;
        cout << "第 " << trial + 1 << " 次试验：" << endl;
        cout << "****************" << endl;
        srand(time(0));
    
        // Guaker 20210521 输入生成
        int32_t input1 = rand() % 2;
        int32_t input2 = rand() % 2;
        cout << "输入 1 = " << input1 << endl;
        cout << "输入 2 = " << input2 << endl;
        int32_t and_ans = input1 * input2;
        int32_t or_ans = input1 + input2 - input1 * input2;
        int32_t nand_ans = 1 - and_ans;
        int32_t xor_ans = input1 + input2 - 2 * input1 * input2;
        int32_t not_ans = 1 - input1;

        // 输入加数
        MKLweSample *enc_in1;
        MKLweSample *enc_in2;
        // 存放结果
        MKLweSample *And_out1, *And_out2;
        MKLweSample *Or_out1, *Or_out2;
        MKLweSample *Nand_out;
        MKLweSample *Xor_out1, *Xor_out2;
        MKLweSample *Not_out1, *Not_out2;

        // Guaker 20210521 加密
        clock_t begin_enc = clock();
        enc_in1 = new_MKLweSample(LWEparams, MKparams);
        enc_in2 = new_MKLweSample(LWEparams, MKparams);
        MKbootsSymEncrypt(enc_in1, input1, MKlwekey);
        MKbootsSymEncrypt(enc_in2, input2, MKlwekey);
        And_out1 = new_MKLweSample(LWEparams, MKparams);
        And_out2 = new_MKLweSample(LWEparams, MKparams);
        Or_out1 = new_MKLweSample(LWEparams, MKparams);
        Or_out2 = new_MKLweSample(LWEparams, MKparams);
        Nand_out = new_MKLweSample(LWEparams, MKparams);
        Xor_out1 = new_MKLweSample(LWEparams, MKparams);
        Xor_out2 = new_MKLweSample(LWEparams, MKparams);
        Not_out1 = new_MKLweSample(LWEparams, MKparams);
        Not_out2 = new_MKLweSample(LWEparams, MKparams);
        clock_t end_enc = clock();
        double time_enc = ((double) end_enc - begin_enc)/CLOCKS_PER_SEC;
        cout << "加密时间： " << time_enc << "秒" << endl;

        // 开始拼接与门运算
        MKAND(And_out1, enc_in1, enc_in2, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        cout << "正确答案：" << and_ans << " 解密：" << MKbootsSymDecrypt(And_out1, MKlwekey) << endl;

        // 开始重写与门运算
        MKbootsAND_FFT_v2m2(And_out2, enc_in1, enc_in2, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        cout << "正确答案：" << and_ans << " 解密：" << MKbootsSymDecrypt(And_out2, MKlwekey) << endl;

        // 开始拼接或门运算
        MKOR(Or_out1, enc_in1, enc_in2, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        cout << "正确答案：" << or_ans << " 解密：" << MKbootsSymDecrypt(Or_out1, MKlwekey) << endl;

        // 开始重写或门运算
        MKbootsOR_FFT_v2m2(Or_out2, enc_in1, enc_in2, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        cout << "正确答案：" << or_ans << " 解密：" << MKbootsSymDecrypt(Or_out2, MKlwekey) << endl;

        // 开始与非门运算
        clock_t begin_MKNAND = clock();
        MKbootsNAND_FFT_v2m2(Nand_out, enc_in1, enc_in2, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        cout << "正确答案：" << nand_ans << " 解密：" << MKbootsSymDecrypt(Nand_out, MKlwekey) << endl;
        clock_t end_MKNAND = clock();
        double time_MKNAND = ((double) end_MKNAND - begin_MKNAND)/CLOCKS_PER_SEC;
        cout << "MKNAND 完成，耗时： " << time_MKNAND << " 秒" << endl;

        // 开始拼接异或门运算
        MKXOR(Xor_out1, enc_in1, enc_in2, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        cout << "正确答案：" << xor_ans << " 解密：" << MKbootsSymDecrypt(Xor_out1, MKlwekey) << endl;

        // 开始重写异或门运算
        MKbootsXOR_FFT_v2m2(Xor_out2, enc_in1, enc_in2, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        cout << "正确答案：" << xor_ans << " 解密：" << MKbootsSymDecrypt(Xor_out2, MKlwekey) << endl;

        // 开始拼接非门运算
        MKNOT(Not_out1, enc_in1, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        cout << "正确答案：" << not_ans << " 解密：" << MKbootsSymDecrypt(Not_out1, MKlwekey) << endl;

        // 开始重写非门运算
        MKlweNegate(Not_out2, enc_in1, MKparams);
        cout << "正确答案：" << not_ans << " 解密：" << MKbootsSymDecrypt(Not_out2, MKlwekey) << endl;
    }


    // 释放密钥空间
    delete_MKLweBootstrappingKeyFFT_v2(MKlweBK_FFT);
    delete_MKLweBootstrappingKey_v2(MKlweBK);
    delete_MKLweKey(MKextractedlwekey);
    delete_MKRLweKey(MKrlwekey);
    delete_MKLweKey(MKlwekey);
    // 释放参数空间
    delete_MKTFHEParams(MKparams);
    delete_TLweParams(RLWEparams);
    delete_LweParams(LWEparams);
    delete_LweParams(extractedLWEparams);

    return;
}

// Guaker 20210526 测试 4 位补码乘法器
void test_fourBitMulter(){
    // 测试轮数
    const int32_t nb_trials = 10;

    // 生成参数 
    static const int32_t k = 1;
    static const double ks_stdev = 2.44e-5;// 2.44e-5; //standard deviation
    static const double bk_stdev = 3.29e-10; // 3.29e-10; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space
    static const int32_t n = 560; //500;            // LWE modulus
    static const int32_t n_extract = 1024;    // LWE extract modulus (used in bootstrapping)
    static const int32_t hLWE = 0;         // HW secret key LWE --> not used
    static const double stdevLWE = 0.012467;      // LWE ciphertexts standard deviation
    static const int32_t Bksbit = 2;       // Base bit key switching
    static const int32_t dks = 8;          // dimension key switching
    static const double stdevKS = ks_stdev; // 2.44e-5;       // KS key standard deviation
    static const int32_t N = 1024;            // RLWE,RGSW modulus
    static const int32_t hRLWE = 0;        // HW secret key RLWE,RGSW --> not used
    static const double stdevRLWEkey = bk_stdev; // 3.29e-10; // 0; // 0.012467;  // RLWE key standard deviation
    static const double stdevRLWE = bk_stdev; // 3.29e-10; // 0; // 0.012467;     // RLWE ciphertexts standard deviation
    static const double stdevRGSW = bk_stdev; // 3.29e-10;     // RGSW ciphertexts standard deviation 
    static const int32_t Bgbit = 9;        // Base bit gadget
    static const int32_t dg = 3;           // dimension gadget
    static const double stdevBK = bk_stdev; // 3.29e-10;       // BK standard deviation
    static const int32_t parties = 2;      // 参与方数量

    const int32_t bitnum = 8;   // 尝试 4 位 乘法器
    int errorNum = 0;           // 统计错误个数
    double avgtime_enc = 0;     // 计算加密的平均时间
    double avgtime_Mult4 = 0;   // 计算乘法器的平均时间
    double avgtime_dec = 0;     // 计算解密的平均时间

    // params
    LweParams *extractedLWEparams = new_LweParams(n_extract, ks_stdev, max_stdev);
    LweParams *LWEparams = new_LweParams(n, ks_stdev, max_stdev);
    TLweParams *RLWEparams = new_TLweParams(N, k, bk_stdev, max_stdev);
    MKTFHEParams *MKparams = new_MKTFHEParams(n, n_extract, hLWE, stdevLWE, Bksbit, dks, stdevKS, N, 
                            hRLWE, stdevRLWEkey, stdevRLWE, stdevRGSW, Bgbit, dg, stdevBK, parties);

    cout << "参数生成完成。" << endl;

    cout << "密钥生成阶段开始。" << endl;
    ofstream OutFile("Test.txt"); //利用构造函数创建txt文本，并且打开该文本
    // Key generation 
    clock_t begin_KG = clock();

    // LWE 密钥       
    MKLweKey* MKlwekey = new_MKLweKey(LWEparams, MKparams);
    MKLweKeyGen(MKlwekey);
    OutFile << "MKLWE 密钥生成完成：MKlwekey" << endl;
    for (int i = 0;i < parties; ++i){
        for (int j = 0;j < n; ++j){
            OutFile << MKlwekey->key[i].key[j];
        }
        OutFile << endl;
    }
    OutFile << endl;

    // RLWE 密钥
    MKRLweKey* MKrlwekey = new_MKRLweKey(RLWEparams, MKparams);
    MKRLweKeyGen(MKrlwekey);
    OutFile << "MKRLWE 密钥生成完成：MKrlwekey" << endl;
    for (int i = 0;i < parties; ++i){
        for (int j = 0;j < N; ++j){
            OutFile << MKrlwekey->key->key->coefs[i*N + j];
        }
        OutFile << endl;
    }

    // LWE key extracted 
    MKLweKey* MKextractedlwekey = new_MKLweKey(extractedLWEparams, MKparams);
    MKtLweExtractKey(MKextractedlwekey, MKrlwekey);
    OutFile << "提取 MKLWE 密钥生成完成：MKextractedlwekey" << endl;
    for (int i = 0;i < parties; ++i){
        for (int j = 0;j < n; ++j){
            OutFile << MKextractedlwekey->key[i].key[j];
        }
        OutFile << endl;
    }
    OutFile << endl;
    OutFile.close();            //关闭Test.txt文件

    // bootstrapping + key switching keys
    MKLweBootstrappingKey_v2* MKlweBK = new_MKLweBootstrappingKey_v2(LWEparams, RLWEparams, MKparams);
    MKlweCreateBootstrappingKey_v2(MKlweBK, MKlwekey, MKrlwekey, MKextractedlwekey, 
                                extractedLWEparams, LWEparams, RLWEparams, MKparams);
    cout << "MKLWE 自举密钥生成完成：MKlweBK" << endl;

    // bootstrapping FFT + key switching keys
    MKLweBootstrappingKeyFFT_v2* MKlweBK_FFT = new_MKLweBootstrappingKeyFFT_v2(MKlweBK, LWEparams, RLWEparams, MKparams);
    cout << "MKLWE FFT 自举密钥生成完成：MKlweBK_FFT" << endl;   

    clock_t end_KG = clock();
    double time_KG = ((double) end_KG - begin_KG)/CLOCKS_PER_SEC;
    cout << "密钥生成阶段完成，耗时： " << time_KG << " 秒" << endl;


    for (int trial = 0; trial < nb_trials; ++trial)
    {
        cout << "****************" << endl;
        cout << "第 " << trial + 1 << " 次试验：" << endl;
        cout << "****************" << endl;
        srand(time(0));
    
        //Buaker 20210526 生成两个输入（原码）
        int32_t mod = pow(2, bitnum);
        int32_t multnum1 = rand() % (mod - 1);
        int32_t multnum2 = rand() % (mod - 1);
        int32_t answer_mult4 = multnum1 * multnum2;
        cout << "原码乘数 1 = " << multnum1 << endl;
        cout << "原码乘数 2 = " << multnum2 << endl;
        cout << "原码正确答案 = " << answer_mult4 << endl;

        // 进制转换：存二进制输入加数
        int32_t multnumarr1[bitnum];
        int32_t multnumarr2[bitnum];

        // 进制转换：存二进制输出结果
        clock_t begin_conv = clock();
        for (int32_t i = 0; i < bitnum; i++){
            multnumarr1[bitnum - 1 - i] = multnum1 % 2;
            multnumarr2[bitnum - 1 - i] = multnum2 % 2;
            multnum1 = multnum1>>1;
            multnum2 = multnum2>>1;
        }
        // 输出二进制转换结果
        for (int32_t i = 0; i < bitnum; i++){
            cout << multnumarr1[bitnum - i - 1];
        }
        cout << endl;
        for (int32_t i = 0; i < bitnum; i++){
            cout << multnumarr2[bitnum - i - 1];
        }
        cout << endl;
        clock_t end_conv = clock();
        double time_conv = ((double) end_conv - begin_conv)/CLOCKS_PER_SEC;
        cout << "进制转换时间： " << time_conv << "秒" << endl;


        // 输入乘数
        MKLweSample *Mult_in1[bitnum];
        MKLweSample *Mult_in2[bitnum];
        // 存放结果
        MKLweSample *Mult_out[bitnum * 2];
        // 加密
        clock_t begin_enc = clock();
        for (int32_t i = 0; i < bitnum; i++){
            Mult_in1[i] = new_MKLweSample(LWEparams, MKparams);
            Mult_in2[i] = new_MKLweSample(LWEparams, MKparams);
            MKbootsSymEncrypt(Mult_in1[i], multnumarr1[i], MKlwekey);
            MKbootsSymEncrypt(Mult_in2[i], multnumarr2[i], MKlwekey);
        }
        for (int32_t i = 0; i < bitnum * 2; i++){
            Mult_out[i] = new_MKLweSample(LWEparams, MKparams);
        }
        clock_t end_enc = clock();
        double time_enc = ((double) end_enc - begin_enc)/CLOCKS_PER_SEC;
        avgtime_enc += time_enc;
        cout << "加密运行时间： " << time_enc << "秒" << endl;

        // 4 位补码乘法
        clock_t begin_Mult4 = clock();
        fourBitMulter_com(Mult_out, Mult_in1, Mult_in1, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        // fourBitMulter_com(Mult_out, Mult_in1, Mult_in2, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, MKlwekey);
        clock_t end_Mult4 = clock();
        double time_Mult4 = ((double) end_Mult4 - begin_Mult4)/CLOCKS_PER_SEC;
        avgtime_Mult4 += time_Mult4;
        cout << "一个 4 位乘法器运行时间： " << time_Mult4 << "秒" << endl;

        // 补码解密
        clock_t begin_dec = clock();
        int32_t multoutnum = 0;
        for (int32_t i = 0; i < bitnum * 2 - 1; i++){
            cout << MKbootsSymDecrypt(Mult_out[i], MKlwekey);
            multoutnum += MKbootsSymDecrypt(Mult_out[i], MKlwekey);
            multoutnum = multoutnum<<1;
        }
        cout << MKbootsSymDecrypt(Mult_out[bitnum * 2 - 1], MKlwekey) << endl;
        multoutnum += MKbootsSymDecrypt(Mult_out[bitnum * 2 - 1], MKlwekey);
        cout << "结果 = " << multoutnum << endl;

        if (multoutnum > (mod / 2) - 1) multoutnum -= mod;

        clock_t end_dec = clock();
        double time_dec = ((double) end_dec - begin_dec)/CLOCKS_PER_SEC;
        avgtime_dec += time_dec;
        cout << "解密运行时间： " << time_dec << "秒" << endl;
    }

    cout << "正确率：" << (nb_trials - errorNum) / nb_trials << endl;
    cout << "加密的平均时间：" << avgtime_enc / nb_trials << endl;
    cout << "乘法器的平均时间：" << avgtime_Mult4 / nb_trials << endl;
    cout << "解密的平均时间：" << avgtime_dec / nb_trials << endl;

    // 释放密钥空间
    delete_MKLweBootstrappingKeyFFT_v2(MKlweBK_FFT);
    delete_MKLweBootstrappingKey_v2(MKlweBK);
    delete_MKLweKey(MKextractedlwekey);
    delete_MKRLweKey(MKrlwekey);
    delete_MKLweKey(MKlwekey);
    // 释放参数空间
    delete_MKTFHEParams(MKparams);
    delete_TLweParams(RLWEparams);
    delete_LweParams(LWEparams);
    delete_LweParams(extractedLWEparams);

    return;
}

int32_t main(int32_t argc, char **argv) {
    test_eightBitAdder_com();
    // test_fourBitMulter();
    return 0;
}