#include "posit_10_1.hpp"
#include <vector>
#include <bitset>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cassert>
using namespace std;
class Posit_10_1; 
Posit_10_1::Scalar is_zero(Posit_10_1 num){
    return (num.m_posit == 0);
}

Posit_10_1::Scalar is_inf(Posit_10_1 num){
    return (num.m_posit == 0x400);
}
Posit_10_1::Posit_10_1():m_size(10){}


Posit_10_1::Posit_10_1(const Posit_10_1& other):m_size(10){
    this->m_posit = other.m_posit;
    this->regime_bit = other.regime_bit;
    this->exponent_bit = other.exponent_bit;
    this->fraction_bit = other.fraction_bit;
    this->sign_bit = other.sign_bit;
}

Posit_10_1::Posit_10_1(const int other, const int sig_bit, const int reg_bit, const int expo_bit, const int frac_bit):m_size(10){
    this->m_posit = other;
    this->regime_bit = reg_bit;
    this->exponent_bit = expo_bit;
    this->fraction_bit = frac_bit;
    this->sign_bit = sig_bit;
}

Posit_10_1& Posit_10_1::operator=(const Posit_10_1& other){
    if (this != &other) {
        m_size = other.m_size;
        this->m_posit = other.m_posit;
        this->sign_bit = other.sign_bit;
        this->regime_bit = other.regime_bit;
        this->exponent_bit = other.exponent_bit;
        this->fraction_bit = other.fraction_bit;
    }
    return *this;
}

Posit_10_1::~Posit_10_1(){}

Posit_10_1::Scalar Posit_10_1::at(Posit_10_1::Index i){
    if ((int)i < 0 || i >= m_size) {
        throw std::out_of_range("Index out of range");
    }
    return std::bitset<10>(this->m_posit)[i];
}

Posit_10_1::Index Posit_10_1::getsize(){
    return m_size;
}

Posit_10_1::Scalar Posit_10_1::Posit_10_1_Get_Sign(){
    return sign_bit;
}

float Posit_10_1::translate_posit_10_1_to_float(){
    if((this->m_posit&0x3ff) == 0){
        return (float)0;
    }
    Posit_10_1::Scalar m_sign = this->Posit_10_1_Get_Sign();      
    float m_reg = 0, m_frac = 0, m_expo = 0;
    int regime_mi = 0, exponent_mi = 0;
    int fraction_bit = this->fraction_bit, exponent_bit = this->exponent_bit, regime_bit = this->regime_bit - 1, current_posit_num = (this->m_posit) ;
    if(m_sign == 1){/*negative*/
        current_posit_num = ((~current_posit_num) + 1);
    }
    current_posit_num &= 0x1ff;
    //cout<<"erjinzhi=="<<std::bitset<15>(current_posit_num)<<endl;
    std::bitset<9> binary_arr(current_posit_num);
    for(int i = 0; i < fraction_bit; i++){
        m_frac += std::pow(2,i - fraction_bit) * binary_arr[i];
    }
    //cout << "m_frac = " << m_frac << endl;
    m_frac += 1;
    for(int i = fraction_bit;i < exponent_bit + fraction_bit; i++){
        exponent_mi += std::pow(2, i -fraction_bit) * binary_arr[i];
    }
    if(binary_arr[exponent_bit+fraction_bit] == 0){/*m-1*/
        regime_mi = regime_bit - 1;
    }
    else if(binary_arr[exponent_bit+fraction_bit] == 1){/*-m*/
        regime_mi = - regime_bit ;
    }
    m_expo = std::pow(2, exponent_mi);
    //std::cout << "expo_mi = " << exponent_mi << endl;
    m_reg = std::pow(useed, regime_mi);
    //std::cout << "m_reg = " << regime_mi << endl;
    float m_val = std::pow((-1),m_sign) * m_reg * m_expo * m_frac;
    return m_val;
}

Posit_10_1 translate_float_to_posit_10_1(const float num){
    if(num == 0){
        return zero;
    }
    int signbit = (num >= 0) ? 0 : 1;
    float tmp = std::abs(num);
    /*binary representation of float*/
    int binarypresentation = *reinterpret_cast<int*>(&tmp); /*binary representation of float so we can use bit operation to simplify our code*/
    /*bit mast*/
 //   std::cout <<"IEEE_REPRESENT ===" <<bitset<32>(binarypresentation) << std::endl;
    int exponentmask = 0x7f800000;
    int mantissamask = 0x007fffff;
    /*sign bit*/
 //   std::cout <<"signbit===" <<bitset<1>(signbit) << std::endl;
    /*real exponent*/
    int ori_exponent = ((binarypresentation & exponentmask) >> 23) -127; 
 //   std::cout << "ori_exponent==="<<bitset<32>(ori_exponent) << std::endl;
 //   printf("ori_expo===%d\n",ori_exponent);
    /*reduce the difficulty in the process of transforming k to m*/
    bool regemi_sign_bit = (ori_exponent >= 0) ? 0 : 1;
    /*mantissa*/
    unsigned int ori_mantissa = binarypresentation & mantissamask;          
 //   std::cout <<"MANTISSA===" <<bitset<23>(ori_mantissa) << std::endl;
    /*exponent*/
    int exponent = ori_exponent & 0x1; /*eg 36-> 4  that is the regime k is then we need to calculate the m...... 4 that is the exponent need to represent*/
    /*regime*/
 //   std::cout <<"EXPONENT===" <<bitset<3>(exponent) << std::endl;
    int regime_bit = std::abs((ori_exponent >> 1)); 
 //   std::cout <<"EXPONENT===" <<bitset<32>((ori_exponent >> 3)) << std::endl;
    int regime_bit_with = 1, exponent_bit = 0, fraction_bit = 0;
    int m, final_posit;
    unsigned int regime = 1;
   // printf("regemi_sign_bit ==%d\n",regemi_sign_bit);
   // printf("regime_bit====%d\n",regime_bit);
    /*calculate regime representation, because that is the only different between posit and IEEE754*/
    if(regemi_sign_bit == 1) {/*add 0*/
        m = regime_bit;
    }
    else if(regemi_sign_bit == 0){/*add 1*/

        m = regime_bit + 1;
        int tmp = m;
        while(tmp != 1){
            regime = (regime << 1) + 1;
            tmp--;
        }
        regime = regime << 1;
    }
        
   // printf("m ==%d\n",m);
    //printf("regemi ==%d\n",regime);
    //std::cout << bitset<32>(regime) << std::endl;
    /*calculate exponent and fraction*/
    int remain_bit = 9 - m - regime_bit_with;
     //   printf("remian_bit ==%d\n",remain_bit);
    if(remain_bit >= 1){
        exponent_bit = 1;
        fraction_bit = remain_bit - exponent_bit;
    }
    else if(remain_bit < 1){
        exponent_bit = remain_bit;
        fraction_bit = 0;
    }
    int reversemask = 0x1ff;
    /*get the final posit number*/
   // printf("exponentfinal ==%d\n",exponent_bit);
   // printf("fractionfinal ==%d\n",fraction_bit);
   // puts("111111111111111111111");
    if(signbit == 0){
        final_posit = (signbit << 9)  
                            + (regime << (remain_bit))
                            + ((exponent >> (1 - exponent_bit)) << fraction_bit)
                                + (ori_mantissa >> (23-fraction_bit));
    }
    else if(signbit == 1) {

        final_posit = (signbit << 9)  
                        + (((regime << (remain_bit))
                            + ((exponent >> (1 - exponent_bit)) << fraction_bit)
                                + (ori_mantissa >> (23-fraction_bit)))^reversemask)+1;
    }
    //std::cout << bitset<16>(((((~(regime<<(32-m)))>>(32-m)) << (remain_bit)))) << std::endl;
    //std::cout << bitset<16>(final_posit) << std::endl;
    //std::cout << bitset<16>((ori_mantissa >> (23-fraction_bit))) << std::endl;
    Posit_10_1 mmm_posit;
    mmm_posit.m_posit = final_posit;
    mmm_posit.sign_bit = signbit;
    mmm_posit.regime_bit = m + regime_bit_with;
    mmm_posit.exponent_bit = exponent_bit;
    mmm_posit.fraction_bit = fraction_bit;

    return mmm_posit;
}

Posit_10_1 Posit_10_1::operator+(const Posit_10_1& other) const{
    Posit_10_1 ltmp = *this;
    Posit_10_1 rtmp = other;
    if(is_zero(ltmp) || is_zero(rtmp)){
        return zero;
    }
    float num1 = ltmp.translate_posit_10_1_to_float();
    float num2 = rtmp.translate_posit_10_1_to_float();
    float num3 = num1 + num2;
    return translate_float_to_posit_10_1(num3);
}
Posit_10_1 Posit_10_1::operator+=(const Posit_10_1& other) {
    Posit_10_1 tmp = *this + other;
    *this = tmp;
    return *this;
}
Posit_10_1 Posit_10_1::operator-(const Posit_10_1& other) const{
    Posit_10_1 ltmp = *this;
    Posit_10_1 rtmp = other;
    if(is_zero(ltmp) || is_zero(rtmp)){
        return zero;
    }
    float num1 = ltmp.translate_posit_10_1_to_float();
    float num2 = rtmp.translate_posit_10_1_to_float();
    float num3 = num1 - num2;
    return translate_float_to_posit_10_1(num3);
}
Posit_10_1 Posit_10_1::operator-=(const Posit_10_1& other) {
    Posit_10_1 tmp = *this - other;
    *this = tmp;
    return *this;
}


Posit_10_1 Posit_10_1::operator*(const Posit_10_1& other) const{
    Posit_10_1 ltmp = *this;
    Posit_10_1 rtmp = other;
    if(is_zero(ltmp) || is_zero(rtmp)){
        return zero;
    }
    float num1 = ltmp.translate_posit_10_1_to_float();
    float num2 = rtmp.translate_posit_10_1_to_float();
    float num3 = num1*num2;
    return translate_float_to_posit_10_1(num3);
}

Posit_10_1 Posit_10_1::operator*=(const Posit_10_1& other) {
    Posit_10_1 tmp = *this * other;
    *this = tmp;
    return *this;
}

Posit_10_1 Posit_10_1::operator/(const Posit_10_1& other) const{
    if(!is_zero(*this) && is_zero(other)){ // the division of any nonzero number by zero is inf
        return inf;
    }
    else if(is_zero(*this) && !is_zero(other)){
        return zero;
    }
    else if(is_zero(*this) && is_zero(other)){
        throw std::invalid_argument("zero cannot be divided by zero.");
    }
    Posit_10_1 ltmp = *this;
    Posit_10_1 rtmp = other;
    if(is_zero(ltmp) || is_zero(rtmp)){
        return zero;
    }
    float num1 = ltmp.translate_posit_10_1_to_float();
    float num2 = rtmp.translate_posit_10_1_to_float();
    float num3 = num1/num2;
    return translate_float_to_posit_10_1(num3);
}

Posit_10_1 Posit_10_1::operator/=(const Posit_10_1& other) {
    Posit_10_1 tmp = *this / other;
    *this = tmp;
    return *this;
}


bool Posit_10_1::operator<(const Posit_10_1& other) const{
    Posit_10_1 ltmp = *this;
    Posit_10_1 rtmp = other;
    float left = ltmp.translate_posit_10_1_to_float();
    float right = rtmp.translate_posit_10_1_to_float();
    return left < right;
}
bool Posit_10_1::operator<=(const Posit_10_1& other) const{
    Posit_10_1 ltmp = *this;
    Posit_10_1 rtmp = other;
    if (ltmp == rtmp) {
        return 1;
    }
    float left = ltmp.translate_posit_10_1_to_float();
    float right = rtmp.translate_posit_10_1_to_float();
    return left <= right;
}
bool Posit_10_1::operator>(const Posit_10_1& other) const{
    Posit_10_1 ltmp = *this;
    Posit_10_1 rtmp = other;
    float left = ltmp.translate_posit_10_1_to_float();
    float right = rtmp.translate_posit_10_1_to_float();
    return left > right;
}
bool Posit_10_1::operator>=(const Posit_10_1& other) const{
    Posit_10_1 ltmp = *this;
    Posit_10_1 rtmp = other;
    if (ltmp == rtmp) {
        return 1;
    }
    float left = ltmp.translate_posit_10_1_to_float();
    float right = rtmp.translate_posit_10_1_to_float();
    return left >= right;
}
bool Posit_10_1::operator==(const Posit_10_1& other) const{
    return this->m_posit == other.m_posit;
}
bool Posit_10_1::operator!=(const Posit_10_1& other) const{
    return this->m_posit != other.m_posit;
}

Posit_10_1& Posit_10_1::operator++(){
    *this += one;
    return *this;
}
Posit_10_1 Posit_10_1::operator++(int){
    Posit_10_1 tmp(*this);
    *this += one;
    return tmp;
}
Posit_10_1& Posit_10_1::operator--(){
    *this -= one;
    return *this;
}
Posit_10_1 Posit_10_1::operator--(int){
    Posit_10_1 tmp(*this);
    *this -= one;
    return tmp;
}

