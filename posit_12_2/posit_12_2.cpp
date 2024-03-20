#include "posit_12_2.hpp"
#include <vector>
#include <bitset>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;
class Posit_12_2; 
Posit_12_2::Scalar is_zero(Posit_12_2 num){
    return (num.m_posit == 0);
}

Posit_12_2::Scalar is_inf(Posit_12_2 num){
    return (num.m_posit == 0x800);
}
Posit_12_2::Posit_12_2():m_size(12){}


Posit_12_2::Posit_12_2(const Posit_12_2& other):m_size(12){
    this->m_posit = other.m_posit;
    this->regime_bit = other.regime_bit;
    this->exponent_bit = other.exponent_bit;
    this->fraction_bit = other.fraction_bit;
    this->sign_bit = other.sign_bit;
}

Posit_12_2::Posit_12_2(const int other, const int sig_bit, const int reg_bit, const int expo_bit, const int frac_bit):m_size(12){
    this->m_posit = other;
    this->regime_bit = reg_bit;
    this->exponent_bit = expo_bit;
    this->fraction_bit = frac_bit;
    this->sign_bit = sig_bit;
}

Posit_12_2& Posit_12_2::operator=(const Posit_12_2& other){
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

Posit_12_2::~Posit_12_2(){}

Posit_12_2::Scalar Posit_12_2::at(Posit_12_2::Index i){
    if ((int)i < 0 || i >= m_size) {
        throw std::out_of_range("Index out of range");
    }
    return std::bitset<12>(this->m_posit)[i];
}

Posit_12_2::Index Posit_12_2::getsize(){
    return m_size;
}

Posit_12_2::Scalar Posit_12_2::Posit_12_2_Get_Sign(){
    return sign_bit;
}

float Posit_12_2::translate_posit_12_2_to_float(){
    if((this->m_posit&0xfff) == 0){
        return (float)0;
    }
    Posit_12_2::Scalar m_sign = this->Posit_12_2_Get_Sign();      
    float m_reg = 0, m_frac = 0, m_expo = 0;
    int regime_mi = 0, exponent_mi = 0;
    int fraction_bit = this->fraction_bit, exponent_bit = this->exponent_bit, regime_bit = this->regime_bit - 1, current_posit_num = (this->m_posit) ;
    if(m_sign == 1){/*negative*/
        current_posit_num = ((~current_posit_num) + 1);
    }
    current_posit_num &= 0x7ff;
    //cout<<"erjinzhi=="<<std::bitset<15>(current_posit_num)<<endl;
    std::bitset<11> binary_arr(current_posit_num);
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

Posit_12_2 translate_float_to_posit_12_2(const float num){
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
    int exponent = ori_exponent & 0x3; /*eg 36-> 9  that is the regime k is then we need to calculate the m...... 0 that is the exponent need to represent*/
    /*regime*/
 //   std::cout <<"EXPONENT===" <<bitset<3>(exponent) << std::endl;
    int regime_bit = std::abs((ori_exponent >> 2)); 
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
    int remain_bit = 11 - m - regime_bit_with;
     //   printf("remian_bit ==%d\n",remain_bit);
    if(remain_bit >= 2){
        exponent_bit = 2;
        fraction_bit = remain_bit - exponent_bit;
    }
    else if(remain_bit < 2){
        exponent_bit = remain_bit;
        fraction_bit = 0;
    }
    int reversemask = 0x7ff;
    /*get the final posit number*/
   // printf("exponentfinal ==%d\n",exponent_bit);
   // printf("fractionfinal ==%d\n",fraction_bit);
   // puts("111111111111111111111");
    if(signbit == 0){
        final_posit = (signbit << 11)  
                            + (regime << (remain_bit))
                            + ((exponent >> (2 - exponent_bit)) << fraction_bit)
                                + (ori_mantissa >> (23-fraction_bit));
    }
    else if(signbit == 1) {

        final_posit = (signbit << 11)  
                        + (((regime << (remain_bit))
                            + ((exponent >> (2 - exponent_bit)) << fraction_bit)
                                + (ori_mantissa >> (23-fraction_bit)))^reversemask)+1;
    }
    //std::cout << bitset<16>(((((~(regime<<(32-m)))>>(32-m)) << (remain_bit)))) << std::endl;
    //std::cout << bitset<16>(final_posit) << std::endl;
    //std::cout << bitset<16>((ori_mantissa >> (23-fraction_bit))) << std::endl;
    Posit_12_2 mmm_posit;
    mmm_posit.m_posit = final_posit;
    mmm_posit.sign_bit = signbit;
    mmm_posit.regime_bit = m + regime_bit_with;
    mmm_posit.exponent_bit = exponent_bit;
    mmm_posit.fraction_bit = fraction_bit;

    return mmm_posit;
}

Posit_12_2 Posit_12_2::operator+(const Posit_12_2& other) const{
    Posit_12_2 ltmp = *this;
    Posit_12_2 rtmp = other;
    if(is_zero(ltmp) || is_zero(rtmp)){
        return zero;
    }
    float num1 = ltmp.translate_posit_12_2_to_float();
    float num2 = rtmp.translate_posit_12_2_to_float();
    float num3 = num1 + num2;
    return translate_float_to_posit_12_2(num3);
}
Posit_12_2 Posit_12_2::operator+=(const Posit_12_2& other) {
    Posit_12_2 tmp = *this + other;
    *this = tmp;
    return *this;
}
Posit_12_2 Posit_12_2::operator-(const Posit_12_2& other) const{
    Posit_12_2 ltmp = *this;
    Posit_12_2 rtmp = other;
    if(is_zero(ltmp) || is_zero(rtmp)){
        return zero;
    }
    float num1 = ltmp.translate_posit_12_2_to_float();
    float num2 = rtmp.translate_posit_12_2_to_float();
    float num3 = num1 - num2;
    return translate_float_to_posit_12_2(num3);
}
Posit_12_2 Posit_12_2::operator-=(const Posit_12_2& other) {
    Posit_12_2 tmp = *this - other;
    *this = tmp;
    return *this;
}


Posit_12_2 Posit_12_2::operator*(const Posit_12_2& other) const{
    Posit_12_2 ltmp = *this;
    Posit_12_2 rtmp = other;
    if(is_zero(ltmp) || is_zero(rtmp)){
        return zero;
    }
    float num1 = ltmp.translate_posit_12_2_to_float();
    float num2 = rtmp.translate_posit_12_2_to_float();
    float num3 = num1*num2;
    return translate_float_to_posit_12_2(num3);
}

Posit_12_2 Posit_12_2::operator*=(const Posit_12_2& other) {
    Posit_12_2 tmp = *this * other;
    *this = tmp;
    return *this;
}

Posit_12_2 Posit_12_2::operator/(const Posit_12_2& other) const{
    if(!is_zero(*this) && is_zero(other)){ // the division of any nonzero number by zero is inf
        return inf;
    }
    else if(is_zero(*this) && !is_zero(other)){
        return zero;
    }
    else if(is_zero(*this) && is_zero(other)){
        throw std::invalid_argument("zero cannot be divided by zero.");
    }
    Posit_12_2 ltmp = *this;
    Posit_12_2 rtmp = other;
    if(is_zero(ltmp) || is_zero(rtmp)){
        return zero;
    }
    float num1 = ltmp.translate_posit_12_2_to_float();
    float num2 = rtmp.translate_posit_12_2_to_float();
    float num3 = num1/num2;
    return translate_float_to_posit_12_2(num3);
}

Posit_12_2 Posit_12_2::operator/=(const Posit_12_2& other) {
    Posit_12_2 tmp = *this / other;
    *this = tmp;
    return *this;
}


bool Posit_12_2::operator<(const Posit_12_2& other) const{
    Posit_12_2 ltmp = *this;
    Posit_12_2 rtmp = other;
    float left = ltmp.translate_posit_12_2_to_float();
    float right = rtmp.translate_posit_12_2_to_float();
    return left < right;
}
bool Posit_12_2::operator<=(const Posit_12_2& other) const{
    Posit_12_2 ltmp = *this;
    Posit_12_2 rtmp = other;
    if (ltmp == rtmp) {
        return 1;
    }
    float left = ltmp.translate_posit_12_2_to_float();
    float right = rtmp.translate_posit_12_2_to_float();
    return left <= right;
}
bool Posit_12_2::operator>(const Posit_12_2& other) const{
    Posit_12_2 ltmp = *this;
    Posit_12_2 rtmp = other;
    float left = ltmp.translate_posit_12_2_to_float();
    float right = rtmp.translate_posit_12_2_to_float();
    return left > right;
}
bool Posit_12_2::operator>=(const Posit_12_2& other) const{
    Posit_12_2 ltmp = *this;
    Posit_12_2 rtmp = other;
    if (ltmp == rtmp) {
        return 1;
    }
    float left = ltmp.translate_posit_12_2_to_float();
    float right = rtmp.translate_posit_12_2_to_float();
    return left >= right;
}
bool Posit_12_2::operator==(const Posit_12_2& other) const{
    return this->m_posit == other.m_posit;
}
bool Posit_12_2::operator!=(const Posit_12_2& other) const{
    return this->m_posit != other.m_posit;
}

Posit_12_2& Posit_12_2::operator++(){
    *this += one;
    return *this;
}
Posit_12_2 Posit_12_2::operator++(int){
    Posit_12_2 tmp(*this);
    *this += one;
    return tmp;
}
Posit_12_2& Posit_12_2::operator--(){
    *this -= one;
    return *this;
}
Posit_12_2 Posit_12_2::operator--(int){
    Posit_12_2 tmp(*this);
    *this -= one;
    return tmp;
}


int main(){
    /*const std::vector<Posit_12_2::Scalar> m_Posit_12_2_vector1 = {0,0,1,0, 1,1,0,1, 1,1,0,0, 0,1,0,1};
    const std::vector<Posit_12_2::Scalar> m_Posit_12_2_vector2 = {0,0,1,0, 1,1,0,1, 1,0,0,1, 0,1,1,0};
    Posit_12_2 m_Posit_12_2_1(0x2DC5, 0, 2, 3, 10);
    Posit_12_2 m_Posit_12_2_2(0x2D96, 0, 2, 3, 10); */
    /*============================================normal test==========================================================================*/
    /*for(int i = 0;i < 16; i++){
        std::cout << m_Posit_12_2_1.at(i) ;
    }
    std::cout << " " << std::endl;
    std::vector<Posit_12_2::Scalar> m_reg = m_Posit_12_2_1.Posit_12_2_Get_Reg();
    for(int i = 0;i < m_reg.size(); i++){
        std::cout << m_reg[i] ;
    }
    std::cout << " " << std::endl;
    std::cout << "sign = " << m_Posit_12_2_1.Posit_12_2_Get_Sign() << std::endl;
    std::cout << "reg_size = " << m_Posit_12_2_1.Posit_12_2_Get_Reg_Size() << std::endl;
    std::cout << "reg_k = " << m_Posit_12_2_1.Posit_12_2_Get_Reg_k() << std::endl;
    std::cout << "reg_val = " << m_Posit_12_2_1.Posit_12_2_Get_Reg_Val() << std::endl;
    std::vector<Posit_12_2::Scalar> m_expo = m_Posit_12_2_1.Posit_12_2_Get_Expo();
    for(int i = 0;i < m_expo.size(); i++){
        std::cout << m_expo[i] ;
    }
    std::cout << " " << std::endl;
    std::cout << "expo_val = " << m_Posit_12_2_1.Posit_12_2_Get_Expo_Val() << std::endl;
    std::vector<Posit_12_2::Scalar> m_frac = m_Posit_12_2_1.Posit_12_2_Get_Frac();
    for(int i = 0;i < m_frac.size(); i++){
        std::cout << m_frac[i] ;
    }
    std::cout << " " << std::endl;
    std::cout << "frac_val = " << m_Posit_12_2_1.Posit_12_2_Get_Frac_Val() << std::endl;
    printf("%.20f\n", m_Posit_12_2_1.Posit_12_2_Get_Frac_Val());
    std::cout << "float_val = " << m_Posit_12_2_1.translate_Posit_12_2_to_float() << std::endl;
    printf("%.30f\n", m_Posit_12_2_1.translate_Posit_12_2_to_float());
    //Posit_12_2 new_posit(translate_float_to_Posit_12_2(0.000003553926944732666015625000));
    Posit_12_2 new_posit(translate_float_to_Posit_12_2(0.0000000000007815970093));
    for(int i = 0;i < new_posit.getsize();i++){
        std::cout << new_posit.at(i) ;
    }
    std::cout << " " << std::endl;
    float error1 = translate_float_to_Posit_12_2(0.0000000000007815970093).translate_Posit_12_2_to_float() - 0.0000000000007815970093;
    std::cout << "error1 = " << error1 << std::endl;
    float num2 = std::pow(256,-6)*std::pow(2,4) *(1+(float)22/32);
    Posit_12_2 error2 = translate_float_to_Posit_12_2(num2);
    for(int i = 0;i < error2.getsize();i++){
        std::cout << error2.at(i) ;
    }
    std::cout << " " << std::endl;*/


    /*======================================multiply test===========================================*/
   /* float a = m_Posit_12_2_1.translate_Posit_12_2_to_float();
    printf("a = %.30f\n",a);
    float b = m_Posit_12_2_2.translate_Posit_12_2_to_float();
    printf("b = %.30f\n",b);
    float c = a /b;
    printf("%.30f\n", c);
    for(int i = 15;i >= 0;i--){
        std::cout << translate_float_to_Posit_12_2(c).at(i);
    }

    
    std::cout << " " << std::endl;
    Posit_12_2 multi_posit = m_Posit_12_2_1 / m_Posit_12_2_2;

    std::cout << "result = " << std::bitset<16>(multi_posit.m_posit) << std::endl;
    std::cout << "result = " << multi_posit.translate_Posit_12_2_to_float() << std::endl;*/
    
    
    
    
    // printf("%.30f\n",translate_float_to_Posit_12_2(0).translate_Posit_12_2_to_float());
    // float a= 2.784336;
    // float b = -0.6576;
    // float c = a * b;
    // Posit_12_2 p1= translate_float_to_Posit_12_2(a);
    // Posit_12_2 p2= translate_float_to_Posit_12_2(b);
    // Posit_12_2 p3 = p1 * p2;
    // float f1 = p1.translate_Posit_12_2_to_float();
    // float f2 = p2.translate_Posit_12_2_to_float();
    // float f3 = f1 * f2;
    // printf("c = %.30f\n",c);
    // std::cout << p1.m_posit << endl;
    // printf("p1 = %.30f\n",translate_float_to_Posit_12_2(a).translate_Posit_12_2_to_float());
    // std::cout << "p1 = " << std::bitset<32>(p1.m_posit) << std::endl;
    // std::cout << "p2 = " << std::bitset<32>(p2.m_posit) << std::endl;
    // std::cout << "p3 = " << std::bitset<32>(p3.m_posit) << std::endl;
    // printf("p3 = %.30f\n",p3.translate_Posit_12_2_to_float());
    // printf("f1 = %.30f\n",f1);
    // printf("f2 = %.30f\n",f2);
    // printf("f3 = %.30f\n",f3);

    float f_value, step = 1e-6;
    int cnt = 0;
    for(f_value = -10.0f; f_value < 10.0f; f_value += step){
        Posit_12_2 p_value = translate_float_to_posit_12_2(f_value);
        float f_value_from_p = p_value.translate_posit_12_2_to_float();
        float delta = std::abs(f_value - f_value_from_p);
        if(delta >= 1e-2){
            puts("===============================");
            printf("Error at %d * %f\n", cnt, step);
            printf("F: %f\n", f_value);
            printf("P: %f\n", f_value_from_p);
            //break;
        }
        cnt++;
    }



    /*for(int i = 0;i < multi_posit.getsize(); i++){
        std::cout << multi_posit.at(i) ;
    }
    std::cout << " " << std::endl;
    for(int i = 0;i < multi_posit.Posit_12_2_Get_Reg().size(); i++){
        std::cout << multi_posit.Posit_12_2_Get_Reg()[i] ;
    }
    std::cout << " " << std::endl;
    for(int i = 0;i < multi_posit.Posit_12_2_Get_Expo().size(); i++){
        std::cout << multi_posit.Posit_12_2_Get_Expo()[i] ;
    }
    std::cout << " " << std::endl;
    for(int i = 0;i < multi_posit.Posit_12_2_Get_Frac().size(); i++){
        std::cout << multi_posit.Posit_12_2_Get_Frac()[i] ;
    }
    std::cout << " " << std::endl;
    printf("%.30f\n", multi_posit.translate_Posit_12_2_to_float());*/

    /*============================abs and pow test===========================*/
    /*printf("%.30f\n", m_Posit_12_2_1.translate_Posit_12_2_to_float());
    for(Posit_12_2::Index i = 0;i < m_Posit_12_2_1.getsize();i++){
        std::cout << m_Posit_12_2_1.at(i);
    }
    std::cout << " " << std::endl;
    Posit_12_2 m_abs_1 = m_Posit_12_2_1.abs();
    for(Posit_12_2::Index i = 0;i < m_abs_1.getsize();i++){
        std::cout << m_abs_1.at(i);
    }
    std::cout << " " << std::endl;
    printf("%.30f\n", m_abs_1.translate_Posit_12_2_to_float());
    printf("%.30f\n", m_Posit_12_2_2.translate_Posit_12_2_to_float());
    Posit_12_2 m_abs_2 = m_Posit_12_2_2.abs();
    for(Posit_12_2::Index i = 0;i < m_abs_2.getsize();i++){
        std::cout << m_abs_2.at(i);
    }
    std::cout << " " << std::endl;
    printf("%.30f\n", m_abs_2.translate_Posit_12_2_to_float());*/


    /*printf("%.30f\n", std::pow(256,14));
    float m_float_1 = std::pow(-3907584,4);
    printf("%.30f\n", m_float_1);
    Posit_12_2 m_pow_1 = m_Posit_12_2_1.Posit_12_2::pow(4);
    for(Posit_12_2::Index i = 0;i < m_pow_1.getsize();i++){
        std::cout << m_pow_1.at(i);
    }
    std::cout << " " << std::endl;
    printf("%.30f\n", m_pow_1.translate_Posit_12_2_to_float());
    printf("%.30f\n",(m_float_1 - m_pow_1.translate_Posit_12_2_to_float())/m_pow_1.translate_Posit_12_2_to_float());
    
    float m_float_2 = std::pow(m_Posit_12_2_2.translate_Posit_12_2_to_float(),4);
    printf("%.30f\n", m_float_2);
    Posit_12_2 m_pow_2 = m_Posit_12_2_2.Posit_12_2::pow(4);
    for(Posit_12_2::Index i = 0;i < m_pow_2.getsize();i++){
        std::cout << m_pow_2.at(i);
    }
    std::cout << " " << std::endl;
    printf("%.30f\n", m_pow_2.translate_Posit_12_2_to_float());
    printf("%.30f\n",(m_float_2 - m_pow_2.translate_Posit_12_2_to_float())/m_pow_2.translate_Posit_12_2_to_float());
    */

    /*
    
    for(int j = 0;j < 10;j++){
        Posit_12_2 m_pow_1 = Three.Posit_12_2::pow(j);
        cout << "j = " << j << endl;
        cout << "std::pow(3,j) = " << std::pow(3,j) << endl;
        for(Posit_12_2::Index i = 0;i < m_pow_1.getsize();i++){
            std::cout << m_pow_1.at(i);
        }
        std::cout << " " << std::endl;
        printf("%.10f\n",m_pow_1.translate_Posit_12_2_to_float());
    }
    puts("=============================");
    
    for(int j = 0;j < 10;j++){
        Posit_12_2 m_pow_1 = nge_Three.Posit_12_2::pow(j);
        cout << "j = " << j << endl;
        cout << "std::pow(-3,j) = " << std::pow(-3,j) << endl;
        for(Posit_12_2::Index i = 0;i < m_pow_1.getsize();i++){
            std::cout << m_pow_1.at(i);
        }
        std::cout << " " << std::endl;
        printf("%.10f\n",m_pow_1.translate_Posit_12_2_to_float());
    }
    puts("=============================");
    
    for(int j = 0;j < 10;j++){
        Posit_12_2 m_pow_1 = Two_point_five.Posit_12_2::pow(j);
        cout << "j = " << j << endl;
        cout << "std::pow(2.5,j) = " << std::pow(2.5,j) << endl;
        for(Posit_12_2::Index i = 0;i < m_pow_1.getsize();i++){
            std::cout << m_pow_1.at(i);
        }
        std::cout << " " << std::endl;
        printf("%.10f\n",m_pow_1.translate_Posit_12_2_to_float());
    }*/

    /*===============================add sub test======================*/
    /*float m_float = m_Posit_12_2_1.translate_Posit_12_2_to_float() +  m_Posit_12_2_2.translate_Posit_12_2_to_float();
    printf("%.30f\n", m_float);
    Posit_12_2 m_add_1 = Three - negative_three;
    for(Posit_12_2::Index i = 0;i < m_add_1.getsize();i++){
        std::cout << m_add_1.at(i);
    }
    std::cout << " " << std::endl;
    printf("%.30f\n", m_add_1.translate_Posit_12_2_to_float());
    printf("%.30f\n", (++m_add_1).translate_Posit_12_2_to_float());
    printf("%.30f\n", (m_add_1++).translate_Posit_12_2_to_float());
    printf("%.30f\n", (m_add_1--).translate_Posit_12_2_to_float());
    printf("%.30f\n", (--m_add_1).translate_Posit_12_2_to_float());*/

    return 0;
}