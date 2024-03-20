#include "float.h"
#include <vector>
#include <bitset>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;
class Float; 
//Float:: Float():m_float(0b0),m_size(32){}
const Float zero({ 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });//0
const Float one({ 0,  0, 1, 1, 1, 1, 1, 1, 1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });//1
const Float two({ 0,  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });//2
const Float epsilon1({ 0,  0, 1, 1, 1, 0, 0, 1, 0,  1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1 });//0.0002
const Float epsilon2({ 0,  0, 1, 1, 0, 1, 0, 0, 0,  1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1 });//0.0000002
long translate_vector_to_long(std::vector<Float::Scalar> vec) {
    Float::Index size = vec.size();
    long result = 0;
    for (Float::Index i = 0; i < size; i++) {
        result += (std::pow(2, size - 1 - i) * vec[i]);
    }
    return result;
}

double translate_vector_to_double(std::vector<Float::Scalar> vec) {
    Float::Index size = vec.size();
    double result = 0;
    for (Float::Index i = 0; i < size; i++) {
        int power = -i - 1;
        result += (std::pow(2, power) * vec[i]);
    }
    return result;
}

Float::Scalar is_zero(Float num) {
    for (Float::Index i = 0; i < num.getsize(); i++) {
        if (num.at(i)) {
            return 0;
        }
    }
    return 1;
}

Float:: Float(Index size) :m_size(size) {
    m_float.resize(m_size);
    
    for (std::size_t i = 0; i < m_size; i++) {
        m_float.push_back(0);
    }

}

Float:: Float(const Float& other) :m_size(other.m_size) {
    m_float.resize(m_size);
    for (std::size_t i = 0; i < m_size; i++) {
        m_float[i] = other.m_float[i];
    }
}

Float:: Float(const std::vector<Float::Scalar> other) :m_size(other.size()) {
    m_float.resize(m_size);
    for (std::size_t i = 0; i < m_size; i++) {
        m_float[i] = other[i];
    }
}

Float& Float::operator=(const Float& other) {
    if (this != &other) {
        {
            std::vector<Scalar> temp;
            m_float.swap(temp);
        }
        m_size = other.m_size;
        m_float.resize(m_size);

        for (std::size_t i = 0; i < m_size; i++) {
            m_float[i] = other.m_float[i];
        }

    }
    return *this;
}

Float::~Float() {
    {
        std::vector<Scalar> temp;
        m_float.swap(temp);
    }
}

Float::Scalar Float::at(Index i) {
    if (i < 0 || i >= m_size) {
        throw std::out_of_range("Index out of range");
    }
    return m_float[i];
}

const Float::Scalar Float::at(Index i) const {
    if (i < 0 || i >= m_size) {
        throw std::out_of_range("Index out of range");
    }
    return m_float[i];
}

Float::Index Float::getsize() {
    return m_size;
}

std::vector<Float::Scalar> Float::block(Index i, Index j) const {
    if (i < 0 || i > m_size || j < 0 || j > m_size) {
        throw std::invalid_argument("Invalid block range: Out of range");
    }
    if (i > j) {
        throw std::invalid_argument("Invalid block range: Index i > j");
    }

    std::vector<Float::Scalar> m_block(j - i + 1);
    for (Index k = i; k <= j; ++k) {
        m_block[k - i] = m_float[k];
    }

    return m_block;
}

Float::Scalar Float::Float_Get_Sign() {
    return m_float[0];
}
std::vector<Float::Scalar> Float::Float_Get_Expo() {
    Index m_expo_size = 0;
    std::vector<Scalar> m_expo;
    switch(m_size){
        case 32: {
            m_expo_size = 8;
            m_expo.resize(8);
            for (Index i = 1; i < 1 + m_expo_size; i++) {
                m_expo[i-1] = m_float[i];
            }
            break;
        }
        case 64: {
            m_expo_size = 11;
            m_expo.resize(11);
            for (Index i = 1; i < 1 + m_expo_size; i++) {
                m_expo[i-1] = m_float[i];
            }
            break;
        }
        case 80: {
            m_expo_size = 15;
            m_expo.resize(15);
            for (Index i = 1; i < 1 + m_expo_size; i++) {
                m_expo[i-1] = m_float[i];
            }
            break;
        }
        default:
            break;
    }
    return m_expo;
}
std::vector<Float::Scalar> Float::Float_Get_Manti() {
    Index m_manti_size = 0;
    std::vector<Scalar> m_manti;
    switch (m_size) {
    case 32: {
        m_manti_size = 23;
        m_manti.resize(23);
        for (Index i = 9; i < 9 + m_manti_size; i++) {
            m_manti[i-9] = m_float[i];
        }
        break;
       }
    case 64: {
        m_manti_size = 52;
        m_manti.resize(52);
            for (Index i = 12; i < 12 + m_manti_size; i++) {
                m_manti[i-12] = m_float[i];
            }
            break;
        }
    case 80: {
            m_manti_size = 64;
            m_manti.resize(64);
            for (Index i = 16; i < 16 + m_manti_size; i++) {
                m_manti[i-16] = m_float[i];
            }
            break;
        }
    default:
        break;
    }
    
    return m_manti;
}

double Float::translate_float_to_double() {
    double result = 0;
    double sign = std::pow(- 1, this->Float_Get_Sign());
    std::vector<bool> my_e = this->Float_Get_Expo();
    std::vector<bool> my_m = this->Float_Get_Manti();
    double e = 0, m = 0;
    switch (m_size) {

    case 32: {
        for (Index i = 0; i < 8; i++) {
            e += my_e[i] * (std::pow(2, (7 - i)));
        }
        ///std::cout << "e =" << e << std::endl;
        e = std::pow(2, int((e - 127)));
        ///std::cout << "e =" << e << std::endl;
        for (int i = 0; i < 23; i++) {
            int power = -i - 1;
            ///std::cout << "my_m[" << i << "]=" << my_m[i] << std::endl;
            ///std::cout << power << std::endl;
            ///std::cout << std::pow(2, power) << std::endl;
            m += my_m[i] * (std::pow(2, power));
            ///std::cout << "m =" << m << std::endl;
        }
        ///std::cout << "m =" << m << std::endl;
        m += 1;
        result = sign * e * m;
        break;
        }
    case 64: {
        for (int i = 0; i < 11; i++) {
            e += my_e[i] * (std::pow(2, (10 - i)));
        }
        e = std::pow(2, int((e - 1023)));
        for (int i = 0; i < 52; i++) {
            int power = -i - 1;
            m += my_m[i] * (std::pow(2, power));
        }
        m += 1;
        result = sign * e * m;
        break;
    }
    default: 
        break;
    }
    return result;

}

Float Float::translate_double_to_float32(const long integ,const float dec) {//integer and decimal
    std::vector<Float::Scalar> result_vector(32);
    long integer = std::abs(integ);
    long remainder = 0;//the remainder of integer part
    float decimal = dec;
    Float::Index expo_size = 0;//size of expo part
    Float::Index manti_size = 0;//size of manti part
    unsigned int e = 0;
    if (integ >= 0) {//sign bit
        result_vector[0] = 0;
    }
    else {
        result_vector[0] = 1;
    }
    std::vector<long> expo_manti;
    
    while (integer > 0) {
        remainder = integer % 2;
        integer /= 2;
        expo_manti.push_back(remainder);
    }
    std::reverse(expo_manti.begin(), expo_manti.end());//get the correct expo part
    expo_size = expo_manti.size();
    //printf("%d", expo_size);
    //expo_manti.push_back('.');
    while (decimal != 0 && manti_size < 64) {
        decimal *= 2;
        //std::cout << decimal << std::endl;
        if (decimal >= 1.0f) {
            
            decimal -= 1;
            expo_manti.push_back(1);
        }
        else {
            expo_manti.push_back(0);
        }
        manti_size += 1;
    }
    while (manti_size < 64) {
        expo_manti.push_back(0);
        manti_size++;
    }
    //std::cout << expo_manti.size() << std::endl;
    //for (int i = 0; i < expo_manti.size(); i++) {
    //    std::cout << "expo_manti[" << i << "]=" << expo_manti[i] << std::endl;
    //}
    if (expo_size > 1) {
        e = expo_size + 126;//expoλ
        std::bitset<8> expo_bit{ e };
        //std::cout << expo_bit << std::endl;
        for (int ii = 1; ii <= 8; ii++) {
            result_vector[ii] = expo_bit[8 - ii];
        }
        for (int jj = 1; jj <= 23; jj++) {
            result_vector[jj + 8] = expo_manti[jj];
        }
    }
    else if(expo_size == 1) {
        e = 127;//expoλ
        std::bitset<8> expo_bit{ e };
        //std::cout << expo_bit << std::endl;
        for (int ii = 1; ii <= 8; ii++) {
            result_vector[ii] = expo_bit[8 - ii];
        }
        for (int jj = 1; jj <= 23; jj++) {
            result_vector[jj + 8] = expo_manti[jj];
        }
        //std::cout << expo_bit << std::endl;
    }
    else if (expo_size == 0) {
        auto it = std::find(expo_manti.begin(), expo_manti.end(), 1);
        expo_size = std::distance(expo_manti.begin(), it);
        //std::cout << expo_size << std::endl;
        e = (-1) * expo_size + 126;
        std::bitset<8> expo_bit{ e };
        for (int ii = 1; ii <= 8; ii++) {
            result_vector[ii] = expo_bit[8 - ii];
        }
        for (int jj = 1; jj <= 23; jj++) {
            result_vector[jj + 8] = expo_manti[jj + expo_size];
        }
        //std::cout << expo_bit << std::endl;
    }
    Float result(result_vector);
    return result;
}

Float Float::translate_double_to_float64(const long integ, const float dec) {
    std::vector<Float::Scalar> result_vector(64);
    long integer = std::abs(integ);
    long remainder = 0;//the remainder of integer part
    float decimal = dec;
    Float::Index expo_size = 0;//size of expo part
    Float::Index manti_size = 0;//size of manti part
    unsigned int e = 0;
    if (integ >= 0) {//sign bit
        result_vector[0] = 0;
    }
    else {
        result_vector[0] = 1;
    }
    std::vector<long> expo_manti;

    while (integer > 0) {
        remainder = integer % 2;
        integer /= 2;
        expo_manti.push_back(remainder);
    }
    std::reverse(expo_manti.begin(), expo_manti.end());//get the correct expo part
    expo_size = expo_manti.size();
    //printf("%d", expo_size);
    //expo_manti.push_back('.');
    while (decimal != 0 && manti_size < 128) {
        decimal *= 2;
        //std::cout << decimal << std::endl;
        if (decimal >= 1.0f) {

            decimal -= 1;
            expo_manti.push_back(1);
        }
        else {
            expo_manti.push_back(0);
        }
        manti_size += 1;
    }
    while (manti_size < 128) {
        expo_manti.push_back(0);
        manti_size++;
    }
    //std::cout << expo_manti.size() << std::endl;
    //for (int i = 0; i < expo_manti.size(); i++) {
    //    std::cout << "expo_manti[" << i << "]=" << expo_manti[i] << std::endl;
    //}
    if (expo_size > 1) {
        e = expo_size + 1022;//expoλ
        std::bitset<11> expo_bit{ e };
        //std::cout << expo_bit << std::endl;
        for (int ii = 1; ii <= 11; ii++) {
            result_vector[ii] = expo_bit[11 - ii];
        }
        for (int jj = 1; jj <= 52; jj++) {
            result_vector[jj + 11] = expo_manti[jj];
        }
    }
    else if (expo_size == 1) {
        e = 1023;//expoλ
        std::bitset<11> expo_bit{ e };
        //std::cout << expo_bit << std::endl;
        for (int ii = 1; ii <= 11; ii++) {
            result_vector[ii] = expo_bit[11 - ii];
        }
        for (int jj = 1; jj <= 52; jj++) {
            result_vector[jj + 11] = expo_manti[jj];
        }
    }
    else if (expo_size == 0) {
        auto it = std::find(expo_manti.begin(), expo_manti.end(), 1);
        expo_size = std::distance(expo_manti.begin(), it);
        //std::cout << expo_size << std::endl;
        e = (-1) * expo_size + 1022;
        std::bitset<11> expo_bit{ e };
        for (int ii = 1; ii <= 11; ii++) {
            result_vector[ii] = expo_bit[11 - ii];
        }
        for (int jj = 1; jj <= 52; jj++) {
            result_vector[jj + 11] = expo_manti[jj + expo_size];
        }
    }
    Float result(result_vector);
    return result;
}



Float Float::operator+(const Float& other) const {
    std::vector<Float::Scalar> result_vector;
    Float ltmp = *this;
    Float rtmp = other;
    std::vector<Float::Scalar> lexpo = ltmp.Float_Get_Expo();
    std::vector<Float::Scalar> rexpo = rtmp.Float_Get_Expo();
    std::vector<Float::Scalar> lmanti = ltmp.Float_Get_Manti();
    std::vector<Float::Scalar> rmanti = rtmp.Float_Get_Manti();
    Float::Scalar lsign = ltmp.Float_Get_Sign();
    Float::Scalar rsign = rtmp.Float_Get_Sign();
    if (lexpo == rexpo && lmanti == rmanti && lsign != rsign) {
        for (Float::Index i = 0; i < this->m_size; i++) {
            result_vector.push_back(0);
        }
        Float result(result_vector);
        return result;
    }
    switch (this->m_size) {
    case 32: {
        
        long l_e = translate_vector_to_long(lexpo) - 127;
        long r_e = translate_vector_to_long(rexpo) - 127;
        double l_m, r_m;
        long delta_e = std::abs(l_e - r_e);
        //std::cout << "delta_e=" << delta_e << std::endl;
        if (l_e < r_e) {
            l_m = (1 + translate_vector_to_double(lmanti)) * std::pow(2, l_e-r_e);
            r_m = (1 + translate_vector_to_double(rmanti));
            l_e += delta_e;
        }
        else {
            l_m = (1 + translate_vector_to_double(lmanti));
            r_m = (1 + translate_vector_to_double(rmanti)) * std::pow(2, r_e-l_e);
            r_e += delta_e;
        }
        //std::cout << "l_e2=" << l_e << std::endl;
        //std::cout << "r_e2=" << r_e << std::endl;
        //std::cout << "l_m2=" << l_m << std::endl;
        //std::cout << "r_m2=" << r_m << std::endl;
        double result_manti = std::pow((-1), lsign) * l_m + std::pow((-1), rsign) * r_m;
        Float::Scalar result_sign = (result_manti < 0);
        result_vector.push_back(result_sign);
        signed int result_expo = std::max(l_e, r_e);
        //std::cout << "result_expo=" << result_expo << std::endl;
        result_manti = std::abs(result_manti);
        int bit_cnt = 0;
        //std::cout << "result_manti=" << result_manti << std::endl;
        while (floor(result_manti) != 1) {//find out how many bits to change,and get final manti
            if (result_manti < 1) {
                result_manti *= 2;
                bit_cnt--;
            }
            else {
                result_manti /= 2;
                bit_cnt++;
            }
        }
        //std::cout << "result_manti2=" << result_manti << std::endl;
        //std::cout << "bit_cnt=" << bit_cnt << std::endl;
        result_manti -= 1;
        unsigned int result_expo2 = result_expo+(bit_cnt + 127);//final exponent
        if (result_expo2 > 255 ) {//check overflow
            throw std::out_of_range("Overflow!");
        }
        if (result_expo2 < 0) {//check underflow
            throw std::out_of_range("Underflow!");
        }
        //std::cout << "result_expo2=" << result_expo2 << std::endl;
        std::bitset<8> final_expo{ result_expo2 };
        for (int ii = 0; ii < 8; ii++) {//push to result (expo part)
            result_vector.push_back(final_expo[7 - ii]);
        }
        Float::Index manti_size = 0;//manti���ֳ���
        while (result_manti != 0 && manti_size < 23) {//push to result (manti part)
            result_manti *= 2;
            if (result_manti >= 1.0f) {

                result_manti -= 1;
                result_vector.push_back(1);
            }
            else {
                result_vector.push_back(0);
            }
            manti_size += 1;
        }
        while (manti_size < 23) {
            result_vector.push_back(0);
            manti_size++;
        }
        break;
    }
    case 64:{
        long l_e = translate_vector_to_long(lexpo) - 1023; /*look expo difference*/
        long r_e = translate_vector_to_long(rexpo) - 1023;
        double l_m, r_m;
        long delta_e = std::abs(l_e - r_e);
        if (l_e < r_e) {
            l_m = (1 + translate_vector_to_double(lmanti)) * std::pow(2, l_e-r_e);  /*align*/
            r_m = (1 + translate_vector_to_double(rmanti));
            l_e += delta_e;
        }
        else {
            l_m = (1 + translate_vector_to_double(lmanti));
            r_m = (1 + translate_vector_to_double(rmanti)) * std::pow(2, r_e-l_e);
            r_e += delta_e;
        }
        double result_manti = std::pow((-1), lsign) * l_m + std::pow((-1), rsign) * r_m;
        Float::Scalar result_sign = (result_manti < 0);
        result_vector.push_back(result_sign);
        signed int result_expo = std::max(l_e, r_e);
        //std::cout << "result_expo=" << result_expo << std::endl;
        result_manti = std::abs(result_manti);
        int bit_cnt = 0;
        while (floor(result_manti) != 1) {//find out how many bits to change,and get final manti
            if (result_manti < 1) {
                result_manti *= 2;
                bit_cnt--;
            }
            else {
                result_manti /= 2;
                bit_cnt++;
            }
        }
        //std::cout << "result_manti2=" << result_manti << std::endl;
        //std::cout << "bit_cnt=" << bit_cnt << std::endl;
        result_manti -= 1;                    //now itis 2
        unsigned int result_expo2 = result_expo+(bit_cnt + 1023);//final exponent
        if (result_expo2 > 2047 ) {//check overflow
            throw std::out_of_range("Overflow!");
        }
        if (result_expo2 < 0) {//check underflow
            throw std::out_of_range("Underflow!");
        }
        //std::cout << "result_expo2=" << result_expo2 << std::endl;
        std::bitset<8> final_expo{ result_expo2 };
        for (int ii = 0; ii < 11; ii++) {//push to result (expo part)
            result_vector.push_back(final_expo[10 - ii]);
        }
        Float::Index manti_size = 0;//size of manti part
        while (result_manti != 0 && manti_size < 52) {//push to result (manti part)
            result_manti *= 2;
            if (result_manti >= 1.0f) {

                result_manti -= 1;
                result_vector.push_back(1);
            }
            else {
                result_vector.push_back(0);
            }
            manti_size += 1;
        }
        while (manti_size < 52) {
            result_vector.push_back(0);
            manti_size++;
        }
        break;
    }
    default:
        break;
    }
    Float result(result_vector);
    return result;
}
Float Float::operator+=(const Float& other) {
    Float tmp = *this + other;
    *this = tmp;
    return *this;
}
Float Float::operator-(const Float& other) const {
    std::vector<Float::Scalar> result_vector;
    Float ltmp = *this;
    Float rtmp = other;
    if (ltmp == rtmp) {
        for (Float::Index i = 0; i < this->m_size; i++) {
            result_vector.push_back(0);
        }
        Float result(result_vector);
        return result;
    }
    std::vector<Float::Scalar> lexpo = ltmp.Float_Get_Expo();
    std::vector<Float::Scalar> rexpo = rtmp.Float_Get_Expo();
    std::vector<Float::Scalar> lmanti = ltmp.Float_Get_Manti();
    std::vector<Float::Scalar> rmanti = rtmp.Float_Get_Manti();
    Float::Scalar lsign = ltmp.Float_Get_Sign();
    Float::Scalar rsign = rtmp.Float_Get_Sign();
    switch (this->m_size) {
    case 32: {
        long l_e = translate_vector_to_long(lexpo) - 127;
        long r_e = translate_vector_to_long(rexpo) - 127;
        double l_m, r_m;
        long delta_e = std::abs(l_e - r_e);
        if (l_e < r_e) {
            l_m = (1 + translate_vector_to_double(lmanti)) * std::pow(2, l_e-r_e);
            r_m = (1 + translate_vector_to_double(rmanti));
            l_e += delta_e;
        }
        else {
            l_m = (1 + translate_vector_to_double(lmanti));
            r_m = (1 + translate_vector_to_double(rmanti)) * std::pow(2, r_e-l_e);
            r_e += delta_e;
        }
        double result_manti = std::pow((-1), lsign) * l_m - std::pow((-1), rsign) * r_m;
        Float::Scalar result_sign = (result_manti < 0);
        result_vector.push_back(result_sign);
        signed int result_expo = std::max(l_e, r_e);
        result_manti = std::abs(result_manti);
        int bit_cnt = 0;
        while (floor(result_manti) != 1) {//find out how many bits to change,and get final manti
            if (result_manti < 1) {
                result_manti *= 2;
                bit_cnt--;
            }
            else {
                result_manti /= 2;
                bit_cnt++;
            }
        }
        result_manti -= 1;
        unsigned int result_expo2 = result_expo + (bit_cnt + 127);//final exponent
        if (result_expo2 > 255 ) {//check overflow
            throw std::out_of_range("Overflow!");
        }
        if (result_expo2 < 0) {//check underflow
            throw std::out_of_range("Underflow!");
        }
        std::bitset<8> final_expo{ result_expo2 };
        for (int ii = 0; ii < 8; ii++) {//push to result (expo part)
            result_vector.push_back(final_expo[7 - ii]);
        }
        Float::Index manti_size = 0;//size of manti part
        while (result_manti != 0 && manti_size < 23) {//push to result (manti part)
            result_manti *= 2;
            if (result_manti >= 1.0f) {

                result_manti -= 1;
                result_vector.push_back(1);
            }
            else {
                result_vector.push_back(0);
            }
            manti_size += 1;
        }
        while (manti_size < 23) {
            result_vector.push_back(0);
            manti_size++;
        }
        break;
    }
    case 64: {
        long l_e = translate_vector_to_long(lexpo) - 1023;
        long r_e = translate_vector_to_long(rexpo) - 1023;
        double l_m, r_m;
        long delta_e = std::abs(l_e - r_e);
        if (l_e < r_e) {
            l_m = (1 + translate_vector_to_double(lmanti)) * std::pow(2, l_e-r_e);
            r_m = (1 + translate_vector_to_double(rmanti));
            l_e += delta_e;
        }
        else {
            l_m = (1 + translate_vector_to_double(lmanti));
            r_m = (1 + translate_vector_to_double(rmanti)) * std::pow(2, r_e-l_e);
            r_e += delta_e;
        }
        double result_manti = std::pow((-1), lsign) * l_m - std::pow((-1), rsign) * r_m;
        Float::Scalar result_sign = (result_manti < 0);
        result_vector.push_back(result_sign);
        signed int result_expo = std::max(l_e, r_e);
        result_manti = std::abs(result_manti);
        int bit_cnt = 0;
        while (floor(result_manti) != 1) {//find out how many bits to change,and get final manti
            if (result_manti < 1) {
                result_manti *= 2;
                bit_cnt--;
            }
            else {
                result_manti /= 2;
                bit_cnt++;
            }
        }
        result_manti -= 1;
        unsigned int result_expo2 = result_expo + (bit_cnt + 1023);//final exponent
        if (result_expo2 > 2047 ) {//check overflow
            throw std::out_of_range("Overflow!");
        }
        if (result_expo2 < 0) {//check underflow
            throw std::out_of_range("Underflow!");
        }
        std::bitset<8> final_expo{ result_expo2 };
        for (int ii = 0; ii < 11; ii++) {//push to result (expo part)
            result_vector.push_back(final_expo[10 - ii]);
        }
        Float::Index manti_size = 0;//manti���ֳ���
        while (result_manti != 0 && manti_size < 52) {//push to result (manti part)
            result_manti *= 2;
            if (result_manti >= 1.0f) {

                result_manti -= 1;
                result_vector.push_back(1);
            }
            else {
                result_vector.push_back(0);
            }
            manti_size += 1;
        }
        while (manti_size < 52) {
            result_vector.push_back(0);
            manti_size++;
        }
        break;
    }
    default:
        break;
    }
    Float result(result_vector);
    return result;
}
Float Float::operator-=(const Float& other) {
    Float tmp = *this - other;
    *this = tmp;
    return *this;
}
Float Float::operator*(const Float& other) const {
    std::vector<Float::Scalar> result_vector;
    Float ltmp = *this;
    Float rtmp = other;
    if (is_zero(ltmp) || is_zero(rtmp)) {
        for (Float::Index i = 0; i < this->m_size; i++) {
            result_vector.push_back(0);
        }
        Float result(result_vector);
        return result;
    }
    std::vector<Float::Scalar> lexpo = ltmp.Float_Get_Expo();
    std::vector<Float::Scalar> rexpo = rtmp.Float_Get_Expo();
    std::vector<Float::Scalar> lmanti = ltmp.Float_Get_Manti();
    std::vector<Float::Scalar> rmanti = rtmp.Float_Get_Manti();
    Float::Scalar lsign = ltmp.Float_Get_Sign();
    Float::Scalar rsign = rtmp.Float_Get_Sign();
    switch (this->m_size) {
    case 32: {
        Float::Scalar result_sign = lsign ^ rsign;
        //std::cout << "result_sign = " << result_sign << std::endl;
        long l_e = translate_vector_to_long(lexpo) - 127;
        long r_e = translate_vector_to_long(rexpo) - 127;
        //std::cout << "l_e = " << l_e << std::endl;
        //std::cout << "r_e = " << r_e << std::endl;
        long result_expo = l_e + r_e;
        //std::cout << "result_expo = " << result_expo << std::endl;
        double l_m = (1 + translate_vector_to_double(lmanti));
        double r_m = (1 + translate_vector_to_double(rmanti));
        //std::cout << "l_m = " << l_m << std::endl;
        //std::cout << "r_m = " << r_m << std::endl;
        double result_manti = l_m * r_m;
        //std::cout << "result_manti = " << result_manti << std::endl;
        int bit_cnt = 0;
        while (floor(result_manti) != 1) {//find out how many bits to change,and get final manti
            result_manti /= 2;
            bit_cnt++;
        }
        //std::cout << "bit_cnt = " << bit_cnt << std::endl;
        result_manti -= 1;
        //std::cout << "result_manti2 = " << result_manti << std::endl;
        unsigned int result_expo2 = result_expo + (bit_cnt + 127);//final exponent
        //std::cout << "result_expo2 = " << result_expo2 << std::endl;
        if (result_expo2 > 255 ) {//check overflow
            throw std::out_of_range("Overflow!");
        }
        if (result_expo2 < 0) {//check underflow
            throw std::out_of_range("Underflow!");
        }
        result_vector.push_back(result_sign);
        std::bitset<8> final_expo{ result_expo2 };
        for (int ii = 0; ii < 8; ii++) {//push to result (expo part)
            result_vector.push_back(final_expo[7 - ii]);
        }
        Float::Index manti_size = 0;//size of manti part
        while (result_manti != 0 && manti_size < 23) {//push to result (manti part)
            result_manti *= 2;
            if (result_manti >= 1.0f) {

                result_manti -= 1;
                result_vector.push_back(1);
            }
            else {
                result_vector.push_back(0);
            }
            manti_size += 1;
        }
        while (manti_size < 23) {
            result_vector.push_back(0);
            manti_size++;
        }
        break;
    }
    case 64: {
        Float::Scalar result_sign = lsign ^ rsign;
        //std::cout << "result_sign = " << result_sign << std::endl;
        long l_e = translate_vector_to_long(lexpo) - 1023;
        long r_e = translate_vector_to_long(rexpo) - 1023;
        //std::cout << "l_e = " << l_e << std::endl;
        //std::cout << "r_e = " << r_e << std::endl;
        long result_expo = l_e + r_e;
        //std::cout << "result_expo = " << result_expo << std::endl;
        double l_m = (1 + translate_vector_to_double(lmanti));
        double r_m = (1 + translate_vector_to_double(rmanti));
        //std::cout << "l_m = " << l_m << std::endl;
        //std::cout << "r_m = " << r_m << std::endl;
        double result_manti = l_m * r_m;
        //std::cout << "result_manti = " << result_manti << std::endl;
        int bit_cnt = 0;
        while (floor(result_manti) != 1) {//find out how many bits to change,and get final manti
            result_manti /= 2;
            bit_cnt++;
        }
        //std::cout << "bit_cnt = " << bit_cnt << std::endl;
        result_manti -= 1;
        //std::cout << "result_manti2 = " << result_manti << std::endl;
        unsigned int result_expo2 = result_expo + (bit_cnt + 1023);//final exponent
        //std::cout << "result_expo2 = " << result_expo2 << std::endl;
        if (result_expo2 > 2047 ) {//check overflow
            throw std::out_of_range("Overflow!");
        }
        if (result_expo2 < 0) {//check underflow
            throw std::out_of_range("Underflow!");
        }
        result_vector.push_back(result_sign);
        std::bitset<8> final_expo{ result_expo2 };
        for (int ii = 0; ii < 11; ii++) {//push to result (expo part)
            result_vector.push_back(final_expo[10 - ii]);
        }
        Float::Index manti_size = 0;//size of manti part
        while (result_manti != 0 && manti_size < 52) {//push to result (manti part)
            result_manti *= 2;
            if (result_manti >= 1.0f) {

                result_manti -= 1;
                result_vector.push_back(1);
            }
            else {
                result_vector.push_back(0);
            }
            manti_size += 1;
        }
        while (manti_size < 52) {
            result_vector.push_back(0);
            manti_size++;
        }
        break;
    }
    default:
        break;
    }
    Float result(result_vector);
    return result;
}
Float Float::operator*=(const Float& other) {
    Float tmp = *this * other;
    *this = tmp;
    return *this;
}
Float Float::operator/(const Float& other) const {
    std::vector<Float::Scalar> result_vector;
    Float ltmp = *this;
    Float rtmp = other;
    if (is_zero(ltmp)) {//0 be divided
        for (Float::Index i = 0; i < this->m_size; i++) {
            result_vector.push_back(0);
        }
        Float result(result_vector);
        return result;
    }
    if (is_zero(rtmp)) {//divide 0
        throw std::invalid_argument("Invalid divisor:cannot be divided by 0");
    }
    std::vector<Float::Scalar> lexpo = ltmp.Float_Get_Expo();
    std::vector<Float::Scalar> rexpo = rtmp.Float_Get_Expo();
    std::vector<Float::Scalar> lmanti = ltmp.Float_Get_Manti();
    std::vector<Float::Scalar> rmanti = rtmp.Float_Get_Manti();
    Float::Scalar lsign = ltmp.Float_Get_Sign();
    Float::Scalar rsign = rtmp.Float_Get_Sign();
    switch (this->m_size) {
    case 32: {
        Float::Scalar result_sign = lsign ^ rsign;
        //std::cout << "result_sign = " << result_sign << std::endl;
        long l_e = translate_vector_to_long(lexpo) - 127;
        long r_e = translate_vector_to_long(rexpo) - 127;
        //std::cout << "l_e = " << l_e << std::endl;
        //std::cout << "r_e = " << r_e << std::endl;
        long result_expo = l_e - r_e;
        //std::cout << "result_expo = " << result_expo << std::endl;
        double l_m = (1 + translate_vector_to_double(lmanti));
        double r_m = (1 + translate_vector_to_double(rmanti));
        //std::cout << "l_m = " << l_m << std::endl;
        //std::cout << "r_m = " << r_m << std::endl;
        double result_manti = l_m / r_m;
        //std::cout << "result_manti = " << result_manti << std::endl;
        int bit_cnt = 0;
        while (floor(result_manti) != 1) {//find out how many bits to change,and get final manti
            if (result_manti < 1) {
                result_manti *= 2;
                bit_cnt--;
            }
            else {
                result_manti /= 2;
                bit_cnt++;
            }
        }
        //std::cout << "bit_cnt = " << bit_cnt << std::endl;
        result_manti -= 1;
        //std::cout << "result_manti2 = " << result_manti << std::endl;
        unsigned int result_expo2 = result_expo + (bit_cnt + 127);//final exponent
        //std::cout << "result_expo2 = " << result_expo2 << std::endl;
        if (result_expo2 > 255) {//check overflow
            throw std::out_of_range("Overflow!");
        }
        if (result_expo2 < 0) {//check underflow
            throw std::out_of_range("Underflow!");
        }
        result_vector.push_back(result_sign);
        std::bitset<8> final_expo{ result_expo2 };
        for (int ii = 0; ii < 8; ii++) {//push to result (expo part)
            result_vector.push_back(final_expo[7 - ii]);
        }
        Float::Index manti_size = 0;//size of manti part
        while (result_manti != 0 && manti_size < 23) {//push to result (manti part)
            result_manti *= 2;
            if (result_manti >= 1.0f) {

                result_manti -= 1;
                result_vector.push_back(1);
            }
            else {
                result_vector.push_back(0);
            }
            manti_size += 1;
        }
        while (manti_size < 23) {
            result_vector.push_back(0);
            manti_size++;
        }
        break;
    }
    case 64: {
        Float::Scalar result_sign = lsign ^ rsign;
        //std::cout << "result_sign = " << result_sign << std::endl;
        long l_e = translate_vector_to_long(lexpo) - 1023;
        long r_e = translate_vector_to_long(rexpo) - 1023;
        //std::cout << "l_e = " << l_e << std::endl;
        //std::cout << "r_e = " << r_e << std::endl;
        long result_expo = l_e - r_e;
        //std::cout << "result_expo = " << result_expo << std::endl;
        double l_m = (1 + translate_vector_to_double(lmanti));
        double r_m = (1 + translate_vector_to_double(rmanti));
        //std::cout << "l_m = " << l_m << std::endl;
        //std::cout << "r_m = " << r_m << std::endl;
        double result_manti = l_m / r_m;
        //std::cout << "result_manti = " << result_manti << std::endl;
        int bit_cnt = 0;
        while (floor(result_manti) != 1) {//find out how many bits to change,and get final manti
            if (result_manti < 1) {
                result_manti *= 2;
                bit_cnt--;
            }
            else {
                result_manti /= 2;
                bit_cnt++;
            }
        }
        //std::cout << "bit_cnt = " << bit_cnt << std::endl;
        result_manti -= 1;
        //std::cout << "result_manti2 = " << result_manti << std::endl;
        unsigned int result_expo2 = result_expo + (bit_cnt + 1023);//final exponent
        //std::cout << "result_expo2 = " << result_expo2 << std::endl;
        if (result_expo2 > 2047) {//check overflow
            throw std::out_of_range("Overflow!");
        }
        if (result_expo2 < 0) {//check underflow
            throw std::out_of_range("Underflow!");
        }
        result_vector.push_back(result_sign);
        std::bitset<8> final_expo{ result_expo2 };
        for (int ii = 0; ii < 11; ii++) {//push to result (expo part)
            result_vector.push_back(final_expo[10 - ii]);
        }
        Float::Index manti_size = 0;//size of manti part
        while (result_manti != 0 && manti_size < 52) {//push to result (manti part)
            result_manti *= 2;
            if (result_manti >= 1.0f) {

                result_manti -= 1;
                result_vector.push_back(1);
            }
            else {
                result_vector.push_back(0);
            }
            manti_size += 1;
        }
        while (manti_size < 52) {
            result_vector.push_back(0);
            manti_size++;
        }
        break;
    }
    default:
        break;
    }
    Float result(result_vector);
    return result;
}
Float Float::operator/=(const Float& other) {
    Float tmp = *this / other;
    *this = tmp;
    return *this;
}
bool Float::operator<(const Float& other) const {
    Float ltmp = *this;
    Float rtmp = other;
    double left = ltmp.translate_float_to_double();
    double right = rtmp.translate_float_to_double();
    return left < right;
}
bool Float::operator<=(const Float& other) const {
    Float ltmp = *this;
    Float rtmp = other;
    if (ltmp == rtmp) {
        return 1;
    }
    double left = ltmp.translate_float_to_double();
    double right = rtmp.translate_float_to_double();
    return left <= right;
}
bool Float::operator>(const Float& other) const {
    Float ltmp = *this;
    Float rtmp = other;
    double left = ltmp.translate_float_to_double();
    double right = rtmp.translate_float_to_double();
    return left > right;
}
bool Float::operator>=(const Float& other) const {
    Float ltmp = *this;
    Float rtmp = other;
    if (ltmp == rtmp) {
        return 1;
    }
    double left = ltmp.translate_float_to_double();
    double right = rtmp.translate_float_to_double();
    return left >= right;
}
bool Float::operator==(const Float& other) const {
    if (this->m_size != other.m_size) {
        return 0;
    }
    Index size = this->m_size;
    for (Index i = 0; i < size; i++) {
        if (this->at(i) != other.at(i)) {
            return 0;
        }
    }
    return 1;
}
bool Float::operator!=(const Float& other) const {
    if (this->m_size != other.m_size) {
        return 1;
    }
    Index size = this->m_size;
    for (Index i = 0; i < size; i++) {
        if (this->at(i) != other.at(i)) {
            return 1;
        }
    }
    return 0;
}

Float& Float::operator++() {
    *this += one;
    return *this;
}
Float Float::operator++(int num) {
    Float tmp(*this);
    *this += one;
    return tmp;
}
Float& Float::operator--() {
    *this -= one;
    return *this;
}
Float Float::operator--(int num) {
    Float tmp(*this);
    *this -= one;
    return tmp;
}

Float Float::abs() {
    std::vector<Float::Scalar> result_vector;
    Float tmp = *this;
    result_vector.push_back(0);
    for (Float::Index ii = 1; ii < tmp.getsize(); ii++) {
        result_vector.push_back(tmp.at(ii));
    }
    Float result(result_vector);
    return result;
}
Float Float::sqrt(){//the smaller the epsilon,the more accurate the result
    if (*this < zero) {
        throw std::invalid_argument("value cannot be a negative number.");
    }
    Float low = zero;
    Float up = (*this < one ? one : *this);
    Float mid = (low + up) / two;
    Float last(32);
    do
    {
        if (mid * mid > *this) {
            up = mid;
        }
        else {
            low = mid;
        }
        last = mid;
        mid = (up + low) / two;
    } while ((mid - last).abs() > epsilon2);
    return mid;
}
Float Float::pow(int power) {
    Float tmp = *this;
    if (power > 0) {
        for (int i = 0; i < power - 1; i++) {
            tmp *= *this;
        }
    }
    else if (power < 0) {
        for (int i = 0; i < -power - 1; i++) {
            tmp *= *this;
        }
        tmp = one / tmp;
    }
    else if (power == 0) {
        tmp = one;
    }
    return tmp;
}
//Float Float::operator<<(const long num) const;
//Float Float::operator>>(const long num) const;

int main() {
    /*
    Float::Index size = 32;
    Float my_float(size);
    std::vector<Float::Scalar> new_float = { 1,0,1,1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
    Float myf(new_float);
    Float::Scalar num = myf.at(0);
    std::cout << num << std::endl;
    std::cout << myf.getsize() << std::endl;
    std::cout << my_float.at(1) << std::endl;
    std::cout << my_float.getsize() << std::endl;
    for (int i = 0; i < 8; i++) {
        std::cout << myf.Float_Get_Expo()[i];
    }
    std::cout << " " << std::endl;
    for (int i = 0; i < 23; i++) {
        std::cout << myf.Float_Get_Manti()[i];
    }
    std::cout << " " << std::endl;
    std::cout << myf.translate_float_to_double() << std::endl;
    */

    /*------------------------------translate_double_to_float test----------------------------*/
    /*
    float integ = 0;
    float dec = 0.64;
    Float m_f(translate_double_to_float64(integ,dec));
    std::cout << m_f.getsize() << std::endl;
    std::cout << m_f.Float_Get_Sign() << std::endl;
    for (int i = 0; i < 11; i++) {
        std::cout << m_f.Float_Get_Expo()[i];
    }
    std::cout << " " << std::endl;
    for (int i = 0; i < 52; i++) {
        std::cout << m_f.Float_Get_Manti()[i];
    }
    std::cout << " " << std::endl;
    for (int i = 0; i < 64; i++) {
        std::cout << m_f.at(i);
    }
    */
    //std::cout << translate_vector_to_long({1,0,1,1}) << std::endl;
    
    /*---------------------------+-test----------------------------*/
    /*std::vector<Float::Scalar> new_float1 = {0, 0,1,1,1,1,0,0,0, 1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    std::vector<Float::Scalar> new_float2 = { 1, 0,1,1,1,1,0,0,1, 1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
    //std::vector<Float::Scalar> new_float2 = { 0, 0,1,1,1,1,1,1,1, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
    Float myf1(new_float1);
    puts("-----1");
    Float myf2(new_float2);
    //Float m_f = myf1 + myf2;
    //myf1 += myf2;
    Float m_f = myf1 - myf2;puts("-----2");
    std::cout << m_f.getsize() << std::endl;
    for (int i = 0; i < m_f.getsize(); i++) {
        std::cout << m_f.at(i) ;
    }
    std::cout << ' ' << std::endl;
    for (int i = 0; i < 8; i++) {
        std::cout << m_f.Float_Get_Expo()[i];
    }
    std::cout << " " << std::endl;
    for (int i = 0; i < 23; i++) {
        std::cout << m_f.Float_Get_Manti()[i];
    }
    std::cout << " " << std::endl;
    double a = m_f.translate_float_to_double();
    printf("%.10f\n", a);
    std::cout << a << std::endl;
    std::cout << (myf1 < myf2) << std::endl;
    std::cout << (myf1 <= myf2) << std::endl;
    std::cout << (myf1 > myf2) << std::endl;
    std::cout << (myf1 >= myf2) << std::endl;
    std::cout << (myf1 == myf2) << std::endl;
    std::cout << (myf1 != myf2) << std::endl;*/
    //std::vector<Float::Scalar> a = { 1,1,1,1,0,0,1 };
    //std::cout << translate_vector_to_double(a) << std::endl;

    /*------------------translate 0 test,output:5.87747e-39----------------------------*/
    /*std::vector<Float::Scalar> new_float3 = {0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    Float myf3(new_float3);
    std::cout << myf3.translate_float_to_double() << std::endl;
    */

    /*------------------------ * / test ---------------------------------------*/
    /*std::vector<Float::Scalar> new_float1 = {0, 0,1,1,1,1,0,0,0, 1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    std::vector<Float::Scalar> new_float2 = { 0, 0,1,1,1,1,0,1,1, 1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
    std::vector<Float::Scalar> new_float3 = { 0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
    Float myf1(new_float1);
    Float myf2(new_float2);
    Float m_f = myf1 / myf2; 
    //myf1 *= myf2;
    std::cout << m_f.getsize() << std::endl;
    for (int i = 0; i < m_f.getsize(); i++) {
        std::cout << m_f.at(i);
    }
    std::cout << ' ' << std::endl;
    for (int i = 0; i < 8; i++) {
        std::cout << m_f.Float_Get_Expo()[i];
    }
    std::cout << " " << std::endl;
    for (int i = 0; i < 23; i++) {
        std::cout << m_f.Float_Get_Manti()[i];
    }
    std::cout << " " << std::endl;
    double a = m_f.translate_float_to_double();
    printf("%.20f\n", a);
    std::cout << a << std::endl;
    */

    /*Float b(translate_double_to_float32(1, 0));
    for (int i = 0; i < b.getsize(); i++) {
        std::cout << b.at(i);
    }*/

    /*--------------------  ++ -- test -------------------------------*/
    /*std::vector<Float::Scalar> new_float1 = {0, 0,1,1,1,1,0,0,0, 1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    Float myf1(new_float1);
    myf1--;
    std::cout << myf1.getsize() << std::endl;
    for (int i = 0; i < myf1.getsize(); i++) {
        std::cout << myf1.at(i);
    }
    std::cout << ' ' << std::endl;
    for (int i = 0; i < 8; i++) {
        std::cout << myf1.Float_Get_Expo()[i];
    }
    std::cout << " " << std::endl;
    for (int i = 0; i < 23; i++) {
        std::cout << myf1.Float_Get_Manti()[i];
    }
    std::cout << " " << std::endl;
    double a = myf1.translate_float_to_double();
    printf("%.20f\n", a);
    std::cout << a << std::endl;*/
    
    /*
    std::vector<Float::Scalar> x = translate_double_to_float32(0, 0.0000002);
    for (int i = 0; i < 32; i++) {
        std::cout << x[i];
    }
    std::cout << " " << std::endl;
    Float X(x);
    std::cout << X.translate_float_to_double() << std::endl;*/

    /*---------------------abs and sqrt test---------------------------*/
    /*std::vector<Float::Scalar> new_float1 = {1, 0,1,1,1,1,0,0,0, 1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    std::vector<Float::Scalar> new_float2 = { 0, 0,1,1,1,1,0,1,1, 1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
    std::vector<Float::Scalar> new_float3 = { 0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
    Float myf1(new_float1);
    Float myf2(new_float2);
    Float a = myf1.abs();
    for (int i = 0; i < a.getsize(); i++) {
        std::cout << a.at(i);
    }
    std::cout << " " << std::endl;
    std::cout << myf2.translate_float_to_double() << std::endl;
    printf("%.20f\n", myf2.translate_float_to_double());
    Float b = myf2.sqrt();
    for (int i = 0; i < b.getsize(); i++) {
        std::cout << b.at(i);
    }
    std::cout << " " << std::endl;
    for (int i = 0; i < 8; i++) {
        std::cout << b.Float_Get_Expo()[i];
    }
    std::cout << " " << std::endl;
    for (int i = 0; i < 23; i++) {
        std::cout << b.Float_Get_Manti()[i];
    }
    std::cout << " " << std::endl;
    std::cout << b.translate_float_to_double() << std::endl;
    printf("%.20f\n", b.translate_float_to_double());*/
    
    /*-----------------------pow test-----------------------------------*/
    std::vector<Float::Scalar> new_float1 = {0, 0,1,1,1,1,0,0,0, 1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    std::vector<Float::Scalar> new_float2 = { 0, 0,1,1,1,1,0,1,1, 1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
    Float myf1(new_float1);
    Float myf2(new_float2);
    Float m_f = myf1.pow(-3);
    std::cout << m_f.getsize() << std::endl;
    for (int i = 0; i < m_f.getsize(); i++) {
        std::cout << m_f.at(i);
    }
    std::cout << ' ' << std::endl;
    for (int i = 0; i < 8; i++) {
        std::cout << m_f.Float_Get_Expo()[i];
    }
    std::cout << " " << std::endl;
    for (int i = 0; i < 23; i++) {
        std::cout << m_f.Float_Get_Manti()[i];
    }
    std::cout << " " << std::endl;
    double a = m_f.translate_float_to_double();
    printf("%.20f\n", a);
    std::cout << a << std::endl;


    return 0;
}