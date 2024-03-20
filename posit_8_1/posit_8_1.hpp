#ifndef POSIT_8_1_H__
#define POSIT_8_1_H__
#include <cstddef>
#include <iostream>
#include <vector>
class Posit_8_1 {
public:
    typedef bool Scalar;
    typedef std::size_t Index;
    Posit_8_1();
    Posit_8_1(const Posit_8_1& other);
    Posit_8_1(const int other, const int sig_bit, const int reg_bit, const int expo_bit, const int frac_bit);
    Posit_8_1& operator=(const Posit_8_1& other);
    ~Posit_8_1();
    Scalar at(Index i);
    Index getsize();
    Posit_8_1(const int other):m_size(8){m_posit = other;}
    Scalar Posit_8_1_Get_Sign();

    float translate_posit_8_1_to_float(); // translate from posit_8_1 to float
    

    Posit_8_1 operator+(const Posit_8_1& other) const;
    Posit_8_1 operator+=(const Posit_8_1& other) ;
    Posit_8_1 operator-(const Posit_8_1& other) const;
    Posit_8_1 operator-=(const Posit_8_1& other) ;
    Posit_8_1 operator*(const Posit_8_1& other) const;
    Posit_8_1 operator*=(const Posit_8_1& other) ;
    Posit_8_1 operator/(const Posit_8_1& other) const;
    Posit_8_1 operator/=(const Posit_8_1& other) ;


    bool operator<(const Posit_8_1& other) const;
    bool operator<=(const Posit_8_1& other) const;
    bool operator>(const Posit_8_1& other) const;
    bool operator>=(const Posit_8_1& other) const;
    bool operator==(const Posit_8_1& other) const;
    bool operator!=(const Posit_8_1& other) const;

    Posit_8_1& operator++();//prefix increment ++()
    Posit_8_1 operator++(int);//postfix increment ()++
    Posit_8_1& operator--();//prefix decrement --()
    Posit_8_1 operator--(int);//postfix decrement ()--

    int m_posit = 0;
    bool sign_bit = 0;
    int regime_bit = 0;
    int exponent_bit = 0;
    int fraction_bit = 0;
private:
    Index m_size;

};
const int useed = 4;
const Posit_8_1 zero(0x0, 0, 7, 0, 0);
const Posit_8_1 one (0x40, 0, 2, 1, 4);
const Posit_8_1 inf (0x80, 1, 7, 0, 0);
Posit_8_1::Scalar is_zero(Posit_8_1 num);
Posit_8_1::Scalar is_inf(Posit_8_1 num);
Posit_8_1 translate_float_to_posit_8_1(const float num); // translate from float to posit_8_1
#endif