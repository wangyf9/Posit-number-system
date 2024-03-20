#ifndef POSIT_10_1_H__
#define POSIT_10_1_H__
#include <cstddef>
#include <iostream>
#include <vector>
class Posit_10_1 {
public:
    typedef bool Scalar;
    typedef std::size_t Index;
    Posit_10_1();
    Posit_10_1(const Posit_10_1& other);
    Posit_10_1(const int other, const int sig_bit, const int reg_bit, const int expo_bit, const int frac_bit);
    Posit_10_1& operator=(const Posit_10_1& other);
    ~Posit_10_1();
    Scalar at(Index i);
    Index getsize();
    Posit_10_1(const int other):m_size(10){m_posit = other;}
    Scalar Posit_10_1_Get_Sign();

    float translate_posit_10_1_to_float(); // translate from posit_10_1 to float
    

    Posit_10_1 operator+(const Posit_10_1& other) const;
    Posit_10_1 operator+=(const Posit_10_1& other) ;
    Posit_10_1 operator-(const Posit_10_1& other) const;
    Posit_10_1 operator-=(const Posit_10_1& other) ;
    Posit_10_1 operator*(const Posit_10_1& other) const;
    Posit_10_1 operator*=(const Posit_10_1& other) ;
    Posit_10_1 operator/(const Posit_10_1& other) const;
    Posit_10_1 operator/=(const Posit_10_1& other) ;


    bool operator<(const Posit_10_1& other) const;
    bool operator<=(const Posit_10_1& other) const;
    bool operator>(const Posit_10_1& other) const;
    bool operator>=(const Posit_10_1& other) const;
    bool operator==(const Posit_10_1& other) const;
    bool operator!=(const Posit_10_1& other) const;

    Posit_10_1& operator++();//prefix increment ++()
    Posit_10_1 operator++(int);//postfix increment ()++
    Posit_10_1& operator--();//prefix decrement --()
    Posit_10_1 operator--(int);//postfix decrement ()--

    int m_posit = 0;
    bool sign_bit = 0;
    int regime_bit = 0;
    int exponent_bit = 0;
    int fraction_bit = 0;
private:
    Index m_size;

};
const int useed = 4;
const Posit_10_1 zero(0x0, 0, 9, 0, 0);
const Posit_10_1 one (0x100, 0, 2, 1, 6);
const Posit_10_1 inf (0x400, 1, 9, 0, 0);
Posit_10_1::Scalar is_zero(Posit_10_1 num);
Posit_10_1::Scalar is_inf(Posit_10_1 num);
Posit_10_1 translate_float_to_posit_10_1(const float num); // translate from float to posit_10_1
#endif