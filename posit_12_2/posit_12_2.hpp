#ifndef POSIT_12_2_H__
#define POSIT_12_2_H__
#include <cstddef>
#include <iostream>
#include <vector>
class Posit_12_2 {
public:
    typedef bool Scalar;
    typedef std::size_t Index;
    Posit_12_2();
    Posit_12_2(const Posit_12_2& other);
    Posit_12_2(const int other, const int sig_bit, const int reg_bit, const int expo_bit, const int frac_bit);
    Posit_12_2& operator=(const Posit_12_2& other);
    ~Posit_12_2();
    Scalar at(Index i);
    Index getsize();
    Posit_12_2(const int other):m_size(12){m_posit = other;}
    Scalar Posit_12_2_Get_Sign();

    float translate_posit_12_2_to_float(); // translate from posit_16_3 to float
    

    Posit_12_2 operator+(const Posit_12_2& other) const;
    Posit_12_2 operator+=(const Posit_12_2& other) ;
    Posit_12_2 operator-(const Posit_12_2& other) const;
    Posit_12_2 operator-=(const Posit_12_2& other) ;
    Posit_12_2 operator*(const Posit_12_2& other) const;
    Posit_12_2 operator*=(const Posit_12_2& other) ;
    Posit_12_2 operator/(const Posit_12_2& other) const;
    Posit_12_2 operator/=(const Posit_12_2& other) ;


    bool operator<(const Posit_12_2& other) const;
    bool operator<=(const Posit_12_2& other) const;
    bool operator>(const Posit_12_2& other) const;
    bool operator>=(const Posit_12_2& other) const;
    bool operator==(const Posit_12_2& other) const;
    bool operator!=(const Posit_12_2& other) const;

    Posit_12_2& operator++();//prefix increment ++()
    Posit_12_2 operator++(int);//postfix increment ()++
    Posit_12_2& operator--();//prefix decrement --()
    Posit_12_2 operator--(int);//postfix decrement ()--

    int m_posit = 0;
    bool sign_bit = 0;
    int regime_bit = 0;
    int exponent_bit = 0;
    int fraction_bit = 0;
private:
    Index m_size;

};
const int useed = 16;
const Posit_12_2 zero(0x0, 0, 11, 0, 0);
const Posit_12_2 one (0x400, 0, 2, 2, 7);
const Posit_12_2 inf (0x800, 1, 11, 0, 0);
Posit_12_2::Scalar is_zero(Posit_12_2 num);
Posit_12_2::Scalar is_inf(Posit_12_2 num);
Posit_12_2 translate_float_to_posit_12_2(const float num); // translate from float to posit_16_3
#endif