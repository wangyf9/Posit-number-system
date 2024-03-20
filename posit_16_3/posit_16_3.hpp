#ifndef POSIT_16_3_H__
#define POSIT_16_3_H__
#include <cstddef>
#include <iostream>
#include <vector>
class Posit_16_3 {
public:
    typedef bool Scalar;
    typedef std::size_t Index;
    Posit_16_3();
    Posit_16_3(const Posit_16_3& other);
    Posit_16_3(const int other, const int sig_bit, const int reg_bit, const int expo_bit, const int frac_bit);
    Posit_16_3& operator=(const Posit_16_3& other);
    ~Posit_16_3();
    Scalar at(Index i);
    Index getsize();
    Posit_16_3(const int other):m_size(16){m_posit = other;}
    Scalar Posit_16_3_Get_Sign();

    float translate_posit_16_3_to_float(); // translate from posit_16_3 to float
    

    Posit_16_3 operator+(const Posit_16_3& other) const;
    Posit_16_3 operator+=(const Posit_16_3& other) ;
    Posit_16_3 operator-(const Posit_16_3& other) const;
    Posit_16_3 operator-=(const Posit_16_3& other) ;
    Posit_16_3 operator*(const Posit_16_3& other) const;
    Posit_16_3 operator*=(const Posit_16_3& other) ;
    Posit_16_3 operator/(const Posit_16_3& other) const;
    Posit_16_3 operator/=(const Posit_16_3& other) ;


    bool operator<(const Posit_16_3& other) const;
    bool operator<=(const Posit_16_3& other) const;
    bool operator>(const Posit_16_3& other) const;
    bool operator>=(const Posit_16_3& other) const;
    bool operator==(const Posit_16_3& other) const;
    bool operator!=(const Posit_16_3& other) const;

    Posit_16_3& operator++();//prefix increment ++()
    Posit_16_3 operator++(int);//postfix increment ()++
    Posit_16_3& operator--();//prefix decrement --()
    Posit_16_3 operator--(int);//postfix decrement ()--

    int m_posit = 0;
    bool sign_bit = 0;
    int regime_bit = 0;
    int exponent_bit = 0;
    int fraction_bit = 0;
private:
    Index m_size;

};
const int useed = 256;
const Posit_16_3 zero(0x0, 0, 15, 0, 0);
const Posit_16_3 one (0x4000, 0, 2, 3, 10);
const Posit_16_3 inf (0x8000, 1, 15, 0, 0);
Posit_16_3::Scalar is_zero(Posit_16_3 num);
Posit_16_3::Scalar is_inf(Posit_16_3 num);
Posit_16_3 translate_float_to_posit_16_3(const float num); // translate from float to posit_16_3
#endif