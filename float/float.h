#ifndef FLOAT_H__
#define FLOAT_H__
#include <cstddef>
#include <iostream>
#include <vector>
class Float {
public:
    typedef bool Scalar;
    typedef std::size_t Index;
    //Float();
    Float(Index size);
    Float(const std::vector<Scalar> other);
    Float(const Float& other);
    Float& operator=(const Float& other);
    ~Float();
    Scalar at(Index i);
    const Scalar at(Index i) const;
    Index getsize();
    std::vector<Scalar> block(Index i, Index j) const; // Block starting from i to j

    double translate_float_to_double(); // translate from float to double
    Float translate_double_to_float32(const long integer, const float num); // translate from double to float32
    Float translate_double_to_float64(const long integer, const float num); // translate from double to float64

    Scalar Float_Get_Sign();
    std::vector<Scalar> Float_Get_Expo() ;
    std::vector<Scalar> Float_Get_Manti() ;

    Float operator+(const Float& other) const;
    Float operator+=(const Float& other) ;
    Float operator-(const Float& other) const;
    Float operator-=(const Float& other) ;
    Float operator*(const Float& other) const;
    Float operator*=(const Float& other) ;
    Float operator/(const Float& other) const;
    Float operator/=(const Float& other) ;
    bool operator<(const Float& other) const;
    bool operator<=(const Float& other) const;
    bool operator>(const Float& other) const;
    bool operator>=(const Float& other) const;
    bool operator==(const Float& other) const;
    bool operator!=(const Float& other) const;

    Float& operator++();//prefix increment ++()
    Float operator++(int num);//postfix increment ()++
    Float& operator--();//prefix decrement --()
    Float operator--(int num);//postfix decrement ()--
    
    Float abs() ;
    Float sqrt() ;
    Float pow(int power);
    
    Float operator<<(const long num) const;
    Float operator>>(const long num) const;


    
private:
    std::vector<Scalar> m_float;
    Index m_size;
};

long translate_vector_to_long(std::vector<Float::Scalar> vec);//used to translate expo part vector
//std::vector<Float::Scalar> teanslate_bits_to_vector32(uint32_t num);
double translate_vector_to_double(std::vector<Float::Scalar> vec);//used to translate manti part vector
Float::Scalar is_zero(Float num);
#endif