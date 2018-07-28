//
// Created by Pavel on 28.07.2018.
//
//#include <initializer_list>
#include <algorithm>
#ifndef C_ARRAY_H
#define C_ARRAY_H

template <class T, size_t size>
class Array
{
private:
    T * A;
public:
    Array();
    Array(const Array<T, size> &);
    Array(Array<T, size> &&);
    Array(T *, size_t massize);
    Array<T, size> & operator=(const Array<T, size> &);
    Array<T, size> & operator=(Array<T, size> &&);
    template<class U, size_t fsize>
    friend Array<U, fsize> operator+(const Array<U, fsize> &,const  Array<U, fsize> &);
    template<class U, size_t fsize>
    friend Array<U, fsize> operator-(const Array<U, fsize> &,const  Array<U, fsize> &);
    template<class U, size_t fsize>
    friend Array<U, fsize> operator*(const Array<U, fsize> &,const  Array<U, fsize> &);
    template<class U, size_t fsize>
    friend Array<U, fsize> operator/(const Array<U, fsize> &,const  Array<U, fsize> &);
    template<class U, size_t fsize>
    friend Array<U, fsize> operator%(const Array<U, fsize> &,const  Array<U, fsize> &);
    Array<T, size> & operator+=(const Array<T, size> &);
    Array<T, size> & operator-=(const Array<T, size> &);
    Array<T, size> & operator*=(const Array<T, size> &);
    Array<T, size> & operator/=(const Array<T, size> &);
    Array<T, size> & operator%=(const Array<T, size> &);
    T & operator[](size_t i) { return A[i];}
    ~Array();

    friend class iterator
    {
    private:
        T * pos;
        iterator(T * ptr) {pos = ptr;}
    public:
        iterator(const iter &) {this->pos = iter.pos;}
        iterator & operator++() {++pos; return *this;}
        iterator operator++(int) {iterator copy(*this); ++pos; return copy;}
        iterator & operator--() {--pos; return *this;}
        iterator operator--(int) {iterator copy(*this); --pos; return copy;}
        iterator operator*() {return *pos;}
        iterator & operator+=(int i) {pos += i; return *this;}
        iterator & operator-=(int i) {pos -= i; return *this;}
        friend iterator operator+(iterator a, int b) {return iterator(a.pos + b);}
        friend iterator operator-(iterator a, int b) {return iterator(a.pos - b);}
        friend int operator-(iterator a, iterator b) {return (int)a.pos - (int)b.pos;}
    };
};

template <class T, size_t size>
Array<T, size>::Array()
{
    A = new T[size];
}

template <class T, size_t size>
Array<T, size>::Array(const Array<T, size> &b)
{
    A = new T[size];
    for (size_t i = 0; i < size; ++i)
        A[i] = b.A[i];
}

template <class T, size_t size>
Array<T, size>::Array(Array<T, size> &&b)
{
    A = b.A;
    //test it well
}

template <class T, size_t size>
Array<T, size>::Array(T * b, size_t massize)
{
    A = new T[size];
    for (size_t i = 0; i < massize; ++i)
        A[i] = b.A[i];
}

template <class T, size_t size>
Array<T, size> & Array<T, size>::operator=(const Array<T, size> &b)
{
    delete[] A;
    A = new T[size];
    for (size_t i = 0; i < size; ++i)
        A[i] = b.A[i];
    return *this;
}

template <class T, size_t size>
Array<T, size> & Array<T, size>::operator=(Array<T, size> &&b)
{
    delete[] A;
    A = b.A;
    return *this;
}

template <class T, size_t size>
Array<T, size> & Array<T, size>::operator+=(const Array<T, size> &b)
{
    for (size_t i = 0; i < size; ++i)
        A[i] += b.A[i];
    return *this;
}

template <class T, size_t size>
Array<T, size> & Array<T, size>::operator-=(const Array<T, size> &b)
{
    for (size_t i = 0; i < size; ++i)
        A[i] -= b.A[i];
    return *this;
}

template <class T, size_t size>
Array<T, size> & Array<T, size>::operator*=(const Array<T, size> &b)
{
    for (size_t i = 0; i < size; ++i)
        A[i] *= b.A[i];
    return *this;
}

template <class T, size_t size>
Array<T, size> & Array<T, size>::operator/=(const Array<T, size> &b)
{
    for (size_t i = 0; i < size; ++i)
        A[i] /= b.A[i];
    return *this;
}

template <class T, size_t size>
Array<T, size> & Array<T, size>::operator%=(const Array<T, size> &b)
{
    for (size_t i = 0; i < size; ++i)
        A[i] %= b.A[i];
    return *this;
}

template <class T, size_t size>
Array<T, size>::~Array()
{
    delete[] A;
}


template<class U, size_t fsize>
Array<U, fsize> operator+(const Array<U, fsize> &a,const  Array<U, fsize> &b)
{
    Array<U, fsize> res;
    for (size_t i = 0; i < fsize; ++i)
        res.A[i] = a.A[i] + b.A[i];
    return res;
}

template<class U, size_t fsize>
Array<U, fsize> operator-(const Array<U, fsize> &a,const  Array<U, fsize> &b)
{
    Array<U, fsize> res;
    for (size_t i = 0; i < fsize; ++i)
        res.A[i] = a.A[i] - b.A[i];
    return res;
}

template<class U, size_t fsize>
Array<U, fsize> operator*(const Array<U, fsize> &a,const  Array<U, fsize> &b)
{
    Array<U, fsize> res;
    for (size_t i = 0; i < fsize; ++i)
        res.A[i] = a.A[i] * b.A[i];
    return res;
}

template<class U, size_t fsize>
Array<U, fsize> operator/(const Array<U, fsize> &a,const  Array<U, fsize> &b)
{
    Array<U, fsize> res;
    for (size_t i = 0; i < fsize; ++i)
        res.A[i] = a.A[i] / b.A[i];
    return res;
}

template<class U, size_t fsize>
Array<U, fsize> operator%(const Array<U, fsize> &a,const  Array<U, fsize> &b)
{
    Array<U, fsize> res;
    for (size_t i = 0; i < fsize; ++i)
        res.A[i] = a.A[i] % b.A[i];
    return res;
}

#endif //C_ARRAY_H
