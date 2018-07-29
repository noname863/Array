//
// Created by Pavel on 28.07.2018.
//
//#include <initializer_list>
#include <algorithm>
#include <iostream>
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
    template<class U, size_t n>
    friend Array<U, n> operator+(const Array<U, n> &,const  Array<U, n> &);
    template<class U, size_t n>
    friend Array<U, n> operator-(const Array<U, n> &,const  Array<U, n> &);
    template<class U, size_t n>
    friend Array<U, n> operator*(const Array<U, n> &,const  Array<U, n> &);
    template<class U, size_t n>
    friend Array<U, n> operator/(const Array<U, n> &,const  Array<U, n> &);
    template<class U, size_t n>
    friend Array<U, n> operator%(const Array<U, n> &,const  Array<U, n> &);
    Array<T, size> & operator+=(const Array<T, size> &);
    Array<T, size> & operator-=(const Array<T, size> &);
    Array<T, size> & operator*=(const Array<T, size> &);
    Array<T, size> & operator/=(const Array<T, size> &);
    Array<T, size> & operator%=(const Array<T, size> &);
    T & operator[](size_t i) { return A[i];}
    template<class U, size_t n>
    friend std::ostream & operator<<(std::ostream & out, const Array<U, n> &b);
    template<class U, size_t n, size_t k, size_t m>
    friend Array<Array<U, m>, n> dot(const Array<Array<U, k>, n> &, const Array<Array<U, m>, k> &);
    template<class U, size_t n, size_t m>
    friend Array<Array<U, m>, n> T(const Array<Array<U, n>, m> &);
    ~Array();

    class iterator
    {
    private:
        T * pos;
        iterator(T * ptr) {pos = ptr;}
    public:
        iterator(const iterator &iter) {this->pos = iter.pos;}
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
/*
template<class U, size_t m, size_t n>
class Array<Array<U, n>, m>
{
    //Array<Array<U, n>, m> T();
    //template <size_t k>
    //Array<Array<U, k>, m> dot(const Array<Array<U, k>, n> &b);
};
*/

template<class U, size_t n, size_t m>
Array<Array<U, m>, n> T(const Array<Array<U, n>, m> &b)
{
    Array<Array<U, m>, n> res;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            res.A[i].A[j] = b.A[j].A[i];
    return res;
}
/*
template<class U, size_t n, size_t m, size_t k>
Array<Array<U, k>, m> Array<Array<U, n>, m>::dot(const Array<Array<U, k>, n> &b)
{
    Array<Array<U, k>, m> res;
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < k; ++j)
        {
            res.A[i].A[j] = 0;
            for (size_t l = 0; l < n; ++l)
                res.A[i].A[j] += (a.A[i].A[l] + b.A[l].A[j]);
        }
    return res;
}
*/
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
    b.A = new T[size];
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
    b.A = new T[size];
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

template <class U, size_t fsize>
std::ostream & operator<<(std::ostream & out, const Array<U, fsize> &b)
{
    out << "[";
    for (size_t i = 0; i < fsize - 1; ++i)
        out << b.A[i] << " ";
    return out << b.A[fsize - 1] << "]" << std::endl;
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

template<class U, size_t n, size_t k, size_t m>
Array<Array<U, m>, n> dot(const Array<Array<U, k>, n> &a, const Array<Array<U, m>, k> &b)
{
    Array<Array<U, m>, n> res;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
        {
            res.A[i].A[j] = 0;
            for (size_t l = 0; l < k; ++l)
                res.A[i].A[j] += a.A[i].A[l] * b.A[l].A[j];
        }
    return res;
}


#endif //C_ARRAY_H
