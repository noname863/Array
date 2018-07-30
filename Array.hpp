//
// Created by Pavel on 28.07.2018.
//
//#include <initializer_list>
#include <algorithm>
#include <iostream>
#ifndef C_ARRAY_H
#define C_ARRAY_H

template <class U, class O>
U & iadd(U &a, const O &b){ return a += b; }
template <class U, class O>
U & isub(U &a, const O &b){ return a -= b; }
template <class U, class O>
U & imul(U &a, const O &b){ return a *= b; }
template <class U, class O>
U & idiv(U &a, const O &b){ return a /= b; }
template <class U, class O>
U add(U &a, const O &b){ return a + b; }
template <class U, class O>
U sub(U &a, const O &b){ return a - b; }
template <class U, class O>
U mul(U &a, const O &b){ return a * b; }
template <class U, class O>
U div(U &a, const O &b){ return a / b; }


template <size_t size, class T>
class Array
{
private:
    T * A;
    template<size_t n, class U>
    friend Array<n, U> operation(const Array<n, U> &, const Array<n, U> &, U (*func)(const U &, const U &));
    template<size_t n, class U, typename O>
    friend Array<n, U> scalar_operation(const Array<n, U> &, const O &, U (*func)(const U &, const O &));
    template<size_t n, class U>
    friend Array<n, U> & i_operation(Array<n, U> &, const Array<n, U> &, U & (*func)(U &, const U &));
    template<size_t n, class U, typename O>
    friend Array<n, U> & i_scalar_operation(Array<n, U> &, const O &, U & (*func)(U &, const O &));
public:
    Array();
    Array(const Array<size, T> &);
    Array(Array<size, T> &&);
    Array(T *, size_t massize);
    Array<size, T> & operator=(const Array<size, T> &);
    Array<size, T> & operator=(Array<size, T> &&);
    template<size_t n, class U>
    friend Array<n, U> operator+(const Array<n, U> &a, const Array<n, U> &b);
    template<size_t n, class U, typename O>
    friend Array<n, U> operator+(const Array<n, U> &a,const O &b);
    template<size_t n, class U>
    friend Array<n, U> operator-(const Array<n, U> &a,const Array<n, U> &b);
    template<size_t n, class U, typename O>
    friend Array<n, U> operator-(const Array<n, U> &a,const O &b);
    template<size_t n, class U>
    friend Array<n, U> operator*(const Array<n, U> &a,const Array<n, U> &b);
    template<size_t n, class U, typename O>
    friend Array<n, U> operator*(const Array<n, U> &a,const O &b);
    template<size_t n, class U>
    friend Array<n, U> operator/(const Array<n, U> &a,const Array<n, U> &b);
    template<size_t n, class U, typename O>
    friend Array<n, U> operator/(const Array<n, U> &a,const O &b);
    Array<size, T> & operator+=(const Array<size, T> &b) { return i_operation(*this, b, iadd); }
    template<typename O>
    Array<size, T> & operator+=(const O &b) { return i_scalar_operation(*this, b, iadd); }
    Array<size, T> & operator-=(const Array<size, T> &b) { return i_operation(*this, b, isub ); }
    template<typename O>
    Array<size, T> & operator-=(const O &b) { return i_scalar_operation(*this, b, isub); }
    Array<size, T> & operator*=(const Array<size, T> &b) { return i_operation(*this, b, imul ); }
    template<typename O>
    Array<size, T> & operator*=(const O &b) { return i_scalar_operation(*this, b, imul); }
    Array<size, T> & operator/=(const Array<size, T> &b) { return i_operation(*this, b, idiv ); }
    template<typename O>
    Array<size, T> & operator/=(const O &b) { return i_scalar_operation(*this, b, idiv); }
    template<typename O>
    Array<size, T> apply_func(void func(O &));
    T & operator[](size_t i) { return A[i];}
    template<size_t n, class U>
    friend std::ostream & operator<<(std::ostream & out, const Array<n, U> &b);
    template<class U, size_t n, size_t k, size_t m>
    friend Array<n, Array<m, U>> dot(const Array<n, Array<k, U>> &, const Array<k, Array<m, U>> &);
    // very later TODO: Strassen algorithm
    template<class U, size_t n, size_t m>
    friend Array<n, Array<m, U>> T(const Array<m, Array<n, U>> &);
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

template<size_t n, class U>
Array<n, U> operator+(const Array<n, U> &a, const Array<n, U> &b)
{
    return operation(a, b, add);
}

template<size_t n, class U, typename O>
Array<n, U> operator+(const Array<n, U> &a, const O &b)
{
    return scalar_operation(a, b, add);
}

template<size_t n, class U>
Array<n, U> operator-(const Array<n, U> &a, const Array<n, U> &b)
{
    return operation(a, b, sub);
}

template<size_t n, class U, typename O>
Array<n, U> operator-(const Array<n, U> &a, const O &b)
{
    return scalar_operation(a, b, sub);
}

template<size_t n, class U>
Array<n, U> operator*(const Array<n, U> &a, const Array<n, U> &b)
{
    return operation(a, b, mul);
}

template<size_t n, class U, typename O>
Array<n, U> operator*(const Array<n, U> &a, const O &b)
{
    return scalar_operation(a, b, mul);
}

template<size_t n, class U>
Array<n, U> operator/(const Array<n, U> &a, const Array<n, U> &b)
{
    return operation(a, b, div);
}

template<size_t n, class U, typename O>
Array<n, U> operator/(const Array<n, U> &a, const O &b)
{
    return scalar_operation(a, b, div);
}

template<class U, size_t n, size_t m>
Array<n, Array<m, U>> T(const Array<m, Array<n, U>> &b)
{
    Array<n, Array<m, U>> res;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            res.A[i].A[j] = b.A[j].A[i];
    return res;
}

template <size_t size, class T>
Array<size, T>::Array()
{
    A = new T[size];
}

template <size_t size, class T>
Array<size, T>::Array(const Array<size, T> &b)
{
    A = new T[size];
    for (size_t i = 0; i < size; ++i)
        A[i] = b.A[i];
}

template <size_t size, class T>
Array<size, T>::Array(Array<size, T> &&b)
{
    A = b.A;
    b.A = NULL;
}

template <size_t size, class T>
Array<size, T>::Array(T * b, size_t massize)
{
    A = new T[size];
    for (size_t i = 0; i < massize; ++i)
        A[i] = b.A[i];
}

template <size_t size, class T>
Array<size, T> & Array<size, T>::operator=(const Array<size, T> &b)
{
    delete[] A;
    A = new T[size];
    for (size_t i = 0; i < size; ++i)
        A[i] = b.A[i];
    return *this;
}

template <size_t size, class T>
Array<size, T> & Array<size, T>::operator=(Array<size, T> &&b)
{
    delete[] A;
    A = b.A;
    b.A = NULL;
    //b.A = new T[size]
    return *this;
}

template<size_t n, class U>
Array<n, U> & i_operation(Array<n, U> &a, const Array<n, U> &b, U & (*func)(U &, const U &))
{
    for (size_t i = 0; i < n; ++i)
        func(a.A[i], b.A[i]);
    return a;
}

template<size_t n, class U, typename O>
Array<n, U> & i_scalar_operation(Array<n, U> &a, const O &b, U & (*func)(U &, const O &))
{
    for (size_t i = 0; i < n; ++i)
        func(a.A[i], b);
    return a;
}

template <size_t n, class U>
std::ostream & operator<<(std::ostream & out, const Array<n, U> &b)
{
    out << "[";
    for (size_t i = 0; i < n - 1; ++i)
        out << b.A[i] << " ";
    return out << b.A[n - 1] << "]" << std::endl;
}

template <size_t size, class T>
Array<size, T>::~Array()
{
    delete[] A;
}

template<size_t n, class U>
Array<n, U> operation(const Array<n, U> &a, const Array<n, U> &b, U (*func)(const U &, const U &))
{
    Array<n, U> res;
    for (size_t i = 0; i < n; ++i)
        res.A[i] = func(a.A[i], b.A[i]);
    return res;
}

template<size_t n, class U, typename O>
Array<n, U> scalar_operation(const Array<n, U> &a, const O &b, U (*func)(const U &, const O &))
{
    Array<n, U> res;
    for (size_t i = 0; i < n; ++i)
        res.A[i] = func(a.A[i], b);
    return res;
}

template<class U, size_t n, size_t k, size_t m>
Array<n, Array<m, U>> dot(const Array<n, Array<k, U>> &a, const Array<k, Array<m, U>> &b)
{
    Array<n, Array<m, U>> res;
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
