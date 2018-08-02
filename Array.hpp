//
// Created by Pavel on 28.07.2018.
//
//#include <initializer_list>
#include <algorithm>
#include <iostream>
#include <type_traits>
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
U add(const U &a, const O &b){ return a + b; }
template <class U, class O>
U sub(const U &a, const O &b){ return a - b; }
template <class U, class O>
U mul(const U &a, const O &b){ return a * b; }
template <class U, class O>
U div(const U &a, const O &b){ return a / b; }


template <size_t size, class T1>
class Array
{
private:
    T1 * A;
    template<size_t n, class U, class O>
    friend Array<n, U> operation(const Array<n, U> &, const Array<n, O> &, U (*func)(const U &, const O &));
    template<size_t n, class U, typename O>
    friend Array<n, U> scalar_operation(const Array<n, U> &, const O &, U (*func)(const U &, const O &));
    template<size_t n, class U, class O>
    friend Array<n, U> & i_operation(Array<n, U> &, const Array<n, O> &, U & (*func)(U &, const O &));
    template<size_t n, class U, typename O>
    friend Array<n, U> & i_scalar_operation(Array<n, U> &, const O &, U & (*func)(U &, const O &));
public:
    constexpr static unsigned int size_ = size;
    Array();
    Array(const Array<size, T1> &);
    Array(Array<size, T1> &&);
    Array(T1 *, size_t massize);
    Array<size, T1> & operator=(const Array<size, T1> &);
    Array<size, T1> & operator=(Array<size, T1> &&);
    template <typename O>
    Array<size, T1> & operator=(const O &);
    template<class O>
    Array<size, T1> operator+(const Array<size, O> &b) const { return operation(*this, b, add); }
    template<typename O>
    Array<size, T1> operator+(const O &b) const { return scalar_operation(*this, b, add); }
    template<class O>
    Array<size, T1> operator-(const Array<size, O> &b) const { return operation(*this, b, sub); }
    template<typename O>
    Array<size, T1> operator-(const O &b) const { return scalar_operation(*this, b, sub); }
    template<class O>
    Array<size, T1> operator*(const Array<size, O> &b) const { return operation(*this, b, mul); }
    template<typename O>
    Array<size, T1> operator*(const O &b) const { return scalar_operation(*this, b, mul); }
    template<class O>
    Array<size, T1> operator/(const Array<size, O> &b) const { return operation(*this, b, div); }
    template<typename O>
    Array<size, T1> operator/(const O &b) const { return scalar_operation(*this, b, div); }
    template <class O>
    Array<size, T1> & operator+=(const Array<size, O> &b) { return i_operation(*this, b, iadd); }
    template<typename O>
    Array<size, T1> & operator+=(const O &b) { return i_scalar_operation(*this, b, iadd); }
    template <class O>
    Array<size, T1> & operator-=(const Array<size, O> &b) { return i_operation(*this, b, isub ); }
    template<typename O>
    Array<size, T1> & operator-=(const O &b) { return i_scalar_operation(*this, b, isub); }
    template <class O>
    Array<size, T1> & operator*=(const Array<size, O> &b) { return i_operation(*this, b, imul ); }
    template<typename O>
    Array<size, T1> & operator*=(const O &b) { return i_scalar_operation(*this, b, imul); }
    template <class O>
    Array<size, T1> & operator/=(const Array<size, O> &b) { return i_operation(*this, b, idiv ); }
    template<typename O>
    Array<size, T1> & operator/=(const O &b) { return i_scalar_operation(*this, b, idiv); }
    template <typename O>
    Array<size, O> apply_func(O (*func)(const T1 &));
    T1 & operator[](size_t i) { return A[i];}
    template<size_t n, class U>
    friend std::ostream & operator<<(std::ostream & out, const Array<n, U> &b);
    template<class U, size_t n, size_t k, size_t m>
    friend Array<n, Array<m, U>> dot(const Array<n, Array<k, U>> &, const Array<k, Array<m, U>> &);
    // very later TODO: Strassen algorithm
    Array<size, Array<1, T1>> T();
    ~Array();

    class iterator
    {
    private:
        T1 * pos;
        explicit iterator(T1 * ptr) {pos = ptr;}
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


template<size_t m, size_t n, class U>
class Array<m, Array<n, U>>
{
    Array<n, U> * A;
    template<size_t k, class G, class O>
    friend Array<k, G> operation(const Array<k, G> &, const Array<k, O> &, G (*func)(const G &, const O &));
    template<size_t k, class G, typename O>
    friend Array<k, G> scalar_operation(const Array<k, G> &, const O &, G (*func)(const G &, const O &));
    template<size_t k, class G, class O>
    friend Array<k, G> & i_operation(Array<k, G> &, const Array<k, O> &, G & (*func)(G &, const O &));
    template<size_t k, class G, typename O>
    friend Array<k, G> & i_scalar_operation(Array<k, G> &, const O &, G & (*func)(G &, const O &));

public:
    constexpr static unsigned int size_ = m;
    Array();
    Array(const Array<m, Array<n, U>> &);
    Array(Array<m, Array<n, U>> &&);
    Array<m, Array<n, U>> & operator=(const Array<m, Array<n, U>> &);
    Array<m, Array<n, U>> & operator=(Array<m, Array<n, U>> &&);
    template <typename O>
    Array<m, Array<n, U>> & operator=(const O &);
    template <class O>
    Array<m, Array<n, U>> operator+(const Array<m, Array<n, O>> &b) const { return operation(*this, b, add); }
    template <typename O>
    Array<m, Array<n, U>> operator+(const O &b) const { return scalar_operation(*this, b, add); }
    template <class O>
    Array<m, Array<n, U>> operator-(const Array<m, Array<n, O>> &b) const { return operation(*this, b, sub); }
    template <typename O>
    Array<m, Array<n, U>> operator-(const O &b) const { return scalar_operation(*this, b, sub); }
    template <class O>
    Array<m, Array<n, U>> operator*(const Array<m, Array<n, O>> &b) const { return operation(*this, b, mul); }
    template <typename O>
    Array<m, Array<n, U>> operator*(const O &b) const { return scalar_operation(*this, b, mul); }
    template <class O>
    Array<m, Array<n, U>> operator/(const Array<m, Array<n, O>> &b) const { return operation(*this, b, div); }
    template <typename O>
    Array<m, Array<n, U>> operator/(const O &b) const { return scalar_operation(*this, b, div); }
    template <class O>
    Array<m, Array<n, U>> & operator+=(const Array<m, Array<n, O>> &b) { return i_operation(*this, b, iadd); }
    template<typename O>
    Array<m, Array<n, U>> & operator+=(const O &b) { return i_scalar_operation(*this, b, iadd); }
    template <class O>
    Array<m, Array<n, U>> & operator-=(const Array<m, Array<n, O>> &b) { return i_operation(*this, b, isub); }
    template<typename O>
    Array<m, Array<n, U>> & operator-=(const O &b) { return i_scalar_operation(*this, b, isub); }
    template <class O>
    Array<m, Array<n, U>> & operator*=(const Array<m, Array<n, O>> &b) { return i_operation(*this, b, imul); }
    template<typename O>
    Array<m, Array<n, U>> & operator*=(const O &b) { return i_scalar_operation(*this, b, imul); }
    template <class O>
    Array<m, Array<n, U>> & operator/=(const Array<m, Array<n, O>> &b) { return i_operation(*this, b, idiv); }
    template<typename O>
    Array<m, Array<n, U>> & operator/=(const O &b) { return i_scalar_operation(*this, b, idiv); }
    Array<n, U> & operator[](size_t i) { return A[i];}
    template <size_t fm, size_t fn, class fU>
    friend std::ostream & operator<<(std::ostream & out, const Array<fm, Array<fn, fU>> &b);
    Array<n, Array<m, U>> T();
    template <typename O>
    Array<m, Array<n, O>> apply_func(O (*func)(const U &));
    template <size_t k>
    Array<m, Array<k, U>> dot(const Array<n, Array<k, U>> &);

};

template<size_t m, size_t n, class U> template <size_t k>
Array<m, Array<k, U>> Array<m, Array<n, U>>::dot(const Array<n, Array<k, U>> &b)
{
    Array<n, Array<m, U>> res;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j) {
            res.A[i][j] = 0;
            for (size_t l = 0; l < k; ++l)
                res.A[i][j] += A[i][l] * b.A[l][j];
        }
    return res;
}

template<size_t fm, size_t fn, class fU>
std::ostream & operator<<(std::ostream & out, const Array<fm, Array<fn, fU>> &b)
{
    out << '[';
    for (size_t i = 0; i < fn - 1; ++i)
        out << b.A[i] << std::endl << ' ';
    return out << b.A[fn - 1] << ']' << std::endl;
}


template<size_t size, class T1>
Array<size, Array<1, T1>> Array<size, T1>::T()
{
    Array<size, Array<1, T1>> res;
    for (size_t i = 0; i < size; ++i)
        res.A[i].A[0] = A[i];
    return res;
}

template<size_t m, size_t n, class U>
Array<n, Array<m, U>> Array<m, Array<n, U>>::T()
{
    Array<n, Array<m, U>> res;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            res.A[i][j] = A[j][i];
    return res;
}

template <size_t size, class T1>
Array<size, T1>::Array()
{
    A = new T1[size];
}

template<size_t m, size_t n, class U>
Array<m, Array<n, U>>::Array()
{
    A = new Array<n, U>[m];
}

template<size_t m, size_t n, class U>
Array<m, Array<n, U>>::Array(const Array<m, Array<n, U>> &b)
{
    A = new Array<n, U>[m];
    for (size_t i = 0; i < m; ++i)
        A[i] = b.A[i];
}

template <size_t size, class T1>
Array<size, T1>::Array(const Array<size, T1> &b)
{
    A = new T1[size];
    for (size_t i = 0; i < size; ++i)
        A[i] = b.A[i];
}

template<size_t m, size_t n, class U>
Array<m, Array<n, U>>::Array(Array<m, Array<n, U>> &&b)
{
    A = b.A;
    b.A = NULL;
}

template <size_t size, class T1>
Array<size, T1>::Array(Array<size, T1> &&b)
{
    A = b.A;
    b.A = NULL;
}

template <size_t size, class T1>
Array<size, T1>::Array(T1 * b, size_t massize)
{
    A = new T1[size];
    for (size_t i = 0; i < massize; ++i)
        A[i] = b[i];
}

template<size_t m, size_t n, class U>
Array<m, Array<n, U>> & Array<m, Array<n, U>>::operator=(const Array<m, Array<n, U>> &b)
{
    delete[] A;
    A = new Array<n, U>[m];
    for (size_t i = 0; i < m; ++i)
        A[i] = b.A[i];
    return *this;
}

template <size_t size, class T1>
Array<size, T1> & Array<size, T1>::operator=(const Array<size, T1> &b)
{
    delete[] A;
    A = new T1[size];
    for (size_t i = 0; i < size; ++i)
        A[i] = b.A[i];
    return *this;
}



template<size_t m, size_t n, class U>
Array<m, Array<n, U>> & Array<m, Array<n, U>>::operator=(Array<m, Array<n, U>> &&b)
{
    delete[] A;
    A = b.A;
    b.A = NULL;
    return *this;
}

template <size_t size, class T1>
Array<size, T1> & Array<size, T1>::operator=(Array<size, T1> &&b)
{
    delete[] A;
    A = b.A;
    b.A = NULL;
    return *this;
}

template <size_t size, class T1> template <typename O>
Array<size, T1> & Array<size, T1>::operator=(const O &b)
{
    delete[] A;
    A = new T1[size];
    for (size_t i = 0; i < size; ++i)
        A[i] = b;
    return *this;
}

template<size_t m, size_t n, class U> template <typename O>
Array<m, Array<n, U>> & Array<m, Array<n, U>>::operator=(const O &b)
{
    delete[] A;
    A = new Array<n, U>[m];
    for (size_t i = 0; i < m; ++i)
        A[i] = b;
    return *this;
}

template <size_t size, class T1> template <typename O>
Array<size, O> Array<size, T1>::apply_func(O (*func)(const T1 &))
{
    Array<size, O> res = *this;
    std::for_each(res.A, res.A + size, func);
    return res;
}

template<size_t m, size_t n, class U> template <typename O>
Array<m, Array<n, O>> Array<m, Array<n, U>>::apply_func(O (*func)(const U &))
{
    Array<m, Array<n, O>> res = *this;
    for (size_t i = 0; i < m; ++i)
        res.A[i] = res.A[i].apply_func(func);
    return res;
}

template<size_t n, class U, class O>
Array<n, U> & i_operation(Array<n, U> &a, const Array<n, O> &b, U & (*func)(U &, const O &))
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

template<size_t n, class U, class O>
Array<n, U> operation(const Array<n, U> &a, const Array<n, O> &b, U (*func)(const U &, const O &))
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

template <size_t n, class U>
std::ostream & operator<<(std::ostream & out, const Array<n, U> &b)
{
    out << "[";
    for (size_t i = 0; i < n - 1; ++i)
        out << b.A[i] << " ";
    return out << b.A[n - 1] << "]";
}

template <size_t size, class T1>
Array<size, T1>::~Array()
{
    delete[] A;
}

#endif //C_ARRAY_H
