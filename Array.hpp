//
// Created by Pavel on 28.07.2018.
//
//#include <initializer_list>
#include <algorithm>
#include <iostream>
#include <type_traits>
#include <tuple>
#include <utility>
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
    constexpr static auto shape = std::make_tuple(size);
    Array();
    Array(const Array<size, T1> &);
    Array(Array<size, T1> &&);
    explicit Array(T1 (*generator)());
    Array(T1 *, size_t massize);
    template <typename O>
    Array(const O &);
    Array<size, T1> & operator=(const Array<size, T1> &);
    Array<size, T1> & operator=(Array<size, T1> &&);
    template <typename O>
    Array<size, T1> & operator=(const O &);
    void generate(T1 (*generator)());
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
    const T1 & operator[](size_t i) const { return A[i];}
    template<size_t n, class U>
    friend std::ostream & operator<<(std::ostream & out, const Array<n, U> &b);
    void write(std::ostream &) const;
    void read(std::istream &);
    Array<size, Array<1, T1>> T();
    ~Array();
    class iterator
        {
        protected:
            T1 * pos;
        public:
            explicit iterator(T1 * ptr) {pos = ptr;}
            iterator(const iterator &iter) {this->pos = iter.pos;}
            iterator & operator++() {++pos; return *this;}
            iterator operator++(int) {iterator copy(*this); ++pos; return copy;}
            iterator & operator--() {--pos; return *this;}
            iterator operator--(int) {iterator copy(*this); --pos; return copy;}
            T1 & operator*() const {return *pos;}
            iterator & operator+=(int i) {pos += i; return *this;}
            iterator & operator-=(int i) {pos -= i; return *this;}
            friend iterator operator+(const iterator &a, int b) {return iterator(a.pos + b);}
            friend iterator operator-(const iterator &a, int b) {return iterator(a.pos - b);}
            friend bool operator==(const iterator &a, const iterator &b) {return a.pos == b.pos;}
            friend bool operator!=(const iterator &a, const iterator &b) {return a.pos != b.pos;}
            friend int operator-(const iterator &a, const iterator &b) {return (int)a.pos - (int)b.pos;}
    };

    class citerator : public iterator
    {
    public:
        explicit citerator(T1 * ptr) : iterator(ptr) {}
        T1 operator*() {return *iterator::pos;}
    };
    iterator begin() {return iterator(A);}
    iterator end() {return iterator(A + size);}
    citerator begin() const { return citerator(A); }
    citerator end() const {return citerator(A + size);}

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
    struct skip_init{};
    Array(skip_init) {}
public:
    constexpr static auto shape = std::tuple_cat(std::make_tuple(m), Array<n, U>::shape);
    Array();
    Array(const Array<m, Array<n, U>> &);
    Array(Array<m, Array<n, U>> &&);
    explicit Array(U (*generator)());
    template <typename O>
    Array(const O &);
    Array<m, Array<n, U>> & operator=(const Array<m, Array<n, U>> &);
    Array<m, Array<n, U>> & operator=(Array<m, Array<n, U>> &&);
    template <typename O>
    Array<m, Array<n, U>> & operator=(const O &);
    void generate(U (*generator)());
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
    template <typename O>
    Array<m, Array<n, U>> & operator+=(const O &b) { return i_scalar_operation(*this, b, iadd); }
    template <class O>
    Array<m, Array<n, U>> & operator-=(const Array<m, Array<n, O>> &b) { return i_operation(*this, b, isub); }
    template <typename O>
    Array<m, Array<n, U>> & operator-=(const O &b) { return i_scalar_operation(*this, b, isub); }
    template <class O>
    Array<m, Array<n, U>> & operator*=(const Array<m, Array<n, O>> &b) { return i_operation(*this, b, imul); }
    template <typename O>
    Array<m, Array<n, U>> & operator*=(const O &b) { return i_scalar_operation(*this, b, imul); }
    template <class O>
    Array<m, Array<n, U>> & operator/=(const Array<m, Array<n, O>> &b) { return i_operation(*this, b, idiv); }
    template <typename O>
    Array<m, Array<n, U>> & operator/=(const O &b) { return i_scalar_operation(*this, b, idiv); }
    Array<n, U> & operator[](size_t i) { return A[i];}
    const Array<n, U> & operator[](size_t i) const { return A[i];}
    template <size_t fm, size_t fn, class fU>
    friend std::ostream & operator<<(std::ostream & out, const Array<fm, Array<fn, fU>> &b);
    Array<n, Array<m, U>> T();
    template <typename O>
    Array<m, Array<n, O>> apply_func(O (*func)(const U &));
    template <size_t k>
    Array<m, Array<k, U>> easy_dot(const Array<n, Array<k, U>> &) const;
    template <class O, size_t l, size_t k, size_t p>
    friend Array<l, Array<p, O>> dot(const Array<l, Array<k, O>> &a, const Array<k, Array<p, O>> &b);
    template <size_t k>
    Array<m, Array<k, U>> dot(const Array<n, Array<k, U>> &b) const;
    template <size_t l, size_t k>
    Array<l, Array<k, U>> submatrix(size_t x1, size_t x2) const;
    void write(std::ostream &) const;
    void read(std::istream &);
    class iterator
    {
    protected:
        Array<n, U> * pos;
    public:
        explicit iterator(Array<n, U> * ptr) {pos = ptr;}
        iterator(const iterator &iter) {this->pos = iter.pos;}
        iterator & operator++() {++pos; return *this;}
        iterator operator++(int) {iterator copy(*this); ++pos; return copy;}
        iterator & operator--() {--pos; return *this;}
        iterator operator--(int) {iterator copy(*this); --pos; return copy;}
        Array<n, U> & operator*() const {return *pos;}
        iterator & operator+=(int i) {pos += i; return *this;}
        iterator & operator-=(int i) {pos -= i; return *this;}
        friend iterator operator+(const iterator &a, int b) {return iterator(a.pos + b);}
        friend iterator operator-(const iterator &a, int b) {return iterator(a.pos - b);}
        friend bool operator==(const iterator &a, const iterator &b) {return a.pos == b.pos;}
        friend bool operator!=(const iterator &a, const iterator &b) {return a.pos != b.pos;}
        friend int operator-(const iterator &a, const iterator &b) {return (int)a.pos - (int)b.pos;}
    };
    class citerator : public iterator
    {
    public:
        explicit citerator(Array<n, U> * ptr) : iterator(ptr) {}
        Array<n, U> operator*() {return *iterator::pos;}
    };
    iterator begin() {return iterator(A);}
    iterator end() {return iterator(A + m);}
    citerator begin() const { return citerator(A); }
    citerator end() const {return citerator(A + m);}
};

template<size_t m, size_t n, class U> template <size_t k>
Array<m, Array<k, U>> Array<m, Array<n, U>>::easy_dot(const Array<n, Array<k, U>> &b) const
{
    Array<n, Array<m, U>> res;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
        {
            res[i][j] = 0;
            for (size_t l = 0; l < k; ++l)
                res[i][j] += this->A[i][l] * b[l][j];
        }
    return res;
}

template<size_t m, size_t n, class U> template <size_t k>
Array<m, Array<k, U>> Array<m, Array<n, U>>::dot(const Array<n, Array<k, U>> &b) const
{
    constexpr bool is_small = ((n * m * k) < (64 * 64 * 64));
    constexpr size_t smaller_n = (n >> 1) + (n & 1);
    constexpr size_t smaller_k = (k >> 1) + (k & 1);
    constexpr size_t smaller_m = (m >> 1) + (m & 1);
    if (is_small)
        return this->easy_dot(b);
    auto a11 = this->template submatrix<smaller_m, smaller_n>(0, 0);
    auto a12 = this->template submatrix<smaller_m, smaller_n>(0, smaller_n);
    auto a21 = this->template submatrix<smaller_m, smaller_n>(smaller_m, 0);
    auto a22 = this->template submatrix<smaller_m, smaller_n>(smaller_m, smaller_n);
    auto b11 = b.template submatrix<smaller_n, smaller_k>(0, 0);
    auto b12 = b.template submatrix<smaller_n, smaller_k>(0, smaller_k);
    auto b21 = b.template submatrix<smaller_n, smaller_k>(smaller_n, 0);
    auto b22 = b.template submatrix<smaller_n, smaller_k>(smaller_n, smaller_k);
    auto p1 = (a11 + a22).dot(b11 + b22);
    auto p2 = (a21 + a22).dot(b11);
    auto p3 = a11.dot(b12 - b22);
    auto p4 = a22.dot(b21 - b11);
    auto p5 = (a11 + a12).dot(b22);
    auto p6 = (a21 - a11).dot(b11 + b12);
    auto p7 = (a12 - a22).dot(b21 + b22);
    Array<m, Array<k, U>> res;
    size_t j;
    size_t i = 0;
    for (; i < smaller_m; ++i)
    {
        for (j = 0; j < smaller_k; ++j)
        {
            res[i][j] = p1[i][j] + p4[i][j] - p5[i][j] + p7[i][j];
            res[i][j + smaller_k] = p3[i][j] + p5[i][j];
            res[i + smaller_m][j] = p2[i][j] + p4[i][j];
            res[i + smaller_m][j + smaller_k] = p1[i][j] - p2[i][j] + p3[i][j] + p6[i][j];
        }
    }
    return res;
}

template<class O, size_t l, size_t k, size_t p>
Array<l, Array<p, O>> dot(const Array<l, Array<k, O>> &a, const Array<k, Array<p, O>> &b)
{
    return a.dot(b);
}

template<size_t fm, size_t fn, class fU>
std::ostream & operator<<(std::ostream & out, const Array<fm, Array<fn, fU>> &b)
{
    out << '[';
    for (size_t i = 0; i < fm - 1; ++i)
        out << b.A[i] << std::endl << ' ';
    return out << b.A[fm - 1] << ']' << std::endl;
}


template<size_t size, class T1>
Array<size, Array<1, T1>> Array<size, T1>::T()
{
    Array<size, Array<1, T1>> res;
    for (size_t i = 0; i < size; ++i)
        res[i][0] = A[i];
    return res;
}

template<size_t m, size_t n, class U>
Array<n, Array<m, U>> Array<m, Array<n, U>>::T()
{
    Array<n, Array<m, U>> res;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            res[i][j] = A[j][i];
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

template <size_t size, class T1>
Array<size, T1>::Array(T1 (*generator)())
{
    A = new T1[size];
    for (size_t i = 0; i < size; ++i)
        A[i] = generator();
}


template<size_t m, size_t n, class U>
Array<m, Array<n, U>>::Array(U (*generator)())
{
    A = new Array<n, U>[m];
    for (size_t i = 0; i < m; ++i)
        A[i].generate(generator);
}

template <size_t size, class T1> template <typename O>
Array<size, T1>::Array(const O &b)
{
    A = new T1[size];
    *this = b;
}

template<size_t m, size_t n, class U> template <typename O>
Array<m, Array<n, U>>::Array(const O &b)
{
    A = new Array<n, U>[m];
    *this = b;
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

template <size_t size, class T1>
void Array<size, T1>::generate(T1 (*generator)())
{
    for (size_t i = 0; i < size; ++i)
        A[i] = generator();
}

template<size_t m, size_t n, class U>
void Array<m, Array<n, U>>::generate(U (*generator)())
{
    for (size_t i = 0; i < m; ++i)
        A[i].generate(generator);
}

template <size_t size, class T1> template <typename O>
Array<size, O> Array<size, T1>::apply_func(O (*func)(const T1 &))
{
    Array<size, O> res;
    for (size_t i = 0; i < size; ++i)
        res[i] = func(A[i]);
    return res;
}

template<size_t m, size_t n, class U> template <typename O>
Array<m, Array<n, O>> Array<m, Array<n, U>>::apply_func(O (*func)(const U &))
{
    Array<m, Array<n, O>> res;
    for (size_t i = 0; i < m; ++i)
        res[i] = A[i].apply_func(func);
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
template <size_t n, class T>
void Array<n, T>::write(std::ostream &out) const
{
    out.write((char*)A, n * sizeof(A[0]));
}

template <size_t n, class T>
void Array<n, T>::read(std::istream &in)
{
    in.read((char *)A, n * sizeof(A[0]));
}

template<size_t m, size_t n, class U>
void Array<m, Array<n, U>>::write(std::ostream &out) const
{
    std::for_each(begin(), end(), [&out](const Array<n, U> &a){ a.write(out); });
}

template<size_t m, size_t n, class U>
void Array<m, Array<n, U>>::read(std::istream &in)
{
    std::for_each(begin(), end(), [&in](Array<n, U> &a){ a.read(in); });
}


template<size_t n, class U, typename O>
Array<n, U> scalar_operation(const Array<n, U> &a, const O &b, U (*func)(const U &, const O &))
{
    Array<n, U> res;
    for (size_t i = 0; i < n; ++i)
        res.A[i] = func(a.A[i], b);
    return res;
}

template<size_t m, size_t n, class U> template<size_t l, size_t k>
Array<l, Array<k, U>> Array<m, Array<n, U>>::submatrix(size_t x1, size_t x2) const
{
    Array<l, Array<k, U>> res;
    size_t i = 0;
    if ((l + x1 <= m) && (k + x2 <= n))
        for (; i < l; ++i)
            res[i] = Array<k, U>((*(U**)(this->A + x1 + i) + x2), k);
    else
    {
        for (; i + x1 < m; ++i)
            res[i] = Array<k, U>((*(U**)(this->A + x1 + i) + x2), n - x2);
        for (; i < l; ++i)
            res[i] = 0;
        for (i = 0; i < l; ++i)
            for (size_t j = n - x2; j < k; ++j)
                res[i][j] = 0;
    }
    return res;
}

template<class U, size_t n, size_t k, size_t m>
Array<n, Array<m, U>> easy_dot(const Array<n, Array<k, U>> &a, const Array<k, Array<m, U>> &b)
{
    Array<n, Array<m, U>> res;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
        {
            res[i][j] = 0;
            for (size_t l = 0; l < k; ++l)
                res[i][j] += a[i][l] * b[l][j];
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
