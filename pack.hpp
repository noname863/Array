//
// Created by Pavel on 07.08.2018.
//
#include <cstring>
#ifndef PACK_HPP
#define PACK_HPP

template <class T, size_t num>
struct pack_get {};

template <class T, size_t num>
struct type_get {};

template <class T, T... values>
struct pack{};

template<typename... Types>
class list
{
public:
    static void init(char * s) {}
};

template<typename... Types>
struct type_pack
{
    constexpr static size_t len = 0;
    template <typename... second_Types>
    using type_cat = type_pack<second_Types...>;
    using to_list = list<Types...>;
};

template<typename T1, typename... Types>
struct type_pack<T1, Types...>
{
    constexpr static size_t len = type_pack<Types...>::len + 1;
    template <typename... second_Types>
    using type_cat = type_pack<T1, Types..., second_Types...>;
    template <typename... second_Types>
    using r_type_cat = type_pack<second_Types..., T1, Types...>;
    using to_list = list<T1, Types...>;
    template <size_t num>
    using get = typename type_get<type_pack<T1, Types...>, num>::value;
};

template <size_t size = 0, typename... Types>
struct sum_of_sizeof
{
    constexpr static size_t value = 0;
};

template <size_t size, typename T1, typename... Types>
struct sum_of_sizeof<size, T1, Types...>
{
    constexpr static size_t value = sizeof(T1) + sum_of_sizeof<size - 1, Types...>::value;
};

template <typename T1, typename... Types>
struct sum_of_sizeof<0, T1, Types...>
{
    constexpr static size_t value = 0;
};

template<typename T1, typename... Types>
class list<T1, Types...>
{
    char c[sum_of_sizeof<type_pack<T1, Types...>::len, T1, Types...>::value];
public:
    static void init(char * s) { new ((T1*)s) T1; list<Types...>::init(s + sizeof(T1));}
    constexpr static size_t len = type_pack<T1, Types...>::len;
    list() {list<T1, Types...>::init(c);}
    template <size_t num>
    using get_type = typename type_pack<T1, Types...>::template get<num>;
    template <size_t num>
    get_type<num> & get() { return *(get_type<num>*)(c + sum_of_sizeof<num, T1, Types...>::value); }
};

template <class T, T first, T... values>
struct pack<T, first, values...>
{
    constexpr static size_t len = pack<T, values...>::len + 1;
    template <size_t num>
    constexpr static T get() { return pack_get<pack<T, first, values...>, num>::value;}
    static T get(size_t i)
    {
        if (i == 0)
            return first;
        else
            return pack<T, values...>::get(i - 1);
    }
};

template <class T>
struct pack<T>
{
    constexpr static size_t len = 0;
    static T get(size_t i) {return T();}
};

template <typename T1, typename... Types, size_t num>
struct type_get<type_pack<T1, Types...>, num>
{
    using value = typename type_get<type_pack<Types...>, num - 1>::value;
};

template <typename T1, typename... Types>
struct type_get<type_pack<T1, Types...>, 0>
{
    using value = T1;
};

template <class T, T first, T... values, size_t num>
struct pack_get<pack<T, first, values...>, num>
{
    constexpr static T value = pack_get<pack<T, values...>, num - 1>::value;
};

template <class T, T first, T... values>
struct pack_get<pack<T, first, values...>, 0>
{
    constexpr static T value = first;
};

#endif //PACK_HPP
