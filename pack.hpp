//
// Created by Pavel on 07.08.2018.
//

#ifndef PACK_HPP
#define PACK_HPP

template <class T, size_t num>
struct pack_get {};

template <class T, T... values>
struct pack{};

template <class T, T first, T... values>
struct pack<T, first, values...>
{
    constexpr static size_t len = pack<T, values...>::len + 1;
    template <size_t num>
    constexpr static T get() { return pack_get<pack<T, first, values...>, num>::value;}
};

template <class T>
struct pack<T>
{
    constexpr static size_t len = 0;
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
