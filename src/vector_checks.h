#pragma once
#include "vector.h"

template<int N, scalar T>
template<int V_SIZE>
bool Vector<N, T>::isOrthogonal(const std::array<Vector<N, T>, V_SIZE>& vectors) {
    for (int i = 0; i < vectors.size(); i++) {
        for (int j = i; j < vectors.size(); j++) {
            if (i == j)
                continue;

            if (vectors[i].componentDot(vectors[j]) != 0) {
                return false;
            }
        }
    }

    return true;
}

template<int N, scalar T>
template<int V_SIZE>
bool Vector<N, T>::isOrthonormal(const std::array<Vector<N, T>, V_SIZE>& vectors) {
    for (int i = 0; i < vectors.size(); i++) {
        if (vectors[i].componentDot(vectors[i]) != 1) {
            return false;
        }

        for (int j = i; j < vectors.size(); j++) {
            if (i == j)
                continue;

            if (vectors[i].componentDot(vectors[j]) != 0) {
                return false;
            }
        }
    }

    return true;
}