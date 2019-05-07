//
// Created by Leandro Ishi Soares de Lima on 23/02/2019.
//

#include "Utils.h"


/**
 * Check if the MSABase accepts the candidateBase (see https://www.bioinformatics.org/sms/iupac.html)
 * @param MSABase - a base of the MSA
 * @param candidateBase - a candidate base (ACGT-)
 * @return
 */
bool Utils::accepts(char MSABase, char candidateBase) {
    static const std::map<char, std::unordered_set<char>> baseToAcceptedBases = {
            {'A', {'A'}},
            {'C', {'C'}},
            {'G', {'G'}},
            {'T', {'T'}},
            {'-', {'-'}},
            {'R', {'A', 'G'}},
            {'Y', {'C', 'T'}},
            {'S', {'C', 'G'}},
            {'W', {'A', 'T'}},
            {'K', {'G', 'T'}},
            {'M', {'A', 'C'}},
            {'B', {'C', 'G', 'T'}},
            {'D', {'A', 'G', 'T'}},
            {'H', {'A', 'C', 'T'}},
            {'V', {'A', 'C', 'G'}},
            {'N', {'A', 'C', 'G', 'T'}}
    };
    return baseToAcceptedBases.at(MSABase).count(candidateBase) > 0;
}
