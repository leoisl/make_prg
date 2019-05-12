//
// Created by Leandro Ishi Soares de Lima on 23/02/2019.
//

#include "Utils.h"


const std::map<char, std::unordered_set<char>> Utils::baseToAcceptedBases = {
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


void Utils::expandSequencesRecursively(uint32_t i, const std::string &seq, std::list<std::string> &allExpandedSequences) {
    if (i < seq.size()) { //base case
        char originalBase = seq[i]; //get the evaluated base
        if (isNotACGT(originalBase)) {
            //need to expand everyone
            std::list<std::string> newExpandedSequences;
            for (const string &seqToBeExpanded : allExpandedSequences) {
                for (char translatedBase : baseToAcceptedBases.at(originalBase))
                    newExpandedSequences.push_back(seqToBeExpanded + translatedBase);
            }

            //move this new sequences to allExpandedSequences
            allExpandedSequences = std::move(newExpandedSequences);
        } else {
            //originalBase is ACGT
            //add to all seqs
            for (string &seqToBeExpanded : allExpandedSequences)
                seqToBeExpanded += originalBase;
        }

        //recursive call
        expandSequencesRecursively(i+1, seq, allExpandedSequences);
    }
}

std::list<std::string> Utils::expandnonACGT(const std::string &seq) {
    std::list<std::string> allExpandedSequences;
    expandSequencesRecursively(0, seq, allExpandedSequences);
    return allExpandedSequences;
}
