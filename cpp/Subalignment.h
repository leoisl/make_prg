//
// Created by Leandro Ishi Soares de Lima on 01/05/2019.
//

#ifndef MAKE_PRG_SUBALIGNMENT_H
#define MAKE_PRG_SUBALIGNMENT_H

#include "includes.h"

/* if begin and end equal, no need for this
class AlignmentSubstring {
public: //TODO: do proper encapsulation here?

    /*
     * Represents an alignment substring (sequence # the alignment come from, begin and end position of the alignment)
     * You need the MSA to retrieve the alignment itself
     sequenceNumber, begin, end;
    AlignmentSubstring() = default;
    AlignmentSubstring(uint32_t sequenceNumber, uint32_t begin, uint32_t end) :
        sequenceNumber{sequenceNumber}, begin{begin}, end{end}
    {}
};
 */

enum IntervalType {MATCH, NONMATCH};

class Interval {

public:
    uint32_t start, end; // (start, end]
    IntervalType intervalType;

    Interval(uint32_t start, uint32_t end, IntervalType intervalType) :
        start{start}, end{end}, intervalType{intervalType}
    {}


};

class SubAlignment {
private:
    /*
     * Represents a sub-alignment - a collection of AlignmentSubstring
     */
    std::vector<uint32_t> sequencesNumbers;
    uint32_t begin, end; //begin and end are always equal to all sub-alignments - TODO: refactor into interval
    const std::vector<std::string> *MSA;


    //builds the consensus string of this subalignment
    //TODO: update this !
    //TODO: replicate from python code
    std::string buildConsensusString() const;

public:
    //Main ctors:
    SubAlignment (uint32_t sequenceNumberLower, uint32_t sequenceNumberUpper, uint32_t begin, uint32_t end, const std::vector<std::string> *MSA) :
        sequencesNumbers{}, begin{begin}, end{end}, MSA{MSA} {
      for (auto i = sequenceNumberLower; i < sequenceNumberUpper; ++i)
        sequencesNumbers.push_back(i);
    }

    SubAlignment (const std::vector<uint32_t> &sequencesNumbers, uint32_t begin, uint32_t end, const std::vector<std::string> *MSA) :
        sequencesNumbers{sequencesNumbers}, begin{begin}, end{end}, MSA{MSA}
    {}

    //This is not nice, but we need a default ctor so that Boost can create nodes
    SubAlignment() = default;

    //default copy ctor/=
    SubAlignment(const SubAlignment &) = default;
    SubAlignment& operator=(const SubAlignment &) = default;

    //move ctor/=
    SubAlignment(SubAlignment && rValueRef) {
      if (this != &rValueRef) {
        this->sequencesNumbers = std::move(rValueRef.sequencesNumbers);
        this->begin = rValueRef.begin;
        this->end = rValueRef.end;
        this->MSA = rValueRef.MSA;
      }
    }

    SubAlignment& operator=(SubAlignment && rValueRef) {
      //TODO: this a copy of the move ctor, refactor
      if (this != &rValueRef) {
        this->sequencesNumbers = std::move(rValueRef.sequencesNumbers);
        this->begin = rValueRef.begin;
        this->end = rValueRef.end;
        this->MSA = rValueRef.MSA;
      }
      return *this;
    }

    inline void addSequence(uint32_t sequenceNumber) {
      sequencesNumbers.push_back(sequenceNumber);
    }

    //get the match and non-match intervals WRT the positions in the global MSA
    std::vector<Interval> getIntervals(uint32_t k) const;

    inline uint32_t getBegin() const { return begin; }
    inline uint32_t getEnd() const { return end; }
    inline const std::vector<uint32_t>& getSequencesNumbers() const { return sequencesNumbers; }

    friend std::ostream& operator<<(std::ostream& os, const SubAlignment& subAlignment) {
      os << "sequencesNumbers: ";
      for (auto sequenceNumber : subAlignment.sequencesNumbers)
        os << sequenceNumber << " ";
      os << std::endl;
      os << "begin: " << subAlignment.begin << std::endl;
      os << "end: " << subAlignment.end << std::endl;
      return os;
    }

    /**
     * Given this subalignment, return the sequences in this alignment AS THEY ARE
     * Do not process anything, just get the sequences and return
     * @return vector with the alignments (strings)
     */
    std::vector<std::string> getSequences() const;

    /**
     * Given this subalignment, return the sequences representing this alignment
     * Deals with N, non-ACGT bases, remove duplicates, etc...
     * @return vector of the representative strings
     */
    std::vector<std::string> getRepresentativeSequences() const;
};

#endif //MAKE_PRG_SUBALIGNMENT_H
