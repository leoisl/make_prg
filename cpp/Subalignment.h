//
// Created by Leandro Ishi Soares de Lima on 01/05/2019.
//

#ifndef MAKE_PRG_SUBALIGNMENT_H
#define MAKE_PRG_SUBALIGNMENT_H

#include "includes.h"
#include "Utils.h"

enum IntervalType {
    UNPROCESSED, //those that do not correspond to any of the classifications below (not still processed) - correspond to subalignments to be broken into subalignments of the types below
    MATCH, // intervals with length >= k, which we manage to get a consensus string -> These should be printed and the result is the consensus string
    NONMATCH, //intervals with length >= k, that we did not manage to get a consensus string
    NONMATCH_MAX_NESTING_LEVEL, //nonmatch that we can't divide further -> These should be printed by getting the unique representative sequences
    TOO_SHORT //non-match or match, it does not matter, these are very short intervals -> These should be printed by getting the unique representative sequences
};
std::ostream &operator<<(std::ostream &os, const IntervalType &intervalType);

class Interval {
    /**
     * Represents an interval and its type
     */
public:
    uint32_t start, end; // (start, end]
    IntervalType intervalType;

    //ctors
    Interval() = default;
    Interval(uint32_t start, uint32_t end, IntervalType intervalType=UNPROCESSED) :
            start{start}, end{end}, intervalType{intervalType} {}
    Interval(const Interval &interval) = default;

    inline uint32_t getLength() const { return end-start; }

    //streams
    friend std::ostream &operator<<(std::ostream &os, const Interval &interval) {
        os << "(" << interval.start << ", " << interval.end << "]: " << interval.intervalType;
        return os;
    }
};

class SubAlignment {
private:
    /*
     * Represents a sub-alignment - a vertical slice - a set of sequences and a horizontal interval
     */
    std::vector<uint32_t> sequencesNumbers; //the sequences in this sub-alignment
    Interval interval; //begin and end are always equal to all sub-alignments
    const std::vector<std::string> *MSA; //allow us to retrieve the sequences of the sub-alignment themselves


    //builds the consensus string of this subalignment
    std::string buildConsensusString() const;

public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //CONSTRUCTORS
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Main ctors:
    SubAlignment(uint32_t sequenceNumberLower, uint32_t sequenceNumberUpper, const Interval &interval,
                 const std::vector<std::string> *MSA) :
            sequencesNumbers{}, interval{interval}, MSA{MSA} {
        for (auto i = sequenceNumberLower; i < sequenceNumberUpper; ++i)
            sequencesNumbers.push_back(i);
    }
    SubAlignment(const std::vector<uint32_t> &sequencesNumbers, const Interval &interval,
                 const std::vector<std::string> *MSA) :
            sequencesNumbers{sequencesNumbers}, interval{interval}, MSA{MSA} {}

    //This is not nice, but we need a default ctor so that Boost can create nodes
    //TODO: hard to deal with this, check after if we have options
    SubAlignment() = default;

    //default copy ctor/=
    SubAlignment(const SubAlignment &) = default;
    SubAlignment &operator=(const SubAlignment &) = default;

    //move ctor/=
    SubAlignment(SubAlignment &&rValueRef) {
        if (this != &rValueRef) {
            this->sequencesNumbers = std::move(rValueRef.sequencesNumbers);
            this->interval = std::move(rValueRef.interval); //nothing changes, but it is a nice practice
            this->MSA = std::move(rValueRef.MSA); //nothing changes, but it is a nice practice
        }
    }
    SubAlignment &operator=(SubAlignment &&rValueRef) {
        //TODO: this a copy of the move ctor, refactor
        if (this != &rValueRef) {
            this->sequencesNumbers = std::move(rValueRef.sequencesNumbers);
            this->interval = std::move(rValueRef.interval); //nothing changes, but it is a nice practice
            this->MSA = std::move(rValueRef.MSA); //nothing changes, but it is a nice practice
        }
        return *this;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //CONSTRUCTORS
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //GETTERS
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inline const std::vector<uint32_t>& getSequencesNumbers() const { return sequencesNumbers; }
    inline const Interval& getInterval() const { return interval; }
    inline Interval& getInterval() { return interval; }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //GETTERS
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //MAIN METHODS
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * add a sequence to this subalignment
     * @param sequenceNumber
     */
    inline void addSequence(uint32_t sequenceNumber) {
        sequencesNumbers.push_back(sequenceNumber);
    }


    /**
     * Return the match and non-match intervals of this subalignment WRT the positions in the global MSA.
     * Consensus sequences longer than k are match intervals and the rest as non-match intervals.
     * @param k - the minimum length to consider a vertical stripe as a match interval
     * @return list of intervals
     */
    std::vector<Interval> getMatchAndNonMatchIntervals(uint32_t k) const;


    /**
     * Given this subalignment, return the sequences in this alignment AS THEY ARE
     * Do not process anything, just get the sequences and return
     * @return vector with the alignments (strings)
     */
    std::vector<std::string> getSequences() const;

    /**
    * 1/ Removes "-" from all alignments
    * 2/ Remove all duplicates
    */
    std::vector<std::string> getRepresentativeSequences() const;



    /**
     * Split this subalignment into several subaligments, where each is a cluster of similar sequences
     * @return a vector of subalignments
     */
    std::vector<SubAlignment> kMeansCluster(uint32_t k) const;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //MAIN METHODS
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //MISC
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend std::ostream &operator<<(std::ostream &os, const SubAlignment &subAlignment) {
        os << "sequencesNumbers: ";
        for (auto sequenceNumber : subAlignment.sequencesNumbers)
            os << sequenceNumber << " ";
        os << std::endl;
        os << "Interval: " << subAlignment.interval << std::endl;
        return os;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //MISC
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};

#endif //MAKE_PRG_SUBALIGNMENT_H
