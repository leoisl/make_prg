//
// Created by Leandro Ishi Soares de Lima on 01/05/2019.
//

#include "Subalignment.h"


std::ostream &operator<<(std::ostream &os, const IntervalType &intervalType) {
    switch (intervalType) {
        case UNCLASSIFIED:
            os << "UNCLASSIFIED";
            break;
        case MATCH:
            os << "MATCH";
            break;
        case NONMATCH:
            os << "NONMATCH";
            break;
        case NONMATCH_MAX_NESTING_LEVEL:
            os << "NONMATCH_MAX_NESTING_LEVEL";
            break;
        case TOO_SHORT:
            os << "TOO_SHORT";
            break;
    }
    return os;
}


/**
 * Creates a consensus string from the aligment represented by this.
 * IUPAC bases results in consensus on that base if possible (see e.g. https://www.bioinformatics.org/sms/iupac.html)
 * - is taken into account and we can have a column full of -.
 * If there is a column where the consensus could be more than one base, we choose it at random
 *
 * @return a consensus of the alignment
 */
std::string SubAlignment::buildConsensusString() const {
    //TODO: issue warning on many non ACGT- bases?
    static const std::string consensusBases{"ACGT-"}; //which bases should we have in the final consensus?

    //generate the consensus string
    std::string consensusString;
    for (size_t j = interval.start; j < interval.end; ++j) {
        //1. Checks which base can be accepted as consensus in this column
        std::string acceptedBases;
        for (char candidateBase : consensusBases) {
            if (std::all_of(sequencesNumbers.begin(), sequencesNumbers.end(),
                    [&j, &candidateBase, this](uint32_t sequenceNumber) {
                        char MSABase = (*(this->MSA))[sequenceNumber][j];
                        try {
                            return Utils::accepts(MSABase, candidateBase);
                        }catch (const std::out_of_range &exception) {
                            BOOST_LOG_TRIVIAL(fatal) << "Unknown base in MSA: " << MSABase;
                            std::exit(1);
                        }
                    }
            )) {
                //everyone accepted the candidate base, add it
                acceptedBases+=candidateBase;
            }
        }

        //2. chooses a random accepted base if there was a consensus
        char acceptedBase = '*'; //assumes no consensus
        if (acceptedBases.size() > 0) {
            //we had a consensus, choose random base from the accepted ones
            acceptedBase = acceptedBases[std::rand() % acceptedBases.size()];
        }

        //3. add the accpted base to the consensus string
        consensusString += acceptedBase;
    }

    return consensusString;
}

/**
 * Return the match and non-match intervals of this subalignment WRT the positions in the global MSA.
 * Consensus sequences longer than k are match intervals and the rest as non-match intervals.
 * @param k - the minimum length to consider a vertical stripe as a match interval
 * @return list of intervals
 */
std::vector<Interval> SubAlignment::getMatchAndNonMatchIntervals(uint32_t k) const {
    std::vector<Interval> intervals; //represent match and non-match intervals

    //check if the interval is too short
    if (interval.getLength() < k) {
        //no reason to continue, let's stop here
        Interval shortInterval {interval};
        shortInterval.intervalType = IntervalType::TOO_SHORT;
        intervals.push_back(shortInterval);
        return intervals;
    }

    //the interval is big enough, divide into MATCH and NONMATCH
    //get consensus
    std::string consensusString = buildConsensusString();
    BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: consensusString = " << consensusString;

    //1. get the basic match and non-match intervals with a finite state machine
    enum State {
        BEGIN, MATCH, NONMATCH
    };

    //vars of the finite state machine: currentState, intervalStart, i, c
    State currentState = BEGIN;
    uint32_t intervalStart;
    for (size_t i=0; i<consensusString.size(); ++i) {
        char c = consensusString[i];
        switch (currentState) {
            case BEGIN:
                switch (c) {
                    //setting up initial state according to first char
                    case '*':
                        currentState = NONMATCH;
                        intervalStart=0;
                        break;
                    default:
                        currentState = MATCH;
                        intervalStart=0;
                        break;
                }
                break;
            case MATCH:
                switch (c) {
                    case '*':
                        //end of match interval
                        //saves match interval
                        intervals.push_back(Interval(intervalStart, i, IntervalType::MATCH));
                        BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: found interval: " << intervals.back();

                        //configures next nonmatch interval
                        currentState = NONMATCH;
                        intervalStart=i;
                        break;
                    default:
                        //nothing to do - i is increased
                        break;
                }
                break;
            case NONMATCH:
                switch (c) {
                    case '*':
                        //nothing to do - i is increased
                        break;
                    default:
                        //end of nonmatch interval
                        //saves nonmatch interval
                        intervals.push_back(Interval(intervalStart, i, IntervalType::NONMATCH));
                        BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: found interval: " << intervals.back();

                        //configures next match interval
                        currentState = MATCH;
                        intervalStart=i;
                        break;
                }
                break;
        }
    }



    //2. Fix the basic intervals
    /*
     * A match can become a non-match if the length is < k
     *
     * A non-match can become a match in some cases (remove ):
     * ----AAAA
     * AAAA----
     * ********
     * If we remove the spaces, then it becomes a match interval
     * */
    for (Interval &interval : intervals) {
        switch (interval.intervalType) {
            case IntervalType::MATCH:
                if (interval.getLength() < k) {
                    //not a match interval anymore
                    BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: match interval " << interval << " is too short to be a match interval - transformed to nonmatch";
                    interval.intervalType = IntervalType::NONMATCH;
                }
                break;
            case IntervalType::NONMATCH:
                //check if removing all spaces we have a consensus
                //TODO: leave this for after, marginal case and will require a good amount of lines
                break;
        }
    }

    //3. Merge consecutive intervals with the same type
    //another finite state machine...
    currentState = BEGIN;
    intervalStart = 0;
    std::vector<Interval> joinedIntervals;
    uint32_t i=0;
    for (const Interval &interval : intervals) {
        switch (currentState) {
            case BEGIN:
                switch (interval.intervalType) {
                    //setting up initial state according to first interval
                    case IntervalType::NONMATCH:
                        currentState = NONMATCH;
                        intervalStart=0;
                        break;
                    case IntervalType::MATCH:
                        currentState = MATCH;
                        intervalStart=0;
                        break;
                }
                break;
            case MATCH:
                switch (interval.intervalType) {
                    case IntervalType::NONMATCH:
                        //end of several match intervals
                        //join all previous match intervals
                        joinedIntervals.push_back(Interval(intervals[intervalStart].start, intervals[i-1].end, IntervalType::MATCH));
                        BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: final interval: " << joinedIntervals.back();

                        //configures next nonmatch interval
                        currentState = NONMATCH;
                        intervalStart=i;
                        break;
                    case IntervalType::MATCH:
                        //nothing to do - i is increased
                        break;
                }
                break;
            case NONMATCH:
                switch (interval.intervalType) {
                    case IntervalType::NONMATCH:
                        //nothing to do - i is increased
                        break;
                    case IntervalType::MATCH:
                        //end of several nonmatch intervals
                        //join all previous nonmatch intervals
                        joinedIntervals.push_back(Interval(intervals[intervalStart].start, intervals[i-1].end, IntervalType::NONMATCH));
                        BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: final interval: " << joinedIntervals.back();

                        //configures next match interval
                        currentState = MATCH;
                        intervalStart=i;
                        break;
                }
                break;
        }
        i++;
    }

    return joinedIntervals;
}

std::vector<std::string> SubAlignment::getRepresentativeSequences() const {
    /**
     * 1/ Removes "-" from all alignments
     * 2/ Remove all duplicates
     */
    //get the sequences
    std::vector<std::string> seqs = getSequences();

    //remove all spaces from all seqs
    for (std::string &seq : seqs)
        boost::erase_all(seq, "-");

    //remove all duplicates now
    auto it = std::unique(seqs.begin(), seqs.end());
    seqs.resize(std::distance(seqs.begin(), it));

    return seqs;
}


std::vector<std::string> SubAlignment::getSequences() const {
    std::vector<std::string> sequences;
    sequences.reserve(sequencesNumbers.size());

    //get the sequences
    for (uint32_t sequenceNumber : sequencesNumbers)
        sequences.push_back(MSA->at(sequenceNumber).substr(interval.start, interval.end - interval.start));

    return sequences;
}