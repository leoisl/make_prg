//
// Created by Leandro Ishi Soares de Lima on 01/05/2019.
//

#include "Subalignment.h"

std::ostream &operator<<(std::ostream &os, const IntervalType &intervalType) {
    switch (intervalType) {
        case UNPROCESSED:
            os << "UNPROCESSED";
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

    //get consensus
    std::string consensusString = buildConsensusString();
    BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: consensusString = " << consensusString;

    //check if the interval is too short
    //if (interval.getLength() < k) { //at first I did like this, but Rachel's condition is below better
    if (boost::erase_all_copy(consensusString, "-").size() < k) { //if len(self.consensus.replace('-', '')) < self.min_match_length:
        //no reason to continue, let's stop here
        Interval shortInterval {interval};
        shortInterval.intervalType = IntervalType::TOO_SHORT;
        BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: Found a TOO SHORT interval: " << shortInterval;

        intervals.push_back(shortInterval);
        return intervals;
    }


    //the interval is big enough, divide into MATCH and NONMATCH
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


/**
 * Visitor that will operate in a model to build the kmer occurance map
 */
//TODO: everytime KSIZE_LIST changes, this should be changed
/*typedef boost::variant <
        BooMap<Kmer<KMER_SPAN(0)>::Type, double>,
        BooMap<Kmer<KMER_SPAN(1)>::Type, double>,
        BooMap<Kmer<KMER_SPAN(2)>::Type, double>,
        BooMap<Kmer<KMER_SPAN(3)>::Type, double>
> BooMapVariant;
class KmerOcurranceVisitor : public boost::static_visitor<BooMapVariant> {
private:
    const std::unordered_map<const std::string *, std::vector<uint32_t>> &seqWithNoSpace2seqNbsBig;
    int k;
public:
    KmerOcurranceVisitor (const std::unordered_map<const std::string *, std::vector<uint32_t>> &seqWithNoSpace2seqNbsBig, int k) :
            seqWithNoSpace2seqNbsBig{seqWithNoSpace2seqNbsBig}, k{k} {}

    template<typename T>
    struct extract_value_type //lets call it extract_value_type
    {
        typedef T value_type;
    };

    template<template<typename> class X, typename T>
    struct extract_value_type<X<T>>   //specialization
    {
        typedef T value_type;
    };


    template<class Model>
    BooMapVariant  operator() (Model &model) const
    {
        BooMap<typename Model::Kmer::value, double> kmerToOccurance; //BooMap to hash
        for (const auto &[seq, dontcare] : seqWithNoSpace2seqNbsBig) {
            for (uint32_t i=0; i <= seq->size()-k; ++i){
                // we get the k-mer
                std::string kmer = seq->substr(i, k);

                //we expand the non-ACGT bases to all possible translations of this k-mer
                std::list<std::string> expandedKmers = Utils::expandnonACGT(kmer);

                //the occurance is 1/nb of expanded kmers (i.e. if we have a N in the kmer, each expanded kmer has an occurence of 0.25)
                double occurance = 1.0 / ((double)expandedKmers.size());

                //now we go through all expanded kmers
                for (const string &expandedKmer : expandedKmers) {
                    //here we have just ACGT, so we can use GATB's 2-bit encoding
                    auto expandedKmerEncoded = model.codeSeed (expandedKmer.c_str(), Data::ASCII);

                    //add to the map
                    kmerToOccurance.add(expandedKmerEncoded, occurance);
                }
            }
        }
        return kmerToOccurance;
    }
};*/


void SubAlignment::kMeansCluster(const std::unordered_map<const std::string *, std::vector<uint32_t>> &seqWithNoSpace2seqNbsBig, int k) const {
    //transform sequences into kmer occurance vectors
    BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::kMeansCluster: transforming sequences into kmer occurance vectors";

    //TODO: this should be done better, but I don't know how yet... we should use boost MPL
    //TODO: everytime KSIZE_LIST changes, this should be changed
    typedef boost::variant <
            Kmer<KMER_SPAN(0)>::ModelDirect,
            Kmer<KMER_SPAN(1)>::ModelDirect,
            Kmer<KMER_SPAN(2)>::ModelDirect,
            Kmer<KMER_SPAN(3)>::ModelDirect
    >  ModelDirectVariant;

    ModelDirectVariant model;
    if (k < KMER_SPAN(0))  {  model = Kmer<KMER_SPAN(0)>::ModelDirect(k); }
    else if (k < KMER_SPAN(1))  {  model = Kmer<KMER_SPAN(1)>::ModelDirect(k); }
    else if (k < KMER_SPAN(2))  {  model = Kmer<KMER_SPAN(2)>::ModelDirect(k); }
    else if (k < KMER_SPAN(3))  {  model = Kmer<KMER_SPAN(3)>::ModelDirect(k); }
    else { throw gatb::core::system::Exception ("Subalignment::kMeansCluster failure because of unhandled kmer size %d", k); }
    //auto kmerToOccurance = boost::apply_visitor (KmerOcurranceVisitor(seqWithNoSpace2seqNbsBig,k),  model);
}


/**
     * Split this subalignment into several subaligments, where each is a cluster of similar sequences
     * @return a vector of subalignments
     */
std::vector<SubAlignment> SubAlignment::kMeansCluster(uint32_t k) const {
    BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::kMeansCluster: clustering " << *this;

    //get the aligments without - with their IDs
    std::unordered_map<uint32_t, std::string> seqNb2seqWithNoSpace;
    //get the sequences without space
    for (uint32_t sequenceNumber : sequencesNumbers)
        seqNb2seqWithNoSpace[sequenceNumber] = boost::erase_all_copy(MSA->at(sequenceNumber).substr(interval.start, interval.end - interval.start), "-");

    //divide the sequences into two sets: tooShort (<k) and big (>=k)
    //also, we will only work on unique sequences, remembering their original numbers
    std::unordered_map<const std::string *, std::vector<uint32_t>> seqWithNoSpace2seqNbsTooShort, seqWithNoSpace2seqNbsBig;
    for (const auto &[seqNb, seqWithNoSpace] : seqNb2seqWithNoSpace) {
        if (seqWithNoSpace.size() < k)
            seqWithNoSpace2seqNbsTooShort[&seqWithNoSpace].push_back(seqNb);
        else
            seqWithNoSpace2seqNbsBig[&seqWithNoSpace].push_back(seqNb);
    }

    //cluster seqWithNoSpace2seqNbsBig
    //??? std::vector<std>
    if (seqWithNoSpace2seqNbsBig.size()>1) {
        //kMeansCluster(seqWithNoSpace2seqNbsBig, k);
    }else {

    }


    //each seqWithNoSpace2seqNbsTooShort becomes a cluster

    //decompress to get the clusters
}


std::vector<std::string> SubAlignment::getSequences() const {
    std::vector<std::string> sequences;
    sequences.reserve(sequencesNumbers.size());

    //get the sequences
    for (uint32_t sequenceNumber : sequencesNumbers)
        sequences.push_back(MSA->at(sequenceNumber).substr(interval.start, interval.end - interval.start));

    return sequences;
}

