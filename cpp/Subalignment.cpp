//
// Created by Leandro Ishi Soares de Lima on 01/05/2019.
//

#include "Subalignment.h"
#include "Utils.h"
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

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
Consensus SubAlignment::buildConsensusString() const {
    //TODO: issue warning on many non ACGT- bases?
    static const std::string consensusBases{"ACGT-"}; //which bases should we have in the final consensus?

    //generate the consensus string
    Consensus consensusString;
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

        //2. checks if we only have one accepted bases - this should be the case, otherwise sequence curation must be done
        if (acceptedBases.size() > 1) {
            BOOST_LOG_TRIVIAL(fatal) << "There are more than 1 candidate base for column " << j << " of the MSA, please redo sequence curation";
            std::exit(1);
        }

        //3. checks if we have a consensus
        char acceptedBase = '*'; //assumes no consensus
        if (acceptedBases.size() == 1) acceptedBase = acceptedBases[0];

        //3. add the accepted base to the consensus string
        consensusString += acceptedBase;
    }

    return consensusString;
}

/**
 * Return the match and non-match intervals of this subalignment WRT the positions in the global MSA.
 * The intervals returned here can be MATCH, NONMATCH or TOO_SHORT
 * Consensus sequences longer than k are match intervals and the rest as non-match intervals.
 * @param k - the minimum length to consider a vertical stripe as a match interval
 * @return list of intervals that can be MATCH, NONMATCH or TOO_SHORT
 */
std::vector<Interval> SubAlignment::getMatchAndNonMatchIntervals(uint32_t k) const {
    std::vector<Interval> intervals; //represent match and non-match intervals

    //get consensus
    Consensus consensusString = buildConsensusString();
    BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: consensusString" << std::endl << consensusString;

    //check if the interval is too short
    if (consensusString.getSizeWithoutSpaces() < k) { //if len(self.consensus.replace('-', '')) < self.min_match_length:
        //no reason to continue, let's stop here
        Interval shortInterval {interval};
        shortInterval.intervalType = IntervalType::TOO_SHORT;
        BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: Found a TOO SHORT interval: " << std::endl << shortInterval;

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
    uint32_t intervalEnd;
    for (size_t i=0; i<consensusString.size(); ++i) {
        char c = consensusString[i];
        switch (currentState) {
            case BEGIN:
                switch (c) {
                    //setting up initial state according to first char
                    case '*':
                        currentState = NONMATCH;
                        intervalStart=interval.start;
                        break;
                    default:
                        currentState = MATCH;
                        intervalStart=interval.start;
                        break;
                }
                break;
            case MATCH:
                switch (c) {
                    case '*':
                        //end of match interval
                        intervalEnd=i+interval.start;

                        //saves match interval
                        intervals.push_back(Interval(intervalStart, intervalEnd, IntervalType::MATCH));
                        BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: found BASIC interval: " << std::endl
                                                 << intervals.back();

                        //configures next nonmatch interval
                        currentState = NONMATCH;
                        intervalStart=intervalEnd;
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
                        intervalEnd=i+interval.start;

                        //saves nonmatch interval
                        intervals.push_back(Interval(intervalStart, intervalEnd, IntervalType::NONMATCH));
                        BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: found BASIC interval: " << std::endl
                                                 << intervals.back();

                        //configures next match interval
                        currentState = MATCH;
                        intervalStart=intervalEnd;
                        break;
                }
                break;
        }
    }
    //add the last interval
    intervalEnd = interval.start + consensusString.size();
    switch (currentState) {
        case MATCH:
            //saves match interval
            intervals.push_back(Interval(intervalStart, intervalEnd, IntervalType::MATCH));
            BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: found BASIC interval: " << std::endl
                                     << intervals.back();
            break;
        case NONMATCH:
            //saves match interval
            intervals.push_back(Interval(intervalStart, intervalEnd, IntervalType::NONMATCH));
            BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: found BASIC interval: " << std::endl
                                     << intervals.back();
            break;
    }



    //2. Fix the basic intervals
    /*
     * A match can become a non-match if the length is < k
     *
     * A non-match can become a match in some cases eg:
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
                    BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: " << std::endl << "MATCH interval " << interval << " is too short to be a match interval - transformed to NONMATCH";
                    interval.intervalType = IntervalType::NONMATCH;
                }
                break;
            case IntervalType::NONMATCH:
                //check if removing all spaces we have a consensus
                //TODO: leave this for after, marginal case and will require a good amount of lines
                BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: " << std::endl << "WARNING: NOT YET IMPLEMENTED";
                break;
        }
    }

    //3. Merge consecutive intervals with the same type
    BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals - MERGING intervals...";
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
                        BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: FINAL interval: " << joinedIntervals.back();

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
                        BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: FINAL interval: " << joinedIntervals.back();

                        //configures next match interval
                        currentState = MATCH;
                        intervalStart=i;
                        break;
                }
                break;
        }
        i++;
    }
    //merge the last interval
    switch (currentState) {
        case MATCH:
            //join all previous match intervals
            joinedIntervals.push_back(Interval(intervals[intervalStart].start, intervals.back().end, IntervalType::MATCH));
            BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: FINAL interval: " << joinedIntervals.back();
            break;
        case NONMATCH:
            //join all previous nonmatch intervals
            joinedIntervals.push_back(Interval(intervals[intervalStart].start, intervals.back().end, IntervalType::NONMATCH));
            BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getMatchAndNonMatchIntervals: FINAL interval: " << joinedIntervals.back();
            break;
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
 * Main clustering algorithm
 * @param seqWithNoSpace2seqNbsBig
 * @param k
 * @return
 */
std::vector< std::vector<const std::string *>> SubAlignment::kMeansCluster(const std::vector<const std::string *> &bigSequences, int k) const {
    //transform sequences into kmer occurance vectors
    BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::kMeansCluster: transforming sequences into kmer occurance vectors";

    //TODO: this should be done better, but I don't know how yet... we should use boost MPL
    //TODO: everytime KSIZE_LIST changes, this should be changed
    typedef boost::variant <
            KMeansClusterer<KMER_SPAN(0)>,
            KMeansClusterer<KMER_SPAN(1)>,
            KMeansClusterer<KMER_SPAN(2)>,
            KMeansClusterer<KMER_SPAN(3)>
    >  KMeansClustererVariant;
    KMeansClustererVariant kMeansClustererVariant;
    if (k < KMER_SPAN(0))
        kMeansClustererVariant = KMeansClusterer<KMER_SPAN(0)>(&bigSequences, k);
    else if (k < KMER_SPAN(1))
        kMeansClustererVariant = KMeansClusterer<KMER_SPAN(1)>(&bigSequences, k);
    else if (k < KMER_SPAN(2))
        kMeansClustererVariant = KMeansClusterer<KMER_SPAN(2)>(&bigSequences, k);
    else if (k < KMER_SPAN(3))
        kMeansClustererVariant = KMeansClusterer<KMER_SPAN(3)>(&bigSequences, k);
    else
        throw gatb::core::system::Exception ("Subalignment::kMeansCluster failure because of unhandled kmer size %d", k);

    boost::apply_visitor(KMeansClustererVisitor(), kMeansClustererVariant);
}


/**
 * Split this subalignment into several subaligments, where each is a cluster of similar sequences
 */
SubAlignment::Clusters SubAlignment::kMeansCluster(uint32_t k) const {
    BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::kMeansCluster: clustering " << *this;

    Clusters clusters;

    //get a random filename
    std::string filename = Utils::random_string(64);
    std::string out_filename = filename+".out";

    //create the file with the MSA to be clustered
    {
        std::ofstream msaFile;
        Utils::openFileForWriting(filename, msaFile);
        for (uint32_t sequenceNumber : sequencesNumbers)
            msaFile << ">" << sequenceNumber << std::endl <<
                    MSA->at(sequenceNumber).substr(interval.start, interval.end - interval.start) << std::endl;
        msaFile.close();
    }

    //cluster the sequences using python script for the moment
    //TODO: fix this: use c++
    {
        std::stringstream ss;
        ss << "../venv/bin/python cluster.py " << filename << " " << k;
        Utils::executeCommand(ss.str());
    }

    //create the clusters
    {
        //read the text output file
        auto clustersAsText = Utils::getVectorStringFromFile(out_filename);


        for (const auto &clusterAsText : clustersAsText) {
            //get the sequences in this cluster
            std::vector<uint32_t> sequencesNumbersInTheCluster;
            uint32_t sequenceNumber;
            stringstream ss;
            ss << clusterAsText;
            while (ss >> sequenceNumber) {
                sequencesNumbersInTheCluster.push_back(sequenceNumber);
            }

            //add the cluster
            Interval newInterval(interval.start, interval.end, IntervalType::UNPROCESSED); //new cluster is unprocessed - we will process is recursively
            clusters.push_back(SubAlignment(sequencesNumbersInTheCluster, newInterval, MSA));
        }
    }

    //clean up
    boost::filesystem::remove(filename);
    boost::filesystem::remove(out_filename);


    return clusters;


    /*
     * TODO: old clustering code, keeping here since we might need it
    //get the aligments without '-' with their IDs
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
    vector<const std::string *> bigSequences;
    for (const auto &[seq, dontcare] : seqWithNoSpace2seqNbsBig)
        bigSequences.push_back(seq);

    if (seqWithNoSpace2seqNbsBig.size()>1) {
        kMeansCluster(bigSequences, k);
    }else {

    }


    //each seqWithNoSpace2seqNbsTooShort becomes a cluster

    //decompress to get the clusters
     */
}


std::vector<std::string> SubAlignment::getSequences() const {
    std::vector<std::string> sequences;
    sequences.reserve(sequencesNumbers.size());

    //get the sequences
    for (uint32_t sequenceNumber : sequencesNumbers)
        sequences.push_back(MSA->at(sequenceNumber).substr(interval.start, interval.end - interval.start));

    return sequences;
}

