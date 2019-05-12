//
// Created by Leandro Ishi Soares de Lima on 01/05/2019.
//

#ifndef MAKE_PRG_SUBALIGNMENT_H
#define MAKE_PRG_SUBALIGNMENT_H

#include "includes.h"
#include "Utils.h"
#include "BooMap.hpp"

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

    //TODO: document
    void kMeansCluster(const std::unordered_map<const std::string *, std::vector<uint32_t>> &seqWithNoSpace2seqNbsBig, int k) const;

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




//Build a map where the key is all present k-mers and the value is their occurance,
//Using GATB's 2-bit encoding (thus the template)
template<size_t span>
class KmerOcurranceBuilder {
private:
    const std::unordered_map<const std::string *, std::vector<uint32_t>> *seqWithNoSpace2seqNbsBig;
    uint32_t k;
    typename Kmer<span>::ModelDirect model;
public:
    KmerOcurranceBuilder() = default;
    KmerOcurranceBuilder(const std::unordered_map<const std::string *, std::vector<uint32_t>> *seqWithNoSpace2seqNbsBig,
                         uint32_t k) :
            seqWithNoSpace2seqNbsBig{seqWithNoSpace2seqNbsBig}, k{k}, model{k} {}
    void build() const
    {
        //TODO: I thought about using BooMap but it is not the good structure here because it works on a fixed set of k-mers
        //TODO: more in general, BooMap, GATB's MPHF, BLight, etc... all work on a fixed set of k-mers
        //TODO: for cheap updates afterwards, a hash map allowing for updates is better
        //BooMap<typename Kmer<span>::Type , double> kmerToOccurance; //BooMap to hash
        std::map<typename Kmer<span>::Type , double> kmerToOccurance; //using normal std::map. TODO: improve?

        for (const auto &[seq, dontcare] : *seqWithNoSpace2seqNbsBig) {
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
                    kmerToOccurance[expandedKmerEncoded.value()] += occurance;
                }
            }
        }

        //Now, cluster
/*
                # cluster sequences using kmeans
                logging.debug("Now cluster:")
                kmeans = KMeans(n_clusters=1, random_state=2).fit(seq_kmer_counts)
                pre_cluster_inertia = kmeans.inertia_

                if pre_cluster_inertia == 0:
                    logging.debug("pre_cluster_intertia is 0!")
                    for key in list(interval_seq_dict.keys()):
                        logging.debug("seq: %s, num_seqs with this seq: %d", key, len(interval_seq_dict[key]))

                cluster_inertia = pre_cluster_inertia
                number_of_clusters = 1
                logging.debug("number of clusters: %d, inertia: %f", number_of_clusters, cluster_inertia)
                while (cluster_inertia > 0
                       and cluster_inertia > pre_cluster_inertia / 2 #we cluster until we reach less than half of the initial cluster inertia
                       and number_of_clusters <= len(interval_seqs)):
                    number_of_clusters += 1
                    kmeans = KMeans(n_clusters=number_of_clusters, random_state=2).fit(seq_kmer_counts)
                    cluster_inertia = kmeans.inertia_
                    logging.debug("number of clusters: %d, inertia: %f", number_of_clusters, cluster_inertia)

                # now extract the equivalence class details from this partition and return
                logging.debug("Extract equivalence classes from this partition")
                if pre_cluster_inertia > 0:
                    equiv_class_ids = list(kmeans.predict(seq_kmer_counts))
                    for i in range(max(equiv_class_ids) + 1):
                        big_return_id_lists.append([])
                    for i, val in enumerate(equiv_class_ids):
                        big_return_id_lists[val].extend(interval_seq_dict[interval_seqs[i]])
                else:
                    logging.debug("default to not clustering")
                    big_return_id_lists = [interval_seq_dict[key] for key in interval_seq_dict.keys()]
            elif len(interval_seqs) == 1:
                big_return_id_lists = [interval_seq_dict[interval_seqs[0]]]

 */
    }
};


/**
 * Visitor that will operate in a model to build the kmer occurance map
 */
class KmerOccuranceBuilderVisitor : public boost::static_visitor<> {
public:
    template<class T>
    void operator()(T& kmerOccuranceBuilder) const {
        kmerOccuranceBuilder.build();
    }
};



#endif //MAKE_PRG_SUBALIGNMENT_H
