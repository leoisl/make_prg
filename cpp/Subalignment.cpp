//
// Created by Leandro Ishi Soares de Lima on 01/05/2019.
//

#include "Subalignment.h"


std::ostream &operator<<(std::ostream &os, const IntervalType &intervalType) {
    switch (intervalType) {
        case NONMATCH:
            os << "NONMATCH";
            break;
        case MATCH:
            os << "MATCH";
            break;
        case UNDEFINED:
            os << "UNDEFINED";
            break;
    }
    return os;
}


/**
 * Creates a consensus string from the aligment represented by this.
 * Non AGCT symbols RYKMSW result in non-consensus
 * N results in consensus at that position.
 * - is taken into account and we can have a column full of -.
 *
 * @return
 */
std::string SubAlignment::buildConsensusString() const {
    //Doubts:
    //TODO: N is considered as any base, resulting in consensus, but RYKMSW directly results in non consensus. Why? Can't we expand RYKMSW? It makes sense since this is how 'N' is treated.
    //TODO: what about the other IUPAC (B,D,H,V)??
    //TODO: what to call in a column full of non-ACGT (e.g. several N and Rs? Random? non-consensus?)

    //Idea:
    //TODO: To remember: - can be in the consensus string, but never in the PRG
    //TODO: Idea to define consensus string:
    /*
     * AAAAAAAR
     * ANNNNNAR
     * AAAAAAAR
     * ARRRAACR
     * ARRGAGCR
     * AAA*A**?
     * Consensus string should have only ACGT*
     */

    //TODO: we need both methods in fact (buildConsensusString() and getRepresentativeSequences())

    //TODO: unit tests needed here - and compare with the python code


    std::string consensusString;
    for (size_t j = interval.start; j < interval.end; ++j) {
        //1. scans the contents of the column
        std::unordered_set<char> uniqueBasesInThisColumn; //'N' is not considered here
        bool containsRYKMSW=false;

        for (uint32_t sequenceNumber : sequencesNumbers) {
            char base = (*MSA)[sequenceNumber][j];
            if (base != 'N')
                uniqueBasesInThisColumn.insert(base);
            if (base=='R' || base=='Y' || base=='K' || base=='M' || base=='S' || base=='W') {
                containsRYKMSW = true;
                break; //no reason to continue
            }
        }

        //2. set the consensus base
        char consensusBase;
        if (containsRYKMSW) { //no consensus
            consensusBase = '*';
        }else if (uniqueBasesInThisColumn.size() == 0) { //all bases in this column are 'N', consensus is 'N'
            consensusBase = 'N'; //TODO: change this? Should we really consensus to 'N' here?
        }else if (uniqueBasesInThisColumn.size() == 1) { //here we have a consensus - everyone agrees on a unique base
            consensusBase = *uniqueBasesInThisColumn.begin();
        }else { //agreement is on 2+ bases, no consensus
            consensusBase = '*';
        }

        //3. add the consensus base to the consensus string
        consensusString += consensusBase;
    }

    return consensusString;
}

/**
 * Return the match and non-match intervals of this subalignment WRT the positions in the global MSA.
 * Consensus sequences longer than k are match intervals and the rest as non-match intervals.
 * @param k - the minimum length to consider a vertical stripe as a match interval
 * @return vector of intervals
 */
std::vector<Interval> SubAlignment::getMatchAndNonMatchIntervals(uint32_t k) const {
    std::vector<Interval> intervals; //represent match and non-match intervals
    uint32_t matchCount = 0;
    uint32_t matchStart = 0;
    uint32_t nonMatchStart = 0;

    std::string consensusString = buildConsensusString();
    BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getIntervals: consensusString = " << consensusString;


    int nbOfNonSpaceBases = std::count_if(consensusString.begin(), consensusString.end(),
                                         [](char c) { return c != '-'; });


    if (nbOfNonSpaceBases < k) { //if len(self.consensus.replace('-', '')) < self.min_match_length:
        /* From Rachel:
         * It makes no sense to classify a fully consensus sequence as non-match just because it is too short.
         */
        if (consensusString.find('*') != std::string::npos) { //if '*' in self.consensus: - tell us if it is a non-match or not

            //it is probably a non-match, but if we can expand to only one seq, then let's treat this as a match
            //TODO: this should be encoded already in the consensus string - we should not care about this here...
            std::unordered_set<std::string> representativeSequences = getRepresentativeSequences(); //interval_seqs = get_interval_seqs(interval_alignment)
            if (representativeSequences.size() > 1) { //if len(interval_seqs) > 1:
                //non-match confirmed
                Interval newIntervalToAdd(this->getInterval().start, this->getInterval().end, NONMATCH);
                BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getIntervals: adding SHORT NON-MATCH whole interval: " << std::endl << newIntervalToAdd;
                intervals.push_back(newIntervalToAdd); //non_match_intervals.append([0, self.length - 1])
            }else {
                //expanded to only one, this is a match
                Interval newIntervalToAdd(this->getInterval().start, this->getInterval().end, MATCH);
                BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getIntervals: adding SHORT MATCH whole interval: " << std::endl << newIntervalToAdd;
                intervals.push_back(newIntervalToAdd); //match_intervals.append([0, self.length - 1])
            }
        }else { //this is for sure a short match interval
            //TODO: these last two elses are identical, refactor?
            Interval newIntervalToAdd(this->getInterval().start, this->getInterval().end, MATCH);
            BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getIntervals: adding SHORT MATCH whole interval: " << std::endl << newIntervalToAdd;
            intervals.push_back(newIntervalToAdd); //match_intervals.append([0, self.length - 1])
        }
    } else {
        /*
         * TODO: finish add the match and non-match intervals

        for (size_t i = 0; i < consensusString.size(); ++i) {
            if (consensusString[i] != '*') {
                size_t j;
                for (j = i + 1; j < consensusString.size() && consensusString[j] != '*'; ++j);
                if (j - i >= k) {
                    //new match interval
                    Interval interval(i + begin, j + begin, MATCH);
                    intervals.push_back(interval);
                } else {
                    //too small non-match interval
                    Interval interval(i + begin, j + begin, NONMATCH);
                    intervals.push_back(interval);
                }
                i = j - 1;
            }
        }
         */
    }

    return intervals;
}


void SubAlignment::expandRYKMSW(const std::string &seq, std::unordered_set<std::string> &representativeSeqs) const {
    static const std::map<char, std::pair<char, char>> translations = {{'R', {'G', 'A'}},
                                                                       {'Y', {'T', 'C'}},
                                                                       {'K', {'G', 'T'}},
                                                                       {'M', {'A', 'C'}},
                                                                       {'S', {'G', 'C'}},
                                                                       {'W', {'A', 'T'}}};

    //get the positions of the bases that are RYKMSW
    std::vector<size_t> posRYKMSW;
    for (size_t pos = 0; pos < seq.size(); ++pos) {
        char c = seq[pos];
        if (translations.find(c) != translations.end()) //c is RYKMSW
            posRYKMSW.push_back(pos);
    }

    //generate a representative seq where the first option is always chosen
    std::string baseRepresentativeSeq(seq);
    for (char &c : baseRepresentativeSeq) {
        if (translations.find(c) != translations.end()) //c is RYKMSW
            c = translations.at(c).first;
    }

    //generate all subsets from posRYKMSW
    //bases on the subset are switched to the second option
    for (auto &&subset : iter::powerset(posRYKMSW)) {
        //switch
        std::string newSeq(baseRepresentativeSeq);
        for (auto &&pos : subset) {
            newSeq[pos] = translations.at(seq[pos]).second; //seq contains the original sequence, with not replacements
        }

        //save newSeq
        representativeSeqs.insert(newSeq);
    }
}


std::unordered_set<std::string> SubAlignment::getRepresentativeSequences() const {
    /**
     * 1/ Removes "-" from all alignments - TODO: this should be kept on the new version
     * 2/ Disregards sequences with forbidden chars (allowed are ['A','C','G','T','R','Y','K','M','S','W']). Note: 'N' is not allowed - TODO: N should be allowed and other IUPAC also
     * 3/ Remove all duplicates
     * 4/ Expands all IUPAC chars
     */
     //TODO: why not expand 'N' also - in the consensus, 'N' is allowed, here is forbidden (make us ignore the whole sequence in fact)
     //TODO: 'N' is allowed in the consensus but disallowed here, this can be a source of bugs...
     //TODO: cosensus string should represent a consensus of the representative sequences...
    static const std::vector<char> allowedBases = {'A', 'C', 'G', 'T', 'R', 'Y', 'K', 'M', 'S',
                                                   'W'}; //static so that we don't initialize this over and over again

    //get the sequences
    std::vector<std::string> seqs = getSequences();

    //remove all spaces from all seqs
    for (std::string &seq : seqs)
        boost::erase_all(seq, "-");


    //filter out seqs with no allowed bases
    {
        //TODO: use memory better here?
        std::vector<std::string> allowedSeqs;
        allowedSeqs.reserve(seqs.size());
        for (const std::string &seq : seqs) {
            if (std::all_of(seq.begin(), seq.end(),
                    //unary predicator that checks if c is an allowed base
                            [](char c) {
                                return std::find(allowedBases.begin(), allowedBases.end(), c) != allowedBases.end();
                            })) {
                allowedSeqs.push_back(seq);
            } else {
                BOOST_LOG_TRIVIAL(warning)
                    << "Disconsidering the following sequence in SubAlignment::getRepresentativeSequences() due to having non-allowed base:"
                    << std::endl << seq;
            }
        }

        //move this vector to seqs
        seqs = std::move(allowedSeqs);
    }

    //remove all duplicates now
    auto it = std::unique(seqs.begin(), seqs.end());
    seqs.resize(std::distance(seqs.begin(), it));

    //expands RYKMSW and saves all new strings to representativeSeqs
    std::unordered_set<std::string> representativeSeqs;
    for (const std::string &seq : seqs)
        expandRYKMSW(seq, representativeSeqs);

    //final check
    BOOST_ASSERT_MSG(representativeSeqs.size() > 0,
                     "Every sequence must have contained an N in this slice - redo sequence curation because this is nonsense"); //keeping Rachel's nice error message

    return representativeSeqs;
}


std::vector<std::string> SubAlignment::getSequences() const {
    std::vector<std::string> sequences;
    sequences.reserve(sequencesNumbers.size());

    //get the sequences
    for (uint32_t sequenceNumber : sequencesNumbers)
        sequences.push_back(MSA->at(sequenceNumber).substr(interval.start, interval.end - interval.start));

    return sequences;
}