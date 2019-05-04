//
// Created by Leandro Ishi Soares de Lima on 01/05/2019.
//

#include "Subalignment.h"


//builds the consensus string of this subalignment
//TODO: update this !
//TODO: replicate from python code
std::string SubAlignment::buildConsensusString() const {
    std::string consensusString;
    for (size_t j = begin; j < end; ++j) {
        char consensusBase = (*MSA)[sequencesNumbers[0]][j];
        for (size_t i = 1; i < sequencesNumbers.size(); ++i) {
            //check
            if (consensusBase != (*MSA)[sequencesNumbers[i]][j]) { //inconsistency, we are done
                consensusBase = '*';
                break; //go to next column
            }
        }
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
        if (consensusString.find('*') != std::string::npos) { //if '*' in self.consensus:
            std::set<std::string> representativeSequences = getRepresentativeSequences(); //interval_seqs = get_interval_seqs(interval_alignment)
            if (representativeSequences.size() > 1) { //if len(interval_seqs) > 1:
                BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getIntervals: adding short non-match whole interval: " << std::endl <<
                                         *this;
                intervals.push_back(Interval);
            }

        }



        /**
                if '*' in self.consensus:
                    interval_alignment = self.alignment[:, 0:self.length]
                    interval_seqs = get_interval_seqs(interval_alignment)
                    if len(interval_seqs) > 1:
                        logging.debug("add short non-match whole interval [%d,%d]" %(0,self.length - 1))
                        non_match_intervals.append([0, self.length - 1])
                    else:
                        logging.debug("add short match whole interval [%d,%d]" %(0,self.length - 1))
                        match_intervals.append([0, self.length - 1])
                else:
                    match_intervals.append([0, self.length - 1])
                    logging.debug("add short match whole interval [%d,%d]" % (0, self.length - 1))

         */
        if (consensusString.find('*') != std::string::npos) { //if '*' in self.consensus:

        }
    } else {
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
    }

    return intervals;


/*

              if len(self.consensus.replace('-', '')) < self.min_match_length:
              # It makes no sense to classify a fully consensus sequence as
              # a non-match just because it is too short.
              if '*' in self.consensus:
              interval_alignment = self.alignment[:, 0:self.length]
              interval_seqs = get_interval_seqs(interval_alignment)
              if len(interval_seqs) > 1:
              logging.debug("add short non-match whole interval [%d,%d]" %(0,self.length - 1))
              non_match_intervals.append([0, self.length - 1])
              else:
              logging.debug("add short match whole interval [%d,%d]" %(0,self.length - 1))
              match_intervals.append([0, self.length - 1])
              else:
              match_intervals.append([0, self.length - 1])
              logging.debug("add short match whole interval [%d,%d]" % (0, self.length - 1))
              else:
              for i in range(self.length):
              letter = self.consensus[i]
              if letter != '*':
              # In a match region.
              if match_count == 0:
              match_start = i
              match_count += 1
              elif match_count > 0:
              # Have reached a non-match. Check if previous match string is long enough to add to match_regions
              match_string = self.consensus[match_start: match_start + match_count].replace('-', '') #Leandro: remove "-" due to subalignments
              match_len = len(match_string)
              logging.debug("have match string %s" % match_string)

              if match_len >= self.min_match_length:
              # if the non_match sequences in the interval are really the same, add a match interval
              interval_alignment = self.alignment[:, non_match_start:match_start + 1] #interval alignment is the subalignment from the last non-match until the first match (the edges in my model)
                        interval_seqs = get_interval_seqs(interval_alignment)
                        if non_match_start < match_start and len(interval_seqs) > 1:
                            non_match_intervals.append([non_match_start, match_start - 1])
                            logging.debug("add non-match interval as have alts [%d,%d]"
                                          % (non_match_start, match_start - 1))
                        elif non_match_start < match_start:
                            match_intervals.append([non_match_start, match_start - 1])
                            logging.debug("add match interval as only one seq [%d,%d]"
                                          % (non_match_start, match_start - 1))
                        match_intervals.append([match_start, match_start + match_count - 1])
                        logging.debug("add match interval to complete step [%d,%d]"
                                      % (match_start, match_start + match_count- 1))
                        non_match_start = i
                    match_count = 0
                    match_start = non_match_start

            # At end add last intervals
            match_string = self.consensus[match_start: match_start + match_count].replace('-', '')
            match_len = len(match_string)
            logging.debug("at end have match string %s" % match_string)
            if 0 < match_len < self.min_match_length:
                logging.debug("have short match region at end, so include it in non-match-region before - "
                              "match count was %d" %match_count)
                match_count = 0
                match_start = non_match_start
                logging.debug("match count is now %d" % match_count)

            if match_count > 0:
                interval_alignment = self.alignment[:, non_match_start:match_start + 1]
            else:
                interval_alignment = self.alignment[:, non_match_start:self.length]
            interval_seqs = get_interval_seqs(interval_alignment)
            if len(interval_seqs) == 1:
                match_intervals.append([non_match_start, self.length - 1])
                logging.debug("add match interval at end as only one seq [%d,%d]" % (non_match_start, self.length - 1))
            elif len(interval_seqs) > 1 and non_match_start < match_start:
                non_match_intervals.append([non_match_start, match_start - 1])
                logging.debug("add non-match interval at end as have alts [%d,%d]" % (non_match_start, match_start - 1))
                match_intervals.append([match_start, self.length - 1])
                logging.debug("add match interval at end [%d,%d]" % (match_start, self.length - 1))
            else:
                non_match_intervals.append([non_match_start, self.length - 1])
                logging.debug("add only non-match interval at end as have alts [%d,%d]" % (non_match_start, self.length - 1))

        # check all stretches of consensus are in an interval, and intervals don't overlap
        for i in range(self.length):
            count_match = 0
            for interval in match_intervals:
                if interval[0] <= i <= interval[1]:
                    count_match += 1
            count_non_match = 0
            for interval in non_match_intervals:
                if interval[0] <= i <= interval[1]:
                    count_non_match += 1

            assert (count_match | count_non_match), "Failed to correctly identify match intervals: position %d " \
                                                    "appeared in both/neither match and non-match intervals" % i
            assert (count_match + count_non_match == 1), "Failed to correctly identify match intervals: position " \
                                                         "%d appeared in %d intervals" % (
                                                             i, count_match + count_non_match)

        return match_intervals, non_match_intervals
*/
}


void SubAlignment::expandRYKMSW(const std::string &seq, std::set<std::string> &representativeSeqs) const {
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


std::set<std::string> SubAlignment::getRepresentativeSequences() const {
    /**
     * 1/ Removes "-" from all alignments
     * 2/ Disregards sequences with forbidden chars (allowed are ['A','C','G','T','R','Y','K','M','S','W']). Note: 'N' is not allowed
     * 3/ Remove all duplicates
     * 4/ Expands all IUPAC chars
     */

/*
    if len(ret_list) == 0: #we enter here if all seqs contain at least one N
        print("Every sequence must have contained an N in this slice - redo sequence curation because this is nonsense")
        assert len(ret_list) > 0
    return list(set(seqs))

 */

    static const std::vector<char> allowedBases = {'A', 'C', 'G', 'T', 'R', 'Y', 'K', 'M', 'S',
                                                   'W'}; //static so that we don't initialize this over and over again

    //get the sequences
    std::vector<std::string> seqs = getSequences();

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

    //remove all spaces from all seqs
    for (std::string &seq : seqs)
        boost::erase_all(seq, "-");

    //remove all duplicates now
    auto it = std::unique(seqs.begin(), seqs.end());
    seqs.resize(std::distance(seqs.begin(), it));

    //expands RYKMSW and saves all new strings to representativeSeqs
    std::set<std::string> representativeSeqs;
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
        sequences.push_back(MSA->at(sequenceNumber).substr(begin, end - begin));

    return sequences;
}