//
// Created by Leandro Ishi Soares de Lima on 01/05/2019.
//

#include "Subalignment.h"


//builds the consensus string of this subalignment
//TODO: update this !
//TODO: replicate from python code
std::string SubAlignment::buildConsensusString() const {
  std::string consensusString;
  for (size_t j=begin; j<end; ++j) {
    char consensusBase=(*MSA)[sequencesNumbers[0]][j];
    for (size_t i=1; i<sequencesNumbers.size(); ++i) {
      //check
      if (consensusBase!=(*MSA)[sequencesNumbers[i]][j]) { //inconsistency, we are done
        consensusBase = '*';
        break; //go to next column
      }
    }
    consensusString+=consensusBase;
  }
  return consensusString;
}


std::vector<Interval> SubAlignment::getIntervals(uint32_t k) const {
  /**
   *  Returns a list of intervals in which we have
   *  consensus sequences longer than k as match intervals and the rest as non-match intervals.
   *  Mostly inspired by https://github.com/rmcolq/make_prg/blob/master/make_prg_from_msa.py#L117
   * @param k : min_match_length
   */
  auto consensusString = buildConsensusString();

  BOOST_LOG_TRIVIAL(debug) << "@SubAlignment::getIntervals: consensusString = " << consensusString;

  std::vector<Interval> intervals;
  uint32_t matchCount = 0;
  uint32_t matchStart = 0;
  uint32_t nonMatchStart = 0;

  int nbOfNonStarBases = std::count_if(consensusString.begin(), consensusString.end(),
                                        [](char c){ return c != '*'; });


  if (nbOfNonStarBases < k) { //if len(self.consensus.replace('-', '')) < self.min_match_length:
    /* From Rachel:
     * It makes no sense to classify a fully consensus sequence as non-match just because it is too short.
     */
        std::vector<std::string> representativeSequences = getRepresentativeSequences();

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
  }
  else {
      for (size_t i=0; i<consensusString.size(); ++i) {
          if (consensusString[i]!='*') {
              size_t j;
              for (j=i+1; j<consensusString.size() && consensusString[j]!='*'; ++j);
              if (j-i>=k){
                  //new match interval
                  Interval interval(i+begin, j+begin, MATCH);
                  intervals.push_back(interval);
              }else {
                  //too small non-match interval
                  Interval interval(i+begin, j+begin, NONMATCH);
                  intervals.push_back(interval);
              }
              i=j-1;
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


std::vector<std::string> SubAlignment::getRepresentativeSequences() const {
    /**
     * 1/ Removes "-" from all alignments
     * 2/ Disregards sequences with forbidden chars (allowed are ['A','C','G','T','R','Y','K','M','S','W']). Note: 'N' is not allowed
     * 3/ Remove all duplicates
     * 4/ Expands all IUPAC chars
     */

/*
 *     allowed = ['A','C','G','T','R','Y','K','M','S','W']
    iupac = {'R': ['G', 'A'], 'Y': ['T', 'C'], 'K': ['G', 'T'], 'M': ['A', 'C'], 'S': ['G', 'C'], 'W': ['A', 'T']}
    seqs = []
    for s in list(remove_duplicates([str(record.seq).replace('-', '').upper() for record in interval_alignment])): #get all alignments, remove '-', upper(), and remove alignment duplicates (alternative: transforming this list in set() to remove duplicates)
        if contains_only(s, allowed): #check if we have only the allowed chars
            new_seqs = [s]
            for letter in iupac.keys(): #expands iupac into all possible strings (see slacks to remove any doubts on this)
                letter_seqs = []
                for t in new_seqs:
                    if letter in t:
                        letter_seqs.append(t.replace(letter, iupac[letter][0])) #so we first replace all Rs by Gs
                        letter_seqs.append(t.replace(letter, iupac[letter][1])) #and then by As, but shouldn't we do all 2^|R| combinations? #TODO: not sure if it is a bug or not, check this with a counter-example
                    else:
                        letter_seqs.append(t)
                new_seqs = letter_seqs
            seqs.extend(new_seqs)
    ret_list = list(set(seqs))
    if len(ret_list) == 0: #we enter here if all seqs contain at least one N
        print("Every sequence must have contained an N in this slice - redo sequence curation because this is nonsense")
        assert len(ret_list) > 0
    return list(set(seqs))

 */

    static const std::vector<char> allowedBases = {'A','C','G','T','R','Y','K','M','S','W'}; //static so that we don't initialize this over and over again
    static const std::map<char, std::pair<char,char>> translations =    {{'R', {'G', 'A'}},
                                                                        {'Y', {'T', 'C'}},
                                                                        {'K', {'G', 'T'}},
                                                                        {'M', {'A', 'C'}},
                                                                        {'S', {'G', 'C'}},
                                                                        {'W', {'A', 'T'}}};

}


std::vector<std::string> SubAlignment::getSequences() const {
    std::vector<std::string> sequences;
    sequences.reserve(sequencesNumbers.size());

    //get the sequences
    for (uint32_t sequenceNumber : sequencesNumbers)
        sequences.push_back(MSA->at(sequenceNumber).substr(begin, end-begin));

    return sequences;
}