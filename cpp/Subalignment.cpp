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


//get the match and non-match intervals WRT the positions in the global MSA
std::vector<Interval> SubAlignment::getIntervals(uint32_t k) const {
  /**
   *  Returns a list of intervals in which we have
   *  consensus sequences longer than k as match intervals and the rest as non-match intervals
   * @param k : min_match_length
   */
  auto consensusString = buildConsensusString();



  std::vector<Interval> intervals;
  uint32_t matchCount = 0;
  uint32_t matchStart = 0;
  uint32_t nonMatchStart = 0;

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
  return intervals;


/*

              logging.debug("consensus: %s" %self.consensus)
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