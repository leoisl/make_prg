//
// Created by Leandro Ishi Soares de Lima on 22/02/2019.
//

#ifndef MAKE_PRG_PRGBOOSTGRAPH_H
#define MAKE_PRG_PRGBOOSTGRAPH_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <string>
#include <vector>
#include <cstdint>
#define NULL_VERTEX_ID UINT32_MAX

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
    std::string buildConsensusString() const {
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

    void addSequence(uint32_t sequenceNumber) {
      sequencesNumbers.push_back(sequenceNumber);
    }

    //get the match and non-match intervals WRT the positions in the global MSA
    std::vector<Interval> getIntervals(uint32_t k) const {
      /**
       * @param k : min_match_length
       */
      auto consensusString = buildConsensusString();
      std::vector<Interval> intervals;
      for (size_t i=0; i<consensusString.size(); ++i) {
        if (consensusString[i]!='*') {
          size_t j;
          for (j=i+1; j<consensusString.size() && consensusString[j]!='*'; ++j);
          if (j-i>=k){
            //new match interval
            Interval interval(i+begin, j+begin, MATCH);
            intervals.push_back(interval);
          }else {
            //new non-match interval
            //new match interval
            Interval interval(i+begin, j+begin, NONMATCH);
            intervals.push_back(interval);
          }
          i=j-1;
        }
      }
      return intervals;
    }

    uint32_t getBegin() const { return begin; }
    uint32_t getEnd() const { return end; }
    const std::vector<uint32_t>& getSequencesNumbers() const { return sequencesNumbers; }


};

//vertex properties
struct VertexInfo {
public:
    //the subalignment represented by this vertex
    SubAlignment subalignments;
    IntervalType intervalType;
    VertexInfo() = default;
    //TODO: remove this, boost always use default constructor
    /*
    VertexInfo (const SubAlignment &subalignments, VertexType, vertexType) :
        subalignments{subalignments}, vertexType{vertexType} {}
    */
};

//edge properties
struct EdgeInfo {
    //empty for now
    EdgeInfo() = default;
};

//we want to add the index properties to vertices and edges
typedef boost::property<boost::vertex_index_t, uint32_t, VertexInfo> VertexProps;
typedef boost::property<boost::edge_index_t, uint32_t, EdgeInfo> EdgeProps;

//declare the Boost graph we want
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, VertexProps, EdgeProps> BoostGraph;

#endif //MAKE_PRG_PRGBOOSTGRAPH_H
