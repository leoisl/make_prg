//
// Created by Leandro Ishi Soares de Lima on 23/02/2019.
//

#include "Utils.h"
#include <pstream.h>
#include <random>
#include <string>

std::string Utils::random_string(std::string::size_type length)
{
    static auto& chrs = "0123456789"
                        "abcdefghijklmnopqrstuvwxyz"
                        "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    thread_local static std::mt19937 rg{std::random_device{}()};
    thread_local static std::uniform_int_distribution<std::string::size_type> pick(0, sizeof(chrs) - 2);

    std::string s;

    s.reserve(length);

    while(length--)
        s += chrs[pick(rg)];

    return s;
}

const std::map<char, std::unordered_set<char>> Utils::baseToAcceptedBases = {
        {'A', {'A'}},
        {'C', {'C'}},
        {'G', {'G'}},
        {'T', {'T'}},
        {'-', {'-'}},
        {'R', {'A', 'G'}},
        {'Y', {'C', 'T'}},
        {'S', {'C', 'G'}},
        {'W', {'A', 'T'}},
        {'K', {'G', 'T'}},
        {'M', {'A', 'C'}},
        {'B', {'C', 'G', 'T'}},
        {'D', {'A', 'G', 'T'}},
        {'H', {'A', 'C', 'T'}},
        {'V', {'A', 'C', 'G'}},
        {'N', {'A', 'C', 'G', 'T'}}
};


void Utils::expandSequencesRecursively(uint32_t i, const std::string &seq, std::list<std::string> &allExpandedSequences) {
    if (i < seq.size()) { //base case
        char originalBase = seq[i]; //get the evaluated base
        if (isNotACGT(originalBase)) {
            //need to expand everyone
            std::list<std::string> newExpandedSequences;
            for (const string &seqToBeExpanded : allExpandedSequences) {
                for (char translatedBase : baseToAcceptedBases.at(originalBase))
                    newExpandedSequences.push_back(seqToBeExpanded + translatedBase);
            }

            //move this new sequences to allExpandedSequences
            allExpandedSequences = std::move(newExpandedSequences);
        } else {
            //originalBase is ACGT
            //add to all seqs
            for (string &seqToBeExpanded : allExpandedSequences)
                seqToBeExpanded += originalBase;
        }

        //recursive call
        expandSequencesRecursively(i+1, seq, allExpandedSequences);
    }
}

std::list<std::string> Utils::expandnonACGT(const std::string &seq) {
    std::list<std::string> allExpandedSequences;
    expandSequencesRecursively(0, seq, allExpandedSequences);
    return allExpandedSequences;
}


void Utils::executeCommand(const string &command, bool verbose, const string &messageIfItFails) {
    // run a process and create a streambuf that reads its stdout and stderr
    if (verbose)
        cerr << "Executing " << command << "..." << endl;

    //create the process
    redi::ipstream proc(command, redi::pstreams::pstdout | redi::pstreams::pstderr);
    string line;

    // read child's stdout
    while (getline(proc.out(), line)) {
        if (verbose)
            cout << line << endl;
    }
    // read child's stderr
    while (getline(proc.err(), line)) {
        if (verbose)
            cerr << line << endl;
    }

    //check exit status
    proc.close();
    if (proc.rdbuf()->exited()) {
        if (proc.rdbuf()->status() != 0) {
            stringstream ss;
            ss << "Error executing " << command << ". Exit status: " << proc.rdbuf()->status() << endl;
            if (messageIfItFails != "")
                ss << "Message: " << messageIfItFails << endl;
            fatalError(ss.str());
        }
        if (verbose)
            cerr << "Executing " << command << " - Done!" << endl;
    }
    else
        fatalError("On executeCommand()");
}