//
// Created by Leandro Ishi Soares de Lima on 23/02/2019.
//

#ifndef MAKE_PRG_UTILS_H
#define MAKE_PRG_UTILS_H

#include "includes.h"

//utility functions (auto explanatory)
class Utils {
public:
    static void fatalError(const std::string &message) {
        std::cerr << std::endl << std::endl << "[FATAL ERROR] " << message << std::endl << std::endl;
        std::cerr.flush();
        std::exit(1);
    }


    static void openFileForReading(const std::string &filePath, std::ifstream &stream) {
        stream.open(filePath);
        if (!stream.is_open()) {
            std::stringstream ss;
            ss << "Error opening file " << filePath;
            Utils::fatalError(ss.str());
        }
    }

    static void openFileForWriting(const std::string &filePath, std::ofstream &stream) {
        stream.open(filePath);
        if (!stream.is_open()) {
            stringstream ss;
            ss << "Error opening file " << filePath;
            fatalError(ss.str());
        }
    }

    static void executeCommand(const string &command, bool verbose=true, const string &messageIfItFails="");


    //Read all strings in the readsFile file and return them as a vector of strings
    static std::vector<std::string> getVectorStringFromFile(const std::string &readsFile) {
        std::vector<std::string> allReadFilesNames;
        std::string tempStr;

        std::ifstream readsFileStream;
        Utils::openFileForReading(readsFile, readsFileStream);
        while (getline(readsFileStream, tempStr)) {
            if (tempStr.size() > 0)
                allReadFilesNames.push_back(tempStr);
        }
        readsFileStream.close();

        return allReadFilesNames;
    }


    //IUPAC dealing
    static const std::map<char, std::unordered_set<char>> baseToAcceptedBases;

    /**
     * Check if the MSABase accepts the candidateBase (see https://www.bioinformatics.org/sms/iupac.html)
     * @param MSABase - a base of the MSA
     * @param candidateBase - a candidate base (ACGT-)
     * @return
     */
    static inline bool accepts(char MSABase, char candidateBase) {
        return baseToAcceptedBases.at(MSABase).count(candidateBase) > 0;
    }

    static inline bool isNotACGT (char c) {
        return c!='A' && c!='C' && c!='G' && c!='T';
    }

    /**
    * Expands all non-ACGT bases on seq, returning all expanded sequences
    */
private: static void expandSequencesRecursively(uint32_t i, const std::string &seq, std::list<std::string> &allExpandedSequences);
public:
    static std::list<std::string> expandnonACGT(const std::string &seq);
};


#endif //MAKE_PRG_UTILS_H
