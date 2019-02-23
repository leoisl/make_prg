//
// Created by Leandro Ishi Soares de Lima on 23/02/2019.
//

#ifndef MAKE_PRG_UTILS_H
#define MAKE_PRG_UTILS_H


#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>


//utility functions (auto explanatory)
class Utils {
public:
    static void fatalError (const std::string &message) {
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

};


#endif //MAKE_PRG_UTILS_H
