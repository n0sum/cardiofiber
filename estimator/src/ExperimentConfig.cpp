#include "ExperimentConfig.h"
#include <fstream>
#include <map>
#include <sstream>

namespace {
    std::string trim(const std::string& s) {
        const auto b = s.find_first_not_of(" \t");
        if (b == std::string::npos) return "";
        const auto e = s.find_last_not_of(" \t");
        return s.substr(b, e - b + 1);
    }

    // "0,2,4" -> {0, 2, 4}
    std::vector<int> parseMapList(const std::string& v) {
        std::vector<int> result;
        std::stringstream ss(v);
        for (std::string tok; std::getline(ss, tok, ',');) result.push_back(std::stoi(tok));
        return result;
    }
}

ExperimentConfig loadExperimentConfig(const std::string& path) {
    std::map<std::string, std::string> kv;
    std::ifstream in(path);
    for (std::string line; std::getline(in, line);) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        const size_t eq = line.find('=');
        kv[trim(line.substr(0, eq))] = trim(line.substr(eq + 1));
    }

    ExperimentConfig c;
    c.experimentId = kv.at("experiment_id");
    c.latPly = kv.at("lat_ply");
    c.groundTruthPly = kv.at("ground_truth_ply");
    c.outputCsv = kv.at("output_csv");
    if (kv.at("maps") == "all") c.allMaps = true;
    else c.maps = parseMapList(kv.at("maps"));
    c.minSingularValue = std::stod(kv.at("mask.min_singular_value"));
    return c;
}
