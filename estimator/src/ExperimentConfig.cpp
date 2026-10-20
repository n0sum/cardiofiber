#include "ExperimentConfig.h"
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>

namespace {
    std::string trim(const std::string& s) {
        const auto b = s.find_first_not_of(" \t");
        if (b == std::string::npos) return "";
        const auto e = s.find_last_not_of(" \t");
        return s.substr(b, e - b + 1);
    }

    std::vector<int> parseMapList(const std::string& v) {
        std::vector<int> result;
        size_t start = 0;
        while (true) {
            const size_t comma = v.find(',', start);
            result.push_back(std::stoi(v.substr(start, comma == std::string::npos ? std::string::npos : comma - start)));
            if (comma == std::string::npos) break;
            start = comma + 1;
        }
        return result;
    }
}

ExperimentConfig loadExperimentConfig(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error(path + ": impossibile aprire il file di esperimento");

    const std::set<std::string> known = {
        "experiment_id", "lat_ply", "ground_truth_ply", "region_cv_csv", "maps",
        "mask.min_singular_value", "mask.collision_angle_deg", "output_csv", "output_mask_csv"};
    std::map<std::string, std::string> kv;

    std::string line;
    for (int lineNo = 1; std::getline(in, line); ++lineNo) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        const size_t eq = line.find('=');
        const std::string key = trim(line.substr(0, eq));
        if (!known.contains(key))
            throw std::runtime_error(path + ":" + std::to_string(lineNo) + ": chiave non riconosciuta: '" + key + "'");
        kv[key] = trim(line.substr(eq + 1));
    }

    for (const auto& k : known)
        if (!kv.contains(k)) throw std::runtime_error(path + ": chiave obbligatoria mancante: '" + k + "'");

    ExperimentConfig c;
    c.sourcePath = path;
    c.experimentId = kv["experiment_id"];
    c.latPly = kv["lat_ply"];
    c.groundTruthPly = kv["ground_truth_ply"];
    c.regionCvCsv = kv["region_cv_csv"];
    c.outputCsv = kv["output_csv"];
    c.outputMaskCsv = kv["output_mask_csv"];
    if (kv["maps"] == "all") c.allMaps = true;
    else c.maps = parseMapList(kv["maps"]);
    c.masks.minSingularValue = std::stod(kv["mask.min_singular_value"]);
    c.masks.collisionAngleDeg = std::stod(kv["mask.collision_angle_deg"]);
    return c;
}
