#include "ExperimentConfig.h"
#include <fstream>
#include <map>
#include <set>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace {
    std::string trim(const std::string& s) {
        const auto b = s.find_first_not_of(" \t");
        if (b == std::string::npos) return "";
        const auto e = s.find_last_not_of(" \t");
        return s.substr(b, e - b + 1);
    }

    std::string where(const std::string& path, int line) {
        return path + ":" + std::to_string(line) + ": ";
    }

    double parseDouble(const std::string& v, const std::string& at) {
        size_t pos = 0;
        double x;
        try { x = std::stod(v, &pos); }
        catch (const std::exception&) { throw std::runtime_error(at + "numero non valido: '" + v + "'"); }
        if (pos != v.size() || !std::isfinite(x)) throw std::runtime_error(at + "numero non valido: '" + v + "'");
        return x;
    }

    std::vector<int> parseMapList(const std::string& v, const std::string& at) {
        std::vector<int> result;
        std::set<int> seen;
        size_t start = 0;
        while (true) {
            const size_t comma = v.find(',', start);
            const std::string item = trim(v.substr(start, comma == std::string::npos ? std::string::npos : comma - start));
            size_t pos = 0;
            int k = -1;
            try { k = std::stoi(item, &pos); }
            catch (const std::exception&) { pos = 0; }
            if (item.empty() || pos != item.size() || k < 0)
                throw std::runtime_error(at + "indice di mappa non valido: '" + item + "'");
            if (!seen.insert(k).second)
                throw std::runtime_error(at + "indice di mappa " + std::to_string(k) + " ripetuto");
            result.push_back(k);
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
        "experiment_id", "description", "lat_ply", "ground_truth_ply", "region_cv_csv",
        "maps", "mask.min_singular_value", "mask.min_area_vector_norm",
        "mask.collision_angle_deg", "output_csv", "output_mask_csv"};
    std::map<std::string, std::pair<std::string, int>> kv;   // chiave -> (valore, riga)

    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        const std::string at = where(path, lineNo);
        const size_t eq = line.find('=');
        if (eq == std::string::npos) throw std::runtime_error(at + "riga senza '=': '" + line + "'");
        const std::string key = trim(line.substr(0, eq));
        const std::string value = trim(line.substr(eq + 1));
        if (!known.contains(key)) throw std::runtime_error(at + "chiave non riconosciuta: '" + key + "'");
        if (value.empty()) throw std::runtime_error(at + "valore vuoto per '" + key + "'");
        if (kv.contains(key))
            throw std::runtime_error(at + "chiave '" + key + "' ripetuta (prima alla riga " +
                                     std::to_string(kv[key].second) + ")");
        kv[key] = {value, lineNo};
    }

    for (const auto& k : known)
        if (k != "description" && !kv.contains(k))
            throw std::runtime_error(where(path, lineNo) + "chiave obbligatoria mancante: '" + k + "'");

    ExperimentConfig c;
    c.sourcePath = path;
    const auto& [id, idLine] = kv["experiment_id"];
    // finisce non quotato in ogni riga del CSV: una virgola o un apice lo romperebbero in silenzio
    for (const char ch : id)
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_' && ch != '-' && ch != '.')
            throw std::runtime_error(where(path, idLine) + "experiment_id '" + id +
                                     "' non valido: ammessi solo lettere, cifre, '_', '-', '.'");
    c.experimentId = id;
    if (kv.contains("description")) c.description = kv["description"].first;
    c.latPly = kv["lat_ply"].first;
    c.groundTruthPly = kv["ground_truth_ply"].first;
    c.regionCvCsv = kv["region_cv_csv"].first;
    c.outputCsv = kv["output_csv"].first;
    c.outputMaskCsv = kv["output_mask_csv"].first;

    const auto& [maps, mapsLine] = kv["maps"];
    c.mapsLine = mapsLine;
    if (maps == "all") c.allMaps = true;
    else c.maps = parseMapList(maps, where(path, mapsLine));

    const auto& [sv, svLine] = kv["mask.min_singular_value"];
    c.masks.minSingularValue = parseDouble(sv, where(path, svLine));
    const auto& [av, avLine] = kv["mask.min_area_vector_norm"];
    c.masks.minAreaVectorNorm = parseDouble(av, where(path, avLine));
    const auto& [ca, caLine] = kv["mask.collision_angle_deg"];
    c.masks.collisionAngleDeg = parseDouble(ca, where(path, caLine));
    return c;
}
