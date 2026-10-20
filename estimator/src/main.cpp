#include <iostream>
#include <filesystem>
#include <stdexcept>
#include "Mesh.h"
#include "ExperimentConfig.h"

// Uso, dalla radice del repository:
//   ./estimator/build/Tesi0 experiments/a.conf [experiments/b.conf ...]
// Un CSV di errori per faccia per esperimento, eseguiti in sequenza.

namespace {
    void runExperiment(const ExperimentConfig& cfg) {
        std::cout << "== esperimento " << cfg.experimentId << "  (" << cfg.sourcePath << ")\n";
        std::cout << "   maps = ";
        if (cfg.allMaps) std::cout << "all";
        else for (size_t i = 0; i < cfg.maps.size(); ++i) std::cout << (i ? "," : "") << cfg.maps[i];
        std::cout << "  mask.min_singular_value = " << cfg.masks.minSingularValue
                  << "  mask.collision_angle_deg = " << cfg.masks.collisionAngleDeg << "\n";

        Mesh m;
        m.importMultiLatPly(cfg.latPly);
        if (!cfg.allMaps) m.selectMaps(cfg.maps);
        m.updateVertexNormals();
        m.computeGradients();
        // maschere e stima prima di caricare la ground truth: non possono leggerla
        m.updateMasks(cfg.masks);
        m.updateAnisotropyEstimates(cfg.masks);
        m.importPlyGroundTruth(cfg.groundTruthPly, Mesh::importRegionCV(cfg.regionCvCsv));

        for (const std::string& p : {cfg.outputCsv, cfg.outputMaskCsv}) {
            const std::filesystem::path out(p);
            if (out.has_parent_path()) std::filesystem::create_directories(out.parent_path());
        }
        m.exportErrorsAsCsv(cfg.outputCsv, cfg.experimentId);
        std::cout << "   scritto " << cfg.outputCsv << "\n";
        m.exportMaskPairsAsCsv(cfg.outputMaskCsv, cfg.experimentId);
        std::cout << "   scritto " << cfg.outputMaskCsv << "\n";
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "uso: " << argv[0] << " experiments/<esperimento>.conf [altri .conf ...]\n"
                     "     da lanciare dalla radice del repository\n";
        return 2;
    }
    try {
        // prima tutti i file: un errore nell'ultimo non deve emergere dopo una run lunga
        std::vector<ExperimentConfig> configs;
        for (int i = 1; i < argc; ++i) configs.push_back(loadExperimentConfig(argv[i]));
        for (const auto& cfg : configs) runExperiment(cfg);
    } catch (const std::exception& e) {
        std::cerr << "errore: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
