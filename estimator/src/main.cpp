#include <iostream>
#include <filesystem>
#include "Mesh.h"
#include "ExperimentConfig.h"

// Uso, dalla radice del repository:
//   ./estimator/build/Tesi0 experiments/a.conf [experiments/b.conf ...]
// Un CSV di errori per faccia per esperimento, in sequenza.
//   ./estimator/build/Tesi0 --ply MESH.ply ERRORS.csv BASE
// I PLY degli errori di un errors.csv (di Tesi0 o di estimator/roney.py) sulla mesh di MESH.ply:
// BASE_angolo.ply, BASE_cv_l.ply, BASE_rapporto.ply.

namespace {
    void runExperiment(const ExperimentConfig& cfg) {
        std::cout << cfg.experimentId << " -> " << cfg.outputCsv << "\n";

        Mesh m;
        m.importMultiLatPly(cfg.latPly);
        if (!cfg.allMaps) m.selectMaps(cfg.maps);
        m.computeGradients();
        // stima prima di caricare la ground truth: non puo' leggerla
        auto est = m.estimateAnisotropy(cfg.minSingularValue);
        m.importPlyGroundTruth(cfg.groundTruthPly);

        std::filesystem::create_directories(std::filesystem::path(cfg.outputCsv).parent_path());
        m.exportErrorsAsCsv(cfg.outputCsv, cfg.experimentId, est);
    }
}

int main(int argc, char** argv) {
    const std::vector<std::string> args(argv + 1, argv + argc);
    if (!args.empty() && args[0] == "--ply") {
        Mesh m;
        m.importMultiLatPly(args.at(1));
        m.exportErrorsAsPly(args.at(2), args.at(3));
        return 0;
    }
    for (const auto& conf : args) runExperiment(loadExperimentConfig(conf));
    return 0;
}
