#include <iostream>
#include <filesystem>
#include <stdexcept>
#include "Mesh.h"
#include "ExperimentConfig.h"

// Uso, dalla radice del repository:
//   ./estimator/build/Tesi0 experiments/a.conf [experiments/b.conf ...]
// Un CSV di errori per faccia per esperimento, eseguiti in sequenza, e il PLY se output_ply non e' none.
//   ./estimator/build/Tesi0 --ply MESH.ply ERRORS.csv OUT.ply
// Solo il PLY degli errori di un errors.csv qualsiasi (es. di estimator/roney.py) sulla mesh di MESH.ply.

namespace {
    void runExperiment(const ExperimentConfig& cfg) {
        std::cout << cfg.experimentId << " -> " << cfg.outputCsv << "\n";

        Mesh m;
        m.importMultiLatPly(cfg.latPly);
        if (!cfg.allMaps) m.selectMaps(cfg.maps);
        m.updateVertexNormals();
        m.computeGradients();
        // stima prima di caricare la ground truth: non puo' leggerla
        auto est = m.estimateAnisotropy(cfg.minSingularValue);
        m.importPlyGroundTruth(cfg.groundTruthPly);

        const std::filesystem::path out(cfg.outputCsv);
        if (out.has_parent_path()) std::filesystem::create_directories(out.parent_path());
        m.exportErrorsAsCsv(cfg.outputCsv, cfg.experimentId, est);
        // dal CSV appena scritto: lo stesso percorso dei PLY di roney.py
        if (cfg.outputPly != "none") m.exportErrorsAsPly(cfg.outputCsv, cfg.outputPly);
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "uso: " << argv[0] << " experiments/<esperimento>.conf [altri .conf ...]\n"
                     "     " << argv[0] << " --ply MESH.ply ERRORS.csv OUT.ply\n"
                     "     da lanciare dalla radice del repository\n";
        return 2;
    }
    try {
        if (std::string(argv[1]) == "--ply") {
            if (argc != 5) throw std::runtime_error("uso: --ply MESH.ply ERRORS.csv OUT.ply");
            Mesh m;
            m.importMultiLatPly(argv[2]);
            m.exportErrorsAsPly(argv[3], argv[4]);
            return 0;
        }
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
