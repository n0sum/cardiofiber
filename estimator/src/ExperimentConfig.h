#pragma once
#include <string>
#include <vector>

// Un esperimento letto da experiments/*.conf (formato in experiments/README.md).
// Percorsi relativi alla radice del repository: l'eseguibile si lancia da la'.
struct ExperimentConfig {
    std::string experimentId;
    std::string latPly, groundTruthPly, outputCsv;
    std::string outputPly;              // base dei PLY degli errori (BASE_angolo.ply, ...), oppure "none"
    bool allMaps = false;               // maps = all
    std::vector<int> maps;              // altrimenti gli indici, nell'ordine dato
    double minSingularValue = 0.0;      // sigma_min assoluta del sistema (fronti quasi paralleli)
};

// Una chiave mancante ferma tutto (std::out_of_range da map::at).
ExperimentConfig loadExperimentConfig(const std::string& path);
