#pragma once
#include <string>
#include <vector>

// Un esperimento letto da experiments/*.conf (formato in experiments/README.md).
// Percorsi relativi alla radice del repository: l'eseguibile si lancia da la'.
struct ExperimentConfig {
    std::string sourcePath;             // il .conf, per i messaggi d'errore
    std::string experimentId;
    std::string latPly, groundTruthPly, outputCsv;
    bool allMaps = false;               // maps = all
    std::vector<int> maps;              // altrimenti gli indici, nell'ordine dato
    double minSingularValue = 0.0;      // sigma_min assoluta del sistema (fronti quasi paralleli)
    double collisionAngleDeg = 0.0;     // coppia esclusa se l'angolo max col gradiente di un vicino > soglia, gradi
};

// Errore (runtime_error) su: file mancante, chiave sconosciuta (con file e riga), chiave obbligatoria mancante.
ExperimentConfig loadExperimentConfig(const std::string& path);
