#pragma once
#include <string>
#include <vector>
#include "Mesh.h"

// Un esperimento letto da experiments/*.conf (formato in experiments/README.md).
// Percorsi relativi alla radice del repository: l'eseguibile si lancia da la'.
struct ExperimentConfig {
    std::string sourcePath;             // il .conf, per i messaggi d'errore
    std::string experimentId;
    std::string description;            // opzionale
    std::string latPly, groundTruthPly, regionCvCsv, outputCsv;
    bool allMaps = false;               // maps = all
    std::vector<int> maps;              // altrimenti gli indici, nell'ordine dato
    int mapsLine = 0;                   // riga di "maps", per l'errore di indice fuori range
    EstimationMasks masks{};
};

// Errore (runtime_error con file e riga) su: file mancante, riga senza '=', valore vuoto, chiave
// sconosciuta o duplicata, chiave obbligatoria mancante, numero o lista di indici non validi,
// experiment_id con caratteri fuori da [A-Za-z0-9_.-].
ExperimentConfig loadExperimentConfig(const std::string& path);
