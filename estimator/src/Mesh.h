//
// Created by gian on 12/19/25.
//

#pragma once
#include <string>
#include <vector>
#include <map>
#include <array>
#include <cstdint>
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>
#include <glm/vec2.hpp>
#include <glm/mat2x2.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_access.hpp>

using namespace glm;

struct Face { int a, b, c; };

// Stima dell'anisotropia di una faccia (stadio 5).
// Base tangente locale, deterministica (dipende solo dall'ordine dei vertici memorizzato):
//   e1 = p_b - p_a,  n = normalize(cross(e1, e2)),  u = normalize(e1),  w = cross(n, u)
// Il tensore (mm^2/ms^2) si stima in questa base; fiber_dir e' riportata in 3D.
struct FaceAnisotropy {
    float CV_L = 0.0f;              // mm/ms
    float CV_T = 0.0f;              // mm/ms
    vec3 fiber_dir = vec3(0);       // unitario, tangente alla faccia (senza verso)
    float cond = 0.0f;              // sigma_max / sigma_min del sistema
    int k_used = 0;                 // righe (mappe) non escluse dalle maschere, entrate nel sistema
    bool estimable = false;         // false: faccia lasciata marcata, campi non significativi
};

// Ground truth per faccia (stadio 6). Usata solo per la valutazione: mai passata allo stimatore.
// region_tag e' solo un'etichetta di provenienza dei valori di CV.
struct FaceGroundTruth {
    vec3 fiber_dir = vec3(0);       // unitario, tangente alla faccia (senza verso)
    float CV_L = 0.0f;              // mm/ms
    float CV_T = 0.0f;              // mm/ms
    int region_tag = 0;
};

// Soglie delle maschere della stima (stadio 5). Vengono dal file di esperimento, non hanno default nel codice.
struct EstimationMasks {
    double minSingularValue;        // sigma_min assoluta di R (fronti quasi paralleli)
    double collisionAngleDeg;       // coppia esclusa se l'angolo max col gradiente di un vicino > soglia, gradi
};

// Velocita' di conduzione vere associate a un tag di regione.
struct RegionCV {
    float CV_L, CV_T;               // mm/ms
};

// Errori di una faccia (stimato vs vero). Se valid == false (faccia non stimabile) err_angle e' NaN.
// Gli errori relativi sulle CV si calcolano in analisi.py da CV_*_est e CV_*_true.
struct FaceError {
    float err_angle = 0.0f;         // gradi, modulo pi: min(|d|, 180 - |d|)
    bool valid = false;
};

class Mesh{
    bool imported = false;
    std::vector<Face> faces;
    std::vector<std::vector<float>> vQualities;
    std::vector<std::string> qualityNames;          // nome del sito di ogni mappa, allineato a vQualities
    std::vector<vec3> vPos;
    // allineato a faces; slot s = vicino oltre lo spigolo v[s]-v[(s+1)%3], v = {a,b,c}; -1 se di bordo
    std::vector<std::array<int,3>> faceNeighbors;

    std::vector<std::vector<vec3>> fTangents;

    // fMask[k][i]: 1 se la coppia (faccia i, mappa k) e' esclusa per collisione, 0 se tenuta; map-major come vQualities
    std::vector<std::vector<uint8_t>> fMask;

    std::vector<vec3> vNorm;

    std::vector<FaceAnisotropy> fAniso;

    std::vector<FaceGroundTruth> fGT;

public:
    Mesh();
    //import export
    void exportAsPly(const std::string &filename) const;
    void exportAllTangentsAsObj(std::string filename) const;
    void importMultiLatPly(const std::string filename);
    // Tiene solo le mappe LAT indicate, nell'ordine dato. Da chiamare prima del calcolo dei gradienti.
    void selectMaps(const std::vector<int>& indices);
    // Ground truth per faccia (elemento face: fx, fy, fz, elem_tag). Richiede la mesh gia' importata.
    // Ogni tag presente nel file deve comparire in cvByRegion, altrimenti errore.
    void importPlyGroundTruth(const std::string& filename, const std::map<int, RegionCV>& cvByRegion);
    // Legge la mappatura region_tag -> (CV_L, CV_T) da CSV (mm/ms; righe vuote, '#' e intestazione saltate).
    static std::map<int, RegionCV> importRegionCV(const std::string& filename);

    // Errori per faccia di una stima rispetto alla ground truth. Senza stato: invocabili piu' volte su stime diverse.
    std::vector<FaceError> computeAnisotropyErrors(const std::vector<FaceAnisotropy>& est) const;
    // CSV con una riga per faccia, nessuna aggregazione. experimentId e' la prima colonna, ripetuta su ogni
    // riga: i CSV di esperimenti diversi si concatenano e si raggruppano senza perdere la provenienza.
    void exportErrorsAsCsv(const std::string& filename, const std::string& experimentId,
                           const std::vector<FaceAnisotropy>& est) const;
    void exportErrorsAsCsv(const std::string& filename, const std::string& experimentId) const {
        exportErrorsAsCsv(filename, experimentId, fAniso);
    }


    //processing
    // Gradiente LAT per faccia, per ogni mappa.
    void computeGradients();

    void updateVertexNormals();

    // Maschere per coppia (faccia, mappa), dai soli LAT, gradienti grezzi e geometria. Calcolate una volta.
    // Richiede computeGradients() gia' chiamato.
    void updateMasks(const EstimationMasks& masks);

    // Stima per faccia (indipendente dalle altre) di CV_L, CV_T, fiber_dir a partire da fTangents,
    // usando solo le righe con fMask == 0. Richiede updateMasks() gia' chiamato.
    void updateAnisotropyEstimates(const EstimationMasks& masks);


private:
    //processing
    std::vector<vec3> areaVectors() const;
    std::vector<vec3> distributeNormalOnVertex(const std::vector<vec3> &fNorm) const;
    vec3 faceGradientOf(int faceIndex, int mapIndex ) const;
    //import export helpers
    void exportTangentsAsObj(std::string filename, int mapIndex) const;


};
