//
// Created by gian on 12/19/25.
//

#pragma once
#include <string>
#include <vector>
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>
#include <glm/vec2.hpp>
#include <glm/mat2x2.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_access.hpp>

using namespace glm;

struct Face { int a, b, c; };

// Stima dell'anisotropia di una faccia.
// Base tangente locale, deterministica (dipende solo dall'ordine dei vertici memorizzato):
//   e1 = p_b - p_a,  n = normalize(cross(e1, e2)),  u = normalize(e1),  w = cross(n, u)
// Il tensore (mm^2/ms^2) si stima in questa base; fiber_dir e' riportata in 3D.
struct FaceAnisotropy {
    float CV_L = 0.0f;              // mm/ms
    float CV_T = 0.0f;              // mm/ms
    vec3 fiber_dir = vec3(0);       // unitario, tangente alla faccia (senza verso)
    float cond = 0.0f;              // sigma_max / sigma_min del sistema
    int k_used = 0;                 // righe (mappe) del sistema: tutte le mappe
    bool estimable = false;         // false: faccia lasciata marcata, campi non significativi
};

// Ground truth per faccia. Usata solo per la valutazione: mai passata allo stimatore.
// region_tag e' il tag della zona di velocita' della faccia: la riga del CSV dello scenario
// (experiments/region_cv/) da cui vengono CV_L e CV_T.
struct FaceGroundTruth {
    vec3 fiber_dir = vec3(0);       // unitario (senza verso)
    float CV_L = 0.0f;              // mm/ms
    float CV_T = 0.0f;              // mm/ms
    int region_tag = 0;
};

// Errori di una faccia (stimato vs vero). Se valid == false (faccia non stimabile) err_angle e' NaN.
// Gli errori relativi sulle CV si calcolano in analisi.py da CV_*_est e CV_*_true.
struct FaceError {
    float err_angle = 0.0f;         // gradi, acos(|a.b|): senza verso
    bool valid = false;
};

class Mesh{
    bool imported = false;
    std::vector<Face> faces;
    std::vector<std::vector<float>> vQualities;
    std::vector<vec3> vPos;

    std::vector<std::vector<vec3>> fTangents;

    std::vector<vec3> vNorm;

    std::vector<FaceGroundTruth> fGT;

public:
    Mesh();
    //import export
    void exportAsPly(const std::string &filename) const;
    void exportAllTangentsAsObj(std::string filename) const;
    void importMultiLatPly(const std::string filename);
    // Tiene solo le mappe LAT indicate, nell'ordine dato. Da chiamare prima del calcolo dei gradienti.
    void selectMaps(const std::vector<int>& indices);
    // Ground truth per faccia (elemento face: fx, fy, fz, CV_L, CV_T, elem_tag), copiata cosi' com'e'.
    // Richiede la mesh gia' importata.
    void importPlyGroundTruth(const std::string& filename);

    // Errori per faccia di una stima rispetto alla ground truth. Senza stato: invocabili piu' volte su stime diverse.
    std::vector<FaceError> computeAnisotropyErrors(const std::vector<FaceAnisotropy>& est) const;
    // CSV con una riga per faccia, nessuna aggregazione. experimentId e' la prima colonna, ripetuta su ogni
    // riga: i CSV di esperimenti diversi si concatenano e si raggruppano senza perdere la provenienza.
    void exportErrorsAsCsv(const std::string& filename, const std::string& experimentId,
                           const std::vector<FaceAnisotropy>& est) const;
    // PLY binario con gli errori per faccia letti da un errors.csv (di exportErrorsAsCsv o di estimator/roney.py)
    // sulla mesh importata: per faccia estimable, err_angle (gradi), err_cv_l e err_rapporto (%, con segno,
    // r = CV_L / CV_T), come load() in analysis/analisi.py; NaN sulle facce non stimabili.
    void exportErrorsAsPly(const std::string& errorsCsv, const std::string& filename) const;


    //processing
    // Gradiente LAT per faccia, per ogni mappa.
    void computeGradients();

    void updateVertexNormals();

    // Stima di CV_L, CV_T, fiber_dir della faccia i a partire da fTangents, una riga per mappa.
    // minSingularValue: sigma_min assoluta del sistema (fronti quasi paralleli).
    // Richiede computeGradients() gia' chiamato.
    FaceAnisotropy estimateFace(int i, double minSingularValue) const;
    // estimateFace per ogni faccia (indipendenti fra loro).
    std::vector<FaceAnisotropy> estimateAnisotropy(double minSingularValue) const;


private:
    //processing
    std::vector<vec3> areaVectors() const;
    std::vector<vec3> distributeNormalOnVertex(const std::vector<vec3> &fNorm) const;
    vec3 faceGradientOf(int faceIndex, int mapIndex ) const;
    //import export helpers
    void exportTangentsAsObj(std::string filename, int mapIndex) const;


};
