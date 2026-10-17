//
// Created by gian on 12/19/25.
//

#pragma once
#include <string>
#include <vector>
#include <map>
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>
#include <glm/vec2.hpp>
#include <glm/mat2x2.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_access.hpp>
#include "Face.h"

using namespace glm;



// Stima dell'anisotropia di una faccia (stadio 5).
// Base tangente locale, deterministica (dipende solo dall'ordine dei vertici memorizzato):
//   e1 = p_b - p_a,  n = normalize(cross(e1, e2)),  u = normalize(e1),  w = cross(n, u)
// M11/M12/M22 sono le componenti del tensore (mm^2/ms^2) in questa base {u, w}.
struct FaceAnisotropy {
    vec3 u = vec3(0), w = vec3(0);
    float CV_L = 0.0f;              // mm/ms
    float CV_T = 0.0f;              // mm/ms
    vec3 fiber_dir = vec3(0);       // unitario, tangente alla faccia (senza verso)
    float M11 = 0.0f, M12 = 0.0f, M22 = 0.0f;
    float cond = 0.0f;              // sigma_max / sigma_min del sistema
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
    double minAreaVectorNorm;       // su |cross(e1,e2)| (= 2|a|), mm^2
};

// Velocita' di conduzione vere associate a un tag di regione.
struct RegionCV {
    float CV_L, CV_T;               // mm/ms
};

// Errori di una faccia (stimato vs vero). Se valid == false (faccia non stimabile) i campi sono NaN.
struct FaceError {
    float err_angle = 0.0f;         // gradi, modulo pi: min(|d|, 180 - |d|)
    float err_CV_L_rel = 0.0f;      // (stimato - vero) / vero
    float err_CV_T_rel = 0.0f;
    float err_ratio_rel = 0.0f;     // errore relativo su CV_L / CV_T
    bool valid = false;
};

class Mesh{
    bool imported = false;
    std::vector<Face> faces;
    std::vector<std::vector<float>> vQualities;
    std::vector<vec3> vPos;
    //std::vector<int> vExtEml;
    //std::vector<int> vEml;
    //std::vector<int> vScar;
    //std::vector<int> vOutOfRange;
    // std::vector<int> etoh;
    // std::vector<int> vtoh;
    //std::vector<Half> halfs;
    //bool halfed;


    std::vector<std::vector<vec3>> fTangents;

    std::vector<vec3> vNorm;

    std::vector<FaceAnisotropy> fAniso;

    std::vector<FaceGroundTruth> fGT;

public:
    Mesh();
    //import export
    //bool importPly(const std::string &filename);
    //bool importPlyDualLat(const std::string &filename);
    //bool importMesh(const std::string &filename);
    void exportAsPly(const std::string &filename) const;
    //void exportAsPly2(const std::string &filename) const;
    //void exportTangents1AsObj(const std::string &filename, float scala , int densify = 0)  const;
    //void exportTangents2AsObj(const std::string &filename, float scala , int densify = 0)  const;
    void exportAllTangentsAsObj(std::string filename) const;
    bool importMultiLatPly(const std::string filename);
    int mapCount() const { return vQualities.size(); }
    // Tiene solo le mappe LAT indicate, nell'ordine dato. Da chiamare prima del calcolo dei gradienti.
    void selectMaps(const std::vector<int>& indices);
    // Ground truth per faccia (elemento face: fx, fy, fz, elem_tag). Richiede la mesh gia' importata.
    // Ogni tag presente nel file deve comparire in cvByRegion, altrimenti errore.
    bool importPlyGroundTruth(const std::string& filename, const std::map<int, RegionCV>& cvByRegion);
    // Legge la mappatura region_tag -> (CV_L, CV_T) da CSV (mm/ms; righe vuote, '#' e intestazione saltate).
    // Errore su file mancante, riga malformata, tag duplicato, CV <= 0, nessuna riga dati.
    static std::map<int, RegionCV> importRegionCV(const std::string& filename);
    const std::vector<FaceGroundTruth>& groundTruth() const { return fGT; }

    // Errori per faccia di una stima rispetto alla ground truth. Senza stato: invocabili piu' volte su stime diverse.
    std::vector<FaceError> computeAnisotropyErrors(const std::vector<FaceAnisotropy>& est) const;
    // CSV con una riga per faccia, nessuna aggregazione. experimentId e' la prima colonna, ripetuta su ogni
    // riga: i CSV di esperimenti diversi si concatenano e si raggruppano senza perdere la provenienza.
    bool exportErrorsAsCsv(const std::string& filename, const std::string& experimentId,
                           const std::vector<FaceAnisotropy>& est) const;
    bool exportErrorsAsCsv(const std::string& filename, const std::string& experimentId) const {
        return exportErrorsAsCsv(filename, experimentId, fAniso);
    }


    //processing
    bool updateAllTangetsAsQualitiesGradients();

    bool updateVertexNormals();

    // Stima per faccia (indipendente dalle altre) di CV_L, CV_T, fiber_dir a partire da fTangents.
    // Richiede updateAllTangetsAsQualitiesGradients() gia' chiamato.
    bool updateAnisotropyEstimates(const EstimationMasks& masks);
    const std::vector<FaceAnisotropy>& anisotropy() const { return fAniso; }


private:
    //processing
    //int farthestPointFrom(int x) const;
    //float euristicDiameter() const;
    std::vector<vec3> areaVectors() const;
    std::vector<vec3> distributeNormalOnVertex(const std::vector<vec3> &fNorm) const;
    vec3 faceGradientOf(int faceIndex, int mapIndex ) const;
    //std::vector<vec3> distributeGradientOnVertex(const std::vector<vec3> &fGrad, bool piInvariant = false, bool projectOnTangentPlane = true) const;
    /* bool hasTangentField(const Face &face) const;

    template<typename T, typename Func, typename... Args>
    static bool filterBy(const Face &face, const std::vector<T> &values, Func f,  Args&&... args){
        return f(values[face.a], values[face.b], values[face.c], std::forward<Args>(args)...);
    }
    template<typename T>
    static bool allEqualsTo(T a, T b, T c, T x) {
        return (a == x && b == x && c == x);
    }
    template<typename T>
    static bool equals(T a, T b, T c) {
        return (a == b && b == c);
    }
    template<typename T>
    T valueAt(int fIndex, vec2 bari, const std::vector<T> &vAttr) const {
        Face face = faces[fIndex];
        return vAttr[face.a] * bari.x + vAttr[face.b] * bari.y + vAttr[face.c] * (1 - bari.x - bari.y);
    } */

    bool updateTangetsAsQualityGradient(int x);
    //vec3 posAt( int fIndex, vec2 bari) const { return valueAt<vec3>(fIndex, bari, vPos); }
    //vec3 tangentAt(const std::vector<vec3> vTan, int fIndex, vec2 bari) const { return valueAt<vec3>(fIndex, bari, vTan); }
    //vec3 normalAt( int fIndex, vec2 bari) const { return normalize(valueAt<vec3>(fIndex, bari, vNorm)); }
    //bool populateHalfs();
    //getters
    //static int ftoh(int f) {return f*3;}
    //import export helpers
    //void exportTangentsAsObj(const std::vector<vec3> vTan, const std::string &filename, float scala , int densify)  const;
    void exportTangentsAsObj(std::string filename, int mapIndex) const;


};

