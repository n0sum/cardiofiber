//
// Created by gian on 12/19/25.
//

#include <happly.h>
#include "Mesh.h"
#include "Face.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <functional>
#include <cctype>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <iomanip>
#include <cstdlib>
#include <stdexcept>

//PROSSIMI TODO
// sistemare normali


//helpers

std::vector<vec3> normalizeVectors(const std::vector<vec3> &vectors) {
    std::vector<vec3> result(vectors.size());
    for (int i = 0; i<vectors.size(); ++i) {
        result[i] = normalize(vectors[i]);
    }
    return result;
}


// Soglie strutturali della stima di anisotropia (stadio 5). Le soglie delle maschere
// (EstimationMasks) vengono dal file di esperimento.
namespace {
    constexpr int    kMinGradients       = 3;       // sotto: sistema sottodeterminato
    constexpr double kMinEigenvalue      = 0.0;     // lambda_2 deve essere > di questo (CV_T^2, mm^2/ms^2)

    struct LsqResult {
        std::array<double,3> m{};
        double sigmaMax = 0.0, sigmaMin = 0.0;
    };

    // Minimi quadrati R m = b con R (K x 3) via SVD one-sided Jacobi: ruota le colonne di R
    // fino a renderle ortogonali, cosi' R V = U S con S = norme delle colonne.
    // m = V S^-1 U^T b = sum_j v_j (u_j . b) / sigma_j. Se sigma_min <= minSingularValue ritorna false.
    bool svdLeastSquares(std::vector<std::array<double,3>> R, const std::vector<double> &b, double minSingularValue, LsqResult &out) {
        const int K = R.size();
        double V[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
        for (int sweep = 0; sweep < 60; ++sweep) {
            bool rotated = false;
            for (int p = 0; p < 2; ++p) for (int q = p + 1; q < 3; ++q) {
                double alpha = 0, beta = 0, gamma = 0;
                for (int i = 0; i < K; ++i) {
                    alpha += R[i][p]*R[i][p];
                    beta  += R[i][q]*R[i][q];
                    gamma += R[i][p]*R[i][q];
                }
                if (std::abs(gamma) <= 1e-15 * std::sqrt(alpha*beta) || gamma == 0.0) continue;
                rotated = true;
                const double zeta = (beta - alpha) / (2.0*gamma);
                const double t = (zeta >= 0 ? 1.0 : -1.0) / (std::abs(zeta) + std::sqrt(1.0 + zeta*zeta));
                const double c = 1.0 / std::sqrt(1.0 + t*t), s = c*t;
                for (int i = 0; i < K; ++i) {
                    const double rp = R[i][p], rq = R[i][q];
                    R[i][p] = c*rp - s*rq;
                    R[i][q] = s*rp + c*rq;
                }
                for (int i = 0; i < 3; ++i) {
                    const double vp = V[i][p], vq = V[i][q];
                    V[i][p] = c*vp - s*vq;
                    V[i][q] = s*vp + c*vq;
                }
            }
            if (!rotated) break;
        }
        double sigma[3];
        for (int j = 0; j < 3; ++j) {
            double n2 = 0;
            for (int i = 0; i < K; ++i) n2 += R[i][j]*R[i][j];
            sigma[j] = std::sqrt(n2);
        }
        out.sigmaMax = std::max({sigma[0], sigma[1], sigma[2]});
        out.sigmaMin = std::min({sigma[0], sigma[1], sigma[2]});
        if (!(out.sigmaMin > minSingularValue)) return false;
        out.m = {0,0,0};
        for (int j = 0; j < 3; ++j) {
            double ub = 0;                                   // u_j . b, con u_j = colonna_j / sigma_j
            for (int i = 0; i < K; ++i) ub += R[i][j]*b[i];
            ub /= sigma[j];
            for (int r = 0; r < 3; ++r) out.m[r] += V[r][j]*ub/sigma[j];
        }
        return true;
    }
}

// metodi processing

// int Mesh::farthestPointFrom(int x) const {
//     int maxInd = x;
//     float max = 0.0;
//     for (int i = 0; i < vPos.size(); ++i){
//         float dist = distance(vPos[x], vPos[i]);
//         if (dist > max) {
//             maxInd = i;
//             max = dist;
//         }
//     }
//     return maxInd;
// }
// float Mesh::euristicDiameter() const {
//     int x = 0;
//     for (int i = 0; i < 3; ++i) {
//         x = farthestPointFrom(x);
//     }
//     return distance(vPos[x],vPos[farthestPointFrom(x)]);
// }

std::vector<vec3> Mesh::areaVectors() const {
    std::vector<vec3>  result(faces.size() );
   for (int i = 0; Face face : faces) {
        vec3 e1 = vPos[face.b] - vPos[face.a];
        vec3 e2 = vPos[face.c] - vPos[face.a];
        result[i++] = cross(e1,e2);
    }
    return result;
}

std::vector<vec3> Mesh::distributeNormalOnVertex(const std::vector<vec3> &fNorm) const {
    std::vector<vec3> result(vPos.size());
    for (int i = 0; i < faces.size(); ++i) {
        result[faces[i].a] += fNorm[i];
        result[faces[i].b] += fNorm[i];
        result[faces[i].c] += fNorm[i];
    }
    return result;
}

vec3 Mesh::faceGradientOf(int faceIndex, int mapIndex) const {
    Face face = faces[faceIndex];
    vec3 e1 = vPos[face.b] - vPos[face.a];
    vec3 e2 = vPos[face.c] - vPos[face.a];
    mat2 A(
        dot(e1,e1),dot(e1,e2),
        dot(e1,e2), dot(e2,e2)
        );
    vec2 k(vQualities[mapIndex][face.b]-vQualities[mapIndex][face.a], vQualities[mapIndex][face.c]-vQualities[mapIndex][face.a]);
    if (std::abs(determinant(A)) < 1e-6f) {
        return vec3(0,0,0);
    }
    mat2 inv = inverse(A);
    vec2 soluzione = inv*k;
    return (soluzione.x * e1 + soluzione.y * e2);
}

// std::vector<vec3> Mesh::distributeGradientOnVertex(const std::vector<vec3> &fGrad, bool piInvariant, bool projectOnTangentPlane) const {
//     const int vNumber = vPos.size();
//     std::vector<vec3> result(vNumber);
//     std::vector<int> count(vNumber);
//
//     if (piInvariant) {
//         for (int i = 0; i < faces.size(); ++i) if (hasTangentField(faces[i])){
//             const Face &f = faces[i];
//             const int verts[3] = {f.a, f.b, f.c};
//             for (int vert : verts) {
//                 vec3 grad = fGrad[i];
//                 if (dot(result[vert],grad) < 0) {
//                     grad = -grad;
//                 }
//                 result[vert] += grad;
//                 ++count[vert];
//             }
//         }
//         for (int i = 0; i < vNumber; ++i) {
//             if (count[i] > 0) result[i] /= count[i];
//         }
//     } else {
//         for (int i = 0; i < faces.size(); ++i) if (hasTangentField(faces[i])){
//             const Face &f = faces[i];
//             result[f.a] += fGrad[i];
//             result[f.b] += fGrad[i];
//             result[f.c] += fGrad[i];
//             ++count[f.a];
//             ++count[f.b];
//             ++count[f.c];
//         }
//         for (int i = 0; i < vNumber; ++i) {
//             if (count[i] > 0) result[i] /= static_cast<float>(count[i]);
//         }
//     }
//
//     if (projectOnTangentPlane) {
//         for (int i = 0; i < vNumber; ++i) {
//             const vec3 &n = vNorm[i];
//             result[i] -= dot(result[i], n) * n;
//         }
//     }
//
//     return result;
// }

bool Mesh::updateVertexNormals() {
    vNorm = normalizeVectors(distributeNormalOnVertex(areaVectors()));
    return true;
}

bool Mesh::updateTangetsAsQualityGradient(int map) {
    if (fTangents.size() != vQualities.size() ){
        fTangents.resize(vQualities.size());
    }
    fTangents[map].resize(faces.size());
    for (int i = 0; i < faces.size(); ++i){
        fTangents[map][i] = faceGradientOf(i,map);
    }
    return true;
}

bool Mesh::updateAnisotropyEstimates(const EstimationMasks& masks) {
    fAniso.assign(faces.size(), FaceAnisotropy());
    const int K = fTangents.size();
    if (K < kMinGradients) return true;     // tutte le facce restano estimable = false

    for (int i = 0; i < faces.size(); ++i) {
        FaceAnisotropy &fa = fAniso[i];
        const Face &face = faces[i];

        // base locale, calcolata una sola volta e usata sia per proiettare sia per ricostruire fiber_dir
        const vec3 e1 = vPos[face.b] - vPos[face.a];
        const vec3 e2 = vPos[face.c] - vPos[face.a];
        const vec3 a2 = cross(e1, e2);
        const float a2len = length(a2);
        const float e1len = length(e1);
        if (!(a2len > masks.minAreaVectorNorm) || !(e1len > 0.0f)) continue;
        const vec3 n = a2 / a2len;
        const vec3 u = e1 / e1len;
        const vec3 w = cross(n, u);
        fa.u = u;
        fa.w = w;

        std::vector<std::array<double,3>> R(K);
        const std::vector<double> b(K, 1.0);
        for (int k = 0; k < K; ++k) {
            const vec3 g = fTangents[k][i];
            const double g1 = dot(g, u), g2 = dot(g, w);
            R[k] = {g1*g1, 2.0*g1*g2, g2*g2};
        }

        LsqResult lsq;
        if (!svdLeastSquares(std::move(R), b, masks.minSingularValue, lsq)) continue;
        fa.cond = static_cast<float>(lsq.sigmaMax / lsq.sigmaMin);

        const double M11 = lsq.m[0], M12 = lsq.m[1], M22 = lsq.m[2];
        const double somma = M11 + M22;
        const double differenza = std::sqrt((M11 - M22)*(M11 - M22) + (2.0*M12)*(2.0*M12));
        const double lambda1 = (somma + differenza) / 2.0;   // CV_L^2
        const double lambda2 = (somma - differenza) / 2.0;   // CV_T^2
        if (!(lambda2 > kMinEigenvalue)) continue;

        const double theta = 0.5 * std::atan2(2.0*M12, M11 - M22);
        fa.M11 = M11; fa.M12 = M12; fa.M22 = M22;
        fa.CV_L = static_cast<float>(std::sqrt(lambda1));
        fa.CV_T = static_cast<float>(std::sqrt(lambda2));
        fa.fiber_dir = static_cast<float>(std::cos(theta))*u + static_cast<float>(std::sin(theta))*w;
        fa.estimable = true;
    }
    return true;
}

bool Mesh::updateAllTangetsAsQualitiesGradients(){
    for(int i = 0; i < vQualities.size(); ++i){
        updateTangetsAsQualityGradient(i);
    }
    return true;
};


// bool Mesh::hasTangentField(const Face &face) const {
//     bool notExtEml = filterBy(face, vExtEml,    allEqualsTo<int>, 0);
//     bool notEml    = filterBy(face, vEml,        allEqualsTo<int>, 0);
//     bool notScar   = filterBy(face, vScar,       allEqualsTo<int>, 0);
//     bool notOor    = filterBy(face, vOutOfRange, allEqualsTo<int>, 0);
//     return (notEml && notExtEml && notScar && notOor);
// }

// bool Mesh::populateHalfs() {
//     struct PairHash {
//         size_t operator()(const std::pair<int, int>& p) const {
//             return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
//         }
//     };
//
//     std::unordered_map<std::pair<int, int>, int, PairHash> mappa;
//
//     for (int i = 0; i < faces.size(); ++i) {
//         int a = faces[i].a;
//         int b = faces[i].b;
//         int c = faces[i].c;
//         halfs.emplace_back(a, i*3 + 1, i);
//         if (mappa.contains({b, a})) {
//             Half half = halfs[mappa[{b,a}]];
//             half.opp = i*3;
//             halfs[i*3].opp = mappa[{b,a}];
//         }
//         halfs.emplace_back(faces[i].b, i*3 + 2, i);
//         halfs.emplace_back(faces[i].c, i*3, i);
//     }
// }

// import export

Mesh::Mesh() = default;

bool Mesh::importMultiLatPly(std::string filename) {
    if (imported) throw std::runtime_error("already imported a mesh in this object");
    imported = true;

    happly::PLYData plyIn(filename);
    auto& v = plyIn.getElement("vertex");

    std::vector<float> x = v.getProperty<float>("x");
    std::vector<float> y = v.getProperty<float>("y");
    std::vector<float> z = v.getProperty<float>("z");

    vPos.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i) vPos[i] = vec3(x[i], y[i], z[i]);

    // vQualities[mapIndex][vertexIndex]
    for (int k = 0; v.hasProperty("lat_" + std::to_string(k)); ++k)
        vQualities.push_back(v.getProperty<float>("lat_" + std::to_string(k)));

    const auto& prov_faces = plyIn.getElement("face").getListProperty<int>("vertex_indices");
    faces.reserve(prov_faces.size());
    for (const auto& f : prov_faces) faces.emplace_back(f[0], f[1], f[2]);

    std::cout << "importPlyLat: " << vPos.size() << " vertici, "
              << faces.size() << " facce, " << vQualities.size() << " mappe\n";
    return true;
}

void Mesh::selectMaps(const std::vector<int>& indices) {
    if (!fTangents.empty()) throw std::runtime_error("selectMaps: gradienti gia' calcolati");
    std::vector<std::vector<float>> selected;
    selected.reserve(indices.size());
    for (int k : indices) {
        if (k < 0 || k >= static_cast<int>(vQualities.size()))
            throw std::runtime_error("selectMaps: mappa " + std::to_string(k) + " fuori range (" +
                                     std::to_string(vQualities.size()) + " mappe)");
        selected.push_back(vQualities[k]);
    }
    vQualities = std::move(selected);
}

bool Mesh::importPlyGroundTruth(const std::string& filename, const std::map<int, RegionCV>& cvByRegion) {
    if (faces.empty()) throw std::runtime_error("importPlyGroundTruth: importare prima la mesh");

    happly::PLYData plyIn(filename);
    auto& fe = plyIn.getElement("face");
    const std::vector<float> fx = fe.getProperty<float>("fx");
    const std::vector<float> fy = fe.getProperty<float>("fy");
    const std::vector<float> fz = fe.getProperty<float>("fz");
    const std::vector<int> tag = fe.getProperty<int>("elem_tag");
    if (fx.size() != faces.size())
        throw std::runtime_error("importPlyGroundTruth: " + std::to_string(fx.size()) + " facce nel file, "
                                 + std::to_string(faces.size()) + " nella mesh");

    std::vector<FaceGroundTruth> gt(faces.size());
    std::map<int, int> tagCount;
    double maxDevDeg = 0.0;
    for (int i = 0; i < faces.size(); ++i) {
        const Face &face = faces[i];
        const vec3 n = normalize(cross(vPos[face.b] - vPos[face.a], vPos[face.c] - vPos[face.a]));
        const vec3 d(fx[i], fy[i], fz[i]);
        const float dlen = length(d);
        const vec3 t = d - dot(d, n) * n;                    // proiezione sul piano della faccia
        const float tlen = length(t);
        if (!(dlen > 0.0f) || !(tlen > 1e-6f * dlen))
            throw std::runtime_error("importPlyGroundTruth: fibra nulla o normale alla faccia " + std::to_string(i));
        maxDevDeg = std::max(maxDevDeg, std::asin(std::min(1.0f, std::abs(dot(d, n)) / dlen)) * 180.0 / M_PI);

        const auto it = cvByRegion.find(tag[i]);
        if (it == cvByRegion.end())
            throw std::runtime_error("importPlyGroundTruth: tag " + std::to_string(tag[i]) + " (faccia "
                                     + std::to_string(i) + ") assente dalla mappatura region_tag -> CV");
        gt[i].fiber_dir = t / tlen;
        gt[i].CV_L = it->second.CV_L;
        gt[i].CV_T = it->second.CV_T;
        gt[i].region_tag = tag[i];
        ++tagCount[tag[i]];
    }
    fGT = std::move(gt);

    std::cout << "importPlyGroundTruth: " << fGT.size() << " facce, deviazione max dal piano "
              << maxDevDeg << " gradi; tag:";
    for (const auto &[t, n] : tagCount) std::cout << " " << t << "=" << n;
    std::cout << "\n";
    return true;
}

std::map<int, RegionCV> Mesh::importRegionCV(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("importRegionCV: impossibile aprire " + filename);

    auto trim = [](std::string s) {
        const char *ws = " \t\r\n";
        const size_t b = s.find_first_not_of(ws);
        if (b == std::string::npos) return std::string();
        return s.substr(b, s.find_last_not_of(ws) - b + 1);
    };
    // intero/float che deve consumare tutto il campo
    auto parseInt = [](const std::string &s, long &v) {
        char *end = nullptr;
        v = std::strtol(s.c_str(), &end, 10);
        return !s.empty() && *end == '\0';
    };
    auto parseFloat = [](const std::string &s, float &v) {
        char *end = nullptr;
        v = std::strtof(s.c_str(), &end);
        return !s.empty() && *end == '\0';
    };

    std::map<int, RegionCV> result;
    std::string line;
    for (int lineNo = 1; std::getline(in, line); ++lineNo) {
        line = trim(line);
        if (line.empty() || line[0] == '#' || line == "region_tag,CV_L,CV_T") continue;

        const std::string where = filename + ":" + std::to_string(lineNo) + ": ";
        std::vector<std::string> f;
        std::stringstream ss(line);
        for (std::string tok; std::getline(ss, tok, ','); ) f.push_back(trim(tok));
        if (!line.empty() && line.back() == ',') f.push_back("");
        long tag; float cvL, cvT;
        if (f.size() != 3 || !parseInt(f[0], tag) || !parseFloat(f[1], cvL) || !parseFloat(f[2], cvT)
            || tag < std::numeric_limits<int>::min() || tag > std::numeric_limits<int>::max())
            throw std::runtime_error("importRegionCV: " + where + "riga malformata: '" + line + "'");
        if (!std::isfinite(cvL) || !std::isfinite(cvT) || cvL <= 0.0f || cvT <= 0.0f)
            throw std::runtime_error("importRegionCV: " + where + "CV non positiva o non finita");
        if (!result.emplace(static_cast<int>(tag), RegionCV{cvL, cvT}).second)
            throw std::runtime_error("importRegionCV: " + where + "region_tag " + std::to_string(tag) + " duplicato");
    }
    if (result.empty()) throw std::runtime_error("importRegionCV: nessuna riga dati in " + filename);
    return result;
}

std::vector<FaceError> Mesh::computeAnisotropyErrors(const std::vector<FaceAnisotropy>& est) const {
    if (est.size() != faces.size() || fGT.size() != faces.size())
        throw std::runtime_error("computeAnisotropyErrors: stima o ground truth non allineate alle facce");
    const float nan = std::numeric_limits<float>::quiet_NaN();
    std::vector<FaceError> err(faces.size());
    for (int i = 0; i < faces.size(); ++i) {
        FaceError &e = err[i];
        e.err_angle = e.err_CV_L_rel = e.err_CV_T_rel = e.err_ratio_rel = nan;
        if (!est[i].estimable) continue;                      // lasciata marcata, non esclusa ne' riempita

        const Face &face = faces[i];
        const vec3 n = normalize(cross(vPos[face.b] - vPos[face.a], vPos[face.c] - vPos[face.a]));
        const vec3 a = normalize(est[i].fiber_dir - dot(est[i].fiber_dir, n) * n);
        const vec3 b = normalize(fGT[i].fiber_dir - dot(fGT[i].fiber_dir, n) * n);
        // angolo nel piano, poi modulo pi: nessuna scelta di orientamento
        const double delta = std::abs(std::atan2(double(dot(n, cross(a, b))), double(dot(a, b))));
        e.err_angle = static_cast<float>(std::min(delta, M_PI - delta) * 180.0 / M_PI);

        const FaceGroundTruth &g = fGT[i];
        e.err_CV_L_rel = (est[i].CV_L - g.CV_L) / g.CV_L;
        e.err_CV_T_rel = (est[i].CV_T - g.CV_T) / g.CV_T;
        const double ratioTrue = double(g.CV_L) / g.CV_T;
        e.err_ratio_rel = static_cast<float>((double(est[i].CV_L) / est[i].CV_T - ratioTrue) / ratioTrue);
        e.valid = true;
    }
    return err;
}

bool Mesh::exportErrorsAsCsv(const std::string& filename, const std::string& experimentId,
                             const std::vector<FaceAnisotropy>& est) const {
    const std::vector<FaceError> err = computeAnisotropyErrors(est);
    const std::vector<vec3> av = areaVectors();                // cross(e1,e2): area = |cross| / 2
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("exportErrorsAsCsv: impossibile aprire " + filename);
    out << std::scientific << std::setprecision(9);
    out << "experiment_id,face_id,area,region_tag,estimable,cond,CV_L_est,CV_T_est,CV_L_true,CV_T_true,"
           "err_angle,err_CV_L_rel,err_CV_T_rel,err_ratio_rel\n";
    for (int i = 0; i < faces.size(); ++i) {
        const FaceAnisotropy &a = est[i];
        const FaceGroundTruth &g = fGT[i];
        const FaceError &e = err[i];
        out << experimentId << ',' << i << ',' << length(av[i]) / 2.0f << ',' << g.region_tag << ',' << (a.estimable ? 1 : 0) << ',';
        if (a.cond > 0.0f) out << a.cond;                       // vuota se non calcolata
        out << ',';
        if (a.estimable) out << a.CV_L << ',' << a.CV_T << ','; else out << ",,";
        out << g.CV_L << ',' << g.CV_T << ',';
        if (e.valid) out << e.err_angle << ',' << e.err_CV_L_rel << ',' << e.err_CV_T_rel << ',' << e.err_ratio_rel;
        else out << ",,,";
        out << '\n';
    }
    return static_cast<bool>(out);
}

// bool Mesh::importMesh(const std::string& filename) {
//
//     if (imported) throw std::runtime_error("already imported a mesh in this object");
//     imported = true;
//
//     auto split = [](const std::string& str) {
//         std::vector<std::string> tokens;
//         std::istringstream iss(str);
//         std::string token;
//
//         while (iss >> token) {
//             tokens.push_back(token);
//         }
//         return tokens;
//     };
//
//
//     // Lettura file
//     std::ifstream file(filename);
//     if (!file.is_open()) {
//         throw std::runtime_error("Impossibile aprire il file: " + filename);
//     }
//
//     std::string line;
//     int count = 0;
//     int mod = 0;
//     while (std::getline(file, line)) {
//         std::vector<std::string> splitted = split(line);
//         if (splitted.size() == 0) {
//             mod = 0;
//             continue;
//         }
//         if (splitted[0]=="0") {
//             count++;
//             mod = count;
//         }
//         switch (mod) {
//             case 0: break;
//             case 1: {
//                 vPos.emplace_back(stof(splitted[2]),stof(splitted[3]),stof(splitted[4]));
//                 break;
//             }
//             case 2: {
//                 if (stoi(splitted[8])!= -1000000) {
//                     faces.emplace_back(stoi(splitted[2]),stoi(splitted[3]),stoi(splitted[4]));
//                 }
//                 break;
//             }
//             case 3: {
//                 vQuality.push_back(stof(splitted[4]));
//                 break;
//             }
//             case 4: {
//                 vEml.push_back(stoi(splitted[2]));
//                 vExtEml.push_back(stoi(splitted[3]));
//                 vScar.push_back(stoi(splitted[4]));
//
//                 break;
//             }
//         }
//
//     }
//     return true;
// }

// bool Mesh::importPly(const std::string &filename){
//
//     if (imported) throw std::runtime_error("already imported a mesh in this object");
//     imported = true;
//
//     happly::PLYData plyIn(filename);
//
//     vQuality = plyIn.getElement("vertex").getProperty<float>("quality");
//     std::vector<float> x = plyIn.getElement("vertex").getProperty<float>("x");
//     std::vector<float> y = plyIn.getElement("vertex").getProperty<float>("y");
//     std::vector<float> z = plyIn.getElement("vertex").getProperty<float>("z");
//     const std::vector<std::vector<int>> prov_faces = plyIn.getElement("face").getListProperty<int>("vertex_indices");
//     std::cout << x.size() << "\n";
//     vPos.resize(x.size());
//     for (int j = 0; j < x.size(); ++j) {
//         vPos[j] = vec3(x[j],y[j],z[j]);
//     }
//
//     faces.reserve(prov_faces.size());
//     for (const auto & prov_face : prov_faces) {
//         faces.emplace_back(Face(prov_face[0],prov_face[1],prov_face[2]));
//     }
//
//
//     return true;
// }

// bool Mesh::importPlyDualLat(const std::string &filename) {
//     if (imported) throw std::runtime_error("already imported a mesh in this object");
//     imported = true;
//
//     happly::PLYData plyIn(filename);
//
//     vQuality  = plyIn.getElement("vertex").getProperty<float>("quality");
//     vQuality2 = plyIn.getElement("vertex").getProperty<float>("quality_transferred");
//
//     std::vector<float> x  = plyIn.getElement("vertex").getProperty<float>("x");
//     std::vector<float> y  = plyIn.getElement("vertex").getProperty<float>("y");
//     std::vector<float> z  = plyIn.getElement("vertex").getProperty<float>("z");
//     std::vector<float> nx = plyIn.getElement("vertex").getProperty<float>("nx");
//     std::vector<float> ny = plyIn.getElement("vertex").getProperty<float>("ny");
//     std::vector<float> nz = plyIn.getElement("vertex").getProperty<float>("nz");
//
//     int n = static_cast<int>(x.size());
//     vPos.resize(n);
//     vNorm.resize(n);
//     for (int i = 0; i < n; ++i) {
//         vPos[i]  = vec3(x[i],  y[i],  z[i]);
//         vNorm[i] = vec3(nx[i], ny[i], nz[i]);
//     }
//
//     // Inizializza i flag di esclusione a 0 (nessuna zona esclusa) —
//     // necessario affinché hasTangentField() funzioni correttamente.
//     vExtEml.assign(n, 0);
//     vEml.assign(n, 0);
//     vScar.assign(n, 0);
//
//     const auto& prov_faces = plyIn.getElement("face").getListProperty<int>("vertex_indices");
//
//     {
//         auto flagF = plyIn.getElement("vertex").getProperty<float>("oor_flag");
//         vOutOfRange.resize(n);
//         for (int i = 0; i < n; ++i)
//             vOutOfRange[i] = (flagF[i] > 0.5f) ? 1 : 0;
//     }
//
//     // Espansione 2-ring: due passate → maschera vertici OOR + 2 anelli di vicini.
//     for (int ring = 0; ring < 2; ++ring) {
//         std::vector<int> expanded = vOutOfRange;
//         for (const auto& pf : prov_faces) {
//             if (vOutOfRange[pf[0]] || vOutOfRange[pf[1]] || vOutOfRange[pf[2]])
//                 expanded[pf[0]] = expanded[pf[1]] = expanded[pf[2]] = 1;
//         }
//         vOutOfRange = std::move(expanded);
//     }
//
//     faces.reserve(prov_faces.size());
//     for (const auto& f : prov_faces) {
//         faces.emplace_back(f[0], f[1], f[2]);
//     }
//
//     std::cout << "importPlyDualLat: " << n << " vertici, "
//               << faces.size() << " facce\n";
//     return true;
// }

void Mesh::exportAsPly(const std::string &filename) const {
    std::vector<float> x(vPos.size()),y(vPos.size()),z(vPos.size()),nx(vPos.size()),ny(vPos.size()),nz(vPos.size());
    for (int i = 0; i < vPos.size(); i++) {
        x[i] = vPos[i].x;
        y[i] = vPos[i].y;
        z[i] = vPos[i].z;
        nx[i] = vNorm[i].x;
        ny[i] = vNorm[i].y;
        nz[i] = vNorm[i].z;
    }
    happly::PLYData plyOut;
    plyOut.addElement("vertex", vPos.size());
    plyOut.getElement("vertex").addProperty<float>("x", x);
    plyOut.getElement("vertex").addProperty<float>("y", y);
    plyOut.getElement("vertex").addProperty<float>("z", z);
    plyOut.getElement("vertex").addProperty<float>("nx", nx);
    plyOut.getElement("vertex").addProperty<float>("ny", ny);
    plyOut.getElement("vertex").addProperty<float>("nz", nz);
    // plyOut.getElement("vertex").addProperty<float>("quality", vQuality);
    std::vector<std::vector<int>> vertex_indices(faces.size());
    for (int i = 0; i < faces.size(); ++i) {
        vertex_indices[i] = {faces[i].a, faces[i].b, faces[i].c};
    }
    plyOut.addElement("face", faces.size());
    plyOut.getElement("face").addListProperty<int>("vertex_indices", vertex_indices);
    plyOut.write(filename, happly::DataFormat::ASCII);
}

// void Mesh::exportAsPly2(const std::string &filename) const {
//     std::vector<float> x(vPos.size()),y(vPos.size()),z(vPos.size()),nx(vPos.size()),ny(vPos.size()),nz(vPos.size());
//     for (int i = 0; i < vPos.size(); i++) {
//         x[i] = vPos[i].x;
//         y[i] = vPos[i].y;
//         z[i] = vPos[i].z;
//         nx[i] = vNorm[i].x;
//         ny[i] = vNorm[i].y;
//         nz[i] = vNorm[i].z;
//     }
//     happly::PLYData plyOut;
//     plyOut.addElement("vertex", vPos.size());
//     plyOut.getElement("vertex").addProperty<float>("x", x);
//     plyOut.getElement("vertex").addProperty<float>("y", y);
//     plyOut.getElement("vertex").addProperty<float>("z", z);
//     plyOut.getElement("vertex").addProperty<float>("nx", nx);
//     plyOut.getElement("vertex").addProperty<float>("ny", ny);
//     plyOut.getElement("vertex").addProperty<float>("nz", nz);
//     plyOut.getElement("vertex").addProperty<float>("quality", vQuality2);
//     std::vector<std::vector<int>> vertex_indices(faces.size());
//     for (int i = 0; i < faces.size(); ++i) {
//         vertex_indices[i] = {faces[i].a, faces[i].b, faces[i].c};
//     }
//     plyOut.addElement("face", faces.size());
//     plyOut.getElement("face").addListProperty<int>("vertex_indices", vertex_indices);
//     plyOut.write(filename, happly::DataFormat::ASCII);
// }

class ObjExporter {
public:
    ObjExporter(const std::string &filename):f(filename){
        f << "OBJ\n";
        vCount = 1;
        nCount = 1;
    }
    ~ObjExporter() {
        if (f.is_open()) f.close();
    }
    void exportVertex(vec3 p) {
        f <<"v " << p.x <<" "<< p.y << ' ' << p.z << '\n';
    }
    void exportNormal(vec3 p) {
        f <<"vn " << p.x <<" "<< p.y << ' ' << p.z << '\n';
    }
    void exportArrow(vec3 p, vec3 v) {
        exportVertex(p);
        exportVertex(p+v);
        exportVertex(p+vec3(0.001, 0,0));
        f <<"f " << vCount++ << ' ' << vCount++ << ' ' << vCount++ << '\n';
    }
    void exportArrow(vec3 p, vec3 v, vec3 n) {
        exportVertex(p);
        exportVertex(p+v);
        exportVertex(p+vec3(0.001, 0,0));
        exportNormal(n);
        f <<"f " << (vCount++) << "//" << nCount << ' ' << vCount++ << "//" << nCount<< ' ' << vCount++<< "//" << nCount << '\n';
        ++nCount;
    }
private:
    std::ofstream f;
    int vCount;
    int nCount;
};

void Mesh::exportTangentsAsObj(std::string filename, int mapIndex) const {
    int nv = faces.size();
    ObjExporter obj(filename);
    for (int i = 0; i < nv; ++i) {
        vec3 v = 1.0f/3.0f * vPos[faces[i].a] + 1.0f/3.0f * vPos[faces[i].b] + 1.0f/3.0f * vPos[faces[i].c];
        obj.exportArrow(v, fTangents[mapIndex][i]);
    }
}

void Mesh::exportAllTangentsAsObj(std::string filename) const{
    for (int i = 0; i < fTangents.size(); ++i){
        int pos = filename.find('.', 0);
        std::string nFilename = filename.substr(0, pos);
        nFilename += std::to_string(i);                          // numero prima del punto
        nFilename += '.';                                // ripristina il punto
        nFilename += filename.substr(pos + 1);
        exportTangentsAsObj(nFilename, i);
    }
}

// void Mesh::exportTangentsAsObj(const std::vector<vec3> vTan, const std::string &filename, float scale, int densify) const {
//     int nf = faces.size();
//     ObjExporter obj(filename);
//     for (int fi = 0; fi < nf; ++fi){
//         for (int i = 0; i < densify; ++i) {
//             for (int j = 0; j < i; ++j) {
//                 vec2 bari(i/float(densify), j/float(densify));
//                 bari += vec2(i/(3.0*densify), j/(3.0*densify));
//                 obj.exportArrow(
//                     posAt(fi,bari),
//                      tangentAt(vTan, fi, bari),
//                      normalAt(fi,bari)
//                     );
//             }
//         }
//     }
// }

// void Mesh::exportTangents1AsObj(const std::string &filename, float scala, int densify) const {
//     if (densify == 0) exportTangentsAsObj(vTangent, filename, scala);
//     else exportTangentsAsObj(vTangent, filename, scala, densify);
// }

// void Mesh::exportTangents2AsObj(const std::string &filename, float scala, int densify) const {
//     if (densify == 0) exportTangentsAsObj(vTangent2, filename, scala);
//     else exportTangentsAsObj(vTangent2, filename, scala, densify);
// }






