//
// Created by gian on 12/19/25.
//

#include <happly.h>
#include "Mesh.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <iomanip>
#include <stdexcept>


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
    mat2 inv = inverse(A);
    vec2 soluzione = inv*k;
    return (soluzione.x * e1 + soluzione.y * e2);
}

void Mesh::updateVertexNormals() {
    vNorm = normalizeVectors(distributeNormalOnVertex(areaVectors()));
}

void Mesh::computeGradients() {
    fTangents.assign(vQualities.size(), std::vector<vec3>(faces.size()));
    for (int k = 0; k < vQualities.size(); ++k)
        for (int i = 0; i < faces.size(); ++i)
            fTangents[k][i] = faceGradientOf(i, k);
}

void Mesh::updateMasks(const EstimationMasks& masks) {
    const int K = fTangents.size();
    if (K != vQualities.size()) throw std::runtime_error("updateMasks: gradienti non calcolati");
    fMask.assign(K, std::vector<uint8_t>(faces.size(), 0));

    for (int k = 0; k < K; ++k) {
        const std::vector<vec3> &g = fTangents[k];
        int nColl = 0;

        for (int i = 0; i < faces.size(); ++i) {
            const Face &face = faces[i];
            const int vi[3] = {face.a, face.b, face.c};

            // collisione: gradiente del vicino aperto nel piano di i attorno allo spigolo condiviso
            const dvec3 gi(g[i]);
            double maxAngle = 0.0;
            for (int s = 0; s < 3; ++s) {
                const int j = faceNeighbors[i][s];
                if (j < 0) continue;
                const int p0 = vi[s], p1 = vi[(s + 1) % 3], pi = vi[(s + 2) % 3];
                const int pj = faces[j].a != p0 && faces[j].a != p1 ? faces[j].a
                             : faces[j].b != p0 && faces[j].b != p1 ? faces[j].b : faces[j].c;
                const dvec3 x0(vPos[p0]);
                const dvec3 e = normalize(dvec3(vPos[p1]) - x0);
                const dvec3 di = dvec3(vPos[pi]) - x0, dj = dvec3(vPos[pj]) - x0;
                const dvec3 bi = normalize(di - dot(di, e) * e);
                const dvec3 bj = normalize(dj - dot(dj, e) * e);
                const dvec3 gj(g[j]);
                const dvec3 gjp = dot(gj, e) * e - dot(gj, bj) * bi;
                const double angle = std::atan2(length(cross(gi, gjp)), dot(gi, gjp)) * 180.0 / M_PI;
                if (std::isnan(angle))
                    throw std::runtime_error("updateMasks: angolo NaN alla faccia " + std::to_string(i) + ", vicino "
                                             + std::to_string(j) + ", mappa " + std::to_string(k) + " (" + qualityNames[k] + ")");
                maxAngle = std::max(maxAngle, angle);
            }
            fMask[k][i] = maxAngle > masks.collisionAngleDeg;
            nColl += fMask[k][i];
        }
        std::cout << "updateMasks: " << qualityNames[k] << ": collisione " << nColl
                  << " su " << faces.size() << " facce\n";
    }
}

void Mesh::updateAnisotropyEstimates(const EstimationMasks& masks) {
    fAniso.assign(faces.size(), FaceAnisotropy());
    const int K = fTangents.size();
    if (fMask.size() != K) throw std::runtime_error("updateAnisotropyEstimates: chiamare prima updateMasks");

    for (int i = 0; i < faces.size(); ++i) {
        FaceAnisotropy &fa = fAniso[i];
        const Face &face = faces[i];

        // solo le righe non escluse, in ordine di mappa
        std::vector<int> used;
        for (int k = 0; k < K; ++k)
            if (!fMask[k][i]) used.push_back(k);
        fa.k_used = used.size();
        if (fa.k_used < kMinGradients) continue;    // sottodeterminato: estimable = false

        // base locale, calcolata una sola volta e usata sia per proiettare sia per ricostruire fiber_dir
        const vec3 e1 = vPos[face.b] - vPos[face.a];
        const vec3 e2 = vPos[face.c] - vPos[face.a];
        const vec3 a2 = cross(e1, e2);
        const vec3 n = a2 / length(a2);
        const vec3 u = e1 / length(e1);
        const vec3 w = cross(n, u);

        std::vector<std::array<double,3>> R(used.size());
        const std::vector<double> b(used.size(), 1.0);
        for (int r = 0; r < used.size(); ++r) {
            const vec3 g = fTangents[used[r]][i];
            const double g1 = dot(g, u), g2 = dot(g, w);
            // embedding di Frobenius (g1^2, sqrt(2) g1 g2, g2^2): una rotazione della base locale
            // agisce come rotazione 3x3, quindi valori singolari e cond non dipendono dalla scelta di u
            R[r] = {g1*g1, std::sqrt(2.0)*g1*g2, g2*g2};
        }

        LsqResult lsq;
        if (!svdLeastSquares(std::move(R), b, masks.minSingularValue, lsq)) continue;
        fa.cond = static_cast<float>(lsq.sigmaMax / lsq.sigmaMin);

        const double M11 = lsq.m[0], M12 = lsq.m[1] / std::sqrt(2.0), M22 = lsq.m[2];
        const double somma = M11 + M22;
        const double differenza = std::sqrt((M11 - M22)*(M11 - M22) + (2.0*M12)*(2.0*M12));
        const double lambda1 = (somma + differenza) / 2.0;   // CV_L^2
        const double lambda2 = (somma - differenza) / 2.0;   // CV_T^2
        if (!(lambda2 > kMinEigenvalue)) continue;

        const double theta = 0.5 * std::atan2(2.0*M12, M11 - M22);
        fa.CV_L = static_cast<float>(std::sqrt(lambda1));
        fa.CV_T = static_cast<float>(std::sqrt(lambda2));
        fa.fiber_dir = static_cast<float>(std::cos(theta))*u + static_cast<float>(std::sin(theta))*w;
        fa.estimable = true;
    }
}

// import export

Mesh::Mesh() = default;

void Mesh::importMultiLatPly(std::string filename) {
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

    // nomi dei siti dai commenti "lat_k NOME": finiscono non quotati nel CSV delle maschere
    for (int k = 0; k < vQualities.size(); ++k) {
        const std::string prefix = "lat_" + std::to_string(k) + " ";
        const auto c = std::find_if(plyIn.comments.begin(), plyIn.comments.end(),
                                    [&](const std::string &s) { return s.rfind(prefix, 0) == 0; });
        if (c == plyIn.comments.end())
            throw std::runtime_error("importMultiLatPly: manca il commento 'lat_" + std::to_string(k) + " NOME'");
        qualityNames.push_back(c->substr(prefix.size()));
    }

    const auto& prov_faces = plyIn.getElement("face").getListProperty<int>("vertex_indices");
    faces.reserve(prov_faces.size());
    for (const auto& f : prov_faces) faces.emplace_back(f[0], f[1], f[2]);

    // adiacenza per spigolo: slot s = spigolo v[s]-v[(s+1)%3], terzo vertice v[(s+2)%3]
    faceNeighbors.assign(faces.size(), {-1, -1, -1});
    std::unordered_map<uint64_t, std::pair<int,int>> firstSeen;   // spigolo -> (faccia, slot) della prima vista
    firstSeen.reserve(faces.size() * 2);
    double edgeLenSum = 0.0;
    for (int i = 0; i < faces.size(); ++i) {
        const int vv[3] = {faces[i].a, faces[i].b, faces[i].c};
        for (int s = 0; s < 3; ++s) {
            const uint32_t p = vv[s], q = vv[(s + 1) % 3];
            const uint64_t key = (uint64_t(std::min(p, q)) << 32) | std::max(p, q);
            const auto [it, inserted] = firstSeen.try_emplace(key, i, s);
            if (inserted) {
                edgeLenSum += distance(vPos[p], vPos[q]);
            } else {
                faceNeighbors[i][s] = it->second.first;
                faceNeighbors[it->second.first][it->second.second] = i;
            }
        }
    }

    std::cout << "importPlyLat: " << vPos.size() << " vertici, "
              << faces.size() << " facce, " << vQualities.size() << " mappe (";
    for (int k = 0; k < qualityNames.size(); ++k) std::cout << (k ? "," : "") << qualityNames[k];
    std::cout << "), lunghezza media spigoli " << edgeLenSum / firstSeen.size() << " mm\n";
}

void Mesh::selectMaps(const std::vector<int>& indices) {
    if (!fTangents.empty()) throw std::runtime_error("selectMaps: gradienti gia' calcolati");
    if (qualityNames.size() != vQualities.size()) throw std::runtime_error("selectMaps: nomi e mappe non allineati");
    std::vector<std::vector<float>> selected;
    std::vector<std::string> selectedNames;
    selected.reserve(indices.size());
    for (int k : indices) {
        if (k < 0 || k >= static_cast<int>(vQualities.size()))
            throw std::runtime_error("selectMaps: mappa " + std::to_string(k) + " fuori range (" +
                                     std::to_string(vQualities.size()) + " mappe)");
        selected.push_back(vQualities[k]);
        selectedNames.push_back(qualityNames[k]);
    }
    vQualities = std::move(selected);
    qualityNames = std::move(selectedNames);
}

void Mesh::importPlyGroundTruth(const std::string& filename, const std::map<int, RegionCV>& cvByRegion) {
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
    for (int i = 0; i < faces.size(); ++i) {
        const Face &face = faces[i];
        const vec3 n = normalize(cross(vPos[face.b] - vPos[face.a], vPos[face.c] - vPos[face.a]));
        const vec3 d(fx[i], fy[i], fz[i]);
        const vec3 t = d - dot(d, n) * n;                    // proiezione sul piano della faccia

        const auto it = cvByRegion.find(tag[i]);
        if (it == cvByRegion.end())
            throw std::runtime_error("importPlyGroundTruth: tag " + std::to_string(tag[i]) + " (faccia "
                                     + std::to_string(i) + ") assente dalla mappatura region_tag -> CV");
        gt[i].fiber_dir = t / length(t);
        gt[i].CV_L = it->second.CV_L;
        gt[i].CV_T = it->second.CV_T;
        gt[i].region_tag = tag[i];
    }
    fGT = std::move(gt);

    std::cout << "importPlyGroundTruth: " << fGT.size() << " facce\n";
}

std::map<int, RegionCV> Mesh::importRegionCV(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("importRegionCV: impossibile aprire " + filename);
    std::map<int, RegionCV> result;
    for (std::string line; std::getline(in, line); ) {
        if (line.empty() || line[0] == '#' || line.rfind("region_tag", 0) == 0) continue;
        std::istringstream ss(line);
        int tag;
        char comma1, comma2;
        RegionCV cv;
        if (!(ss >> tag >> comma1 >> cv.CV_L >> comma2 >> cv.CV_T))
            throw std::runtime_error("importRegionCV: " + filename + ": riga malformata: '" + line + "'");
        result[tag] = cv;
    }
    return result;
}

std::vector<FaceError> Mesh::computeAnisotropyErrors(const std::vector<FaceAnisotropy>& est) const {
    if (est.size() != faces.size() || fGT.size() != faces.size())
        throw std::runtime_error("computeAnisotropyErrors: stima o ground truth non allineate alle facce");
    const float nan = std::numeric_limits<float>::quiet_NaN();
    std::vector<FaceError> err(faces.size());
    for (int i = 0; i < faces.size(); ++i) {
        FaceError &e = err[i];
        e.err_angle = nan;
        if (!est[i].estimable) continue;                      // lasciata marcata, non esclusa ne' riempita

        const Face &face = faces[i];
        const vec3 n = normalize(cross(vPos[face.b] - vPos[face.a], vPos[face.c] - vPos[face.a]));
        // entrambe gia' nel piano: la stimata e' costruita nella base {u, w}, la vera proiettata all'import
        const vec3 a = est[i].fiber_dir;
        const vec3 b = fGT[i].fiber_dir;
        // angolo nel piano, poi modulo pi: nessuna scelta di orientamento
        const double delta = std::abs(std::atan2(double(dot(n, cross(a, b))), double(dot(a, b))));
        e.err_angle = static_cast<float>(std::min(delta, M_PI - delta) * 180.0 / M_PI);
        e.valid = true;
    }
    return err;
}

void Mesh::exportErrorsAsCsv(const std::string& filename, const std::string& experimentId,
                             const std::vector<FaceAnisotropy>& est) const {
    const std::vector<FaceError> err = computeAnisotropyErrors(est);
    const std::vector<vec3> av = areaVectors();                // cross(e1,e2): area = |cross| / 2
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("exportErrorsAsCsv: impossibile aprire " + filename);
    out << std::scientific << std::setprecision(9);
    out << "experiment_id,face_id,area,region_tag,estimable,k_used,cond,CV_L_est,CV_T_est,CV_L_true,CV_T_true,"
           "err_angle\n";
    for (int i = 0; i < faces.size(); ++i) {
        const FaceAnisotropy &a = est[i];
        const FaceGroundTruth &g = fGT[i];
        const FaceError &e = err[i];
        out << experimentId << ',' << i << ',' << length(av[i]) / 2.0f << ',' << g.region_tag << ',' << (a.estimable ? 1 : 0) << ','
            << a.k_used << ',';
        if (a.cond > 0.0f) out << a.cond;                       // vuota se non calcolata
        out << ',';
        if (a.estimable) out << a.CV_L << ',' << a.CV_T << ','; else out << ",,";
        out << g.CV_L << ',' << g.CV_T << ',';
        if (e.valid) out << e.err_angle;
        out << '\n';
    }
    out.close();
    if (!out) throw std::runtime_error("exportErrorsAsCsv: errore di scrittura su " + filename);
}

void Mesh::exportMaskPairsAsCsv(const std::string& filename, const std::string& experimentId) const {
    const int K = fTangents.size();
    if (fMask.size() != K || qualityNames.size() != K || fGT.size() != faces.size())
        throw std::runtime_error("exportMaskPairsAsCsv: maschere, nomi o ground truth non allineati");
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("exportMaskPairsAsCsv: impossibile aprire " + filename);
    out << std::scientific << std::setprecision(5);                // 6 cifre significative
    out << "experiment_id,face_id,map_name,masked,r\n";
    for (int k = 0; k < K; ++k) {
        for (int i = 0; i < faces.size(); ++i) {
            // r = g^T M_true g - 1, M_true = CV_L^2 f f^T + CV_T^2 t t^T: solo valutazione, mai nella stima
            const Face &face = faces[i];
            const dvec3 n = normalize(dvec3(cross(vPos[face.b] - vPos[face.a], vPos[face.c] - vPos[face.a])));
            const dvec3 f(fGT[i].fiber_dir);                       // gia' proiettata e normalizzata all'import
            const dvec3 t = cross(n, f);
            const dvec3 g(fTangents[k][i]);
            const double gf = dot(g, f), gt = dot(g, t);
            const double cvL = fGT[i].CV_L, cvT = fGT[i].CV_T;
            const double r = cvL*cvL*gf*gf + cvT*cvT*gt*gt - 1.0;
            out << experimentId << ',' << i << ',' << qualityNames[k] << ',' << int(fMask[k][i]) << ',' << r << '\n';
        }
    }
    out.close();
    if (!out) throw std::runtime_error("exportMaskPairsAsCsv: errore di scrittura su " + filename);
}

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
    std::vector<std::vector<int>> vertex_indices(faces.size());
    for (int i = 0; i < faces.size(); ++i) {
        vertex_indices[i] = {faces[i].a, faces[i].b, faces[i].c};
    }
    plyOut.addElement("face", faces.size());
    plyOut.getElement("face").addListProperty<int>("vertex_indices", vertex_indices);
    plyOut.write(filename, happly::DataFormat::ASCII);
}

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
