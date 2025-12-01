#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <algorithm>
#include <regex>

using namespace std;


struct Grid {
    double* data = nullptr;   // matriz lineal Ny * Nx
    double* x = nullptr;      // eje X (longitud)
    double* y = nullptr;      // eje Y (latitud)
    size_t Nx = 0;
    size_t Ny = 0;
    double resolution = 0.0;
    double lat_min = 0.0, lat_max = 0.0;
    double lon_min = 0.0, lon_max = 0.0;

    // Constructor vacío
    Grid() = default;

    // Constructor con tamaño
    Grid(size_t ny, size_t nx) : Nx(nx), Ny(ny) {
        data = new double[Ny * Nx];
        x = new double[Nx];
        y = new double[Ny];
        std::fill(data, data + Ny * Nx, std::numeric_limits<double>::quiet_NaN());
    }

    // Destructor
    ~Grid() {
        delete[] data;
        delete[] x;
        delete[] y;
    }

    // Copia profunda
    Grid(const Grid& other) {
        Nx = other.Nx;
        Ny = other.Ny;
        resolution = other.resolution;
        lat_min = other.lat_min;
        lat_max = other.lat_max;
        lon_min = other.lon_min;
        lon_max = other.lon_max;

        x = new double[Nx];
        y = new double[Ny];
        data = new double[Ny * Nx];
        std::copy(other.x, other.x + Nx, x);
        std::copy(other.y, other.y + Ny, y);
        std::copy(other.data, other.data + Ny * Nx, data);
    }

    // Asignación profunda
    Grid& operator=(const Grid& other) {
        if (this == &other) return *this;
        delete[] data;
        delete[] x;
        delete[] y;
        Nx = other.Nx;
        Ny = other.Ny;
        resolution = other.resolution;
        lat_min = other.lat_min;
        lat_max = other.lat_max;
        lon_min = other.lon_min;
        lon_max = other.lon_max;

        x = new double[Nx];
        y = new double[Ny];
        data = new double[Ny * Nx];
        std::copy(other.x, other.x + Nx, x);
        std::copy(other.y, other.y + Ny, y);
        std::copy(other.data, other.data + Ny * Nx, data);
        return *this;
    }
};


struct Config {
    double cel = 0.01;
    double rel = 8.0;
    double ratio = 10.0;
    vector<vector<vector<string>>> nombres; 
    // niveles → grupos → archivos
};

using Nivel = vector<Grid>;
using MultiNivel = vector<Nivel>;


void linspace(double* vec, size_t n, double min, double max);
void imprimirLimites(const Grid &g);
Grid leerBat(const string &nombreArchivo, double celda);
void rellenar(Grid &g);
void exportarCSV(const Grid &g, const string &nombre);
Grid fusionarGrids(const vector<Grid> &grids);
void adaptarGrid(const Grid &coarse, Grid &fine, double ratio);
Grid fusionarMallas(const Nivel &mallasNivel, const Nivel &mallasSuperior, double ratio);
Grid procesoNivel(const string &entrada, const string &salida, double celda);
Grid leerCSV(const string &nombreArchivo);
string trim(const string& s);
vector<string> split(const string& s, char delim);
Config leer_config_avanzado(const string& archivo);



int main() {
    try {
        Config cfg = leer_config_avanzado("config0.txt");

        double cel = cfg.cel;
        double rel = cfg.rel;
        double ratio = cfg.ratio;

        vector<vector<Grid>> niveles;
        double resol_actual = cel;

        for (size_t nivel = 0; nivel < cfg.nombres.size(); ++nivel) {
            cout << "\n===============================" << endl;
            cout << "   🧩 Procesando Nivel " << nivel << endl;
            cout << "===============================" << endl;

            vector<Grid> resultado_nivel;

            for (const auto& grupo : cfg.nombres[nivel]) {
                vector<Grid> mallasGrupo;
                for (const auto& archivo : grupo) {
                    cout << "📄 Leyendo: " << archivo << endl;
                    Grid g = leerBat(archivo, resol_actual);
                    rellenar(g);
                    mallasGrupo.push_back(std::move(g));
                }

                if (nivel == 0) {
                    resultado_nivel.push_back(fusionarGrids(mallasGrupo));
                } else {
                    resultado_nivel.push_back(
                        fusionarMallas(mallasGrupo, niveles.back(), ratio)
                    );
                }
            }

            niveles.push_back(std::move(resultado_nivel));
            resol_actual /= rel;
        }

        // Exportar todos los niveles
        for (size_t n = 0; n < niveles.size(); ++n) {
            for (size_t i = 0; i < niveles[n].size(); ++i) {
                string nombre = "nivel" + to_string(n) + "_" + to_string(i);
                exportarCSV(niveles[n][i], nombre);
            }
        }

        cout << "\n✅ Proceso completado correctamente.\n";

    } catch (const exception& e) {
        cerr << "❌ Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}



// FUNCIONES


void linspace(double* vec, size_t n, double min, double max) {
    if (n <= 1) return;
    double step = (max - min) / (n - 1);
    for (size_t i = 0; i < n; ++i)
        vec[i] = min + i * step;
}

void imprimirLimites(const Grid &g) {
    cout << fixed << setprecision(6);
    cout << "   Latitud:  [" << g.lat_min << ", " << g.lat_max << "]\n";
    cout << "   Longitud: [" << g.lon_min << ", " << g.lon_max << "]\n";
    cout << "   Resolución: " << g.resolution << "°\n";
    cout << "   Dimensiones: " << g.Ny << " x " << g.Nx << "\n";
}

Grid leerBat(const string &nombreArchivo, double celda) {
    ifstream file(nombreArchivo);
    if (!file.is_open())
        throw runtime_error("No se pudo abrir el archivo: " + nombreArchivo);

    double lon, lat, prof;
    double lon_min = numeric_limits<double>::max(), lon_max = numeric_limits<double>::lowest();
    double lat_min = numeric_limits<double>::max(), lat_max = numeric_limits<double>::lowest();

    vector<tuple<double, double, double>> puntos;

    while (file >> lon >> lat >> prof) {
        lon_min = min(lon_min, lon);
        lon_max = max(lon_max, lon);
        lat_min = min(lat_min, lat);
        lat_max = max(lat_max, lat);
        puntos.emplace_back(lon, lat, prof);
    }

    file.close();

    lon_min = floor(lon_min / celda) * celda;
    lon_max = ceil(lon_max / celda) * celda;
    lat_min = floor(lat_min / celda) * celda;
    lat_max = ceil(lat_max / celda) * celda;

    size_t Nx = static_cast<size_t>((lon_max - lon_min) / celda) + 1;
    size_t Ny = static_cast<size_t>((lat_max - lat_min) / celda) + 1;

    Grid g(Ny, Nx);
    g.resolution = celda;
    g.lat_min = lat_min;
    g.lat_max = lat_max;
    g.lon_min = lon_min;
    g.lon_max = lon_max;

    linspace(g.y, Ny, lat_min, lat_max);
    linspace(g.x, Nx, lon_min, lon_max);

    for (auto &[lo, la, pr] : puntos) {
        int i = round((la - lat_min) / celda);
        int j = round((lo - lon_min) / celda);
        if (i >= 0 && i < (int)Ny && j >= 0 && j < (int)Nx)
            g.data[i * Nx + j] = pr;
    }

    imprimirLimites(g);
    return g;
}

void rellenar(Grid &g) {
    const int MAX_ITER = 10;
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        bool cambios = false;
        vector<double> copia(g.data, g.data + g.Nx * g.Ny);

        for (size_t i = 0; i < g.Ny; ++i) {
            for (size_t j = 0; j < g.Nx; ++j) {
                double &val = g.data[i * g.Nx + j];
                if (isnan(val)) {
                    double suma = 0; int n = 0;
                    if (i > 0 && !isnan(g.data[(i-1)*g.Nx + j])) { suma += g.data[(i-1)*g.Nx + j]; n++; }
                    if (i < g.Ny-1 && !isnan(g.data[(i+1)*g.Nx + j])) { suma += g.data[(i+1)*g.Nx + j]; n++; }
                    if (j > 0 && !isnan(g.data[i*g.Nx + (j-1)])) { suma += g.data[i*g.Nx + (j-1)]; n++; }
                    if (j < g.Nx-1 && !isnan(g.data[i*g.Nx + (j+1)])) { suma += g.data[i*g.Nx + (j+1)]; n++; }
                    if (n > 0) { copia[i * g.Nx + j] = suma / n; cambios = true; }
                }
            }
        }

        copy(copia.begin(), copia.end(), g.data);
        if (!cambios) break;
    }
    cout << "✅ Interpolación completada.\n";
}

void exportarCSV(const Grid &g, const string &nombre) {
    string ruta = "topobatimetrias/" + nombre + ".csv";
    ofstream file(ruta);
    if (!file.is_open()) {
        cerr << "No se pudo escribir: " << ruta << "\n";
        return;
    }

    file << fixed << setprecision(6);
    file << "0";
    for (size_t j = 0; j < g.Nx; ++j) file << "," << g.x[j];
    file << "\n";

    for (size_t i = 0; i < g.Ny; ++i) {
        file << g.y[i];
        for (size_t j = 0; j < g.Nx; ++j)
            file << "," << g.data[i * g.Nx + j];
        file << "\n";
    }

    file.close();
    cout << "  CSV exportado: " << ruta << endl;
}


Grid leerCSV(const string &nombreArchivo) {
    ifstream file(nombreArchivo);
    if (!file.is_open()) throw runtime_error("No se pudo abrir el CSV: " + nombreArchivo);

    string line;
    getline(file, line);
    stringstream ss(line);
    string tmp;

    // leer longitudes (x)
    vector<double> lon;
    getline(ss, tmp, ','); // saltar cabecera
    while (getline(ss, tmp, ',')) lon.push_back(stod(tmp));

    vector<double> lat;
    vector<double> datos_flat;

    // leer cada fila
    while (getline(file, line)) {
        stringstream ls(line);
        string token;
        getline(ls, token, ',');
        double la = stod(token);
        lat.push_back(la);

        size_t col = 0;
        while (getline(ls, token, ',')) {
            double val = stod(token);
            datos_flat.push_back(val);
            col++;
        }
    }

    size_t Nx = lon.size();
    size_t Ny = lat.size();

    Grid g(Ny, Nx);
    g.resolution = (Nx > 1) ? (lon[1] - lon[0]) : 0;
    g.lat_min = lat.front();
    g.lat_max = lat.back();
    g.lon_min = lon.front();
    g.lon_max = lon.back();

    // Copiar ejes
    for (size_t i = 0; i < Nx; ++i) g.x[i] = lon[i];
    for (size_t i = 0; i < Ny; ++i) g.y[i] = lat[i];

    // Copiar datos
    std::copy(datos_flat.begin(), datos_flat.end(), g.data);

    return g;
}


Grid fusionarGrids(const vector<Grid> &grids) {
    if (grids.empty()) throw runtime_error("fusionarGrids: sin mallas para fusionar");

    double lat_min = numeric_limits<double>::max(), lat_max = numeric_limits<double>::lowest();
    double lon_min = numeric_limits<double>::max(), lon_max = numeric_limits<double>::lowest();
    double res = grids[0].resolution;

    for (const auto &g : grids) {
        lat_min = min(lat_min, g.lat_min);
        lat_max = max(lat_max, g.lat_max);
        lon_min = min(lon_min, g.lon_min);
        lon_max = max(lon_max, g.lon_max);
    }

    size_t Nx = static_cast<size_t>((lon_max - lon_min) / res) + 1;
    size_t Ny = static_cast<size_t>((lat_max - lat_min) / res) + 1;

    Grid out(Ny, Nx);
    out.resolution = res;
    out.lat_min = lat_min;
    out.lat_max = lat_max;
    out.lon_min = lon_min;
    out.lon_max = lon_max;

    linspace(out.x, Nx, lon_min, lon_max);
    linspace(out.y, Ny, lat_min, lat_max);

    // Inicializar con NaN
    fill(out.data, out.data + Ny * Nx, numeric_limits<double>::quiet_NaN());

    for (const auto &g : grids) {
        for (size_t i = 0; i < g.Ny; ++i) {
            for (size_t j = 0; j < g.Nx; ++j) {
                double val = g.data[i * g.Nx + j];
                if (!isnan(val)) {
                    size_t i2 = round((g.y[i] - lat_min) / res);
                    size_t j2 = round((g.x[j] - lon_min) / res);
                    if (i2 < Ny && j2 < Nx) out.data[i2 * Nx + j2] = val;
                }
            }
        }
    }

    cout << "✅ Mallas fusionadas correctamente (" << grids.size() << " submallas)\n";
    return out;
}

// Ajusta una malla fina para que no sea demasiado pequeña respecto a la gruesa
void adaptarGrid(const Grid &coarse, Grid &fine, double ratio) {
    // Tamaño de la malla fina
    double lat_fine_size = fine.lat_max - fine.lat_min;
    double lon_fine_size = fine.lon_max - fine.lon_min;

    // Tamaño mínimo permitido según malla gruesa y ratio
    double lat_min_size = (coarse.lat_max - coarse.lat_min) / ratio;
    double lon_min_size = (coarse.lon_max - coarse.lon_min) / ratio;

    // Verificar si es necesario ampliar
    if (lat_fine_size >= lat_min_size && lon_fine_size >= lon_min_size) {
        cout << "✅ Malla ya cumple el tamaño mínimo, no se amplía\n";
        return;
    }

    // Nueva extensión de la malla fina
    double new_lat_size = max(lat_fine_size, lat_min_size);
    double new_lon_size = max(lon_fine_size, lon_min_size);

    // Calcular nuevo número de celdas (manteniendo resolución original)
    size_t Ny_new = static_cast<size_t>(round(new_lat_size / fine.resolution)) + 1;
    size_t Nx_new = static_cast<size_t>(round(new_lon_size / fine.resolution)) + 1;

    Grid enlarged(Ny_new, Nx_new);
    enlarged.resolution = fine.resolution;

    // Centrar la malla fina dentro de la malla ampliada
    double lat_center = (fine.lat_min + fine.lat_max) / 2.0;
    double lon_center = (fine.lon_min + fine.lon_max) / 2.0;

    enlarged.lat_min = lat_center - new_lat_size / 2.0;
    enlarged.lat_max = lat_center + new_lat_size / 2.0;
    enlarged.lon_min = lon_center - new_lon_size / 2.0;
    enlarged.lon_max = lon_center + new_lon_size / 2.0;

    linspace(enlarged.y, Ny_new, enlarged.lat_min, enlarged.lat_max);
    linspace(enlarged.x, Nx_new, enlarged.lon_min, enlarged.lon_max);

    // Calcular el desplazamiento en celdas
    int i_offset = static_cast<int>(round((fine.lat_min - enlarged.lat_min) / fine.resolution));
    int j_offset = static_cast<int>(round((fine.lon_min - enlarged.lon_min) / fine.resolution));

    // Copiar datos de la malla fina a la ampliada sin distorsión
    for (size_t i = 0; i < fine.Ny; ++i) {
        for (size_t j = 0; j < fine.Nx; ++j) {
            size_t i_new = i + i_offset;
            size_t j_new = j + j_offset;
            if (i_new < Ny_new && j_new < Nx_new) {
                enlarged.data[i_new * Nx_new + j_new] = fine.data[i * fine.Nx + j];
            }
        }
    }

    // Rellenar los bordes con datos de la malla gruesa donde haya NaN
    for (size_t i = 0; i < Ny_new; ++i) {
        for (size_t j = 0; j < Nx_new; ++j) {
            double &val = enlarged.data[i * Nx_new + j];
            if (isnan(val)) {
                double lat = enlarged.y[i];
                double lon = enlarged.x[j];
                if (lat >= coarse.lat_min && lat <= coarse.lat_max &&
                    lon >= coarse.lon_min && lon <= coarse.lon_max) {
                    size_t si = static_cast<size_t>((lat - coarse.lat_min) / coarse.resolution);
                    size_t sj = static_cast<size_t>((lon - coarse.lon_min) / coarse.resolution);
                    if (si < coarse.Ny && sj < coarse.Nx) {
                        double valSup = coarse.data[si * coarse.Nx + sj];
                        if (!isnan(valSup)) val = valSup;
                    }
                }
            }
        }
    }

    fine = enlarged; // Sustituir la malla fina por la ampliada
    cout << "✅ Malla adaptada correctamente: nueva dimensión [" << Ny_new << " x " << Nx_new << "]\n";
}


// Fusiona mallas del nivel actual y rellena con los datos del nivel superior
Grid fusionarMallas(const Nivel &mallasNivel, const Nivel &mallasSuperior, double ratio) {
    if (mallasNivel.empty())
        throw runtime_error("fusionarMallas: sin mallas del nivel actual");

    cout << "\n🔷 Fusionando " << mallasNivel.size() << " mallas de este nivel...\n";

    // 1️⃣ Fusionar todas las mallas del nivel actual
    Grid fusionado = fusionarGrids(mallasNivel);

    // 2️⃣ Adaptar respecto al nivel superior
    for (const auto &sup : mallasSuperior)
        adaptarGrid(sup, fusionado, ratio);

    // 3️⃣ Rellenar huecos con datos del nivel superior
    for (const auto &sup : mallasSuperior) {
        for (size_t i = 0; i < fusionado.Ny; ++i) {
            for (size_t j = 0; j < fusionado.Nx; ++j) {
                double &val = fusionado.data[i * fusionado.Nx + j];
                if (isnan(val)) {
                    double lat = fusionado.y[i];
                    double lon = fusionado.x[j];
                    if (lat >= sup.lat_min && lat <= sup.lat_max &&
                        lon >= sup.lon_min && lon <= sup.lon_max) {

                        size_t si = static_cast<size_t>((lat - sup.lat_min) / sup.resolution);
                        size_t sj = static_cast<size_t>((lon - sup.lon_min) / sup.resolution);
                        if (si < sup.Ny && sj < sup.Nx) {
                            double valSup = sup.data[si * sup.Nx + sj];
                            if (!isnan(valSup))
                                val = valSup;
                        }
                    }
                }
            }
        }
    }

    rellenar(fusionado);
    cout << "✅ Nivel fusionado y rellenado con datos superiores\n";
    return fusionado;
}

// Función para eliminar espacios en blanco al inicio y final de una cadena
string trim(const string& s) {
    size_t first = s.find_first_not_of(" \t\r\n");
    if (first == string::npos) return "";
    size_t last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

// Función para dividir una cadena por un delimitador
vector<string> split(const string& s, char delim) {
    vector<string> result;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        item = trim(item);
        if (!item.empty()) {
            result.push_back(item);
        }
    }
    return result;
}

Config leer_config_avanzado(const string& archivo) {
    Config cfg;
    ifstream file(archivo);
    if (!file.is_open())
        throw runtime_error("No se pudo abrir el archivo de configuración: " + archivo);

    string line;
    bool in_nombres = false;
    bool in_grupo = false;

    vector<vector<string>> nivel_actual;  // grupos del nivel actual
    vector<string> grupo_actual;          // archivos del grupo actual

    regex cel_re("^cel\\s*=\\s*([0-9.]+)");
    regex rel_re("^rel\\s*=\\s*([0-9.]+)");
    regex ratio_re("^ratio\\s*=\\s*([0-9.]+)");
    regex nombres_re("^nombres\\s*=");

    while (getline(file, line)) {
        line = regex_replace(line, regex("^\\s+|\\s+$"), ""); // trim
        if (line.empty() || line[0] == '#') {
            // Detectar "# Nivel X"
            if (line.rfind("# Nivel", 0) == 0) {
                if (!nivel_actual.empty()) {
                    cfg.nombres.push_back(nivel_actual);
                    nivel_actual.clear();
                }
            }
            continue;
        }

        smatch m;
        if (regex_search(line, m, cel_re)) { cfg.cel = stod(m[1]); continue; }
        if (regex_search(line, m, rel_re)) { cfg.rel = stod(m[1]); continue; }
        if (regex_search(line, m, ratio_re)) { cfg.ratio = stod(m[1]); continue; }
        if (regex_search(line, m, nombres_re)) { in_nombres = true; continue; }

        if (in_nombres) {
            if (line == "{") {
                in_grupo = true;
                grupo_actual.clear();
                continue;
            } else if (line == "}") {
                if (!grupo_actual.empty()) {
                    nivel_actual.push_back(grupo_actual);
                    grupo_actual.clear();
                }
                in_grupo = false;
                continue;
            }

            if (in_grupo) {
                vector<string> archivos;
                stringstream ss(line);
                string token;
                while (getline(ss, token, ',')) {
                    token = regex_replace(token, regex("^\\s+|\\s+$"), "");
                    if (!token.empty())
                        archivos.push_back(token);
                }
                grupo_actual.insert(grupo_actual.end(), archivos.begin(), archivos.end());
            }
        }
    }

    // ⚠️ Añadir el último nivel que quedó pendiente
    if (!grupo_actual.empty())
        nivel_actual.push_back(grupo_actual);
    if (!nivel_actual.empty())
        cfg.nombres.push_back(nivel_actual);

    file.close();

    // Mostrar resumen
    cout << "Cel: " << cfg.cel << "  Rel: " << cfg.rel << "  Ratio: " << cfg.ratio << endl;
    cout << "  Niveles: " << cfg.nombres.size() << endl;

    for (size_t i = 0; i < cfg.nombres.size(); ++i) {
        cout << "   Nivel " << i << " con " << cfg.nombres[i].size() << " grupos" << endl;
        for (size_t g = 0; g < cfg.nombres[i].size(); ++g) {
            cout << "     Grupo " << g << ": ";
            for (auto& n : cfg.nombres[i][g]) cout << n << " ";
            cout << endl;
        }
    }

    return cfg;
}
