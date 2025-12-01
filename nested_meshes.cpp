#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <sys/stat.h>  
#include <algorithm>
#include <regex>

using namespace std;


struct Data {
    double* data = nullptr;   // matriz lineal Ny * Nx
    double* x = nullptr;      // eje X
    double* y = nullptr;      // eje Y
    size_t Nx = 0;
    size_t Ny = 0;
};

struct SubBlock {
    double* data = nullptr;
    double* x = nullptr;
    double* y = nullptr;
    size_t bloque = 0;
    size_t i0 = 0;
    size_t j0 = 0;
};

struct Config {
    vector<int> numporlevel;
    vector<int> ratios;
    vector<vector<string>> nombres;
    vector<vector<int>> bloque;
};


Data devolvercsv(const string& filename);
Data comprimir(const Data& datos, size_t bloque);
vector<SubBlock> extraerSubbloquesDesdeRed(const Data& grid, const Data& red, size_t bloque);
vector<SubBlock> refinar(const vector<SubBlock>& bloques, const string& nombre2, double factor_refinamiento)  ;
void linspace(double* arr, size_t N, double min, double max);
double* redimensionar(const double* mat, size_t tam, double rel, size_t& nuevo_tam);
Data pegar_refinados(const vector<SubBlock>& bloques_refinados, const Data& datos_superior);
void exportar_csv(const Data& datos, const string& filename, bool print);
void exportar_subbloques_csv(const vector<SubBlock>& bloques, const string& carpeta);
static double coord_tol(double range) ;
void liberarData(Data& d);
void proceso(int numlevel, const vector<int>& numporlevel, const vector<int>& ratios, const vector<vector<string>>& nombres,const vector<vector<int>>& bloques_por_level);
string trim(const string& s);
vector<string> split(const string& s, char delim) ;
Config leer_config(const string& archivo);

int main() {
    Config cfg = leer_config("config.txt");
    if (cfg.nombres.empty()) {
        cerr << "Error: no se cargó la configuración correctamente." << endl;
        return 1;
    }
    cout << "Configuración cargada ✅" << endl;
    cout << "Niveles: " << cfg.numporlevel.size() << endl;
    proceso(cfg.numporlevel.size(), cfg.numporlevel, cfg.ratios, cfg.nombres, cfg.bloque);
    return 0;
}



Data devolvercsv(const string& filename) {
    ifstream file(filename);
    Data result;

    if (!file.is_open()) {
        cerr << "❌ Error: no se pudo abrir el archivo " << filename << endl;
        return result;
    }

    vector<vector<double>> temp;
    string line;

    // --- Leer todas las líneas ---
    while (getline(file, line)) {
        if (line.empty()) continue;

        stringstream ss(line);
        string campo;
        vector<double> fila;

        while (getline(ss, campo, ',')) {
            if (campo.empty()) {
                fila.push_back(0.0);
            } else {
                try {
                    fila.push_back(stod(campo));
                } catch (...) {
                    fila.push_back(0.0);
                }
            }
        }

        temp.push_back(fila);
    }

    file.close();

    if (temp.empty() || temp.size() < 2 || temp[0].size() < 2) {
        cerr << "❌ Formato CSV inválido o vacío: " << filename << endl;
        return result;
    }

    // --- Extraer dimensiones ---
    result.Nx = temp[0].size() - 1; // primera fila tiene Nx+1 valores (el 0 inicial)
    result.Ny = temp.size() - 1;    // primera columna tiene Ny+1 valores (el 0 inicial)

    // --- Reservar memoria ---
    result.x = new double[result.Nx];
    result.y = new double[result.Ny];
    result.data = new double[result.Nx * result.Ny];

    // --- Leer ejes X ---
    for (size_t j = 0; j < result.Nx; ++j)
        result.x[j] = temp[0][j + 1];

    // --- Leer ejes Y y matriz de datos ---
    for (size_t i = 0; i < result.Ny; ++i) {
        result.y[i] = temp[i + 1][0];
        for (size_t j = 0; j < result.Nx; ++j) {
            result.data[i * result.Nx + j] = temp[i + 1][j + 1];
        }
    }

    return result;
}

Data comprimir(const Data& datos, size_t bloque) {
    // Ajustar a múltiplos del bloque
    size_t ny_c = datos.Ny - (datos.Ny % bloque);
    size_t nx_c = datos.Nx - (datos.Nx % bloque);

    size_t ny_red = ny_c / bloque;
    size_t nx_red = nx_c / bloque;

    Data datos_red;
    datos_red.Nx = nx_red;
    datos_red.Ny = ny_red;

    // Reservar memoria lineal
    datos_red.data = new double[ny_red * nx_red];
    datos_red.x = new double[nx_red];
    datos_red.y = new double[ny_red];

    // Compresión principal
    for (size_t i = 0; i < ny_red; ++i) {
        for (size_t j = 0; j < nx_red; ++j) {
            double suma = 0.0;
            bool nom = false;

            double ref = datos.data[(i * bloque) * datos.Nx + (j * bloque)];
            int signo_ref = (ref > 0) - (ref < 0);

            // Bloque compacto
            const size_t base_i = i * bloque;
            const size_t base_j = j * bloque;

            for (size_t bi = 0; bi < bloque; ++bi) {
                const size_t idx_row = (base_i + bi) * datos.Nx + base_j;
                for (size_t bj = 0; bj < bloque; ++bj) {
                    double val = datos.data[idx_row + bj];
                    int signo_val = (val > 0) - (val < 0);
                    if (signo_ref != signo_val)
                        nom = true;
                    suma += val;
                }
            }

            datos_red.data[i * nx_red + j] =
                nom ? -1.0 : ((suma / (bloque * bloque)) > 0 ? 1.0 : 0.0);
        }
    }

    // Reducción de coordenadas
    for (size_t i = 0; i < ny_red; ++i)
        datos_red.y[i] = datos.y[i * bloque + bloque / 2];

    for (size_t j = 0; j < nx_red; ++j)
        datos_red.x[j] = datos.x[j * bloque + bloque / 2];

    return datos_red;
}



vector<SubBlock> extraerSubbloquesDesdeRed(const Data& grid, const Data& red, size_t bloque) {
    vector<SubBlock> bloques;

    if (!grid.data || !red.data) {
        cerr << "❌ Error: estructuras sin datos." << endl;
        return bloques;
    }

    for (size_t i = 0; i < red.Ny; ++i) {
        for (size_t j = 0; j < red.Nx; ++j) {
            double val = red.data[i * red.Nx + j];
            if (val == -1.0) { // zona costera detectada
                SubBlock sb;
                sb.bloque = bloque;
                sb.i0 = i;
                sb.j0 = j;

                sb.data = new double[bloque * bloque];
                sb.x = new double[bloque];
                sb.y = new double[bloque];

                for (size_t bi = 0; bi < bloque; ++bi) {
                    memcpy(sb.data + bi * bloque,
                           grid.data + (sb.i0*bloque + bi) * grid.Nx + sb.j0*bloque,
                           bloque * sizeof(double));
                }

                for (size_t bi = 0; bi < bloque; ++bi)
                    sb.y[bi] = grid.y[sb.i0*bloque + bi];
                for (size_t bj = 0; bj < bloque; ++bj)
                    sb.x[bj] = grid.x[sb.j0*bloque + bj];

                bloques.push_back(sb);
            }
        }
    }

    cout << "[extraerSubbloquesDesdeRed] Subbloques detectados: " << bloques.size() << endl;
    return bloques;
}

void linspace(double* arr, size_t N, double min, double max) {
    if (N == 0) return;
    if (N == 1) {
        arr[0] = min;
        return;
    }
    double step = (max - min) / (N - 1);
    for (size_t i = 0; i < N; ++i)
        arr[i] = min + i * step;
}

double* redimensionar(const double* mat, size_t tam, double rel, size_t& nuevo_tam) {
    if (tam < 2) {
        nuevo_tam = tam;
        double* res = new double[nuevo_tam];
        if (tam == 1) res[0] = mat[0];
        return res;
    }

    double dx = (mat[tam - 1] - mat[0]) / (tam - 1);

    int it = static_cast<int>(log(rel) / log(2));
    double des = 0.0;
    for (int i = 0; i < abs(it); ++i)
        des += 1.0 / pow(2.0, i);

    nuevo_tam = static_cast<size_t>(tam * rel);
    double* mat_mod = new double[nuevo_tam];

    double min = mat[0] - dx * des / 4.0;
    double max = mat[tam - 1] + dx * des / 4.0;

    linspace(mat_mod, nuevo_tam, min, max);

    return mat_mod;
}



vector<SubBlock> refinar(const vector<SubBlock>& bloques, const string& nombre2, double factor_refinamiento) {
    vector<SubBlock> bloques_refinados;

    // Abrir CSV con datos finos
    Data datos_finos = devolvercsv(nombre2);
    if (!datos_finos.data) {
        cerr << "Error al abrir archivo de referencia: " << nombre2 << endl;
        return bloques_refinados;
    }

    for (const auto& blk : bloques) {
        bool refinable = false;

        // Verificar si las coordenadas del sub-bloque están dentro del grid fino
        double x_min = blk.x[0];
        double x_max = blk.x[blk.bloque - 1];
        double y_min = blk.y[0];
        double y_max = blk.y[blk.bloque - 1];

        if (x_min >= datos_finos.x[0] && x_max <= datos_finos.x[datos_finos.Nx - 1] &&
            y_min >= datos_finos.y[0] && y_max <= datos_finos.y[datos_finos.Ny - 1]) {
            refinable = true;
        }

        if (refinable) {
            SubBlock blk_ref;
            blk_ref.i0 = blk.i0;
            blk_ref.j0 = blk.j0;
        
            size_t nuevo_tam_x, nuevo_tam_y;
            blk_ref.bloque = static_cast<size_t>(blk.bloque * factor_refinamiento);
            //cout << endl << "el tamaño del bloque es" << blk_ref.bloque;
            // Redimensionar coordenadas usando el factor de refinamiento
            vector<double> x_temp(blk.bloque), y_temp(blk.bloque);
            for (size_t i = 0; i < blk.bloque; i++) {
                x_temp[i] = blk.x[i];
                y_temp[i] = blk.y[i];
            }
        
            double* x_red = redimensionar(x_temp.data(), blk.bloque, factor_refinamiento, nuevo_tam_x);
            double* y_red = redimensionar(y_temp.data(), blk.bloque, factor_refinamiento, nuevo_tam_y);
        
            blk_ref.x = new double[nuevo_tam_x];
            blk_ref.y = new double[nuevo_tam_y];
        
            memcpy(blk_ref.x, x_red, nuevo_tam_x * sizeof(double));
            memcpy(blk_ref.y, y_red, nuevo_tam_y * sizeof(double));
        
            blk_ref.data = new double[nuevo_tam_x * nuevo_tam_y];
        
            // Rellenar datos desde datos_finos
            for (size_t i = 0; i < nuevo_tam_y; i++) {
                size_t idx_y = 0;
                while (idx_y < datos_finos.Ny && datos_finos.y[idx_y] < blk_ref.y[i]) idx_y++;
                if (idx_y >= datos_finos.Ny) idx_y = datos_finos.Ny - 1;
        
                for (size_t j = 0; j < nuevo_tam_x; j++) {
                    size_t idx_x = 0;
                    while (idx_x < datos_finos.Nx && datos_finos.x[idx_x] < blk_ref.x[j]) idx_x++;
                    if (idx_x >= datos_finos.Nx) idx_x = datos_finos.Nx - 1;
        
                    blk_ref.data[i * nuevo_tam_x + j] = datos_finos.data[idx_y * datos_finos.Nx + idx_x];
                }
            }
        
            bloques_refinados.push_back(blk_ref);
        
            // Liberar arrays temporales
            delete[] x_red;
            delete[] y_red;
        }
    }

    // Liberar memoria del CSV fino
    delete[] datos_finos.data;
    delete[] datos_finos.x;
    delete[] datos_finos.y;

    return bloques_refinados;
}

// tolerancia para considerar coordenadas iguales (relativa)
static double coord_tol(double range) {
    return max(1e-12, range * 1e-10);
}

Data pegar_refinados(const vector<SubBlock>& bloques_refinados, const Data& datos_superior)
{
    Data result;

    if (bloques_refinados.empty()) {
        // No hay refinamientos → copia directa del nivel superior
        result.Nx = datos_superior.Nx;
        result.Ny = datos_superior.Ny;
        result.x = new double[result.Nx];
        result.y = new double[result.Ny];
        result.data = new double[result.Nx * result.Ny];
        memcpy(result.x, datos_superior.x, result.Nx * sizeof(double));
        memcpy(result.y, datos_superior.y, result.Ny * sizeof(double));
        memcpy(result.data, datos_superior.data, result.Nx * result.Ny * sizeof(double));
        return result;
    }

    size_t i_min = bloques_refinados.front().i0;
    size_t i_max = bloques_refinados.front().i0;
    size_t j_min = bloques_refinados.front().j0;
    size_t j_max = bloques_refinados.front().j0;

    double x_min = 1e30, x_max = -1e30, y_min = 1e30, y_max = -1e30;
    
    for (const auto& blk : bloques_refinados) {
        i_min = min(i_min, blk.i0);
        i_max = max(i_max, blk.i0);
        j_min = min(j_min, blk.j0);
        j_max = max(j_max, blk.j0);
        if (!blk.x || !blk.y) continue;
        x_min = min(x_min, blk.x[0]);
        x_max = max(x_max, blk.x[blk.bloque - 1]);
        y_min = min(y_min, blk.y[0]);
        y_max = max(y_max, blk.y[blk.bloque - 1]);
    }

    size_t bloque = bloques_refinados.front().bloque;

    size_t total_bloques_y = (i_max - i_min + 1);
    size_t total_bloques_x = (j_max - j_min + 1);

    size_t Ny_final = total_bloques_y * bloque;
    size_t Nx_final = total_bloques_x * bloque;

    result.Nx = Nx_final;
    result.Ny = Ny_final;
    result.x = new double[Nx_final];
    result.y = new double[Ny_final];
    result.data = new double[Nx_final * Ny_final];

    linspace(result.x, Nx_final, x_min, x_max);
    linspace(result.y, Ny_final, y_min, y_max);

    // Inicializar datos con interpolación del nivel superior
    for (size_t i = 0; i < Ny_final; ++i) {
        for (size_t j = 0; j < Nx_final; ++j) {
            // Buscar índice más cercano en datos_superior
            size_t iy = lower_bound(datos_superior.y, datos_superior.y + datos_superior.Ny, result.y[i]) - datos_superior.y;
            size_t ix = lower_bound(datos_superior.x, datos_superior.x + datos_superior.Nx, result.x[j]) - datos_superior.x;

            if (iy >= datos_superior.Ny) iy = datos_superior.Ny - 1;
            if (ix >= datos_superior.Nx) ix = datos_superior.Nx - 1;

            result.data[i * Nx_final + j] = datos_superior.data[iy * datos_superior.Nx + ix];
        }
    }


    // Pegar los datos refinados sobre el fondo superior
    for (const auto& blk : bloques_refinados) {
        if (!blk.data) continue;

        size_t i_global = (blk.i0 - i_min) * bloque;
        size_t j_global = (blk.j0 - j_min) * bloque;

        for (size_t bi = 0; bi < bloque; ++bi) {
            size_t row = i_global + bi;
            memcpy(result.data + row * Nx_final + j_global,
                   blk.data + bi * bloque,
                   bloque * sizeof(double));
        }
    }

    cout << "✅ Pegado completado (con datos superiores rellenando huecos)." << endl;
    return result;
}


void liberarData(Data& d) {
    delete[] d.data;
    delete[] d.x;
    delete[] d.y;
    d.data = nullptr;
    d.x = nullptr;
    d.y = nullptr;
    d.Nx = d.Ny = 0;
}

void proceso(int numlevel,
             const vector<int>& numporlevel,
             const vector<int>& ratios,
             const vector<vector<string>>& nombres,
             const vector<vector<int>>& bloques_por_level)
{
    if (numlevel < 2) {
        cerr << "⚠️  Se necesita al menos 2 niveles para anidar." << endl;
        return;
    }

    for (int i = 0; i < numlevel - 1; ++i) {

        cout << "\n=== 🔄 Procesando paso de nivel " << i
             << " → " << i + 1 << " ===" << endl;

        int ratio = ratios[i];
        

        // Iterar sobre todos los subdominios del nivel siguiente
        for (int j = 0; j < numporlevel[i+1]; ++j) {
            
            
            int bloque = bloques_por_level[i][j];
            
            
            cout << "📦 Usando bloque = " << bloque << endl;

            
            // Archivo base y refinado
            string nombre_base = nombres[i][0];        // nivel actual (por ej. nivel0.csv)
            string nombre_ref  = nombres[i+1][j];      // archivo de nivel refinado
            string salida      = nombres[i+1][j];      // sobrescribir mismo archivo

            cout << "🧩 Refinando subbloque " << j + 1
                 << " con ratio " << ratio
                 << "\n   Base: " << nombre_base
                 << "\n   Fino: " << nombre_ref
                 << "\n   → Guardando en: " << salida << endl;

            // === Cargar nivel base ===
            Data grid = devolvercsv(nombre_base);
            if (!grid.data) {
                cerr << "❌ Error al cargar archivo base: " << nombre_base << endl;
                continue;
            }

            // === Crear red y subbloques ===
            Data red = comprimir(grid, bloque);
            vector<SubBlock> bloques = extraerSubbloquesDesdeRed(grid, red, bloque);
            
            vector<SubBlock> refinados = refinar(bloques, nombre_ref, ratio);
            //exportar_subbloques_csv(refinados, "topobatimetrias");

            // === Pegar refinamientos y exportar ===
            
            Data malla = pegar_refinados(refinados, grid);
            exportar_csv(malla, salida, false); // sobrescribe el mismo archivo

            cout << "✅ Guardado refinado en " << salida << endl;

            // === Liberar memoria ===
            for (auto& b : bloques) {
                delete[] b.data;
                delete[] b.x;
                delete[] b.y;
            }
            for (auto& b : refinados) {
                delete[] b.data;
                delete[] b.x;
                delete[] b.y;
            }

            liberarData(grid);
            liberarData(red);
            liberarData(malla);
        }
    }

    cout << "\n🏁 Proceso de anidamiento completado correctamente.\n" << endl;
}


void exportar_csv(const Data& datos, const string& filename, bool print) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "No se pudo abrir el archivo para escribir: " << filename << endl;
        return;
    }

    // Precisión
    file << fixed << setprecision(7);

    // Primera fila con ejes X
    file << "0";
    for (size_t j = 0; j < datos.Nx; ++j)
        file << ',' << datos.x[j];
    file << '\n';

    // Filas con datos y eje Y
    for (size_t i = 0; i < datos.Ny; ++i) {
        file << datos.y[i];
        for (size_t j = 0; j < datos.Nx; ++j)
            file << ',' << datos.data[i * datos.Nx + j];
        file << '\n';
    }

    file.close();
    if (print) cout << "Datos exportados exitosamente a " << filename << endl;
}


void exportar_subbloques_csv(const vector<SubBlock>& bloques, const string& carpeta) {
    // Crear carpeta si no existe
    struct stat info;
    if (stat(carpeta.c_str(), &info) != 0) {
        cout << "📁 Creando carpeta: " << carpeta << endl;
        system(("mkdir -p " + carpeta).c_str());
    }

    // Exportar cada subbloque
    for (size_t idx = 0; idx < bloques.size(); ++idx) {
        const SubBlock& blk = bloques[idx];
        string filename = carpeta + "/subbloque_" + to_string(idx + 1) + ".csv";

        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "❌ No se pudo abrir " << filename << " para escribir." << endl;
            continue;
        }

        file << fixed << setprecision(10);

        // Cabecera con coordenadas X
        file << "0";
        for (size_t j = 0; j < blk.bloque; ++j)
            file << ',' << blk.x[j];
        file << '\n';

        // Filas con Y y datos
        for (size_t i = 0; i < blk.bloque; ++i) {
            file << blk.y[i];
            for (size_t j = 0; j < blk.bloque; ++j)
                file << ',' << blk.data[i * blk.bloque + j];
            file << '\n';
        }

        file.close();
        cout << "✅ Subbloque " << idx + 1 << " exportado a " << filename << endl;
    }

    cout << "🟢 Total de subbloques exportados: " << bloques.size() << endl;
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

// Función para leer el archivo de configuración
Config leer_config(const string& archivo) {
    Config cfg;
    ifstream file(archivo);
    string line;

    if (!file.is_open()) {
        cerr << "Error: No se pudo abrir el archivo '" << archivo
             << "'. Verifica que exista y tenga permisos de lectura." << endl;
        return cfg;
    }

    regex numporlevel_regex("^numporlevel\\s*=\\s*([0-9,]+)");
    regex ratios_regex("^ratios\\s*=\\s*([0-9,]+)");
    regex nombres_start_regex("^nombres\\s*=\\s*$");
    regex bloque_start_regex("^bloque\\s*=\\s*$");

    bool in_nombres_block = false;
    bool in_bloque_block = false;

    while (getline(file, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        smatch match;

        // --- numporlevel ---
        if (regex_search(line, match, numporlevel_regex) && match.size() > 1) {
            vector<string> nums = split(match[1].str(), ',');
            for (const auto& num : nums) {
                try { cfg.numporlevel.push_back(stoi(num)); }
                catch (...) { cerr << "Error al parsear numporlevel: " << num << endl; }
            }
        }

        // --- ratios ---
        else if (regex_search(line, match, ratios_regex) && match.size() > 1) {
            vector<string> nums = split(match[1].str(), ',');
            for (const auto& num : nums) {
                try { cfg.ratios.push_back(stoi(num)); }
                catch (...) { cerr << "Error al parsear ratios: " << num << endl; }
            }
        }

        // --- nombres ---
        else if (regex_search(line, match, nombres_start_regex)) {
            in_nombres_block = true;
            in_bloque_block = false;
            continue;
        } 
        else if (in_nombres_block) {
            if (line.find('=') != string::npos) { // siguiente bloque detectado
                in_nombres_block = false;
                // seguimos, no hacemos continue aquí
            } else {
                vector<string> archivos = split(line, ',');
                vector<string> nombres_limpios;
                for (auto& archivo : archivos) {
                    archivo = trim(archivo);
                    if (!archivo.empty()) nombres_limpios.push_back(archivo);
                }
                if (!nombres_limpios.empty())
                    cfg.nombres.push_back(nombres_limpios);
                continue;
            }
        }

        // --- bloque ---
        if (regex_search(line, match, bloque_start_regex)) {
            in_bloque_block = true;
            in_nombres_block = false;
            continue;
        } 
        else if (in_bloque_block) {
            if (line.find('=') != string::npos) { // final del bloque
                in_bloque_block = false;
            } else {
                vector<string> nums = split(line, ',');
                vector<int> bloques_nivel;
                for (const auto& num : nums) {
                    try { bloques_nivel.push_back(stoi(num)); }
                    catch (...) { cerr << "Error al parsear bloque: " << num << endl; }
                }
                if (!bloques_nivel.empty())
                    cfg.bloque.push_back(bloques_nivel);
                continue;
            }
        }
    }

    file.close();

    // === Validación básica ===
    if (cfg.numporlevel.size() != cfg.nombres.size()) {
        cerr << "Error: numporlevel y nombres tienen distinto número de niveles." << endl;
    }
    if (!cfg.bloque.empty() && cfg.bloque.size() != cfg.nombres.size()) {
        cerr << "Advertencia: bloque no tiene la misma cantidad de niveles que nombres." << endl;
    }

    // === Mostrar resumen ===
    cout << "\n✅ Configuración cargada correctamente:" << endl;
    cout << "  Niveles: " << cfg.numporlevel.size() << endl;
    cout << "  numporlevel: ";
    for (int n : cfg.numporlevel) cout << n << " ";
    cout << "\n  ratios: ";
    for (int r : cfg.ratios) cout << r << " ";
    cout << "\n  nombres:" << endl;
    for (const auto& nivel : cfg.nombres) {
        cout << "   ";
        for (const auto& n : nivel) cout << n << " ";
        cout << endl;
    }
    cout << "  bloque:" << endl;
    for (const auto& nivel : cfg.bloque) {
        cout << "   ";
        for (int b : nivel) cout << b << " ";
        cout << endl;
    }

    return cfg;
}
